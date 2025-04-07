use crate::error::{Error, Mp3Error};
use crate::tables::{LFS_TABLE};
use crate::types::*;
use bitstream_io::{BigEndian, BitReader};
use byteorder::ReadBytesExt;
use std::io::Read;

pub fn read_frame_header<R: Read>(mut data: R) -> Result<FrameHeader, Error> {
    if data.read_u8()? != 0xff {
        panic!("Invalid frame header1");
    }

    let byte = data.read_u8()?;
    if byte & 0b1110_0000 != 0b1110_0000 {
        panic!("Invalid frame header2");
    }

    let version = match byte & 0b0001_1000 {
        0b00_000 => MpegVersion::Mpeg2_5,
        _ => panic!("Invalid MPEG version"),
    };

    let layer = match byte & 0b110 {
        0b010 => MpegLayer::Layer3,
        _ => panic!("Invalid MPEG layer"),
    };

    // CRC is ignored for now.
    let crc = byte & 1 == 0;

    let mut bytes = [0u8; 2];
    data.read_exact(&mut bytes)?;

    let is_version2 = version == MpegVersion::Mpeg2 || version == MpegVersion::Mpeg2_5;
    let bitrate = match (bytes[0] & 0b1111_0000, is_version2) {
        (0b0001_0000, true) => BitRate::Kbps8,
        _ => panic!("Invalid bitrate"),
    };

    let sample_rate = match (bytes[0] & 0b0000_1100, version) {
        (0b10_00, MpegVersion::Mpeg2_5) => SampleRate::Hz8000,
        _ => panic!("Invalid sample rate"),
    };
    let sample_rate_table = ((bytes[0] & 0b0000_1100) >> 2) as usize
        + match version {
            MpegVersion::Mpeg1 => 0,
            MpegVersion::Mpeg2 => 3,
            MpegVersion::Mpeg2_5 => 6,
        };

    let padding = bytes[0] & 0b10 != 0;

    let channels = match bytes[1] & 0b11_000000 {
        0b11_000000 => Channels::Mono,
        _ => panic!("Invalid channel mode"),
    };

    let copyright = bytes[1] & 0b1000 != 0;
    let original = bytes[1] & 0b100 != 0;
    let emphasis = match bytes[1] & 0b11 {
        0b00 => Emphasis::None,
        _ => panic!("Invalid emphasis"),
    };

    if crc {
        // Skip CRC for now.
        data.read_u8()?;
        data.read_u8()?;
    }

    let bits_per_sample = match version {
        MpegVersion::Mpeg2_5 => 72,
        _ => panic!("Invalid bits per sample"),
    };
    let data_size = (bits_per_sample * bitrate.bps() / sample_rate.hz()
        + if padding { 1 } else { 0 }
        - if crc { 2 } else { 0 }
        - 4) as usize;

    // Skip framesize?
    // Skip ancillary data...?

    Ok(FrameHeader {
        version,
        layer,
        crc,
        bitrate,
        sample_rate,
        padding,
        channels,
        copyright,
        original,
        emphasis,

        sample_rate_table,
        data_size,
    })
}

fn read_side_info<R: Read>(mut data: R, header: &FrameHeader) -> Result<SideInfo, Error> {
    let mut info: SideInfo = Default::default();
    let mut bytes = [0u8; 32];
    let size = header.side_data_len();
    data.read_exact(&mut bytes[..size])?;

    let mut reader = BitReader::endian(&bytes[..], BigEndian);

    if header.version == MpegVersion::Mpeg1 {
        info.main_data_begin = reader.read(9)?;

        // Skip private bits.
        if header.channels == Channels::Mono {
            reader.skip(5)?;
        } else {
            reader.skip(3)?;
        }

        for scfsi in &mut info.scfsi[..header.channels.num_channels()] {
            for band in scfsi.iter_mut() {
                *band = reader.read_bit()?;
            }
        }
    } else {
        info.main_data_begin = reader.read(8)?;

        // Skip private bits.
        if header.channels == Channels::Mono {
            reader.skip(1)?;
        } else {
            reader.skip(2)?;
        }
    }

    for granule in &mut info.granules[..header.num_granules()] {
        *granule = read_granule_side_info(&header, &mut reader)?;
    }

    Ok(info)
}

fn read_granule_side_info<R: Read>(
    header: &FrameHeader,
    reader: &mut BitReader<R, BigEndian>,
) -> Result<GranuleSideInfo, Error> {
    let mut info: GranuleSideInfo = Default::default();
    for channel_side_info in &mut info.channels[0..header.channels.num_channels()] {
        *channel_side_info = read_granule_channel_side_info(header, reader)?;
    }
    Ok(info)
}

fn read_granule_channel_side_info<R: Read>(
    header: &FrameHeader,
    reader: &mut BitReader<R, BigEndian>,
) -> Result<GranuleChannelSideInfo, Error> {
    let mut info: GranuleChannelSideInfo = Default::default();

    info.part2_3_length = reader.read(12)?;
    info.big_values = reader.read(9)?;
    if info.big_values > 288 {
        return Err(Error::Mp3Error(Mp3Error::InvalidData("big_values > 288")));
    }
    info.global_gain = reader.read(8)?;
    let scalefac_compress_len = if header.version == MpegVersion::Mpeg1 {
        4
    } else {
        9
    };
    info.scalefac_compress = reader.read(scalefac_compress_len)?;

    let window_switching = reader.read_bit()?;
    if window_switching {
        let block_type_id = reader.read::<u8>(2)?;
        let mixed_block = reader.read_bit()?;
        for region in &mut info.table_select[..2] {
            *region = reader.read(5)?;
        }

        let mut subblock_gain = [0f32; 3];
        for gain in &mut subblock_gain {
            *gain = reader.read::<u8>(3)?.into();
        }
        info.subblock_gain = subblock_gain;

        info.block_type = match block_type_id {
            0b00 => {
                // Block type 00 is only if window switching is off
                return Err(Error::Mp3Error(Mp3Error::InvalidData(
                    "Forbidden block type",
                )));
            }
            0b01 => BlockType::Start,
            0b10 => {
                if mixed_block {
                    BlockType::Mixed
                } else {
                    BlockType::Short
                }
            }
            0b11 => BlockType::End,
            _ => unreachable!(),
        };

        // Mixed blocks are always marked as short.
        assert!(!mixed_block || info.block_type == BlockType::Short);

        info.region0_count = if info.block_type == BlockType::Short {
            8
        } else {
            7
        };
        info.region1_count = 20 - info.region0_count;
    } else {
        info.block_type = BlockType::Long;

        for region in &mut info.table_select {
            *region = reader.read(5)?;
        }

        info.region0_count = reader.read(4)?;
        info.region1_count = reader.read(3)?;
    }

    info.preflag = if header.version == MpegVersion::Mpeg1 {
        reader.read_bit()?
    } else {
        info.scalefac_compress >= 500
    };

    info.scalefac_scale = reader.read_bit()?; // .5f * (1f + frame.ReadBits(1));
    info.count1table_select = reader.read_bit()?;

    Ok(info)
}

fn read_logical_frame_data<'a, R: Read>(
    decoder: &'a mut DecoderState,
    mut reader: R,
    header: &FrameHeader,
    side_info: &SideInfo,
) -> Result<&'a [u8], Error> {
    let side_info_size = header.side_data_len();
    let main_data_size = header.data_size - side_info_size;

    // Copy main_data_begin bytes from the previous frame(s).
    let main_data_begin = side_info.main_data_begin as usize;
    let prev_start = decoder.frame_buffer_len - main_data_begin;
    for i in 0..main_data_begin {
        decoder.frame_buffer[i] = decoder.frame_buffer[prev_start + i];
    }
    decoder.frame_buffer_len = main_data_begin + main_data_size;
    reader.read_exact(&mut decoder.frame_buffer[main_data_begin..decoder.frame_buffer_len])?;

    Ok(&decoder.frame_buffer[0..decoder.frame_buffer_len])
}

fn read_main_data<R: Read>(
    reader: &mut BitReader<R, BigEndian>,
    header: &FrameHeader,
    side_info: &SideInfo,
) -> Result<MainData, Error> {
    let mut data: MainData = Default::default();

    for g in 0..header.num_granules() {
        for c in 0..header.channels.num_channels() {
            let bits_read = if header.version == MpegVersion::Mpeg1 {
                panic!("Mpeg1 not supported");
            } else {
                read_lfs_scale_factors(
                    reader,
                    c == 1 && header.is_intensity_stereo(),
                    &side_info.granules[g].channels[c],
                    &mut data.granules[g].channels[c],
                )?
            };

            let huffman_len =
                u32::from(side_info.granules[g].channels[c].part2_3_length) - bits_read;
            data.granules[g].channels[c].count1 = crate::huffman::read_huffman(
                reader,
                header,
                &side_info.granules[g].channels[c],
                huffman_len,
                &mut data.granules[g].channels[c].samples,
            )?;
        }
    }

    // TODO(Herschel): Skip ancillary data.
    Ok(data)
}

fn read_lfs_scale_factors<R: Read>(
    reader: &mut BitReader<R, BigEndian>,
    intensity_stereo_channel: bool,
    channel_info: &GranuleChannelSideInfo,
    channel_data: &mut MainDataChannel,
) -> Result<u32, Error> {
    let mut bits_read = 0;

    let lfs_table = if intensity_stereo_channel {
        panic!("Intensity stereo not supported");
    } else {
        &LFS_TABLE
    };
    let lfs_table = match channel_info.block_type {
        BlockType::Short => &lfs_table[1],
        BlockType::Mixed => &lfs_table[2],
        _ => &lfs_table[0],
    };

    let (scale_lens, lfs_table) = if intensity_stereo_channel {
        let sfc = u32::from(channel_info.scalefac_compress / 2);
        match sfc {
            0..=179 => ([sfc / 36, (sfc % 36) / 6, sfc % 6, 0], &lfs_table[0]),
            180..=243 => (
                [
                    ((sfc - 180) % 64) / 16,
                    ((sfc - 180) % 16) / 4,
                    (sfc - 180) % 4,
                    0,
                ],
                &lfs_table[1],
            ),
            244..=255 => ([(sfc - 244) / 3, (sfc - 244) % 3, 0, 0], &lfs_table[2]),
            _ => unreachable!(),
        }
    } else {
        let sfc = u32::from(channel_info.scalefac_compress);
        match sfc {
            0..=399 => (
                [sfc / 80, (sfc / 16) % 5, (sfc % 16) / 4, sfc & 3],
                &lfs_table[0],
            ),
            400..=499 => (
                [(sfc - 400) / 20, ((sfc - 400) / 4) % 5, (sfc - 400) % 4, 0],
                &lfs_table[1],
            ),
            500..=512 => ([(sfc - 500) / 3, (sfc - 500) % 3, 0, 0], &lfs_table[2]),
            _ => unreachable!(),
        }
    };

    // TODO(Herschel): We could avoid using this intermediate buffer.
    // Write an iterator for reading scalefacs and/or write an iterator
    // through scalefac_s/l for the block type.
    let mut scalefacs = [0u8; 54];
    let mut i = 0;
    for (&len, &num_blocks) in scale_lens[..].iter().zip(lfs_table.iter()) {
        assert!(len <= 8);
        if len > 0 {
            for _ in 0..num_blocks {
                scalefacs[i] = reader.read(len)?;
                bits_read += len;
                i += 1;
            }
        } else {
            i += num_blocks;
        }
    }

    i = 0;
    if channel_info.block_type == BlockType::Short || channel_info.block_type == BlockType::Mixed {
        let short_start = if channel_info.block_type == BlockType::Mixed {
            for sfb in 0..8 {
                channel_data.scalefac_l[sfb] = scalefacs[i];
                i += 1;
            }
            3
        } else {
            0
        };

        for sfb in short_start..12 {
            for window in 0..3 {
                channel_data.scalefac_s[sfb][window] = scalefacs[i];
                i += 1;
            }
        }
    } else {
        for sfb in 0..21 {
            channel_data.scalefac_l[sfb] = scalefacs[i];
            i += 1;
        }
    }

    Ok(bits_read)
}

pub fn process_frame<R: Read>(
    decoder: &mut DecoderState,
    mut reader: R,
    header: &FrameHeader,
) -> Result<(usize, [[f32; 1152]; 2]), Error> {
    let side_info = read_side_info(&mut reader, header)?;
    let data_buffer = read_logical_frame_data(decoder, &mut reader, header, &side_info)?;

    let mut reader = BitReader::endian(data_buffer, BigEndian);
    let mut main_data = read_main_data(&mut reader, header, &side_info)?;

    let mut out_samples = [[0f32; 1152]; 2];
    let num_samples = decode_frame(
        decoder,
        header,
        &side_info,
        &mut main_data,
        &mut out_samples,
    )?;

    Ok((num_samples, out_samples))
}

fn decode_frame(
    decoder: &mut DecoderState,
    header: &FrameHeader,
    side_info: &SideInfo,
    main_data: &mut MainData,
    out_samples: &mut [[f32; 1152]; 2],
) -> Result<usize, Error> {
    use crate::{requantize, synthesis};

    if header.channels == Channels::Mono {
        for gr in 0..header.num_granules() {
            let side_info = &side_info.granules[gr].channels[0];
            let main_data = &mut main_data.granules[gr].channels[0];

            requantize::requantize(header, side_info, main_data);
            requantize::reorder(header, side_info, main_data);
            synthesis::antialias(side_info, &mut main_data.samples);
            synthesis::hybrid_synthesis(
                side_info.block_type,
                &mut decoder.store[0],
                &mut main_data.samples,
            );
            synthesis::frequency_inversion(&mut main_data.samples);
            synthesis::subband_synthesis(
                &main_data.samples,
                &mut decoder.sbs_v_vec[0],
                &mut out_samples[0][gr * 576..(gr + 1) * 576],
            );
        }

        out_samples[1] = out_samples[0];
    } else {
        panic!("Stereo not supported");
    }
    Ok(header.num_granules() * 576)
}
