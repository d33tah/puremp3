//! An MP3 decoder implemented in pure Rust.
//!
//! Supports MPEG-1, MPEG-2, and MPEG-2.5 Layer III streams.
//! Layers I and II are currently unsupported.
//!
//! # Example
//!
//! ```
//! let data = std::fs::read("tests/vectors/MonoCBR192.mp3").expect("Could not open file");
//! let (header, samples) = puremp3::read_mp3(&data[..]).expect("Invalid MP3");
//! for (left, right) in samples {
//!     // Operate on samples here
//! }
//! ```

///! Error types related to MP3 decoding.
use std::{fmt, io};
use bitstream_io::{BigEndian, BitReader};
use byteorder::ReadBytesExt;
use std::io::Read;


/// Error that can be raised during MP3 decoding.
#[derive(Debug)]
pub enum Error {
    /// An error during the MP3 decoding process.
    Mp3Error(Mp3Error),

    // An IO error reading the underlying stream.
    IoError(io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Mp3Error(e) => write!(f, "MP3 Error: {}", e),
            Error::IoError(e) => write!(f, "IO Error: {}", e),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Mp3Error(e) => Some(e),
            Error::IoError(e) => Some(e),
        }
    }
}

#[derive(Debug)]
pub enum Mp3Error {
    /// Invalid or unknown data was encountered when reading the stream.
    InvalidData(&'static str),

    /// An unsupported MP3 feature is used in this MP3 stream.
    Unsupported(&'static str),
}

impl fmt::Display for Mp3Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Mp3Error::InvalidData(s) => write!(f, "Invalid data: {}", s),
            Mp3Error::Unsupported(s) => write!(f, "Unsupported: {}", s),
        }
    }
}

impl std::error::Error for Mp3Error {}

impl From<Mp3Error> for Error {
    fn from(error: Mp3Error) -> Self {
        Error::Mp3Error(error)
    }
}

impl From<io::Error> for Error {
    fn from(error: io::Error) -> Self {
        Error::IoError(error)
    }
}

/// The maximum number of channels supported in an MP3.
pub const MAX_CHANNELS: usize = 2;

/// The maximum number of granules in an MP3 frame.AsMut
///
/// Depends on the MPEG version.
pub(crate) const MAX_GRANULES: usize = 2;

/// Header of an MP3 frame.
///
/// Contains info about the format of the audio samples.
#[derive(Debug, Clone)]
pub struct FrameHeader {
    /// The MPEG standard used in encoding this frame.
    pub version: MpegVersion,

    /// The MPEG layer of the frame.
    ///
    /// Currently only MPEG Layer III is supported.
    pub layer: MpegLayer,

    /// Whether the frame contains a CRC checksum.
    pub crc: bool,

    /// The bitrate of this frame.
    pub bitrate: BitRate,

    /// The sample rate of this frame.
    pub sample_rate: SampleRate,

    /// Whether the frame has an extra padding bit.
    pub padding: bool,

    /// The channel mode of this frame.
    pub channels: Channels,

    /// Whether this frame is under copyright.
    pub copyright: bool,

    /// Whether this frame contains original data or a copy.
    pub original: bool,

    /// The emphasis of this frame.
    pub emphasis: Emphasis,

    pub(crate) sample_rate_table: usize,
    pub(crate) data_size: usize,
}

impl FrameHeader {
    pub(crate) fn side_data_len(&self) -> usize {
        match self.layer {
            MpegLayer::Layer3 => {
                if self.channels == Channels::Mono && self.version != MpegVersion::Mpeg1 {
                    9
                } else if self.channels != Channels::Mono && self.version == MpegVersion::Mpeg1 {
                    32
                } else {
                    17
                }
            }
            _ => unimplemented!(),
        }
    }

    pub(crate) fn num_granules(&self) -> usize {
        if self.version == MpegVersion::Mpeg1 {
            2
        } else {
            1
        }
    }

    pub(crate) fn is_intensity_stereo(&self) -> bool {
        if let Channels::JointStereo {
            intensity_stereo: true,
            ..
        } = self.channels
        {
            true
        } else {
            false
        }
    }
}

/// The version of the MPEG standard used in encoding audio.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
#[allow(clippy::enum_variant_names)]
pub enum MpegVersion {
    /// MPEG-1 (ISO/IEC 11172-3)
    Mpeg1,

    /// MPEG-2 (ISO/IEC 13818-3)
    Mpeg2,

    /// MPEG-2.5
    Mpeg2_5,
}

/// The MPEG Layer used in encoding audio.
///
/// Higher layers provide better compression, but are more complex.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
#[allow(clippy::enum_variant_names)]
pub enum MpegLayer {
    /// MPEG Layer I
    Layer1,

    /// MPEG Layer II
    Layer2,

    /// MPEG Layer III
    Layer3,
}

/// The channel mode
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Channels {
    /// One audio channel.
    Mono,

    /// Two unrelated audio channels (e.g. for different languages).
    DualMono,

    /// Stereo.
    Stereo,

    /// Joint stereo. Improves compression by utilizing the correlation
    /// in stereo channels.
    JointStereo {
        intensity_stereo: bool,
        mid_side_stereo: bool,
    },
}

impl Channels {
    /// The number of audio channels.
    pub fn num_channels(self) -> usize {
        match self {
            Channels::Mono => 1,
            _ => 2,
        }
    }
}

/// The bit rate of an MP3 stream.
///
/// MP3 supports specific bitrates, depending on the MPEG version and layer.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum BitRate {
    Kbps8,
    Kbps16,
    Kbps24,
    Kbps32,
    Kbps40,
    Kbps48,
    Kbps56,
    Kbps64,
    Kbps80,
    Kbps96,
    Kbps112,
    Kbps128,
    Kbps144,
    Kbps160,
    Kbps192,
    Kbps224,
    Kbps256,
    Kbps320,
}

impl BitRate {
    /// Returns the bit rate in bits per second as a `u32`.
    pub fn bps(self) -> u32 {
        match self {
            BitRate::Kbps8 => 8_000,
            BitRate::Kbps16 => 16_000,
            BitRate::Kbps24 => 24_000,
            BitRate::Kbps32 => 32_000,
            BitRate::Kbps40 => 40_000,
            BitRate::Kbps48 => 48_000,
            BitRate::Kbps56 => 56_000,
            BitRate::Kbps64 => 64_000,
            BitRate::Kbps80 => 80_000,
            BitRate::Kbps96 => 96_000,
            BitRate::Kbps112 => 112_000,
            BitRate::Kbps128 => 128_000,
            BitRate::Kbps144 => 144_000,
            BitRate::Kbps160 => 160_000,
            BitRate::Kbps192 => 192_000,
            BitRate::Kbps224 => 224_000,
            BitRate::Kbps256 => 256_000,
            BitRate::Kbps320 => 320_000,
        }
    }
}

/// The sample rate of an MP3 stream.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum SampleRate {
    Hz8000,
    Hz11025,
    Hz12000,
    Hz16000,
    Hz22050,
    Hz24000,
    Hz32000,
    Hz44100,
    Hz48000,
}

impl SampleRate {
    /// Returns the sample rate in hertz as a `u32`.
    pub fn hz(self) -> u32 {
        match self {
            SampleRate::Hz8000 => 8_000,
            SampleRate::Hz11025 => 11_025,
            SampleRate::Hz12000 => 12_000,
            SampleRate::Hz16000 => 16_000,
            SampleRate::Hz22050 => 22_050,
            SampleRate::Hz24000 => 24_000,
            SampleRate::Hz32000 => 32_000,
            SampleRate::Hz44100 => 44_100,
            SampleRate::Hz48000 => 48_000,
        }
    }
}

/// Emphasis used in encoding an MP3 audio stream.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Emphasis {
    None,
    FiftyFifteen,
    CcitJ17,
}

// Internal types
pub struct DecoderState {
    pub frame_buffer: [u8; 4096],
    pub frame_buffer_len: usize,
    pub store: [[[f32; 18]; 32]; 2],
    pub sbs_v_vec: [[f32; 1024]; 2],
}

impl Default for DecoderState {
    fn default() -> Self {
        Self::new()
    }
}

impl DecoderState {
    pub fn new() -> Self {
        DecoderState {
            frame_buffer: [0; 4096],
            frame_buffer_len: 0,
            store: [[[0f32; 18]; 32]; 2],
            sbs_v_vec: [[0f32; 1024]; 2],
        }
    }
}

#[derive(Debug, Default)]
pub struct SideInfo {
    pub main_data_begin: u16,
    pub granules: [GranuleSideInfo; 2],
}

#[derive(Debug, Default)]
pub struct GranuleSideInfo {
    pub channels: [GranuleChannelSideInfo; 2],
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
#[derive(Default)]
pub enum BlockType {
    #[default]
    Long,
    Short,
    Mixed,
    Start,
    End,
}


#[derive(Debug, Default)]
pub struct GranuleChannelSideInfo {
    pub part2_3_length: u16,
    pub big_values: u16,
    pub global_gain: u8,
    pub scalefac_compress: u16,
    pub block_type: BlockType,
    pub subblock_gain: [f32; 3],

    pub table_select: [u8; 3],
    pub region0_count: u8,
    pub region1_count: u8,
    pub preflag: bool,
    pub scalefac_scale: bool,
    pub count1table_select: bool,
}

#[derive(Debug, Default)]
pub struct MainData {
    pub granules: [MainDataGranule; MAX_GRANULES],
}

#[derive(Debug, Default)]
pub struct MainDataGranule {
    pub channels: [MainDataChannel; MAX_CHANNELS],
}

pub struct MainDataChannel {
    pub scalefac_l: [u8; 22],
    pub scalefac_s: [[u8; 3]; 13],
    pub count1: u32, // TODO(Herschel): What's the actual size of this?
    pub samples: [f32; 576],
}

impl Default for MainDataChannel {
    fn default() -> Self {
        Self {
            scalefac_l: Default::default(),
            scalefac_s: Default::default(),
            count1: Default::default(),
            samples: [Default::default(); 576],
        }
    }
}

impl std::fmt::Debug for MainDataChannel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "MainDataChannel")
    }
}

// (l, s)
const SCALE_FACTOR_BAND_INDICES: [([u32; 23], [u32; 14]); 9] = [
    // MPEG-1
    (
        [
            0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 52, 62, 74, 90, 110, 134, 162, 196, 238, 288, 342,
            418, 576,
        ],
        [0, 4, 8, 12, 16, 22, 30, 40, 52, 66, 84, 106, 136, 192],
    ),
    (
        [
            0, 4, 8, 12, 16, 20, 24, 30, 36, 42, 50, 60, 72, 88, 106, 128, 156, 190, 230, 276, 330,
            384, 576,
        ],
        [0, 4, 8, 12, 16, 22, 28, 38, 50, 64, 80, 100, 126, 192],
    ),
    (
        [
            0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 54, 66, 82, 102, 126, 156, 194, 240, 296, 364,
            448, 550, 576,
        ],
        [0, 4, 8, 12, 16, 22, 30, 42, 58, 78, 104, 138, 180, 192],
    ),
    // MPEG-2
    (
        [
            0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 52, 62, 74, 90, 110, 134, 162, 196, 238, 288, 342,
            418, 576,
        ],
        [0, 4, 8, 12, 16, 22, 30, 40, 52, 66, 84, 106, 136, 192],
    ),
    (
        [
            0, 4, 8, 12, 16, 20, 24, 30, 36, 42, 50, 60, 72, 88, 106, 128, 156, 190, 230, 276, 330,
            384, 576,
        ],
        [0, 4, 8, 12, 16, 22, 28, 38, 50, 64, 80, 100, 126, 192],
    ),
    (
        [
            0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 54, 66, 82, 102, 126, 156, 194, 240, 296, 364,
            448, 550, 576,
        ],
        [0, 4, 8, 12, 16, 22, 30, 42, 58, 78, 104, 138, 180, 192],
    ),
    // MPEG-2.5
    (
        [
            0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396,
            464, 522, 576,
        ],
        [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192],
    ),
    (
        [
            0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396,
            464, 522, 576,
        ],
        [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192],
    ),
    (
        [
            0, 12, 24, 36, 48, 60, 72, 88, 108, 132, 160, 192, 232, 280, 336, 400, 476, 566, 568,
            570, 572, 574, 576,
        ],
        [0, 8, 16, 24, 36, 52, 72, 96, 124, 160, 162, 164, 166, 192],
    ),
];

const LFS_TABLE: [[[usize; 4]; 3]; 3] = [
    [[6, 5, 5, 5], [6, 5, 7, 3], [11, 10, 0, 0]],
    [[9, 9, 9, 9], [9, 9, 12, 6], [18, 18, 0, 0]],
    [[6, 9, 9, 9], [6, 9, 12, 6], [15, 18, 0, 0]],
];

#[derive(Debug)]
pub struct HuffmanTable {
    pub data: &'static [u16],
    pub linbits: usize,
    pub quads: bool,
}

const HUFFMAN_TABLES: [HuffmanTable; 34] = [
    // Table 0
    HuffmanTable {
        data: &[],
        linbits: 0,
        quads: false,
    },
    // Table 1
    HuffmanTable {
        data: &[0x0201, 0x0000, 0x0201, 0x0010, 0x0201, 0x0001, 0x0011],
        linbits: 0,
        quads: false,
    },
    // Table 2
    HuffmanTable {
        data: &[
            0x0201, 0x0000, 0x0401, 0x0201, 0x0010, 0x0001, 0x0201, 0x0011, 0x0401, 0x0201, 0x0020,
            0x0021, 0x0201, 0x0012, 0x0201, 0x0002, 0x0022,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 3
    HuffmanTable {
        data: &[
            0x0401, 0x0201, 0x0000, 0x0001, 0x0201, 0x0011, 0x0201, 0x0010, 0x0401, 0x0201, 0x0020,
            0x0021, 0x0201, 0x0012, 0x0201, 0x0002, 0x0022,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 4
    HuffmanTable {
        data: &[],
        linbits: 0,
        quads: false,
    },
    // Table 5
    HuffmanTable {
        data: &[
            0x0201, 0x0000, 0x0401, 0x0201, 0x0010, 0x0001, 0x0201, 0x0011, 0x0801, 0x0401, 0x0201,
            0x0020, 0x0002, 0x0201, 0x0021, 0x0012, 0x0801, 0x0401, 0x0201, 0x0022, 0x0030, 0x0201,
            0x0003, 0x0013, 0x0201, 0x0031, 0x0201, 0x0032, 0x0201, 0x0023, 0x0033,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 6
    HuffmanTable {
        data: &[
            0x0601, 0x0401, 0x0201, 0x0000, 0x0010, 0x0011, 0x0601, 0x0201, 0x0001, 0x0201, 0x0020,
            0x0021, 0x0601, 0x0201, 0x0012, 0x0201, 0x0002, 0x0022, 0x0401, 0x0201, 0x0031, 0x0013,
            0x0401, 0x0201, 0x0030, 0x0032, 0x0201, 0x0023, 0x0201, 0x0003, 0x0033,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 7
    HuffmanTable {
        data: &[
            0x0201, 0x0000, 0x0401, 0x0201, 0x0010, 0x0001, 0x0801, 0x0201, 0x0011, 0x0401, 0x0201,
            0x0020, 0x0002, 0x0021, 0x1201, 0x0601, 0x0201, 0x0012, 0x0201, 0x0022, 0x0030, 0x0401,
            0x0201, 0x0031, 0x0013, 0x0401, 0x0201, 0x0003, 0x0032, 0x0201, 0x0023, 0x0004, 0x0a01,
            0x0401, 0x0201, 0x0040, 0x0041, 0x0201, 0x0014, 0x0201, 0x0042, 0x0024, 0x0c01, 0x0601,
            0x0401, 0x0201, 0x0033, 0x0043, 0x0050, 0x0401, 0x0201, 0x0034, 0x0005, 0x0051, 0x0601,
            0x0201, 0x0015, 0x0201, 0x0052, 0x0025, 0x0401, 0x0201, 0x0044, 0x0035, 0x0401, 0x0201,
            0x0053, 0x0054, 0x0201, 0x0045, 0x0055,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 8
    HuffmanTable {
        data: &[
            0x0601, 0x0201, 0x0000, 0x0201, 0x0010, 0x0001, 0x0201, 0x0011, 0x0401, 0x0201, 0x0021,
            0x0012, 0x0e01, 0x0401, 0x0201, 0x0020, 0x0002, 0x0201, 0x0022, 0x0401, 0x0201, 0x0030,
            0x0003, 0x0201, 0x0031, 0x0013, 0x0e01, 0x0801, 0x0401, 0x0201, 0x0032, 0x0023, 0x0201,
            0x0040, 0x0004, 0x0201, 0x0041, 0x0201, 0x0014, 0x0042, 0x0c01, 0x0601, 0x0201, 0x0024,
            0x0201, 0x0033, 0x0050, 0x0401, 0x0201, 0x0043, 0x0034, 0x0051, 0x0601, 0x0201, 0x0015,
            0x0201, 0x0005, 0x0052, 0x0601, 0x0201, 0x0025, 0x0201, 0x0044, 0x0035, 0x0201, 0x0053,
            0x0201, 0x0045, 0x0201, 0x0054, 0x0055,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 9
    HuffmanTable {
        data: &[
            0x0801, 0x0401, 0x0201, 0x0000, 0x0010, 0x0201, 0x0001, 0x0011, 0x0a01, 0x0401, 0x0201,
            0x0020, 0x0021, 0x0201, 0x0012, 0x0201, 0x0002, 0x0022, 0x0c01, 0x0601, 0x0401, 0x0201,
            0x0030, 0x0003, 0x0031, 0x0201, 0x0013, 0x0201, 0x0032, 0x0023, 0x0c01, 0x0401, 0x0201,
            0x0041, 0x0014, 0x0401, 0x0201, 0x0040, 0x0033, 0x0201, 0x0042, 0x0024, 0x0a01, 0x0601,
            0x0401, 0x0201, 0x0004, 0x0050, 0x0043, 0x0201, 0x0034, 0x0051, 0x0801, 0x0401, 0x0201,
            0x0015, 0x0052, 0x0201, 0x0025, 0x0044, 0x0601, 0x0401, 0x0201, 0x0005, 0x0054, 0x0053,
            0x0201, 0x0035, 0x0201, 0x0045, 0x0055,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 10
    HuffmanTable {
        data: &[
            0x0201, 0x0000, 0x0401, 0x0201, 0x0010, 0x0001, 0x0a01, 0x0201, 0x0011, 0x0401, 0x0201,
            0x0020, 0x0002, 0x0201, 0x0021, 0x0012, 0x1c01, 0x0801, 0x0401, 0x0201, 0x0022, 0x0030,
            0x0201, 0x0031, 0x0013, 0x0801, 0x0401, 0x0201, 0x0003, 0x0032, 0x0201, 0x0023, 0x0040,
            0x0401, 0x0201, 0x0041, 0x0014, 0x0401, 0x0201, 0x0004, 0x0033, 0x0201, 0x0042, 0x0024,
            0x1c01, 0x0a01, 0x0601, 0x0401, 0x0201, 0x0050, 0x0005, 0x0060, 0x0201, 0x0061, 0x0016,
            0x0c01, 0x0601, 0x0401, 0x0201, 0x0043, 0x0034, 0x0051, 0x0201, 0x0015, 0x0201, 0x0052,
            0x0025, 0x0401, 0x0201, 0x0026, 0x0036, 0x0071, 0x1401, 0x0801, 0x0201, 0x0017, 0x0401,
            0x0201, 0x0044, 0x0053, 0x0006, 0x0601, 0x0401, 0x0201, 0x0035, 0x0045, 0x0062, 0x0201,
            0x0070, 0x0201, 0x0007, 0x0064, 0x0e01, 0x0401, 0x0201, 0x0072, 0x0027, 0x0601, 0x0201,
            0x0063, 0x0201, 0x0054, 0x0055, 0x0201, 0x0046, 0x0073, 0x0801, 0x0401, 0x0201, 0x0037,
            0x0065, 0x0201, 0x0056, 0x0074, 0x0601, 0x0201, 0x0047, 0x0201, 0x0066, 0x0075, 0x0401,
            0x0201, 0x0057, 0x0076, 0x0201, 0x0067, 0x0077,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 11
    HuffmanTable {
        data: &[
            0x0601, 0x0201, 0x0000, 0x0201, 0x0010, 0x0001, 0x0801, 0x0201, 0x0011, 0x0401, 0x0201,
            0x0020, 0x0002, 0x0012, 0x1801, 0x0801, 0x0201, 0x0021, 0x0201, 0x0022, 0x0201, 0x0030,
            0x0003, 0x0401, 0x0201, 0x0031, 0x0013, 0x0401, 0x0201, 0x0032, 0x0023, 0x0401, 0x0201,
            0x0040, 0x0004, 0x0201, 0x0041, 0x0014, 0x1e01, 0x1001, 0x0a01, 0x0401, 0x0201, 0x0042,
            0x0024, 0x0401, 0x0201, 0x0033, 0x0043, 0x0050, 0x0401, 0x0201, 0x0034, 0x0051, 0x0061,
            0x0601, 0x0201, 0x0016, 0x0201, 0x0006, 0x0026, 0x0201, 0x0062, 0x0201, 0x0015, 0x0201,
            0x0005, 0x0052, 0x1001, 0x0a01, 0x0601, 0x0401, 0x0201, 0x0025, 0x0044, 0x0060, 0x0201,
            0x0063, 0x0036, 0x0401, 0x0201, 0x0070, 0x0017, 0x0071, 0x1001, 0x0601, 0x0401, 0x0201,
            0x0007, 0x0064, 0x0072, 0x0201, 0x0027, 0x0401, 0x0201, 0x0053, 0x0035, 0x0201, 0x0054,
            0x0045, 0x0a01, 0x0401, 0x0201, 0x0046, 0x0073, 0x0201, 0x0037, 0x0201, 0x0065, 0x0056,
            0x0a01, 0x0601, 0x0401, 0x0201, 0x0055, 0x0057, 0x0074, 0x0201, 0x0047, 0x0066, 0x0401,
            0x0201, 0x0075, 0x0076, 0x0201, 0x0067, 0x0077,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 12
    HuffmanTable {
        data: &[
            0x0c01, 0x0401, 0x0201, 0x0010, 0x0001, 0x0201, 0x0011, 0x0201, 0x0000, 0x0201, 0x0020,
            0x0002, 0x1001, 0x0401, 0x0201, 0x0021, 0x0012, 0x0401, 0x0201, 0x0022, 0x0031, 0x0201,
            0x0013, 0x0201, 0x0030, 0x0201, 0x0003, 0x0040, 0x1a01, 0x0801, 0x0401, 0x0201, 0x0032,
            0x0023, 0x0201, 0x0041, 0x0033, 0x0a01, 0x0401, 0x0201, 0x0014, 0x0042, 0x0201, 0x0024,
            0x0201, 0x0004, 0x0050, 0x0401, 0x0201, 0x0043, 0x0034, 0x0201, 0x0051, 0x0015, 0x1c01,
            0x0e01, 0x0801, 0x0401, 0x0201, 0x0052, 0x0025, 0x0201, 0x0053, 0x0035, 0x0401, 0x0201,
            0x0060, 0x0016, 0x0061, 0x0401, 0x0201, 0x0062, 0x0026, 0x0601, 0x0401, 0x0201, 0x0005,
            0x0006, 0x0044, 0x0201, 0x0054, 0x0045, 0x1201, 0x0a01, 0x0401, 0x0201, 0x0063, 0x0036,
            0x0401, 0x0201, 0x0070, 0x0007, 0x0071, 0x0401, 0x0201, 0x0017, 0x0064, 0x0201, 0x0046,
            0x0072, 0x0a01, 0x0601, 0x0201, 0x0027, 0x0201, 0x0055, 0x0073, 0x0201, 0x0037, 0x0056,
            0x0801, 0x0401, 0x0201, 0x0065, 0x0074, 0x0201, 0x0047, 0x0066, 0x0401, 0x0201, 0x0075,
            0x0057, 0x0201, 0x0076, 0x0201, 0x0067, 0x0077,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 13
    HuffmanTable {
        data: &[
            0x0201, 0x0000, 0x0601, 0x0201, 0x0010, 0x0201, 0x0001, 0x0011, 0x1c01, 0x0801, 0x0401,
            0x0201, 0x0020, 0x0002, 0x0201, 0x0021, 0x0012, 0x0801, 0x0401, 0x0201, 0x0022, 0x0030,
            0x0201, 0x0003, 0x0031, 0x0601, 0x0201, 0x0013, 0x0201, 0x0032, 0x0023, 0x0401, 0x0201,
            0x0040, 0x0004, 0x0041, 0x4601, 0x1c01, 0x0e01, 0x0601, 0x0201, 0x0014, 0x0201, 0x0033,
            0x0042, 0x0401, 0x0201, 0x0024, 0x0050, 0x0201, 0x0043, 0x0034, 0x0401, 0x0201, 0x0051,
            0x0015, 0x0401, 0x0201, 0x0005, 0x0052, 0x0201, 0x0025, 0x0201, 0x0044, 0x0053, 0x0e01,
            0x0801, 0x0401, 0x0201, 0x0060, 0x0006, 0x0201, 0x0061, 0x0016, 0x0401, 0x0201, 0x0080,
            0x0008, 0x0081, 0x1001, 0x0801, 0x0401, 0x0201, 0x0035, 0x0062, 0x0201, 0x0026, 0x0054,
            0x0401, 0x0201, 0x0045, 0x0063, 0x0201, 0x0036, 0x0070, 0x0601, 0x0401, 0x0201, 0x0007,
            0x0055, 0x0071, 0x0201, 0x0017, 0x0201, 0x0027, 0x0037, 0x4801, 0x1801, 0x0c01, 0x0401,
            0x0201, 0x0018, 0x0082, 0x0201, 0x0028, 0x0401, 0x0201, 0x0064, 0x0046, 0x0072, 0x0801,
            0x0401, 0x0201, 0x0084, 0x0048, 0x0201, 0x0090, 0x0009, 0x0201, 0x0091, 0x0019, 0x1801,
            0x0e01, 0x0801, 0x0401, 0x0201, 0x0073, 0x0065, 0x0201, 0x0056, 0x0074, 0x0401, 0x0201,
            0x0047, 0x0066, 0x0083, 0x0601, 0x0201, 0x0038, 0x0201, 0x0075, 0x0057, 0x0201, 0x0092,
            0x0029, 0x0e01, 0x0801, 0x0401, 0x0201, 0x0067, 0x0085, 0x0201, 0x0058, 0x0039, 0x0201,
            0x0093, 0x0201, 0x0049, 0x0086, 0x0601, 0x0201, 0x00a0, 0x0201, 0x0068, 0x000a, 0x0201,
            0x00a1, 0x001a, 0x4401, 0x1801, 0x0c01, 0x0401, 0x0201, 0x00a2, 0x002a, 0x0401, 0x0201,
            0x0095, 0x0059, 0x0201, 0x00a3, 0x003a, 0x0801, 0x0401, 0x0201, 0x004a, 0x0096, 0x0201,
            0x00b0, 0x000b, 0x0201, 0x00b1, 0x001b, 0x1401, 0x0801, 0x0201, 0x00b2, 0x0401, 0x0201,
            0x0076, 0x0077, 0x0094, 0x0601, 0x0401, 0x0201, 0x0087, 0x0078, 0x00a4, 0x0401, 0x0201,
            0x0069, 0x00a5, 0x002b, 0x0c01, 0x0601, 0x0401, 0x0201, 0x005a, 0x0088, 0x00b3, 0x0201,
            0x003b, 0x0201, 0x0079, 0x00a6, 0x0601, 0x0401, 0x0201, 0x006a, 0x00b4, 0x00c0, 0x0401,
            0x0201, 0x000c, 0x0098, 0x00c1, 0x3c01, 0x1601, 0x0a01, 0x0601, 0x0201, 0x001c, 0x0201,
            0x0089, 0x00b5, 0x0201, 0x005b, 0x00c2, 0x0401, 0x0201, 0x002c, 0x003c, 0x0401, 0x0201,
            0x00b6, 0x006b, 0x0201, 0x00c4, 0x004c, 0x1001, 0x0801, 0x0401, 0x0201, 0x00a8, 0x008a,
            0x0201, 0x00d0, 0x000d, 0x0201, 0x00d1, 0x0201, 0x004b, 0x0201, 0x0097, 0x00a7, 0x0c01,
            0x0601, 0x0201, 0x00c3, 0x0201, 0x007a, 0x0099, 0x0401, 0x0201, 0x00c5, 0x005c, 0x00b7,
            0x0401, 0x0201, 0x001d, 0x00d2, 0x0201, 0x002d, 0x0201, 0x007b, 0x00d3, 0x3401, 0x1c01,
            0x0c01, 0x0401, 0x0201, 0x003d, 0x00c6, 0x0401, 0x0201, 0x006c, 0x00a9, 0x0201, 0x009a,
            0x00d4, 0x0801, 0x0401, 0x0201, 0x00b8, 0x008b, 0x0201, 0x004d, 0x00c7, 0x0401, 0x0201,
            0x007c, 0x00d5, 0x0201, 0x005d, 0x00e0, 0x0a01, 0x0401, 0x0201, 0x00e1, 0x001e, 0x0401,
            0x0201, 0x000e, 0x002e, 0x00e2, 0x0801, 0x0401, 0x0201, 0x00e3, 0x006d, 0x0201, 0x008c,
            0x00e4, 0x0401, 0x0201, 0x00e5, 0x00ba, 0x00f0, 0x2601, 0x1001, 0x0401, 0x0201, 0x00f1,
            0x001f, 0x0601, 0x0401, 0x0201, 0x00aa, 0x009b, 0x00b9, 0x0201, 0x003e, 0x0201, 0x00d6,
            0x00c8, 0x0c01, 0x0601, 0x0201, 0x004e, 0x0201, 0x00d7, 0x007d, 0x0201, 0x00ab, 0x0201,
            0x005e, 0x00c9, 0x0601, 0x0201, 0x000f, 0x0201, 0x009c, 0x006e, 0x0201, 0x00f2, 0x002f,
            0x2001, 0x1001, 0x0601, 0x0401, 0x0201, 0x00d8, 0x008d, 0x003f, 0x0601, 0x0201, 0x00f3,
            0x0201, 0x00e6, 0x00ca, 0x0201, 0x00f4, 0x004f, 0x0801, 0x0401, 0x0201, 0x00bb, 0x00ac,
            0x0201, 0x00e7, 0x00f5, 0x0401, 0x0201, 0x00d9, 0x009d, 0x0201, 0x005f, 0x00e8, 0x1e01,
            0x0c01, 0x0601, 0x0201, 0x006f, 0x0201, 0x00f6, 0x00cb, 0x0401, 0x0201, 0x00bc, 0x00ad,
            0x00da, 0x0801, 0x0201, 0x00f7, 0x0401, 0x0201, 0x007e, 0x007f, 0x008e, 0x0601, 0x0401,
            0x0201, 0x009e, 0x00ae, 0x00cc, 0x0201, 0x00f8, 0x008f, 0x1201, 0x0801, 0x0401, 0x0201,
            0x00db, 0x00bd, 0x0201, 0x00ea, 0x00f9, 0x0401, 0x0201, 0x009f, 0x00eb, 0x0201, 0x00be,
            0x0201, 0x00cd, 0x00fa, 0x0e01, 0x0401, 0x0201, 0x00dd, 0x00ec, 0x0601, 0x0401, 0x0201,
            0x00e9, 0x00af, 0x00dc, 0x0201, 0x00ce, 0x00fb, 0x0801, 0x0401, 0x0201, 0x00bf, 0x00de,
            0x0201, 0x00cf, 0x00ee, 0x0401, 0x0201, 0x00df, 0x00ef, 0x0201, 0x00ff, 0x0201, 0x00ed,
            0x0201, 0x00fd, 0x0201, 0x00fc, 0x00fe,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 14
    HuffmanTable {
        data: &[],
        linbits: 0,
        quads: false,
    },
    // Table 15
    HuffmanTable {
        data: &[
            0x1001, 0x0601, 0x0201, 0x0000, 0x0201, 0x0010, 0x0001, 0x0201, 0x0011, 0x0401, 0x0201,
            0x0020, 0x0002, 0x0201, 0x0021, 0x0012, 0x3201, 0x1001, 0x0601, 0x0201, 0x0022, 0x0201,
            0x0030, 0x0031, 0x0601, 0x0201, 0x0013, 0x0201, 0x0003, 0x0040, 0x0201, 0x0032, 0x0023,
            0x0e01, 0x0601, 0x0401, 0x0201, 0x0004, 0x0014, 0x0041, 0x0401, 0x0201, 0x0033, 0x0042,
            0x0201, 0x0024, 0x0043, 0x0a01, 0x0601, 0x0201, 0x0034, 0x0201, 0x0050, 0x0005, 0x0201,
            0x0051, 0x0015, 0x0401, 0x0201, 0x0052, 0x0025, 0x0401, 0x0201, 0x0044, 0x0053, 0x0061,
            0x5a01, 0x2401, 0x1201, 0x0a01, 0x0601, 0x0201, 0x0035, 0x0201, 0x0060, 0x0006, 0x0201,
            0x0016, 0x0062, 0x0401, 0x0201, 0x0026, 0x0054, 0x0201, 0x0045, 0x0063, 0x0a01, 0x0601,
            0x0201, 0x0036, 0x0201, 0x0070, 0x0007, 0x0201, 0x0071, 0x0055, 0x0401, 0x0201, 0x0017,
            0x0064, 0x0201, 0x0072, 0x0027, 0x1801, 0x1001, 0x0801, 0x0401, 0x0201, 0x0046, 0x0073,
            0x0201, 0x0037, 0x0065, 0x0401, 0x0201, 0x0056, 0x0080, 0x0201, 0x0008, 0x0074, 0x0401,
            0x0201, 0x0081, 0x0018, 0x0201, 0x0082, 0x0028, 0x1001, 0x0801, 0x0401, 0x0201, 0x0047,
            0x0066, 0x0201, 0x0083, 0x0038, 0x0401, 0x0201, 0x0075, 0x0057, 0x0201, 0x0084, 0x0048,
            0x0601, 0x0401, 0x0201, 0x0090, 0x0019, 0x0091, 0x0401, 0x0201, 0x0092, 0x0076, 0x0201,
            0x0067, 0x0029, 0x5c01, 0x2401, 0x1201, 0x0a01, 0x0401, 0x0201, 0x0085, 0x0058, 0x0401,
            0x0201, 0x0009, 0x0077, 0x0093, 0x0401, 0x0201, 0x0039, 0x0094, 0x0201, 0x0049, 0x0086,
            0x0a01, 0x0601, 0x0201, 0x0068, 0x0201, 0x00a0, 0x000a, 0x0201, 0x00a1, 0x001a, 0x0401,
            0x0201, 0x00a2, 0x002a, 0x0201, 0x0095, 0x0059, 0x1a01, 0x0e01, 0x0601, 0x0201, 0x00a3,
            0x0201, 0x003a, 0x0087, 0x0401, 0x0201, 0x0078, 0x00a4, 0x0201, 0x004a, 0x0096, 0x0601,
            0x0401, 0x0201, 0x0069, 0x00b0, 0x00b1, 0x0401, 0x0201, 0x001b, 0x00a5, 0x00b2, 0x0e01,
            0x0801, 0x0401, 0x0201, 0x005a, 0x002b, 0x0201, 0x0088, 0x0097, 0x0201, 0x00b3, 0x0201,
            0x0079, 0x003b, 0x0801, 0x0401, 0x0201, 0x006a, 0x00b4, 0x0201, 0x004b, 0x00c1, 0x0401,
            0x0201, 0x0098, 0x0089, 0x0201, 0x001c, 0x00b5, 0x5001, 0x2201, 0x1001, 0x0601, 0x0401,
            0x0201, 0x005b, 0x002c, 0x00c2, 0x0601, 0x0401, 0x0201, 0x000b, 0x00c0, 0x00a6, 0x0201,
            0x00a7, 0x007a, 0x0a01, 0x0401, 0x0201, 0x00c3, 0x003c, 0x0401, 0x0201, 0x000c, 0x0099,
            0x00b6, 0x0401, 0x0201, 0x006b, 0x00c4, 0x0201, 0x004c, 0x00a8, 0x1401, 0x0a01, 0x0401,
            0x0201, 0x008a, 0x00c5, 0x0401, 0x0201, 0x00d0, 0x005c, 0x00d1, 0x0401, 0x0201, 0x00b7,
            0x007b, 0x0201, 0x001d, 0x0201, 0x000d, 0x002d, 0x0c01, 0x0401, 0x0201, 0x00d2, 0x00d3,
            0x0401, 0x0201, 0x003d, 0x00c6, 0x0201, 0x006c, 0x00a9, 0x0601, 0x0401, 0x0201, 0x009a,
            0x00b8, 0x00d4, 0x0401, 0x0201, 0x008b, 0x004d, 0x0201, 0x00c7, 0x007c, 0x4401, 0x2201,
            0x1201, 0x0a01, 0x0401, 0x0201, 0x00d5, 0x005d, 0x0401, 0x0201, 0x00e0, 0x000e, 0x00e1,
            0x0401, 0x0201, 0x001e, 0x00e2, 0x0201, 0x00aa, 0x002e, 0x0801, 0x0401, 0x0201, 0x00b9,
            0x009b, 0x0201, 0x00e3, 0x00d6, 0x0401, 0x0201, 0x006d, 0x003e, 0x0201, 0x00c8, 0x008c,
            0x1001, 0x0801, 0x0401, 0x0201, 0x00e4, 0x004e, 0x0201, 0x00d7, 0x007d, 0x0401, 0x0201,
            0x00e5, 0x00ba, 0x0201, 0x00ab, 0x005e, 0x0801, 0x0401, 0x0201, 0x00c9, 0x009c, 0x0201,
            0x00f1, 0x001f, 0x0601, 0x0401, 0x0201, 0x00f0, 0x006e, 0x00f2, 0x0201, 0x002f, 0x00e6,
            0x2601, 0x1201, 0x0801, 0x0401, 0x0201, 0x00d8, 0x00f3, 0x0201, 0x003f, 0x00f4, 0x0601,
            0x0201, 0x004f, 0x0201, 0x008d, 0x00d9, 0x0201, 0x00bb, 0x00ca, 0x0801, 0x0401, 0x0201,
            0x00ac, 0x00e7, 0x0201, 0x007e, 0x00f5, 0x0801, 0x0401, 0x0201, 0x009d, 0x005f, 0x0201,
            0x00e8, 0x008e, 0x0201, 0x00f6, 0x00cb, 0x2201, 0x1201, 0x0a01, 0x0601, 0x0401, 0x0201,
            0x000f, 0x00ae, 0x006f, 0x0201, 0x00bc, 0x00da, 0x0401, 0x0201, 0x00ad, 0x00f7, 0x0201,
            0x007f, 0x00e9, 0x0801, 0x0401, 0x0201, 0x009e, 0x00cc, 0x0201, 0x00f8, 0x008f, 0x0401,
            0x0201, 0x00db, 0x00bd, 0x0201, 0x00ea, 0x00f9, 0x1001, 0x0801, 0x0401, 0x0201, 0x009f,
            0x00dc, 0x0201, 0x00cd, 0x00eb, 0x0401, 0x0201, 0x00be, 0x00fa, 0x0201, 0x00af, 0x00dd,
            0x0e01, 0x0601, 0x0401, 0x0201, 0x00ec, 0x00ce, 0x00fb, 0x0401, 0x0201, 0x00bf, 0x00ed,
            0x0201, 0x00de, 0x00fc, 0x0601, 0x0401, 0x0201, 0x00cf, 0x00fd, 0x00ee, 0x0401, 0x0201,
            0x00df, 0x00fe, 0x0201, 0x00ef, 0x00ff,
        ],
        linbits: 0,
        quads: false,
    },
    // Table 16
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 1,
        quads: false,
    },
    // Table 17
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 2,
        quads: false,
    },
    // Table 18
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 3,
        quads: false,
    },
    // Table 19
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 4,
        quads: false,
    },
    // Table 20
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 6,
        quads: false,
    },
    // Table 21
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 8,
        quads: false,
    },
    // Table 22
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 10,
        quads: false,
    },
    // Table 23
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_1,
        linbits: 13,
        quads: false,
    },
    // Table 24
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 4,
        quads: false,
    },
    // Table 25
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 5,
        quads: false,
    },
    // Table 26
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 6,
        quads: false,
    },
    // Table 27
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 7,
        quads: false,
    },
    // Table 28
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 8,
        quads: false,
    },
    // Table 29
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 9,
        quads: false,
    },
    // Table 30
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 11,
        quads: false,
    },
    // Table 31
    HuffmanTable {
        data: &HUFFMAN_TABLES_LIN_2,
        linbits: 13,
        quads: false,
    },
    // Table 32
    HuffmanTable {
        data: &[
            0x0201, 0x0000, 0x0801, 0x0401, 0x0201, 0x0008, 0x0004, 0x0201, 0x0001, 0x0002, 0x0801,
            0x0401, 0x0201, 0x000c, 0x000a, 0x0201, 0x0003, 0x0006, 0x0601, 0x0201, 0x0009, 0x0201,
            0x0005, 0x0007, 0x0401, 0x0201, 0x000e, 0x000d, 0x0201, 0x000f, 0x000b,
        ],
        linbits: 0,
        quads: true,
    },
    // Table 33
    HuffmanTable {
        // TODO (Herschel): Check this one
        data: &[
            0x1001, 0x0801, 0x0401, 0x0201, 0x0000, 0x0001, 0x0201, 0x0002, 0x0003, 0x0401, 0x0201,
            0x0004, 0x0005, 0x0201, 0x0006, 0x0007, 0x0801, 0x0401, 0x0201, 0x0008, 0x0009, 0x0201,
            0x000a, 0x000b, 0x0401, 0x0201, 0x000c, 0x000d, 0x0201, 0x000e, 0x000f,
        ],
        linbits: 0,
        quads: true,
    },
];

const HUFFMAN_TABLES_LIN_1: [u16; 511] = [
    0x0201, 0x0000, 0x0601, 0x0201, 0x0010, 0x0201, 0x0001, 0x0011, 0x2a01, 0x0801, 0x0401, 0x0201,
    0x0020, 0x0002, 0x0201, 0x0021, 0x0012, 0x0a01, 0x0601, 0x0201, 0x0022, 0x0201, 0x0030, 0x0003,
    0x0201, 0x0031, 0x0013, 0x0a01, 0x0401, 0x0201, 0x0032, 0x0023, 0x0401, 0x0201, 0x0040, 0x0004,
    0x0041, 0x0601, 0x0201, 0x0014, 0x0201, 0x0033, 0x0042, 0x0401, 0x0201, 0x0024, 0x0050, 0x0201,
    0x0043, 0x0034, 0x8a01, 0x2801, 0x1001, 0x0601, 0x0401, 0x0201, 0x0005, 0x0015, 0x0051, 0x0401,
    0x0201, 0x0052, 0x0025, 0x0401, 0x0201, 0x0044, 0x0035, 0x0053, 0x0a01, 0x0601, 0x0401, 0x0201,
    0x0060, 0x0006, 0x0061, 0x0201, 0x0016, 0x0062, 0x0801, 0x0401, 0x0201, 0x0026, 0x0054, 0x0201,
    0x0045, 0x0063, 0x0401, 0x0201, 0x0036, 0x0070, 0x0071, 0x2801, 0x1201, 0x0801, 0x0201, 0x0017,
    0x0201, 0x0007, 0x0201, 0x0055, 0x0064, 0x0401, 0x0201, 0x0072, 0x0027, 0x0401, 0x0201, 0x0046,
    0x0065, 0x0073, 0x0a01, 0x0601, 0x0201, 0x0037, 0x0201, 0x0056, 0x0008, 0x0201, 0x0080, 0x0081,
    0x0601, 0x0201, 0x0018, 0x0201, 0x0074, 0x0047, 0x0201, 0x0082, 0x0201, 0x0028, 0x0066, 0x1801,
    0x0e01, 0x0801, 0x0401, 0x0201, 0x0083, 0x0038, 0x0201, 0x0075, 0x0084, 0x0401, 0x0201, 0x0048,
    0x0090, 0x0091, 0x0601, 0x0201, 0x0019, 0x0201, 0x0009, 0x0076, 0x0201, 0x0092, 0x0029, 0x0e01,
    0x0801, 0x0401, 0x0201, 0x0085, 0x0058, 0x0201, 0x0093, 0x0039, 0x0401, 0x0201, 0x00a0, 0x000a,
    0x001a, 0x0801, 0x0201, 0x00a2, 0x0201, 0x0067, 0x0201, 0x0057, 0x0049, 0x0601, 0x0201, 0x0094,
    0x0201, 0x0077, 0x0086, 0x0201, 0x00a1, 0x0201, 0x0068, 0x0095, 0xdc01, 0x7e01, 0x3201, 0x1a01,
    0x0c01, 0x0601, 0x0201, 0x002a, 0x0201, 0x0059, 0x003a, 0x0201, 0x00a3, 0x0201, 0x0087, 0x0078,
    0x0801, 0x0401, 0x0201, 0x00a4, 0x004a, 0x0201, 0x0096, 0x0069, 0x0401, 0x0201, 0x00b0, 0x000b,
    0x00b1, 0x0a01, 0x0401, 0x0201, 0x001b, 0x00b2, 0x0201, 0x002b, 0x0201, 0x00a5, 0x005a, 0x0601,
    0x0201, 0x00b3, 0x0201, 0x00a6, 0x006a, 0x0401, 0x0201, 0x00b4, 0x004b, 0x0201, 0x000c, 0x00c1,
    0x1e01, 0x0e01, 0x0601, 0x0401, 0x0201, 0x00b5, 0x00c2, 0x002c, 0x0401, 0x0201, 0x00a7, 0x00c3,
    0x0201, 0x006b, 0x00c4, 0x0801, 0x0201, 0x001d, 0x0401, 0x0201, 0x0088, 0x0097, 0x003b, 0x0401,
    0x0201, 0x00d1, 0x00d2, 0x0201, 0x002d, 0x00d3, 0x1201, 0x0601, 0x0401, 0x0201, 0x001e, 0x002e,
    0x00e2, 0x0601, 0x0401, 0x0201, 0x0079, 0x0098, 0x00c0, 0x0201, 0x001c, 0x0201, 0x0089, 0x005b,
    0x0e01, 0x0601, 0x0201, 0x003c, 0x0201, 0x007a, 0x00b6, 0x0401, 0x0201, 0x004c, 0x0099, 0x0201,
    0x00a8, 0x008a, 0x0601, 0x0201, 0x000d, 0x0201, 0x00c5, 0x005c, 0x0401, 0x0201, 0x003d, 0x00c6,
    0x0201, 0x006c, 0x009a, 0x5801, 0x5601, 0x2401, 0x1001, 0x0801, 0x0401, 0x0201, 0x008b, 0x004d,
    0x0201, 0x00c7, 0x007c, 0x0401, 0x0201, 0x00d5, 0x005d, 0x0201, 0x00e0, 0x000e, 0x0801, 0x0201,
    0x00e3, 0x0401, 0x0201, 0x00d0, 0x00b7, 0x007b, 0x0601, 0x0401, 0x0201, 0x00a9, 0x00b8, 0x00d4,
    0x0201, 0x00e1, 0x0201, 0x00aa, 0x00b9, 0x1801, 0x0a01, 0x0601, 0x0401, 0x0201, 0x009b, 0x00d6,
    0x006d, 0x0201, 0x003e, 0x00c8, 0x0601, 0x0401, 0x0201, 0x008c, 0x00e4, 0x004e, 0x0401, 0x0201,
    0x00d7, 0x00e5, 0x0201, 0x00ba, 0x00ab, 0x0c01, 0x0401, 0x0201, 0x009c, 0x00e6, 0x0401, 0x0201,
    0x006e, 0x00d8, 0x0201, 0x008d, 0x00bb, 0x0801, 0x0401, 0x0201, 0x00e7, 0x009d, 0x0201, 0x00e8,
    0x008e, 0x0401, 0x0201, 0x00cb, 0x00bc, 0x009e, 0x00f1, 0x0201, 0x001f, 0x0201, 0x000f, 0x002f,
    0x4201, 0x3801, 0x0201, 0x00f2, 0x3401, 0x3201, 0x1401, 0x0801, 0x0201, 0x00bd, 0x0201, 0x005e,
    0x0201, 0x007d, 0x00c9, 0x0601, 0x0201, 0x00ca, 0x0201, 0x00ac, 0x007e, 0x0401, 0x0201, 0x00da,
    0x00ad, 0x00cc, 0x0a01, 0x0601, 0x0201, 0x00ae, 0x0201, 0x00db, 0x00dc, 0x0201, 0x00cd, 0x00be,
    0x0601, 0x0401, 0x0201, 0x00eb, 0x00ed, 0x00ee, 0x0601, 0x0401, 0x0201, 0x00d9, 0x00ea, 0x00e9,
    0x0201, 0x00de, 0x0401, 0x0201, 0x00dd, 0x00ec, 0x00ce, 0x003f, 0x00f0, 0x0401, 0x0201, 0x00f3,
    0x00f4, 0x0201, 0x004f, 0x0201, 0x00f5, 0x005f, 0x0a01, 0x0201, 0x00ff, 0x0401, 0x0201, 0x00f6,
    0x006f, 0x0201, 0x00f7, 0x007f, 0x0c01, 0x0601, 0x0201, 0x008f, 0x0201, 0x00f8, 0x00f9, 0x0401,
    0x0201, 0x009f, 0x00fa, 0x00af, 0x0801, 0x0401, 0x0201, 0x00fb, 0x00bf, 0x0201, 0x00fc, 0x00cf,
    0x0401, 0x0201, 0x00fd, 0x00df, 0x0201, 0x00fe, 0x00ef,
];

const HUFFMAN_TABLES_LIN_2: [u16; 512] = [
    0x3c01, 0x0801, 0x0401, 0x0201, 0x0000, 0x0010, 0x0201, 0x0001, 0x0011, 0x0e01, 0x0601, 0x0401,
    0x0201, 0x0020, 0x0002, 0x0021, 0x0201, 0x0012, 0x0201, 0x0022, 0x0201, 0x0030, 0x0003, 0x0e01,
    0x0401, 0x0201, 0x0031, 0x0013, 0x0401, 0x0201, 0x0032, 0x0023, 0x0401, 0x0201, 0x0040, 0x0004,
    0x0041, 0x0801, 0x0401, 0x0201, 0x0014, 0x0033, 0x0201, 0x0042, 0x0024, 0x0601, 0x0401, 0x0201,
    0x0043, 0x0034, 0x0051, 0x0601, 0x0401, 0x0201, 0x0050, 0x0005, 0x0015, 0x0201, 0x0052, 0x0025,
    0xfa01, 0x6201, 0x2201, 0x1201, 0x0a01, 0x0401, 0x0201, 0x0044, 0x0053, 0x0201, 0x0035, 0x0201,
    0x0060, 0x0006, 0x0401, 0x0201, 0x0061, 0x0016, 0x0201, 0x0062, 0x0026, 0x0801, 0x0401, 0x0201,
    0x0054, 0x0045, 0x0201, 0x0063, 0x0036, 0x0401, 0x0201, 0x0071, 0x0055, 0x0201, 0x0064, 0x0046,
    0x2001, 0x0e01, 0x0601, 0x0201, 0x0072, 0x0201, 0x0027, 0x0037, 0x0201, 0x0073, 0x0401, 0x0201,
    0x0070, 0x0007, 0x0017, 0x0a01, 0x0401, 0x0201, 0x0065, 0x0056, 0x0401, 0x0201, 0x0080, 0x0008,
    0x0081, 0x0401, 0x0201, 0x0074, 0x0047, 0x0201, 0x0018, 0x0082, 0x1001, 0x0801, 0x0401, 0x0201,
    0x0028, 0x0066, 0x0201, 0x0083, 0x0038, 0x0401, 0x0201, 0x0075, 0x0057, 0x0201, 0x0084, 0x0048,
    0x0801, 0x0401, 0x0201, 0x0091, 0x0019, 0x0201, 0x0092, 0x0076, 0x0401, 0x0201, 0x0067, 0x0029,
    0x0201, 0x0085, 0x0058, 0x5c01, 0x2201, 0x1001, 0x0801, 0x0401, 0x0201, 0x0093, 0x0039, 0x0201,
    0x0094, 0x0049, 0x0401, 0x0201, 0x0077, 0x0086, 0x0201, 0x0068, 0x00a1, 0x0801, 0x0401, 0x0201,
    0x00a2, 0x002a, 0x0201, 0x0095, 0x0059, 0x0401, 0x0201, 0x00a3, 0x003a, 0x0201, 0x0087, 0x0201,
    0x0078, 0x004a, 0x1601, 0x0c01, 0x0401, 0x0201, 0x00a4, 0x0096, 0x0401, 0x0201, 0x0069, 0x00b1,
    0x0201, 0x001b, 0x00a5, 0x0601, 0x0201, 0x00b2, 0x0201, 0x005a, 0x002b, 0x0201, 0x0088, 0x00b3,
    0x1001, 0x0a01, 0x0601, 0x0201, 0x0090, 0x0201, 0x0009, 0x00a0, 0x0201, 0x0097, 0x0079, 0x0401,
    0x0201, 0x00a6, 0x006a, 0x00b4, 0x0c01, 0x0601, 0x0201, 0x001a, 0x0201, 0x000a, 0x00b0, 0x0201,
    0x003b, 0x0201, 0x000b, 0x00c0, 0x0401, 0x0201, 0x004b, 0x00c1, 0x0201, 0x0098, 0x0089, 0x4301,
    0x2201, 0x1001, 0x0801, 0x0401, 0x0201, 0x001c, 0x00b5, 0x0201, 0x005b, 0x00c2, 0x0401, 0x0201,
    0x002c, 0x00a7, 0x0201, 0x007a, 0x00c3, 0x0a01, 0x0601, 0x0201, 0x003c, 0x0201, 0x000c, 0x00d0,
    0x0201, 0x00b6, 0x006b, 0x0401, 0x0201, 0x00c4, 0x004c, 0x0201, 0x0099, 0x00a8, 0x1001, 0x0801,
    0x0401, 0x0201, 0x008a, 0x00c5, 0x0201, 0x005c, 0x00d1, 0x0401, 0x0201, 0x00b7, 0x007b, 0x0201,
    0x001d, 0x00d2, 0x0901, 0x0401, 0x0201, 0x002d, 0x00d3, 0x0201, 0x003d, 0x00c6, 0x55fa, 0x0401,
    0x0201, 0x006c, 0x00a9, 0x0201, 0x009a, 0x00d4, 0x2001, 0x1001, 0x0801, 0x0401, 0x0201, 0x00b8,
    0x008b, 0x0201, 0x004d, 0x00c7, 0x0401, 0x0201, 0x007c, 0x00d5, 0x0201, 0x005d, 0x00e1, 0x0801,
    0x0401, 0x0201, 0x001e, 0x00e2, 0x0201, 0x00aa, 0x00b9, 0x0401, 0x0201, 0x009b, 0x00e3, 0x0201,
    0x00d6, 0x006d, 0x1401, 0x0a01, 0x0601, 0x0201, 0x003e, 0x0201, 0x002e, 0x004e, 0x0201, 0x00c8,
    0x008c, 0x0401, 0x0201, 0x00e4, 0x00d7, 0x0401, 0x0201, 0x007d, 0x00ab, 0x00e5, 0x0a01, 0x0401,
    0x0201, 0x00ba, 0x005e, 0x0201, 0x00c9, 0x0201, 0x009c, 0x006e, 0x0801, 0x0201, 0x00e6, 0x0201,
    0x000d, 0x0201, 0x00e0, 0x000e, 0x0401, 0x0201, 0x00d8, 0x008d, 0x0201, 0x00bb, 0x00ca, 0x4a01,
    0x0201, 0x00ff, 0x4001, 0x3a01, 0x2001, 0x1001, 0x0801, 0x0401, 0x0201, 0x00ac, 0x00e7, 0x0201,
    0x007e, 0x00d9, 0x0401, 0x0201, 0x009d, 0x00e8, 0x0201, 0x008e, 0x00cb, 0x0801, 0x0401, 0x0201,
    0x00bc, 0x00da, 0x0201, 0x00ad, 0x00e9, 0x0401, 0x0201, 0x009e, 0x00cc, 0x0201, 0x00db, 0x00bd,
    0x1001, 0x0801, 0x0401, 0x0201, 0x00ea, 0x00ae, 0x0201, 0x00dc, 0x00cd, 0x0401, 0x0201, 0x00eb,
    0x00be, 0x0201, 0x00dd, 0x00ec, 0x0801, 0x0401, 0x0201, 0x00ce, 0x00ed, 0x0201, 0x00de, 0x00ee,
    0x000f, 0x0401, 0x0201, 0x00f0, 0x001f, 0x00f1, 0x0401, 0x0201, 0x00f2, 0x002f, 0x0201, 0x00f3,
    0x003f, 0x1201, 0x0801, 0x0401, 0x0201, 0x00f4, 0x004f, 0x0201, 0x00f5, 0x005f, 0x0401, 0x0201,
    0x00f6, 0x006f, 0x0201, 0x00f7, 0x0201, 0x007f, 0x008f, 0x0a01, 0x0401, 0x0201, 0x00f8, 0x00f9,
    0x0401, 0x0201, 0x009f, 0x00af, 0x00fa, 0x0801, 0x0401, 0x0201, 0x00fb, 0x00bf, 0x0201, 0x00fc,
    0x00cf, 0x0401, 0x0201, 0x00fd, 0x00df, 0x0201, 0x00fe, 0x00ef,
];

#[allow(clippy::unreadable_literal)]
#[allow(clippy::excessive_precision)]
const SYNTH_DTBL: [f32; 512] = [
    0.000000000,
    -0.000015259,
    -0.000015259,
    -0.000015259,
    -0.000015259,
    -0.000015259,
    -0.000015259,
    -0.000030518,
    -0.000030518,
    -0.000030518,
    -0.000030518,
    -0.000045776,
    -0.000045776,
    -0.000061035,
    -0.000061035,
    -0.000076294,
    -0.000076294,
    -0.000091553,
    -0.000106812,
    -0.000106812,
    -0.000122070,
    -0.000137329,
    -0.000152588,
    -0.000167847,
    -0.000198364,
    -0.000213623,
    -0.000244141,
    -0.000259399,
    -0.000289917,
    -0.000320435,
    -0.000366211,
    -0.000396729,
    -0.000442505,
    -0.000473022,
    -0.000534058,
    -0.000579834,
    -0.000625610,
    -0.000686646,
    -0.000747681,
    -0.000808716,
    -0.000885010,
    -0.000961304,
    -0.001037598,
    -0.001113892,
    -0.001205444,
    -0.001296997,
    -0.001388550,
    -0.001480103,
    -0.001586914,
    -0.001693726,
    -0.001785278,
    -0.001907349,
    -0.002014160,
    -0.002120972,
    -0.002243042,
    -0.002349854,
    -0.002456665,
    -0.002578735,
    -0.002685547,
    -0.002792358,
    -0.002899170,
    -0.002990723,
    -0.003082275,
    -0.003173828,
    0.003250122,
    0.003326416,
    0.003387451,
    0.003433228,
    0.003463745,
    0.003479004,
    0.003479004,
    0.003463745,
    0.003417969,
    0.003372192,
    0.003280640,
    0.003173828,
    0.003051758,
    0.002883911,
    0.002700806,
    0.002487183,
    0.002227783,
    0.001937866,
    0.001617432,
    0.001266479,
    0.000869751,
    0.000442505,
    -0.000030518,
    -0.000549316,
    -0.001098633,
    -0.001693726,
    -0.002334595,
    -0.003005981,
    -0.003723145,
    -0.004486084,
    -0.005294800,
    -0.006118774,
    -0.007003784,
    -0.007919312,
    -0.008865356,
    -0.009841919,
    -0.010848999,
    -0.011886597,
    -0.012939453,
    -0.014022827,
    -0.015121460,
    -0.016235352,
    -0.017349243,
    -0.018463135,
    -0.019577026,
    -0.020690918,
    -0.021789551,
    -0.022857666,
    -0.023910522,
    -0.024932861,
    -0.025909424,
    -0.026840210,
    -0.027725220,
    -0.028533936,
    -0.029281616,
    -0.029937744,
    -0.030532837,
    -0.031005859,
    -0.031387329,
    -0.031661987,
    -0.031814575,
    -0.031845093,
    -0.031738281,
    -0.031478882,
    0.031082153,
    0.030517578,
    0.029785156,
    0.028884888,
    0.027801514,
    0.026535034,
    0.025085449,
    0.023422241,
    0.021575928,
    0.019531250,
    0.017257690,
    0.014801025,
    0.012115479,
    0.009231567,
    0.006134033,
    0.002822876,
    -0.000686646,
    -0.004394531,
    -0.008316040,
    -0.012420654,
    -0.016708374,
    -0.021179199,
    -0.025817871,
    -0.030609131,
    -0.035552979,
    -0.040634155,
    -0.045837402,
    -0.051132202,
    -0.056533813,
    -0.061996460,
    -0.067520142,
    -0.073059082,
    -0.078628540,
    -0.084182739,
    -0.089706421,
    -0.095169067,
    -0.100540161,
    -0.105819702,
    -0.110946655,
    -0.115921021,
    -0.120697021,
    -0.125259399,
    -0.129562378,
    -0.133590698,
    -0.137298584,
    -0.140670776,
    -0.143676758,
    -0.146255493,
    -0.148422241,
    -0.150115967,
    -0.151306152,
    -0.151962280,
    -0.152069092,
    -0.151596069,
    -0.150497437,
    -0.148773193,
    -0.146362305,
    -0.143264771,
    -0.139450073,
    -0.134887695,
    -0.129577637,
    -0.123474121,
    -0.116577148,
    -0.108856201,
    0.100311279,
    0.090927124,
    0.080688477,
    0.069595337,
    0.057617188,
    0.044784546,
    0.031082153,
    0.016510010,
    0.001068115,
    -0.015228271,
    -0.032379150,
    -0.050354004,
    -0.069168091,
    -0.088775635,
    -0.109161377,
    -0.130310059,
    -0.152206421,
    -0.174789429,
    -0.198059082,
    -0.221984863,
    -0.246505737,
    -0.271591187,
    -0.297210693,
    -0.323318481,
    -0.349868774,
    -0.376800537,
    -0.404083252,
    -0.431655884,
    -0.459472656,
    -0.487472534,
    -0.515609741,
    -0.543823242,
    -0.572036743,
    -0.600219727,
    -0.628295898,
    -0.656219482,
    -0.683914185,
    -0.711318970,
    -0.738372803,
    -0.765029907,
    -0.791213989,
    -0.816864014,
    -0.841949463,
    -0.866363525,
    -0.890090942,
    -0.913055420,
    -0.935195923,
    -0.956481934,
    -0.976852417,
    -0.996246338,
    -1.014617920,
    -1.031936646,
    -1.048156738,
    -1.063217163,
    -1.077117920,
    -1.089782715,
    -1.101211548,
    -1.111373901,
    -1.120223999,
    -1.127746582,
    -1.133926392,
    -1.138763428,
    -1.142211914,
    -1.144287109,
    1.144989014,
    1.144287109,
    1.142211914,
    1.138763428,
    1.133926392,
    1.127746582,
    1.120223999,
    1.111373901,
    1.101211548,
    1.089782715,
    1.077117920,
    1.063217163,
    1.048156738,
    1.031936646,
    1.014617920,
    0.996246338,
    0.976852417,
    0.956481934,
    0.935195923,
    0.913055420,
    0.890090942,
    0.866363525,
    0.841949463,
    0.816864014,
    0.791213989,
    0.765029907,
    0.738372803,
    0.711318970,
    0.683914185,
    0.656219482,
    0.628295898,
    0.600219727,
    0.572036743,
    0.543823242,
    0.515609741,
    0.487472534,
    0.459472656,
    0.431655884,
    0.404083252,
    0.376800537,
    0.349868774,
    0.323318481,
    0.297210693,
    0.271591187,
    0.246505737,
    0.221984863,
    0.198059082,
    0.174789429,
    0.152206421,
    0.130310059,
    0.109161377,
    0.088775635,
    0.069168091,
    0.050354004,
    0.032379150,
    0.015228271,
    -0.001068115,
    -0.016510010,
    -0.031082153,
    -0.044784546,
    -0.057617188,
    -0.069595337,
    -0.080688477,
    -0.090927124,
    0.100311279,
    0.108856201,
    0.116577148,
    0.123474121,
    0.129577637,
    0.134887695,
    0.139450073,
    0.143264771,
    0.146362305,
    0.148773193,
    0.150497437,
    0.151596069,
    0.152069092,
    0.151962280,
    0.151306152,
    0.150115967,
    0.148422241,
    0.146255493,
    0.143676758,
    0.140670776,
    0.137298584,
    0.133590698,
    0.129562378,
    0.125259399,
    0.120697021,
    0.115921021,
    0.110946655,
    0.105819702,
    0.100540161,
    0.095169067,
    0.089706421,
    0.084182739,
    0.078628540,
    0.073059082,
    0.067520142,
    0.061996460,
    0.056533813,
    0.051132202,
    0.045837402,
    0.040634155,
    0.035552979,
    0.030609131,
    0.025817871,
    0.021179199,
    0.016708374,
    0.012420654,
    0.008316040,
    0.004394531,
    0.000686646,
    -0.002822876,
    -0.006134033,
    -0.009231567,
    -0.012115479,
    -0.014801025,
    -0.017257690,
    -0.019531250,
    -0.021575928,
    -0.023422241,
    -0.025085449,
    -0.026535034,
    -0.027801514,
    -0.028884888,
    -0.029785156,
    -0.030517578,
    0.031082153,
    0.031478882,
    0.031738281,
    0.031845093,
    0.031814575,
    0.031661987,
    0.031387329,
    0.031005859,
    0.030532837,
    0.029937744,
    0.029281616,
    0.028533936,
    0.027725220,
    0.026840210,
    0.025909424,
    0.024932861,
    0.023910522,
    0.022857666,
    0.021789551,
    0.020690918,
    0.019577026,
    0.018463135,
    0.017349243,
    0.016235352,
    0.015121460,
    0.014022827,
    0.012939453,
    0.011886597,
    0.010848999,
    0.009841919,
    0.008865356,
    0.007919312,
    0.007003784,
    0.006118774,
    0.005294800,
    0.004486084,
    0.003723145,
    0.003005981,
    0.002334595,
    0.001693726,
    0.001098633,
    0.000549316,
    0.000030518,
    -0.000442505,
    -0.000869751,
    -0.001266479,
    -0.001617432,
    -0.001937866,
    -0.002227783,
    -0.002487183,
    -0.002700806,
    -0.002883911,
    -0.003051758,
    -0.003173828,
    -0.003280640,
    -0.003372192,
    -0.003417969,
    -0.003463745,
    -0.003479004,
    -0.003479004,
    -0.003463745,
    -0.003433228,
    -0.003387451,
    -0.003326416,
    0.003250122,
    0.003173828,
    0.003082275,
    0.002990723,
    0.002899170,
    0.002792358,
    0.002685547,
    0.002578735,
    0.002456665,
    0.002349854,
    0.002243042,
    0.002120972,
    0.002014160,
    0.001907349,
    0.001785278,
    0.001693726,
    0.001586914,
    0.001480103,
    0.001388550,
    0.001296997,
    0.001205444,
    0.001113892,
    0.001037598,
    0.000961304,
    0.000885010,
    0.000808716,
    0.000747681,
    0.000686646,
    0.000625610,
    0.000579834,
    0.000534058,
    0.000473022,
    0.000442505,
    0.000396729,
    0.000366211,
    0.000320435,
    0.000289917,
    0.000259399,
    0.000244141,
    0.000213623,
    0.000198364,
    0.000167847,
    0.000152588,
    0.000137329,
    0.000122070,
    0.000106812,
    0.000106812,
    0.000091553,
    0.000076294,
    0.000076294,
    0.000061035,
    0.000061035,
    0.000045776,
    0.000045776,
    0.000030518,
    0.000030518,
    0.000030518,
    0.000030518,
    0.000015259,
    0.000015259,
    0.000015259,
    0.000015259,
    0.000015259,
    0.000015259,
];

#[allow(clippy::unreadable_literal)]
#[allow(clippy::excessive_precision)]
#[allow(clippy::approx_constant)]
const SBS_N_WIN: [[f32; 32]; 64] = [
    [
        0.7071067811865476,
        -0.7071067811865475,
        -0.7071067811865477,
        0.7071067811865474,
        0.7071067811865477,
        -0.7071067811865467,
        -0.7071067811865471,
        0.7071067811865466,
        0.7071067811865472,
        -0.7071067811865465,
        -0.7071067811865474,
        0.7071067811865464,
        0.7071067811865475,
        -0.7071067811865464,
        -0.7071067811865476,
        0.7071067811865462,
        0.7071067811865476,
        -0.7071067811865461,
        -0.7071067811865477,
        0.707106781186546,
        0.7071067811865503,
        -0.707106781186546,
        -0.7071067811865479,
        0.7071067811865483,
        0.7071067811865505,
        -0.7071067811865458,
        -0.707106781186548,
        0.7071067811865482,
        0.7071067811865507,
        -0.7071067811865456,
        -0.7071067811865482,
        0.707106781186548,
    ],
    [
        0.6715589548470183,
        -0.8032075314806448,
        -0.5141027441932218,
        0.9039892931234431,
        0.33688985339222005,
        -0.9700312531945441,
        -0.14673047445536166,
        0.9987954562051724,
        -0.04906767432741729,
        -0.9891765099647811,
        0.24298017990326243,
        0.9415440651830208,
        -0.42755509343028003,
        -0.8577286100002726,
        0.5956993044924337,
        0.7409511253549602,
        -0.7409511253549589,
        -0.5956993044924354,
        0.8577286100002715,
        0.4275550934302851,
        -0.9415440651830214,
        -0.24298017990326443,
        0.9891765099647806,
        0.04906767432742292,
        -0.9987954562051724,
        0.1467304744553596,
        0.9700312531945451,
        -0.3368898533922206,
        -0.903989293123444,
        0.5141027441932186,
        0.8032075314806442,
        -0.6715589548470177,
    ],
    [
        0.6343932841636455,
        -0.8819212643483549,
        -0.29028467725446244,
        0.9951847266721969,
        -0.09801714032955997,
        -0.9569403357322087,
        0.47139673682599736,
        0.7730104533627377,
        -0.773010453362737,
        -0.4713967368259983,
        0.9569403357322089,
        0.09801714032956282,
        -0.995184726672197,
        0.2902846772544622,
        0.8819212643483563,
        -0.6343932841636443,
        -0.6343932841636459,
        0.8819212643483553,
        0.29028467725446433,
        -0.9951847266721968,
        0.09801714032956063,
        0.9569403357322086,
        -0.47139673682599326,
        -0.7730104533627394,
        0.7730104533627351,
        0.4713967368259993,
        -0.9569403357322086,
        -0.09801714032956038,
        0.9951847266721968,
        -0.2902846772544578,
        -0.8819212643483568,
        0.6343932841636434,
    ],
    [
        0.5956993044924335,
        -0.9415440651830207,
        -0.04906767432741803,
        0.9700312531945441,
        -0.5141027441932214,
        -0.6715589548470181,
        0.903989293123443,
        0.1467304744553618,
        -0.9891765099647811,
        0.42755509343028014,
        0.7409511253549601,
        -0.8577286100002717,
        -0.24298017990326395,
        0.9987954562051724,
        -0.33688985339221794,
        -0.8032075314806458,
        0.8032075314806444,
        0.33688985339222016,
        -0.9987954562051723,
        0.24298017990326515,
        0.8577286100002729,
        -0.7409511253549561,
        -0.4275550934302822,
        0.9891765099647806,
        -0.146730474455363,
        -0.903989293123444,
        0.6715589548470151,
        0.5141027441932219,
        -0.9700312531945433,
        0.04906767432741926,
        0.9415440651830214,
        -0.5956993044924298,
    ],
    [
        0.5555702330196023,
        -0.9807852804032304,
        0.1950903220161283,
        0.8314696123025455,
        -0.8314696123025451,
        -0.19509032201612803,
        0.9807852804032307,
        -0.5555702330196015,
        -0.5555702330196026,
        0.9807852804032304,
        -0.19509032201612858,
        -0.8314696123025449,
        0.8314696123025438,
        0.19509032201613036,
        -0.9807852804032308,
        0.5555702330196011,
        0.5555702330196061,
        -0.9807852804032297,
        0.19509032201612447,
        0.8314696123025471,
        -0.8314696123025435,
        -0.19509032201613097,
        0.9807852804032309,
        -0.5555702330196005,
        -0.5555702330196036,
        0.9807852804032302,
        -0.19509032201612736,
        -0.8314696123025456,
        0.8314696123025451,
        0.19509032201612808,
        -0.9807852804032303,
        0.555570233019603,
    ],
    [
        0.5141027441932217,
        -0.9987954562051724,
        0.42755509343028214,
        0.5956993044924332,
        -0.989176509964781,
        0.3368898533922202,
        0.6715589548470182,
        -0.9700312531945441,
        0.24298017990326243,
        0.7409511253549601,
        -0.9415440651830203,
        0.14673047445536033,
        0.8032075314806457,
        -0.9039892931234428,
        0.04906767432741668,
        0.8577286100002728,
        -0.8577286100002696,
        -0.04906767432741925,
        0.9039892931234453,
        -0.8032075314806442,
        -0.1467304744553664,
        0.9415440651830211,
        -0.740951125354956,
        -0.24298017990326493,
        0.9700312531945451,
        -0.6715589548470177,
        -0.3368898533922243,
        0.9891765099647811,
        -0.5956993044924298,
        -0.42755509343028286,
        0.9987954562051726,
        -0.5141027441932149,
    ],
    [
        0.4713967368259978,
        -0.9951847266721969,
        0.6343932841636456,
        0.29028467725446255,
        -0.9569403357322087,
        0.773010453362737,
        0.09801714032956081,
        -0.8819212643483562,
        0.8819212643483555,
        -0.09801714032956124,
        -0.7730104533627368,
        0.9569403357322088,
        -0.2902846772544621,
        -0.6343932841636459,
        0.9951847266721968,
        -0.4713967368259935,
        -0.47139673682599587,
        0.9951847266721967,
        -0.6343932841636466,
        -0.2902846772544613,
        0.9569403357322086,
        -0.7730104533627373,
        -0.09801714032956038,
        0.881921264348355,
        -0.8819212643483548,
        0.0980171403295599,
        0.7730104533627375,
        -0.9569403357322085,
        0.29028467725446083,
        0.634393284163647,
        -0.9951847266721959,
        0.47139673682598915,
    ],
    [
        0.4275550934302822,
        -0.970031253194544,
        0.803207531480645,
        -0.04906767432741754,
        -0.7409511253549599,
        0.989176509964781,
        -0.5141027441932212,
        -0.33688985339221955,
        0.9415440651830208,
        -0.8577286100002717,
        0.14673047445536033,
        0.6715589548470199,
        -0.9987954562051723,
        0.5956993044924335,
        0.24298017990326776,
        -0.9039892931234438,
        0.9039892931234441,
        -0.2429801799032616,
        -0.5956993044924329,
        0.9987954562051726,
        -0.6715589548470179,
        -0.14673047445536663,
        0.8577286100002731,
        -0.941544065183021,
        0.33688985339221694,
        0.5141027441932221,
        -0.9891765099647817,
        0.740951125354958,
        0.04906767432741681,
        -0.8032075314806467,
        0.9700312531945422,
        -0.42755509343028464,
    ],
    [
        0.38268343236508984,
        -0.9238795325112868,
        0.9238795325112865,
        -0.3826834323650899,
        -0.38268343236509056,
        0.9238795325112867,
        -0.9238795325112864,
        0.38268343236508956,
        0.3826834323650909,
        -0.9238795325112876,
        0.9238795325112868,
        -0.3826834323650892,
        -0.3826834323650912,
        0.9238795325112877,
        -0.9238795325112854,
        0.38268343236508556,
        0.3826834323650883,
        -0.9238795325112865,
        0.9238795325112866,
        -0.3826834323650885,
        -0.38268343236509195,
        0.9238795325112881,
        -0.9238795325112851,
        0.3826834323650849,
        0.382683432365089,
        -0.9238795325112868,
        0.9238795325112863,
        -0.3826834323650813,
        -0.3826834323650926,
        0.9238795325112856,
        -0.9238795325112849,
        0.3826834323650908,
    ],
    [
        0.33688985339222005,
        -0.8577286100002721,
        0.9891765099647809,
        -0.6715589548470177,
        0.04906767432741742,
        0.5956993044924335,
        -0.9700312531945443,
        0.9039892931234429,
        -0.42755509343028003,
        -0.24298017990326395,
        0.8032075314806457,
        -0.9987954562051723,
        0.7409511253549588,
        -0.1467304744553635,
        -0.5141027441932244,
        0.941544065183021,
        -0.9415440651830213,
        0.5141027441932188,
        0.146730474455363,
        -0.7409511253549584,
        0.9987954562051726,
        -0.8032075314806439,
        0.24298017990326445,
        0.42755509343028597,
        -0.9039892931234442,
        0.9700312531945441,
        -0.5956993044924352,
        -0.04906767432742757,
        0.6715589548470239,
        -0.9891765099647817,
        0.8577286100002706,
        -0.33688985339221944,
    ],
    [
        0.29028467725446233,
        -0.7730104533627371,
        0.9951847266721969,
        -0.8819212643483548,
        0.47139673682599736,
        0.09801714032956081,
        -0.6343932841636456,
        0.9569403357322094,
        -0.9569403357322089,
        0.6343932841636444,
        -0.09801714032956099,
        -0.47139673682599875,
        0.8819212643483564,
        -0.9951847266721968,
        0.7730104533627352,
        -0.2902846772544582,
        -0.2902846772544613,
        0.7730104533627373,
        -0.9951847266721972,
        0.8819212643483533,
        -0.4713967368259928,
        -0.09801714032956063,
        0.6343932841636468,
        -0.9569403357322098,
        0.9569403357322074,
        -0.6343932841636404,
        0.09801714032955235,
        0.47139673682599387,
        -0.8819212643483538,
        0.9951847266721969,
        -0.7730104533627365,
        0.2902846772544601,
    ],
    [
        0.24298017990326398,
        -0.6715589548470187,
        0.9415440651830209,
        -0.989176509964781,
        0.8032075314806448,
        -0.4275550934302818,
        -0.04906767432741852,
        0.5141027441932239,
        -0.8577286100002726,
        0.9987954562051724,
        -0.9039892931234428,
        0.5956993044924335,
        -0.1467304744553635,
        -0.3368898533922236,
        0.7409511253549605,
        -0.9700312531945442,
        0.9700312531945442,
        -0.740951125354956,
        0.33688985339221716,
        0.14673047445536325,
        -0.5956993044924334,
        0.9039892931234457,
        -0.9987954562051722,
        0.8577286100002709,
        -0.5141027441932149,
        0.049067674327418764,
        0.4275550934302864,
        -0.8032075314806426,
        0.9891765099647812,
        -0.9415440651830184,
        0.6715589548470194,
        -0.24298017990325993,
    ],
    [
        0.19509032201612833,
        -0.5555702330196022,
        0.8314696123025455,
        -0.9807852804032307,
        0.9807852804032304,
        -0.831469612302545,
        0.5555702330196015,
        -0.19509032201612858,
        -0.19509032201613025,
        0.5555702330196028,
        -0.831469612302545,
        0.9807852804032309,
        -0.9807852804032297,
        0.8314696123025456,
        -0.5555702330196007,
        0.19509032201612425,
        0.1950903220161276,
        -0.5555702330196036,
        0.8314696123025475,
        -0.9807852804032303,
        0.9807852804032301,
        -0.8314696123025431,
        0.555570233019603,
        -0.1950903220161269,
        -0.19509032201612497,
        0.5555702330196073,
        -0.8314696123025459,
        0.9807852804032298,
        -0.9807852804032293,
        0.8314696123025446,
        -0.5555702330196052,
        0.19509032201612256,
    ],
    [
        0.14673047445536175,
        -0.4275550934302825,
        0.6715589548470188,
        -0.8577286100002723,
        0.9700312531945443,
        -0.9987954562051724,
        0.9415440651830204,
        -0.8032075314806446,
        0.5956993044924337,
        -0.33688985339221794,
        0.04906767432741668,
        0.24298017990326776,
        -0.5141027441932244,
        0.7409511253549605,
        -0.9039892931234439,
        0.989176509964781,
        -0.989176509964781,
        0.9039892931234439,
        -0.7409511253549558,
        0.5141027441932183,
        -0.24298017990326087,
        -0.04906767432742023,
        0.33688985339222133,
        -0.5956993044924337,
        0.8032075314806446,
        -0.9415440651830204,
        0.9987954562051723,
        -0.9700312531945448,
        0.8577286100002668,
        -0.6715589548470114,
        0.4275550934302745,
        -0.14673047445535428,
    ],
    [
        0.09801714032956077,
        -0.29028467725446244,
        0.471396736825998,
        -0.6343932841636454,
        0.7730104533627377,
        -0.8819212643483562,
        0.9569403357322094,
        -0.995184726672197,
        0.9951847266721968,
        -0.9569403357322088,
        0.8819212643483553,
        -0.7730104533627375,
        0.6343932841636439,
        -0.47139673682599326,
        0.2902846772544615,
        -0.09801714032955673,
        -0.09801714032956038,
        0.29028467725446505,
        -0.4713967368259965,
        0.6343932841636468,
        -0.7730104533627399,
        0.8819212643483553,
        -0.9569403357322078,
        0.9951847266721968,
        -0.9951847266721966,
        0.9569403357322073,
        -0.881921264348351,
        0.7730104533627387,
        -0.6343932841636453,
        0.47139673682599476,
        -0.29028467725445634,
        0.09801714032955137,
    ],
    [
        0.049067674327418126,
        -0.1467304744553623,
        0.24298017990326423,
        -0.336889853392221,
        0.4275550934302828,
        -0.5141027441932238,
        0.595699304492435,
        -0.6715589548470199,
        0.7409511253549602,
        -0.8032075314806458,
        0.8577286100002728,
        -0.9039892931234438,
        0.941544065183021,
        -0.9700312531945442,
        0.989176509964781,
        -0.9987954562051724,
        0.9987954562051724,
        -0.989176509964781,
        0.9700312531945441,
        -0.941544065183021,
        0.9039892931234437,
        -0.8577286100002726,
        0.8032075314806414,
        -0.7409511253549601,
        0.6715589548470144,
        -0.5956993044924349,
        0.5141027441932174,
        -0.4275550934302842,
        0.3368898533922158,
        -0.2429801799032666,
        0.14673047445535767,
        -0.049067674327421214,
    ],
    [
        6.123233995736766e-17,
        -1.8369701987210297e-16,
        3.061616997868383e-16,
        -4.286263797015736e-16,
        5.51091059616309e-16,
        -2.4499125789312946e-15,
        -9.803364199544708e-16,
        -2.6948419387607653e-15,
        -7.354070601250002e-16,
        -2.939771298590236e-15,
        -4.904777002955296e-16,
        -3.1847006584197066e-15,
        -2.45548340466059e-16,
        -3.4296300182491773e-15,
        -6.189806365883577e-19,
        -3.674559378078648e-15,
        2.443103791928823e-16,
        -3.919488737908119e-15,
        4.892397390223529e-16,
        -4.164418097737589e-15,
        7.839596456452825e-15,
        -4.40934745756706e-15,
        9.790984586812941e-16,
        2.4511505402044715e-15,
        8.329455176111767e-15,
        -4.899206177226001e-15,
        1.4689571783402355e-15,
        1.96129182054553e-15,
        8.819313895770708e-15,
        -5.389064896884942e-15,
        1.9588158979991767e-15,
        1.471433100886589e-15,
    ],
    [
        -0.04906767432741801,
        0.14673047445536194,
        -0.2429801799032628,
        0.3368898533922202,
        -0.4275550934302818,
        0.5141027441932227,
        -0.5956993044924338,
        0.6715589548470184,
        -0.7409511253549589,
        0.8032075314806444,
        -0.8577286100002696,
        0.9039892931234441,
        -0.9415440651830213,
        0.9700312531945442,
        -0.989176509964781,
        0.9987954562051724,
        -0.9987954562051724,
        0.9891765099647811,
        -0.9700312531945443,
        0.9415440651830214,
        -0.9039892931234473,
        0.8577286100002698,
        -0.8032075314806426,
        0.7409511253549568,
        -0.6715589548470162,
        0.5956993044924314,
        -0.51410274419322,
        0.42755509343028064,
        -0.336889853392219,
        0.24298017990326326,
        -0.14673047445536155,
        0.04906767432741827,
    ],
    [
        -0.09801714032956065,
        0.29028467725446205,
        -0.4713967368259975,
        0.6343932841636447,
        -0.773010453362737,
        0.8819212643483555,
        -0.9569403357322089,
        0.9951847266721968,
        -0.995184726672197,
        0.9569403357322095,
        -0.8819212643483565,
        0.7730104533627371,
        -0.634393284163649,
        0.4713967368259993,
        -0.2902846772544615,
        0.09801714032956405,
        0.0980171403295599,
        -0.29028467725445756,
        0.47139673682599564,
        -0.6343932841636404,
        0.773010453362739,
        -0.8819212643483545,
        0.9569403357322073,
        -0.9951847266721959,
        0.9951847266721969,
        -0.9569403357322102,
        0.8819212643483592,
        -0.7730104533627362,
        0.6343932841636479,
        -0.47139673682600425,
        0.2902846772544601,
        -0.09801714032956259,
    ],
    [
        -0.14673047445536164,
        0.42755509343028214,
        -0.6715589548470177,
        0.8577286100002719,
        -0.9700312531945441,
        0.9987954562051724,
        -0.9415440651830209,
        0.8032075314806457,
        -0.5956993044924354,
        0.33688985339222016,
        -0.04906767432741925,
        -0.2429801799032616,
        0.5141027441932188,
        -0.740951125354956,
        0.9039892931234439,
        -0.989176509964781,
        0.9891765099647811,
        -0.9039892931234442,
        0.7409511253549612,
        -0.5141027441932254,
        0.2429801799032623,
        0.04906767432741142,
        -0.33688985339221944,
        0.5956993044924263,
        -0.8032075314806432,
        0.9415440651830218,
        -0.9987954562051722,
        0.9700312531945438,
        -0.8577286100002759,
        0.6715589548470194,
        -0.42755509343029086,
        0.14673047445536544,
    ],
    [
        -0.1950903220161282,
        0.5555702330196018,
        -0.8314696123025451,
        0.9807852804032304,
        -0.9807852804032307,
        0.8314696123025448,
        -0.5555702330196027,
        0.19509032201613036,
        0.19509032201612822,
        -0.555570233019601,
        0.8314696123025456,
        -0.9807852804032295,
        0.9807852804032309,
        -0.8314696123025455,
        0.5555702330196066,
        -0.19509032201613144,
        -0.19509032201612714,
        0.555570233019603,
        -0.8314696123025429,
        0.9807852804032301,
        -0.9807852804032304,
        0.831469612302544,
        -0.5555702330196105,
        0.19509032201613602,
        0.19509032201612256,
        -0.5555702330195991,
        0.8314696123025443,
        -0.9807852804032305,
        0.98078528040323,
        -0.8314696123025506,
        0.5555702330196085,
        -0.1950903220161336,
    ],
    [
        -0.24298017990326387,
        0.6715589548470183,
        -0.9415440651830205,
        0.9891765099647811,
        -0.8032075314806455,
        0.4275550934302814,
        0.049067674327416926,
        -0.5141027441932223,
        0.8577286100002715,
        -0.9987954562051723,
        0.9039892931234453,
        -0.5956993044924329,
        0.146730474455363,
        0.33688985339221716,
        -0.7409511253549558,
        0.9700312531945441,
        -0.9700312531945443,
        0.7409511253549612,
        -0.33688985339221805,
        -0.14673047445536205,
        0.5956993044924321,
        -0.9039892931234419,
        0.9987954562051726,
        -0.8577286100002757,
        0.5141027441932292,
        -0.04906767432742855,
        -0.42755509343028375,
        0.8032075314806449,
        -0.9891765099647807,
        0.941544065183022,
        -0.6715589548470224,
        0.24298017990327087,
    ],
    [
        -0.29028467725446216,
        0.7730104533627367,
        -0.9951847266721969,
        0.8819212643483553,
        -0.4713967368259983,
        -0.09801714032956124,
        0.6343932841636444,
        -0.9569403357322088,
        0.9569403357322095,
        -0.6343932841636489,
        0.09801714032956356,
        0.47139673682599625,
        -0.8819212643483549,
        0.9951847266721968,
        -0.7730104533627398,
        0.2902846772544653,
        0.29028467725446083,
        -0.7730104533627368,
        0.9951847266721963,
        -0.8819212643483538,
        0.47139673682600036,
        0.09801714032955186,
        -0.6343932841636453,
        0.9569403357322072,
        -0.9569403357322082,
        0.6343932841636479,
        -0.09801714032956942,
        -0.47139673682599736,
        0.8819212643483522,
        -0.9951847266721966,
        0.773010453362739,
        -0.2902846772544709,
    ],
    [
        -0.33688985339221994,
        0.857728610000272,
        -0.9891765099647811,
        0.6715589548470182,
        -0.04906767432741852,
        -0.5956993044924338,
        0.9700312531945435,
        -0.9039892931234437,
        0.4275550934302851,
        0.24298017990326515,
        -0.8032075314806442,
        0.9987954562051726,
        -0.7409511253549584,
        0.14673047445536325,
        0.5141027441932183,
        -0.941544065183021,
        0.9415440651830214,
        -0.5141027441932254,
        -0.14673047445536205,
        0.7409511253549529,
        -0.9987954562051722,
        0.8032075314806449,
        -0.24298017990327322,
        -0.4275550934302776,
        0.9039892931234431,
        -0.9700312531945464,
        0.5956993044924377,
        0.0490676743274173,
        -0.6715589548470108,
        0.9891765099647801,
        -0.8577286100002726,
        0.33688985339223004,
    ],
    [
        -0.3826834323650897,
        0.9238795325112865,
        -0.9238795325112867,
        0.38268343236509067,
        0.38268343236508956,
        -0.923879532511287,
        0.9238795325112876,
        -0.3826834323650912,
        -0.382683432365089,
        0.9238795325112867,
        -0.9238795325112865,
        0.3826834323650885,
        0.3826834323650851,
        -0.9238795325112851,
        0.9238795325112881,
        -0.3826834323650924,
        -0.3826834323650813,
        0.9238795325112835,
        -0.9238795325112897,
        0.3826834323650962,
        0.382683432365084,
        -0.9238795325112846,
        0.9238795325112886,
        -0.3826834323650935,
        -0.38268343236508673,
        0.9238795325112857,
        -0.9238795325112874,
        0.3826834323650908,
        0.38268343236508945,
        -0.9238795325112868,
        0.9238795325112863,
        -0.38268343236508806,
    ],
    [
        -0.42755509343028186,
        0.970031253194544,
        -0.8032075314806454,
        0.0490676743274184,
        0.7409511253549591,
        -0.9891765099647809,
        0.5141027441932241,
        0.33688985339221783,
        -0.9415440651830214,
        0.8577286100002729,
        -0.1467304744553664,
        -0.6715589548470179,
        0.9987954562051726,
        -0.5956993044924334,
        -0.24298017990326087,
        0.9039892931234437,
        -0.9039892931234473,
        0.2429801799032623,
        0.5956993044924321,
        -0.9987954562051722,
        0.6715589548470242,
        0.14673047445536494,
        -0.8577286100002721,
        0.9415440651830218,
        -0.33688985339222594,
        -0.5141027441932137,
        0.9891765099647812,
        -0.7409511253549601,
        -0.04906767432741338,
        0.8032075314806403,
        -0.9700312531945466,
        0.427555093430282,
    ],
    [
        -0.4713967368259977,
        0.9951847266721969,
        -0.6343932841636454,
        -0.29028467725446255,
        0.9569403357322089,
        -0.7730104533627368,
        -0.09801714032956099,
        0.8819212643483553,
        -0.8819212643483565,
        0.09801714032956356,
        0.7730104533627351,
        -0.9569403357322097,
        0.29028467725446505,
        0.6343932841636434,
        -0.9951847266721972,
        0.4713967368259999,
        0.47139673682598915,
        -0.9951847266721966,
        0.6343932841636528,
        0.2902846772544601,
        -0.9569403357322062,
        0.7730104533627383,
        0.09801714032955137,
        -0.881921264348354,
        0.8819212643483594,
        -0.09801714032956259,
        -0.7730104533627312,
        0.9569403357322094,
        -0.2902846772544709,
        -0.6343932841636442,
        0.9951847266721977,
        -0.47139673682601163,
    ],
    [
        -0.5141027441932217,
        0.9987954562051724,
        -0.4275550934302827,
        -0.5956993044924326,
        0.9891765099647809,
        -0.33688985339221983,
        -0.6715589548470184,
        0.9700312531945441,
        -0.24298017990326443,
        -0.7409511253549561,
        0.9415440651830211,
        -0.14673047445536663,
        -0.8032075314806439,
        0.9039892931234457,
        -0.04906767432742023,
        -0.8577286100002726,
        0.8577286100002698,
        0.04906767432741142,
        -0.9039892931234419,
        0.8032075314806449,
        0.14673047445536494,
        -0.9415440651830181,
        0.7409511253549621,
        0.2429801799032628,
        -0.9700312531945445,
        0.6715589548470249,
        0.33688985339221483,
        -0.9891765099647807,
        0.5956993044924325,
        0.42755509343027315,
        -0.9987954562051727,
        0.5141027441932245,
    ],
    [
        -0.555570233019602,
        0.9807852804032304,
        -0.19509032201612803,
        -0.831469612302545,
        0.8314696123025448,
        0.19509032201612844,
        -0.9807852804032303,
        0.5555702330196061,
        0.5555702330196038,
        -0.9807852804032302,
        0.1950903220161276,
        0.8314696123025452,
        -0.8314696123025456,
        -0.19509032201612714,
        0.9807852804032301,
        -0.5555702330196102,
        -0.5555702330196056,
        0.9807852804032298,
        -0.19509032201612544,
        -0.8314696123025465,
        0.8314696123025443,
        0.1950903220161293,
        -0.9807852804032305,
        0.5555702330196024,
        0.5555702330196015,
        -0.9807852804032308,
        0.19509032201613025,
        0.8314696123025438,
        -0.831469612302547,
        -0.19509032201612447,
        0.9807852804032268,
        -0.5555702330196183,
    ],
    [
        -0.5956993044924334,
        0.9415440651830209,
        0.04906767432741742,
        -0.9700312531945441,
        0.5141027441932239,
        0.6715589548470184,
        -0.9039892931234437,
        -0.1467304744553635,
        0.9891765099647806,
        -0.4275550934302822,
        -0.740951125354956,
        0.8577286100002731,
        0.24298017990326445,
        -0.9987954562051722,
        0.33688985339222133,
        0.8032075314806414,
        -0.8032075314806426,
        -0.33688985339221944,
        0.9987954562051726,
        -0.24298017990327322,
        -0.8577286100002721,
        0.7409511253549621,
        0.42755509343027404,
        -0.9891765099647809,
        0.14673047445536544,
        0.9039892931234398,
        -0.6715589548470173,
        -0.5141027441932191,
        0.9700312531945459,
        -0.04906767432741582,
        -0.9415440651830153,
        0.5956993044924388,
    ],
    [
        -0.6343932841636454,
        0.881921264348355,
        0.29028467725446266,
        -0.9951847266721969,
        0.09801714032956282,
        0.9569403357322088,
        -0.47139673682599875,
        -0.7730104533627375,
        0.7730104533627371,
        0.47139673682599625,
        -0.9569403357322097,
        -0.09801714032955648,
        0.9951847266721964,
        -0.290284677254462,
        -0.8819212643483513,
        0.6343932841636472,
        0.6343932841636483,
        -0.8819212643483573,
        -0.2902846772544634,
        0.9951847266721976,
        -0.09801714032956209,
        -0.956940335732206,
        0.47139673682600125,
        0.7730104533627381,
        -0.7730104533627412,
        -0.4713967368259969,
        0.9569403357322115,
        0.09801714032955722,
        -0.9951847266721972,
        0.2902846772544681,
        0.8819212643483483,
        -0.6343932841636412,
    ],
    [
        -0.6715589548470184,
        0.8032075314806453,
        0.5141027441932213,
        -0.9039892931234434,
        -0.33688985339221816,
        0.970031253194544,
        0.1467304744553601,
        -0.9987954562051724,
        0.04906767432742292,
        0.9891765099647806,
        -0.24298017990326493,
        -0.941544065183021,
        0.42755509343028597,
        0.8577286100002709,
        -0.5956993044924337,
        -0.7409511253549601,
        0.7409511253549568,
        0.5956993044924263,
        -0.8577286100002757,
        -0.4275550934302776,
        0.9415440651830218,
        0.2429801799032628,
        -0.9891765099647809,
        -0.04906767432742072,
        0.998795456205172,
        -0.1467304744553693,
        -0.9700312531945426,
        0.3368898533922236,
        0.9039892931234365,
        -0.5141027441932339,
        -0.8032075314806376,
        0.671558954847026,
    ],
    [
        -0.7071067811865475,
        0.7071067811865477,
        0.7071067811865466,
        -0.7071067811865474,
        -0.7071067811865464,
        0.7071067811865476,
        0.707106781186546,
        -0.7071067811865479,
        -0.7071067811865458,
        0.7071067811865507,
        0.707106781186548,
        -0.7071067811865483,
        -0.7071067811865452,
        0.7071067811865511,
        0.7071067811865425,
        -0.7071067811865539,
        -0.7071067811865498,
        0.7071067811865467,
        0.707106781186547,
        -0.7071067811865495,
        -0.7071067811865442,
        0.7071067811865522,
        0.7071067811865415,
        -0.707106781186555,
        -0.7071067811865487,
        0.7071067811865477,
        0.707106781186546,
        -0.7071067811865606,
        -0.7071067811865432,
        0.7071067811865432,
        0.7071067811865405,
        -0.707106781186546,
    ],
    [
        -0.7409511253549589,
        0.5956993044924332,
        0.8577286100002719,
        -0.4275550934302813,
        -0.9415440651830203,
        0.2429801799032641,
        0.9891765099647806,
        -0.04906767432741925,
        -0.9987954562051724,
        -0.146730474455363,
        0.9700312531945451,
        0.33688985339221694,
        -0.9039892931234442,
        -0.5141027441932149,
        0.8032075314806446,
        0.6715589548470144,
        -0.6715589548470162,
        -0.8032075314806432,
        0.5141027441932292,
        0.9039892931234431,
        -0.33688985339222594,
        -0.9700312531945445,
        0.14673047445536544,
        0.998795456205172,
        0.04906767432741681,
        -0.989176509964782,
        -0.24298017990326515,
        0.9415440651830271,
        0.4275550934302855,
        -0.8577286100002731,
        -0.595699304492427,
        0.7409511253549683,
    ],
    [
        -0.773010453362737,
        0.471396736825998,
        0.9569403357322085,
        -0.0980171403295627,
        -0.995184726672197,
        -0.2902846772544621,
        0.8819212643483564,
        0.6343932841636439,
        -0.634393284163649,
        -0.8819212643483549,
        0.29028467725446505,
        0.9951847266721964,
        0.09801714032955966,
        -0.9569403357322078,
        -0.47139673682599215,
        0.7730104533627381,
        0.7730104533627387,
        -0.47139673682600386,
        -0.9569403357322082,
        0.09801714032955867,
        0.9951847266721977,
        0.29028467725445917,
        -0.8819212643483545,
        -0.6343932841636388,
        0.6343932841636487,
        0.8819212643483552,
        -0.2902846772544578,
        -0.995184726672195,
        -0.098017140329546,
        0.9569403357322118,
        0.4713967368259926,
        -0.7730104533627378,
    ],
    [
        -0.8032075314806448,
        0.33688985339222005,
        0.9987954562051724,
        0.24298017990326243,
        -0.8577286100002726,
        -0.7409511253549589,
        0.4275550934302851,
        0.9891765099647806,
        0.1467304744553596,
        -0.903989293123444,
        -0.6715589548470177,
        0.5141027441932221,
        0.9700312531945441,
        0.049067674327418764,
        -0.9415440651830204,
        -0.5956993044924349,
        0.5956993044924314,
        0.9415440651830218,
        -0.04906767432742855,
        -0.9700312531945464,
        -0.5141027441932137,
        0.6715589548470249,
        0.9039892931234398,
        -0.1467304744553693,
        -0.989176509964782,
        -0.42755509343027626,
        0.7409511253549631,
        0.857728610000262,
        -0.24298017990326848,
        -0.9987954562051733,
        -0.3368898533922167,
        0.8032075314806552,
    ],
    [
        -0.8314696123025453,
        0.19509032201612878,
        0.9807852804032307,
        0.5555702330196015,
        -0.5555702330196027,
        -0.9807852804032303,
        -0.19509032201612808,
        0.8314696123025471,
        0.8314696123025455,
        -0.19509032201613122,
        -0.9807852804032303,
        -0.5555702330196002,
        0.5555702330196071,
        0.9807852804032301,
        0.19509032201612303,
        -0.83146961230255,
        -0.8314696123025465,
        0.19509032201612928,
        0.9807852804032313,
        0.5555702330195958,
        -0.5555702330196114,
        -0.9807852804032304,
        -0.19509032201612497,
        0.8314696123025489,
        0.8314696123025397,
        -0.1950903220161413,
        -0.9807852804032337,
        -0.5555702330196093,
        0.5555702330195978,
        0.9807852804032309,
        0.1950903220161269,
        -0.8314696123025479,
    ],
    [
        -0.857728610000272,
        0.049067674327418154,
        0.9039892931234434,
        0.8032075314806447,
        -0.14673047445536216,
        -0.9415440651830209,
        -0.7409511253549563,
        0.242980179903268,
        0.9700312531945451,
        0.6715589548470151,
        -0.3368898533922243,
        -0.9891765099647817,
        -0.5956993044924352,
        0.4275550934302864,
        0.9987954562051723,
        0.5141027441932174,
        -0.51410274419322,
        -0.9987954562051722,
        -0.42755509343028375,
        0.5956993044924377,
        0.9891765099647812,
        0.33688985339221483,
        -0.6715589548470173,
        -0.9700312531945426,
        -0.24298017990326515,
        0.7409511253549631,
        0.9415440651830211,
        0.1467304744553698,
        -0.8032075314806528,
        -0.9039892931234407,
        -0.049067674327418764,
        0.8577286100002681,
    ],
    [
        -0.8819212643483549,
        -0.09801714032955997,
        0.7730104533627377,
        0.9569403357322089,
        0.2902846772544622,
        -0.6343932841636459,
        -0.9951847266721968,
        -0.47139673682599326,
        0.4713967368259993,
        0.9951847266721968,
        0.6343932841636434,
        -0.290284677254462,
        -0.9569403357322078,
        -0.7730104533627321,
        0.09801714032956502,
        0.8819212643483556,
        0.8819212643483558,
        0.09801714032955137,
        -0.7730104533627409,
        -0.9569403357322079,
        -0.29028467725446244,
        0.634393284163654,
        0.9951847266721962,
        0.4713967368259935,
        -0.47139673682601163,
        -0.9951847266721967,
        -0.634393284163638,
        0.29028467725445495,
        0.9569403357322098,
        0.7730104533627278,
        -0.0980171403295577,
        -0.8819212643483589,
    ],
    [
        -0.9039892931234433,
        -0.2429801799032628,
        0.5956993044924335,
        0.9987954562051724,
        0.6715589548470184,
        -0.14673047445536241,
        -0.857728610000271,
        -0.9415440651830213,
        -0.3368898533922206,
        0.5141027441932219,
        0.9891765099647811,
        0.740951125354958,
        -0.04906767432742757,
        -0.8032075314806426,
        -0.9700312531945448,
        -0.4275550934302842,
        0.42755509343028064,
        0.9700312531945438,
        0.8032075314806449,
        0.0490676743274173,
        -0.7409511253549601,
        -0.9891765099647807,
        -0.5141027441932191,
        0.3368898533922236,
        0.9415440651830271,
        0.857728610000262,
        0.1467304744553698,
        -0.6715589548470129,
        -0.9987954562051727,
        -0.595699304492438,
        0.24298017990325896,
        0.9039892931234415,
    ],
    [
        -0.9238795325112867,
        -0.3826834323650899,
        0.38268343236509067,
        0.9238795325112875,
        0.9238795325112868,
        0.3826834323650891,
        -0.38268343236509145,
        -0.9238795325112865,
        -0.9238795325112852,
        -0.3826834323650883,
        0.382683432365089,
        0.9238795325112882,
        0.9238795325112835,
        0.3826834323650908,
        -0.38268343236509306,
        -0.9238795325112898,
        -0.9238795325112873,
        -0.38268343236508673,
        0.3826834323650971,
        0.9238795325112862,
        0.9238795325112856,
        0.3826834323650826,
        -0.38268343236508806,
        -0.9238795325112878,
        -0.9238795325112893,
        -0.38268343236507857,
        0.38268343236509217,
        0.9238795325112841,
        0.9238795325112822,
        0.3826834323650876,
        -0.38268343236508306,
        -0.9238795325112912,
    ],
    [
        -0.9415440651830207,
        -0.5141027441932214,
        0.1467304744553618,
        0.7409511253549601,
        0.9987954562051724,
        0.8032075314806444,
        0.24298017990326515,
        -0.4275550934302822,
        -0.903989293123444,
        -0.9700312531945433,
        -0.5956993044924298,
        0.04906767432741681,
        0.6715589548470239,
        0.9891765099647812,
        0.8577286100002668,
        0.3368898533922158,
        -0.336889853392219,
        -0.8577286100002759,
        -0.9891765099647807,
        -0.6715589548470108,
        -0.04906767432741338,
        0.5956993044924325,
        0.9700312531945459,
        0.9039892931234365,
        0.4275550934302855,
        -0.24298017990326848,
        -0.8032075314806528,
        -0.9987954562051727,
        -0.7409511253549578,
        -0.14673047445535137,
        0.5141027441932381,
        0.9415440651830205,
    ],
    [
        -0.9569403357322088,
        -0.6343932841636448,
        -0.09801714032955972,
        0.4713967368259984,
        0.8819212643483563,
        0.9951847266721968,
        0.7730104533627352,
        0.2902846772544615,
        -0.2902846772544615,
        -0.7730104533627398,
        -0.9951847266721972,
        -0.8819212643483513,
        -0.47139673682599215,
        0.09801714032956502,
        0.6343932841636476,
        0.9569403357322092,
        0.9569403357322092,
        0.6343932841636476,
        0.09801714032955088,
        -0.4713967368260047,
        -0.8819212643483579,
        -0.9951847266721965,
        -0.7730104533627352,
        -0.2902846772544615,
        0.2902846772544615,
        0.7730104533627352,
        0.9951847266721965,
        0.8819212643483579,
        0.47139673682597966,
        -0.09801714032957916,
        -0.6343932841636586,
        -0.9569403357322133,
    ],
    [
        -0.970031253194544,
        -0.7409511253549593,
        -0.3368898533922201,
        0.14673047445536203,
        0.5956993044924352,
        0.9039892931234438,
        0.9987954562051724,
        0.8577286100002712,
        0.5141027441932186,
        0.04906767432741926,
        -0.42755509343028286,
        -0.8032075314806467,
        -0.9891765099647817,
        -0.9415440651830184,
        -0.6715589548470114,
        -0.2429801799032666,
        0.24298017990326326,
        0.6715589548470194,
        0.941544065183022,
        0.9891765099647801,
        0.8032075314806403,
        0.42755509343027315,
        -0.04906767432741582,
        -0.5141027441932339,
        -0.8577286100002731,
        -0.9987954562051733,
        -0.9039892931234407,
        -0.595699304492438,
        -0.14673047445535137,
        0.33688985339221855,
        0.740951125354969,
        0.9700312531945446,
    ],
    [
        -0.9807852804032304,
        -0.8314696123025451,
        -0.5555702330196015,
        -0.19509032201612858,
        0.19509032201613036,
        0.5555702330196061,
        0.8314696123025471,
        0.9807852804032309,
        0.9807852804032302,
        0.8314696123025451,
        0.555570233019603,
        0.19509032201613025,
        -0.19509032201613216,
        -0.5555702330196105,
        -0.8314696123025462,
        -0.980785280403232,
        -0.9807852804032305,
        -0.8314696123025421,
        -0.5555702330196044,
        -0.19509032201612497,
        0.19509032201613746,
        0.5555702330196032,
        0.8314696123025413,
        0.9807852804032302,
        0.9807852804032294,
        0.8314696123025391,
        0.5555702330195881,
        0.1950903220161336,
        -0.1950903220161288,
        -0.5555702330196076,
        -0.8314696123025522,
        -0.9807852804032341,
    ],
    [
        -0.989176509964781,
        -0.9039892931234431,
        -0.7409511253549592,
        -0.5141027441932225,
        -0.24298017990326207,
        0.04906767432742268,
        0.3368898533922204,
        0.5956993044924359,
        0.8032075314806442,
        0.9415440651830214,
        0.9987954562051726,
        0.9700312531945422,
        0.8577286100002706,
        0.6715589548470194,
        0.4275550934302745,
        0.14673047445535767,
        -0.14673047445536155,
        -0.42755509343029086,
        -0.6715589548470224,
        -0.8577286100002726,
        -0.9700312531945466,
        -0.9987954562051727,
        -0.9415440651830153,
        -0.8032075314806376,
        -0.595699304492427,
        -0.3368898533922167,
        -0.049067674327418764,
        0.24298017990325896,
        0.5141027441932381,
        0.740951125354969,
        0.9039892931234478,
        0.9891765099647819,
    ],
    [
        -0.9951847266721968,
        -0.9569403357322085,
        -0.8819212643483547,
        -0.7730104533627357,
        -0.6343932841636443,
        -0.4713967368259935,
        -0.2902846772544582,
        -0.09801714032955673,
        0.09801714032956405,
        0.2902846772544653,
        0.4713967368259999,
        0.6343932841636472,
        0.7730104533627381,
        0.8819212643483556,
        0.9569403357322092,
        0.9951847266721969,
        0.9951847266721969,
        0.9569403357322089,
        0.8819212643483554,
        0.7730104533627378,
        0.6343932841636468,
        0.47139673682599953,
        0.29028467725445123,
        0.09801714032956356,
        -0.09801714032957136,
        -0.2902846772544587,
        -0.4713967368260064,
        -0.6343932841636418,
        -0.7730104533627428,
        -0.8819212643483524,
        -0.9569403357322113,
        -0.9951847266721963,
    ],
    [
        -0.9987954562051724,
        -0.989176509964781,
        -0.9700312531945441,
        -0.9415440651830203,
        -0.9039892931234428,
        -0.8577286100002696,
        -0.8032075314806442,
        -0.740951125354956,
        -0.6715589548470177,
        -0.5956993044924298,
        -0.5141027441932149,
        -0.42755509343028464,
        -0.33688985339221944,
        -0.24298017990325993,
        -0.14673047445535428,
        -0.049067674327421214,
        0.04906767432741827,
        0.14673047445536544,
        0.24298017990327087,
        0.33688985339223004,
        0.427555093430282,
        0.5141027441932245,
        0.5956993044924388,
        0.671558954847026,
        0.7409511253549683,
        0.8032075314806552,
        0.8577286100002681,
        0.9039892931234415,
        0.9415440651830205,
        0.9700312531945446,
        0.9891765099647819,
        0.9987954562051728,
    ],
    [
        -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0,
    ],
    [
        -0.9987954562051724,
        -0.9891765099647811,
        -0.9700312531945443,
        -0.9415440651830209,
        -0.9039892931234437,
        -0.857728610000271,
        -0.803207531480646,
        -0.7409511253549584,
        -0.6715589548470208,
        -0.5956993044924335,
        -0.5141027441932254,
        -0.4275550934302833,
        -0.33688985339221855,
        -0.24298017990327322,
        -0.14673047445536833,
        -0.0490676743274217,
        0.0490676743274173,
        0.14673047445536397,
        0.24298017990325518,
        0.3368898533922144,
        0.42755509343027936,
        0.5141027441932094,
        0.5956993044924357,
        0.6715589548470122,
        0.7409511253549459,
        0.8032075314806435,
        0.857728610000265,
        0.9039892931234449,
        0.9415440651830181,
        0.9700312531945394,
        0.9891765099647807,
        0.9987954562051717,
    ],
    [
        -0.9951847266721969,
        -0.9569403357322087,
        -0.8819212643483562,
        -0.7730104533627368,
        -0.6343932841636459,
        -0.47139673682599587,
        -0.2902846772544613,
        -0.09801714032956038,
        0.0980171403295599,
        0.29028467725446083,
        0.47139673682598915,
        0.6343932841636483,
        0.7730104533627387,
        0.8819212643483558,
        0.9569403357322092,
        0.9951847266721969,
        0.9951847266721969,
        0.9569403357322094,
        0.8819212643483564,
        0.7730104533627393,
        0.63439328416366,
        0.47139673682599,
        0.29028467725445495,
        0.0980171403295538,
        -0.09801714032956649,
        -0.29028467725446716,
        -0.47139673682600125,
        -0.6343932841636479,
        -0.7730104533627383,
        -0.8819212643483556,
        -0.9569403357322089,
        -0.9951847266721968,
    ],
    [
        -0.989176509964781,
        -0.9039892931234433,
        -0.74095112535496,
        -0.514102744193224,
        -0.2429801799032642,
        0.049067674327419986,
        0.3368898533922174,
        0.5956993044924329,
        0.8032075314806417,
        0.9415440651830198,
        0.9987954562051724,
        0.9700312531945453,
        0.8577286100002701,
        0.6715589548470191,
        0.4275550934302873,
        0.1467304744553722,
        -0.14673047445536058,
        -0.4275550934302767,
        -0.6715589548470104,
        -0.857728610000264,
        -0.9700312531945425,
        -0.9987954562051722,
        -0.9415440651830261,
        -0.8032075314806487,
        -0.5956993044924309,
        -0.33688985339223515,
        -0.049067674327424635,
        0.2429801799032666,
        0.5141027441932078,
        0.7409511253549544,
        0.9039892931234444,
        0.9891765099647786,
    ],
    [
        -0.9807852804032304,
        -0.8314696123025456,
        -0.5555702330196026,
        -0.19509032201613025,
        0.19509032201612822,
        0.5555702330196038,
        0.8314696123025455,
        0.9807852804032302,
        0.980785280403231,
        0.8314696123025477,
        0.5555702330196073,
        0.1950903220161288,
        -0.1950903220161192,
        -0.5555702330195991,
        -0.8314696123025462,
        -0.9807852804032291,
        -0.9807852804032308,
        -0.8314696123025509,
        -0.5555702330196061,
        -0.1950903220161413,
        0.19509032201613458,
        0.5555702330196003,
        0.8314696123025391,
        0.9807852804032267,
        0.9807852804032304,
        0.83146961230255,
        0.5555702330196166,
        0.19509032201612592,
        -0.1950903220161221,
        -0.5555702330195897,
        -0.8314696123025479,
        -0.9807852804032297,
    ],
    [
        -0.970031253194544,
        -0.7409511253549599,
        -0.33688985339221955,
        0.14673047445536033,
        0.5956993044924335,
        0.9039892931234441,
        0.9987954562051726,
        0.8577286100002731,
        0.5141027441932221,
        0.04906767432741681,
        -0.42755509343028464,
        -0.8032075314806391,
        -0.9891765099647798,
        -0.9415440651830229,
        -0.671558954847022,
        -0.24298017990326706,
        0.2429801799032623,
        0.6715589548470184,
        0.9415440651830214,
        0.9891765099647826,
        0.8032075314806505,
        0.4275550934302891,
        -0.049067674327411916,
        -0.5141027441932179,
        -0.8577286100002706,
        -0.9987954562051723,
        -0.9039892931234431,
        -0.5956993044924317,
        -0.14673047445535817,
        0.336889853392225,
        0.7409511253549447,
        0.9700312531945392,
    ],
    [
        -0.9569403357322089,
        -0.6343932841636454,
        -0.0980171403295627,
        0.4713967368259969,
        0.8819212643483553,
        0.9951847266721967,
        0.7730104533627373,
        0.29028467725446505,
        -0.29028467725445756,
        -0.7730104533627368,
        -0.9951847266721966,
        -0.8819212643483573,
        -0.47139673682600386,
        0.09801714032955137,
        0.6343932841636476,
        0.9569403357322089,
        0.9569403357322094,
        0.6343932841636487,
        0.09801714032956697,
        -0.47139673682599,
        -0.8819212643483566,
        -0.9951847266721981,
        -0.7730104533627378,
        -0.2902846772544793,
        0.29028467725445684,
        0.7730104533627409,
        0.9951847266721959,
        0.8819212643483543,
        0.47139673682601074,
        -0.0980171403295577,
        -0.6343932841636305,
        -0.9569403357322067,
    ],
    [
        -0.9415440651830208,
        -0.514102744193222,
        0.14673047445536058,
        0.7409511253549589,
        0.9987954562051723,
        0.8032075314806439,
        0.24298017990326823,
        -0.4275550934302789,
        -0.9039892931234422,
        -0.9700312531945444,
        -0.5956993044924396,
        0.04906767432741828,
        0.671558954847014,
        0.9891765099647812,
        0.8577286100002741,
        0.3368898533922296,
        -0.33688985339221805,
        -0.8577286100002678,
        -0.989176509964781,
        -0.6715589548470231,
        -0.049067674327430505,
        0.5956993044924184,
        0.9700312531945449,
        0.9039892931234444,
        0.42755509343028997,
        -0.24298017990324947,
        -0.8032075314806324,
        -0.9987954562051723,
        -0.7409511253549624,
        -0.1467304744553727,
        0.514102744193207,
        0.9415440651830225,
    ],
    [
        -0.9238795325112868,
        -0.38268343236509056,
        0.38268343236508956,
        0.9238795325112868,
        0.9238795325112877,
        0.3826834323650883,
        -0.3826834323650885,
        -0.9238795325112851,
        -0.9238795325112868,
        -0.3826834323650926,
        0.3826834323650908,
        0.9238795325112833,
        0.9238795325112886,
        0.38268343236509034,
        -0.3826834323650799,
        -0.9238795325112843,
        -0.9238795325112876,
        -0.38268343236508806,
        0.3826834323650822,
        0.9238795325112852,
        0.9238795325112867,
        0.3826834323650858,
        -0.38268343236507135,
        -0.9238795325112807,
        -0.9238795325112912,
        -0.38268343236509667,
        0.38268343236508673,
        0.9238795325112871,
        0.9238795325112849,
        0.38268343236510755,
        -0.38268343236507585,
        -0.9238795325112826,
    ],
    [
        -0.9039892931234434,
        -0.24298017990326348,
        0.5956993044924339,
        0.9987954562051723,
        0.6715589548470174,
        -0.14673047445536325,
        -0.8577286100002693,
        -0.9415440651830225,
        -0.33688985339222455,
        0.5141027441932179,
        0.9891765099647803,
        0.7409511253549618,
        -0.04906767432741436,
        -0.8032075314806429,
        -0.9700312531945448,
        -0.42755509343028464,
        0.4275550934302798,
        0.9700312531945434,
        0.8032075314806546,
        0.04906767432741974,
        -0.7409511253549486,
        -0.9891765099647811,
        -0.5141027441932347,
        0.33688985339221944,
        0.9415440651830159,
        0.8577286100002721,
        0.1467304744553756,
        -0.6715589548470188,
        -0.9987954562051731,
        -0.5956993044924325,
        0.24298017990325138,
        0.903989293123444,
    ],
    [
        -0.881921264348355,
        -0.09801714032956069,
        0.7730104533627358,
        0.9569403357322095,
        0.29028467725446433,
        -0.6343932841636466,
        -0.9951847266721972,
        -0.4713967368259965,
        0.47139673682599564,
        0.9951847266721963,
        0.6343932841636528,
        -0.2902846772544634,
        -0.9569403357322082,
        -0.7730104533627409,
        0.09801714032955088,
        0.8819212643483554,
        0.8819212643483564,
        0.09801714032956697,
        -0.7730104533627397,
        -0.9569403357322087,
        -0.2902846772544653,
        0.6343932841636404,
        0.9951847266721979,
        0.4713967368260099,
        -0.4713967368259822,
        -0.9951847266721948,
        -0.6343932841636426,
        0.29028467725446244,
        0.9569403357322078,
        0.7730104533627414,
        -0.0980171403295499,
        -0.8819212643483483,
    ],
    [
        -0.8577286100002721,
        0.04906767432741742,
        0.9039892931234429,
        0.8032075314806457,
        -0.1467304744553635,
        -0.9415440651830213,
        -0.7409511253549584,
        0.24298017990326445,
        0.9700312531945441,
        0.6715589548470239,
        -0.33688985339221944,
        -0.9891765099647798,
        -0.5956993044924345,
        0.42755509343027404,
        0.9987954562051723,
        0.5141027441932301,
        -0.5141027441932191,
        -0.998795456205173,
        -0.4275550934302855,
        0.5956993044924357,
        0.9891765099647838,
        0.33688985339223143,
        -0.6715589548470144,
        -0.9700312531945436,
        -0.2429801799032837,
        0.7409511253549499,
        0.9415440651830231,
        0.14673047445536203,
        -0.8032075314806318,
        -0.9039892931234499,
        -0.04906767432742659,
        0.8577286100002711,
    ],
    [
        -0.8314696123025455,
        0.1950903220161272,
        0.9807852804032304,
        0.5555702330196028,
        -0.555570233019601,
        -0.9807852804032302,
        -0.19509032201613122,
        0.8314696123025451,
        0.8314696123025477,
        -0.19509032201611967,
        -0.9807852804032293,
        -0.5555702330196048,
        0.555570233019602,
        0.98078528040323,
        0.195090322016137,
        -0.8314696123025419,
        -0.831469612302547,
        0.19509032201612786,
        0.9807852804032281,
        0.5555702330195978,
        -0.5555702330195971,
        -0.9807852804032339,
        -0.1950903220161288,
        0.8314696123025386,
        0.8314696123025425,
        -0.1950903220161221,
        -0.980785280403227,
        -0.5555702330196027,
        0.5555702330195922,
        0.9807852804032294,
        0.19509032201613458,
        -0.8314696123025354,
    ],
    [
        -0.8032075314806449,
        0.3368898533922202,
        0.9987954562051724,
        0.2429801799032641,
        -0.8577286100002696,
        -0.7409511253549583,
        0.4275550934302822,
        0.9891765099647811,
        0.14673047445537077,
        -0.903989293123445,
        -0.6715589548470162,
        0.5141027441932233,
        0.9700312531945438,
        0.04906767432741827,
        -0.9415440651830204,
        -0.5956993044924352,
        0.5956993044924306,
        0.9415440651830271,
        -0.0490676743274266,
        -0.9700312531945459,
        -0.5141027441932162,
        0.6715589548470224,
        0.9039892931234415,
        -0.14673047445536494,
        -0.9891765099647812,
        -0.4275550934302811,
        0.7409511253549591,
        0.8577286100002726,
        -0.24298017990326182,
        -0.9987954562051722,
        -0.33688985339222405,
        0.8032075314806417,
    ],
    [
        -0.7730104533627371,
        0.47139673682599736,
        0.9569403357322094,
        -0.09801714032956099,
        -0.9951847266721968,
        -0.2902846772544613,
        0.8819212643483533,
        0.6343932841636468,
        -0.6343932841636404,
        -0.8819212643483538,
        0.2902846772544601,
        0.9951847266721976,
        0.09801714032955867,
        -0.9569403357322079,
        -0.4713967368260047,
        0.7730104533627378,
        0.7730104533627393,
        -0.47139673682599,
        -0.9569403357322087,
        0.0980171403295421,
        0.9951847266721959,
        0.29028467725446244,
        -0.881921264348346,
        -0.6343932841636533,
        0.6343932841636449,
        0.8819212643483644,
        -0.2902846772544521,
        -0.995184726672197,
        -0.09801714032958111,
        0.9569403357322056,
        0.47139673682599953,
        -0.7730104533627234,
    ],
    [
        -0.7409511253549591,
        0.5956993044924327,
        0.8577286100002725,
        -0.4275550934302798,
        -0.9415440651830222,
        0.24298017990326493,
        0.9891765099647811,
        -0.049067674327415586,
        -0.9987954562051722,
        -0.14673047445536058,
        0.9700312531945421,
        0.3368898533922222,
        -0.9039892931234447,
        -0.5141027441932267,
        0.8032075314806446,
        0.6715589548470253,
        -0.6715589548470154,
        -0.803207531480644,
        0.5141027441932153,
        0.9039892931234503,
        -0.33688985339222316,
        -0.9700312531945453,
        0.1467304744553475,
        0.9987954562051722,
        0.0490676743274217,
        -0.9891765099647791,
        -0.24298017990328463,
        0.9415440651830201,
        0.42755509343029174,
        -0.857728610000262,
        -0.5956993044924334,
        0.7409511253549532,
    ],
];

#[allow(clippy::unreadable_literal)]
const COS_N36: [[f32; 36]; 18] = [
    [
        0.6755902,
        0.6087614,
        0.53729963,
        0.4617486,
        0.38268343,
        0.3007058,
        0.21643962,
        0.13052619,
        0.043619387,
        -0.043619387,
        -0.13052619,
        -0.21643962,
        -0.3007058,
        -0.38268343,
        -0.4617486,
        -0.53729963,
        -0.6087614,
        -0.6755902,
        -0.7372773,
        -0.7933533,
        -0.8433914,
        -0.8870108,
        -0.9238795,
        -0.95371693,
        -0.976296,
        -0.9914449,
        -0.99904823,
        -0.99904823,
        -0.9914449,
        -0.976296,
        -0.95371693,
        -0.9238795,
        -0.8870108,
        -0.8433914,
        -0.7933533,
        -0.7372773,
    ],
    [
        -0.7933533,
        -0.9238795,
        -0.9914449,
        -0.9914449,
        -0.9238795,
        -0.7933533,
        -0.6087614,
        -0.38268343,
        -0.13052619,
        0.13052619,
        0.38268343,
        0.6087614,
        0.7933533,
        0.9238795,
        0.9914449,
        0.9914449,
        0.9238795,
        0.7933533,
        0.6087614,
        0.38268343,
        0.13052619,
        -0.13052619,
        -0.38268343,
        -0.6087614,
        -0.7933533,
        -0.9238795,
        -0.9914449,
        -0.9914449,
        -0.9238795,
        -0.7933533,
        -0.6087614,
        -0.38268343,
        -0.13052619,
        0.13052619,
        0.38268343,
        0.6087614,
    ],
    [
        -0.53729963,
        -0.13052619,
        0.3007058,
        0.6755902,
        0.9238795,
        0.99904823,
        0.8870108,
        0.6087614,
        0.21643962,
        -0.21643962,
        -0.6087614,
        -0.8870108,
        -0.99904823,
        -0.9238795,
        -0.6755902,
        -0.3007058,
        0.13052619,
        0.53729963,
        0.8433914,
        0.9914449,
        0.95371693,
        0.7372773,
        0.38268343,
        -0.043619387,
        -0.4617486,
        -0.7933533,
        -0.976296,
        -0.976296,
        -0.7933533,
        -0.4617486,
        -0.043619387,
        0.38268343,
        0.7372773,
        0.95371693,
        0.9914449,
        0.8433914,
    ],
    [
        0.8870108,
        0.9914449,
        0.7372773,
        0.21643962,
        -0.38268343,
        -0.8433914,
        -0.99904823,
        -0.7933533,
        -0.3007058,
        0.3007058,
        0.7933533,
        0.99904823,
        0.8433914,
        0.38268343,
        -0.21643962,
        -0.7372773,
        -0.9914449,
        -0.8870108,
        -0.4617486,
        0.13052619,
        0.6755902,
        0.976296,
        0.9238795,
        0.53729963,
        -0.043619387,
        -0.6087614,
        -0.95371693,
        -0.95371693,
        -0.6087614,
        -0.043619387,
        0.53729963,
        0.9238795,
        0.976296,
        0.6755902,
        0.13052619,
        -0.4617486,
    ],
    [
        0.38268343,
        -0.38268343,
        -0.9238795,
        -0.9238795,
        -0.38268343,
        0.38268343,
        0.9238795,
        0.9238795,
        0.38268343,
        -0.38268343,
        -0.9238795,
        -0.9238795,
        -0.38268343,
        0.38268343,
        0.9238795,
        0.9238795,
        0.38268343,
        -0.38268343,
        -0.9238795,
        -0.9238795,
        -0.38268343,
        0.38268343,
        0.9238795,
        0.9238795,
        0.38268343,
        -0.38268343,
        -0.9238795,
        -0.9238795,
        -0.38268343,
        0.38268343,
        0.9238795,
        0.9238795,
        0.38268343,
        -0.38268343,
        -0.9238795,
        -0.9238795,
    ],
    [
        -0.95371693,
        -0.7933533,
        0.043619387,
        0.8433914,
        0.9238795,
        0.21643962,
        -0.6755902,
        -0.9914449,
        -0.4617486,
        0.4617486,
        0.9914449,
        0.6755902,
        -0.21643962,
        -0.9238795,
        -0.8433914,
        -0.043619387,
        0.7933533,
        0.95371693,
        0.3007058,
        -0.6087614,
        -0.99904823,
        -0.53729963,
        0.38268343,
        0.976296,
        0.7372773,
        -0.13052619,
        -0.8870108,
        -0.8870108,
        -0.13052619,
        0.7372773,
        0.976296,
        0.38268343,
        -0.53729963,
        -0.99904823,
        -0.6087614,
        0.3007058,
    ],
    [
        -0.21643962,
        0.7933533,
        0.8870108,
        -0.043619387,
        -0.9238795,
        -0.7372773,
        0.3007058,
        0.9914449,
        0.53729963,
        -0.53729963,
        -0.9914449,
        -0.3007058,
        0.7372773,
        0.9238795,
        0.043619387,
        -0.8870108,
        -0.7933533,
        0.21643962,
        0.976296,
        0.6087614,
        -0.4617486,
        -0.99904823,
        -0.38268343,
        0.6755902,
        0.95371693,
        0.13052619,
        -0.8433914,
        -0.8433914,
        0.13052619,
        0.95371693,
        0.6755902,
        -0.38268343,
        -0.99904823,
        -0.4617486,
        0.6087614,
        0.976296,
    ],
    [
        0.9914449,
        0.38268343,
        -0.7933533,
        -0.7933533,
        0.38268343,
        0.9914449,
        0.13052619,
        -0.9238795,
        -0.6087614,
        0.6087614,
        0.9238795,
        -0.13052619,
        -0.9914449,
        -0.38268343,
        0.7933533,
        0.7933533,
        -0.38268343,
        -0.9914449,
        -0.13052619,
        0.9238795,
        0.6087614,
        -0.6087614,
        -0.9238795,
        0.13052619,
        0.9914449,
        0.38268343,
        -0.7933533,
        -0.7933533,
        0.38268343,
        0.9914449,
        0.13052619,
        -0.9238795,
        -0.6087614,
        0.6087614,
        0.9238795,
        -0.13052619,
    ],
    [
        0.043619387,
        -0.9914449,
        -0.21643962,
        0.95371693,
        0.38268343,
        -0.8870108,
        -0.53729963,
        0.7933533,
        0.6755902,
        -0.6755902,
        -0.7933533,
        0.53729963,
        0.8870108,
        -0.38268343,
        -0.95371693,
        0.21643962,
        0.9914449,
        -0.043619387,
        -0.99904823,
        -0.13052619,
        0.976296,
        0.3007058,
        -0.9238795,
        -0.4617486,
        0.8433914,
        0.6087614,
        -0.7372773,
        -0.7372773,
        0.6087614,
        0.8433914,
        -0.4617486,
        -0.9238795,
        0.3007058,
        0.976296,
        -0.13052619,
        -0.99904823,
    ],
    [
        -0.99904823,
        0.13052619,
        0.976296,
        -0.3007058,
        -0.9238795,
        0.4617486,
        0.8433914,
        -0.6087614,
        -0.7372773,
        0.7372773,
        0.6087614,
        -0.8433914,
        -0.4617486,
        0.9238795,
        0.3007058,
        -0.976296,
        -0.13052619,
        0.99904823,
        -0.043619387,
        -0.9914449,
        0.21643962,
        0.95371693,
        -0.38268343,
        -0.8870108,
        0.53729963,
        0.7933533,
        -0.6755902,
        -0.6755902,
        0.7933533,
        0.53729963,
        -0.8870108,
        -0.38268343,
        0.95371693,
        0.21643962,
        -0.9914449,
        -0.043619387,
    ],
    [
        0.13052619,
        0.9238795,
        -0.6087614,
        -0.6087614,
        0.9238795,
        0.13052619,
        -0.9914449,
        0.38268343,
        0.7933533,
        -0.7933533,
        -0.38268343,
        0.9914449,
        -0.13052619,
        -0.9238795,
        0.6087614,
        0.6087614,
        -0.9238795,
        -0.13052619,
        0.9914449,
        -0.38268343,
        -0.7933533,
        0.7933533,
        0.38268343,
        -0.9914449,
        0.13052619,
        0.9238795,
        -0.6087614,
        -0.6087614,
        0.9238795,
        0.13052619,
        -0.9914449,
        0.38268343,
        0.7933533,
        -0.7933533,
        -0.38268343,
        0.9914449,
    ],
    [
        0.976296,
        -0.6087614,
        -0.4617486,
        0.99904823,
        -0.38268343,
        -0.6755902,
        0.95371693,
        -0.13052619,
        -0.8433914,
        0.8433914,
        0.13052619,
        -0.95371693,
        0.6755902,
        0.38268343,
        -0.99904823,
        0.4617486,
        0.6087614,
        -0.976296,
        0.21643962,
        0.7933533,
        -0.8870108,
        -0.043619387,
        0.9238795,
        -0.7372773,
        -0.3007058,
        0.9914449,
        -0.53729963,
        -0.53729963,
        0.9914449,
        -0.3007058,
        -0.7372773,
        0.9238795,
        -0.043619387,
        -0.8870108,
        0.7933533,
        0.21643962,
    ],
    [
        -0.3007058,
        -0.6087614,
        0.99904823,
        -0.53729963,
        -0.38268343,
        0.976296,
        -0.7372773,
        -0.13052619,
        0.8870108,
        -0.8870108,
        0.13052619,
        0.7372773,
        -0.976296,
        0.38268343,
        0.53729963,
        -0.99904823,
        0.6087614,
        0.3007058,
        -0.95371693,
        0.7933533,
        0.043619387,
        -0.8433914,
        0.9238795,
        -0.21643962,
        -0.6755902,
        0.9914449,
        -0.4617486,
        -0.4617486,
        0.9914449,
        -0.6755902,
        -0.21643962,
        0.9238795,
        -0.8433914,
        0.043619387,
        0.7933533,
        -0.95371693,
    ],
    [
        -0.9238795,
        0.9238795,
        -0.38268343,
        -0.38268343,
        0.9238795,
        -0.9238795,
        0.38268343,
        0.38268343,
        -0.9238795,
        0.9238795,
        -0.38268343,
        -0.38268343,
        0.9238795,
        -0.9238795,
        0.38268343,
        0.38268343,
        -0.9238795,
        0.9238795,
        -0.38268343,
        -0.38268343,
        0.9238795,
        -0.9238795,
        0.38268343,
        0.38268343,
        -0.9238795,
        0.9238795,
        -0.38268343,
        -0.38268343,
        0.9238795,
        -0.9238795,
        0.38268343,
        0.38268343,
        -0.9238795,
        0.9238795,
        -0.38268343,
        -0.38268343,
    ],
    [
        0.4617486,
        0.13052619,
        -0.6755902,
        0.976296,
        -0.9238795,
        0.53729963,
        0.043619387,
        -0.6087614,
        0.95371693,
        -0.95371693,
        0.6087614,
        -0.043619387,
        -0.53729963,
        0.9238795,
        -0.976296,
        0.6755902,
        -0.13052619,
        -0.4617486,
        0.8870108,
        -0.9914449,
        0.7372773,
        -0.21643962,
        -0.38268343,
        0.8433914,
        -0.99904823,
        0.7933533,
        -0.3007058,
        -0.3007058,
        0.7933533,
        -0.99904823,
        0.8433914,
        -0.38268343,
        -0.21643962,
        0.7372773,
        -0.9914449,
        0.8870108,
    ],
    [
        0.8433914,
        -0.9914449,
        0.95371693,
        -0.7372773,
        0.38268343,
        0.043619387,
        -0.4617486,
        0.7933533,
        -0.976296,
        0.976296,
        -0.7933533,
        0.4617486,
        -0.043619387,
        -0.38268343,
        0.7372773,
        -0.95371693,
        0.9914449,
        -0.8433914,
        0.53729963,
        -0.13052619,
        -0.3007058,
        0.6755902,
        -0.9238795,
        0.99904823,
        -0.8870108,
        0.6087614,
        -0.21643962,
        -0.21643962,
        0.6087614,
        -0.8870108,
        0.99904823,
        -0.9238795,
        0.6755902,
        -0.3007058,
        -0.13052619,
        0.53729963,
    ],
    [
        -0.6087614,
        0.38268343,
        -0.13052619,
        -0.13052619,
        0.38268343,
        -0.6087614,
        0.7933533,
        -0.9238795,
        0.9914449,
        -0.9914449,
        0.9238795,
        -0.7933533,
        0.6087614,
        -0.38268343,
        0.13052619,
        0.13052619,
        -0.38268343,
        0.6087614,
        -0.7933533,
        0.9238795,
        -0.9914449,
        0.9914449,
        -0.9238795,
        0.7933533,
        -0.6087614,
        0.38268343,
        -0.13052619,
        -0.13052619,
        0.38268343,
        -0.6087614,
        0.7933533,
        -0.9238795,
        0.9914449,
        -0.9914449,
        0.9238795,
        -0.7933533,
    ],
    [
        -0.7372773,
        0.7933533,
        -0.8433914,
        0.8870108,
        -0.9238795,
        0.95371693,
        -0.976296,
        0.9914449,
        -0.99904823,
        0.99904823,
        -0.9914449,
        0.976296,
        -0.95371693,
        0.9238795,
        -0.8870108,
        0.8433914,
        -0.7933533,
        0.7372773,
        -0.6755902,
        0.6087614,
        -0.53729963,
        0.4617486,
        -0.38268343,
        0.3007058,
        -0.21643962,
        0.13052619,
        -0.043619387,
        -0.043619387,
        0.13052619,
        -0.21643962,
        0.3007058,
        -0.38268343,
        0.4617486,
        -0.53729963,
        0.6087614,
        -0.6755902,
    ],
];

#[allow(clippy::unreadable_literal)]
const COS_N12: [[f32; 12]; 6] = [
    [
        0.6087614,
        0.38268343,
        0.13052619,
        -0.13052619,
        -0.38268343,
        -0.6087614,
        -0.7933533,
        -0.9238795,
        -0.9914449,
        -0.9914449,
        -0.9238795,
        -0.7933533,
    ],
    [
        -0.9238795,
        -0.9238795,
        -0.38268343,
        0.38268343,
        0.9238795,
        0.9238795,
        0.38268343,
        -0.38268343,
        -0.9238795,
        -0.9238795,
        -0.38268343,
        0.38268343,
    ],
    [
        -0.13052619,
        0.9238795,
        0.6087614,
        -0.6087614,
        -0.9238795,
        0.13052619,
        0.9914449,
        0.38268343,
        -0.7933533,
        -0.7933533,
        0.38268343,
        0.9914449,
    ],
    [
        0.9914449,
        -0.38268343,
        -0.7933533,
        0.7933533,
        0.38268343,
        -0.9914449,
        0.13052619,
        0.9238795,
        -0.6087614,
        -0.6087614,
        0.9238795,
        0.13052619,
    ],
    [
        -0.38268343,
        -0.38268343,
        0.9238795,
        -0.9238795,
        0.38268343,
        0.38268343,
        -0.9238795,
        0.9238795,
        -0.38268343,
        -0.38268343,
        0.9238795,
        -0.9238795,
    ],
    [
        -0.7933533,
        0.9238795,
        -0.9914449,
        0.9914449,
        -0.9238795,
        0.7933533,
        -0.6087614,
        0.38268343,
        -0.13052619,
        -0.13052619,
        0.38268343,
        -0.6087614,
    ],
];

#[allow(clippy::unreadable_literal)]
const IMDCT_WIN: [[f32; 36]; 4] = [
    [
        0.043619387,
        0.13052619,
        0.21643962,
        0.3007058,
        0.38268343,
        0.4617486,
        0.53729963,
        0.6087614,
        0.6755902,
        0.7372773,
        0.7933533,
        0.8433914,
        0.8870108,
        0.9238795,
        0.95371693,
        0.976296,
        0.9914449,
        0.99904823,
        0.99904823,
        0.9914449,
        0.976296,
        0.95371693,
        0.9238795,
        0.8870108,
        0.8433914,
        0.7933533,
        0.7372773,
        0.6755902,
        0.6087614,
        0.53729963,
        0.4617486,
        0.38268343,
        0.3007058,
        0.21643962,
        0.13052619,
        0.043619387,
    ],
    [
        0.043619387,
        0.13052619,
        0.21643962,
        0.3007058,
        0.38268343,
        0.4617486,
        0.53729963,
        0.6087614,
        0.6755902,
        0.7372773,
        0.7933533,
        0.8433914,
        0.8870108,
        0.9238795,
        0.95371693,
        0.976296,
        0.9914449,
        0.99904823,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        0.9914449,
        0.9238795,
        0.7933533,
        0.6087614,
        0.38268343,
        0.13052619,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.13052619, 0.38268343, 0.6087614, 0.7933533, 0.9238795, 0.9914449, 0.9914449, 0.9238795,
        0.7933533, 0.6087614, 0.38268343, 0.13052619, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.13052619,
        0.38268343,
        0.6087614,
        0.7933533,
        0.9238795,
        0.9914449,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        0.99904823,
        0.9914449,
        0.976296,
        0.95371693,
        0.9238795,
        0.8870108,
        0.8433914,
        0.7933533,
        0.7372773,
        0.6755902,
        0.6087614,
        0.53729963,
        0.4617486,
        0.38268343,
        0.3007058,
        0.21643962,
        0.13052619,
        0.043619387,
    ],
];

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
        panic!("Mpeg1 not supported");
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
        *granule = read_granule_side_info(header, &mut reader)?;
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
            data.granules[g].channels[c].count1 = read_huffman(
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

    let (scale_lens, lfs_table) = {
        let sfc = u32::from(channel_info.scalefac_compress);
        match sfc {
            0..=399 => (
                [sfc / 80, (sfc / 16) % 5, (sfc % 16) / 4, sfc & 3],
                &lfs_table[0],
            ),
            _ => panic!("Invalid scalefac_compress"),
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
            panic!("Mixed block type not supported");
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

    if header.channels == Channels::Mono {
        for gr in 0..header.num_granules() {
            let side_info = &side_info.granules[gr].channels[0];
            let main_data = &mut main_data.granules[gr].channels[0];

            requantize(header, side_info, main_data);
            reorder(header, side_info, main_data);
            antialias(side_info, &mut main_data.samples);
            hybrid_synthesis(
                side_info.block_type,
                &mut decoder.store[0],
                &mut main_data.samples,
            );
            frequency_inversion(&mut main_data.samples);
            subband_synthesis(
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



#[allow(clippy::unreadable_literal)]
pub fn antialias(side_info: &GranuleChannelSideInfo, samples: &mut [f32; 576]) {
    const CS: [f32; 8] = [
        0.857493, 0.881742, 0.949629, 0.983315, 0.995518, 0.999161, 0.999899, 0.999993,
    ];
    const CA: [f32; 8] = [
        -0.514496, -0.471732, -0.313377, -0.181913, -0.094574, -0.040966, -0.014199, -0.003700,
    ];

    let sblim = if side_info.block_type == BlockType::Short {
        // No anti-aliasing done for short blocks.
        return;
    } else if side_info.block_type == BlockType::Mixed {
        2
    } else {
        32
    };

    for sb in 1..sblim {
        for i in 0..8 {
            let li = 18 * sb - 1 - i;
            let ui = 18 * sb + i;
            let lb = samples[li] * CS[i] - samples[ui] * CA[i];
            let ub = samples[ui] * CS[i] + samples[li] * CA[i];
            samples[li] = lb;
            samples[ui] = ub;
        }
    }
}

pub(crate) fn hybrid_synthesis(
    block_type: BlockType,
    store: &mut [[f32; 18]; 32],
    samples: &mut [f32; 576],
) {
    for sb in 0..32 {
        let block_type = match block_type {
            BlockType::Long => 0,
            BlockType::Start => 1,
            BlockType::Short => 2,
            BlockType::Mixed => {
                if sb < 2 {
                    0
                } else {
                    2
                }
            }
            BlockType::End => 3,
        };

        let out = imdct_win(block_type, &samples[sb * 18..sb * 18 + 18]);
        for i in 0..18 {
            samples[sb * 18 + i] = out[i] + store[sb][i];
            store[sb][i] = out[i + 18];
        }
    }
}

fn imdct_win(block_type: usize, samples: &[f32]) -> [f32; 36] {
    let mut out = [0f32; 36];
    let imdct_table = &IMDCT_WIN[block_type];
    if block_type == 2 {
        for i in 0..3 {
            for p in 0..12 {
                let mut sum = 0.0;
                for m in 0..6 {
                    sum += samples[i + 3 * m] * COS_N12[m][p];
                }
                out[6 * i + p + 6] = sum * imdct_table[p];
            }
        }
    } else {
        for p in 0..36 {
            let mut sum = 0.0;
            for m in 0..18 {
                sum += samples[m] * COS_N36[m][p];
            }

            out[p] = sum * imdct_table[p];
        }
    }
    out
}

pub fn frequency_inversion(samples: &mut [f32; 576]) {
    for sb in (1..32).step_by(2) {
        for i in (1..18).step_by(2) {
            let n = sb * 18 + i;
            samples[n] = -samples[n];
        }
    }
}

pub fn subband_synthesis(samples: &[f32; 576], v_vec: &mut [f32; 1024], out: &mut [f32]) {
    let mut s_vec = [0f32; 32];
    let mut u_vec = [0f32; 512];

    for ss in 0..18 {
        for i in (64..=1023).rev() {
            v_vec[i] = v_vec[i - 64];
        }

        for i in 0..32 {
            s_vec[i] = samples[i * 18 + ss];
        }

        for (i, row) in SBS_N_WIN.iter().enumerate() {
            let mut sum = 0.0;
            for (j, &sbs_n_win) in row.iter().enumerate() {
                sum += sbs_n_win * s_vec[j];
            }
            v_vec[i] = sum;
        }

        for i in 0..8 {
            for j in 0..32 {
                let i6 = i << 6;
                let i7 = i << 7;

                u_vec[i6 + j] = v_vec[i7 + j];
                u_vec[i6 + j + 32] = v_vec[i7 + j + 96];
            }
        }

        for i in 0..512 {
            u_vec[i] *= SYNTH_DTBL[i];
        }

        for i in 0..32 {
            let mut sum = 0.0;
            for j in 0..16 {
                sum += u_vec[(j << 5) + i];
            }
            out[(32 * ss) + i] = sum;
        }
    }
}


pub fn requantize(
    header: &FrameHeader,
    side_info: &GranuleChannelSideInfo,
    main_data: &mut MainDataChannel,
) {
    if side_info.block_type == BlockType::Short || side_info.block_type == BlockType::Mixed {
        if side_info.block_type == BlockType::Mixed {
            let mut sfb = 0;
            let mut next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].0[sfb + 1];
            for i in 0..36 {
                if i == next_sfb {
                    sfb += 1;
                    next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].0[sfb + 1];
                }

                requantize_long(side_info, i as usize, sfb, main_data);
            }

            sfb = 3;
            next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1] * 3;
            let mut window_len = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1]
                - SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb];

            let mut i = 36;
            while i < main_data.count1 {
                if i == next_sfb {
                    assert!(sfb < 14);
                    sfb += 1;
                    next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1] * 3;
                    window_len = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1]
                        - SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb];
                }

                for win in 0..3 {
                    for _ in 0..window_len {
                        requantize_short(
                            side_info,
                            i as usize,
                            sfb,
                            win,
                            &side_info.subblock_gain[..],
                            main_data,
                        );
                        i += 1;
                    }
                }
            }
        } else {
            // Data only contains short blocks.
            let mut sfb = 0;
            let mut next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1] * 3;
            let mut window_len = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1]
                - SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb];

            let mut i = 0;
            while i < main_data.count1 {
                if i == next_sfb {
                    assert!(sfb < 14);
                    sfb += 1;
                    next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1] * 3;
                    window_len = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb + 1]
                        - SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1[sfb];
                }

                for win in 0..3 {
                    for _ in 0..window_len {
                        requantize_short(
                            side_info,
                            i as usize,
                            sfb,
                            win,
                            &side_info.subblock_gain[..],
                            main_data,
                        );
                        i += 1;
                    }
                }
            }
        }
    } else {
        // Data contains only long blocks.
        let mut sfb = 0;
        let mut next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].0[sfb + 1];

        for i in 0..main_data.count1 {
            if i == next_sfb {
                assert!(sfb < 23);
                sfb += 1;
                next_sfb = SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].0[sfb + 1];
            }

            requantize_long(side_info, i as usize, sfb, main_data);
        }
    }
}

// Requnaitze subband using long blocks.
fn requantize_long(
    side_info: &GranuleChannelSideInfo,
    pos: usize,
    sfb: usize,
    main_data: &mut MainDataChannel,
) {
    const PRE_TAB: [f32; 22] = [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0,
        3.0, 3.0, 2.0, 0.0,
    ];

    assert!(pos < 576);
    let sf_mult = if side_info.scalefac_scale { 1.0 } else { 0.5 };
    let pf_x_pt = if side_info.preflag { PRE_TAB[sfb] } else { 0.0 };
    let tmp1 = f64::powf(
        2.0,
        -sf_mult * (f64::from(main_data.scalefac_l[sfb]) + f64::from(pf_x_pt)),
    );
    let tmp2 = f64::powf(2.0, 0.25 * (f64::from(side_info.global_gain) - 210.0));
    let tmp3 = if main_data.samples[pos] < 0.0 {
        -requantize_pow_43(-main_data.samples[pos])
    } else {
        requantize_pow_43(main_data.samples[pos])
    };

    main_data.samples[pos] = (tmp1 * tmp2 * f64::from(tmp3)) as f32;
}

// Requanitze short block subband.
fn requantize_short(
    side_info: &GranuleChannelSideInfo,
    pos: usize,
    sfb: usize,
    window: usize,
    subblock_gain: &[f32],
    data: &mut MainDataChannel,
) {
    assert!(pos < 576);
    let sf_mult = if side_info.scalefac_scale { 1.0 } else { 0.5 };
    let tmp1 = f64::powf(2.0, -sf_mult * f64::from(data.scalefac_s[sfb][window]));
    let tmp2 = f64::powf(
        2.0,
        0.25 * (f64::from(side_info.global_gain) - 210.0 - 8.0 * f64::from(subblock_gain[window])),
    );
    let tmp3 = if data.samples[pos] < 0.0 {
        -requantize_pow_43(-data.samples[pos])
    } else {
        requantize_pow_43(data.samples[pos])
    };
    data.samples[pos] = (tmp1 * tmp2 * f64::from(tmp3)) as f32;
}

fn requantize_pow_43(sample: f32) -> f32 {
    f32::powf(f32::trunc(sample), 4.0 / 3.0)
}

pub fn reorder(
    header: &FrameHeader,
    side_info: &GranuleChannelSideInfo,
    main_data: &mut MainDataChannel,
) {
    let mut reorder_buffer = [0f32; 576];

    let band_indices = &SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].1;
    if side_info.block_type == BlockType::Short || side_info.block_type == BlockType::Mixed {
        let mut sfb = if side_info.block_type == BlockType::Mixed {
            3
        } else {
            0
        };
        let mut next_sfb = band_indices[sfb + 1] * 3;
        let mut window_len = (band_indices[sfb + 1] - band_indices[sfb]) as usize;
        let mut i = if sfb == 0 { 0 } else { 36 };
        while i < 576 {
            if i == next_sfb {
                for (j, &val) in reorder_buffer[0..3 * window_len].iter().enumerate() {
                    main_data.samples[3 * band_indices[sfb] as usize + j] = val;
                }

                if i >= main_data.count1 {
                    return;
                }

                sfb += 1;
                next_sfb = band_indices[sfb + 1] * 3;
                window_len = (band_indices[sfb + 1] - band_indices[sfb]) as usize;
            }

            for win in 0..3 {
                for j in 0..window_len {
                    reorder_buffer[j * 3 + win] = main_data.samples[i as usize];
                    i += 1;
                }
            }
        }
    }
}

pub fn read_huffman<R: Read>(
    reader: &mut BitReader<R, BigEndian>,
    header: &FrameHeader,
    side_info: &GranuleChannelSideInfo,
    len: u32,
    samples: &mut [f32; 576],
) -> Result<u32, Error> {
    if len == 0 {
        for sample in samples.iter_mut() {
            *sample = 0.0;
        }
        return Ok(0);
    }

    let mut bits_read = 0;
    // ? let bit_pos_end = part_2_start + side_info.part2_3_length - 1;

    let (region1_start, region2_start) =
        if side_info.block_type == BlockType::Short || side_info.block_type == BlockType::Mixed {
            (36, 576)
        } else {
            (
                SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].0
                    [side_info.region0_count as usize + 1],
                SCALE_FACTOR_BAND_INDICES[header.sample_rate_table].0
                    [side_info.region0_count as usize + side_info.region1_count as usize + 2],
            )
        };

    // Read big_values.
    let mut is_pos: usize = 0;
    let is_len = side_info.big_values as usize * 2;
    let mut state: HuffmanState = Default::default();
    while is_pos < is_len {
        let table_num = if is_pos < region1_start as usize {
            side_info.table_select[0]
        } else if is_pos < region2_start as usize {
            side_info.table_select[1]
        } else {
            side_info.table_select[2]
        };

        let huffman_table = &HUFFMAN_TABLES[table_num as usize];
        // TODO(Herschel): Is state an inout parameter or just output?
        bits_read += huffman_decode(reader, huffman_table, &mut state)?;

        samples[is_pos] = state.x as f32;
        is_pos += 1;
        samples[is_pos] = state.y as f32;
        is_pos += 1;
    }

    // Read small values until is_pos is 576
    let table_num = if side_info.count1table_select { 33 } else { 32 };
    let huffman_table = &HUFFMAN_TABLES[table_num];
    is_pos = is_len;
    while is_pos <= 572 && bits_read < len as usize {
        bits_read += huffman_decode(reader, huffman_table, &mut state)?;
        samples[is_pos] = state.v as f32;
        is_pos += 1;
        samples[is_pos] = state.w as f32;
        is_pos += 1;
        samples[is_pos] = state.x as f32;
        is_pos += 1;
        samples[is_pos] = state.y as f32;
        is_pos += 1;
    }

    if bits_read < len as usize {
        reader.skip(len - bits_read as u32)?;
    } else if bits_read > len as usize {
        is_pos -= 4;
    }

    for sample in &mut samples[is_pos..576] {
        *sample = 0.0;
    }

    Ok(is_pos as u32)
}

#[derive(Debug, Default)]
struct HuffmanState {
    x: i32,
    y: i32,
    v: i32,
    w: i32,
}
fn huffman_decode<R: Read>(
    reader: &mut BitReader<R, BigEndian>,
    huffman_table: &HuffmanTable,
    state: &mut HuffmanState,
) -> Result<usize, Error> {
    let mut point = 0;
    let mut bits_left = 32;
    let mut bits_read = 0;
    if !huffman_table.data.is_empty() {
        loop {
            if huffman_table.data[point] & 0xff00 == 0 {
                state.x = ((huffman_table.data[point] >> 4) & 0xf).into();
                state.y = (huffman_table.data[point] & 0xf).into();
                break;
            }

            bits_read += 1;
            if reader.read_bit()? {
                while (huffman_table.data[point] & 0xff) >= 250 {
                    point += (huffman_table.data[point] & 0xff) as usize;
                }
                point += (huffman_table.data[point] & 0xff) as usize;
            } else {
                while (huffman_table.data[point] >> 8) >= 250 {
                    point += (huffman_table.data[point] >> 8) as usize;
                }
                point += (huffman_table.data[point] >> 8) as usize;
            }

            bits_left -= 1;
            if bits_left <= 0 || point >= huffman_table.data.len() {
                break;
            }
        }

        if huffman_table.quads {
            state.v = (state.y >> 3) & 1;
            state.w = (state.y >> 2) & 1;
            state.x = (state.y >> 1) & 1;
            state.y &= 1;

            if state.v > 0 {
                bits_read += 1;
                if reader.read_bit()? {
                    state.v = -state.v;
                }
            }
            if state.w > 0 {
                bits_read += 1;
                if reader.read_bit()? {
                    state.w = -state.w;
                }
            }
            if state.x > 0 {
                bits_read += 1;
                if reader.read_bit()? {
                    state.x = -state.x;
                }
            }
            if state.y > 0 {
                bits_read += 1;
                if reader.read_bit()? {
                    state.y = -state.y;
                }
            }
        } else {
            if huffman_table.linbits > 0 && state.x == 15 {
                bits_read += huffman_table.linbits;
                // TODO(Herschel): u32?
                state.x += reader.read::<u32>(huffman_table.linbits as u32)? as i32;
            }

            if state.x > 0 {
                bits_read += 1;
                if reader.read_bit()? {
                    state.x = -state.x;
                }
            }

            if huffman_table.linbits > 0 && state.y == 15 {
                bits_read += huffman_table.linbits;
                state.y += reader.read::<u32>(huffman_table.linbits as u32)? as i32;
            }

            if state.y > 0 {
                bits_read += 1;
                if reader.read_bit()? {
                    state.y = -state.y;
                }
            }
        }
    } else {
        *state = Default::default();
    }
    Ok(bits_read)
}

/// Convenience method to decode an MP3.
/// Returns the first frame header found in the MP3, and an `Iterator` that
/// yields MP3 `Sample`s`.
///
/// Each `Sample` represents one left and right sample at the sample rate of
/// the MP3. Any invalid data is ignored. The iterator will provide `Sample`s
/// until there is no more data, or an error occurs.
///
/// If you need to handle errors or changes in the format mid-stream, use
/// `Mp3Decoder` driectly.
pub fn read_mp3<R: Read>(
    reader: R,
) -> Result<(FrameHeader, impl Iterator<Item = (f32, f32)>), Error> {
    let mut decoder = Mp3Decoder::new(reader);
    let mut frame = decoder.next_frame()?;
    let header = frame.header.clone();
    let mut i = 0;
    let iter = std::iter::from_fn(move || {
        if i >= frame.num_samples {
            i = 0;
            frame = if let Ok(frame) = decoder.next_frame() {
                frame
            } else {
                return None;
            }
        }
        let sample = (frame.samples[0][i], frame.samples[1][i]);
        i += 1;
        Some(sample)
    });
    Ok((header, iter))
}

/// Decodes MP3 streams.
pub struct Mp3Decoder<R: Read> {
    reader: R,
    state: DecoderState,
}

impl<R: Read> Mp3Decoder<R> {
    /// Creates a new `MP3Decoder` from the given reader.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            state: DecoderState::new(),
        }
    }

    /// Gets a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        &self.reader
    }

    /// Gets a mutable reference to the underlying reader.
    ///
    /// It is inadvisable to directly read from the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.reader
    }

    /// Unwraps the `Mp3Decoder`, returning the underlying reader.
    pub fn into_inner(self) -> R {
        self.reader
    }

    /// Returns an `Iterator` that yields MP3 `Frame`s.
    ///
    /// Each `Frame` contains header information and the decoded samples.
    /// Any invalid data is skipped. The iterator will provide `Frame`s until
    /// there is no more valid MP3 data or an error occurs.
    ///
    /// If you wish to inspect any errors, Use `next_frame` instead.
    pub fn frames(mut self) -> impl Iterator<Item = Frame> {
        std::iter::from_fn(move || self.next_frame().ok())
    }

    /// Decodes the next MP3 `Frame` in the stream.
    ///
    /// Data is read until a valid `Frame` is found. Invalid data is skipped.
    /// Other errors are returned.
    pub fn next_frame(&mut self) -> Result<Frame, Error> {
        let header;
        loop {
            match read_frame_header(&mut self.reader) {
                Ok(frame_header) => {
                    header = frame_header;
                    break;
                }
                Err(Error::Mp3Error(Mp3Error::InvalidData(_))) => (),
                Err(e) => return Err(e),
            }
        }

        let (num_samples, samples) =
            process_frame(&mut self.state, &mut self.reader, &header)?;

        Ok(Frame {
            header,
            samples,
            num_samples,
        })
    }
}

/// A frame of MP3 data.
///
/// Each frame contains a header describing the format of the data, and the decoded
/// samples. An MP3 frame contains either 576 or 1152 samples (depending on the
/// format).
pub struct Frame {
    /// The header of this MP3 frame.
    pub header: FrameHeader,

    /// The decoded MP3 samples for the left and right channels.
    /// Each sample is in the range of [-1.0, 1.0].
    /// Only the first `num_samples` entries will contain valid data.
    /// For mono streams, the data will be duplicated to the left and right
    /// channels.
    pub samples: [[f32; 1152]; 2],

    /// The number of samples in the `samples` array.
    /// This will be either 576 or 1152 samples depending on the
    /// format of the MP3.
    pub num_samples: usize,
}


fn main() {
    let fpath = std::env::args().nth(1).expect("Please provide a file path");

    let data = std::fs::read(fpath).expect("Could not open file");
    let (_header, samples) = read_mp3(&data[..]).expect("Invalid MP3");

    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: 8000,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer =
        hound::WavWriter::create("output.wav", spec).expect("Failed to create WAV writer");

    for (left, _right) in samples {
        let left_i16 = (left * i16::MAX as f32) as i16;
        writer
            .write_sample(left_i16)
            .expect("Failed to write left sample");
    }
}
