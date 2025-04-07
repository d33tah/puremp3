fn main() {

    let fpath = std::env::args().nth(1).expect("Please provide a file path");

    let data = std::fs::read(fpath).expect("Could not open file");
    let (_header, samples) = puremp3::read_mp3(&data[..]).expect("Invalid MP3");


    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: 8000,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create("output.wav", spec).expect("Failed to create WAV writer");


    for (left, _right) in samples {
        let left_i16 = (left * i16::MAX as f32) as i16;
        writer.write_sample(left_i16).expect("Failed to write left sample");
    }
}

