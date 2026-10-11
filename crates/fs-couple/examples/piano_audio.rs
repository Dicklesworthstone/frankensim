//! Export the piano's unquantized receiver trace without rerunning mechanics.
use fs_couple::pcm_wav::high_resolution::{Encoding, encode_pressure_wav};
use std::io::{BufRead, Read, Write};

const MONO: &str = "sample,time_s,pressure_left_pa";
const STEREO: &str = "sample,time_s,pressure_left_pa,pressure_right_pa";
const MAX_VALUES: usize = 32 * 1024 * 1024;
const MAX_LINE_BYTES: u64 = 4096;
const USAGE: &str = "piano_audio --input receivers.csv --output piano.wav
    --sample-rate 8000..192000 --full-scale-pa positive-Pa --format pcm24|float32
Input is grand_piano --receiver-pressure-csv, NOT a previously quantized WAV.
The original clock, channel order and one common pressure scale are preserved.
PCM24 refuses clipping. Float32 retains finite excursions beyond +/-1, but a
player may clip them; explicitly attenuate before fixed-point playback.
No normalization, dither, resampling, EQ or reconstructed overtones. All input
is validated before exclusive output creation. Maximum 33554432 channel values.";

#[derive(Debug)]
struct Options { input: String, output: String, rate: u32, full_scale: f64, encoding: Encoding }
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut seen = std::collections::BTreeSet::new();
        let (mut input, mut output, mut rate, mut full_scale, mut encoding) = (None, None, None, None, None);
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
            let value = args.next().filter(|v| !v.is_empty() && !v.starts_with("--"))
                .ok_or_else(|| format!("missing value for {flag}"))?;
            match flag.as_str() {
                "--input" => input = Some(value.clone()),
                "--output" => output = Some(value.clone()),
                "--sample-rate" => rate = Some(value.parse::<u32>().map_err(|_| "invalid sample rate")?),
                "--full-scale-pa" => full_scale = Some(value.parse::<f64>().map_err(|_| "invalid full scale")?),
                "--format" => encoding = Some(match value.as_str() {
                    "pcm24" => Encoding::Pcm24, "float32" => Encoding::Float32,
                    _ => return Err("format must be pcm24 or float32".into()),
                }),
                _ => return Err(format!("unknown option {flag}")),
            }
        }
        let rate = rate.ok_or("--sample-rate is required; the clock is never inferred")?;
        let full_scale = full_scale.ok_or("--full-scale-pa is required; no peak normalization")?;
        if !(8000..=192000).contains(&rate) || !full_scale.is_finite() || full_scale <= 0.0 {
            return Err("invalid clock or finite positive pressure full scale".into());
        }
        Ok(Self { input: input.ok_or("--input is required")?, output: output.ok_or("--output is required")?,
            rate, full_scale, encoding: encoding.ok_or("--format is required")? })
    }
}

fn row(reader: &mut impl BufRead, line: &mut String) -> Result<usize, String> {
    line.clear();
    let size = reader.take(MAX_LINE_BYTES + 1).read_line(line).map_err(|e| e.to_string())?;
    if size as u64 > MAX_LINE_BYTES { return Err("receiver CSV row exceeds 4096 bytes".into()); }
    Ok(size)
}

fn read_pressure(mut reader: impl BufRead, rate: u32) -> Result<(Vec<f64>, u16), String> {
    if !(8000..=192000).contains(&rate) { return Err("invalid receiver clock".into()); }
    let mut line = String::new();
    row(&mut reader, &mut line)?;
    let channels = match line.trim() {
        MONO => 1u16, STEREO => 2u16,
        _ => return Err("expected the unquantized mono/stereo receiver-pressure CSV header".into()),
    };
    let mut pressure = Vec::new();
    let mut frame = 0usize;
    while row(&mut reader, &mut line)? != 0 {
        let mut fields = line.trim().split(',').map(str::trim);
        let index = fields.next().and_then(|v| v.parse::<usize>().ok());
        if index != Some(frame) { return Err(format!("receiver frame {frame}: noncontiguous sample index")); }
        let time = fields.next().and_then(|v| v.parse::<f64>().ok())
            .ok_or_else(|| format!("receiver frame {frame}: invalid time"))?;
        let expected = frame as f64 / f64::from(rate);
        let tolerance = 16.0 * f64::EPSILON * expected.abs().max(1.0 / f64::from(rate));
        if !time.is_finite() || (time - expected).abs() > tolerance {
            return Err(format!("receiver frame {frame}: timestamp does not match --sample-rate"));
        }
        if pressure.len() > MAX_VALUES - usize::from(channels) {
            return Err("receiver history exceeds the 33554432 channel-value memory budget".into());
        }
        for _ in 0..channels {
            let value = fields.next().and_then(|v| v.parse::<f64>().ok())
                .filter(|v| v.is_finite()).ok_or_else(|| format!("receiver frame {frame}: invalid pressure"))?;
            pressure.push(value);
        }
        if fields.next().is_some() { return Err(format!("receiver frame {frame}: extra channel/field")); }
        frame += 1;
    }
    if pressure.is_empty() { return Err("receiver trace has no frames".into()); }
    Ok((pressure, channels))
}

fn run(options: &Options) -> Result<(), String> {
    let input = std::fs::File::open(&options.input).map_err(|e| format!("{}: {e}", options.input))?;
    let (pressure, channels) = read_pressure(std::io::BufReader::new(input), options.rate)?;
    let wav = encode_pressure_wav(&pressure, options.rate, channels, options.full_scale, options.encoding)
        .map_err(|e| e.to_string())?;
    // No destination is opened before the entire source and encoding pass.
    // create_new also refuses existing files, aliases and dangling symlinks.
    let mut output = std::fs::OpenOptions::new().write(true).create_new(true).open(&options.output)
        .map_err(|e| format!("{}: {e}", options.output))?;
    output.write_all(&wav).and_then(|()| output.sync_all()).map_err(|e| e.to_string())?;
    println!("{} frames, {channels} channels, {} Hz, {:?}, {} Pa per digital full scale: {}",
        pressure.len() / usize::from(channels), options.rate, options.encoding, options.full_scale, options.output);
    if options.encoding == Encoding::Float32 && pressure.iter().any(|p| p.abs() > options.full_scale) {
        println!("Float samples exceed +/-1 without clipping; attenuate explicitly before fixed-point playback.");
    }
    Ok(())
}
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.as_slice() == ["--help"] || args.as_slice() == ["-h"] { println!("{USAGE}"); return; }
    if let Err(error) = Options::parse(&args).and_then(|options| run(&options)) {
        eprintln!("piano_audio: {error}\n{USAGE}"); std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn source(stereo: bool) -> String {
        let mut csv = format!("{}\n", if stereo { STEREO } else { MONO });
        for frame in 0..3 {
            csv.push_str(&format!("{frame},{:.17e},1e-6{}\n", frame as f64 / 48_000.0,
                if stereo { ",-2e-6" } else { "" }));
        }
        csv
    }
    #[test]
    fn original_receiver_rows_preserve_quiet_pressure_and_stereo_order_in_both_formats() {
        for stereo in [false, true] {
            let (p, channels) = read_pressure(source(stereo).as_bytes(), 48_000).unwrap();
            assert_eq!(channels, if stereo { 2 } else { 1 });
            for frame in p.chunks_exact(usize::from(channels)) {
                assert_eq!(frame[0], 1e-6); if stereo { assert_eq!(frame[1], -2e-6); }
            }
            for encoding in [Encoding::Pcm24, Encoding::Float32] {
                let wav = encode_pressure_wav(&p, 48_000, channels, 1.0, encoding).unwrap();
                let start = if encoding == Encoding::Pcm24 { 68 } else { 80 };
                assert!(wav[start..].iter().any(|b| *b != 0));
            }
        }
    }
    #[test]
    fn altered_clocks_gaps_channels_and_nonfinite_pressure_refuse() {
        assert!(read_pressure(source(false).as_bytes(), 44_100).is_err());
        for bad in [String::new(), format!("{MONO}\n"), source(false).replace("1,", "5,"),
            source(false).replace("1e-6", "NaN"), source(false).replace("1e-6", "1e-6,0"),
            source(true).replace(",-2e-6", ""), format!("{MONO}\n0,1,1e-6\n"),
            format!("{MONO}\n{}", "x".repeat(4097))] {
            assert!(read_pressure(bad.as_bytes(), 48_000).is_err(), "accepted malformed input");
        }
    }
    fn parse(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }
    #[test]
    fn explicit_options_refuse_defaults_duplicates_and_unknown_formats() {
        let args = ["--input", "p.csv", "--output", "p.wav", "--sample-rate", "48000",
            "--full-scale-pa", "2", "--format", "pcm24"];
        let options = parse(&args).unwrap(); assert_eq!(options.encoding, Encoding::Pcm24);
        for index in [0, 2, 4, 6, 8] {
            let mut missing = args.to_vec(); missing.drain(index..index+2);
            assert!(parse(&missing).is_err());
        }
        for (index, value) in [(5,"0"), (7,"NaN"), (9,"pcm16")] {
            let mut bad = args; bad[index] = value; assert!(parse(&bad).is_err());
        }
        let mut duplicate = args.to_vec(); duplicate.extend(["--format", "float32"]);
        assert!(parse(&duplicate).is_err());
    }
    #[test]
    fn fresh_output_runs_once_and_rejects_overwrite_and_invalid_input_before_creation() {
        let directory = std::env::temp_dir().join(format!("piano-audio-{}-{}", std::process::id(),
            std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        std::fs::create_dir(&directory).unwrap();
        let input = directory.join("input.csv"); let output = directory.join("output.wav");
        std::fs::write(&input, source(true)).unwrap();
        let mut options = Options { input: input.to_str().unwrap().into(), output: output.to_str().unwrap().into(),
            rate: 48_000, full_scale: 1.0, encoding: Encoding::Float32 };
        run(&options).unwrap(); let bytes = std::fs::read(&output).unwrap();
        assert!(run(&options).is_err()); assert_eq!(bytes, std::fs::read(&output).unwrap());
        let fresh = directory.join("invalid.wav"); options.output = fresh.to_str().unwrap().into();
        options.rate = 44_100; assert!(run(&options).is_err()); assert!(!fresh.exists());
    }
}
