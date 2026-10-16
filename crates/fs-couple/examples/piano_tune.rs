//! Prepare physical octave-stretched string tensions for the existing renderer.
//! No sound generation or mechanical model is duplicated here.
#[allow(dead_code)]
#[path = "grand_piano/geometry.rs"]
mod geometry;
#[allow(dead_code)]
#[path = "grand_piano/steinway_scale.rs"]
mod steinway_scale;
#[path = "grand_piano/stretch_tuning.rs"]
mod stretch_tuning;

use std::io::Write;
use stretch_tuning::partials::{OctaveWeights, Profile, OctaveReport};

const USAGE: &str = "piano_tune [--preset steinway-d | --scale strings.csv]
    [--concert-pitch 430..450] [--bass-octave-weights w21,w42,w63,w84]
    [--treble-octave-weights w21,w42,w63,w84] [--tuning-band max_partial_hz]
    --output tuned.csv
Default source: the explicitly estimated demonstration scale. Default A4: 440 Hz.
Preserves length, mass, EI, unison spread and contact geometry; changes tension.
C4..B4 remains equal tempered. Defaults match lower second and upper first partials.
Explicit register weights blend 2:1, 4:2, 6:3 and 8:4 partial matches by minimizing
weighted squared beat/n in physical Hz. Weights must be nonnegative, not all zero.
The optional tuning band excludes neighbour target partials above its ceiling;
an octave with no admitted weighted target refuses. This is not a measured curve.
Existing output files refuse. Render the result with grand_piano --scale tuned.csv
without --concert-pitch or --stretch-tuning, which would replace this preparation.";

#[derive(Debug)]
struct Options {
    scale: Option<String>, preset: bool, concert_a: f64, output: String,
    tuning: Profile,
}
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut seen = std::collections::BTreeSet::new();
        let (mut scale, mut preset, mut concert_a, mut output) = (None, false, 440.0, None);
        let mut tuning = Profile::default();
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
            let value = args.next().filter(|v| !v.starts_with("--"))
                .ok_or_else(|| format!("missing value for {flag}"))?;
            match flag.as_str() {
                "--scale" => scale = Some(value.clone()),
                "--preset" if value == "steinway-d" => preset = true,
                "--preset" => return Err("only --preset steinway-d is supported".into()),
                "--concert-pitch" => concert_a = value.parse::<f64>()
                    .map_err(|_| "invalid concert pitch")?,
                "--output" => output = Some(value.clone()),
                "--bass-octave-weights" => tuning.bass = OctaveWeights::parse(value)?,
                "--treble-octave-weights" => tuning.treble = OctaveWeights::parse(value)?,
                "--tuning-band" => tuning.max_partial_hz = Some(value.parse::<f64>()
                    .map_err(|_| "invalid tuning band")?),
                _ => return Err(format!("unknown option {flag}")),
            }
        }
        if scale.is_some() && preset { return Err("choose --scale or --preset, not both".into()); }
        if !concert_a.is_finite() || !(430.0..=450.0).contains(&concert_a) {
            return Err("concert pitch must be finite and in 430..450 Hz".into());
        }
        if tuning.max_partial_hz.is_some_and(|hz| !hz.is_finite() || hz <= 0.0) {
            return Err("tuning band must be finite and positive".into());
        }
        let output = output.filter(|s| !s.is_empty()).ok_or("--output is required")?;
        Ok(Self { scale, preset, concert_a, output, tuning })
    }
    fn prepare(&self) -> Result<(Vec<geometry::Course>, Vec<OctaveReport>), String> {
        let mut courses = match &self.scale {
            Some(path) => geometry::read_scale(&std::fs::read_to_string(path)
                .map_err(|e| format!("read {path}: {e}"))?)?,
            None if self.preset => steinway_scale::courses()?,
            None => geometry::demonstration_scale()?,
        };
        let report = stretch_tuning::partials::tune(&mut courses, self.concert_a, self.tuning)?;
        Ok((courses, report))
    }
}

fn run(options: &Options) -> Result<(), String> {
    // Complete all physical validation before creating the destination. Exclusive
    // creation also refuses input/output aliases, existing files and symlinks.
    let (courses, report) = options.prepare()?;
    let text = geometry::write_scale(&courses);
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true)
        .open(&options.output).map_err(|e| format!("create {}: {e}", options.output))?;
    file.write_all(text.as_bytes()).and_then(|()| file.sync_all())
        .map_err(|e| format!("write {}: {e}", options.output))?;
    println!("Prepared {} courses at A4={} Hz; bass weights {:?}, treble weights {:?}: {}",
        courses.len(), options.concert_a, options.tuning.bass.values(),
        options.tuning.treble.values(), options.output);
    println!("Checked {} octave pairs against the declared physical partial objective; target band {:?} Hz",
        report.len(), options.tuning.max_partial_hz);
    if let Some(worst) = report.iter().max_by(|a,b| a.weighted_rms_hz.total_cmp(&b.weighted_rms_hz)) {
        println!("Largest weighted RMS beat/n: {:.9e} Hz at MIDI {}:{}; signed beats {:?} Hz",
            worst.weighted_rms_hz, worst.lower_midi, worst.lower_midi+12, worst.signed_beats_hz);
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.as_slice() == ["--help"] || args.as_slice() == ["-h"] {
        println!("{USAGE}");
        return;
    }
    if let Err(error) = Options::parse(&args).and_then(|options| run(&options)) {
        eprintln!("piano_tune: {error}\n{USAGE}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parse(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }
    #[test]
    fn source_and_reference_selection_are_explicit() {
        let demo = parse(&["--output", "demo.csv"]).unwrap();
        assert!(!demo.preset && demo.scale.is_none());
        assert_eq!(demo.concert_a, 440.0);
        let preset = parse(&["--preset", "steinway-d", "--concert-pitch", "442",
            "--output", "tuned.csv"]).unwrap();
        assert!(preset.preset);
        assert_eq!(preset.concert_a, 442.0);
        assert_eq!(parse(&["--scale", "measured.csv", "--output", "out.csv"])
            .unwrap().scale.as_deref(), Some("measured.csv"));
    }
    #[test]
    fn ambiguous_invalid_and_incomplete_options_refuse() {
        for args in [vec![], vec!["--output"], vec!["--output", ""],
            vec!["--output", "a", "--output", "b"],
            vec!["--output", "a", "--preset", "other"],
            vec!["--output", "a", "--preset", "steinway-d", "--scale", "b"],
            vec!["--output", "a", "--concert-pitch", "NaN"],
            vec!["--output", "a", "--concert-pitch", "451"],
            vec!["--output", "a", "--unknown", "b"]] {
            assert!(parse(&args).is_err(), "{args:?}");
        }
    }
    #[test]
    fn both_scales_match_octaves_and_roundtrip_as_renderer_input() {
        for preset in [false, true] {
            let options = Options { preset, scale: None, concert_a: 442.0, output: String::new(), tuning: Profile::default() };
            let (courses, _) = options.prepare().unwrap();
            assert_eq!(geometry::read_scale(&geometry::write_scale(&courses)).unwrap(), courses);
            assert!((courses[48].partial_hz(1, courses[48].tension_n) / 442.0 - 1.0).abs() < 1e-11);
            for lower in 0..76 {
                let a = courses[lower]; let b = courses[lower + 12];
                assert!((a.partial_hz(2, a.tension_n) / b.partial_hz(1, b.tension_n) - 1.0).abs() < 1e-11);
            }
        }
    }
    #[test]
    fn existing_output_is_not_truncated_or_overwritten() {
        let dir = std::env::temp_dir().join(format!("piano-tune-{}-{}", std::process::id(),
            std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        std::fs::create_dir(&dir).unwrap();
        let path = dir.join("retained.csv");
        std::fs::write(&path, "retain me").unwrap();
        let options = Options { scale: None, preset: false, concert_a: 440.0,
            output: path.to_str().unwrap().to_owned(), tuning: Profile::default() };
        assert!(run(&options).is_err());
        assert_eq!(std::fs::read_to_string(path).unwrap(), "retain me");
    }
}

#[cfg(test)]
#[path = "grand_piano/partial_tuning_cli_tests.rs"]
mod partial_tuning_cli_tests;
