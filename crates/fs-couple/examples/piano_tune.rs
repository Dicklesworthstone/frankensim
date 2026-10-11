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

const USAGE: &str = "piano_tune [--preset steinway-d | --scale strings.csv]
    [--concert-pitch 430..450] --output tuned.csv
Default source: the explicitly estimated demonstration scale. Default A4: 440 Hz.
Preserves length, mass, EI, unison spread and contact geometry; changes tension.
C4..B4 remains equal tempered. Outside it, lower second and upper first partials
match. This is a declared 2:1 octave policy, not a measured tuning curve.
Existing output files refuse. Render the result with grand_piano --scale tuned.csv
without --concert-pitch, which would replace the stretch with equal temperament.";

#[derive(Debug)]
struct Options { scale: Option<String>, preset: bool, concert_a: f64, output: String }
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut seen = std::collections::BTreeSet::new();
        let (mut scale, mut preset, mut concert_a, mut output) = (None, false, 440.0, None);
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
                _ => return Err(format!("unknown option {flag}")),
            }
        }
        if scale.is_some() && preset { return Err("choose --scale or --preset, not both".into()); }
        if !concert_a.is_finite() || !(430.0..=450.0).contains(&concert_a) {
            return Err("concert pitch must be finite and in 430..450 Hz".into());
        }
        let output = output.filter(|s| !s.is_empty()).ok_or("--output is required")?;
        Ok(Self { scale, preset, concert_a, output })
    }
    fn prepare(&self) -> Result<Vec<geometry::Course>, String> {
        let mut courses = match &self.scale {
            Some(path) => geometry::read_scale(&std::fs::read_to_string(path)
                .map_err(|e| format!("read {path}: {e}"))?)?,
            None if self.preset => steinway_scale::courses()?,
            None => geometry::demonstration_scale()?,
        };
        stretch_tuning::tune(&mut courses, self.concert_a)?;
        Ok(courses)
    }
}

fn run(options: &Options) -> Result<(), String> {
    // Complete all physical validation before creating the destination. Exclusive
    // creation also refuses input/output aliases, existing files and symlinks.
    let courses = options.prepare()?;
    let text = geometry::write_scale(&courses);
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true)
        .open(&options.output).map_err(|e| format!("create {}: {e}", options.output))?;
    file.write_all(text.as_bytes()).and_then(|()| file.sync_all())
        .map_err(|e| format!("write {}: {e}", options.output))?;
    println!("Prepared {} courses at A4={} Hz with physical 2:1 octave stretch: {}",
        courses.len(), options.concert_a, options.output);
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
            let options = Options { preset, scale: None, concert_a: 442.0, output: String::new() };
            let courses = options.prepare().unwrap();
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
            output: path.to_str().unwrap().to_owned() };
        assert!(run(&options).is_err());
        assert_eq!(std::fs::read_to_string(path).unwrap(), "retain me");
    }
}
