use super::*;

fn parse(args: &[&str]) -> Result<Options, String> {
    Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
}

#[test]
fn explicit_bass_and_treble_policies_reach_the_actual_scale() {
    for source in [vec![], vec!["--preset", "steinway-d"]] {
        let mut args = source;
        args.extend(["--bass-octave-weights", "0,1,1,0", "--treble-octave-weights", "1,0,0,0",
            "--concert-pitch", "442", "--tuning-band", "20000", "--output", "tuned.csv"]);
        let options = parse(&args).unwrap();
        let (courses, report) = options.prepare().unwrap();
        assert_eq!(courses.len(), 88); assert_eq!(report.len(), 76);
        assert_eq!(options.tuning.bass.values(), [0.,1.,1.,0.]);
        assert_eq!(options.tuning.treble, OctaveWeights::TWO_ONE);
        assert_eq!(geometry::read_scale(&geometry::write_scale(&courses)).unwrap(), courses);
        for row in report {
            if row.lower_midi < 60 {
                assert!(row.signed_beats_hz[0].is_none());
                assert!(row.signed_beats_hz[1].is_some() && row.signed_beats_hz[2].is_some());
            } else {
                assert!(row.signed_beats_hz[0].unwrap().abs() < 1e-8);
                assert!(row.signed_beats_hz[1..].iter().all(Option::is_none));
            }
        }
        assert!((courses[48].partial_hz(1,courses[48].tension_n)-442.).abs() < 1e-9);
    }
}

#[test]
fn ambiguous_nonphysical_and_empty_policy_options_refuse() {
    for tail in [vec!["--bass-octave-weights", "0,0,0,0"],
        vec!["--treble-octave-weights", "1,2,3"],
        vec!["--bass-octave-weights", "1,-1,0,0"],
        vec!["--treble-octave-weights", "NaN,0,1,0"],
        vec!["--tuning-band", "0"], vec!["--tuning-band", "NaN"],
        vec!["--tuning-band", "inf"],
        vec!["--bass-octave-weights", "1,0,0,0", "--bass-octave-weights", "0,1,0,0"]] {
        let mut args = vec!["--output", "tuned.csv"]; args.extend(tail);
        assert!(parse(&args).is_err(), "{args:?}");
    }
}

#[test]
fn default_cli_preparation_preserves_the_original_tensions_exactly() {
    let options = parse(&["--output", "tuned.csv"]).unwrap();
    let (actual, _) = options.prepare().unwrap();
    let mut expected = geometry::demonstration_scale().unwrap();
    stretch_tuning::tune(&mut expected, 440.).unwrap();
    assert_eq!(actual, expected);
}

#[test]
fn failed_tuning_creates_no_output_and_success_exports_the_prepared_physics() {
    let dir = std::env::temp_dir().join(format!("partial-tuning-{}-{}", std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    std::fs::create_dir(&dir).unwrap();
    let input = dir.join("source.csv"); let output = dir.join("tuned.csv");
    let source = geometry::demonstration_scale().unwrap();
    std::fs::write(&input, geometry::write_scale(&source)).unwrap();
    let mut options = Options { scale: Some(input.to_str().unwrap().into()), preset: false,
        concert_a: 440., output: output.to_str().unwrap().into(),
        tuning: Profile { max_partial_hz: Some(1.), ..Profile::default() } };
    assert!(run(&options).is_err()); assert!(!output.exists());
    options.tuning = Profile { bass: OctaveWeights::new([0.,1.,1.,0.]).unwrap(), ..Profile::default() };
    let (expected, _) = options.prepare().unwrap();
    run(&options).unwrap();
    let actual = geometry::read_scale(&std::fs::read_to_string(&output).unwrap()).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(std::fs::read_to_string(&input).unwrap(), geometry::write_scale(&source));
    assert!(run(&options).is_err()); // Exclusive creation still protects existing files.
    assert_eq!(geometry::read_scale(&std::fs::read_to_string(&output).unwrap()).unwrap(), expected);
}
