//! Authored geometry and mechanics for a new stress-constrained study.
use super::*;
use fs_cutfem::CutSdf;
use fs_marquee::level_set_csv::read_field;
use std::path::PathBuf;

#[derive(Default)]
pub(super) struct Input {
    field: Option<PathBuf>,
    load: Option<f64>,
    band: Option<f64>,
    youngs: Option<f64>,
    poisson: Option<f64>,
}

/// Strip only the authored-problem flags. The existing checkpoint parser still
/// rejects unknown options; resume never calls this parser or admits overrides.
pub(super) fn options(args: &[String]) -> Result<(Vec<String>, Input), Box<dyn Error>> {
    let mut input = Input::default();
    let mut rest = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        if matches!(flag, "--initial-field" | "--load" | "--load-band" | "--youngs" | "--poisson") {
            let value = args.get(i + 1).ok_or_else(|| format!("{flag} requires a value"))?;
            if flag == "--initial-field" {
                if input.field.is_some() || value.is_empty() || value.starts_with("--") {
                    return Err("--initial-field requires one nonempty path and may not be repeated".into());
                }
                input.field = Some(PathBuf::from(value));
            } else {
                let number: f64 = value.parse()?;
                let valid = number.is_finite() && match flag {
                    "--load" | "--youngs" => number > 0.0,
                    "--load-band" => number > 0.0 && number <= 0.5,
                    _ => number > -1.0 && number < 0.5,
                };
                if !valid { return Err(format!("invalid {flag} value").into()); }
                let slot = match flag {
                    "--load" => &mut input.load,
                    "--load-band" => &mut input.band,
                    "--youngs" => &mut input.youngs,
                    _ => &mut input.poisson,
                };
                if slot.replace(number).is_some() { return Err(format!("repeated {flag}").into()); }
            }
            i += 2;
        } else {
            rest.push(args[i].clone());
            i += 1;
        }
    }
    Ok((rest, input))
}

impl Input {
    pub(super) fn prepare(&self, settings: &mut OptimizeSettings) -> Result<(GridSdf, Cantilever), Box<dyn Error>> {
        if let Some(youngs) = self.youngs { settings.youngs = youngs; }
        if let Some(poisson) = self.poisson { settings.poisson = poisson; }
        let n = 1usize << settings.level;
        let geometry = match &self.field {
            Some(path) => read_field(path, n)?,
            None => GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.35),
        };
        let fixture = Cantilever {
            load: self.load.unwrap_or(1.0), band: self.band.unwrap_or(0.125),
        };
        // Canonical checkpoint admission checks the actual material law and
        // optimizer settings before projection or a PDE solve. No rival law.
        let _ = fs_topols::OptimizeCheckpoint::new(geometry.clone(), fixture, *settings)?;
        Ok((geometry, fixture))
    }
}

/// A fixed-load study cannot erase a load by deleting the material beneath it.
/// The two boundary traces are frozen after this check. The bilinear enclosure
/// covers the whole band, including gaps missed by endpoint/midpoint probes.
pub(super) fn require_load_support(geometry: &GridSdf, fixture: Cantilever) -> Result<(), Box<dyn Error>> {
    if !(fixture.band.is_finite() && fixture.band > 0.0 && fixture.band <= 0.5) {
        return Err("load band half-width must be in (0,0.5]".into());
    }
    let enclosure = geometry.enclose([1.0, 0.5 - fixture.band], [1.0, 0.5 + fixture.band]);
    if !(enclosure.lo().is_finite() && enclosure.hi().is_finite() && enclosure.hi() < 0.0) {
        return Err("the complete declared load band must lie inside material; absent, cut or uncertain load support cannot be discarded".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn authored_flags_are_checked_and_do_not_consume_checkpoint_options() {
        let args = ["out", "1e12", "--load", "2", "--youngs", "3", "--poisson", "0.2",
            "--load-band", "0.125", "--checkpoint", "--pause-after", "1"].map(str::to_string);
        let (remaining, input) = options(&args).unwrap();
        assert_eq!(remaining, ["out", "1e12", "--checkpoint", "--pause-after", "1"]);
        assert_eq!((input.load, input.youngs, input.poisson, input.band),
            (Some(2.0), Some(3.0), Some(0.2), Some(0.125)));
        for bad in [vec!["--load"], vec!["--youngs", "NaN"], vec!["--load", "0"],
            vec!["--load-band", "0"], vec!["--load-band", "0.6"],
            vec!["--poisson", "0.5"], vec!["--youngs", "2", "--youngs", "3"],
            vec!["--initial-field", "a", "--initial-field", "b"]]
        {
            assert!(options(&bad.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
        }
    }

    #[test]
    fn entire_loaded_trace_is_required_not_just_three_probes() {
        let fixture = Cantilever { load: 1.0, band: 0.375 };
        let mut geometry = GridSdf::from_fn(8, &|_, _| -1.0);
        assert!(require_load_support(&geometry, fixture).is_ok());
        *geometry.node_mut(8, 2) = 1.0;
        for y in [0.125, 0.5, 0.875] { assert!(geometry.value_at([1.0, y]) < 0.0); }
        assert!(require_load_support(&geometry, fixture).is_err());
        assert!(require_load_support(&GridSdf::from_fn(8, &|x, _| x - 0.8), fixture).is_err());
    }

    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!("fs-stress-input-{}-{}", std::process::id(),
                NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
    }
    fn arguments(output: &Path, source: &Path) -> Vec<String> {
        let mut args = vec![output.to_string_lossy().into_owned()];
        args.extend(["1e12", "3", "2", "0.6", "8", "0", "300", "--initial-field"].map(str::to_string));
        args.push(source.to_string_lossy().into_owned());
        args.extend(["--load", "2", "--youngs", "2", "--poisson", "0.3", "--load-band", "0.125",
            "--checkpoint"].map(str::to_string));
        args
    }
    fn invoke(args: &[String]) -> Result<u8, Box<dyn Error>> {
        super::super::run_at(args, Instant::now())
    }

    #[test]
    fn authored_study_accepts_and_resumes_without_the_original_input_file() {
        let scratch = Scratch::new();
        let source = scratch.0.join("source.csv");
        let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        field(&source, &geometry).unwrap();
        let original = std::fs::read(&source).unwrap();
        let full = scratch.0.join("full");
        let full_exit = invoke(&arguments(&full, &source)).unwrap();
        assert!(matches!(full_exit, 0 | 11));
        let first = scratch.0.join("first");
        let mut args = arguments(&first, &source);
        args.extend(["--pause-after", "1"].map(str::to_string));
        assert_eq!(invoke(&args).unwrap(), 6, "must accept a real authored update before pausing");
        assert_eq!(std::fs::read(first.join("input-level-set.csv")).unwrap(), original);
        assert_eq!(std::fs::read(&source).unwrap(), original);
        // Recovery must depend on the exact checkpoint, not re-import the CSV.
        std::fs::rename(&source, source.with_extension("retained-source")).unwrap();
        let second = scratch.0.join("second");
        let checkpoint = first.join("checkpoint-0001.fscp");
        let bytes = std::fs::read(&checkpoint).unwrap();
        let args = vec!["--resume".into(), checkpoint.to_string_lossy().into_owned(),
            second.to_string_lossy().into_owned(), "--wall-seconds".into(), "300".into()];
        assert_eq!(invoke(&args).unwrap(), full_exit);
        assert_eq!(std::fs::read(checkpoint).unwrap(), bytes);
        assert_eq!(std::fs::read(full.join("level-set.csv")).unwrap(),
            std::fs::read(second.join("level-set.csv")).unwrap());
        for name in ["trajectory.jsonl", "attempts.jsonl", "stress-checks.jsonl"] {
            let mut joined = std::fs::read(first.join(name)).unwrap();
            joined.extend(std::fs::read(second.join(name)).unwrap());
            assert_eq!(std::fs::read(full.join(name)).unwrap(), joined, "{name}");
        }
        let summary = std::fs::read_to_string(second.join("summary.json")).unwrap();
        assert!(summary.contains(&format!("\"youngs\":{:.17e}", 2.0)));
        assert!(summary.contains(&format!("\"traction_y\":{:.17e}", -2.0)));
        assert!(summary.contains("\"initial_field\":null"));
    }

    #[test]
    fn authored_geometry_and_material_match_independent_baseline_physics() {
        let scratch = Scratch::new();
        let source = scratch.0.join("sloped.csv");
        let geometry = GridSdf::from_fn(8, &|x, y| (y - 0.5).abs() - (0.33 + 0.04 * x));
        field(&source, &geometry).unwrap();
        let output = scratch.0.join("baseline");
        let mut args = vec![output.to_string_lossy().into_owned()];
        args.extend(["1e12", "3", "1", "0.6", "8", "0", "300", "--initial-field"].map(str::to_string));
        args.push(source.to_string_lossy().into_owned());
        args.extend(["--load", "2", "--youngs", "3", "--poisson", "0.2", "--load-band", "0.1",
            "--pause-after", "0"].map(str::to_string));
        assert_eq!(invoke(&args).unwrap(), 6);
        let baseline = read_field(&output.join("baseline-level-set.csv"), 8).unwrap();
        for j in 0..=8 { for i in [0, 8] {
            assert_eq!(baseline.node(i, j).to_bits(), geometry.node(i, j).to_bits());
        }}
        let settings = OptimizeSettings { level: 3, youngs: 3.0, poisson: 0.2, ..OptimizeSettings::default() };
        let expected = fs_topols::evaluate_sampled_stress(&baseline,
            Cantilever { load: 2.0, band: 0.1 }, settings).unwrap();
        let summary = std::fs::read_to_string(output.join("summary.json")).unwrap();
        assert!(summary.contains(&format!("\"baseline\":{}", stress_json(&expected))));
        assert_eq!(std::fs::read(&source).unwrap(), std::fs::read(output.join("input-level-set.csv")).unwrap());
        let wrong = fs_topols::evaluate_sampled_stress(&baseline,
            Cantilever { load: 1.0, band: 0.125 }, OptimizeSettings { level: 3, ..OptimizeSettings::default() }).unwrap();
        assert_ne!(wrong.compliance.to_bits(), expected.compliance.to_bits(), "test must detect ignored inputs");
    }

    #[test]
    fn invalid_csv_or_lost_load_refuses_before_study_publication() {
        let scratch = Scratch::new();
        for (index, geometry) in [GridSdf::from_fn(8, &|x, _| x - 0.7),
            GridSdf::from_fn(8, &|_, y| y - 0.5)].into_iter().enumerate()
        {
            let source = scratch.0.join(format!("lost-load-{index}.csv"));
            field(&source, &geometry).unwrap();
            let output = scratch.0.join(format!("refused-{index}"));
            let error = invoke(&arguments(&output, &source)).unwrap_err();
            assert!(error.to_string().contains("load support"), "{error}");
            assert!(!output.exists());
        }
        let source = scratch.0.join("bad.csv");
        std::fs::write(&source, "x_normalized,y_normalized,phi_normalized\n0.125,0,-1\n").unwrap();
        let output = scratch.0.join("bad-study");
        assert!(invoke(&arguments(&output, &source)).is_err());
        assert!(!output.exists());
    }

    #[test]
    fn resume_refuses_authored_problem_overrides_before_reading_files() {
        let scratch = Scratch::new();
        for (index, flag) in ["--initial-field", "--load", "--load-band", "--youngs", "--poisson"].iter().enumerate() {
            let output = scratch.0.join(format!("refused-{index}"));
            let args = vec!["--resume".into(), "missing.fscp".into(), output.to_string_lossy().into_owned(),
                "--wall-seconds".into(), "300".into(), flag.to_string(), "1".into()];
            let error = invoke(&args).unwrap_err();
            assert!(error.to_string().contains("immutable"), "{error}");
            assert!(!output.exists());
        }
    }
}
