//! Product path for feasible-baseline, hard-area-constrained elasticity descent.
use super::{json_string, parse, writer};
use fs_topols::projected::{
    ProjectedAttempt, ProjectedOptimizer, ProjectedProgress, ProjectedSettings,
};
use fs_topols::volume::VolumeProjectionSettings;
use fs_topols::{Cantilever, GridSdf, OptimizeSettings};
use std::collections::BTreeSet;
use std::error::Error;
use std::io::Write;
use std::path::{Path, PathBuf};

const USAGE: &str = "usage: fs-marquee-elasticity --projected OUTPUT_DIR [LEVEL=4] [ITERATIONS=30] [VOLFRAC=0.45] [MAX_CANDIDATES=6] [--initial-field CSV] [--youngs E] [--poisson NU] [--load TRACTION] [--load-band HALF_WIDTH]";

/// The geometry and material declarations are inputs to the actual PDE, not
/// merely report labels. All coordinates and loads retain the normalized model.
struct Options {
    output_dir: PathBuf,
    initial_field: Option<PathBuf>,
    level: u32,
    iterations: usize,
    volfrac: f64,
    max_candidates: usize,
    youngs: f64,
    poisson: f64,
    load: f64,
    load_band: f64,
}

impl Options {
    fn parse(args: &[String]) -> Result<Self, Box<dyn Error>> {
        let defaults = OptimizeSettings::default();
        let mut options = Self {
            output_dir: PathBuf::new(), initial_field: None, level: 4,
            iterations: 30, volfrac: 0.45, max_candidates: 6,
            youngs: defaults.youngs, poisson: defaults.poisson,
            load: 1.0, load_band: 0.125,
        };
        let mut positional = Vec::new();
        let mut seen = BTreeSet::new();
        let mut index = 0;
        while index < args.len() {
            let name = args[index].as_str();
            if name.starts_with("--") {
                if !matches!(name, "--initial-field" | "--youngs" | "--poisson" | "--load" | "--load-band") {
                    return Err(format!("unknown projected option {name}; {USAGE}").into());
                }
                if !seen.insert(name) {
                    return Err(format!("duplicate projected option {name}").into());
                }
                index += 1;
                let value = args.get(index).filter(|value| !value.starts_with("--"))
                    .ok_or_else(|| format!("missing value for {name}"))?;
                match name {
                    "--initial-field" => options.initial_field = Some(PathBuf::from(value)),
                    "--youngs" => options.youngs = value.parse()?,
                    "--poisson" => options.poisson = value.parse()?,
                    "--load" => options.load = value.parse()?,
                    "--load-band" => options.load_band = value.parse()?,
                    _ => unreachable!("option name checked above"),
                }
            } else {
                positional.push(args[index].clone());
            }
            index += 1;
        }
        if positional.is_empty() || positional.len() > 5 { return Err(USAGE.into()); }
        options.output_dir = PathBuf::from(&positional[0]);
        options.level = parse(&positional, 1, options.level)?;
        options.iterations = parse(&positional, 2, options.iterations)?;
        options.volfrac = parse(&positional, 3, options.volfrac)?;
        options.max_candidates = parse(&positional, 4, options.max_candidates)?;
        if !(2..=7).contains(&options.level) || !(1..=200).contains(&options.iterations)
            || !(options.volfrac.is_finite() && options.volfrac > 0.001 && options.volfrac < 0.999)
            || !(1..=16).contains(&options.max_candidates)
        {
            return Err("projected mode requires LEVEL in [2,7], ITERATIONS in [1,200], VOLFRAC in (0.001,0.999), and MAX_CANDIDATES in [1,16]".into());
        }
        if !(options.youngs.is_finite() && options.youngs > 0.0
            && options.poisson.is_finite() && options.poisson > -1.0 && options.poisson <= 1.0 / 3.0
            && options.load.is_finite() && options.load > 0.0
            && options.load_band.is_finite() && options.load_band > 0.0 && options.load_band <= 0.5)
        {
            return Err("projected material/load requires finite E > 0, -1 < NU <= 1/3, TRACTION > 0, and 0 < HALF_WIDTH <= 0.5".into());
        }
        Ok(options)
    }

    fn geometry(&self) -> Result<GridSdf, Box<dyn Error>> {
        let n = 1usize << self.level;
        match &self.initial_field {
            Some(path) => fs_marquee::level_set_csv::read_field(path, n),
            None => Ok(GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.42)),
        }
    }

    fn settings(&self) -> OptimizeSettings {
        OptimizeSettings {
            level: self.level, volfrac: self.volfrac, iterations: self.iterations,
            move_cells: 0.2, nucleation_period: 4, hole_radius_cells: 1.5,
            youngs: self.youngs, poisson: self.poisson, ..OptimizeSettings::default()
        }
    }

    fn fixture(&self) -> Cantilever { Cantilever { load: self.load, band: self.load_band } }
}

fn field(path: &Path, phi: &GridSdf) -> Result<(), Box<dyn Error>> {
    let mut output = writer(path)?;
    writeln!(output, "x_normalized,y_normalized,phi_normalized")?;
    for j in 0..=phi.n() {
        for i in 0..=phi.n() {
            let [x, y] = phi.pos(i, j);
            writeln!(output, "{x:.17e},{y:.17e},{:.17e}", phi.node(i, j))?;
        }
    }
    output.flush()?;
    Ok(())
}

fn attempt_rows(
    output: &mut impl Write,
    iteration: usize,
    attempts: &[ProjectedAttempt],
) -> std::io::Result<()> {
    for attempt in attempts {
        let state = attempt.state.map_or_else(|| "null".to_string(), |state| format!(
            "{{\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\"}}",
            state.compliance, state.volume, state.snapshot,
        ));
        let projection = attempt.projection.map_or_else(|| "null".to_string(), |projection| format!(
            "{{\"shift\":{:.17e},\"volume\":{:.17e},\"evaluations\":{}}}",
            projection.shift, projection.volume, projection.evaluations,
        ));
        let refusal = attempt.refusal.as_deref().map_or_else(|| "null".to_string(), json_string);
        writeln!(output,
            "{{\"iter\":{iteration},\"candidate\":{},\"move_cells\":{:.17e},\"projection\":{projection},\"state\":{state},\"refusal\":{refusal}}}",
            attempt.index, attempt.move_cells,
        )?;
    }
    Ok(())
}

pub(super) fn run(args: &[String]) -> Result<u8, Box<dyn Error>> {
    let options = Options::parse(args)?;
    let output_dir = &options.output_dir;
    if output_dir.try_exists()? {
        return Err("output directory already exists; refusing to overwrite it".into());
    }
    let Options { level, iterations, volfrac, max_candidates, .. } = options;
    let geometry = options.geometry()?;
    let n = geometry.n();
    // Hold the actual imported supported/loaded boundary traces fixed. No
    // material can disappear from a loaded edge merely to reduce external work.
    let fixed: Vec<_> = geometry.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % (n + 1) == 0 || index % (n + 1) == n).collect();
    let settings = options.settings();
    let projection = VolumeProjectionSettings {
        target: volfrac, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    };
    let controls = ProjectedSettings { max_candidates, ..ProjectedSettings::default() };
    let mut optimizer = ProjectedOptimizer::new(geometry,
        options.fixture(), settings, fixed, projection, controls)?;
    let baseline = optimizer.baseline();
    std::fs::create_dir(output_dir)?;
    field(&output_dir.join("baseline-level-set.csv"), optimizer.checkpoint().geometry())?;
    let mut trajectory = writer(&output_dir.join("trajectory.jsonl"))?;
    let mut attempts = writer(&output_dir.join("attempts.jsonl"))?;
    let mut proposals = writer(&output_dir.join("unprojected-proposals.jsonl"))?;
    let mut candidate_count = 0usize;
    let status = loop {
        let iteration = optimizer.checkpoint().next_iteration();
        match optimizer.advance_one()? {
            ProjectedProgress::Accepted(step) => {
                candidate_count += step.attempts.len();
                attempt_rows(&mut attempts, step.iteration, &step.attempts)?;
                writeln!(trajectory,
                    "{{\"iter\":{},\"previous_compliance\":{:.17e},\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\",\"projection_shift\":{:.17e},\"projection_evaluations\":{},\"authority\":\"estimated\"}}",
                    step.iteration, step.previous.compliance, step.state.compliance,
                    step.state.volume, step.state.snapshot, step.projection.shift,
                    step.projection.evaluations,
                )?;
                for row in &step.proposal.rows { writeln!(proposals, "{row}")?; }
                trajectory.flush()?;
                attempts.flush()?;
                proposals.flush()?;
            }
            ProjectedProgress::IterationLimit => break "iteration_limit",
            ProjectedProgress::NoDescent(rows) => {
                candidate_count += rows.len();
                attempt_rows(&mut attempts, iteration, &rows)?;
                break "no_descent";
            }
        }
    };
    trajectory.flush()?;
    attempts.flush()?;
    proposals.flush()?;
    field(&output_dir.join("level-set.csv"), optimizer.checkpoint().geometry())?;
    let current = optimizer.current();
    let reduction = if baseline.compliance > 0.0 {
        (baseline.compliance - current.compliance) / baseline.compliance
    } else { 0.0 };
    let initial_field = options.initial_field.as_ref().map_or_else(|| "null".to_string(),
        |path| json_string(&path.to_string_lossy()));
    let summary = format!(
        concat!(
            "{{\"schema\":\"fs-marquee-projected-v1\",\"model\":\"normalized_unit_square_plane_strain_cantilever\",",
            "\"authority\":\"estimated\",\"status\":\"{}\",\"level\":{},\"requested_updates\":{},",
            "\"accepted_updates\":{},\"candidate_count\":{},\"candidate_budget_per_update\":{},",
            "\"volume_target\":{:.17e},\"volume_tolerance\":{:.17e},",
            "\"baseline_compliance\":{:.17e},\"baseline_volume\":{:.17e},",
            "\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\",",
            "\"relative_reduction_from_feasible_baseline\":{:.17e},",
            "\"initial_field\":{},\"youngs\":{:.17e},\"poisson\":{:.17e},",
            "\"load\":{:.17e},\"load_band\":{:.17e},\"fixed_boundaries\":[\"left\",\"right\"],",
            "\"claims\":{{\"converged\":false,\"global_optimum\":false,\"physical_validation\":false,\"certified_continuum_volume\":false}}}}"
        ),
        status, level, iterations, optimizer.checkpoint().next_iteration(), candidate_count,
        max_candidates, volfrac, projection.tolerance, baseline.compliance, baseline.volume,
        current.compliance, current.volume, current.snapshot, reduction,
        initial_field, options.youngs, options.poisson, options.load, options.load_band,
    );
    // Create the completion marker only after every field and trace export.
    let mut completion = writer(&output_dir.join("summary.json"))?;
    writeln!(completion, "{summary}")?;
    completion.flush()?;
    println!("{summary}");
    Ok(if status == "iteration_limit" { 0 } else { 11 })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|v| (*v).into()).collect() }

    #[test]
    fn legacy_defaults_and_custom_physics_reach_solver_inputs() {
        let original = Options::parse(&args(&["out"])).unwrap();
        assert_eq!(original.settings().youngs, OptimizeSettings::default().youngs);
        assert_eq!(original.fixture().load, 1.0);
        assert!(original.initial_field.is_none());
        let custom = Options::parse(&args(&["out", "3", "2", "0.6", "4",
            "--initial-field", "beam.csv", "--youngs", "2.5", "--poisson", "-0.2",
            "--load", "3", "--load-band", "0.2"])).unwrap();
        assert_eq!(custom.settings().level, 3);
        assert_eq!(custom.settings().iterations, 2);
        assert_eq!(custom.settings().volfrac, 0.6);
        assert_eq!(custom.settings().youngs, 2.5);
        assert_eq!(custom.settings().poisson, -0.2);
        assert_eq!(custom.fixture().load, 3.0);
        assert_eq!(custom.fixture().band, 0.2);
        assert_eq!(custom.initial_field, Some(PathBuf::from("beam.csv")));
    }

    #[test]
    fn malformed_or_ambiguous_inputs_are_refused() {
        for values in [
            vec!["out", "--unknown", "1"], vec!["out", "--initial-field"],
            vec!["out", "--load", "1", "--load", "2"],
            vec!["out", "--load", "--youngs", "2"],
            vec!["out", "--youngs", "NaN"], vec!["out", "--youngs", "0"],
            vec!["out", "--poisson", "0.4"], vec!["out", "--poisson", "-1"],
            vec!["out", "--load", "inf"], vec!["out", "--load", "-1"],
            vec!["out", "--load-band", "0"], vec!["out", "--load-band", "0.6"],
        ] {
            assert!(Options::parse(&args(&values)).is_err(), "{values:?}");
        }
    }

    #[test]
    fn projected_options_refuse_before_creating_output() {
        let path = std::env::temp_dir().join(format!("frankensim-projected-refusal-{}", std::process::id()));
        assert!(!path.exists(), "test output must not preexist");
        let args = vec![path.to_string_lossy().into_owned(), "40".into()];
        assert!(run(&args).is_err());
        assert!(!path.exists());
    }

    #[test]
    fn imported_geometry_is_used_by_the_actual_projected_solver() {
        let root = std::env::temp_dir().join(format!("frankensim-projected-input-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        let input = root.join("input.csv");
        let output = root.join("run");
        let geometry = GridSdf::from_fn(8, &|x, y| ((y - 0.5).abs() - 0.42) * (1.0 + 0.25 * x));
        field(&input, &geometry).unwrap();
        let arguments = vec![output.to_string_lossy().into_owned(), "3".into(), "1".into(),
            "0.6".into(), "1".into(), "--initial-field".into(), input.to_string_lossy().into_owned(),
            "--youngs".into(), "2".into(), "--load".into(), "1.5".into()];
        assert!(matches!(run(&arguments).unwrap(), 0 | 11));
        let baseline = fs_marquee::level_set_csv::read_field(&output.join("baseline-level-set.csv"), 8).unwrap();
        for j in 0..=8 {
            for i in [0, 8] {
                assert_eq!(baseline.node(i, j).to_bits(), geometry.node(i, j).to_bits());
            }
        }
        let summary = std::fs::read_to_string(output.join("summary.json")).unwrap();
        assert!(summary.contains(&format!("\"youngs\":{:.17e}", 2.0)));
        assert!(summary.contains(&format!("\"load\":{:.17e}", 1.5)));
        assert!(summary.contains(&json_string(&input.to_string_lossy())));
    }
}
