//! Runnable SI-unit, two-load, fixed-mesh elasticity topology study.
//! See examples/elastic-topology/README.md. No external runtime dependencies.
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};
use std::time::Instant;

use fs_topopt::pipeline::LoadCase;
use fs_topopt::{
    DensityElasticity, DensityFilter, DesignPipeline, EvaluationStop, MultiLoadOcOptions,
    MultiLoadOcReport, MultiLoadOcTermination, SimpParams, SolveBudget, SolveControl,
    controlled_multi_load_optimality_criteria,
};

const HELP: &str = "Usage: elastic_topology --output NEW_DIRECTORY [options]
  --cells N          Kuhn subdivisions per axis, 1..8 (default 2)
  --iterations N     Accepted-update budget, 0..200 (default 30)
  --volume F         Projected material fraction, 0.01..1 (default 0.5)
  --length M         Cantilever length in meters (default 0.2)
  --width M          Width in meters (default 0.1)
  --height M         Height in meters (default 0.1)
  --youngs PA        Solid Young's modulus in pascals (default 7e10)
  --poisson NU       Poisson ratio, -0.9..0.49 (default 0.3)
  --filter-radius M  Helmholtz radius in meters (default 0.015)
  --load-y N         Signed total y-directed end load in newtons (default 1)
  --load-z N         Signed total z-directed end load in newtons (default -1)
  --weight-y W       y-case weight, 0..1; z weight is 1-W (default 0.5)
  --seconds S        Optional wall-time stop, including baseline preparation
  --linear-iterations N        Per-solve Krylov cap, 0..50000 (default 50000)
  --total-linear-iterations N  Whole-run Krylov cap, 0..50000000 (default 2000000)
Loads are independent cases, each shared equally among x=length end nodes.
Coordinates/displacements use m; force uses N; compliance uses N*m.
This exports a density field, not a binary solid or certified optimum.
Solves poll cancellation every 32 Krylov iterations, not at fixed wall intervals.";

#[derive(Debug, Clone)]
struct Config {
    output: PathBuf,
    cells: usize,
    iterations: usize,
    volume: f64,
    length: f64,
    width: f64,
    height: f64,
    youngs: f64,
    poisson: f64,
    radius: f64,
    load_y: f64,
    load_z: f64,
    weight_y: f64,
    seconds: Option<f64>,
    linear_iterations: usize,
    total_linear_iterations: usize,
}

fn parse(args: impl IntoIterator<Item = String>) -> Result<Option<Config>, String> {
    let mut config = Config {
        output: PathBuf::new(), cells: 2, iterations: 30, volume: 0.5,
        length: 0.2, width: 0.1, height: 0.1, youngs: 7e10, poisson: 0.3,
        radius: 0.015, load_y: 1.0, load_z: -1.0, weight_y: 0.5, seconds: None,
        linear_iterations: 50_000, total_linear_iterations: 2_000_000,
    };
    let mut args = args.into_iter();
    let mut seen = std::collections::BTreeSet::new();
    while let Some(flag) = args.next() {
        if flag == "--help" || flag == "-h" { return Ok(None); }
        if !seen.insert(flag.clone()) { return Err(format!("duplicate option: {flag}")); }
        let value = args.next().ok_or_else(|| format!("missing value for {flag}"))?;
        let number = || value.parse::<f64>().map_err(|_| format!("invalid number for {flag}: {value}"));
        match flag.as_str() {
            "--output" => config.output = PathBuf::from(&value),
            "--cells" => config.cells = value.parse().map_err(|_| "--cells requires an integer")?,
            "--iterations" => config.iterations = value.parse().map_err(|_| "--iterations requires an integer")?,
            "--volume" => config.volume = number()?,
            "--length" => config.length = number()?,
            "--width" => config.width = number()?,
            "--height" => config.height = number()?,
            "--youngs" => config.youngs = number()?,
            "--poisson" => config.poisson = number()?,
            "--filter-radius" => config.radius = number()?,
            "--load-y" => config.load_y = number()?,
            "--load-z" => config.load_z = number()?,
            "--weight-y" => config.weight_y = number()?,
            "--seconds" => config.seconds = Some(number()?),
            "--linear-iterations" => config.linear_iterations = value.parse()
                .map_err(|_| "--linear-iterations requires an integer")?,
            "--total-linear-iterations" => config.total_linear_iterations = value.parse()
                .map_err(|_| "--total-linear-iterations requires an integer")?,
            _ => return Err(format!("unknown option: {flag}")),
        }
    }
    if config.output.as_os_str().is_empty() { return Err("--output NEW_DIRECTORY is required".into()); }
    if !(1..=8).contains(&config.cells) || config.iterations > 200
        || config.linear_iterations > 50_000 || config.total_linear_iterations > 50_000_000
    {
        return Err("resource admission requires 1..8 subdivisions, at most 200 updates, 50000 per-solve and 50000000 total Krylov iterations".into());
    }
    let inside = |x: f64, lo: f64, hi: f64| x.is_finite() && (lo..=hi).contains(&x);
    if !inside(config.volume, 0.01, 1.0)
        || !inside(config.youngs, 1e3, 1e13)
        || !inside(config.poisson, -0.9, 0.49)
        || !inside(config.weight_y, 0.0, 1.0)
        || [config.length, config.width, config.height].iter().any(|&x| !inside(x, 1e-4, 1e3))
        || !inside(config.radius, 0.0, config.length.max(config.width).max(config.height))
        || [config.load_y, config.load_z].iter().any(|&x| !inside(x, -1e9, 1e9))
        || config.seconds.is_some_and(|s| !inside(s, 1e-3, 86_400.0))
    {
        return Err("non-finite or unsupported physical input: dimensions 1e-4..1e3 m, E 1e3..1e13 Pa, loads +/-1e9 N, time 1e-3..86400 s; see --help for other bounds".into());
    }
    if !(config.weight_y > 0.0 && config.load_y != 0.0
        || config.weight_y < 1.0 && config.load_z != 0.0)
    {
        return Err("at least one positive-weight load must be nonzero".into());
    }
    Ok(Some(config))
}

fn physical_volume(
    pipeline: &DesignPipeline, rho: &[f64], volumes: &[f64], control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    let total: f64 = volumes.iter().sum();
    Ok(pipeline.try_forward(rho, control)?.1.iter().zip(volumes).map(|(r, v)| r * (v / total)).sum())
}

fn uniform_baseline(
    pipeline: &DesignPipeline, volumes: &[f64], cap: f64, control: &mut SolveControl<'_>,
) -> Result<Vec<f64>, EvaluationStop> {
    // Raw density is not projected volume, even for a uniform field when
    // cap != eta. Invert the ACTUAL filter/projection, retaining feasibility.
    let mut lower = vec![1e-3; volumes.len()];
    assert!(physical_volume(pipeline, &lower, volumes, control)? <= cap, "density floor exceeds volume cap");
    let (mut lo, mut hi) = (1e-3, 1.0);
    for _ in 0..56 {
        let mid = 0.5 * (lo + hi);
        let candidate = vec![mid; volumes.len()];
        if physical_volume(pipeline, &candidate, volumes, control)? <= cap {
            lo = mid;
            lower = candidate;
        } else { hi = mid; }
    }
    Ok(lower)
}

fn fresh_writer(path: &Path) -> io::Result<BufWriter<File>> {
    Ok(BufWriter::new(File::options().write(true).create_new(true).open(path)?))
}

fn vtk(
    out: &mut impl Write,
    positions: &[[f64; 3]],
    cells: &[[u32; 4]],
    report: &MultiLoadOcReport,
) -> io::Result<()> {
    writeln!(out, "# vtk DataFile Version 3.0\nFrankenSim SI fixed-mesh SIMP; not a binary solid\nASCII\nDATASET UNSTRUCTURED_GRID")?;
    writeln!(out, "POINTS {} double", positions.len())?;
    for p in positions { writeln!(out, "{:.17e} {:.17e} {:.17e}", p[0], p[1], p[2])?; }
    writeln!(out, "CELLS {} {}", cells.len(), 5 * cells.len())?;
    for t in cells { writeln!(out, "4 {} {} {} {}", t[0], t[1], t[2], t[3])?; }
    writeln!(out, "CELL_TYPES {}", cells.len())?;
    for _ in cells { writeln!(out, "10")?; }
    writeln!(out, "POINT_DATA {}", positions.len())?;
    for (name, displacement) in ["displacement_y_case_m", "displacement_z_case_m"]
        .iter().zip(&report.displacements)
    {
        writeln!(out, "VECTORS {name} double")?;
        for u in displacement.chunks_exact(3) {
            writeln!(out, "{:.17e} {:.17e} {:.17e}", u[0], u[1], u[2])?;
        }
    }
    writeln!(out, "CELL_DATA {}", cells.len())?;
    for (name, field) in [("raw_density", &report.rho), ("projected_density", &report.projected_rho)] {
        writeln!(out, "SCALARS {name} double 1\nLOOKUP_TABLE default")?;
        for value in field { writeln!(out, "{value:.17e}")?; }
    }
    Ok(())
}

fn terminal(reason: MultiLoadOcTermination) -> &'static str {
    match reason {
        MultiLoadOcTermination::IterationBudget => "iteration-budget",
        MultiLoadOcTermination::DesignChange => "design-change-threshold",
        MultiLoadOcTermination::Cancelled => "wall-time-stop",
        MultiLoadOcTermination::LinearBudget => "linear-work-budget",
        MultiLoadOcTermination::NumericalFailure => "numerical-failure",
        MultiLoadOcTermination::NoAcceptableStep => "no-acceptable-step",
    }
}

fn export(
    config: &Config,
    positions: &[[f64; 3]],
    cells: &[[u32; 4]],
    report: &MultiLoadOcReport,
    options: MultiLoadOcOptions,
) -> io::Result<()> {
    if report.history.is_empty() || report.rho.len() != cells.len()
        || report.projected_rho.len() != cells.len() || report.displacements.len() != 2
        || report.displacements.iter().any(|u| u.len() != 3 * positions.len())
        || report.rho.iter().chain(&report.projected_rho)
            .chain(report.displacements.iter().flatten()).any(|v| !v.is_finite())
        || report.history.iter().any(|row| row.case_compliances.len() != 2
            || !row.compliance.is_finite() || !row.volume_fraction.is_finite()
            || !row.max_change.is_finite()
            || row.case_compliances.iter().any(|c| !c.is_finite()))
        || report.history.first().is_none_or(|row| row.compliance <= 0.0)
    {
        return Err(io::Error::other("no complete accepted two-load state to export"));
    }
    // The caller must supply a NEW directory with an existing parent. No
    // artifact is overwritten. An I/O failure can leave a partial directory,
    // but cannot emit the success message or a final summary.
    fs::create_dir(&config.output)?;
    let mut mesh = fresh_writer(&config.output.join("design.vtk"))?;
    vtk(&mut mesh, positions, cells, report)?;
    mesh.flush()?;
    let mut csv = fresh_writer(&config.output.join("iterations.csv"))?;
    writeln!(csv, "iteration,weighted_compliance_N_m,y_compliance_N_m,z_compliance_N_m,projected_volume_fraction,max_raw_density_change")?;
    for row in &report.history {
        writeln!(csv, "{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            row.iteration, row.compliance, row.case_compliances[0], row.case_compliances[1],
            row.volume_fraction, row.max_change)?;
    }
    csv.flush()?;
    // Write summary LAST. It contains every modeling control, but is not a
    // ledger package, a checkpoint/resume file, or a certificate.
    let first = &report.history[0];
    let last = report.history.last().expect("nonempty history checked above");
    let mut summary = fresh_writer(&config.output.join("summary.txt"))?;
    writeln!(summary, "FrankenSim {} fixed-mesh two-load elasticity topology study", fs_topopt::VERSION)?;
    writeln!(summary, "terminal={}\naccepted_updates={}\ntetrahedra={}\nsubdivisions_per_axis={}",
        terminal(report.termination), last.iteration, cells.len(), config.cells)?;
    writeln!(summary, "length_m={:.17e}\nwidth_m={:.17e}\nheight_m={:.17e}\nyoungs_Pa={:.17e}\npoisson={:.17e}\nfilter_radius_m={:.17e}",
        config.length, config.width, config.height, config.youngs, config.poisson, config.radius)?;
    writeln!(summary, "y_end_force_N={:.17e}\nz_end_force_N={:.17e}\ny_case_weight={:.17e}\nz_case_weight={:.17e}",
        config.load_y, config.load_z, config.weight_y, 1.0 - config.weight_y)?;
    writeln!(summary, "projection_beta=2\nprojection_eta=0.5\nsimp_penal=3\nrelative_void_modulus=1e-6\nraw_density_floor=1e-3")?;
    writeln!(summary, "volume_cap={:.17e}\nmove_limit={:.17e}\nmax_updates={}\nchange_tolerance={:.17e}\nvolume_tolerance={:.17e}\nmax_backtracks={}\nwall_time_seconds={:?}",
        options.volume_fraction, options.move_limit, options.max_iterations,
        options.change_tolerance, options.volume_tolerance, options.max_backtracks, config.seconds)?;
    writeln!(summary, "per_solve_iteration_budget={}\ntotal_linear_iteration_budget={}\nlinear_solves={}\nlinear_iterations={}\nevaluation_stop={:?}",
        config.linear_iterations, config.total_linear_iterations, report.work.linear_solves,
        report.work.linear_iterations, report.evaluation_stop)?;
    writeln!(summary, "initial_compliance_N_m={:.17e}\nfinal_compliance_N_m={:.17e}\nfinal_projected_volume_fraction={:.17e}\nrelative_compliance_reduction={:.17e}",
        first.compliance, last.compliance, last.volume_fraction,
        (first.compliance - last.compliance) / first.compliance)?;
    writeln!(summary, "deterministic_no_rng=true\nartifacts=design.vtk,iterations.csv\nAll end nodes share each total force equally; loads are NOT simultaneous.\nNo binary-solid, stationarity, mesh-convergence, stress-safety, experimental-validation, or free-boundary claim.\nWall/work limits cover baseline and optimization; polls bound Krylov batches, not milliseconds.\nCG convergence uses recursive residual estimates, not certified Euclidean error bounds.\nNo ledger integration or resumable checkpoint is supplied by this example.")?;
    summary.flush()
}

fn main() -> Result<(), Box<dyn Error>> {
    let Some(config) = parse(std::env::args().skip(1)).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))? else {
        println!("{HELP}");
        return Ok(());
    };
    if config.output.try_exists()? {
        return Err(io::Error::new(io::ErrorKind::AlreadyExists, "output directory must be new").into());
    }
    // Bounds are admitted before mesh/filter/elasticity allocations.
    let (complex, mut positions) = fs_feec::kuhn_cube(config.cells);
    for p in &mut positions { p[0] *= config.length; p[1] *= config.width; p[2] *= config.height; }
    let mut elasticity = DensityElasticity::new(&complex, &positions, config.youngs,
        config.poisson, &|p| p[0] == 0.0);
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, config.radius),
        params: SimpParams::default(),
    };
    let volumes: Vec<f64> = fs_feec::element_geometry(&complex, &positions)
        .vol_signed.iter().map(|v| v.abs()).collect();
    let end: Vec<usize> = positions.iter().enumerate()
        .filter(|(_, p)| (p[0] - config.length).abs() <= config.length * 1e-12)
        .map(|(i, _)| i).collect();
    let mut fy = vec![0.0; elasticity.n()];
    let mut fz = vec![0.0; elasticity.n()];
    for &i in &end {
        fy[3 * i + 1] = config.load_y / end.len() as f64;
        fz[3 * i + 2] = config.load_z / end.len() as f64;
    }
    let options = MultiLoadOcOptions {
        volume_fraction: config.volume, max_iterations: config.iterations,
        ..MultiLoadOcOptions::default()
    };
    // One control covers baseline inversion, accepted designs and rejected
    // trials. Matrix assembly and output I/O remain outside this wall budget.
    let start = Instant::now();
    let mut callback = |_| {
        if config.seconds.is_some_and(|limit| start.elapsed().as_secs_f64() >= limit) {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget {
        per_solve_iterations: config.linear_iterations,
        total_iterations: config.total_linear_iterations,
    }, &mut callback);
    let rho0 = uniform_baseline(&pipeline, &volumes, config.volume, &mut control)?;
    let report = controlled_multi_load_optimality_criteria(&pipeline, &mut elasticity, &[
        LoadCase { force: &fy, weight: config.weight_y },
        LoadCase { force: &fz, weight: 1.0 - config.weight_y },
    ], &rho0, &volumes, options, &mut control);
    if report.history.is_empty() {
        if let Some(stop) = &report.evaluation_stop { return Err(stop.clone().into()); }
    }
    export(&config, &positions, &complex.tets, &report, options)?;
    let final_row = report.history.last().expect("export requires a solved state");
    println!("Retained {}: terminal={}, updates={}, compliance={:.9e} N*m, projected volume={:.9}, Krylov iterations={}",
        config.output.display(), terminal(report.termination), final_row.iteration,
        final_row.compliance, final_row.volume_fraction, report.work.linear_iterations);
    if report.termination == MultiLoadOcTermination::NumericalFailure {
        return Err(io::Error::other(format!("numerical failure; accepted prefix was exported: {:?}", report.evaluation_stop)).into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_owned).collect() }

    #[test]
    fn rejects_nonfinite_inputs_and_oversized_mesh_before_allocation() {
        for bad in ["--youngs NaN", "--volume inf", "--cells 9", "--iterations 201", "--seconds -1",
            "--linear-iterations 50001", "--total-linear-iterations 50000001", "--linear-iterations -1"] {
            assert!(parse(args(&format!("--output test-output {bad}"))).is_err());
        }
    }

    #[test]
    fn refuses_silent_unknown_or_duplicate_controls() {
        assert!(parse(args("--output x --volum 0.4")).is_err());
        assert!(parse(args("--output x --volume 0.4 --volume 0.6")).is_err());
        assert!(parse(args("--output x --load-y 0 --load-z 0")).is_err());
    }

    #[test]
    fn load_settings_and_volume_are_actual_inputs() {
        let c = parse(args("--output x --load-y -2 --load-z 3 --weight-y 0.25 --volume 0.4"))
            .unwrap().unwrap();
        assert_eq!((c.load_y, c.load_z, c.weight_y, c.volume), (-2.0, 3.0, 0.25, 0.4));
    }

    #[test]
    fn linear_budgets_are_explicit_and_zero_is_not_silently_defaulted() {
        let c = parse(args("--output x --linear-iterations 0 --total-linear-iterations 1234"))
            .unwrap().unwrap();
        assert_eq!((c.linear_iterations, c.total_linear_iterations), (0, 1234));
    }
}
