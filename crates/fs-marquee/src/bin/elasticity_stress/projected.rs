//! User-facing same-material, sampled-stress-constrained design search.
use super::{json_string, parse, writer};
use fs_topols::projected::{ProjectedAttempt, ProjectedOptimizer, ProjectedProgress, ProjectedSettings};
use fs_topols::volume::VolumeProjectionSettings;
use fs_topols::{
    Cantilever, GridSdf, OptimizeSettings, ProjectedStressOptimizer,
    SampledStressEvaluation, SampledStressLimit,
};
use std::error::Error;
use std::io::Write;
use std::ops::ControlFlow;
use std::path::Path;
use std::time::{Duration, Instant};

#[path = "projected/checkpoint.rs"]
mod checkpoint;
#[path = "projected/input.rs"]
mod input;
#[path = "projected/regions.rs"]
mod regions;

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
    output.get_ref().sync_all()?;
    Ok(())
}

fn stress_json(state: &SampledStressEvaluation) -> String {
    format!(
        "{{\"compliance\":{:.17e},\"volume\":{:.17e},\"sampled_max_von_mises\":{:.17e},\"max_location\":[{:.17e},{:.17e}],\"sample_count\":{},\"snapshot\":\"{:#018x}\"}}",
        state.compliance, state.volume, state.sampled_max_von_mises,
        state.max_location[0], state.max_location[1], state.sample_count, state.snapshot,
    )
}

fn attempt_rows(output: &mut impl Write, iteration: usize, rows: &[ProjectedAttempt]) -> std::io::Result<()> {
    for row in rows {
        let refusal = row.refusal.as_deref().map_or_else(|| "null".to_string(), json_string);
        let state = row.state.map_or_else(|| "null".to_string(), |state| format!(
            "{{\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\"}}",
            state.compliance, state.volume, state.snapshot,
        ));
        writeln!(output,
            "{{\"iter\":{iteration},\"candidate\":{},\"move_cells\":{:.17e},\"state\":{state},\"refusal\":{refusal}}}",
            row.index, row.move_cells,
        )?;
    }
    Ok(())
}

pub(super) fn run(args: &[String]) -> Result<u8, Box<dyn Error>> {
    run_at(args, Instant::now())
}

fn poll_wall(started: Instant, wall_seconds: u64) -> ControlFlow<()> {
    if started.elapsed() >= Duration::from_secs(wall_seconds) { ControlFlow::Break(()) }
    else { ControlFlow::Continue(()) }
}

fn stopped_before_study(phase: &str) -> u8 {
    eprintln!("wall budget exhausted during {phase}; no study published");
    6
}

// An explicit start also permits deadline regressions without sleeping or
// adding test-only flags to the user-facing command.
fn run_at(args: &[String], started: Instant) -> Result<u8, Box<dyn Error>> {
    if args.first().is_some_and(|arg| arg == "--resume") {
        return checkpoint::resume(&args[1..], started);
    }
    let (args, region_path) = regions::options(args)?;
    let (args, input_options) = input::options(&args)?;
    let (args, checkpoint_options) = checkpoint::options(&args)?;
    if args.len() < 2 || args.len() > 8 {
        return Err("usage: fs-marquee-elasticity-stress --projected OUTPUT_DIR MAX_SAMPLED_VON_MISES [LEVEL=3] [ITERATIONS=30] [VOLFRAC=0.6] [MAX_CANDIDATES=8] [ABS_STRESS_TOL=0] [WALL_SECONDS=300] [--checkpoint] [--pause-after N] [--initial-field FIELD.csv] [--load F] [--load-band HALF_WIDTH] [--youngs E] [--poisson NU] [--design-regions REGIONS.csv]".into());
    }
    let output_dir = Path::new(&args[0]);
    if output_dir.try_exists()? {
        return Err("output directory already exists; refusing to overwrite it".into());
    }
    let limit = SampledStressLimit::new(args[1].parse()?, parse(&args, 6, 0.0)?)?;
    let level: u32 = parse(&args, 2, 3)?;
    let iterations: usize = parse(&args, 3, 30)?;
    let volfrac: f64 = parse(&args, 4, 0.6)?;
    let max_candidates: usize = parse(&args, 5, 8)?;
    let wall_seconds: u64 = parse(&args, 7, 300)?;
    if !(2..=7).contains(&level) || !(1..=200).contains(&iterations)
        || !(volfrac.is_finite() && volfrac > 0.001 && volfrac < 0.999)
        || !(1..=16).contains(&max_candidates) || !(1..=3600).contains(&wall_seconds)
    {
        return Err("projected stress mode requires LEVEL in [2,7], ITERATIONS in [1,200], VOLFRAC in (0.001,0.999), MAX_CANDIDATES in [1,16], WALL_SECONDS in [1,3600]".into());
    }
    if poll_wall(started, wall_seconds).is_break() {
        return Ok(stopped_before_study("initialization"));
    }
    let executable = if checkpoint_options.enabled { Some(checkpoint::executable()?) } else { None };
    let n = 1usize << level;
    let mut settings = OptimizeSettings {
        level, iterations, volfrac, move_cells: 0.1,
        nucleation_period: 0, ..OptimizeSettings::default()
    };
    let (source, fixture) = input_options.prepare(&mut settings)?;
    // Preserve both support and load traces during projection. In particular,
    // disappearing loaded-edge material cannot masquerade as reduced compliance.
    let fixed: Vec<_> = source.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % (n + 1) == 0 || index % (n + 1) == n).collect();
    let authored = match region_path.as_deref() {
        Some(path) => match regions::load_controlled(Path::new(path), &source, &fixed,
            |_| poll_wall(started, wall_seconds),
        )? {
            ControlFlow::Continue(regions) => Some(regions),
            ControlFlow::Break(()) => return Ok(stopped_before_study("design-region authoring")),
        },
        None => None,
    };
    let (geometry, fixed) = match authored.as_ref() {
        Some(regions) => (regions.prepared().geometry.clone(), regions.prepared().fixed_nodes.clone()),
        None => (source.clone(), fixed),
    };
    input::require_load_support(&geometry, fixture)?;
    let projection = VolumeProjectionSettings {
        target: volfrac, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    };
    let controls = ProjectedSettings { max_candidates, ..ProjectedSettings::default() };
    let policy = checkpoint::Policy { fixed: fixed.clone(), projection, controls };
    let area_optimizer = match ProjectedOptimizer::new_controlled(&geometry,
        fixture, settings, fixed, projection, controls,
        |_| poll_wall(started, wall_seconds),
    )? {
        ControlFlow::Continue(optimizer) => optimizer,
        ControlFlow::Break(()) => return Ok(stopped_before_study("area-feasible initialization")),
    };
    let optimizer = match ProjectedStressOptimizer::new_controlled(&area_optimizer, limit,
        |_| poll_wall(started, wall_seconds),
    )? {
        ControlFlow::Continue(optimizer) => optimizer,
        ControlFlow::Break(()) => return Ok(stopped_before_study("baseline stress admission")),
    };
    run_optimizer_with_input(output_dir, optimizer, policy, checkpoint_options, executable,
        started, wall_seconds, false, Some(&source), authored.as_ref())
}

fn run_optimizer(
    output_dir: &Path,
    optimizer: ProjectedStressOptimizer,
    policy: checkpoint::Policy,
    options: checkpoint::Options,
    executable: Option<String>,
    started: Instant,
    wall_seconds: u64,
    resumed: bool,
) -> Result<u8, Box<dyn Error>> {
    run_optimizer_with_input(output_dir, optimizer, policy, options, executable,
        started, wall_seconds, resumed, None, None)
}

fn run_optimizer_with_input(
    output_dir: &Path,
    mut optimizer: ProjectedStressOptimizer,
    policy: checkpoint::Policy,
    options: checkpoint::Options,
    executable: Option<String>,
    started: Instant,
    wall_seconds: u64,
    resumed: bool,
    source: Option<&GridSdf>,
    authored: Option<&regions::Regions>,
) -> Result<u8, Box<dyn Error>> {
    if poll_wall(started, wall_seconds).is_break() {
        return Ok(stopped_before_study("study publication"));
    }
    let wall_budget = Duration::from_secs(wall_seconds);
    let baseline = optimizer.baseline().clone();
    let start_iteration = optimizer.checkpoint().next_iteration();
    let settings = optimizer.checkpoint().settings();
    let fixture = optimizer.checkpoint().fixture();
    let limit = optimizer.limit();
    std::fs::create_dir(output_dir)?;
    if let Some(source) = source { field(&output_dir.join("input-level-set.csv"), source)?; }
    let design_regions = match authored {
        Some(regions) => regions.export(output_dir)?,
        None => "null".to_string(),
    };
    field(&output_dir.join("baseline-level-set.csv"), optimizer.checkpoint().geometry())?;
    let mut last_checkpoint = None;
    if let Some(executable) = executable.as_deref() {
        let name = format!("checkpoint-{start_iteration:04}.fscp");
        checkpoint::save(&output_dir.join(&name), &optimizer, &policy, executable)?;
        last_checkpoint = Some(name);
    }
    let mut trajectory = writer(&output_dir.join("trajectory.jsonl"))?;
    let mut attempts = writer(&output_dir.join("attempts.jsonl"))?;
    let mut stress_checks = writer(&output_dir.join("stress-checks.jsonl"))?;
    let mut candidate_count = 0usize;
    let status = loop {
        let iteration = optimizer.checkpoint().next_iteration();
        if optimizer.checkpoint().is_complete() { break "iteration_limit"; }
        if options.pause_after.is_some_and(|count| iteration - start_iteration >= count) {
            break "paused";
        }
        let update = match optimizer.advance_one_controlled(|_| {
            if started.elapsed() >= wall_budget { ControlFlow::Break(()) }
            else { ControlFlow::Continue(()) }
        })? {
            ControlFlow::Continue(update) => update,
            ControlFlow::Break(()) => break "wall_budget",
        };
        for check in &update.stress_checks {
            let state = check.evaluation.as_ref().map_or_else(|| "null".to_string(), stress_json);
            let refusal = check.refusal.as_deref().map_or_else(|| "null".to_string(), json_string);
            writeln!(stress_checks,
                "{{\"iter\":{iteration},\"candidate\":{},\"state\":{state},\"refusal\":{refusal}}}", check.index,
            )?;
        }
        match update.progress {
            ProjectedProgress::Accepted(step) => {
                candidate_count += step.attempts.len();
                attempt_rows(&mut attempts, iteration, &step.attempts)?;
                let state = stress_json(optimizer.current());
                let accepted_field = format!("accepted-{:04}.csv", iteration + 1);
                field(&output_dir.join(&accepted_field), optimizer.checkpoint().geometry())?;
                if let Some(executable) = executable.as_deref() {
                    let name = format!("checkpoint-{:04}.fscp", iteration + 1);
                    checkpoint::save(&output_dir.join(&name), &optimizer, &policy, executable)?;
                    last_checkpoint = Some(name);
                }
                writeln!(trajectory,
                    "{{\"iter\":{iteration},\"previous_compliance\":{:.17e},\"state\":{state},\"field\":{},\"ell\":{:.17e},\"authority\":\"estimated\"}}",
                    step.previous.compliance, json_string(&accepted_field), optimizer.checkpoint().ell(),
                )?;
                // Every accepted field and checkpoint is flushed before more
                // physics. A partial checkpoint fails its length/digest on load.
                trajectory.flush()?;
                trajectory.get_ref().sync_all()?;
                attempts.flush()?;
                stress_checks.flush()?;
            }
            ProjectedProgress::IterationLimit => break "iteration_limit",
            ProjectedProgress::NoDescent(rows) => {
                candidate_count += rows.len();
                attempt_rows(&mut attempts, iteration, &rows)?;
                break "no_feasible_descent";
            }
        }
    };
    trajectory.flush()?;
    attempts.flush()?;
    stress_checks.flush()?;
    field(&output_dir.join("level-set.csv"), optimizer.checkpoint().geometry())?;
    let current = optimizer.current();
    let reduction = if baseline.compliance > 0.0 {
        (baseline.compliance - current.compliance) / baseline.compliance
    } else { 0.0 };
    let checkpoint = last_checkpoint.as_deref().map_or_else(|| "null".to_string(), json_string);
    let input_field = if source.is_some() { "\"input-level-set.csv\"" } else { "null" };
    let summary = format!(concat!(
        "{{\"schema\":\"fs-marquee-projected-stress-v2\",",
        "\"model\":\"normalized_unit_square_plane_strain_cantilever\",\"authority\":\"estimated\",",
        "\"status\":\"{}\",\"level\":{},\"requested_updates\":{},\"accepted_updates\":{},",
        "\"candidate_count\":{},\"candidate_budget_per_update\":{},\"wall_budget_seconds\":{},",
        "\"volume_target\":{:.17e},\"volume_tolerance\":{:.17e},",
        "\"sampled_von_mises_limit\":{:.17e},\"sampled_von_mises_tolerance\":{:.17e},",
        "\"baseline\":{},\"current\":{},\"relative_reduction_from_feasible_baseline\":{:.17e},",
        "\"resumed\":{},\"start_iteration\":{},\"segment_accepted_updates\":{},",
        "\"baseline_scope\":\"current_segment\",\"checkpoint\":{},",
        "\"fixed_boundaries\":[\"left\",\"right\"],",
        "\"fixed_node_count\":{},\"design_regions\":{},",
        "\"material\":{{\"youngs\":{:.17e},\"poisson\":{:.17e}}},",
        "\"load\":{{\"traction_y\":{:.17e},\"band_half_width\":{:.17e}}},\"initial_field\":{},",
        "\"claims\":{{\"converged\":false,\"global_optimum\":false,\"physical_validation\":false,",
        "\"stress_adjoint_kkt\":false,\"continuous_max_stress\":false,\"certified_continuum_volume\":false}}}}"
    ), status, settings.level, settings.iterations, optimizer.checkpoint().next_iteration(), candidate_count,
        policy.controls.max_candidates, wall_seconds, settings.volfrac, policy.projection.tolerance,
        limit.max_von_mises, limit.absolute_tolerance, stress_json(&baseline), stress_json(current), reduction,
        resumed, start_iteration, optimizer.checkpoint().next_iteration() - start_iteration, checkpoint,
        policy.fixed.len(), design_regions,
        settings.youngs, settings.poisson, -fixture.load, fixture.band, input_field,
    );
    let mut completion = writer(&output_dir.join("summary.json"))?;
    writeln!(completion, "{summary}")?;
    completion.flush()?;
    completion.get_ref().sync_all()?;
    println!("{summary}");
    Ok(match status { "iteration_limit" => 0, "wall_budget" | "paused" => 6, _ => 11 })
}

#[cfg(test)]
#[path = "projected/cancellation_tests.rs"]
mod cancellation_tests;
