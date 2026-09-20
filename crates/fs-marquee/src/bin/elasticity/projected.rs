//! Product path for feasible-baseline, hard-area-constrained elasticity descent.
use super::{json_string, parse, writer};
use fs_topols::projected::{
    ProjectedAttempt, ProjectedOptimizer, ProjectedProgress, ProjectedSettings,
};
use fs_topols::volume::VolumeProjectionSettings;
use fs_topols::{Cantilever, GridSdf, OptimizeSettings};
use std::error::Error;
use std::io::Write;
use std::path::Path;

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
    if args.is_empty() || args.len() > 5 {
        return Err("usage: fs-marquee-elasticity --projected OUTPUT_DIR [LEVEL=4] [ITERATIONS=30] [VOLFRAC=0.45] [MAX_CANDIDATES=6]".into());
    }
    let output_dir = Path::new(&args[0]);
    if output_dir.try_exists()? {
        return Err("output directory already exists; refusing to overwrite it".into());
    }
    let level: u32 = parse(args, 1, 4)?;
    let iterations: usize = parse(args, 2, 30)?;
    let volfrac: f64 = parse(args, 3, 0.45)?;
    let max_candidates: usize = parse(args, 4, 6)?;
    if !(2..=7).contains(&level) || !(1..=200).contains(&iterations)
        || !(volfrac.is_finite() && volfrac > 0.001 && volfrac < 0.999)
        || !(1..=16).contains(&max_candidates)
    {
        return Err("projected mode requires LEVEL in [2,7], ITERATIONS in [1,200], VOLFRAC in (0.001,0.999), and MAX_CANDIDATES in [1,16]".into());
    }
    let n = 1usize << level;
    let geometry = GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.42);
    // Hold both supported boundary traces fixed. No loaded-edge material can
    // disappear merely to reduce the external work or satisfy the area budget.
    let fixed: Vec<_> = geometry.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % (n + 1) == 0 || index % (n + 1) == n).collect();
    let settings = OptimizeSettings {
        level, volfrac, iterations, move_cells: 0.2,
        nucleation_period: 4, hole_radius_cells: 1.5, ..OptimizeSettings::default()
    };
    let projection = VolumeProjectionSettings {
        target: volfrac, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    };
    let controls = ProjectedSettings { max_candidates, ..ProjectedSettings::default() };
    let mut optimizer = ProjectedOptimizer::new(geometry,
        Cantilever { load: 1.0, band: 0.125 }, settings, fixed, projection, controls)?;
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
                // Do not buffer the entire expensive study until completion.
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
    let summary = format!(
        concat!(
            "{{\"schema\":\"fs-marquee-projected-v1\",\"model\":\"normalized_unit_square_plane_strain_cantilever\",",
            "\"authority\":\"estimated\",\"status\":\"{}\",\"level\":{},\"requested_updates\":{},",
            "\"accepted_updates\":{},\"candidate_count\":{},\"candidate_budget_per_update\":{},",
            "\"volume_target\":{:.17e},\"volume_tolerance\":{:.17e},",
            "\"baseline_compliance\":{:.17e},\"baseline_volume\":{:.17e},",
            "\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\",",
            "\"relative_reduction_from_feasible_baseline\":{:.17e},",
            "\"fixed_boundaries\":[\"left\",\"right\"],",
            "\"claims\":{{\"converged\":false,\"global_optimum\":false,\"physical_validation\":false,\"certified_continuum_volume\":false}}}}"
        ),
        status, level, iterations, optimizer.checkpoint().next_iteration(), candidate_count,
        max_candidates, volfrac, projection.tolerance, baseline.compliance, baseline.volume,
        current.compliance, current.volume, current.snapshot, reduction,
    );
    // Summary is the completion marker and is created only AFTER every field
    // and trace export succeeded. A failed write never prints a success result.
    let mut completion = writer(&output_dir.join("summary.json"))?;
    writeln!(completion, "{summary}")?;
    completion.flush()?;
    println!("{summary}");
    Ok(if status == "iteration_limit" { 0 } else { 11 })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projected_options_refuse_before_creating_output() {
        let path = std::env::temp_dir().join(format!("frankensim-projected-refusal-{}", std::process::id()));
        assert!(!path.exists(), "test output must not preexist");
        let args = vec![path.to_string_lossy().into_owned(), "40".into()];
        assert!(run(&args).is_err());
        assert!(!path.exists());
    }
}
