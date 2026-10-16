//! Real projected topology updates with local analysis refinement before admission.
//! Run with --strict to deliberately exhaust the endpoint-resolution allowance.
use std::ops::ControlFlow;
use fs_topols::{DesignBoxEdge, GridSdf, OptimizeSettings, RobustAggregate,
    RobustLoadCase, SampledStressLimit};
use fs_topols::adaptive_dwr::{AdaptiveDwrPolicy, gain::AdaptiveGainPolicy};
use fs_topols::robust_descent::{MultiLoadProjectedOptimizer,
    MultiLoadProjectedProgress, MultiLoadProjectedSettings};
use fs_topols::volume::VolumeProjectionSettings;

fn run(strict: bool) -> Result<(), fs_cutfem::CutFemError> {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = geometry.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    let loads = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0)?,
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0)?,
    ];
    let settings = OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6,
        move_cells: 0.1, nucleation_period: 0, ..OptimizeSettings::default() };
    let policy = AdaptiveGainPolicy { refinement: AdaptiveDwrPolicy {
        max_level: 5, max_rounds: 2, max_cells: 512,
        absolute_tolerance: if strict { 1e-30 } else { 1e6 },
        relative_tolerance: 0.0, mark_fraction: 0.25 },
        max_error_fraction: 0.99, area_tolerance: 0.05 };
    let mut owner = MultiLoadProjectedOptimizer::new(geometry, &loads, settings,
        RobustAggregate::WeightedSum, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 3, max_solves: 256,
            ..MultiLoadProjectedSettings::default() })?
        .with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0)?)?;
    println!("authority=Estimated; adaptive gain screens are not continuum or stress certificates");
    println!("initial_same_area_compliance_j={:.17e}", owner.current().objective);
    loop {
        let ControlFlow::Continue(update) = owner.advance_one_adaptive_dwr_polling(policy, 1,
            |_| ControlFlow::<()>::Continue(()))? else { unreachable!("no cancellation requested") };
        for check in &update.candidates {
            if let Some(report) = &check.report {
                for (round, row) in report.rows.iter().enumerate() {
                    println!("candidate={} round={} cells={} gain_j={:.8e} enriched_gain_j={:.8e} correction_change_j={:.8e} gain_change_j={:.8e} allowance_j={:.8e} areas_match={} admitted={}",
                        check.index, round, row.baseline.cells, row.coarse_gain, row.enriched_gain,
                        row.correction_change, row.observed_gain_change, row.gain_allowance,
                        row.areas_match, row.admissible);
                }
                println!("comparison_stop={} comparison_solves={}", report.stop.as_str(), report.solves);
            } else { println!("comparison_refusal={:?}", check.refusal); }
        }
        println!("accepted_updates={} lifetime_solve_debit={} current_compliance_j={:.17e}",
            owner.next_iteration(), owner.solves_started(), owner.current().objective);
        match update.progress {
            MultiLoadProjectedProgress::Accepted(_) => continue,
            MultiLoadProjectedProgress::IterationLimit => println!("stop=iteration-limit-not-optimality"),
            MultiLoadProjectedProgress::NoDescent(_) => println!("stop=no-admissible-candidate"),
            MultiLoadProjectedProgress::SolveBudget(_) => println!("stop=solve-budget"),
        }
        break;
    }
    Ok(())
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let strict = match args.as_slice() {
        [] => false,
        [arg] if arg == "--strict" => true,
        _ => { eprintln!("usage: adaptive_descent [--strict]"); std::process::exit(2); }
    };
    if let Err(error) = run(strict) { eprintln!("adaptive descent refused: {error}"); std::process::exit(1); }
}
