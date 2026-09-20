//! Opt-in box-constrained route through the reusable EquilibriumStudy owner.
use super::{DesignControl, DesignEvaluation, EquilibriumDesign};
use fs_ascent::{EquilibriumStudy, SqpRunReport, SqpStop};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::DesignWork;
use fs_exec::CancelGate;

type Error = Box<dyn std::error::Error>;

struct FitOutcome {
    initial: f64,
    report: SqpRunReport,
    audited: DesignEvaluation,
    work: DesignWork,
}

fn fit(problem: &EquilibriumDesign, iterations: usize, evaluations: usize, gate: &CancelGate)
    -> Result<FitOutcome, Error>
{
    if !(2..=4096).contains(&evaluations) || iterations > 512 {
        return Err("SQP fit needs 2..=4096 evaluations and 0..=512 iterations".into());
    }
    let cases = evaluations.checked_mul(problem.load_cases().len())
        .ok_or("SQP fit case budget overflow")?;
    let mut control = DesignControl::new(evaluations, cases);
    let point = vec![0.0; problem.variables().len()];
    let cap = point.len().checked_mul(3).ok_or("SQP fit dimension overflow")?;
    let mut study = EquilibriumStudy::new(problem, &point, &mut control, cap, gate)?;
    let initial = study.accepted().value;
    // Preserve the original example's independently recomputed final report.
    // The optimizer's hard ceiling leaves one complete family for this audit.
    let report = study.run(1e-8, iterations, evaluations - 1, gate)?;
    let accepted = study.accepted().clone();
    let final_point = study.optimizer().point().to_vec();
    drop(study);
    let audited = problem.evaluate(&final_point, &mut control, gate)?;
    if audited != accepted {
        return Err("final physical re-solve disagrees with the accepted SQP evidence".into());
    }
    Ok(FitOutcome { initial, report, audited, work: control.work() })
}

pub(super) fn run(problem: &EquilibriumDesign, source: &str, iterations: usize,
    evaluations: usize, gate: &CancelGate) -> Result<(), Error>
{
    let result = fit(problem, iterations, evaluations, gate)?;
    let report = &result.report;
    let audited = &result.audited;
    let mut maximum_error = 0.0_f64;
    let mut adjoint_residual = 0.0_f64;
    for (case, experiment) in audited.cases.iter().zip(problem.load_cases()) {
        adjoint_residual = adjoint_residual.max(case.adjoint_relative_residual);
        for (prediction, target) in case.observations_m.iter().zip(&experiment.targets) {
            maximum_error = maximum_error.max((prediction - target.target_m).abs());
        }
    }
    let raw_gradient = audited.gradient.iter().map(|g| g.abs()).fold(0.0, f64::max);
    let predictions: Vec<_> = audited.cases.iter().map(|case| &case.observations_m).collect();
    let kkt = &report.solution.kkt;
    println!("{{\"scope\":\"local-box-constrained-static-inverse-design\",\"solver\":\"sqp\",\"data\":\"{source}\",\"stop\":\"{:?}\",\"converged\":{},\"kkt_within_tolerance\":{},\"initial_objective\":{:.17e},\"objective\":{:.17e},\"support_n_m\":{:.17e},\"contact_n_m2\":{:.17e},\"gap_m\":{:.17e},\"maximum_observation_error_m\":{maximum_error:.17e},\"adjoint_residual\":{adjoint_residual:.17e},\"raw_gradient_norm\":{raw_gradient:.17e},\"kkt_stationarity\":{:.17e},\"kkt_feasibility\":{:.17e},\"kkt_dual_feasibility\":{:.17e},\"kkt_complementarity\":{:.17e},\"bound_multipliers\":{:?},\"predictions_m\":{:?},\"iterations\":{},\"evaluations_including_audit\":{},\"case_solves\":{}}}",
        report.stop, report.stop == SqpStop::Converged, report.solution.converged,
        result.initial, audited.value, audited.physical_parameters[0], audited.physical_parameters[1],
        audited.physical_parameters[2], kkt.stationarity, kkt.feasibility, kkt.dual_feasibility,
        kkt.complementarity, report.solution.nu, predictions, report.solution.iters,
        result.work.evaluations, result.work.case_solves);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::problem;

    #[test]
    fn sqp_recovers_all_three_contact_parameters_with_a_fresh_physical_audit() {
        let gate = CancelGate::new();
        let problem = problem::build(&problem::synthetic_data(), 0.0002, &gate).unwrap();
        let result = fit(&problem, 256, 1024, &gate).unwrap();
        assert_eq!(result.report.stop, SqpStop::Converged, "{:?}", result.report);
        assert!(result.report.solution.kkt.within_tolerance(1e-8));
        assert!(result.audited.value < 1e-12);
        for (actual, expected) in result.audited.physical_parameters.iter().zip([600.0, 1.8e8, 0.0002]) {
            assert!((actual - expected).abs() < 1e-4 * expected);
        }
        assert_eq!(result.work.evaluations, result.report.solution.evals + 1);
        assert_eq!(result.audited.cases.len(), 3);
    }

    #[test]
    fn inconsistent_contact_data_reaches_a_bound_with_nonzero_raw_gradient() {
        let gate = CancelGate::new();
        // Disclosed independent quadratic targets require K=8e8; the admitted
        // family has K<=4e8. A constrained optimum need not have zero raw g.
        let data = problem::parse_data(include_bytes!("bounded-observations.txt")).unwrap();
        let problem = problem::build(&data, 0.0002, &gate).unwrap();
        let result = fit(&problem, 256, 1024, &gate).unwrap();
        assert_eq!(result.report.stop, SqpStop::Converged, "{:?}", result.report);
        assert!(result.report.solution.kkt.within_tolerance(1e-8));
        assert!((result.audited.physical_parameters[1] / 4e8 - 1.0).abs() < 1e-7);
        assert!(result.audited.gradient[1] < -1e-6);
        assert!(result.report.solution.nu[3] > 1e-6);
        assert!(result.audited.value > 1e-5); // the observations are not exactly fit
        assert!(result.audited.value < result.initial);
    }

    #[test]
    fn budget_stop_still_recomputes_every_case_without_claiming_convergence() {
        let gate = CancelGate::new();
        let problem = problem::build(&problem::synthetic_data(), 0.0002, &gate).unwrap();
        let result = fit(&problem, 128, 2, &gate).unwrap();
        assert_eq!(result.report.stop, SqpStop::EvaluationLimit);
        assert!(!result.report.solution.converged);
        assert_eq!(result.work, DesignWork { evaluations: 2, case_solves: 6 });
        assert_eq!(result.initial.to_bits(), result.audited.value.to_bits());
        assert_eq!(result.report.solution.iters, 0);
    }

    #[test]
    fn physical_refusal_is_not_laundered_into_a_finite_fit_report() {
        let gate = CancelGate::new();
        let mut data = problem::synthetic_data();
        data[1][0] = 1e9;
        let problem = problem::build(&data, 0.0002, &gate).unwrap();
        assert!(fit(&problem, 128, 256, &gate).is_err());
        let cancelled = CancelGate::new(); cancelled.request();
        assert!(fit(&problem, 128, 256, &cancelled).is_err());
    }
}
