//! Physical inverse design through existing L-BFGS or opt-in bounded SQP.
//! This example fixes one supported-mass test rig, not an arbitrary assembly
//! importer. The reusable EquilibriumDesign evaluator accepts other networks.
mod problem;
#[cfg(feature = "equilibrium-design")]
mod bounded;
use fs_ascent::{LbfgsError, LbfgsReport, LbfgsState, StopReason, StopRule};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignControl, DesignError, DesignEvaluation, EquilibriumDesign,
};
use fs_exec::CancelGate;
use std::io::Read;

// No new optimization algorithm: these fields bind the immutable problem and
// cumulative physical work to the existing accepted-iterate checkpoint.
#[derive(Clone)]
struct FitSession<'a> {
    problem: &'a EquilibriumDesign,
    state: LbfgsState,
    control: DesignControl,
    maximum_evaluations: usize,
    rejected_domains: usize,
}
impl<'a> FitSession<'a> {
    fn new(problem: &'a EquilibriumDesign, maximum_evaluations: usize, gate: &CancelGate)
        -> Result<Self, LbfgsError<DesignError>>
    {
        if !(2..=4096).contains(&maximum_evaluations) {
            return Err(LbfgsError::InvalidInput("fit needs 2..=4096 total evaluations, including final re-solve"));
        }
        let maximum_cases = maximum_evaluations.checked_mul(problem.load_cases().len())
            .ok_or(LbfgsError::InvalidInput("fit case budget overflow"))?;
        let mut control = DesignControl::new(maximum_evaluations, maximum_cases);
        let state = LbfgsState::try_new(&vec![0.0; problem.variables().len()], 8, &mut |x| {
            problem.evaluate(x, &mut control, gate).map(|r| (r.value, r.gradient))
        })?;
        Ok(Self { problem, state, control, maximum_evaluations, rejected_domains: 0 })
    }
    fn run(&mut self, additional_iterations: usize, gate: &CancelGate)
        -> Result<LbfgsReport, LbfgsError<DesignError>>
    {
        if gate.is_requested() { return Err(LbfgsError::Evaluation(DesignError::Cancelled)); }
        let remaining = self.maximum_evaluations.saturating_sub(self.control.work().evaluations).saturating_sub(1);
        // Always reserve a complete final objective/gradient re-solve; cached
        // optimizer values alone are not the final physical design report.
        let ceiling = self.state.evals.checked_add(remaining)
            .ok_or(LbfgsError::InvalidInput("fit evaluation clock overflow"))?;
        let problem = self.problem;
        let control = &mut self.control;
        let rejected = &mut self.rejected_domains;
        let mut fg = |x: &[f64]| match problem.evaluate(x, control, gate) {
            Ok(r) => Ok((r.value, r.gradient)),
            Err(DesignError::OutsideBounds { .. }) => {
                *rejected += 1;
                // The existing Wolfe owner treats +infinity as an unavailable
                // domain point. This placeholder gradient is NEVER accepted.
                Ok((f64::INFINITY, vec![0.0; x.len()]))
            }
            Err(error) => Err(error), // physics, cancellation, budgets stay typed
        };
        let rule = StopRule::GradNorm(1e-8);
        let mut report = self.state.try_run(&mut fg, &rule, 0, ceiling)?;
        for _ in 0..additional_iterations {
            if report.reason != StopReason::IterationCap { break; }
            if gate.is_requested() { return Err(LbfgsError::Evaluation(DesignError::Cancelled)); }
            report = self.state.try_run(&mut fg, &rule, 1, ceiling)?;
        }
        Ok(report)
    }
    fn audit(&mut self, gate: &CancelGate) -> Result<DesignEvaluation, DesignError> {
        self.problem.evaluate(&self.state.x, &mut self.control, gate)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Solver { Lbfgs, Sqp }

#[derive(Debug, PartialEq)]
struct Options { data: Option<String>, iterations: usize, evaluations: usize, scale_m: f64, solver: Solver }
fn options(args: &[String]) -> Result<Options, String> {
    let mut result = Options { data: None, iterations: 128, evaluations: 256, scale_m: 0.0002, solver: Solver::Lbfgs };
    let mut seen = Vec::new();
    let mut args = args.iter();
    while let Some(flag) = args.next() {
        if seen.contains(flag) { return Err(format!("duplicate option {flag}")); }
        seen.push(flag.clone());
        let value = args.next().ok_or_else(|| format!("missing value for {flag}"))?;
        match flag.as_str() {
            "--data" => result.data = Some(value.clone()),
            "--solver" => result.solver = match value.as_str() {
                "lbfgs" => Solver::Lbfgs,
                "sqp" => Solver::Sqp,
                _ => return Err("--solver must be lbfgs or sqp".into()),
            },
            "--iterations" => result.iterations = value.parse::<usize>().ok().filter(|x| *x <= 512)
                .ok_or("--iterations must be in 0..=512")?,
            "--evaluations" => result.evaluations = value.parse::<usize>().ok().filter(|x| (2..=4096).contains(x))
                .ok_or("--evaluations must be in 2..=4096")?,
            "--scale-m" => result.scale_m = value.parse::<f64>().ok().filter(|x| x.is_finite() && *x > 0.0)
                .ok_or("--scale-m must be finite and positive")?,
            _ => return Err(format!("unknown option {flag}")),
        }
    }
    Ok(result)
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("contact_fit [--solver lbfgs|sqp] [--data observations.txt] [--iterations 128] [--evaluations 256] [--scale-m 0.0002]\nEach data row: load_N mass_displacement_m receiver_displacement_m. Without --data, use disclosed synthetic targets. SQP requires --features equilibrium-design; L-BFGS remains the default.");
        return Ok(());
    }
    let options = options(&args)?;
    #[cfg(not(feature = "equilibrium-design"))]
    if options.solver == Solver::Sqp {
        return Err("--solver sqp requires cargo run --features equilibrium-design".into());
    }
    let (data, source) = if let Some(path) = &options.data {
        let mut bytes = Vec::new();
        std::fs::File::open(path)?.take(8193).read_to_end(&mut bytes)?;
        (problem::parse_data(&bytes)?, "provided-displacements")
    } else { (problem::synthetic_data(), "synthetic-quadratic-oracle") };
    let gate = CancelGate::new();
    let problem = problem::build(&data, options.scale_m, &gate)?;
    #[cfg(feature = "equilibrium-design")]
    if options.solver == Solver::Sqp {
        return bounded::run(&problem, source, options.iterations, options.evaluations, &gate);
    }
    let mut fit = FitSession::new(&problem, options.evaluations, &gate)?;
    let initial = fit.state.f;
    let report = fit.run(options.iterations, &gate)?;
    let audited = fit.audit(&gate)?;
    let mut maximum_error = 0.0_f64;
    let mut adjoint_residual = 0.0_f64;
    for (case, experiment) in audited.cases.iter().zip(problem.load_cases()) {
        adjoint_residual = adjoint_residual.max(case.adjoint_relative_residual);
        for (prediction, target) in case.observations_m.iter().zip(&experiment.targets) {
            maximum_error = maximum_error.max((prediction-target.target_m).abs());
        }
    }
    let work = fit.control.work();
    println!("{{\"scope\":\"local-static-inverse-design\",\"data\":\"{source}\",\"stop\":\"{:?}\",\"initial_objective\":{initial:.17e},\"objective\":{:.17e},\"support_n_m\":{:.17e},\"contact_n_m2\":{:.17e},\"gap_m\":{:.17e},\"maximum_observation_error_m\":{maximum_error:.17e},\"adjoint_residual\":{adjoint_residual:.17e},\"gradient_norm\":{:.17e},\"iterations\":{},\"evaluations_including_audit\":{},\"case_solves\":{},\"rejected_domains\":{}}}",
        report.reason, audited.value, audited.physical_parameters[0], audited.physical_parameters[1],
        audited.physical_parameters[2], report.grad_norm, report.iters, work.evaluations, work.case_solves, fit.rejected_domains);
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("contact fit refused: {error}"); std::process::exit(1); }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rig() -> EquilibriumDesign { problem::build(&problem::synthetic_data(),0.0002,&CancelGate::new()).unwrap() }

    #[test]
    fn existing_lbfgs_recovers_three_shared_parameters_from_independent_experiments() {
        let p=rig();let gate=CancelGate::new();let mut fit=FitSession::new(&p,256,&gate).unwrap();
        let initial=fit.state.f;let report=fit.run(128,&gate).unwrap();let audit=fit.audit(&gate).unwrap();
        assert!(audit.value<1e-12 && audit.value<initial*1e-8,"{:?}: {}",report.reason,audit.value);
        for (actual,expected) in audit.physical_parameters.iter().zip([600.0,1.8e8,0.0002]) {
            assert!((actual-expected).abs()<1e-4*expected,"{actual} != {expected}");
        }
        assert!(fit.control.work().evaluations<=256);
        assert!(fit.control.work().case_solves<=768);
        assert!(audit.cases.iter().all(|c|c.adjoint_relative_residual<=1e-10));
    }

    #[test]
    fn accepted_iteration_checkpoints_and_cancellation_resume_the_same_fit() {
        let p=rig();let gate=CancelGate::new();let mut whole=FitSession::new(&p,256,&gate).unwrap();
        let mut split=whole.clone();whole.run(80,&gate).unwrap();split.run(7,&gate).unwrap();
        let snapshot=split.clone();let cancelled=CancelGate::new();cancelled.request();
        assert!(matches!(split.run(1,&cancelled),Err(LbfgsError::Evaluation(DesignError::Cancelled))));
        assert_eq!(split.state.x,snapshot.state.x);assert_eq!(split.state.history,snapshot.state.history);
        assert_eq!(split.control.work(),snapshot.control.work());
        split.run(73,&gate).unwrap();
        assert_eq!(whole.state.x,split.state.x);assert_eq!(whole.state.g,split.state.g);
        assert_eq!(whole.state.history,split.state.history);assert_eq!(whole.control.work(),split.control.work());
        assert_eq!(whole.rejected_domains,split.rejected_domains);
    }

    #[test]
    fn exhausted_budget_retains_a_finite_design_and_reserves_the_final_physical_audit() {
        let p=rig();let gate=CancelGate::new();let mut fit=FitSession::new(&p,2,&gate).unwrap();
        let before=fit.state.x.clone();let report=fit.run(128,&gate).unwrap();
        assert_eq!(report.reason,StopReason::Budget);assert_eq!(before,fit.state.x);
        let audit=fit.audit(&gate).unwrap();assert!(audit.value.is_finite());
        assert_eq!(fit.control.work().evaluations,2);assert_eq!(fit.control.work().case_solves,6);
        assert!(matches!(fit.audit(&gate),Err(DesignError::Budget {..})));
    }

    #[test]
    fn solver_selection_is_explicit_and_preserves_the_default() {
        assert_eq!(options(&[]).unwrap().solver, Solver::Lbfgs);
        let args = vec!["--solver".into(), "sqp".into()];
        assert_eq!(options(&args).unwrap().solver, Solver::Sqp);
        let repeated = vec!["--solver".into(), "sqp".into(), "--solver".into(), "lbfgs".into()];
        assert!(options(&repeated).is_err());
    }

    #[test]
    fn external_observations_and_options_require_complete_finite_units() {
        assert_eq!(problem::parse_data(b"0.8 0.0003 0.00006\n1.6 0.0004 0.00013\n").unwrap().len(),2);
        for input in [b"".as_slice(),b"1 2",b"1 2 3 4",b"1 NaN 3",b"1 2 3\n\n"] {
            assert!(problem::parse_data(input).is_err());
        }
        assert!(problem::parse_data(&vec![b' ';8193]).is_err());
        for values in [vec!["--evaluations","1"],vec!["--scale-m","NaN"],vec!["--unknown","x"],
            vec!["--data","a","--data","b"],vec!["--iterations"],vec!["--solver","unknown"]] {
            assert!(options(&values.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
    }
}
