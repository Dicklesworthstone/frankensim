use fs_uq::{CorrelationModel, ImportanceConfig, ImportanceExecution, ImportanceFitConfig,
    ParameterUncertainty, PilotImportanceConfig, PilotImportanceExecution, PilotImportanceStage,
    PropagationMethod, UqPlan, UqStatus};

fn plan(n: usize, threshold: f64) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, n)
        .with_parameter(ParameterUncertainty::gaussian("load", 0.0, 1.0, "N"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(threshold);
    plan.seed = 83;
    plan
}
fn config(rounds: usize, samples: usize) -> PilotImportanceConfig {
    PilotImportanceConfig {
        pilot_rounds: rounds, samples_per_pilot: samples,
        initial_proposal: ImportanceConfig { germ_shift: vec![0.0], defensive_fraction: 0.2 },
        fit: ImportanceFitConfig { elite_fraction: 0.1, relaxation: 1.0, maximum_shift_norm: 12.0 },
    }
}
fn response(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0]) }

#[test]
fn learns_a_five_sigma_proposal_and_reproduces_independent_production_exactly() {
    let plan = plan(4*512+8192, 5.0);
    let mut study = PilotImportanceExecution::new(&plan, config(4,512)).unwrap();
    let report = study.advance(usize::MAX, || false, response);
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.stage, PilotImportanceStage::Complete);
    assert_eq!(report.model_calls_evaluated, plan.budget_max_samples);
    assert_eq!(report.completed_pilots.len(), 4);
    assert_eq!(study.plan(), &plan);
    let production = study.production().unwrap();
    assert_eq!(production.observations().len(), 8192);
    assert!(production.config().germ_shift[0] > 4.5);
    let target = 2.866515718791933e-7;
    let estimate = report.estimate.unwrap();
    assert!((estimate.failure_probability / target - 1.0).abs() < 0.2);
    assert!(report.production.as_ref().unwrap().proposal_failures > 2000);
    // Every round can be reconstructed from its own immutable physical law,
    // key and sampled proposal; earlier pilot responses are not pooled.
    for round in study.pilot_history() {
        let mut p = plan.clone(); p.seed = round.seed; p.budget_max_samples = 512;
        let mut reference = ImportanceExecution::new(&p, round.sampled_proposal.clone()).unwrap();
        reference.advance(512, || false, response);
        assert_eq!(reference.fit_proposal(study.config().fit, || false).unwrap(), round.fit);
        assert_ne!(round.seed, production.plan().seed);
    }
    let mut reference = ImportanceExecution::new(production.plan(), production.config().clone()).unwrap();
    let reference_report = reference.advance(8192, || false, response);
    assert_eq!(reference_report, production.report());
    assert_eq!(reference.observations(), production.observations());
    assert_eq!(report.estimate, reference_report.estimate);
}

#[test]
fn every_call_fit_transition_and_publication_checkpoint_resumes_without_repaying_physics() {
    let plan = plan(32, 2.0);
    let cfg = config(2,8);
    let mut full = PilotImportanceExecution::new(&plan, cfg.clone()).unwrap();
    let mut expected_inputs = Vec::new();
    let mut polls = 0;
    let expected = full.advance(usize::MAX, || { polls += 1; false }, |x| {
        expected_inputs.push(x.to_vec()); response(x)
    });
    assert_eq!(expected.status, UqStatus::Complete);
    for boundary in 1..=polls {
        let mut run = PilotImportanceExecution::new(&plan, cfg.clone()).unwrap();
        let mut visited = 0;
        let mut actual_inputs = Vec::new();
        let paused = run.advance(usize::MAX, || { visited += 1; visited == boundary }, |x| {
            actual_inputs.push(x.to_vec()); response(x)
        });
        assert_eq!(paused.status, UqStatus::Cancelled, "checkpoint {boundary}");
        assert_eq!(visited, boundary);
        assert_eq!(paused.model_calls_evaluated, actual_inputs.len());
        let mut resumed = run.clone();
        let result = resumed.advance(usize::MAX, || false, |x| {
            actual_inputs.push(x.to_vec()); response(x)
        });
        assert_eq!(result, expected, "checkpoint {boundary}");
        assert_eq!(actual_inputs, expected_inputs, "checkpoint {boundary}");
    }
}

#[test]
fn one_call_chunks_and_unfinished_callbacks_keep_the_same_phase_design() {
    let plan = plan(32, 1.0);
    let cfg = config(2,8);
    let mut full = PilotImportanceExecution::new(&plan, cfg.clone()).unwrap();
    let mut inputs = Vec::new();
    let expected = full.advance(32, || false, |x| { inputs.push(x.to_vec()); response(x) });
    let mut run = PilotImportanceExecution::new(&plan, cfg).unwrap();
    assert_eq!(run.advance(0, || panic!("zero"), |_| -> Result<f64,&str> { panic!("zero") }), run.report());
    for (index, expected_input) in inputs.iter().enumerate() {
        let paused = run.advance_interruptible(1, || false, |x| {
            assert_eq!(x, expected_input.as_slice()); Ok::<_, &str>(None)
        });
        assert_eq!(paused.model_calls_evaluated, index);
        assert!(paused.estimate.is_none());
        let completed = run.advance(1, || false, |x| {
            assert_eq!(x, expected_input.as_slice()); response(x)
        });
        assert_eq!(completed.model_calls_evaluated, index+1);
    }
    assert_eq!(run.report(), expected);
    assert_eq!(run.advance(1, || panic!("terminal"), |_| -> Result<f64,&str> { panic!("terminal") }), expected);
}

#[test]
fn physical_failures_are_terminal_in_both_pilot_and_production() {
    for paid in [3, 16, 20] {
        for bad in [Err("physical model refused"), Ok(f64::NAN)] {
            let mut run = PilotImportanceExecution::new(&plan(32, 2.0), config(2,8)).unwrap();
            run.advance(paid, || false, response);
            let refused = run.advance(1, || false, |_| bad);
            assert_eq!(refused.status, UqStatus::Refused);
            assert_eq!(refused.model_calls_evaluated, paid+1);
            assert!(refused.estimate.is_none());
            assert!(refused.rejection_reason.is_some());
            assert_eq!(run.advance(32, || panic!("terminal"), |_| -> Result<f64,&str> { panic!("terminal") }), refused);
        }
    }
}

#[test]
fn singular_laws_flat_responses_and_wrapping_keys_do_not_fabricate_risk_certainty() {
    let mut p = plan(64, 1.0); p.seed = u64::MAX;
    p.parameters.push(ParameterUncertainty::gaussian("copy", 0.0, 1.0, "N"));
    p.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0,1.0];2] };
    let mut cfg = config(2,16); cfg.initial_proposal.germ_shift = vec![0.0;2];
    let mut study = PilotImportanceExecution::new(&p, cfg).unwrap();
    let report = study.advance(64, || false, |x| {
        assert_eq!(x[0].to_bits(), x[1].to_bits()); Ok::<_, &str>(1.0)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.estimate.unwrap().failure_probability, 0.0);
    assert!(report.estimate.unwrap().failure_effective_samples.is_none());
    let mut keys: Vec<_> = study.pilot_history().iter().map(|r| r.seed).collect();
    keys.push(report.current_seed); keys.sort_unstable(); keys.dedup(); assert_eq!(keys.len(), 3);
    let confidence = study.production().unwrap().assess_failure_probability(0.05, 1e-9).unwrap().unwrap();
    assert!(!confidence.converged && confidence.hi > 0.0);
}

#[test]
fn invalid_counts_laws_and_fit_policies_refuse_before_physical_work() {
    let p = plan(32, 2.0);
    for (rounds,samples) in [(0,8),(17,2),(2,1),(2,16),(usize::MAX,usize::MAX)] {
        assert!(PilotImportanceExecution::new(&p, config(rounds,samples)).is_err());
    }
    assert!(PilotImportanceExecution::new(&plan(1_000_001,2.0), config(2,8)).is_err());
    let mut cfg = config(2,8); cfg.fit.relaxation = 0.0;
    assert!(PilotImportanceExecution::new(&p, cfg).is_err());
    let mut wrong = p.clone(); wrong.parameters[0] = ParameterUncertainty::uniform("load",0.0,1.0,"N");
    assert!(PilotImportanceExecution::new(&wrong, config(2,8)).is_err());
    wrong = p; wrong.method = PropagationMethod::PolynomialChaos;
    assert!(PilotImportanceExecution::new(&wrong, config(2,8)).is_err());
}
