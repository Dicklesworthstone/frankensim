use fs_uq::{CorrelationModel, ImportanceConfig, ImportanceExecution, ParameterUncertainty,
    PropagationMethod, UqPlan, UqStatus};
use fs_uq::product_importance::proposal::{ImportanceFitConfig, ImportanceFitError};

fn plan(n: usize, threshold: f64) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, n)
        .with_parameter(ParameterUncertainty::gaussian("load", 0.0, 1.0, "N"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(threshold);
    plan.seed = 83;
    plan
}
fn policy() -> ImportanceFitConfig {
    ImportanceFitConfig { elite_fraction: 0.1, relaxation: 1.0, maximum_shift_norm: 12.0 }
}
fn run(n: usize, threshold: f64, shift: f64) -> ImportanceExecution {
    let mut run = ImportanceExecution::new(&plan(n, threshold), ImportanceConfig {
        germ_shift: vec![shift], defensive_fraction: 0.2,
    }).unwrap();
    assert_eq!(run.advance(n, || false, |x| Ok::<_, &str>(x[0])).status, UqStatus::Complete);
    run
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 1e-11 * (1.0+a.abs().max(b.abs())), "{a} {b}"); }

#[test]
fn shifted_pilot_fit_matches_independent_likelihood_weighted_reconstruction() {
    let run = run(1024, 1.0, 2.0);
    let before = run.report();
    let observations = run.observations().to_vec();
    let fit = run.fit_proposal(policy(), || false).unwrap();
    assert!(fit.targets_failure_event);
    let selected: Vec<_> = observations.iter().filter(|v| v.qoi > 1.0).collect();
    let weights: f64 = selected.iter().map(|v| v.likelihood_weight).sum();
    let expected = selected.iter().map(|v| v.likelihood_weight * v.qoi).sum::<f64>() / weights;
    close(fit.proposal.germ_shift[0], expected);
    close(fit.effective_samples, weights*weights /
        selected.iter().map(|v| v.likelihood_weight.powi(2)).sum::<f64>());
    let unweighted = selected.iter().map(|v| v.qoi).sum::<f64>() / selected.len() as f64;
    assert!((unweighted - expected).abs() > 0.1);
    assert_eq!(fit.selected, selected.len());
    assert_eq!(run.report(), before);
    assert_eq!(run.observations(), observations);
}

#[test]
fn no_failures_still_learns_toward_the_upper_tail_and_reports_smoothing_and_cap() {
    let run = run(256, 100.0, 0.0);
    let mut ranked: Vec<_> = run.observations().iter().map(|v| v.qoi).collect();
    ranked.sort_by(f64::total_cmp);
    let cutoff = ranked[230]; // ceil(0.1 * 256) = 26.
    let selected: Vec<_> = ranked.into_iter().filter(|&x| x >= cutoff).collect();
    let expected = selected.iter().sum::<f64>() / selected.len() as f64;
    let fit = run.fit_proposal(policy(), || false).unwrap();
    assert_eq!(fit.observed_failures, 0);
    assert!(!fit.targets_failure_event);
    assert_eq!(fit.cutoff, cutoff);
    close(fit.proposal.germ_shift[0], expected);
    assert!(expected > 1.0);
    let fit = run.fit_proposal(ImportanceFitConfig {
        relaxation: 0.5, maximum_shift_norm: 0.2, ..policy()
    }, || false).unwrap();
    close(fit.norm_before_limit, 0.5 * expected);
    close(fit.proposal.germ_shift[0], 0.2);
    assert!(fit.norm_limited);
    assert_eq!(fit.proposal.defensive_fraction, 0.2);
}

#[test]
fn ties_at_the_failure_ceiling_are_not_strict_failures() {
    let mut pilot = ImportanceExecution::new(&plan(64, 1.0), ImportanceConfig {
        germ_shift: vec![0.0], defensive_fraction: 0.2,
    }).unwrap();
    let mut germs = Vec::new();
    pilot.advance(64, || false, |x| {
        germs.push(x[0]);
        Ok::<_, &str>(1.0)
    });
    let fit = pilot.fit_proposal(policy(), || false).unwrap();
    assert_eq!(fit.selected, 64);
    assert_eq!(fit.observed_failures, 0);
    assert!(!fit.targets_failure_event);
    close(fit.proposal.germ_shift[0], germs.iter().sum::<f64>() / 64.0);
}

#[test]
fn singular_physical_maps_replay_the_hidden_germs_without_inversion() {
    let mut p = plan(256, 0.0);
    p.parameters.push(ParameterUncertainty::gaussian("copy", 0.0, 1.0, "N"));
    p.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0,1.0];2] };
    let mut pilot = ImportanceExecution::new(&p, ImportanceConfig {
        germ_shift: vec![1.0, 0.7], defensive_fraction: 0.2,
    }).unwrap();
    pilot.advance(256, || false, |x| {
        assert_eq!(x[0].to_bits(), x[1].to_bits());
        Ok::<_, &str>(x[0])
    });
    let fit = pilot.fit_proposal(policy(), || false).unwrap();
    assert_eq!(fit.proposal.germ_shift.len(), 2);
    assert!(fit.proposal.germ_shift.iter().all(|x| x.is_finite()));
    // Both original germs remain identifiable by sampler address, even though
    // the second coordinate cannot be reconstructed from either physical input.
    let mut next = p.clone(); next.seed += 1;
    assert!(ImportanceExecution::new(&next, fit.proposal).is_ok());
}

#[test]
fn incomplete_failed_and_invalid_pilots_refuse_without_mutation() {
    let mut pilot = ImportanceExecution::new(&plan(16, 1.0), ImportanceConfig {
        germ_shift: vec![0.0], defensive_fraction: 0.2,
    }).unwrap();
    pilot.advance(3, || false, |x| Ok::<_, &str>(x[0]));
    assert_eq!(pilot.fit_proposal(policy(), || false), Err(ImportanceFitError::IncompletePilot));
    pilot.advance(1, || false, |_| Err::<f64, _>("solver refused"));
    assert_eq!(pilot.fit_proposal(policy(), || false), Err(ImportanceFitError::IncompletePilot));
    let pilot = run(16, 1.0, 0.0);
    for invalid in [
        ImportanceFitConfig { elite_fraction: 0.0, ..policy() },
        ImportanceFitConfig { elite_fraction: f64::NAN, ..policy() },
        ImportanceFitConfig { relaxation: 1.1, ..policy() },
        ImportanceFitConfig { maximum_shift_norm: f64::INFINITY, ..policy() },
    ] {
        assert_eq!(pilot.fit_proposal(invalid, || false), Err(ImportanceFitError::InvalidConfig));
    }
}

#[test]
fn every_cancellation_boundary_preserves_the_completed_pilot_and_retries_identically() {
    let pilot = run(16, 2.0, 1.0);
    let before = pilot.report();
    let mut polls = 0;
    let expected = pilot.fit_proposal(policy(), || { polls += 1; false }).unwrap();
    for boundary in 1..=polls {
        let mut visited = 0;
        assert_eq!(pilot.fit_proposal(policy(), || { visited += 1; visited == boundary }),
            Err(ImportanceFitError::Cancelled));
        assert_eq!(visited, boundary);
    }
    assert_eq!(pilot.report(), before);
    assert_eq!(pilot.fit_proposal(policy(), || false).unwrap(), expected);
}
