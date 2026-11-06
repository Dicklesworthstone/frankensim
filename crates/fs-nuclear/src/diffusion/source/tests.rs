//! G1 continuous-source/age laws and G3 cache, partition and zero-source checks.

use std::sync::Arc;

use crate::FissionError;
use crate::diffusion::*;
use crate::irradiation::NeutronGroups;

fn data() -> DiffusionData {
    DiffusionData {
        data_id: "manufactured-continuous-source".into(),
        groups: Arc::new(
            NeutronGroups::new("manufactured-source-grid".into(), vec![0.0, 1.0]).unwrap(),
        ),
        speeds_m_s: vec![2.0],
        cell_volumes_m3: vec![4.0],
        delayed_groups: vec![],
        reactions: vec![],
        links: vec![],
        leakage: vec![],
    }
}

fn reaction(sigma: f64, kind: DiffusionReactionKind) -> DiffusionReaction {
    DiffusionReaction {
        data_id: "manufactured-source-reaction".into(),
        cell: 0,
        incident_group: 0,
        macroscopic_cross_section_m_inv: sigma,
        kind,
    }
}

fn state(data: DiffusionData, neutrons: f64, precursors: Vec<f64>) -> DiffusionState {
    DiffusionState::new(
        Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
        vec![neutrons],
        precursors,
    )
    .unwrap()
}

fn budget() -> DiffusionBudget {
    DiffusionBudget {
        limits: DiffusionLimits::default(),
        population_tolerance: 1.0e-7,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn balances(step: &DiffusionStep) {
    close(step.neutron_balance_residual, 0.0);
    close(step.precursor_balance_residual, 0.0);
    for &residual in &step.population_balance_residuals {
        close(residual, 0.0);
    }
    assert!(step.population_roundoff <= budget().population_tolerance);
}

#[test]
fn continuous_source_drives_subcritical_multiplication_with_actual_source_age_fluence() {
    let mut supplied = data();
    supplied.reactions = vec![
        reaction(
            0.1,
            DiffusionReactionKind::Fission {
                prompt_yields: vec![2.0],
                delayed_yields: vec![],
            },
        ),
        reaction(0.3, DiffusionReactionKind::Capture),
    ];
    let initial = state(supplied, 0.0, vec![]);
    let duration = 3.0_f64;
    let source = 5.0;
    let net_loss = 0.4_f64;
    let fraction = -(-net_loss * duration).exp_m1();
    let integral = source / net_loss * (duration - fraction / net_loss);
    let (next, step) = initial
        .trial_with_sources(duration, &[source], &[], budget(), || Ok(()))
        .unwrap();
    close(next.neutron_counts()[0], source / net_loss * fraction);
    close(step.injected_neutrons, source * duration);
    close(step.injected_neutron_counts[0], source * duration);
    close(step.fissions, 0.2 * integral);
    close(step.captures, 0.6 * integral);
    close(step.prompt_neutrons, 2.0 * step.fissions);
    close(step.integrated_fluence_m2[0], 0.5 * integral);
    close(
        step.mean_fluxes[0].values_m2_s()[0],
        0.5 * integral / duration,
    );
    assert!(next.neutron_counts()[0] > source * duration * (-net_loss * duration).exp());
    assert_eq!(initial.neutron_counts(), &[0.0]);
    balances(&step);
    let (half, first) = initial
        .trial_with_sources(1.5, &[source], &[], budget(), || Ok(()))
        .unwrap();
    let (split, second) = half
        .trial_with_sources(1.5, &[source], &[], budget(), || Ok(()))
        .unwrap();
    close(split.neutron_counts()[0], next.neutron_counts()[0]);
    close(
        first.integrated_fluence_m2[0] + second.integrated_fluence_m2[0],
        step.integrated_fluence_m2[0],
    );
    close(first.fissions + second.fissions, step.fissions);
}

#[test]
fn uniformly_injected_precursors_can_decay_and_irradiate_in_the_same_interval() {
    let mut supplied = data();
    supplied.delayed_groups.push(DelayedGroup {
        data_id: "manufactured-source-family".into(),
        decay_constant_s_inv: 0.5,
        spectrum: vec![1.0],
    });
    let initial = state(supplied, 0.0, vec![0.0]);
    let (next, step) = initial
        .trial_with_sources(2.0, &[0.0], &[6.0], budget(), || Ok(()))
        .unwrap();
    let lost = -(-1.0_f64).exp_m1();
    let precursor_integral = 12.0 * (2.0 - 2.0 * lost);
    let neutron_integral = 12.0 - precursor_integral;
    close(next.precursor_counts()[0], 12.0 * lost);
    close(next.neutron_counts()[0], 12.0 * (1.0 - lost));
    close(step.injected_precursors, 12.0);
    assert_eq!(step.injected_precursor_counts, vec![12.0]);
    close(step.precursor_decays, 0.5 * precursor_integral);
    close(step.delayed_neutrons, step.precursor_decays);
    close(step.integrated_fluence_m2[0], 0.5 * neutron_integral);
    assert!(step.integrated_fluence_m2[0] > 0.0);
    balances(&step);
}

#[test]
fn free_injection_has_half_interval_mean_age_and_zero_duration_is_identity() {
    let initial = state(data(), 4.0, vec![]);
    let (next, step) = initial
        .trial_with_sources(3.0, &[6.0], &[], budget(), || Ok(()))
        .unwrap();
    close(next.neutron_counts()[0], 22.0);
    close(
        step.integrated_fluence_m2[0],
        0.5 * (4.0 * 3.0 + 6.0 * 9.0 / 2.0),
    );
    close(step.injected_neutrons, 18.0);
    let (identity, zero) = initial
        .trial_with_sources(0.0, &[6.0], &[], budget(), || Ok(()))
        .unwrap();
    assert_eq!(identity.neutron_counts(), initial.neutron_counts());
    assert_eq!(zero.injected_neutron_counts, vec![0.0]);
    assert_eq!(zero.integrated_fluence_m2, vec![0.0]);
    close(zero.mean_fluxes[0].values_m2_s()[0], 2.0);
    balances(&step);
    balances(&zero);
}

#[test]
fn source_rates_do_not_invalidate_preparation_and_shutdown_reuses_the_same_operator() {
    let mut supplied = data();
    supplied
        .reactions
        .push(reaction(0.25, DiffusionReactionKind::Capture));
    let initial = state(supplied.clone(), 2.0, vec![]);
    let (closed, _) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let (first, _) = closed
        .trial_with_sources(1.0, &[4.0], &[], budget(), || Ok(()))
        .unwrap();
    assert!(!Arc::ptr_eq(
        closed.prepared.as_ref().unwrap(),
        first.prepared.as_ref().unwrap()
    ));
    let (second, _) = first
        .trial_with_sources(1.0, &[1.0], &[], budget(), || Ok(()))
        .unwrap();
    assert!(Arc::ptr_eq(
        first.prepared.as_ref().unwrap(),
        second.prepared.as_ref().unwrap()
    ));
    close(
        second.neutron_counts()[0],
        first.neutron_counts()[0] * (-0.5_f64).exp() + 2.0 * -(-0.5_f64).exp_m1(),
    );
    let (off, off_step) = second.trial(1.0, budget(), || Ok(())).unwrap();
    assert!(Arc::ptr_eq(
        second.prepared.as_ref().unwrap(),
        off.prepared.as_ref().unwrap()
    ));
    let (explicit, explicit_step) = second
        .trial_with_sources(1.0, &[0.0], &[], budget(), || Ok(()))
        .unwrap();
    assert_eq!(off.neutron_counts(), explicit.neutron_counts());
    assert_eq!(off_step, explicit_step);
    supplied.reactions[0].macroscopic_cross_section_m_inv = 0.5;
    let changed = Arc::new(
        first
            .model()
            .rebuild_with_data(supplied, DiffusionLimits::default(), || Ok(()))
            .unwrap(),
    );
    let rebound = first.rebind_model(changed).unwrap();
    assert!(rebound.prepared.is_none());
    let (rebound, _) = rebound
        .trial_with_sources(1.0, &[1.0], &[], budget(), || Ok(()))
        .unwrap();
    close(
        rebound.neutron_counts()[0],
        first.neutron_counts()[0] * (-1.0_f64).exp() + -(-1.0_f64).exp_m1(),
    );
}

#[test]
fn source_shape_range_account_and_checkpoint_refusals_leave_all_histories_unchanged() {
    let initial = state(data(), 4.0, vec![]);
    for rates in [vec![], vec![1.0, 2.0], vec![-1.0], vec![f64::INFINITY]] {
        assert!(
            initial
                .trial_with_sources(1.0, &rates, &[], budget(), || Ok(()))
                .is_err()
        );
    }
    assert!(
        initial
            .trial_with_sources(0.0, &[1.0], &[0.0], budget(), || Ok(()))
            .is_err()
    );
    assert!(
        initial
            .trial_with_sources(2.0, &[f64::MAX], &[], budget(), || Ok(()))
            .is_err()
    );
    assert!(
        initial
            .trial_with_sources(1.0e-300, &[1.0e-300], &[], budget(), || Ok(()))
            .is_err()
    );
    let large = state(data(), 2.0_f64.powi(60), vec![]);
    let exact_budget = DiffusionBudget {
        population_tolerance: 0.0,
        ..budget()
    };
    assert!(
        large
            .trial_with_sources(1.0, &[1.0], &[], exact_budget, || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    initial
        .trial_with_sources(1.0, &[3.0], &[], budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut candidate = initial.clone();
        let mut poll = 0;
        let result = candidate.advance_with_sources(1.0, &[3.0], &[], budget(), || {
            poll += 1;
            if poll == stop {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert_eq!(
            result.unwrap_err(),
            DiffusionError::Fission(FissionError::Cancelled)
        );
        assert_eq!(candidate.neutron_counts(), initial.neutron_counts());
        assert_eq!(candidate.precursor_counts(), initial.precursor_counts());
        assert!(candidate.prepared.is_none());
    }
}
