//! G0/G3/G4 checks on the native closed molecular entropy admission.
use super::*;
use crate::multicomponent::MAX_SPECIES;

fn closed_field(boundaries: Boundaries, formation: f64) -> MixtureEuler1d {
    field(gas(2, formation, false), 16, boundaries, |x| {
        let wave = (2.0 * core::f64::consts::PI * x).sin();
        let a = 0.5 + 0.2 * wave;
        (vec![a, 1.0 - a], 0.15 * wave, 2.0 + 0.2 * wave)
    })
}

// For these original NASA9 fixtures M_i=R, p0=1, and cp_i/R is constant.
// Thus S/V = sum_i rho_i [cv_i ln(T) - ln(rho_i)], including ideal mixing.
fn physical_entropy(flow: &MixtureEuler1d, cv: &[f64]) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let cell = flow.cell(i).unwrap();
            let t = flow.gas().primitive(cell, || false).unwrap().temperature;
            flow.mesh().volumes()[i]
                * cell
                    .partial_densities
                    .iter()
                    .zip(cv)
                    .filter(|(rho, _)| **rho > 0.0)
                    .map(|(rho, capacity)| rho * (capacity * t.ln() - rho.ln()))
                    .sum::<f64>()
        })
        .sum()
}

#[test]
fn native_closed_trajectory_accounts_for_heat_mixing_and_viscous_conversion() {
    for boundaries in [
        Boundaries::Periodic,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
    ] {
        let mut reference: Option<(Vec<f64>, MolecularEntropyReport)> = None;
        for formation in [0.0, 1e18] {
            let mut native = closed_field(boundaries, formation);
            let mut admitted = closed_field(boundaries, formation);
            let before = physical_entropy(&admitted, &[2.5, 2.5]);
            let mut ordinary = operator(&native, 0.03, 0.01, 0.05, 0.03);
            let mut audited = operator(&admitted, 0.03, 0.01, 0.05, 0.03);
            let controls = DissipationControl {
                maximum_step_s: 0.005,
                ..DissipationControl::default()
            };
            let native_report = ordinary
                .advance_for(&mut native, Time::new(0.03), controls, || false)
                .unwrap();
            let (report, entropy) = audited
                .advance_for_with_entropy(
                    &mut admitted,
                    Time::new(0.03),
                    controls,
                    MolecularEntropyControl::default(),
                    || false,
                )
                .unwrap();
            assert_eq!(admitted.states, native.states);
            assert_eq!(report, native_report);
            assert_eq!(admitted.time().value(), 0.0);
            assert_eq!(report.boundary_change.energy, 0.0);
            assert!(entropy.entropy_change_j_k > 0.0);
            assert!(entropy.minimum_step_entropy_change_j_k >= -1e-10);
            assert!(entropy.maximum_relative_species_residual <= 1e-12);
            assert_eq!(entropy.rejected_entropy_steps, 0);
            let expected = physical_entropy(&admitted, &[2.5, 2.5]) - before;
            assert!((entropy.source_entropy_change_j_k - expected).abs() < 3e-11);
            assert!(
                (entropy.source_entropy_change_j_k
                    - entropy.entropy_change_j_k
                    - entropy.reference_entropy_residual_j_k)
                    .abs()
                    < 2e-15
            );
            if let Some((states, account)) = &reference {
                assert_eq!(
                    &admitted.states, states,
                    "source formation datum changed dynamics"
                );
                assert_eq!(&entropy, account, "source formation datum changed entropy");
            } else {
                assert!(report.absolute_energy_residual_j.abs() < 2e-12);
                reference = Some((admitted.states.clone(), entropy));
            }
        }
    }
}

fn negative_rate_field() -> MixtureEuler1d {
    let gas = MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.0, 0.0),
            source("B", R, 6.5, 0.0, 0.0),
        ]),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    field(
        gas,
        2,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        |x| {
            if x < 0.5 {
                (vec![0.1, 0.9], 0.0, 1.0)
            } else {
                (vec![0.9, 0.1], 0.0, 10.0)
            }
        },
    )
}

#[test]
fn energy_conservation_does_not_hide_a_negative_molecular_entropy_rate() {
    let controls = DissipationControl {
        maximum_step_s: 1e-4,
        maximum_halvings: 3,
        ..DissipationControl::default()
    };
    let mut native = negative_rate_field();
    let initial_entropy = physical_entropy(&native, &[2.5, 5.5]);
    let mut ordinary = operator(&native, 0.0, 0.0, 0.0, 0.01);
    let report = ordinary
        .advance_for(&mut native, Time::new(1e-4), controls, || false)
        .unwrap();
    // Independent face calculation: j_A=-0.016, j_B=+0.016 kg/(m² s),
    // q_h=0.264 W/m², Tdot=[-0.0830769230769231,-0.1542857142857142] K/s.
    // Both cell volumes are 0.35 m³. Arithmetic face enthalpy is conservative
    // but its original thermal-plus-mixing entropy rate is negative.
    let measured_rate = (physical_entropy(&native, &[2.5, 5.5]) - initial_entropy) / 1e-4;
    assert!((measured_rate - (-0.039735310343068805)).abs() < 2e-6);
    assert!(report.absolute_energy_residual_j.abs() < 2e-13);

    let mut refused = negative_rate_field();
    let before = refused.states.clone();
    let mut op = operator(&refused, 0.0, 0.0, 0.0, 0.01);
    let error = op
        .advance_for_with_entropy(
            &mut refused,
            Time::new(1e-4),
            controls,
            MolecularEntropyControl {
                absolute_tolerance_j_k: 1e-14,
                ..MolecularEntropyControl::default()
            },
            || false,
        )
        .unwrap_err();
    match error {
        DissipationError::StepFailure { last, .. } => {
            assert!(
                matches!(*last, DissipationError::EntropyDecrease { change_j_k, .. }
                if change_j_k < -1e-14)
            );
        }
        other => panic!("expected a bounded entropy refusal, got {other:?}"),
    }
    assert_eq!(refused.states, before);
    assert_eq!(refused.time().value(), 0.0);
}

#[test]
fn complete_interval_uses_one_entropy_allowance_without_publishing_partial_work() {
    let controls = DissipationControl {
        maximum_step_s: 1e-4,
        ..DissipationControl::default()
    };
    let allowance = MolecularEntropyControl {
        absolute_tolerance_j_k: 5e-6,
        ..MolecularEntropyControl::default()
    };
    // Each approximately -3.97e-6 J/K substep fits the allowance independently.
    let mut split = negative_rate_field();
    let mut split_op = operator(&split, 0.0, 0.0, 0.0, 0.01);
    for _ in 0..2 {
        let (report, entropy) = split_op
            .advance_for_with_entropy(&mut split, Time::new(1e-4), controls, allowance, || false)
            .unwrap();
        assert_eq!(report.accepted_steps, 1);
        assert_eq!(report.rejected_steps, 0);
        assert!(entropy.entropy_change_j_k < 0.0);
        assert!(entropy.entropy_change_j_k >= -allowance.absolute_tolerance_j_k);
    }
    let mut whole = negative_rate_field();
    let before = whole.states.clone();
    let mut whole_op = operator(&whole, 0.0, 0.0, 0.0, 0.01);
    let error = whole_op
        .advance_for_with_entropy(&mut whole, Time::new(2e-4), controls, allowance, || false)
        .unwrap_err();
    // A direct final refusal, after the internal candidates passed: the
    // whole-duration allowance was not multiplied by accepted_steps.
    assert!(
        matches!(error, DissipationError::EntropyDecrease { change_j_k, .. }
        if change_j_k < -allowance.absolute_tolerance_j_k)
    );
    assert_eq!(whole.states, before);
    assert_eq!(whole.time().value(), 0.0);
}

#[test]
fn weak_mixing_and_absent_species_keep_their_actual_entropy_increment() {
    let g = gas(2, 0.0, false);
    let reference = [0.5; MAX_SPECIES];
    let change = |old: &[f64], new: &[f64]| {
        super::super::entropy::cell_change(&g, &reference, 2.0, old, new, 2.0, 2.0, &mut || false)
            .unwrap()
            .to_f64()
    };
    let epsilon = 2.0_f64.powi(-30);
    let actual = change(&[0.5 + epsilon, 0.5 - epsilon], &[0.5, 0.5]);
    // The exact binary mixing increment starts with 2 epsilon²; subtracting
    // independently rounded absolute entropy densities loses this signal.
    let expected = 2.0 * epsilon * epsilon;
    assert!(actual > 0.0);
    assert!((actual / expected - 1.0).abs() < 1e-9);
    let pure_to_mixed = change(&[1.0, 0.0], &[0.5, 0.5]) + change(&[0.0, 1.0], &[0.5, 0.5]);
    assert!((pure_to_mixed - 2.0 * core::f64::consts::LN_2).abs() < 2e-14);
}

#[test]
fn late_entropy_cancellation_is_atomic_and_workspaces_remain_reusable() {
    let mut reference = closed_field(Boundaries::Periodic, 0.0);
    let mut reference_op = operator(&reference, 0.03, 0.01, 0.05, 0.03);
    let controls = DissipationControl {
        maximum_step_s: 0.005,
        ..DissipationControl::default()
    };
    let mut polls = 0;
    reference_op
        .advance_for_with_entropy(
            &mut reference,
            Time::new(0.03),
            controls,
            MolecularEntropyControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = closed_field(Boundaries::Periodic, 0.0);
    let before = flow.states.clone();
    let mut op = operator(&flow, 0.03, 0.01, 0.05, 0.03);
    let mut visits = 0;
    let error = op
        .advance_for_with_entropy(
            &mut flow,
            Time::new(0.03),
            controls,
            MolecularEntropyControl::default(),
            || {
                visits += 1;
                visits >= polls
            },
        )
        .unwrap_err();
    assert_eq!(error, DissipationError::Cancelled);
    assert_eq!(flow.states, before);
    assert_eq!(flow.time().value(), 0.0);
    op.advance_for_with_entropy(
        &mut flow,
        Time::new(0.03),
        controls,
        MolecularEntropyControl::default(),
        || false,
    )
    .unwrap();
    for (actual, expected) in flow.states.iter().zip(&reference.states) {
        assert!((actual - expected).abs() < 1e-11);
    }
    let state = flow.states.clone();
    let (empty, entropy) = op
        .advance_for_with_entropy(
            &mut flow,
            Time::new(0.0),
            controls,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(empty.accepted_steps, 0);
    assert_eq!(entropy, MolecularEntropyReport::default());
    assert_eq!(flow.states, state);
}
