use super::super::super::tests::{Case, near};
use super::*;
use crate::miscibility::multicomponent::adaptive::MulticomponentMiscibilityTimeControl;
use crate::miscibility::multicomponent::advection::MulticomponentMiscibilityAdvectionFace;

fn time(step: f64, temperature: f64) -> MulticomponentNonisothermalTimeControl {
    MulticomponentNonisothermalTimeControl {
        composition: MulticomponentMiscibilityTimeControl {
            initial_step_s: step,
            minimum_step_s: 1e-9,
            maximum_step_s: step,
            absolute_mole_fraction_tolerance: 1e-7,
            relative_mole_fraction_tolerance: 1e-6,
            maximum_attempts: 10000,
        },
        absolute_temperature_tolerance_k: temperature,
        relative_temperature_tolerance: 0.0,
    }
}
fn flow(field: &MulticomponentNonisothermalField, rate: f64) -> MulticomponentMiscibilityAdvection {
    field
        .network()
        .prepare_advection(
            vec![MulticomponentMiscibilityAdvectionFace {
                forward_m3_s: rate,
                reverse_m3_s: rate,
            }],
            1e-12,
            || false,
        )
        .unwrap()
}

#[test]
fn temperature_error_and_single_donor_cap_prevent_subcycling_alias() {
    let initial = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        volumes: vec![0.001; 2],
        ..Case::default()
    }
    .field();
    let flow = flow(&initial, 0.01);
    let advection = MulticomponentNonisothermalAdvectionControls::default();
    let controls = MulticomponentNonisothermalSolve::default();
    let candidate = initial
        .advance_transport_candidate(&flow, 0.2, advection, controls, time(1.0, 1e-3), || false)
        .unwrap();
    let contrast = 150.0 * (-4.0_f64).exp();
    near(
        candidate.field().temperatures()[0].value(),
        450.0 - contrast,
        0.08,
    );
    near(
        candidate.field().temperatures()[1].value(),
        450.0 + contrast,
        0.08,
    );
    assert!(candidate.advance().accepted_intervals > 4);
    assert!(candidate.advance().rejected_intervals > 0);
    assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
    assert!(candidate.advance().maximum_advection_courant <= advection.maximum_courant);
    assert!(candidate.advance().suggested_next_step_s <= 0.05);
    assert_eq!(
        candidate.field().component_moles(),
        initial.component_moles()
    );
    assert!(
        candidate
            .advance()
            .diffused_mol
            .iter()
            .flatten()
            .all(|&v| v == 0.0)
    );
    assert!(
        candidate
            .advance()
            .diffused_energy_j
            .iter()
            .all(|&v| v == 0.0)
    );
    assert!(candidate.advance().energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert!(candidate.advance().entropy_change_j_per_k > 0.0);
}

#[test]
fn one_coupled_interval_matches_two_full_kernels_and_accepted_receipts() {
    let initial = Case {
        conductivity: 2.0,
        heat: [400.0, -100.0, 0.0],
        ..Case::default()
    }
    .field();
    let flow = flow(&initial, 0.0002);
    let advection = MulticomponentNonisothermalAdvectionControls::default();
    let controls = MulticomponentNonisothermalSolve::default();
    let mut time = time(0.002, 100.0);
    time.composition.absolute_mole_fraction_tolerance = 1.0;
    let candidate = initial
        .advance_transport_candidate(&flow, 0.002, advection, controls, time, || false)
        .unwrap();
    let first = initial
        .transport_candidate(&flow, 0.001, advection, controls, || false)
        .unwrap();
    let second = first
        .field()
        .transport_candidate(&flow, 0.001, advection, controls, || false)
        .unwrap();
    assert_eq!(candidate.field(), second.field());
    assert_eq!(candidate.field().elapsed_s(), 0.002);
    assert_eq!(candidate.advance().attempted_intervals, 1);
    assert_eq!(candidate.advance().completed_transport_steps, 3);
    assert_eq!(candidate.advance().completed_advection_substeps, 3);
    assert_eq!(candidate.advance().completed_maxwell_stefan_evaluations, 3);
    assert!(candidate.advance().completed_caloric_iterations > 0);
    for s in 0..3 {
        near(
            candidate.advance().advected_mol[0][s],
            first.step().advection.transferred_mol[0][s]
                + second.step().advection.transferred_mol[0][s],
            1e-18,
        );
        near(
            candidate.advance().diffused_mol[0][s],
            first.step().diffusion.transferred_mol[0][s]
                + second.step().diffusion.transferred_mol[0][s],
            1e-18,
        );
        near(
            candidate.advance().transferred_mol[0][s],
            first.step().transferred_mol[0][s] + second.step().transferred_mol[0][s],
            1e-18,
        );
        near(
            candidate.field().component_moles()[s] - initial.component_moles()[s],
            -candidate.advance().transferred_mol[0][s],
            2e-15,
        );
    }
    near(
        candidate.advance().transferred_energy_j[0],
        candidate.advance().advected_energy_j[0] + candidate.advance().diffused_energy_j[0],
        1e-13,
    );
}

#[test]
fn rejected_coupled_trials_and_final_cancellation_leave_original_state_unchanged() {
    let initial = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        volumes: vec![0.001; 2],
        ..Case::default()
    }
    .field();
    let before = initial.clone();
    let flow = flow(&initial, 0.001);
    let advection = MulticomponentNonisothermalAdvectionControls::default();
    let controls = MulticomponentNonisothermalSolve::default();
    let mut limited = time(1.0, 1e-9);
    limited.composition.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_transport_candidate(&flow, 1.0, advection, controls, limited, || false),
        Err(
            MulticomponentNonisothermalTransportAdvanceError::WorkLimit {
                accepted_intervals: 0,
                ..
            }
        )
    ));
    assert_eq!(initial, before);
    let bounds = time(0.002, 0.001);
    let mut polls = 0usize;
    let expected = initial
        .advance_transport_candidate(&flow, 0.002, advection, controls, bounds, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0usize;
    assert!(matches!(
        initial.advance_transport_candidate(&flow, 0.002, advection, controls, bounds, || {
            calls += 1;
            calls == polls
        }),
        Err(MulticomponentNonisothermalTransportAdvanceError::Solver(
            MulticomponentNonisothermalError::Cancelled
        ))
    ));
    assert_eq!(initial, before);
    assert_eq!(
        initial
            .advance_transport_candidate(&flow, 0.002, advection, controls, bounds, || false)
            .unwrap(),
        expected
    );
}

#[test]
fn short_terminal_window_keeps_next_suggestion_inside_cfl_despite_larger_minimum() {
    let initial = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        temperatures: vec![300.0; 2],
        volumes: vec![0.001; 2],
        ..Case::default()
    }
    .field();
    let flow = flow(&initial, 0.01);
    let mut bounds = time(1.0, 1e-6);
    bounds.composition.minimum_step_s = 1.0;
    let candidate = initial
        .advance_transport_candidate(
            &flow,
            0.01,
            MulticomponentNonisothermalAdvectionControls::default(),
            MulticomponentNonisothermalSolve::default(),
            bounds,
            || false,
        )
        .unwrap();
    assert_eq!(candidate.field().elapsed_s(), 0.01);
    assert_eq!(
        candidate.field().component_moles(),
        initial.component_moles()
    );
    assert_eq!(candidate.field().temperatures(), initial.temperatures());
    assert_eq!(candidate.advance().accepted_intervals, 1);
    assert!(candidate.advance().suggested_next_step_s > 0.0);
    assert!(candidate.advance().suggested_next_step_s <= 0.05);
}
