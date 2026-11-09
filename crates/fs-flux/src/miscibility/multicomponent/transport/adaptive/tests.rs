use super::*;
use crate::miscibility::multicomponent::advection::MulticomponentMiscibilityAdvectionFace;
use crate::miscibility::multicomponent::tests::{face, field, near, network, uniform_matrix};

fn fixture(diffusivity: f64) -> (MulticomponentMiscibilityField, MulticomponentMiscibilityAdvection) {
    let net = network(
        &[0.018, 0.046, 0.072], vec![0.0; 9], &[0.001; 2],
        vec![face(0, 1)], vec![uniform_matrix(3, diffusivity)],
    );
    let initial = field(&net, &[&[0.6, 0.3, 0.1], &[0.2, 0.3, 0.5]]);
    let flow = net.prepare_advection(
        vec![MulticomponentMiscibilityAdvectionFace {
            forward_m3_s: 0.001, reverse_m3_s: 0.001,
        }], 1e-12, || false,
    ).unwrap();
    (initial, flow)
}
fn time(step: f64) -> MulticomponentMiscibilityTimeControl {
    MulticomponentMiscibilityTimeControl {
        initial_step_s: step, minimum_step_s: 1e-9, maximum_step_s: step,
        absolute_mole_fraction_tolerance: 1e-7,
        relative_mole_fraction_tolerance: 1e-7,
        maximum_attempts: 20_000,
    }
}

#[test]
fn single_donor_trials_resolve_mixing_even_when_cfl_subcycles_would_alias() {
    let (initial, flow) = fixture(1e-20);
    let advection = MulticomponentMiscibilityAdvectionControls {
        maximum_courant: 0.25,
        ..MulticomponentMiscibilityAdvectionControls::default()
    };
    // Independent subcycling would make a 0.5-second coarse path and its two
    // half paths use exactly the same two 0.25-second donor steps.
    let candidate = initial.advance_transport_candidate(
        &flow, 0.5, advection, MulticomponentMiscibilitySolve::default(), time(0.5), || false,
    ).unwrap();
    let report = candidate.advance();
    assert!(report.rejected_intervals > 0);
    assert!(report.accepted_intervals > 2);
    assert!(report.maximum_accepted_error_ratio > 0.0);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert!(report.maximum_accepted_step_s < 0.25);
    assert!(report.maximum_advection_courant <= 0.125);
    assert_eq!(report.completed_advection_substeps, report.completed_transport_steps);
    assert_eq!(candidate.field().elapsed_s(), 0.5);
    // Independent compartment ODE: x_A,left=.4+.2*exp(-2Qt/V).
    // Molecular transport stays enabled but is negligible at this accuracy.
    let expected = 0.4 + 0.2 * (-1.0_f64).exp();
    near(candidate.field().mole_fractions()[0], expected, 8e-5);
    for species in 0..3 {
        near(report.component_residual_mol[species], 0.0, 1e-14);
    }
    assert!(report.free_energy_change_j < 0.0);
}

#[test]
fn coupled_adaptive_receipts_retain_only_accepted_advection_and_diffusion() {
    let (initial, flow) = fixture(1e-4);
    let advection = MulticomponentMiscibilityAdvectionControls::default();
    let controls = MulticomponentMiscibilitySolve::default();
    let duration = 0.08;
    let first = initial.transport_candidate(&flow, 0.5 * duration, advection, controls, || false).unwrap();
    let second = first.field().transport_candidate(&flow, 0.5 * duration, advection, controls, || false).unwrap();
    let mut accuracy = time(duration);
    accuracy.absolute_mole_fraction_tolerance = 1.0;
    accuracy.maximum_attempts = 1;
    let candidate = initial.advance_transport_candidate(
        &flow, duration, advection, controls, accuracy, || false,
    ).unwrap();
    let replay = initial.advance_transport_candidate(
        &flow, duration, advection, controls, accuracy, || false,
    ).unwrap();
    assert_eq!(candidate, replay);
    assert_eq!(candidate.field(), second.field());
    assert_eq!(candidate.advance().accepted_intervals, 1);
    assert_eq!(candidate.advance().completed_transport_steps, 3);
    assert_eq!(candidate.advance().completed_maxwell_stefan_evaluations, 3);
    for species in 0..3 {
        let a = first.step().advection.transferred_mol[0][species]
            + second.step().advection.transferred_mol[0][species];
        let d = first.step().diffusion.transferred_mol[0][species]
            + second.step().diffusion.transferred_mol[0][species];
        near(candidate.advance().advected_mol[0][species], a, 2e-17);
        near(candidate.advance().diffused_mol[0][species], d, 2e-17);
        near(candidate.advance().transferred_mol[0][species], a + d, 2e-17);
        near(initial.component_moles()[species] - candidate.field().component_moles()[species],
            candidate.advance().transferred_mol[0][species], 2e-16);
    }
    near(candidate.advance().stage_energy_residual_j, 0.0, 1e-11);
}

#[test]
fn incompatible_cfl_accuracy_work_and_cancellation_refuse_atomically() {
    let (initial, flow) = fixture(1e-4);
    let snapshot = initial.clone();
    let advection = MulticomponentMiscibilityAdvectionControls::default();
    let controls = MulticomponentMiscibilitySolve::default();
    let mut accuracy = time(1.0);
    accuracy.minimum_step_s = 1.0;
    assert!(matches!(
        initial.advance_transport_candidate(&flow, 1.0, advection, controls, accuracy, || false),
        Err(MulticomponentMiscibilityAdvanceError::Solver(MulticomponentMiscibilityError::Invalid(_)))
    ));
    accuracy = time(0.1);
    accuracy.absolute_mole_fraction_tolerance = 1e-16;
    accuracy.relative_mole_fraction_tolerance = 1e-16;
    accuracy.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_transport_candidate(&flow, 0.1, advection, controls, accuracy, || false),
        Err(MulticomponentMiscibilityAdvanceError::WorkLimit { .. })
    ));
    accuracy.minimum_step_s = 0.1;
    assert!(matches!(
        initial.advance_transport_candidate(&flow, 0.1, advection, controls, accuracy, || false),
        Err(MulticomponentMiscibilityAdvanceError::TemporalAccuracy { .. })
    ));
    let mut polls = 0usize;
    let result = initial.advance_transport_candidate(&flow, 0.1, advection, controls, time(0.1), || {
        polls += 1;
        polls > 80
    });
    assert!(matches!(result,
        Err(MulticomponentMiscibilityAdvanceError::Solver(MulticomponentMiscibilityError::Cancelled))
    ));
    assert_eq!(initial, snapshot);
}

#[test]
fn final_remainder_and_nonzero_start_time_preserve_the_requested_endpoint() {
    let (initial, flow) = fixture(1e-4);
    let controls = MulticomponentMiscibilitySolve::default();
    let mut advection = MulticomponentMiscibilityAdvectionControls::default();
    advection.maximum_courant = 0.9;
    let first = initial.transport_candidate(&flow, 0.02, advection, controls, || false).unwrap();
    let mut accuracy = time(0.06);
    accuracy.minimum_step_s = 0.05;
    accuracy.absolute_mole_fraction_tolerance = 1.0;
    let result = first.field().advance_transport_candidate(
        &flow, 0.07, advection, controls, accuracy, || false,
    ).unwrap();
    assert_eq!(result.field().elapsed_s(), first.field().elapsed_s() + 0.07);
    assert_eq!(result.advance().accepted_intervals, 2);
    assert!(result.advance().minimum_accepted_step_s < accuracy.minimum_step_s);
    near(result.advance().stage_energy_residual_j, 0.0, 1e-11);

    // A short final interval may lie below a minimum that exceeds the CFL cap.
    // The next-step receipt must nevertheless remain within that physical cap.
    let mut short = time(1.0);
    short.minimum_step_s = 0.6;
    short.absolute_mole_fraction_tolerance = 1.0;
    let result = initial.advance_transport_candidate(
        &flow, 0.1, MulticomponentMiscibilityAdvectionControls::default(),
        controls, short, || false,
    ).unwrap();
    assert!(result.advance().suggested_next_step_s <= 0.5);
    assert_eq!(result.field().elapsed_s(), 0.1);
}
