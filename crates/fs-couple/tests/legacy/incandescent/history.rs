//! Held source transitions exercise the native enthalpy owner without a reset.

use super::{environment, fixture, input, source_fixture};
use fs_conduction::ConductionError;
use fs_couple::incandescent::history::{
    IncandescentDriveHistory, IncandescentHistoryAccuracy, IncandescentHistoryControls,
    IncandescentHistoryError, advance_history, advance_history_adaptive,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{Interp, TimeSignal};

const CURRENT: Dims = Dims([0, 0, 0, 0, 1, 0]);
const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);

fn controls(step: f64, maximum_steps: usize) -> IncandescentHistoryControls {
    IncandescentHistoryControls {
        start_s: 0.0,
        maximum_step_s: step,
        maximum_steps,
    }
}

#[test]
fn g1_drive_jump_ends_the_incoming_interval_then_cools_without_resetting_storage() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let drive = TimeSignal::Table {
        times: vec![0.0, 0.3],
        values: vec![2.0, 0.0],
        dims: CURRENT,
        interp: Interp::Hold,
    };
    let mut bath = environment();
    bath.convection_w_per_m2_k = 10.0; // hA = 0.1 W/K; C = 1 J/K.
    let original = input(&curve, 0.0, 1.0);
    let history = advance_history(
        &node,
        &bath,
        original,
        IncandescentDriveHistory::Current(&drive),
        controls(1.0, 10),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(history.len(), 2);
    assert_eq!((history[0].start_s, history[0].end_s), (0.0, 0.3));
    assert_eq!((history[1].start_s, history[1].end_s), (0.3, 1.0));
    assert_eq!(history[0].step.electrical.current_a, 2.0);
    assert_eq!(history[1].step.electrical.joule_power_w, 0.0);
    // Incoming source P=T/75 persists through the first endpoint. The next
    // interval retains that stored heat while the zero-drive resistor cools.
    let hot = (300.0 + 0.3 * 30.0) / (1.0 + 0.3 * (0.1 - 1.0 / 75.0));
    let cooled = (hot + 0.7 * 30.0) / (1.0 + 0.7 * 0.1);
    assert!((history[0].step.thermal.state.temperature_k() - hot).abs() < 1e-8);
    assert!((history[1].step.thermal.state.temperature_k() - cooled).abs() < 1e-8);
    assert!(cooled < hot && cooled > 300.0);
    let supplied: f64 = history
        .iter()
        .map(|sample| {
            sample.step.thermal.additional_internal_heat_j
                + sample.step.thermal.endpoint.convection_heat_j
                + sample.step.thermal.endpoint.radiation_heat_j
        })
        .sum();
    let final_h = history
        .last()
        .unwrap()
        .step
        .thermal
        .state
        .specific_enthalpy_j_kg();
    assert!((original.mass_kg * final_h - supplied).abs() < 1e-9);
    assert_eq!(original.initial.temperature_k(), 300.0);
}

#[test]
fn g1_prescribed_heat_is_spread_once_on_the_original_absolute_clock() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let zero = TimeSignal::Constant(QtyAny::new(0.0, VOLTAGE));
    let mut original = input(&curve, 0.0, 1.0);
    original.internal_heat_j = 10.0;
    let history = advance_history(
        &node,
        &environment(),
        original,
        IncandescentDriveHistory::Voltage(&zero),
        IncandescentHistoryControls {
            start_s: 2.0,
            ..controls(0.25, 10)
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(history.len(), 4);
    assert_eq!(history[0].start_s, 2.0);
    assert_eq!(history.last().unwrap().end_s, 3.0);
    for sample in &history {
        assert_eq!(sample.step.thermal.additional_internal_heat_j, 0.0);
        assert!((sample.step.thermal.endpoint.stored_energy_change_j - 2.5).abs() < 1e-9);
    }
    assert!((history.last().unwrap().step.thermal.state.temperature_k() - 310.0).abs() < 1e-8);
}

#[test]
fn g1_native_feedback_admission_and_temporal_resolution_are_both_active() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let drive = TimeSignal::Constant(QtyAny::new(2.0, CURRENT));
    let run = |duration, maximum_step| {
        advance_history(
            &node,
            &environment(),
            input(&curve, 0.0, duration),
            IncandescentDriveHistory::Current(&drive),
            controls(maximum_step, 100),
            || Ok(()),
        )
        .unwrap()
    };
    // The unbounded requested thermal step would violate current-control
    // feedback. Native admission splits it below 37.5 s without a root retry.
    let bounded = run(50.0, 50.0);
    assert_eq!(bounded.len(), 2);
    assert!(bounded.iter().all(|row| row.end_s - row.start_s < 37.5));
    assert!((bounded.last().unwrap().step.thermal.state.temperature_k() - 720.0).abs() < 1e-7);
    // Independent continuous-time solution for the manufactured C=1, P=T/75
    // law. A smaller declared step reduces the native first-order time error.
    let exact = 300.0 * fs_math::det::exp(10.0 / 75.0);
    let coarse = run(10.0, 1.0);
    let fine = run(10.0, 0.5);
    let coarse_error = coarse.last().unwrap().step.thermal.state.temperature_k() - exact;
    let fine_error = fine.last().unwrap().step.thermal.state.temperature_k() - exact;
    assert!(fine_error > 0.0 && fine_error < coarse_error);
    assert!((1.9..2.2).contains(&(coarse_error / fine_error)));
}

#[test]
fn g1_history_uses_the_original_phase_dependent_source_inside_a_latent_plateau() {
    let (node, curve) = source_fixture(ElectricalPropertyBasis::Resistivity, true);
    let drive = TimeSignal::Constant(QtyAny::new(2.0, CURRENT));
    let history = advance_history(
        &node,
        &environment(),
        input(&curve, 4000.0, 0.1),
        IncandescentDriveHistory::Current(&drive),
        controls(0.05, 10),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(history.len(), 2);
    // q(f)=1+2f, f=(h-3000)/2000, I=2 and mass=0.1 kg:
    // x=h-3000 follows x_next=(x_previous+40*dt)/(1-0.04*dt).
    let x1 = (1000.0 + 40.0 * 0.05) / (1.0 - 0.04 * 0.05);
    let x2 = (x1 + 40.0 * 0.05) / (1.0 - 0.04 * 0.05);
    let final_state = history.last().unwrap().step.thermal.state;
    assert_eq!(final_state.temperature_k(), 600.0);
    assert!((final_state.specific_enthalpy_j_kg() - (3000.0 + x2)).abs() < 1e-8);
    assert!(final_state.liquid_mass_fraction() > 0.5);
    assert!(history[1].step.electrical.joule_power_w > history[0].step.electrical.joule_power_w);
}

#[test]
fn g0_history_refusals_and_cancellation_leave_the_accepted_input_unchanged() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let drive = TimeSignal::Constant(QtyAny::new(2.0, CURRENT));
    let original = input(&curve, 0.0, 1.0);
    let budget = advance_history(
        &node,
        &environment(),
        original,
        IncandescentDriveHistory::Current(&drive),
        controls(0.25, 1),
        || Ok(()),
    );
    assert!(matches!(budget, Err(IncandescentHistoryError::Invalid(_))));
    let mut calls = 0;
    let cancelled = advance_history(
        &node,
        &environment(),
        original,
        IncandescentDriveHistory::Current(&drive),
        controls(0.01, 200),
        || {
            calls += 1;
            if calls == 200 {
                Err(ConductionError::Cancelled {
                    stage: "history-test",
                    at: calls,
                })
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(IncandescentHistoryError::Native(_))
    ));
    assert_eq!(calls, 200);
    assert_eq!(original.initial.temperature_k(), 300.0);
    let wrong_units = TimeSignal::Constant(QtyAny::new(2.0, VOLTAGE));
    assert!(
        advance_history(
            &node,
            &environment(),
            original,
            IncandescentDriveHistory::Current(&wrong_units),
            controls(0.1, 20),
            || Ok(())
        )
        .is_err()
    );
    let linear = TimeSignal::Table {
        times: vec![0.0, 1.0],
        values: vec![0.0, 2.0],
        dims: CURRENT,
        interp: Interp::Linear,
    };
    assert!(
        advance_history(
            &node,
            &environment(),
            original,
            IncandescentDriveHistory::Current(&linear),
            controls(0.1, 20),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        advance_history(
            &node,
            &environment(),
            original,
            IncandescentDriveHistory::Current(&drive),
            IncandescentHistoryControls {
                start_s: 1e20,
                ..controls(0.1, 20)
            },
            || Ok(())
        )
        .is_err()
    );
    // The endpoint rounds to start + 4, which would otherwise deposit 4/3 of
    // the prescribed heat when dividing it by the requested three seconds.
    let mut rounded = original;
    rounded.duration_s = 3.0;
    rounded.internal_heat_j = 9.0;
    assert!(matches!(
        advance_history(
            &node,
            &environment(),
            rounded,
            IncandescentDriveHistory::Current(&drive),
            IncandescentHistoryControls {
                start_s: 1e16,
                ..controls(3.0, 20)
            },
            || Ok(())
        ),
        Err(IncandescentHistoryError::Invalid(_))
    ));
}

#[test]
fn g1_adaptive_history_refines_real_time_error_and_preserves_native_work_at_a_jump() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let drive = TimeSignal::Table {
        times: vec![0.0, 10.0],
        values: vec![2.0, 0.0],
        dims: CURRENT,
        interp: Interp::Hold,
    };
    let bath = environment();
    let original = input(&curve, 0.0, 20.0);
    let run = |tolerance| {
        advance_history_adaptive(
            &node,
            &bath,
            original,
            IncandescentDriveHistory::Current(&drive),
            controls(20.0, 6000),
            IncandescentHistoryAccuracy {
                maximum_local_enthalpy_error_j_kg: tolerance,
            },
            || Ok(()),
        )
        .unwrap()
    };
    let loose = run(10.0);
    let fine = run(0.01);
    // C=1 J/K, I=2 A, R=T/300 ohm until 10 s; afterwards stored heat is held.
    let exact = 300.0 * fs_math::det::exp(10.0 / 75.0);
    let loose_error = loose
        .steps
        .last()
        .unwrap()
        .step
        .thermal
        .state
        .temperature_k()
        - exact;
    let fine_error = fine
        .steps
        .last()
        .unwrap()
        .step
        .thermal
        .state
        .temperature_k()
        - exact;
    assert!(fine_error > 0.0 && fine_error < loose_error && fine_error < 0.1);
    assert!(fine.native_solves > loose.native_solves);
    assert!(fine.rejected_intervals > 0);
    assert!(fine.maximum_accepted_error_j_kg <= 0.01);
    assert_eq!(
        fine.native_solves,
        3 * (fine.steps.len() / 2 + fine.rejected_intervals)
    );
    assert!(fine.native_solves <= 6000);
    assert_eq!(fine.steps.first().unwrap().start_s, 0.0);
    assert_eq!(fine.steps.last().unwrap().end_s, 20.0);
    assert!(
        fine.steps
            .windows(2)
            .all(|rows| rows[0].end_s == rows[1].start_s)
    );
    let incoming = fine.steps.iter().find(|row| row.end_s == 10.0).unwrap();
    assert_eq!(incoming.step.electrical.current_a, 2.0);
    assert!(
        fine.steps
            .iter()
            .filter(|row| row.start_s >= 10.0)
            .all(|row| row.step.electrical.joule_power_w == 0.0)
    );
    let heat: f64 = fine
        .steps
        .iter()
        .map(|row| {
            row.step.thermal.additional_internal_heat_j
                + row.step.thermal.endpoint.convection_heat_j
                + row.step.thermal.endpoint.radiation_heat_j
        })
        .sum();
    let stored = original.mass_kg
        * fine
            .steps
            .last()
            .unwrap()
            .step
            .thermal
            .state
            .specific_enthalpy_j_kg();
    let native_budget = bath
        .maximum_thermal_residual_j
        .min(original.mass_kg * bath.enthalpy_tolerance_j_kg);
    assert!((stored - heat).abs() <= fine.steps.len() as f64 * native_budget + 1e-10);

    // An odd number of clock ulps gives unequal represented half intervals.
    // Their native heat still sums to the prescribed 6 J; the coarse trial's
    // duplicate proposal must not enter the published energy account.
    let zero = TimeSignal::Constant(QtyAny::new(0.0, CURRENT));
    let mut prescribed = input(&curve, 0.0, 6.0);
    prescribed.internal_heat_j = 6.0;
    let rounded = advance_history_adaptive(
        &node,
        &bath,
        prescribed,
        IncandescentDriveHistory::Current(&zero),
        IncandescentHistoryControls {
            start_s: 1e16,
            ..controls(6.0, 3)
        },
        IncandescentHistoryAccuracy {
            maximum_local_enthalpy_error_j_kg: 0.01,
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(rounded.native_solves, 3);
    assert_eq!(rounded.steps.len(), 2);
    assert_ne!(
        rounded.steps[0].end_s - rounded.steps[0].start_s,
        rounded.steps[1].end_s - rounded.steps[1].start_s
    );
    let deposited: f64 = rounded
        .steps
        .iter()
        .map(|row| row.step.thermal.endpoint.stored_energy_change_j)
        .sum();
    assert!((deposited - 6.0).abs() < 1e-9);
}

#[test]
fn g1_adaptive_error_resolves_latent_phase_without_using_the_enthalpy_datum() {
    let (node, curve) = source_fixture(ElectricalPropertyBasis::Resistivity, true);
    let drive = TimeSignal::Constant(QtyAny::new(2.0, CURRENT));
    let offset = 1e6;
    let shifted = EquilibriumEnthalpyPhaseCurve::try_new(
        curve.material_card_identity(),
        curve
            .knots()
            .iter()
            .map(|knot| {
                let mut shifted = *knot;
                shifted.specific_enthalpy_j_kg += offset;
                shifted
            })
            .collect(),
    )
    .unwrap();
    let run = |chart: &EquilibriumEnthalpyPhaseCurve, initial_h| {
        advance_history_adaptive(
            &node,
            &environment(),
            input(chart, initial_h, 1.0),
            IncandescentDriveHistory::Current(&drive),
            controls(1.0, 1000),
            IncandescentHistoryAccuracy {
                maximum_local_enthalpy_error_j_kg: 0.01,
            },
            || Ok(()),
        )
        .unwrap()
    };
    let history = run(&curve, 4000.0);
    let translated = run(&shifted, 4000.0 + offset);
    assert!(history.rejected_intervals > 0 && history.steps.len() > 2);
    assert!(
        history
            .steps
            .iter()
            .all(|row| row.step.thermal.state.temperature_k() == 600.0)
    );
    assert!(history.maximum_accepted_error_j_kg <= 0.01);
    let state = history.steps.last().unwrap().step.thermal.state;
    // Within the plateau, dh/dt = 40 + .04(h-3000), with h(0)=4000.
    let exact = 2000.0 + 2000.0 * fs_math::det::exp(0.04);
    assert!((state.specific_enthalpy_j_kg() - exact).abs() < 0.1);
    assert!(state.liquid_mass_fraction() > 0.5);
    assert_eq!(translated.native_solves, history.native_solves);
    let translated_state = translated.steps.last().unwrap().step.thermal.state;
    assert!(
        (translated_state.specific_enthalpy_j_kg() - offset - state.specific_enthalpy_j_kg()).abs()
            < 1e-7
    );
    assert!((translated_state.liquid_mass_fraction() - state.liquid_mass_fraction()).abs() < 1e-9);
}

#[test]
fn g0_adaptive_trials_charge_all_work_and_refuse_unresolved_accuracy_atomically() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let drive = TimeSignal::Constant(QtyAny::new(2.0, CURRENT));
    let original = input(&curve, 0.0, 10.0);
    let accuracy = IncandescentHistoryAccuracy {
        maximum_local_enthalpy_error_j_kg: 0.01,
    };
    // Two solves cannot complete one trial. Three solves complete a rejected
    // trial but cannot publish its discarded coarse or half-step proposals.
    for budget in [2, 3] {
        assert!(matches!(
            advance_history_adaptive(
                &node,
                &environment(),
                original,
                IncandescentDriveHistory::Current(&drive),
                controls(10.0, budget),
                accuracy,
                || Ok(()),
            ),
            Err(IncandescentHistoryError::Invalid(_))
        ));
    }
    assert!(matches!(
        advance_history_adaptive(
            &node,
            &environment(),
            original,
            IncandescentDriveHistory::Current(&drive),
            controls(10.0, 1000),
            IncandescentHistoryAccuracy {
                maximum_local_enthalpy_error_j_kg: 1e-12
            },
            || Ok(()),
        ),
        Err(IncandescentHistoryError::Invalid(_))
    ));
    let mut clock_input = original;
    clock_input.duration_s = 2.0;
    assert!(matches!(
        advance_history_adaptive(
            &node,
            &environment(),
            clock_input,
            IncandescentDriveHistory::Current(&drive),
            IncandescentHistoryControls {
                start_s: 1e16,
                ..controls(2.0, 100)
            },
            accuracy,
            || Ok(()),
        ),
        Err(IncandescentHistoryError::Invalid(_))
    ));
    let mut calls = 0;
    let cancelled = advance_history_adaptive(
        &node,
        &environment(),
        original,
        IncandescentDriveHistory::Current(&drive),
        controls(10.0, 1000),
        accuracy,
        || {
            calls += 1;
            if calls == 200 {
                Err(ConductionError::Cancelled {
                    stage: "adaptive-history-test",
                    at: calls,
                })
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(IncandescentHistoryError::Native(_))
    ));
    assert_eq!(calls, 200);
    assert_eq!(original.initial.specific_enthalpy_j_kg(), 0.0);
    assert_eq!(original.initial.temperature_k(), 300.0);
}
