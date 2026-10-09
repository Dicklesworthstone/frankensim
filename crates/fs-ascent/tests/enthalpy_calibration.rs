//! G1/G3: complete spatial latent-heat trajectory gradients and heater recovery.

#[allow(dead_code)]
#[path = "../examples/enthalpy_calibration.rs"]
mod calibration;

use calibration::{EnthalpyCalibration, START, TRUTH, fit};
use fs_ascent::SqpStop;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

#[test]
fn full_temperature_history_gradient_matches_forward_differences() {
    with_cx(|cx| {
        let experiment = EnthalpyCalibration::synthetic(cx).unwrap();
        let midpoint = std::array::from_fn(|j| START[j].midpoint(TRUTH[j]));
        for point in [START, midpoint] {
            let evaluation = experiment.evaluate(cx, &point).unwrap();
            let forward_loss = experiment.forward_loss(cx, &point).unwrap();
            assert!((evaluation.loss_k2 - forward_loss).abs() < 1e-12);
            for j in 0..2 {
                let (mut plus, mut minus) = (point, point);
                let delta = 5e-5;
                plus[j] += delta;
                minus[j] -= delta;
                let fd = (experiment.forward_loss(cx, &plus).unwrap()
                    - experiment.forward_loss(cx, &minus).unwrap())
                    / (2.0 * delta);
                let gradient = evaluation.gradient[j];
                assert!(
                    (gradient - fd).abs() < 2e-6 * (1.0 + fd.abs()),
                    "point={point:?}, pulse={j}, adjoint={gradient}, forward FD={fd}"
                );
            }
        }
    });
}

#[test]
fn bounded_sqp_recovers_both_pulses_through_a_spatial_phase_change() {
    with_cx(|cx| {
        let experiment = EnthalpyCalibration::synthetic(cx).unwrap();
        let initial = experiment.forward_loss(cx, &START).unwrap();
        let (state, report) = fit(cx, &experiment).unwrap();
        assert_eq!(report.stop, SqpStop::Converged, "{report:?}");
        assert!(state.sample().f < initial * 1e-8);
        for (&actual, expected) in state.point().iter().zip(TRUTH) {
            assert!(
                (actual - expected).abs() < 1e-3,
                "recovered {actual} W/m3, synthetic target {expected} W/m3"
            );
        }
        let final_state = experiment.evaluate(cx, state.point()).unwrap();
        assert!(final_state.final_liquid_fraction.iter().any(|f| *f > 0.99));
        assert!(final_state.final_liquid_fraction.iter().any(|f| *f < 0.01));
        assert!(
            final_state
                .final_liquid_fraction
                .iter()
                .any(|f| *f > 0.01 && *f < 0.99)
        );
        assert!(final_state.max_step_energy_residual_j < 1e-9);
        assert!((final_state.loss_k2 - state.sample().f).abs() < 1e-12);
        println!(
            "heater recovery {:?}, initial loss={initial:.9e}, final loss={:.9e}, iterations={}, evaluations={}, max energy defect={:.9e} J",
            state.point(),
            state.sample().f,
            state.iterations(),
            state.evaluations(),
            final_state.max_step_energy_residual_j
        );
    });
}

#[test]
fn radiative_heater_recovery_uses_complete_temperature_history_gradients() {
    with_cx(|cx| {
        let experiment = EnthalpyCalibration::synthetic_with_ambient_radiation(cx).unwrap();
        let insulated = EnthalpyCalibration::synthetic(cx).unwrap();
        let radiative_history = experiment.forward_temperatures(cx, &TRUTH).unwrap();
        let insulated_history = insulated.forward_temperatures(cx, &TRUTH).unwrap();
        let largest_change = radiative_history
            .iter()
            .flatten()
            .zip(insulated_history.iter().flatten())
            .map(|(radiative, insulated)| (radiative - insulated).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            largest_change > 0.01,
            "the surface law must affect the solved field"
        );

        let midpoint = std::array::from_fn(|j| START[j].midpoint(TRUTH[j]));
        for point in [START, midpoint] {
            let evaluation = experiment.evaluate(cx, &point).unwrap();
            let forward_loss = experiment.forward_loss(cx, &point).unwrap();
            assert!((evaluation.loss_k2 - forward_loss).abs() < 1e-11);
            for pulse in 0..2 {
                let (mut plus, mut minus) = (point, point);
                let delta = 5e-5;
                plus[pulse] += delta;
                minus[pulse] -= delta;
                let difference = (experiment.forward_loss(cx, &plus).unwrap()
                    - experiment.forward_loss(cx, &minus).unwrap())
                    / (2.0 * delta);
                assert!(
                    (evaluation.gradient[pulse] - difference).abs()
                        < 2e-6 * (1.0 + difference.abs()),
                    "radiative point={point:?}, pulse={pulse}, adjoint={}, forward FD={difference}",
                    evaluation.gradient[pulse]
                );
            }
        }

        let initial = experiment.forward_loss(cx, &START).unwrap();
        let (state, report) = fit(cx, &experiment).unwrap();
        assert_eq!(report.stop, SqpStop::Converged, "{report:?}");
        assert!(state.sample().f < initial * 1e-8);
        for (&actual, expected) in state.point().iter().zip(TRUTH) {
            assert!(
                (actual - expected).abs() < 1e-3,
                "radiative recovery {actual} W/m3, synthetic target {expected} W/m3"
            );
        }
        let final_state = experiment.evaluate(cx, state.point()).unwrap();
        assert!(final_state.max_step_energy_residual_j < 1e-9);
        assert!(final_state.final_liquid_fraction.iter().any(|f| *f > 0.99));
        assert!(final_state.final_liquid_fraction.iter().any(|f| *f < 0.01));
        assert!(
            final_state
                .final_liquid_fraction
                .iter()
                .any(|f| *f > 0.01 && *f < 0.99)
        );
        assert!((final_state.loss_k2 - state.sample().f).abs() < 1e-11);
        println!(
            "radiative heater recovery {:?}, initial loss={initial:.9e}, final loss={:.9e}, largest thermal field change={largest_change:.9e} K",
            state.point(),
            state.sample().f
        );
    });
}
