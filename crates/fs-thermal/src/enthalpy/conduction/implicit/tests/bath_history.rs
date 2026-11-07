//! Physical-time bath forcing through the original heat and entropy owners.
use super::super::boundary::{
    BathTemperatureKnot, BathTemperatureSchedule, GrayAmbientBoundary, GrayAmbientControls,
    GrayTwoBathBoundary,
};
use super::*;

fn history(samples: &[(f64, f64)]) -> BathTemperatureSchedule {
    BathTemperatureSchedule::new(
        samples
            .iter()
            .map(|&(time_s, temperature_k)| BathTemperatureKnot {
                time_s,
                temperature_k,
            })
            .collect(),
        samples.len(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn ramped_reservoir_uses_actual_endpoint_and_both_stage_entropy_temperatures() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[400.0], Some(300.0), vec![]);
    let bath = history(&[(0.0, 300.0), (0.5, 340.0), (1.0, 320.0), (1.5, 320.0)]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_bath_temperature_histories(&thermal, &[(CellId(0), &bath)], &[], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();

    let (next, report) = chain
        .advance_to(&thermal, 0.5, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), 380.0, 1e-10);
    near(report.energy.external_heat_j, -20.0, 1e-10);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        -20.0 / 340.0,
        1e-12,
    );
    assert!(!report.entropy.has_unclassified_sources);
    assert_eq!(
        next.cells()[&CellId(0)].boundary,
        thermal.cells()[&CellId(0)].boundary
    );

    // Independent scalar SDIRK solution for C=G=1 and Tb(t)=300+80t.
    let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
    let dt = 0.5;
    let bath1 = 300.0 + 80.0 * gamma * dt;
    let t1 = (400.0 + gamma * dt * bath1) / (1.0 + gamma * dt);
    let q1 = bath1 - t1;
    let t2 = (400.0 + (1.0 - gamma) * dt * q1 + gamma * dt * 340.0) / (1.0 + gamma * dt);
    let q2 = 340.0 - t2;
    let (next, report) = chain
        .advance_sdirk2_to(&thermal, dt, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), t2, 1e-10);
    near(report.energy.external_heat_j, t2 - 400.0, 1e-10);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        dt * ((1.0 - gamma) * q1 / bath1 + gamma * q2 / 340.0),
        1e-12,
    );
    assert!(!report.entropy.has_unclassified_sources);
    assert!(report.entropy.conduction_generation_j_k > 0.0);

    // Exact affine-forced lumped solution, independently integrated over each
    // bath segment; the native face has unit half-cell Fourier conductance.
    let mut exact = 400.0_f64;
    for (duration, start, slope) in [(0.5_f64, 300.0, 80.0), (0.5, 340.0, -40.0), (0.5, 320.0, 0.0)] {
        exact = start + slope * duration - slope + (exact - start + slope) * (-duration).exp();
    }
    let mut errors = Vec::new();
    for dt in [0.125, 0.0625, 0.03125] {
        let (next, report) = chain
            .advance_sdirk2_to(&thermal, 1.5, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        errors.push((next.temperature_k(CellId(0)).unwrap() - exact).abs());
        near(
            report.energy.external_heat_j,
            next.temperature_k(CellId(0)).unwrap() - 400.0,
            1e-8,
        );
        near(report.energy.residual_j, 0.0, 1e-8);
    }
    assert!(errors[0] / errors[1] > 3.7, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.7, "{errors:?}");

    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_bath_temperature_histories(&thermal, &[(CellId(0), &bath)], &[], || Ok(()))
    .unwrap();
    for (next, report) in [
        chain
            .advance_sdirk2_to(&thermal, 1.5, controls(2.0), &mut workspace, || Ok(()))
            .unwrap(),
        graph
            .advance_sdirk2_to(&thermal, 1.5, controls(2.0), &mut workspace, || Ok(()))
            .unwrap(),
    ] {
        assert_eq!(report.energy.substeps, 3);
        assert!(next.temperature_k(CellId(0)).unwrap() < 400.0);
        near(report.energy.residual_j, 0.0, 1e-8);
    }
    let adaptive = adaptive::AdaptiveImplicitConductionControls {
        step: controls(2.0),
        initial_dt_s: 2.0,
        enthalpy_absolute_tolerance_j_kg: 1e-5,
        temperature_absolute_tolerance_k: 1e-5,
        relative_tolerance: 0.0,
        ..Default::default()
    };
    let (next, report) = chain
        .advance_sdirk2_adaptive_to(&thermal, 1.5, adaptive, &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), exact, 1e-3);
    assert!(report.maximum_accepted_normalized_error <= 1.0);
    assert!(report.accuracy_rejections > 0);
    near(report.solve.energy.residual_j, 0.0, 1e-8);
}

#[test]
fn scheduled_baths_keep_native_nonlinear_and_resolved_radiative_tangents() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let (mut thermal, conduction) = setup(&law, &[400.0], Some(300.0), vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.5;
    let bath = history(&[(0.0, 300.0), (1.0, 400.0)]);
    let gray = GrayAmbientBoundary::new(&law, 0.8, 9.0, GrayAmbientControls::default()).unwrap();
    let paired =
        GrayTwoBathBoundary::new(&law, 0.8, 9.0, 450.0, GrayAmbientControls::default()).unwrap();
    for variant in 0..3 {
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap();
        let chain = match variant {
            1 => chain
                .with_boundary_laws(&thermal, &[(CellId(0), &gray)], &[], || Ok(()))
                .unwrap(),
            2 => chain
                .with_boundary_exchanges(&thermal, &[(CellId(0), &paired)], &[], || Ok(()))
                .unwrap(),
            _ => chain,
        }
        .with_bath_temperature_histories(&thermal, &[(CellId(0), &bath)], &[], || Ok(()))
        .unwrap();
        let dt = 0.25;
        let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
        let energy = [3.0];
        let mut report = empty_report();
        let trial = chain
            .evaluate(
                &step,
                &energy,
                &BTreeMap::new(),
                &[],
                true,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        let exchange = trial.payload().reservoir_exchange_states[&CellId(0)];
        near(exchange.baths[0].temperature_k, 325.0, 0.0);
        if variant == 2 {
            near(exchange.baths[1].temperature_k, 450.0, 0.0);
            assert!(exchange.baths[0].outward_power_w > 0.0);
            assert!(exchange.baths[1].outward_power_w < 0.0);
        }
        near(
            trial.entropy_audit().unwrap().reservoir_entropy_input_j_k,
            -dt * exchange
                .baths
                .iter()
                .map(|bath| bath.outward_power_w / bath.temperature_k)
                .sum::<f64>(),
            1e-12,
        );
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        chain
            .jacobian(
                &thermal,
                trial.payload(),
                dt,
                &mut workspace,
                &mut || Ok(()),
            )
            .unwrap();
        let mut residuals = Vec::new();
        for perturbation in [-1e-3, 1e-3] {
            let trial = chain
                .evaluate(
                    &step,
                    &[energy[0] + perturbation],
                    &BTreeMap::new(),
                    &[],
                    false,
                    controls(1.0),
                    &mut report,
                    &mut || Ok(()),
                )
                .unwrap();
            residuals.push(trial.residual_j()[0]);
        }
        near(
            workspace.diagonal[0],
            (residuals[1] - residuals[0]) / 0.002,
            2e-8,
        );
    }
}

#[test]
fn native_robin_history_removes_latent_heat_without_a_temperature_surrogate() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    let (mut thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    thermal
        .cells
        .get_mut(&CellId(0))
        .unwrap()
        .boundary
        .htc_w_m2_k = 0.5;
    let bath = history(&[(0.0, 200.0), (1.0, 250.0)]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_bath_temperature_histories(&thermal, &[], &[(CellId(0), &bath)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    // A_exposed=2, htc=.5, T_cell=300 throughout the latent plateau:
    // integral Q_in dt = integral(200+50t-300)dt = -75 J.
    near(next.cells()[&CellId(0)].specific_h_j_kg, 325.0, 1e-10);
    near(next.temperature_k(CellId(0)).unwrap(), 300.0, 0.0);
    near(report.energy.external_heat_j, -75.0, 1e-10);
    near(
        report.entropy.thermal_storage_change_j_k,
        -75.0 / 300.0,
        1e-12,
    );
    assert!(!report.entropy.has_unclassified_sources);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(
        next.cells()[&CellId(0)].boundary,
        thermal.cells()[&CellId(0)].boundary
    );
}

#[test]
fn bath_history_refusals_preserve_original_state_and_event_ownership() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let mut changed = boundary(Some(300.0));
    changed.reservoir.as_mut().unwrap().temperature_k = 350.0;
    let (thermal, conduction) = setup(
        &law,
        &[400.0],
        Some(300.0),
        vec![TimedEvent {
            time_s: 0.5,
            event: Event::Boundary(CellId(0), changed),
        }],
    );
    let bath = history(&[(0.0, 300.0), (1.0, 400.0)]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_bath_temperature_histories(&thermal, &[(CellId(0), &bath)], &[], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    assert!(matches!(
        chain.advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Invalid(
            "native boundary event conflicts with prescribed bath-temperature history",
        )),
    ));
    assert!(matches!(
        chain.advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Invalid(
            "native boundary event conflicts with prescribed bath-temperature history",
        )),
    ));
    near(thermal.time_s(), 0.0, 0.0);
    near(thermal.temperature_k(CellId(0)).unwrap(), 400.0, 0.0);

    let (thermal, conduction) = setup(&law, &[400.0], Some(300.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_bath_temperature_histories(&thermal, &[(CellId(0), &bath)], &[], || Ok(()))
            .unwrap();
    assert!(
        chain
            .advance_sdirk2_to(&thermal, 1.5, controls(1.0), &mut workspace, || Ok(()),)
            .is_err()
    );
    assert!(
        chain
            .advance_to(&thermal, 0.25, controls(1.0), &mut workspace, || Err(
                DepositionError::Cancelled
            ),)
            .is_err()
    );
    near(thermal.time_s(), 0.0, 0.0);
    near(thermal.temperature_k(CellId(0)).unwrap(), 400.0, 0.0);
    assert_eq!(bath.next_transition_time_s(0.0), Some(1.0));
    assert_eq!(bath.next_transition_time_s(1.0), None);
    near(bath.temperature_at(0.25).unwrap(), 325.0, 0.0);
    for time in [f64::NAN, -0.01, 1.0_f64.next_up()] {
        assert!(bath.temperature_at(time).is_err());
    }
    assert!(BathTemperatureSchedule::new(vec![], 2, || Ok(())).is_err());
    assert!(BathTemperatureSchedule::new(bath.knots().to_vec(), 1, || Ok(())).is_err());
    assert!(
        BathTemperatureSchedule::new(bath.knots().to_vec(), 2, || Err(DepositionError::Cancelled),)
            .is_err()
    );
    for samples in [
        [(0.0, 300.0), (0.0, 400.0)],
        [(1.0, 300.0), (0.0, 400.0)],
        [(0.0, 0.0), (1.0, 300.0)],
        [(0.0, 300.0), (1.0, f64::NAN)],
        [(0.0, 300.0), (f64::INFINITY, 400.0)],
    ] {
        assert!(
            BathTemperatureSchedule::new(
                samples
                    .into_iter()
                    .map(|(time_s, temperature_k)| BathTemperatureKnot {
                        time_s,
                        temperature_k,
                    })
                    .collect(),
                2,
                || Ok(()),
            )
            .is_err()
        );
    }
}
