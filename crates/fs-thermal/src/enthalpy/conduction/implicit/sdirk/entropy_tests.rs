//! Actual native stages must not destroy the entropy of an isolated body pair.
use super::*;

const CAPACITIES: [f64; 2] = [15.228_830_82, 0.681_110_36];
const INITIAL_K: [f64; 2] = [230.484_864, 2339.693_649];
const CONDUCTANCE_W_K: f64 = 2.074_837_42;
const INTERVAL_S: f64 = 6.019_507_901;

#[test]
fn passive_sdirk2_retries_positive_temperature_entropy_loss() {
    let law = Law {
        coordinate_bounds: [0.1, 4000.0],
        ..Law::constant(1.0, CONDUCTANCE_W_K)
    };
    let (thermal, observer) = network(&law, &INITIAL_K, &CAPACITIES, &[(0, 1, 1.0)], None, vec![]);
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    // The independent exact SDIRK2 stage solution stays inside the original
    // positive caloric domain and conserves energy, yet loses entropy. Domain,
    // Newton, and joule gates alone therefore cannot reject this physical error.
    let mean = (CAPACITIES[0] * INITIAL_K[0] + CAPACITIES[1] * INITIAL_K[1])
        / CAPACITIES.iter().sum::<f64>();
    let decay_dt = INTERVAL_S * CONDUCTANCE_W_K * (1.0 / CAPACITIES[0] + 1.0 / CAPACITIES[1]);
    let raw = INITIAL_K.map(|t| mean + (t - mean) * amplification(decay_dt));
    assert!(raw.iter().all(|&t| t > law.coordinate_bounds[0]));
    let raw_entropy =
        CAPACITIES[0] * (raw[0] / INITIAL_K[0]).ln() + CAPACITIES[1] * (raw[1] / INITIAL_K[1]).ln();
    assert!(raw_entropy < -0.17);

    let (next, report) = solver
        .advance_sdirk2_to(
            &thermal,
            INTERVAL_S,
            controls(INTERVAL_S),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert!(
        report.rejected_intervals > 0,
        "accepted entropy change: {:?}",
        report.entropy
    );
    assert!(report.energy.substeps > 1);
    assert!(report.entropy.thermal_storage_change_j_k >= -1e-10);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(next.time_s(), INTERVAL_S);
}

#[test]
fn sparse_passive_sdirk2_entropy_admission_ignores_caloric_reference_offsets() {
    let mut endpoints = Vec::new();
    for datum in [0.0, 1e12] {
        let law = Law {
            datum,
            coordinate_bounds: [0.1, 4000.0],
            ..Law::constant(1.0, CONDUCTANCE_W_K)
        };
        let (thermal, observer) =
            network(&law, &INITIAL_K, &CAPACITIES, &[(0, 1, 1.0)], None, vec![]);
        let solver = ImplicitConductionGraph::new(
            &thermal,
            &observer,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        // These joule tolerances resolve the actual stored coordinate even
        // after its 10^12 J/kg reference shift. The entropy allowance must not
        // grow with that arbitrary datum and hide the 0.1757 J/K destruction.
        let (next, report) = solver
            .advance_sdirk2_to(
                &thermal,
                INTERVAL_S,
                ImplicitConductionControls {
                    cell_tolerance_j: 0.003,
                    total_tolerance_j: 0.01,
                    ..controls(INTERVAL_S)
                },
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        assert!(report.rejected_intervals > 0);
        assert!(report.entropy.thermal_storage_change_j_k > 3.5);
        assert!(report.energy.residual_j.abs() <= 0.01);
        endpoints.push([
            next.temperature_k(id(0)).unwrap(),
            next.temperature_k(id(1)).unwrap(),
        ]);
    }
    for i in 0..2 {
        near(endpoints[0][i], endpoints[1][i], 0.01);
    }
}

#[test]
fn passive_sdirk2_checks_cell_plus_actual_bath_entropy() {
    let law = Law {
        coordinate_bounds: [0.1, 4000.0],
        ..Law::constant(1.0, 1.0)
    };
    let dt = 6.3;
    // Both native bath channels have the same independently soluble 1 W/K
    // response. Ordinary cooling loses CELL entropy; the actual bath must be
    // included before testing the second law.
    for reservoir in [false, true] {
        let (mut thermal, observer) = network(
            &law,
            &[1800.0],
            &[1.0],
            &[],
            reservoir.then_some(300.0),
            vec![],
        );
        if !reservoir {
            thermal.cells.get_mut(&id(0)).unwrap().boundary.htc_w_m2_k = 0.25;
        }
        let solver = chain(&thermal, &observer, &law);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let raw = 300.0 + 1500.0 * amplification(dt);
        assert!(raw > law.coordinate_bounds[0]);
        assert!((raw / 1800.0).ln() - (raw - 1800.0) / 300.0 < -0.93);

        let (next, report) = solver
            .advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        assert!(report.rejected_intervals > 0);
        assert!(report.energy.substeps > 1);
        assert!(report.entropy.thermal_storage_change_j_k < 0.0);
        assert!(report.entropy.reservoir_entropy_input_j_k < 0.0);
        assert!(
            report.entropy.thermal_storage_change_j_k - report.entropy.reservoir_entropy_input_j_k
                >= -1e-10
        );
        near(
            report.energy.external_heat_j,
            next.temperature_k(id(0)).unwrap() - 1800.0,
            1e-8,
        );
        near(
            report.entropy.reservoir_entropy_input_j_k,
            report.energy.external_heat_j / 300.0,
            1e-12,
        );
    }
}

#[test]
fn passive_sdirk2_entropy_retry_exhaustion_publishes_no_state() {
    let law = Law {
        coordinate_bounds: [0.1, 4000.0],
        ..Law::constant(1.0, 1.0)
    };
    let (thermal, observer) = network(&law, &[1800.0], &[1.0], &[], Some(300.0), vec![]);
    let solver = chain(&thermal, &observer, &law);
    let before = thermal.cells().clone();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let result = solver.advance_sdirk2_to(
        &thermal,
        6.3,
        ImplicitConductionControls {
            maximum_halvings: 0,
            ..controls(6.3)
        },
        &mut workspace,
        || Ok(()),
    );
    assert!(
        matches!(result, Err(ImplicitConductionError::EntropyDecrease {
        entropy_change_j_k, allowance_j_k,
    }) if entropy_change_j_k < -0.93 && allowance_j_k < 1e-8)
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);

    let result = solver.advance_sdirk2_adaptive_to(
        &thermal,
        6.3,
        AdaptiveImplicitConductionControls {
            step: ImplicitConductionControls {
                maximum_halvings: 0,
                ..controls(6.3)
            },
            initial_dt_s: 6.3,
            enthalpy_absolute_tolerance_j_kg: 10_000.0,
            temperature_absolute_tolerance_k: 10_000.0,
            relative_tolerance: 0.0,
            ..Default::default()
        },
        &mut workspace,
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::EntropyDecrease { .. },
        ))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    // Scratch may have changed; the same original state remains continuable.
    let (next, _) = solver
        .advance_sdirk2_to(&thermal, 0.25, controls(0.25), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(next.time_s(), 0.25);
}

#[test]
fn adaptive_sdirk2_retries_entropy_loss_with_loose_accuracy_tolerances() {
    let law = Law {
        coordinate_bounds: [0.1, 4000.0],
        ..Law::constant(1.0, 1.0)
    };
    let (thermal, observer) = network(&law, &[1800.0], &[1.0], &[], Some(300.0), vec![]);
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = solver
        .advance_sdirk2_adaptive_to(
            &thermal,
            6.3,
            AdaptiveImplicitConductionControls {
                step: controls(6.3),
                initial_dt_s: 6.3,
                enthalpy_absolute_tolerance_j_kg: 10_000.0,
                temperature_absolute_tolerance_k: 10_000.0,
                relative_tolerance: 0.0,
                ..Default::default()
            },
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert!(report.solve.rejected_intervals > 0);
    assert_eq!(report.accuracy_rejections, 0);
    assert!(
        report.solve.entropy.thermal_storage_change_j_k
            - report.solve.entropy.reservoir_entropy_input_j_k
            >= -1e-10
    );
    near(
        report.solve.energy.external_heat_j,
        next.temperature_k(id(0)).unwrap() - 1800.0,
        1e-8,
    );
    assert_eq!(next.time_s(), 6.3);
}

#[test]
fn source_and_classified_stream_cooling_keep_their_original_entropy_scope() {
    let law = Law::constant(1.0, 0.0);
    let (mut thermal, observer) = network(&law, &[600.0], &[1.0], &[], None, vec![]);
    thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = -100.0;
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (_, report) = solver
        .advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    assert!(report.entropy.thermal_storage_change_j_k < 0.0);
    assert!(report.entropy.has_unclassified_sources);
    assert_eq!(report.rejected_intervals, 0);

    thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = 0.0;
    let solver = chain(&thermal, &observer, &law)
        .with_implicit_advection(&thermal, &throughflow(1, 1.0, 300.0), || Ok(()))
        .unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(()))
        .unwrap();
    let (_, report) = solver
        .advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    assert!(report.entropy.thermal_storage_change_j_k < 0.0);
    assert!(!report.entropy.has_unclassified_sources);
    assert!(report.entropy.advection.is_some());
    assert_eq!(report.rejected_intervals, 0);
}
