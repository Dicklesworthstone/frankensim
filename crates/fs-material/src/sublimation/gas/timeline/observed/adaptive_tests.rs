//! G0/G1/G3/G4: adaptive deposition with original physical heat observers.
use super::*;

fn adaptive_controls(cp: f64, tolerance: f64) -> GasTimelineAdaptiveSolve {
    GasTimelineAdaptiveSolve {
        absolute_mass_tolerance_kg: 1e-11,
        absolute_temperature_tolerance_k: tolerance,
        absolute_specific_enthalpy_tolerance_j_kg: cp * tolerance,
        relative_tolerance: 0.0,
        minimum_dt_s: 1e-9,
        maximum_trials: 10000,
    }
}

#[test]
fn empty_observer_preserves_adaptive_physics_and_work_bitwise() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let plain = crate::sublimation::gas::timeline::adaptive::advance_deposition_adaptive(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        adaptive_controls(cp, 1e-5),
        || Ok(()),
    )
    .unwrap();
    let (observed, ()) = advance_deposition_adaptive_with_gas(
        &thermal,
        &surfaces,
        None,
        &(),
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        adaptive_controls(cp, 1e-5),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(observed.endpoint().gas(), plain.endpoint().gas());
    assert_eq!(
        observed.endpoint().thermal().cells(),
        plain.endpoint().thermal().cells()
    );
    assert_eq!(observed.endpoint().audit(), plain.endpoint().audit());
    assert_eq!(observed.audit(), plain.audit());
}

#[test]
fn accepted_fine_history_counts_source_heat_once_across_rejections_and_transitions() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let original = thermal.cells().clone();
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let pulse = Pulse::default();
    let (candidate, accepted) = advance_deposition_adaptive_with_gas(
        &thermal,
        &surfaces,
        None,
        &pulse,
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        adaptive_controls(cp, 1e-5),
        || Ok(()),
    )
    .unwrap();
    let physical = candidate.endpoint().audit();
    assert!(candidate.audit().rejected_trials > 0);
    assert!(candidate.audit().maximum_accepted_error_ratio <= 1.0);
    assert!(candidate.audit().thermal_substeps_evaluated > physical.thermal.substeps);
    assert_eq!(accepted.begin_calls, 1);
    assert_eq!(accepted.accepted_steps, physical.thermal.substeps);
    fixture::near(accepted.accepted_duration_s, 0.004, 1e-15);
    // Exact integral across the declared 2 ms transition: 10 W then -5 W.
    fixture::near(accepted.interval_heat_j, 0.01, 1e-14);
    fixture::near(physical.thermal.external_heat_j, 0.01, 1e-14);
    assert!(physical.transferred_mass_kg < 0.0);
    fixture::near(physical.mass_residual_kg, 0.0, 1e-14);
    fixture::near(physical.energy_residual_j, 0.0, 1e-7);
    assert_eq!(pulse, Pulse::default());
    assert_eq!(thermal.cells(), &original);
}

#[test]
fn adaptive_growing_layer_conduction_converges_to_resolved_original_transport() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let distance = MASS / (2.0 * DENSITY * AREA);
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &Conductivity)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(1),
            half_distance_m: distance,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_planar_phase_layers(
        &thermal,
        &[PlanarPhaseLayer {
            cell_id: CellId(1),
            footprint_area_m2: AREA,
            density_kg_m3: DENSITY,
        }],
        || Ok(()),
    )
    .unwrap();
    let (reference, _) = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        None,
        &conduction,
        &gas,
        &catalog,
        0.02,
        heat(0.02 / 4096.0),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let run = |tolerance| {
        advance_deposition_adaptive_with_gas(
            &thermal,
            &surfaces,
            None,
            &conduction,
            &gas,
            &catalog,
            0.02,
            heat(0.01),
            controls(),
            adaptive_controls(cp, tolerance),
            || Ok(()),
        )
        .unwrap()
    };
    let (loose, _) = run(0.01);
    let (tight, accepted_conduction) = run(1e-6);
    let expected_t = reference.thermal().temperature_k(CellId(1)).unwrap();
    let errors = [&loose, &tight].map(|candidate| {
        (candidate
            .endpoint()
            .thermal()
            .temperature_k(CellId(1))
            .unwrap()
            - expected_t)
            .abs()
    });
    assert!(
        errors[1] < 0.2 * errors[0],
        "temperature errors: {errors:?}"
    );
    assert!(errors[1] < 0.002, "temperature errors: {errors:?}");
    assert!(tight.audit().accepted_intervals > loose.audit().accepted_intervals);
    let endpoint = tight.endpoint();
    let mass = endpoint.thermal().cells()[&CellId(1)].mass_kg;
    assert!(mass > MASS);
    let last_distance = accepted_conduction
        .last_planar_half_distance_m(CellId(1))
        .unwrap();
    assert!(last_distance > distance && last_distance <= mass / (2.0 * DENSITY * AREA));
    fixture::near(endpoint.audit().mass_residual_kg, 0.0, 1e-14);
    fixture::near(endpoint.audit().energy_residual_j, 0.0, 1e-7);
    fixture::near(
        conduction.last_planar_half_distance_m(CellId(1)).unwrap(),
        distance,
        1e-16,
    );
}

#[test]
fn late_observer_refusal_and_cancellation_publish_no_adaptive_history() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let original = thermal.cells().clone();
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let original_gas = gas.clone();
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let refusing = Pulse {
        fail_after: Some(5),
        ..Pulse::default()
    };
    let result = advance_deposition_adaptive_with_gas(
        &thermal,
        &surfaces,
        None,
        &refusing,
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        adaptive_controls(cp, 1e-5),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(GasTimelineAdaptiveObservedError::Observer(
            SourceError::Rejected
        ))
    ));
    assert_eq!(refusing.accepted_steps, 0);
    assert_eq!(refusing.begin_calls, 0);
    let pulse = Pulse::default();
    let mut calls = 0;
    let accepted = advance_deposition_adaptive_with_gas(
        &thermal,
        &surfaces,
        None,
        &pulse,
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        adaptive_controls(cp, 1e-5),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut replay_calls = 0;
    let cancelled = advance_deposition_adaptive_with_gas(
        &thermal,
        &surfaces,
        None,
        &pulse,
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        adaptive_controls(cp, 1e-5),
        || {
            replay_calls += 1;
            if replay_calls == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(GasTimelineAdaptiveObservedError::Adaptive(
            GasTimelineAdaptiveError::Timeline(GasTimelineError::Thermal(
                DepositionError::Cancelled
            ))
        ))
    ));
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(gas, original_gas);
    assert_eq!(pulse, Pulse::default());
    assert_eq!(accepted.1.begin_calls, 1);
}
