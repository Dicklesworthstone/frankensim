//! Actual caloric temperatures retain the complete loaded native rolling cycle.

use super::{DT, binding, controls as native_controls, events, fixture, thermal};
use crate::{angular_temperature as temperature, close};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::dissipation::{
        JkrRollingDissipationHeat, temperature::JkrRollingTemperatureHeatError as Error,
    },
};
use fs_exec::CancelGate;
use fs_thermal::enthalpy::{CellId, DepositionError, StepConfig};
use fs_tribo::{HeatPartition, WorkLedger};

fn ledgers(heat: &JkrRollingDissipationHeat<'_>) -> [WorkLedger; 3] {
    [
        *heat.normal_release_ledger(),
        *heat.rolling_release_ledger(),
        *heat.rolling_loss_ledger(),
    ]
}

#[test]
fn actual_cells_refresh_loaded_rolling_without_reheating_prior_losses() {
    let material = temperature::material_fixture();
    let normal = temperature::normal(&material);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, false);
    let caloric = temperature::curves(&normal);
    let mut thermal = thermal(&caloric, 1.0);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, state.source_state().clone(), &caloric, &thermal),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let gate = CancelGate::new();
    // All three ledgers are nonzero before temperature continuation.
    for delta in [0.0, normal.minimum_attached_indentation_m(), 0.0] {
        heat.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            heat.stamp(),
            delta,
            DT,
            native_controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    }
    let accounts = ledgers(&heat);
    assert!(
        accounts
            .iter()
            .all(|ledger| ledger.dissipated_work_j() > 0.0)
    );
    let remembered = heat.state().attached().unwrap().1;
    let prior_loss = heat.state().accumulated_rolling_loss_j();
    let stamp = heat.stamp();
    let bodies = f.bodies;
    temperature::warm_cells(&mut thermal);
    let cells = thermal.cells().clone();
    let time = thermal.time_s();
    let report = heat
        .apply_native_temperature_continuation(
            &mut state,
            &thermal,
            stamp.clone(),
            material.cards(),
            315.0,
            temperature::controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    close(report.temperatures.surface_a_k, 310.0);
    close(report.temperatures.surface_b_k, 305.0);
    assert_eq!(report.temperatures.interface_k, 315.0);
    assert_eq!(
        report.temperatures.surface_a_k,
        thermal.temperature_k(CellId(0)).unwrap()
    );
    assert_eq!(
        report.temperatures.surface_b_k,
        thermal.temperature_k(CellId(1)).unwrap()
    );
    assert_eq!(
        report.parameter_work_source,
        temperature::controls().parameter_work_source
    );
    assert_ne!(report.temperature.normal_parameter_work_j, 0.0);
    assert!(report.temperature.rolling_parameter_work_j < 0.0);
    assert_ne!(
        report.temperature.original_model_identity,
        report.temperature.refreshed_model_identity
    );
    let (contact, memory) = heat.state().attached().unwrap();
    close(contact.law().parameters().work_of_adhesion_j_m2, 0.0275);
    close(
        contact.law().parameters().reduced_modulus_pa,
        ((1.0 - 0.25_f64.powi(2)) * (1.0 / 9e6 + 1.0 / 19e6)).recip(),
    );
    assert_eq!(memory.elastic_shift_m(), remembered.elastic_shift_m());
    assert_eq!(
        memory.accumulated_angle_rad(),
        remembered.accumulated_angle_rad()
    );
    assert_eq!(heat.state().accumulated_rolling_loss_j(), prior_loss);
    assert_eq!(ledgers(&heat), accounts);
    assert_eq!(heat.state(), state.source_state());
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), time);
    assert_eq!(f.bodies, bodies);
    assert!(matches!(
        heat.propose_native_temperature_continuation(
            &state,
            &thermal,
            stamp,
            material.cards(),
            315.0,
            temperature::controls(),
            &gate,
            || Ok(()),
        ),
        Err(Error::HistoryMismatch)
    ));
    // The updated original source remains usable by the native heat owner,
    // with the same held geometry and the existing loss partitions.
    let next = heat
        .propose_native(
            &f.bodies,
            &f.properties,
            &state,
            f.binding,
            f.frame,
            &thermal,
            heat.stamp(),
            0.0,
            DT,
            native_controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.history().state(), next.native_state().source_state());
    assert!(next.report().mechanics.mechanical.is_some());
    assert_eq!(
        next.report().heat.release_outside_destination.as_deref(),
        Some("resolved elastic waves")
    );
    assert_eq!(
        next.report().heat.rolling_outside_destination.as_deref(),
        Some("rolling microstructural rearrangement")
    );
}

#[test]
fn rolling_temperature_work_clock_cards_and_final_cancellation_refuse_atomically() {
    let material = temperature::material_fixture();
    let normal = temperature::normal(&material);
    let events = events(&normal);
    let (f, mut state) = fixture(&events, false);
    let caloric = temperature::curves(&normal);
    let mut thermal = thermal(&caloric, 1.0);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, state.source_state().clone(), &caloric, &thermal),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    temperature::warm_cells(&mut thermal);
    let stamp = heat.stamp();
    let original = state.clone();
    let cells = thermal.cells().clone();
    let bodies = f.bodies;
    let accounts = ledgers(&heat);
    let gate = CancelGate::new();
    for (name, bound) in [
        ("", 1e-6),
        ("driver", -1.0),
        ("driver", f64::INFINITY),
        ("driver", f64::NAN),
    ] {
        let mut invalid = temperature::controls();
        invalid.parameter_work_source = name;
        invalid.maximum_absolute_parameter_work_j = bound;
        assert!(matches!(
            heat.apply_native_temperature_continuation(
                &mut state,
                &thermal,
                stamp.clone(),
                material.cards(),
                315.0,
                invalid,
                &gate,
                || Ok(()),
            ),
            Err(Error::Invalid(_))
        ));
    }
    let (_, mut foreign) = fixture(heat.source(), true);
    let foreign_before = foreign.clone();
    assert!(matches!(
        heat.apply_native_temperature_continuation(
            &mut foreign,
            &thermal,
            stamp.clone(),
            material.cards(),
            315.0,
            temperature::controls(),
            &gate,
            || Ok(()),
        ),
        Err(Error::HistoryMismatch)
    ));
    assert_eq!(foreign, foreign_before);
    let mut limits = temperature::controls();
    limits.maximum_absolute_parameter_work_j = 0.0;
    assert!(matches!(
        heat.apply_native_temperature_continuation(
            &mut state,
            &thermal,
            stamp.clone(),
            material.cards(),
            315.0,
            limits,
            &gate,
            || Ok(()),
        ),
        Err(Error::ParameterWorkLimit { .. })
    ));
    let cards = material.cards();
    assert!(matches!(
        heat.apply_native_temperature_continuation(
            &mut state,
            &thermal,
            stamp.clone(),
            [cards[1], cards[0]],
            315.0,
            temperature::controls(),
            &gate,
            || Ok(()),
        ),
        Err(Error::Temperature(_))
    ));
    let later = thermal
        .advance_to(
            0.01,
            StepConfig {
                max_dt_s: 0.01,
                energy_tolerance_j: 1e-20,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    assert!(matches!(
        heat.apply_native_temperature_continuation(
            &mut state,
            &later,
            stamp.clone(),
            material.cards(),
            315.0,
            temperature::controls(),
            &gate,
            || Ok(()),
        ),
        Err(Error::HistoryMismatch)
    ));
    let mut checks = 0;
    assert!(matches!(
        heat.apply_native_temperature_continuation(
            &mut state,
            &thermal,
            stamp.clone(),
            material.cards(),
            315.0,
            temperature::controls(),
            &gate,
            || {
                checks += 1;
                if checks == 3 {
                    gate.request();
                }
                Ok(())
            },
        ),
        Err(Error::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(checks, 3); // Final publication checkpoint, after complete staging.
    assert_eq!(state, original);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), accounts);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(f.bodies, bodies);
}
