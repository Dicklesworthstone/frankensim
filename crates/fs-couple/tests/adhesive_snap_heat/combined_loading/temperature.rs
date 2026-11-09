//! Prescribed thermal work retains the actual combined release/receiver history.
#[path = "../../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod original;

use super::{Segment, binding, controls, curves, driver, mixed_source, thermal};
use fs_contact::normal_patch::adhesive::{
    SourceJkrSphereContact,
    mixed_mode::{SourceJkrMixedMode, combined::loading::JkrCombinedLoadingState},
    temperature::{AdhesiveTemperatures, SourceJkrThermalRefresh},
};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError, combined_loading::JkrCombinedHeatError,
};
use fs_exec::CancelGate;
use fs_thermal::enthalpy::{DepositionError, StepConfig};

fn source(data: &original::Fixture) -> SourceJkrMixedMode {
    let template = mixed_source();
    let normal = template.normal_source();
    SourceJkrMixedMode::new(
        SourceJkrSphereContact::new(
            data.normal.clone(),
            normal.effective_radius_m(),
            normal.interaction_range_m(),
            normal.envelope(),
        )
        .unwrap(),
        template.envelope(),
    )
    .unwrap()
}

fn refresh(
    source: &SourceJkrMixedMode,
    state: JkrCombinedLoadingState,
    data: &original::Fixture,
    temperatures: [f64; 3],
) -> (SourceJkrThermalRefresh, SourceJkrMixedMode) {
    let normal = source.normal_source();
    let indentation = state.coordinates()[0];
    let witness = if state.contact().is_some() {
        normal.evaluate_attached(indentation).unwrap().state()
    } else {
        normal.detached_state(indentation).unwrap()
    };
    let update = normal
        .refresh_at_held_indentation(
            witness,
            data.cards(),
            AdhesiveTemperatures {
                surface_a_k: temperatures[0],
                surface_b_k: temperatures[1],
                interface_k: temperatures[2],
            },
        )
        .unwrap();
    let successor = SourceJkrMixedMode::new(update.source().clone(), source.envelope()).unwrap();
    (update, successor)
}

#[test]
fn loaded_temperature_work_preserves_receivers_and_resumes_without_depositing_heat() {
    let data = original::Fixture::with_options(10.0, 0.05, true, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.2, &gate).unwrap();
    let caloric = curves(source.normal_source(), 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let heat = binding(loading, initial, &caloric, &thermal, &gate);
    let old_stamp = heat.stamp();
    // The original heat owner can transport independently while the prescribed
    // mechanical coordinates remain held. The successor pins the new clock.
    thermal = thermal
        .advance_to(
            0.1,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-20,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let cells = thermal.cells().clone();
    let (refresh, successor) = refresh(&source, initial, &data, [302.0, 301.0, 302.0]);
    let continuation = loading
        .continue_at_temperature(initial, &refresh, &successor, &gate)
        .unwrap();
    assert_ne!(continuation.report().parameter_work_j, 0.0);
    assert_ne!(
        continuation.report().previous_limit,
        continuation.report().refreshed_limit
    );
    let mut count = 0;
    assert!(matches!(
        heat.propose_temperature_continuation(&thermal, old_stamp, &continuation, || {
            count += 1;
            if count == 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(JkrCombinedHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(heat.stamp(), old_stamp);
    assert_eq!(heat.state(), initial);
    assert_eq!(thermal.cells(), &cells);
    let mut continued = heat
        .propose_temperature_continuation(&thermal, old_stamp, &continuation, || Ok(()))
        .unwrap();
    assert_eq!(continued.state(), continuation.state());
    assert_eq!(continued.ledger(), heat.ledger());
    assert_eq!(continued.state().steps(), initial.steps());
    assert_eq!(continued.state().coordinates(), initial.coordinates());
    assert_eq!(thermal.cells(), &cells);
    assert!(matches!(
        continued.propose(
            &thermal,
            old_stamp,
            &[Segment::LoadFactor(0.3)],
            controls(),
            &gate,
            || Ok(())
        ),
        Err(JkrCombinedHeatError::Heat(
            JkrSnapHeatError::HistoryMismatch
        ))
    ));
    let resumed = continued
        .apply(
            &mut thermal,
            continued.stamp(),
            &[Segment::LoadFactor(0.3)],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert!(resumed.prescribed_work_j > 0.0);
    assert_eq!(resumed.released_energy_j, 0.0);
    assert_eq!(resumed.thermal.stored_energy_change_j, 0.0);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.1);
}

#[test]
fn previously_heated_detachment_survives_loss_of_every_attached_equilibrium() {
    let data = original::Fixture::with_options(10.0, 0.05, true, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -2.9, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let caloric = curves(source.normal_source(), 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
    heat.apply(
        &mut thermal,
        heat.stamp(),
        &[
            Segment::LoadFactor(2.0 * loading.limit().unwrap().load_factor),
            Segment::DetachAtFold,
        ],
        controls(),
        &gate,
        || Ok(()),
    )
    .unwrap();
    let old = heat.state();
    let ledger = *heat.ledger();
    let cells = thermal.cells().clone();
    assert!(ledger.surface_a_heat_j() > 0.0 && ledger.surface_b_heat_j() > 0.0);
    let (refresh, successor) = refresh(&source, old, &data, [320.0; 3]);
    assert!(
        successor
            .normal_source()
            .evaluate_attached(old.coordinates()[0])
            .is_err()
    );
    let continuation = loading
        .continue_at_temperature(old, &refresh, &successor, &gate)
        .unwrap();
    assert_eq!(continuation.report().parameter_work_j, 0.0);
    let mut continued = heat
        .propose_temperature_continuation(&thermal, heat.stamp(), &continuation, || Ok(()))
        .unwrap();
    assert_eq!(*continued.ledger(), ledger);
    assert_eq!(
        continued.state().released_energy_j(),
        old.released_energy_j()
    );
    assert_eq!(continued.state().contact(), None);
    assert_eq!(thermal.cells(), &cells);
    continued
        .apply(
            &mut thermal,
            continued.stamp(),
            &[Segment::LoadFactor(0.0), Segment::LoadFactor(-1.0)],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(continued.state().contact(), None);
    assert_eq!(
        continued.state().released_energy_j(),
        old.released_energy_j()
    );
    assert_eq!(*continued.ledger(), ledger);
    assert_eq!(thermal.cells(), &cells);
}
