//! Original material temperature and heat with finite radioactive neutron stock.

use super::*;
use crate::burnup::transport::{
    BurnupTransportLimits, BurnupTransportState, NuclideOutflow, NuclideTransfer, NuclideTransport,
};
use crate::diffusion::burnup::decay::BurnupDecayNeutronEmission;

fn radioactive() -> Arc<FusionNeutronActivationState> {
    let original = activation_state(&[[100.0, 20.0, 0.0], [80.0, 0.0, 0.0]]);
    let emissions = original
        .inventories()
        .iter()
        .map(|inventory| {
            Arc::new(
                BurnupDecayNeutronEmission::new(
                    "manufactured-radioactive-product".into(),
                    Arc::clone(inventory.model()),
                    vec![vec![0.0], vec![0.125], vec![0.0]],
                )
                .unwrap(),
            )
        })
        .collect();
    let spatial = original
        .spatial()
        .clone()
        .with_decay_neutrons(emissions)
        .unwrap();
    Arc::new(FusionNeutronActivationState::new(Arc::clone(original.model()), spatial).unwrap())
}

#[test]
fn g1_internal_emission_heats_original_material_cells_without_external_neutrons() {
    let initial = material_cells();
    let state = radioactive();
    let sources = FusionNeutronActivationCellSource::bind_without_external_source(
        Arc::clone(&state),
        Time::new(0.0),
        bindings(),
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let (warm, accepted, _) = heating
        .advance_to(&initial, 0.5, controls(0.05), || Ok(()))
        .unwrap();
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    let last = accepted.last_step(CellId(20)).unwrap().activation();
    assert_eq!(last.external_injected_neutrons, 0.0);
    assert!(last.decay_neutrons > 0.0);
    for (cell_index, (id, initial_k, capacity)) in
        [(20, 300.0, 8.0), (10, 450.0, 6.0)].into_iter().enumerate()
    {
        let events = left.state().inventories()[cell_index].cumulative_reaction_counts();
        let c = events[0];
        let d = events[1];
        let heat = 5.0 * c + 2.0 * d;
        close(
            warm.temperature_k(CellId(id)).unwrap(),
            initial_k + heat / capacity,
            1e-8,
        );
        let exports = left.state().cell_released(cell_index).unwrap();
        close(exports.radiation_j, 3.5 * c + 0.5 * d, 1e-8);
        close(exports.escaped_j, 0.5 * (c + d), 1e-8);
        assert_eq!(
            warm.cells()[&CellId(id)].mass_kg,
            initial.cells()[&CellId(id)].mass_kg
        );
    }
    assert!(left.state().inventories()[1].populations()[0] < 80.0);
    close(
        warm.energy_j() - initial.energy_j(),
        accepted.cumulative_released().unwrap().deposited_j,
        1e-8,
    );
}

#[test]
fn g3_g4_decay_continues_after_external_history_and_rolls_back_a_cancelled_material_trial() {
    let initial = material_cells();
    let state = radioactive();
    let sources = FusionNeutronActivationCellSource::bind_with_history(
        Arc::clone(&state),
        Time::new(0.0),
        bindings(),
        history(&state),
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let (middle, checkpoint, _) = heating
        .advance_to(&initial, 0.5, controls(0.05), || Ok(()))
        .unwrap();
    let (warm, final_state, _) = checkpoint
        .advance_to(&middle, 0.75, controls(0.05), || Ok(()))
        .unwrap();
    let report = final_state.last_step(CellId(20)).unwrap().activation();
    assert_eq!(report.external_injected_neutrons, 0.0);
    assert!(report.decay_neutrons > 0.0 && report.released.deposited_j > 0.0);
    let accepted = Arc::clone(checkpoint.source(CellId(20)).unwrap().state());
    let mut polls = 0;
    checkpoint
        .advance_to(&middle, 0.75, controls(0.05), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut remaining = polls / 2;
    assert!(
        checkpoint
            .advance_to(&middle, 0.75, controls(0.05), || {
                if remaining == 0 {
                    Err(DepositionError::Cancelled)
                } else {
                    remaining -= 1;
                    Ok(())
                }
            })
            .is_err()
    );
    assert!(Arc::ptr_eq(
        checkpoint.source(CellId(20)).unwrap().state(),
        &accepted
    ));
    assert_eq!(middle.time_s(), 0.5);
    let (retry, _, _) = checkpoint
        .advance_to(&middle, 0.75, controls(0.05), || Ok(()))
        .unwrap();
    close(retry.energy_j(), warm.energy_j(), 1e-8);
}

#[test]
fn g1_g3_g4_transported_emitters_heat_their_actual_material_and_outflow_has_no_future_source() {
    // With no capture targets, the original radioactive product obeys
    // A'=-1.2 A, B'=0.6 A-0.6 B. Both cells decay at 0.4/s and lose
    // emitters to the exterior at 0.2/s; only A transfers stock at 0.6/s.
    // Stable daughters and already released heat remain in their event cell.
    let original = activation_state(&[[0.0, 20.0, 0.0], [0.0, 0.0, 0.0]]);
    let limits = BurnupTransportLimits::default();
    let transport = BurnupTransportState::new(
        original.inventories().to_vec(),
        original.neutrons().model().cell_volumes_m3().to_vec(),
        NuclideTransport {
            transfers: vec![NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![0.0, 0.6, 0.0],
            }],
            outflows: (0..2)
                .map(|compartment| NuclideOutflow {
                    compartment,
                    rates_per_s: vec![0.0, 0.2, 0.0],
                })
                .collect(),
        },
        limits,
        || Ok(()),
    )
    .unwrap();
    let emissions = original
        .inventories()
        .iter()
        .map(|inventory| {
            Arc::new(
                BurnupDecayNeutronEmission::new(
                    "manufactured-transported-radioactive-product".into(),
                    Arc::clone(inventory.model()),
                    vec![vec![0.0], vec![0.125], vec![0.0]],
                )
                .unwrap(),
            )
        })
        .collect();
    let spatial = DiffusionBurnupState::new_with_transport(
        original.neutrons().clone(),
        transport,
        original.spatial().bindings().to_vec(),
        limits,
    )
    .unwrap()
    .with_decay_neutrons(emissions)
    .unwrap();
    let state =
        Arc::new(FusionNeutronActivationState::new(Arc::clone(original.model()), spatial).unwrap());
    let initial = material_cells();
    let sources = FusionNeutronActivationCellSource::bind_without_external_source(
        Arc::clone(&state),
        Time::new(0.0),
        bindings(),
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let (middle, checkpoint, _) = heating
        .advance_to(&initial, 0.5, controls(0.25), || Ok(()))
        .unwrap();
    let mut polls = 0;
    let (warm, accepted, _) = checkpoint
        .advance_to(&middle, 1.5, controls(0.25), || {
            polls += 1;
            Ok(())
        })
        .unwrap();

    for (time, material, heating) in [(0.5_f64, &middle, &checkpoint), (1.5, &warm, &accepted)] {
        let shared = heating.source(CellId(20)).unwrap().state();
        assert!(Arc::ptr_eq(
            shared,
            heating.source(CellId(10)).unwrap().state()
        ));
        let a = 20.0 * (-1.2 * time).exp();
        let retained = 20.0 * (-0.6 * time).exp();
        let decays = (20.0 - retained) * (2.0 / 3.0);
        let local_decays = [(20.0 - a) / 3.0, decays - (20.0 - a) / 3.0];
        let movement = shared.spatial().transport().unwrap();
        close(
            movement.cumulative_transfer_counts()[0][1],
            (20.0 - a) / 2.0,
            1e-8,
        );
        for (cell, (id, initial_k, capacity, parent)) in
            [(20, 300.0, 8.0, a), (10, 450.0, 6.0, retained - a)]
                .into_iter()
                .enumerate()
        {
            let d = local_decays[cell];
            let inventory = &shared.inventories()[cell];
            assert!(Arc::ptr_eq(
                inventory.model(),
                original.inventories()[cell].model()
            ));
            close(inventory.populations()[1], parent, 1e-8);
            close(inventory.populations()[2], d, 1e-8);
            close(inventory.cumulative_reaction_counts()[1], d, 1e-8);
            close(movement.cumulative_outflow_counts()[cell][1], d / 2.0, 1e-8);
            close(
                material.temperature_k(CellId(id)).unwrap(),
                initial_k + 2.0 * d / capacity,
                1e-8,
            );
            assert_eq!(
                material.cells()[&CellId(id)].mass_kg,
                initial.cells()[&CellId(id)].mass_kg
            );
            let release = shared.cell_released(cell).unwrap();
            close(release.deposited_j, 2.0 * d, 1e-8);
            close(release.radiation_j, 0.5 * d, 1e-8);
            close(release.escaped_j, 0.5 * d, 1e-8);
        }
        assert!(local_decays[1] > 0.0 && retained - a > 0.0);
        close(retained + 1.5 * decays, 20.0, 1e-8);
        close(
            shared.neutrons().neutron_counts().iter().sum(),
            0.125 * decays,
            1e-8,
        );
        close(shared.kinetic_energy_j().unwrap(), 0.5 * decays, 1e-8);
        let release = heating.cumulative_released().unwrap();
        close(material.energy_j() - initial.energy_j(), 2.0 * decays, 1e-8);
        close(
            release.total_j() + shared.kinetic_energy_j().unwrap(),
            3.5 * decays,
            1e-8,
        );
        assert_eq!(
            heating
                .last_step(CellId(20))
                .unwrap()
                .activation()
                .external_injected_neutrons,
            0.0
        );
    }

    // The accepted checkpoint still owns the finite stock and transport ledger
    // if a later material candidate is cancelled after doing part of its work.
    let accepted_state = Arc::clone(checkpoint.source(CellId(20)).unwrap().state());
    assert!(polls > 2);
    let mut remaining = polls / 2;
    assert!(
        checkpoint
            .advance_to(&middle, 1.5, controls(0.25), || {
                if remaining == 0 {
                    Err(DepositionError::Cancelled)
                } else {
                    remaining -= 1;
                    Ok(())
                }
            })
            .is_err()
    );
    assert!(Arc::ptr_eq(
        checkpoint.source(CellId(20)).unwrap().state(),
        &accepted_state
    ));
    assert_eq!(middle.time_s(), 0.5);
    assert_eq!(state.inventories()[0].populations()[1], 20.0);
    assert_eq!(
        state
            .spatial()
            .transport()
            .unwrap()
            .cumulative_transfer_counts()[0][1],
        0.0
    );
}
