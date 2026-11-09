//! Original material temperature and heat with finite radioactive neutron stock.

use super::*;
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
