//! Evaluated MF1/MF5 prompt SF neutrons and MF8 product decay on native material.
//! The shared manufactured records test source semantics, not isotope calibration.

use super::*;
use crate::endf::decay::model::burnup::neutrons::material_fixture;

const EV_J: f64 = 1.602_176_634e-19;

fn macroscopic_budget() -> FusionNeutronActivationBudget {
    let mut budget = activation_budget();
    // Extensive count allowances scale with the 1e18-nucleus specimen. Energy
    // stays in joules with the original small absolute calorimetric allowance.
    budget.spatial.diffusion.population_tolerance = 1e8;
    budget.spatial.burnup.nuclei_tolerance = 1e8;
    budget
}

#[test]
fn g1_g3_g4_evaluated_sf_and_product_neutrons_heat_original_material_once() {
    let initial_nuclei = 1e18;
    let (_evaluated, state) = material_fixture(initial_nuclei);
    let state = Arc::new(state);
    let initial = thermal(
        vec![cell(20, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let sources = FusionNeutronActivationCellSource::bind_without_external_source(
        Arc::clone(&state),
        Time::new(0.0),
        vec![(CellId(20), MaterialId(1))],
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, macroscopic_budget()).unwrap();
    let (middle, checkpoint, _) = heating
        .advance_to(&initial, 0.5, controls(0.125), || Ok(()))
        .unwrap();
    let mut polls = 0;
    let (warm, accepted, _) = checkpoint
        .advance_to(&middle, 2.0, controls(0.125), || {
            polls += 1;
            Ok(())
        })
        .unwrap();

    for (time, material, heating) in [(0.5_f64, &middle, &checkpoint), (2.0, &warm, &accepted)] {
        let shared = heating.source(CellId(20)).unwrap().state();
        let inventory = &shared.inventories()[0];
        let populations = inventory.populations();
        let lambda = core::f64::consts::LN_2 / 2.0;
        let capture_rate = 0.3;
        let surviving_parent = initial_nuclei * (-lambda * time).exp();
        let parent_decays = initial_nuclei * -(-lambda * time).exp_m1();
        let fragment =
            0.25 * initial_nuclei * ((-lambda * time).exp() - (-2.0 * lambda * time).exp());
        let fragment_decays = 0.25 * parent_decays - fragment;
        // Solve n' = lambda*P0*exp(-lambda*t)
        //           - .5*lambda*P0*exp(-2*lambda*t) - .3*n independently.
        let neutrons = lambda
            * initial_nuclei
            * (((-lambda * time).exp() - (-capture_rate * time).exp()) / (capture_rate - lambda)
                - 0.5 * ((-2.0 * lambda * time).exp() - (-capture_rate * time).exp())
                    / (capture_rate - 2.0 * lambda));
        let captures = 0.5 * parent_decays + fragment_decays - neutrons;
        assert!(captures > 0.0 && fragment_decays > 0.0);
        for (&actual, expected) in populations.iter().zip([
            surviving_parent,
            0.75 * parent_decays,
            fragment,
            fragment_decays,
            0.25 * parent_decays,
        ]) {
            close(actual, expected, 1e7);
        }
        for &count in shared.neutrons().neutron_counts() {
            close(count, 0.5 * neutrons, 1e7);
        }
        let kinetic = shared.cell_kinetic_released(0).unwrap();
        close(kinetic.deposited_j, 0.5 * captures * EV_J, 1e-9);
        close(kinetic.radiation_j, 1.5 * captures * EV_J, 1e-9);
        let heat = (4.0 * parent_decays + fragment_decays + 0.5 * captures) * EV_J;
        let radiation = (4.0 * parent_decays + fragment_decays + 1.5 * captures) * EV_J;
        let released = shared.cell_released(0).unwrap();
        close(released.deposited_j, heat, 1e-9);
        close(released.radiation_j, radiation, 1e-9);
        assert_eq!(released.escaped_j, 0.0);
        close(
            shared.kinetic_energy_j().unwrap(),
            2.0 * neutrons * EV_J,
            1e-9,
        );
        close(
            released.total_j() + shared.kinetic_energy_j().unwrap(),
            (9.0 * parent_decays + 4.0 * fragment_decays) * EV_J,
            1e-9,
        );
        close(material.energy_j() - initial.energy_j(), heat, 1e-8);
        close(
            material.temperature_k(CellId(20)).unwrap(),
            300.0 + heat / 8.0,
            1e-8,
        );
        assert_eq!(material.cells()[&CellId(20)].mass_kg, 2.0);
        assert!(Arc::ptr_eq(
            inventory.model(),
            state.inventories()[0].model()
        ));
        let last = heating.last_step(CellId(20)).unwrap().activation();
        assert_eq!(last.external_injected_neutrons, 0.0);
        let events = &last.spatial.burnup[0].reaction_counts;
        close(last.decay_neutrons, 2.0 * events[0] + events[2], 1e7);
        assert_eq!(last.kinetic.diffusion, last.spatial.neutrons);
    }

    let before = Arc::clone(checkpoint.source(CellId(20)).unwrap().state());
    assert!(polls > 2);
    let mut remaining = polls - 2;
    assert!(
        checkpoint
            .advance_to(&middle, 2.0, controls(0.125), || {
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
        &before
    ));
    assert_eq!(middle.time_s(), 0.5);
    assert_eq!(state.inventories()[0].populations()[0], initial_nuclei);
    let (retry, _, _) = checkpoint
        .advance_to(&middle, 2.0, controls(0.125), || Ok(()))
        .unwrap();
    close(retry.energy_j(), warm.energy_j(), 1e-9);
}
