//! Evaluated MF1/MF5 prompt SF neutrons and MF8 product decay on native material.
//! The shared manufactured records test source semantics, not isotope calibration.

use super::*;
use crate::endf::decay::model::burnup::neutrons::{
    EvaluatedBurnupNeutrons, material_fixture, material_fixture_evaporation,
    material_fixture_general_evaporation, material_fixture_maxwellian,
};

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
    verify_material_history(material_fixture, 2.0);
}

#[test]
fn g1_g3_g4_lf5_spectrum_heats_original_material_with_the_same_decay_history() {
    verify_material_history(material_fixture_general_evaporation, 2.0);
}

#[test]
fn g1_g3_g4_analytic_spectra_heat_material_from_their_original_energy_moments() {
    // Independent 70-digit incomplete-gamma means for theta=2 eV, R=3 eV.
    verify_material_history(material_fixture_maxwellian, 1.4794210037284057);
    verify_material_history(material_fixture_evaporation, 1.7292098607415109);
}

fn verify_material_history(
    fixture: fn(f64) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState),
    prompt_mean_ev: f64,
) {
    let initial_nuclei = 1e18;
    let (_evaluated, state) = fixture(initial_nuclei);
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
        // Separate prompt SF births from the fragment's delayed births in
        // the independent analytic solution. Their neutron means differ.
        let first =
            ((-lambda * time).exp() - (-capture_rate * time).exp()) / (capture_rate - lambda);
        let second = ((-2.0 * lambda * time).exp() - (-capture_rate * time).exp())
            / (capture_rate - 2.0 * lambda);
        let prompt = 0.5 * lambda * initial_nuclei * first;
        let delayed = 0.5 * lambda * initial_nuclei * (first - second);
        let neutrons = prompt + delayed;
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
        close(
            shared.neutrons().neutron_counts()[0],
            prompt * (3.0 - prompt_mean_ev) / 2.0 + delayed / 2.0,
            1e7,
        );
        close(
            shared.neutrons().neutron_counts()[1],
            prompt * (prompt_mean_ev - 1.0) / 2.0 + delayed / 2.0,
            1e7,
        );
        let emitted_ev = 0.5 * prompt_mean_ev * parent_decays + 2.0 * fragment_decays;
        let retained_ev = prompt_mean_ev * prompt + 2.0 * delayed;
        let captured_ev = emitted_ev - retained_ev;
        let kinetic = shared.cell_kinetic_released(0).unwrap();
        close(kinetic.deposited_j, 0.25 * captured_ev * EV_J, 1e-9);
        close(kinetic.radiation_j, 0.75 * captured_ev * EV_J, 1e-9);
        let heat = (4.0 * parent_decays + fragment_decays + 0.25 * captured_ev) * EV_J;
        let radiation =
            (5.0 * parent_decays + 3.0 * fragment_decays - emitted_ev + 0.75 * captured_ev) * EV_J;
        let released = shared.cell_released(0).unwrap();
        close(released.deposited_j, heat, 1e-9);
        close(released.radiation_j, radiation, 1e-9);
        assert_eq!(released.escaped_j, 0.0);
        close(shared.kinetic_energy_j().unwrap(), retained_ev * EV_J, 1e-9);
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
