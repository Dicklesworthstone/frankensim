//! Finite delayed-family storage and release in the original material transaction.

use super::*;
use crate::fusion_neutrons::activation::{material_precursor_fixture, precursor_expected_history};

#[test]
fn g1_g3_g4_delayed_family_energy_reaches_material_at_its_actual_capture_time() {
    let state = Arc::new(material_precursor_fixture(100.0, 4.0).unwrap());
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
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
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
    for (time, material, heating) in [(0.5, &middle, &checkpoint), (2.0, &warm, &accepted)] {
        let native = heating.source(CellId(20)).unwrap().state();
        let (parent, decays, precursors, neutrons) = precursor_expected_history(100.0, time);
        let captures = 1.5 * decays - precursors - neutrons;
        let heat = decays + 0.5 * captures;
        close(native.inventories()[0].populations()[0], parent, 1e-8);
        close(native.neutrons().precursor_counts()[0], precursors, 1e-8);
        close(native.neutrons().neutron_counts()[0], neutrons, 1e-8);
        close(native.precursor_energy_j().unwrap(), 2.0 * precursors, 1e-8);
        let released = native.cell_released(0).unwrap();
        close(released.deposited_j, heat, 1e-8);
        close(released.radiation_j, decays + 1.5 * captures, 1e-8);
        close(released.escaped_j, 0.5 * decays, 1e-8);
        close(
            native.stored_energy_j().unwrap() + released.total_j(),
            5.5 * decays,
            1e-8,
        );
        close(material.energy_j() - initial.energy_j(), heat, 1e-8);
        close(
            material.temperature_k(CellId(20)).unwrap(),
            300.0 + heat / 8.0,
            1e-8,
        );
        assert_eq!(material.cells()[&CellId(20)].mass_kg, 2.0);
        let step = heating.last_step(CellId(20)).unwrap().activation();
        assert_eq!(step.external_injected_neutrons, 0.0);
        assert_eq!(step.external_injected_precursors, 0.0);
        assert!(step.reserved_decay_precursor_energy_j > 0.0);
        assert!(step.kinetic.precursor_conversion_energy_j > 0.0);
    }
    let before = Arc::clone(checkpoint.source(CellId(20)).unwrap().state());
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
    assert_eq!(state.inventories()[0].populations()[0], 100.0);
    let (retry, _, _) = checkpoint
        .advance_to(&middle, 2.0, controls(0.125), || Ok(()))
        .unwrap();
    close(retry.energy_j(), warm.energy_j(), 1e-9);
}

// Independent closed-form recurrence for the existing interval-mean finite
// fission closure. Continuous Bateman refinement is tested by the native owner;
// this material check follows exactly the accepted fixed-duration intervals.
fn induced_interval_history(initial_parent: f64, steps: usize, duration: f64) -> [f64; 3] {
    let mut parent = initial_parent;
    let mut precursor = 0.0;
    let mut slow = 0.0;
    let precursor_survival = (-duration).exp();
    let neutron_survival = (-0.3 * duration).exp();
    let residence = (precursor_survival - neutron_survival) / (0.3 - 1.0);
    for _ in 0..steps {
        let next_parent = parent * (-0.2 * duration).exp();
        let family_rate = 0.5 * (parent - next_parent) / duration;
        slow = slow * neutron_survival
            + precursor * residence
            + family_rate * ((1.0 - neutron_survival) / 0.3 - residence);
        precursor = precursor * precursor_survival + family_rate * (1.0 - precursor_survival);
        parent = next_parent;
    }
    [parent, precursor, slow]
}

#[test]
fn g1_g3_g4_finite_induced_fission_reserves_energy_and_heats_original_material_once() {
    use crate::fusion_neutrons::activation::material_induced_precursor_fixture;

    let initial_parent = 100.0;
    let state = Arc::new(material_induced_precursor_fixture(initial_parent));
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
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
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

    for (steps, material, owner) in [(4, &middle, &checkpoint), (16, &warm, &accepted)] {
        let native = owner.source(CellId(20)).unwrap().state();
        let [parent, family, slow] = induced_interval_history(initial_parent, steps, 0.125);
        let fissions = initial_parent - parent;
        let captures = 0.5 * fissions - family - slow;
        let heat = fissions + 0.5 * captures;
        let radiation = 3.0 * fissions + 0.5 * captures;
        close(native.inventories()[0].populations()[0], parent, 1e-8);
        assert_eq!(&native.inventories()[0].populations()[1..], &[0.0; 3]);
        close(native.neutrons().neutron_counts()[1], 10.0, 1e-8);
        close(native.neutrons().neutron_counts()[0], slow, 1e-8);
        close(native.neutrons().precursor_counts()[0], family, 1e-8);
        close(native.precursor_energy_j().unwrap(), family, 1e-8);
        close(
            native.cumulative_induced_precursor_energy_j().unwrap(),
            0.5 * fissions,
            1e-8,
        );
        let released = native.cell_released(0).unwrap();
        close(released.deposited_j, heat, 1e-8);
        close(released.radiation_j, radiation, 1e-8);
        assert_eq!(released.escaped_j, 0.0);
        close(
            native.stored_energy_j().unwrap() + released.total_j(),
            20.0 + 4.5 * fissions,
            1e-8,
        );
        // The resumed thermal audit begins at the accepted middle checkpoint.
        let prior_fissions = if steps > 4 {
            checkpoint.audit().fissions
        } else {
            0.0
        };
        close(prior_fissions + owner.audit().fissions, fissions, 1e-8);
        close(material.energy_j() - initial.energy_j(), heat, 1e-8);
        close(
            material.temperature_k(CellId(20)).unwrap(),
            300.0 + heat / 8.0,
            1e-8,
        );
        assert_eq!(material.cells()[&CellId(20)].mass_kg, 2.0);
        let last = owner.last_step(CellId(20)).unwrap().activation();
        assert_eq!(last.external_injected_neutrons, 0.0);
        assert_eq!(last.external_injected_precursors, 0.0);
        assert!(last.kinetic.induced_precursor_energy_j > 0.0);
        close(
            last.nuclear_generated_energy_j,
            last.kinetic.nuclear_reaction_energy_j + last.kinetic.induced_precursor_energy_j,
            1e-9,
        );
    }

    let rebound_sources = FusionNeutronActivationCellSource::bind_without_external_source(
        Arc::clone(checkpoint.source(CellId(20)).unwrap().state()),
        Time::new(middle.time_s()),
        vec![(CellId(20), MaterialId(1))],
    )
    .unwrap();
    let rebound =
        FusionNeutronActivationHeating::new(&middle, rebound_sources, activation_budget()).unwrap();
    let (resumed, _, _) = rebound
        .advance_to(&middle, 2.0, controls(0.125), || Ok(()))
        .unwrap();
    close(resumed.energy_j(), warm.energy_j(), 1e-9);
    assert!(polls > 2);
    let before = Arc::clone(checkpoint.source(CellId(20)).unwrap().state());
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
    assert_eq!(state.inventories()[0].populations()[0], initial_parent);
    assert_eq!(middle.time_s(), 0.5);
    let (retry, _, _) = checkpoint
        .advance_to(&middle, 2.0, controls(0.125), || Ok(()))
        .unwrap();
    close(retry.energy_j(), warm.energy_j(), 1e-9);
}
