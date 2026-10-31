use super::*;
use crate::osmosis::OsmoticPressureVolume;
use fs_material::osmosis::IdealDiluteSolution;

fn bulk(concentrations: &[f64], volume: f64, pressure: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(
            298.15,
            18e-6,
            volume,
            concentrations.iter().map(|c| c * volume).collect(),
        )
        .unwrap(),
        OsmoticPressureVolume::new(volume, pressure, 0.0).unwrap(),
        volume * 1e-6,
    )
    .unwrap()
}
fn membrane(lp: f64, sigma: f64, permeability: f64, count: usize) -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(lp).unwrap(),
        vec![
            (
                ReflectionCoefficient::new(sigma).unwrap(),
                SolutePermeance::new(permeability).unwrap()
            );
            count
        ],
    )
    .unwrap()
}
fn film(cells: usize, concentration: f64) -> TransientFilm {
    TransientFilm::new(1e-4, cells, vec![1e-9], vec![concentration; cells]).unwrap()
}
fn retaining_pair(cells: usize, volume: f64, right_pressure: f64) -> TransientOsmoticPair {
    TransientOsmoticPair::new(
        [
            bulk(&[50.0], volume, 101_325.0),
            bulk(&[100.0], volume, right_pressure),
        ],
        1e-4,
        membrane(2e-11, 1.0, 0.0, 1),
        Some(film(cells, 50.0)),
        Some(film(cells, 100.0)),
        1 << 24,
    )
    .unwrap()
}
fn close(left: f64, right: f64, relative: f64, absolute: f64) {
    assert!(
        (left - right).abs() <= absolute + relative * left.abs().max(right.abs()),
        "{left:.16e} != {right:.16e}"
    );
}
fn snapshot(pair: &TransientOsmoticPair) -> Vec<u64> {
    let mut values = vec![pair.time_s().to_bits()];
    for b in pair.bulks() {
        values.push(b.solution().solvent_volume_m3().to_bits());
        values.extend(b.solution().solute_moles().iter().map(|n| n.to_bits()));
    }
    for side in 0..2 {
        for cell in 0..pair.film_cell_count(side).unwrap() {
            values.extend(
                pair.film_solute_moles(side, cell)
                    .unwrap()
                    .iter()
                    .map(|n| n.to_bits()),
            );
        }
    }
    values
}

#[test]
fn explicit_geometry_layout_and_workspace_admission() {
    assert!(TransientFilm::new(0.0, 2, vec![1e-9], vec![1.0; 2]).is_err());
    assert!(TransientFilm::new(1e-4, 0, vec![1e-9], vec![]).is_err());
    assert!(TransientFilm::new(1e-4, 2, vec![0.0], vec![1.0; 2]).is_err());
    assert!(TransientFilm::new(1e-4, 2, vec![1e-9], vec![1.0]).is_err());
    let pair = retaining_pair(4, 1e-6, 101_325.0);
    assert_eq!(pair.film_cell_count(0), Some(4));
    assert_eq!(pair.film_cell_count(2), None);
    assert!(pair.film_solute_moles(1, 4).is_none());
    close(pair.total_solvent_volume_m3(), 2.02e-6, 1e-15, 0.0);
    let mut total = [0.0];
    pair.total_solute_moles_into(&mut total).unwrap();
    close(total[0], 151.5e-6, 1e-15, 0.0);
    assert!(pair.total_solute_moles_into(&mut []).is_err());
    let refused = TransientOsmoticPair::new(
        [bulk(&[50.0], 1e-6, 0.0), bulk(&[100.0], 1e-6, 0.0)],
        1e-4,
        membrane(2e-11, 1.0, 0.0, 1),
        Some(film(4, 50.0)),
        Some(film(4, 100.0)),
        pair.workspace_bytes() - 1,
    );
    assert!(matches!(
        refused,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Workspace { .. }
        ))
    ));
}

#[test]
fn pure_water_counts_film_storage_and_exact_pressure_work() {
    let empty_film = || TransientFilm::new(1e-3, 3, vec![], vec![]).unwrap();
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[], 1e-3, 200_000.0), bulk(&[], 1e-3, 100_000.0)],
        1e-3,
        membrane(1e-10, 0.0, 0.0, 0),
        Some(empty_film()),
        Some(empty_film()),
        1 << 20,
    )
    .unwrap();
    let initial_water = pair.total_solvent_volume_m3();
    close(initial_water, 0.002002, 1e-15, 0.0);
    let report = pair
        .step(10.0, TransientStepControls::default(), || false)
        .unwrap();
    close(report.solvent_transfer_m3, 1e-7, 1e-12, 1e-20);
    close(report.pressure_work_into_solution_j, 0.01, 1e-12, 1e-16);
    close(pair.total_solvent_volume_m3(), initial_water, 1e-15, 0.0);
    close(pair.film_cell_volume_m3(0).unwrap(), 1e-6 / 3.0, 1e-15, 0.0);
    assert_eq!(report.species_count, 0);
}

#[test]
fn no_film_limit_has_exact_two_compartment_backward_euler_diffusion() {
    let mut pair = TransientOsmoticPair::new(
        [
            bulk(&[0.0], 1e-6, 101_325.0),
            bulk(&[100.0], 1e-6, 101_325.0),
        ],
        1e-4,
        membrane(0.0, 0.0, 2e-6, 1),
        None,
        None,
        1 << 20,
    )
    .unwrap();
    let report = pair
        .step(1000.0, TransientStepControls::default(), || false)
        .unwrap();
    let expected_left = 0.5e-4 * (1.0 - 1.0 / 1.4);
    close(
        pair.bulks()[0].solution().solute_moles()[0],
        expected_left,
        1e-14,
        1e-20,
    );
    assert_eq!(report.solvent_transfer_m3, 0.0);
    assert_eq!(pair.film_cell_count(0), Some(0));
    close(report.solute_balance_residual_mol[0], 0.0, 0.0, 1e-19);
}

#[test]
fn retained_film_inventory_diffuses_into_finite_bulk_without_leaking_membrane() {
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[0.0], 1e-6, 0.0), bulk(&[0.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 1.0, 0.0, 1),
        Some(film(8, 100.0)),
        Some(film(8, 0.0)),
        1 << 20,
    )
    .unwrap();
    let mut before = [0.0];
    pair.total_solute_moles_into(&mut before).unwrap();
    close(before[0], 1e-6, 1e-15, 0.0);
    let report = pair
        .step(1000.0, TransientStepControls::default(), || false)
        .unwrap();
    let mut after = [0.0];
    pair.total_solute_moles_into(&mut after).unwrap();
    close(after[0], before[0], 1e-11, 1e-18);
    assert!(pair.bulks()[0].solution().solute_moles()[0] > 0.9e-6);
    assert_eq!(pair.bulks()[1].solution().solute_moles()[0], 0.0);
    for cell in 0..8 {
        assert_eq!(pair.film_solute_moles(1, cell).unwrap()[0], 0.0);
    }
    assert_eq!(report.solvent_transfer_m3, 0.0);
}

#[test]
fn permeating_species_with_distinct_diffusivities_conserve_all_bulk_and_film_moles() {
    let law = MembraneTransport::new(
        HydraulicPermeance::new(2e-11).unwrap(),
        vec![
            (
                ReflectionCoefficient::new(0.7).unwrap(),
                SolutePermeance::new(2e-7).unwrap(),
            ),
            (
                ReflectionCoefficient::new(0.1).unwrap(),
                SolutePermeance::new(4e-7).unwrap(),
            ),
        ],
    )
    .unwrap();
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[80.0, 0.0], 1e-6, 0.0), bulk(&[0.0, 30.0], 1e-6, 0.0)],
        1e-4,
        law,
        Some(TransientFilm::new(1e-4, 4, vec![1e-9, 5e-10], [80.0, 0.0].repeat(4)).unwrap()),
        Some(TransientFilm::new(2e-4, 3, vec![1e-9, 5e-10], [0.0, 30.0].repeat(3)).unwrap()),
        1 << 20,
    )
    .unwrap();
    let mut before = [0.0; 2];
    pair.total_solute_moles_into(&mut before).unwrap();
    let water = pair.total_solvent_volume_m3();
    let report = pair
        .step(100.0, TransientStepControls::default(), || false)
        .unwrap();
    let mut after = [0.0; 2];
    pair.total_solute_moles_into(&mut after).unwrap();
    for i in 0..2 {
        close(before[i], after[i], 1e-12, 1e-19);
    }
    assert!(pair.bulks()[0].solution().solute_moles()[1] > 0.0);
    assert!(pair.bulks()[1].solution().solute_moles()[0] > 0.0);
    assert_eq!(report.species_count, 2);
    close(pair.total_solvent_volume_m3(), water, 1e-14, 0.0);
    for side in 0..2 {
        for cell in 0..pair.film_cell_count(side).unwrap() {
            assert!(
                pair.film_solute_moles(side, cell)
                    .unwrap()
                    .iter()
                    .all(|n| n.is_finite() && *n >= 0.0)
            );
        }
    }
}

#[test]
fn nonideal_owner_is_refused_by_neutral_ideal_film_transport() {
    use fs_material::osmosis::{AqueousSolution, nonideal::VirialSolutionModel};
    let model = std::sync::Arc::new(
        VirialSolutionModel::new(298.15, 18e-6, 1.0, vec![200.0], vec![0.002]).unwrap(),
    );
    let left = OsmoticCompartment::new(
        AqueousSolution::nonideal(model, 1e-6, vec![50e-6]).unwrap(),
        OsmoticPressureVolume::new(1e-6, 0.0, 0.0).unwrap(),
        1e-9,
    )
    .unwrap();
    assert!(
        TransientOsmoticPair::new(
            [left, bulk(&[100.0], 1e-6, 0.0)],
            1e-4,
            membrane(2e-11, 1.0, 0.0, 1),
            Some(film(4, 50.0)),
            Some(film(4, 100.0)),
            1 << 20,
        )
        .is_err()
    );
}

#[test]
fn transient_profile_relaxes_to_original_steady_film_law_in_both_directions() {
    for pressure_right in [101_325.0, 501_325.0] {
        // Large but finite bulk volumes separate film relaxation from changing
        // bulk composition. Their water and moles still enter every balance.
        let mut pair = retaining_pair(8, 1.0, pressure_right);
        let law = PolarizedMembraneTransport::new(
            membrane(2e-11, 1.0, 0.0, 1),
            vec![Some(
                FilmPermeance::from_diffusion_layer(1e-9, 1e-4).unwrap(),
            )],
            vec![Some(
                FilmPermeance::from_diffusion_layer(1e-9, 1e-4).unwrap(),
            )],
        )
        .unwrap();
        let mut total_before = [0.0];
        pair.total_solute_moles_into(&mut total_before).unwrap();
        let mut report = pair
            .step(100.0, TransientStepControls::default(), || false)
            .unwrap();
        for _ in 0..9 {
            report = pair
                .step(100.0, TransientStepControls::default(), || false)
                .unwrap();
        }
        let bulks = pair.bulks();
        let cl = bulks[0].solution().solute_moles()[0] / bulks[0].solution().solvent_volume_m3();
        let cr = bulks[1].solution().solute_moles()[0] / bulks[1].solution().solvent_volume_m3();
        let (fl, fr) = law
            .coefficients(0, report.volume_flux_m_s)
            .unwrap()
            .membrane_face_concentrations(cl, cr)
            .unwrap();
        close(report.membrane_left_concentration_mol_m3[0], fl, 2e-6, 1e-8);
        close(
            report.membrane_right_concentration_mol_m3[0],
            fr,
            2e-6,
            1e-8,
        );
        let mut total_after = [0.0];
        pair.total_solute_moles_into(&mut total_after).unwrap();
        close(total_before[0], total_after[0], 2e-13, 1e-15);
        if pressure_right > 200_000.0 {
            assert!(report.volume_flux_m_s < 0.0);
            assert!(report.pressure_work_into_solution_j > 0.0);
            assert!(fl < cl && fr > cr);
        } else {
            assert!(report.volume_flux_m_s > 0.0);
            assert!(fl > cl && fr < cr);
        }
    }
}

fn uptake(cells: usize, step: f64, duration: f64) -> f64 {
    let mut pair = retaining_pair(cells, 1e-6, 101_325.0);
    let controls = TransientStepControls {
        absolute_water_tolerance_m3: 1e-22,
        relative_root_tolerance: 1e-11,
        ..TransientStepControls::default()
    };
    for _ in 0..(duration / step).round() as usize {
        pair.step(step, controls, || false).unwrap();
    }
    pair.bulks()[1].solution().solvent_volume_m3() - 1e-6
}

#[test]
fn time_and_space_refinement_resolve_film_startup() {
    let reference_time = uptake(8, 0.03125, 4.0);
    let coarse_time = uptake(8, 0.5, 4.0);
    let fine_time = uptake(8, 0.25, 4.0);
    assert!((fine_time - reference_time).abs() < 0.7 * (coarse_time - reference_time).abs());
    let reference_space = uptake(32, 0.03125, 4.0);
    let coarse_space = uptake(2, 0.03125, 4.0);
    let fine_space = uptake(4, 0.03125, 4.0);
    assert!((fine_space - reference_space).abs() < 0.7 * (coarse_space - reference_space).abs());
}

#[test]
fn steep_water_law_accepts_a_resolved_transfer_bracket_without_hiding_its_residual() {
    // Deliberately extreme numerical stiffness: the source law remains intact,
    // while the root's transfer-location tolerance can outresolve its residual.
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[50.0], 1e-6, 0.0), bulk(&[100.0], 1e-6, 0.0)],
        1e-4,
        membrane(1e6, 1.0, 0.0, 1),
        Some(film(8, 50.0)),
        Some(film(8, 100.0)),
        1 << 20,
    )
    .unwrap();
    let report = pair
        .step(
            10.0,
            TransientStepControls {
                absolute_water_tolerance_m3: 1e-24,
                relative_root_tolerance: 1e-13,
                ..TransientStepControls::default()
            },
            || false,
        )
        .unwrap();
    assert!(report.solvent_transfer_m3 > 0.0);
    assert!(report.water_root_residual_m3.is_finite());
    assert!(report.water_root_residual_m3.abs() > 1e-24 + 1e-13 * report.solvent_transfer_m3);
    close(report.solute_balance_residual_mol[0], 0.0, 0.0, 1e-18);
}

#[test]
fn cancellation_and_root_exhaustion_preserve_all_finite_inventories_and_clock() {
    let mut pair = retaining_pair(8, 1e-6, 101_325.0);
    let before = snapshot(&pair);
    let mut polls = 0;
    let result = pair.step(10.0, TransientStepControls::default(), || {
        polls += 1;
        polls == 10
    });
    assert!(matches!(
        result,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Cancelled
        ))
    ));
    assert_eq!(snapshot(&pair), before);
    let result = pair.step(
        10.0,
        TransientStepControls {
            maximum_root_evaluations: 2,
            ..TransientStepControls::default()
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::NonlinearConvergence
        ))
    ));
    assert_eq!(snapshot(&pair), before);
    pair.step(10.0, TransientStepControls::default(), || false)
        .unwrap();
    assert_eq!(pair.time_s(), 10.0);
}
