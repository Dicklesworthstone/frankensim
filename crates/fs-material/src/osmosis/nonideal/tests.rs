use super::*;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}"
    );
}

fn interacting() -> VirialSolutionModel {
    VirialSolutionModel::new(
        298.15,
        18e-6,
        1.0,
        vec![200.0, 100.0],
        vec![0.001, 0.0002, 0.0002, -0.002],
    )
    .unwrap()
}

#[test]
fn g1_zero_interactions_recover_original_dilute_water_law() {
    let law = VirialSolutionModel::new(298.15, 18e-6, 1.0, vec![500.0; 2], vec![0.0; 4]).unwrap();
    for c in [[0.0, 0.0], [100.0, 50.0], [0.0, 300.0]] {
        let sample = law.evaluate(&c).unwrap();
        let ideal = super::super::ideal_dilute_osmotic_pressure_pa(298.15, &c).unwrap();
        close(sample.osmotic_pressure_pa, ideal, 1e-14);
        close(sample.osmotic_stiffness_pa, ideal, 1e-14);
        close(
            sample.water_activity,
            super::super::water_activity_from_osmotic_pressure(298.15, ideal, 18e-6).unwrap(),
            1e-14,
        );
    }
    let water = law.evaluate(&[0.0, 0.0]).unwrap();
    assert_eq!(water.osmotic_pressure_pa, 0.0);
    assert_eq!(water.water_activity, 1.0);
    assert_eq!(
        law.mixing_free_energy_density_j_m3(&[0.0, 0.0]).unwrap(),
        0.0
    );
}

#[test]
fn g1_cross_interactions_change_both_osmosis_and_solute_activities() {
    let law = interacting();
    // Bc=(0.11,-0.08), c^T Bc=7 mol/m³. The two opposite signs
    // exercise repulsive and attractive interactions on one stable branch.
    let sample = law.evaluate(&[100.0, 50.0]).unwrap();
    close(sample.osmotic_pressure_pa, R * 298.15 * 153.5, 1e-14);
    close(sample.osmotic_stiffness_pa, R * 298.15 * 157.0, 1e-14);
    close(sample.log_water_activity, -18e-6 * 153.5, 1e-14);
    let mut coefficients = [0.0; 2];
    law.log_activity_coefficients_into(&[100.0, 50.0], &mut coefficients)
        .unwrap();
    close(coefficients[0], 0.11, 1e-14);
    close(coefficients[1], -0.08, 1e-14);
    let mut chemical = [0.0; 2];
    law.chemical_potentials_into(&[100.0, 50.0], &mut chemical)
        .unwrap();
    close(chemical[0], R * 298.15 * (100.0_f64.ln() + 0.11), 1e-14);
    close(chemical[1], R * 298.15 * (50.0_f64.ln() - 0.08), 1e-14);
}

#[test]
fn g1_free_energy_volume_derivative_is_negative_osmotic_pressure() {
    let law = interacting();
    let amounts = [0.1, 0.05];
    let volume = 1e-3;
    let delta = 1e-8;
    let energy = |v: f64| {
        v * law
            .mixing_free_energy_density_j_m3(&[amounts[0] / v, amounts[1] / v])
            .unwrap()
    };
    let pressure = -(energy(volume + delta) - energy(volume - delta)) / (2.0 * delta);
    close(
        pressure,
        law.evaluate_inventories(volume, &amounts)
            .unwrap()
            .osmotic_pressure_pa,
        2e-9,
    );
    let high = law
        .evaluate_inventories(volume + delta, &amounts)
        .unwrap()
        .osmotic_pressure_pa;
    let low = law
        .evaluate_inventories(volume - delta, &amounts)
        .unwrap()
        .osmotic_pressure_pa;
    let stiffness = -(high - low) * volume / (2.0 * delta);
    close(
        stiffness,
        law.evaluate_inventories(volume, &amounts)
            .unwrap()
            .osmotic_stiffness_pa,
        2e-9,
    );
}

#[test]
fn g1_chemical_potential_is_solute_free_energy_derivative() {
    let law = interacting();
    let c = [100.0, 50.0];
    let mut chemical = [0.0; 2];
    law.chemical_potentials_into(&c, &mut chemical).unwrap();
    for i in 0..2 {
        let mut high = c;
        let mut low = c;
        let delta = 1e-3;
        high[i] += delta;
        low[i] -= delta;
        let slope = (law.mixing_free_energy_density_j_m3(&high).unwrap()
            - law.mixing_free_energy_density_j_m3(&low).unwrap())
            / (2.0 * delta);
        close(slope, chemical[i], 2e-9);
    }
}

#[test]
fn g1_matched_jacobian_satisfies_gibbs_duhem_and_positive_curvature() {
    let law = interacting();
    let c = [100.0, 50.0];
    let mut hessian = [0.0; 4];
    law.chemical_potential_jacobian_into(&c, &mut hessian)
        .unwrap();
    assert_eq!(hessian[1], hessian[2]);
    assert!(hessian[0] > 0.0 && hessian[3] > 0.0);
    assert!(hessian[0] * hessian[3] > hessian[1] * hessian[2]);
    for i in 0..2 {
        let mut high = c;
        let mut low = c;
        let delta = 1e-3;
        high[i] += delta;
        low[i] -= delta;
        let gradient = (law.evaluate(&high).unwrap().osmotic_pressure_pa
            - law.evaluate(&low).unwrap().osmotic_pressure_pa)
            / (2.0 * delta);
        close(gradient, c[0] * hessian[i] + c[1] * hessian[2 + i], 2e-10);
    }
}

#[test]
fn g0_stability_is_a_coupled_domain_condition() {
    // Positive diagonal entries alone do not make this Hessian positive.
    assert!(matches!(
        VirialSolutionModel::new(
            300.0,
            18e-6,
            1.0,
            vec![1000.0; 2],
            vec![0.0, 0.002, 0.002, 0.0]
        ),
        Err(NonidealSolutionError::UnstableDomain { .. })
    ));
    assert!(matches!(
        VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![1000.0], vec![-0.002]),
        Err(NonidealSolutionError::UnstableDomain { .. })
    ));
    assert!(VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![500.0], vec![-0.001]).is_ok());
    assert!(
        VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![1.0], vec![-1.0 + f64::EPSILON]).is_err()
    );
}

#[test]
fn g0_source_box_is_exact_and_never_extrapolated() {
    let law = interacting();
    assert!(law.evaluate(&[200.0, 100.0]).is_ok());
    assert!(matches!(
        law.evaluate(&[200.0_f64.next_up(), 0.0]),
        Err(NonidealSolutionError::OutsideDomain { species: 0, .. })
    ));
    assert!(law.evaluate(&[f64::NAN, 0.0]).is_err());
    assert!(law.evaluate(&[-1.0, 0.0]).is_err());
    assert!(law.evaluate(&[1.0]).is_err());
    assert!(law.evaluate_inventories(0.0, &[1.0, 1.0]).is_err());
    assert!(law.evaluate_inventories(1.0, &[-1.0, 1.0]).is_err());
    assert!(law.evaluate_inventories(1e-4, &[0.1, 0.05]).is_err());
}

#[test]
fn g0_singular_or_invalid_chemical_outputs_are_failure_atomic() {
    let law = interacting();
    let mut chemical = [123.0, 456.0];
    assert!(matches!(
        law.chemical_potentials_into(&[1.0, 0.0], &mut chemical),
        Err(NonidealSolutionError::ZeroConcentration { species: 1 })
    ));
    assert_eq!(chemical, [123.0, 456.0]);
    assert!(
        law.log_activity_coefficients_into(&[1.0, 101.0], &mut chemical)
            .is_err()
    );
    assert_eq!(chemical, [123.0, 456.0]);
    let mut hessian = [789.0; 4];
    assert!(
        law.chemical_potential_jacobian_into(&[0.0, 1.0], &mut hessian)
            .is_err()
    );
    assert_eq!(hessian, [789.0; 4]);
    let mut short = [123.0];
    assert!(
        law.log_activity_coefficients_into(&[1.0, 1.0], &mut short)
            .is_err()
    );
    assert_eq!(short, [123.0]);
    // Only the logarithmic/reciprocal APIs are singular: absent solutes are
    // otherwise genuine supported physical endpoints.
    law.log_activity_coefficients_into(&[0.0, 0.0], &mut chemical)
        .unwrap();
    assert_eq!(chemical, [0.0, 0.0]);
}

#[test]
fn g3_species_permutation_preserves_water_thermodynamics() {
    let original = interacting();
    let permuted = VirialSolutionModel::new(
        298.15,
        18e-6,
        1.0,
        vec![100.0, 200.0],
        vec![-0.002, 0.0002, 0.0002, 0.001],
    )
    .unwrap();
    assert_eq!(
        original.evaluate(&[100.0, 50.0]).unwrap(),
        permuted.evaluate(&[50.0, 100.0]).unwrap()
    );
}

#[test]
fn g0_malformed_models_and_unrepresentable_activity_refuse() {
    assert!(
        VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![1.0, 1.0], vec![0.0, 0.1, 0.2, 0.0])
            .is_err()
    );
    assert!(VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![1.0], vec![]).is_err());
    assert!(VirialSolutionModel::new(0.0, 18e-6, 1.0, vec![1.0], vec![0.0]).is_err());
    assert!(VirialSolutionModel::new(300.0, 18e-6, 0.0, vec![1.0], vec![0.0]).is_err());
    assert!(
        VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![1.0; 65], vec![0.0; 65 * 65]).is_err()
    );
    let law = VirialSolutionModel::new(300.0, 1.0, 1.0, vec![1000.0], vec![0.0]).unwrap();
    assert!(matches!(
        law.evaluate(&[1000.0]),
        Err(NonidealSolutionError::Arithmetic(_))
    ));
    let pure = VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![], vec![]).unwrap();
    assert_eq!(pure.evaluate(&[]).unwrap().water_activity, 1.0);
    assert_eq!(pure.heap_bytes(), 0);
}

#[test]
fn g3_reference_changes_chemical_zero_without_changing_water() {
    let first = interacting();
    let second = VirialSolutionModel::new(
        298.15,
        18e-6,
        1000.0,
        first.maximum_concentrations_mol_m3().to_vec(),
        first.interaction_matrix_m3_mol().to_vec(),
    )
    .unwrap();
    let c = [100.0, 50.0];
    assert_eq!(first.evaluate(&c).unwrap(), second.evaluate(&c).unwrap());
    let mut a = [0.0; 2];
    let mut b = [0.0; 2];
    first.chemical_potentials_into(&c, &mut a).unwrap();
    second.chemical_potentials_into(&c, &mut b).unwrap();
    for i in 0..2 {
        close(a[i] - b[i], R * 298.15 * 1000.0_f64.ln(), 1e-14);
    }
}

#[test]
fn g1_full_species_capacity_retains_all_cross_interactions() {
    // A positive rank-one interaction exercises the full 64-species owner,
    // including the existing blocked Cholesky's native GEMM trailing update.
    let n = MAX_OSMOTIC_SOLUTES;
    let law =
        VirialSolutionModel::new(300.0, 18e-6, 1.0, vec![100.0; n], vec![1e-5; n * n]).unwrap();
    let sample = law.evaluate(&[10.0; MAX_OSMOTIC_SOLUTES]).unwrap();
    let total = 10.0 * n as f64;
    close(
        sample.osmotic_pressure_pa,
        R * 300.0 * (total + 0.5e-5 * total * total),
        1e-14,
    );
    assert_eq!(law.heap_bytes(), (n + n * n) * core::mem::size_of::<f64>());
}
