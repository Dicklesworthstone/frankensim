use super::*;

fn model() -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(298.15, 18e-6, 1.0, vec![200.0], vec![0.002]).unwrap())
}

#[test]
fn g0_ideal_conversion_preserves_exact_legacy_behavior() {
    let original = IdealDiluteSolution::new(298.15, 18e-6, 1e-3, vec![0.1, 0.05]).unwrap();
    let solution = AqueousSolution::from(original.clone());
    assert!(solution.is_ideal());
    assert_eq!(solution.solute_moles(), original.solute_moles());
    assert_eq!(
        solution.osmotic_pressure_pa().unwrap(),
        original.osmotic_pressure_pa().unwrap()
    );
    assert_eq!(
        solution.water_activity().unwrap(),
        original.water_activity().unwrap()
    );
    assert_eq!(solution.heap_bytes(), 2 * core::mem::size_of::<f64>());
    assert!(solution.nonideal_model().is_none());
}

#[test]
fn g1_nonideal_dilution_retains_moles_and_refreshes_actual_activity() {
    let source = model();
    let solution = AqueousSolution::nonideal(source.clone(), 1e-3, vec![0.1]).unwrap();
    assert!(!solution.is_ideal());
    let expected = source.evaluate(&[100.0]).unwrap();
    assert_eq!(
        solution.osmotic_pressure_pa().unwrap(),
        expected.osmotic_pressure_pa
    );
    assert_eq!(solution.water_activity().unwrap(), expected.water_activity);
    let diluted = solution.with_solvent_volume_m3(2e-3).unwrap();
    assert_eq!(diluted.solute_moles(), &[0.1]);
    assert!(diluted.water_activity().unwrap() > solution.water_activity().unwrap());
    assert_eq!(
        diluted.osmotic_pressure_pa().unwrap(),
        source.evaluate(&[50.0]).unwrap().osmotic_pressure_pa
    );
    assert!(Arc::ptr_eq(
        solution.model.as_ref().unwrap(),
        diluted.model.as_ref().unwrap()
    ));
    assert_eq!(solution.temperature_k(), source.temperature_k());
    assert_eq!(
        solution.solvent_molar_volume_m3_mol(),
        source.solvent_molar_volume_m3_mol()
    );
}

#[test]
fn g0_domain_refusal_does_not_mutate_finite_inventory() {
    let mut solution = AqueousSolution::nonideal(model(), 1e-3, vec![0.1]).unwrap();
    let before = solution.clone();
    assert!(solution.replace_inventories(1e-4, &[0.1]).is_err());
    assert_eq!(solution, before);
    assert!(solution.replace_inventories(2e-3, &[-0.1]).is_err());
    assert_eq!(solution, before);
    assert!(solution.replace_inventories(2e-3, &[0.1, 0.2]).is_err());
    assert_eq!(solution, before);
    solution.replace_inventories(2e-3, &[0.05]).unwrap();
    assert_eq!(solution.solute_moles(), &[0.05]);
    let mut concentrations = [0.0];
    solution.concentrations_into(&mut concentrations).unwrap();
    assert_eq!(concentrations, [25.0]);
}

#[test]
fn g0_trials_are_read_only_and_count_shared_model_storage() {
    let source = model();
    let solution = AqueousSolution::nonideal(source.clone(), 1e-3, vec![0.1]).unwrap();
    let before = solution.clone();
    assert_eq!(
        solution.osmotic_pressure_at(2e-3, &[0.1]).unwrap(),
        source.evaluate(&[50.0]).unwrap().osmotic_pressure_pa
    );
    assert!(solution.osmotic_pressure_at(1e-5, &[0.1]).is_err());
    assert_eq!(solution, before);
    assert!(
        solution.heap_bytes()
            >= core::mem::size_of::<f64>()
                + core::mem::size_of::<VirialSolutionModel>()
                + source.heap_bytes()
    );
    assert_eq!(solution.clone().heap_bytes(), solution.heap_bytes());
}

#[test]
fn g0_nonideal_construction_refuses_unsupported_or_malformed_inventory() {
    assert!(AqueousSolution::nonideal(model(), 1e-5, vec![0.1]).is_err());
    assert!(AqueousSolution::nonideal(model(), 0.0, vec![0.1]).is_err());
    assert!(AqueousSolution::nonideal(model(), 1e-3, vec![]).is_err());
    let pure = AqueousSolution::nonideal(model(), 1e-3, vec![0.0]).unwrap();
    assert_eq!(pure.water_activity().unwrap(), 1.0);
    assert_eq!(pure.osmotic_pressure_pa().unwrap(), 0.0);
}
