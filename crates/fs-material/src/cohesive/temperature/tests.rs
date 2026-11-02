use super::*;

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a:e} != {b:e}");
}
fn law(k: f64, strength: f64, energy: f64) -> BilinearCohesive {
    BilinearCohesive::new(k, strength, energy, 20.0, 10.0).unwrap()
}

#[test]
fn heating_separates_parameter_work_from_damage_release_at_fixed_jump() {
    let cold = law(10.0, 5.0, 5.0);
    let hot = law(4.0, 1.0, 0.5);
    let loaded = cold
        .bonded_point()
        .trial([0.5, 0.0, 0.0], 1e-10)
        .unwrap()
        .point;
    let step = loaded.transfer_temperature(&hot, 1e-10).unwrap();
    // Before heating: psi=1.25. New K at d=0 gives psi=0.5.
    // Hot envelope at g=.5 has d=2/3, psi=1/6 and t=2/3.
    near(step.parameter_work_j_m2, -0.75);
    near(step.point.damage(), 2.0 / 3.0);
    near(step.damage_dissipated_j_m2, 1.0 / 3.0);
    near(step.point.traction_pa()[0], 2.0 / 3.0);
    near(step.point.dissipated_j_m2(), 1.0 / 3.0);
    near(step.defect_j_m2, 0.0);
    assert_eq!(step.point.jump_m(), loaded.jump_m());
    // Later separation dissipates only the new envelope's remaining work;
    // the already accumulated fracture loss is not rescaled to the new Gc.
    let failed = step.point.trial([1.0, 0.0, 0.0], 1e-10).unwrap();
    near(failed.energy.dissipated, 1.0 / 3.0);
    near(failed.point.dissipated_j_m2(), 2.0 / 3.0);
    near(failed.energy.defect, 0.0);
}

#[test]
fn cooling_and_toughness_changes_keep_damage_and_cumulative_fracture() {
    let a = law(10.0, 5.0, 5.0);
    let b = law(20.0, 20.0, 40.0);
    let damaged = a
        .bonded_point()
        .trial([1.0, 0.0, 0.0], 1e-10)
        .unwrap()
        .point;
    let loss = damaged.dissipated_j_m2();
    let cooled = damaged.transfer_temperature(&b, 1e-10).unwrap();
    assert!(cooled.point.damage() >= damaged.damage());
    near(cooled.point.damage(), damaged.damage());
    near(cooled.point.dissipated_j_m2(), loss);
    near(cooled.damage_dissipated_j_m2, 0.0);
    near(cooled.parameter_work_j_m2, damaged.stored_energy_j_m2());
    let unloaded = cooled.point.trial([0.2, -0.1, 0.1], 1e-10).unwrap();
    near(unloaded.energy.dissipated, 0.0);
    let returned = unloaded.point.transfer_temperature(&a, 1e-10).unwrap();
    assert!(returned.point.damage() >= damaged.damage());
    near(returned.point.dissipated_j_m2(), loss);
    let (continued, _) = returned
        .point
        .trial_with_mean_traction([2.0, 0.0, 0.0], 1e-10)
        .unwrap();
    near(continued.point.dissipated_j_m2(), a.fracture_energy);
    assert!(continued.point.failed());
}

#[test]
fn failed_bond_keeps_compression_and_never_rebonds() {
    let a = law(10.0, 5.0, 5.0);
    let b = BilinearCohesive::new(20.0, 20.0, 40.0, 40.0, 10.0).unwrap();
    let failed = a
        .bonded_point()
        .trial([2.0, 0.0, 0.0], 1e-10)
        .unwrap()
        .point;
    let closed = failed.trial([-0.5, 0.0, 0.0], 1e-10).unwrap().point;
    let refreshed = closed.transfer_temperature(&b, 1e-10).unwrap();
    assert!(refreshed.point.failed());
    near(refreshed.point.dissipated_j_m2(), 5.0);
    near(refreshed.parameter_work_j_m2, 2.5);
    near(refreshed.damage_dissipated_j_m2, 0.0);
    near(refreshed.point.traction_pa()[0], -20.0);
    let same = refreshed.point.transfer_temperature(&b, 1e-10).unwrap();
    assert_eq!(same.point, refreshed.point);
    assert_eq!(same.parameter_work_j_m2, 0.0);
}

#[test]
fn temperature_transfer_refuses_mixed_history_and_out_of_domain_without_mutation() {
    use crate::cohesive::mixed_mode::MixedModeParameters;
    let common = law(10.0, 5.0, 5.0);
    let mixed = BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 10.0,
            tangential_stiffness_pa_m: 10.0,
            normal_strength_pa: 5.0,
            tangential_strength_pa: 5.0,
            mode_i_fracture_energy_j_m2: 5.0,
            mode_ii_fracture_energy_j_m2: 10.0,
            bk_exponent: 1.5,
        },
        20.0,
        10.0,
    )
    .unwrap();
    let p = common
        .bonded_point()
        .trial([1.0, 0.0, 0.0], 1e-10)
        .unwrap()
        .point;
    let saved = p;
    assert!(p.transfer_temperature(&mixed, 1e-10).is_err());
    assert!(
        mixed
            .bonded_point()
            .transfer_temperature(&common, 1e-10)
            .is_err()
    );
    let narrow = BilinearCohesive::new(10.0, 1.0, 0.1, 20.0, 0.5).unwrap();
    assert_eq!(
        p.transfer_temperature(&narrow, 1e-10),
        Err(CohesiveError::SeparationLimit)
    );
    assert!(p.transfer_temperature(&common, f64::NAN).is_err());
    assert_eq!(p, saved);
}
