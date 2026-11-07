//! Independent surface-energy increments under retained, dominant preloads.

use super::mixed_mode::MixedModeParameters;
use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn mixed(compression: f64) -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 10.0,
            tangential_stiffness_pa_m: 40.0,
            normal_strength_pa: 5.0,
            tangential_strength_pa: 10.0,
            mode_i_fracture_energy_j_m2: 5.0,
            mode_ii_fracture_energy_j_m2: 20.0,
            bk_exponent: 1.6,
        },
        compression,
        10.0,
    )
    .unwrap()
}

#[test]
fn g1_compressive_preload_does_not_erase_elastic_shear_fracture_or_unloading() {
    let common = BilinearCohesive::new(10.0, 5.0, 5.0, 1e24, 10.0).unwrap();
    let mixed = mixed(1e24);
    for (law, elastic_gap, soft_gap, soft_storage, fracture, failure) in [
        (&common, 0.25, 1.0, 5.0 / 3.0, 5.0 / 3.0, 2.0),
        (&mixed, 0.125, 1.0, 4.0, 4.0, 4.0),
    ] {
        let compressed = law
            .bonded_point()
            .trial([-1.0, 0.0, 0.0], 1e-10)
            .unwrap()
            .point;
        let (elastic, mean) = compressed
            .trial_with_mean_traction([-1.0, elastic_gap, 0.0], 1e-10)
            .unwrap();
        // Both scalar totals round to the same dominant compression inventory.
        assert_eq!(
            elastic.point.stored_energy_j_m2(),
            compressed.stored_energy_j_m2()
        );
        near(elastic.energy.work, 0.3125, 1e-12);
        near(elastic.energy.stored_change, 0.3125, 1e-12);
        near(mean[1] * elastic_gap, 0.3125, 1e-12);
        assert_eq!(elastic.energy.dissipated, 0.0);
        let (soft, mean) = elastic
            .point
            .trial_with_mean_traction([-1.0, soft_gap, 0.0], 1e-10)
            .unwrap();
        near(soft.energy.stored_change, soft_storage - 0.3125, 1e-12);
        near(soft.energy.dissipated, fracture, 1e-12);
        near(soft.energy.work, soft_storage - 0.3125 + fracture, 1e-12);
        near(mean[1] * (soft_gap - elastic_gap), soft.energy.work, 1e-12);
        let unload = soft.point.trial([-1.0, 0.0, 0.0], 1e-10).unwrap();
        near(unload.energy.stored_change, -soft_storage, 1e-12);
        near(unload.energy.work, -soft_storage, 1e-12);
        assert_eq!(unload.energy.dissipated, 0.0);
        assert_eq!(unload.point.damage(), soft.point.damage());
        let failed = unload.point.trial([-1.0, failure, 0.0], 1e-10).unwrap();
        assert!(failed.point.failed());
        near(failed.energy.stored_change, 0.0, 1e-12);
        near(
            elastic.energy.work + soft.energy.work + unload.energy.work + failed.energy.work,
            fracture + failed.energy.dissipated,
            1e-11,
        );
        let detached = failed.point.trial([-1.0, 1e100, -1e100], 1e-10).unwrap();
        assert_eq!(detached.energy.work, 0.0);
        assert_eq!(detached.energy.stored_change, 0.0);
        assert_eq!(detached.energy.dissipated, 0.0);
    }
}

#[test]
fn g1_small_orthogonal_motion_retains_work_beside_unchanged_active_storage() {
    let common = BilinearCohesive::new(10.0, 100.0, 1000.0, 20.0, 30.0).unwrap();
    let mixed = BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 10.0,
            tangential_stiffness_pa_m: 40.0,
            normal_strength_pa: 100.0,
            tangential_strength_pa: 400.0,
            mode_i_fracture_energy_j_m2: 1000.0,
            mode_ii_fracture_energy_j_m2: 4000.0,
            bk_exponent: 1.6,
        },
        20.0,
        30.0,
    )
    .unwrap();
    for (law, shear_stiffness) in [(&common, 10.0), (&mixed, 40.0)] {
        let loaded = law
            .bonded_point()
            .trial([0.0, 1.0, 0.0], 1e-12)
            .unwrap()
            .point;
        let delta = 1e-9;
        let expected = 0.5 * shear_stiffness * delta * delta;
        let (step, mean) = loaded
            .trial_with_mean_traction([0.0, 1.0, delta], 1e-30)
            .unwrap();
        assert_eq!(step.point.stored_energy_j_m2(), loaded.stored_energy_j_m2());
        near(step.energy.work, expected, 1e-30);
        near(step.energy.stored_change, expected, 1e-30);
        near(mean[2] * delta, expected, 1e-30);
        let (back, _) = step
            .point
            .trial_with_mean_traction(loaded.jump_m(), 1e-30)
            .unwrap();
        near(back.energy.work, -expected, 1e-30);
        near(back.energy.stored_change, -expected, 1e-30);
        assert_eq!(back.point.damage(), 0.0);
    }
}

#[test]
fn g1_new_fracture_is_not_a_difference_of_accumulated_loss_inventories() {
    let law = BilinearCohesive::new(10.0, 5.0, 5.0, 20.0, 10.0).unwrap();
    let loaded = law
        .bonded_point()
        .trial([0.75, 0.0, 0.0], 1e-12)
        .unwrap()
        .point;
    // This independent history anchor is the same field used by thermal law
    // transfers; previously spent energy cannot suppress a later small loss.
    let point = CohesivePoint {
        dissipated_reference: 1e24,
        ..loaded
    };
    let step = point.trial([1.25, 0.0, 0.0], 1e-12).unwrap();
    assert_eq!(point.dissipated_j_m2(), step.point.dissipated_j_m2());
    near(step.energy.dissipated, 5.0 / 3.0, 1e-12);
    near(step.energy.work, 5.0 / 3.0, 1e-12);
    near(step.energy.stored_change, 0.0, 1e-12);

    let law = mixed(20.0);
    let loaded = law
        .bonded_point()
        .trial([0.0, 1.5, 0.0], 1e-11)
        .unwrap()
        .point;
    let point = CohesivePoint {
        dissipated_reference: 1e24,
        ..loaded
    };
    let step = point.trial([0.0, 2.5, 0.0], 1e-11).unwrap();
    assert_eq!(point.dissipated_j_m2(), step.point.dissipated_j_m2());
    near(step.energy.dissipated, 20.0 / 3.75, 1e-11);
    near(step.energy.stored_change, 0.0, 1e-11);
    near(step.energy.work, 20.0 / 3.75, 1e-11);
}

#[test]
fn g1_thermal_adhesive_release_stays_separate_from_held_compression() {
    let cold = BilinearCohesive::new(10.0, 5.0, 5.0, 1e24, 10.0).unwrap();
    let hot = BilinearCohesive::new(4.0, 1.0, 0.5, 1e24, 10.0).unwrap();
    let old = cold
        .bonded_point()
        .trial([-1.0, 0.5, 0.0], 1e-12)
        .unwrap()
        .point;
    let update = old.transfer_temperature(&hot, 1e-12).unwrap();
    assert_eq!(old.stored_energy_j_m2(), update.point.stored_energy_j_m2());
    near(update.parameter_work_j_m2, -0.75, 1e-12);
    near(update.damage_dissipated_j_m2, 1.0 / 3.0, 1e-12);
    near(update.stored_change_j_m2, -13.0 / 12.0, 1e-12);
    near(update.point.dissipated_j_m2(), 1.0 / 3.0, 1e-12);
    let continued = update.point.trial([-1.0, 1.0, 0.0], 1e-12).unwrap();
    near(continued.energy.dissipated, 1.0 / 3.0, 1e-12);
    assert!(continued.point.failed());
}

#[test]
fn g1_physical_viscosity_keeps_its_own_loss_under_compressive_preload() {
    let law = BilinearCohesive::new(10.0, 5.0, 5.0, 1e24, 10.0)
        .unwrap()
        .with_viscosity(CohesiveViscosity::new(2.0, 3.0, 10.0).unwrap());
    let compressed = law
        .bonded_point()
        .trial_over_time([-1.0, 0.0, 0.0], 1.0, 1e-12)
        .unwrap()
        .equilibrium
        .point;
    let (step, mean) = compressed
        .trial_over_time_with_mean_traction([-1.0, 0.25, 0.0], 0.5, 1e-12)
        .unwrap();
    near(step.equilibrium.energy.stored_change, 0.3125, 1e-12);
    near(step.viscous_dissipated_j_m2, 3.0 * 0.5 * 0.25, 1e-12);
    near(mean[1] * 0.25, 0.3125 + 0.375, 1e-12);
    assert_eq!(step.equilibrium.energy.dissipated, 0.0);
}

#[test]
fn g1_subnormal_surviving_fraction_retains_representable_thermal_work() {
    let cold = BilinearCohesive::new(1e300, 2e-8, 1e-8, 20.0, 1.0).unwrap();
    let hot = BilinearCohesive::new(2e300, 4e-8, 2e-8, 20.0, 1.0).unwrap();
    let opening = 1.0_f64.next_down().next_down();
    let old = cold
        .bonded_point()
        .trial([opening, 0.0, 0.0], 1e-20)
        .unwrap()
        .point;
    let intact = cold.intact(old.maximum);
    assert_eq!(intact.to_bits(), 1);
    assert!(!old.failed());
    let update = old.transfer_temperature(&hot, 1e-35).unwrap();
    let expected = (0.5 * (intact * 1e300) * opening) * opening;
    assert!(expected > 0.0);
    near(update.parameter_work_j_m2, expected, 1e-38);
    near(update.stored_change_j_m2, expected, 1e-38);
    assert_eq!(update.damage_dissipated_j_m2, 0.0);
    assert_eq!(hot.intact(update.point.maximum).to_bits(), intact.to_bits());
}
