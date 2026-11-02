use super::*;
use crate::cohesive::mixed_mode::MixedModeParameters;

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11, "{a:e} != {b:e}");
}

fn cold_parameters() -> MixedModeParameters {
    MixedModeParameters {
        normal_stiffness_pa_m: 4.0,
        tangential_stiffness_pa_m: 9.0,
        normal_strength_pa: 4.0,
        tangential_strength_pa: 6.0,
        mode_i_fracture_energy_j_m2: 10.0,
        mode_ii_fracture_energy_j_m2: 20.0,
        bk_exponent: 2.0,
    }
}

fn hot_parameters() -> MixedModeParameters {
    MixedModeParameters {
        normal_stiffness_pa_m: 1.0,
        tangential_stiffness_pa_m: 4.0,
        normal_strength_pa: 0.5,
        tangential_strength_pa: 1.0,
        mode_i_fracture_energy_j_m2: 0.5,
        mode_ii_fracture_energy_j_m2: 1.0,
        bk_exponent: 2.0,
    }
}

fn mixed(p: MixedModeParameters) -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(p, 20.0, 50.0).unwrap()
}

#[test]
fn thermal_strength_crossing_initiates_on_the_held_physical_mode() {
    let cold = mixed(cold_parameters());
    let hot = mixed(MixedModeParameters {
        normal_stiffness_pa_m: 2.0,
        tangential_stiffness_pa_m: 3.0,
        normal_strength_pa: 1.0,
        tangential_strength_pa: 1.5,
        mode_i_fracture_energy_j_m2: 1.0,
        mode_ii_fracture_energy_j_m2: 3.0,
        bk_exponent: 1.5,
    });
    let loaded = cold
        .bonded_point()
        .trial([0.4, 0.3, -0.2], 1e-10)
        .unwrap()
        .point;
    assert!(loaded.mixed_mode_history().is_none());
    let transfer = loaded.transfer_temperature(&hot, 1e-10).unwrap();
    let history = transfer.point.mixed_mode_history().unwrap();
    let beta: f64 = 0.39 / 0.71;
    let radius = 0.71_f64.sqrt();
    let onset = radius / 1.16_f64.sqrt();
    let toughness = 1.0 + 2.0 * beta.powf(1.5);
    let failure = 2.0 * toughness / onset;
    let intact = onset / radius * (failure - radius) / (failure - onset);
    near(history.energy_shear_fraction, beta);
    near(history.envelope_shear_fraction, beta);
    near(history.onset_sqrt_energy, onset);
    near(history.fracture_energy_j_m2, toughness);
    for i in 0..3 {
        near(
            history.onset_direction[i],
            loaded.jump_m()[i] / 0.29_f64.sqrt(),
        );
    }
    near(transfer.point.damage(), 1.0 - intact);
    near(transfer.parameter_work_j_m2, -0.55);
    near(transfer.damage_dissipated_j_m2, 0.355 * (1.0 - intact));
    near(transfer.stored_change_j_m2, 0.355 * intact - 0.905);
    near(
        transfer.point.dissipated_j_m2(),
        transfer.damage_dissipated_j_m2,
    );
    near(transfer.defect_j_m2, 0.0);
    let below = cold
        .bonded_point()
        .trial([0.02, -0.01, 0.0], 1e-10)
        .unwrap()
        .point;
    let unchanged = below.transfer_temperature(&hot, 1e-10).unwrap();
    assert!(unchanged.point.mixed_mode_history().is_none());
    assert_eq!(unchanged.damage_dissipated_j_m2, 0.0);
}

#[test]
fn anisotropic_temperature_cycles_keep_first_direction_and_original_beta() {
    let cold = mixed(cold_parameters());
    let changed = mixed(MixedModeParameters {
        normal_stiffness_pa_m: 16.0,
        tangential_stiffness_pa_m: 1.0,
        normal_strength_pa: 10.0,
        tangential_strength_pa: 2.0,
        bk_exponent: 1.5,
        ..cold_parameters()
    });
    let loaded = cold
        .bonded_point()
        .trial([1.2, 0.8, -0.6], 1e-10)
        .unwrap()
        .point;
    let first = loaded.mixed_mode_history().unwrap();
    let rotated = loaded.trial([0.15, -0.2, 0.1], 1e-10).unwrap().point;
    let intact = 1.0 - rotated.damage();
    let loss = rotated.dissipated_j_m2();
    let transfer = rotated.transfer_temperature(&changed, 1e-10).unwrap();
    let h = transfer.point.mixed_mode_history().unwrap();
    let beta: f64 = 1.0 / 24.04;
    assert_eq!(h.onset_direction, first.onset_direction);
    assert_eq!(h.energy_shear_fraction, first.energy_shear_fraction);
    near(h.envelope_shear_fraction, beta);
    near(h.onset_sqrt_energy, 24.04_f64.sqrt() / 3.9364_f64.sqrt());
    near(h.fracture_energy_j_m2, 10.0 + 10.0 * beta.powf(1.5));
    assert!((h.envelope_shear_fraction - 0.05 / 0.41).abs() > 0.05);
    near(transfer.point.damage(), rotated.damage());
    near(transfer.point.dissipated_j_m2(), loss);
    near(transfer.damage_dissipated_j_m2, 0.0);
    near(transfer.parameter_work_j_m2, -0.065 * intact);
    let cooled = transfer.point.transfer_temperature(&cold, 1e-10).unwrap();
    let back = cooled.point.mixed_mode_history().unwrap();
    assert_eq!(back.onset_direction, first.onset_direction);
    assert_eq!(back.energy_shear_fraction, first.energy_shear_fraction);
    near(back.envelope_shear_fraction, first.energy_shear_fraction);
    near(cooled.point.damage(), rotated.damage());
    near(cooled.point.dissipated_j_m2(), loss);
    near(
        cooled.parameter_work_j_m2 + transfer.parameter_work_j_m2,
        0.0,
    );
    let again = cooled.point.transfer_temperature(&cold, 1e-10).unwrap();
    assert_eq!(again.point, cooled.point);
    assert_eq!(again.damage_dissipated_j_m2, 0.0);
}

#[test]
fn mixed_heating_release_and_rotating_reload_have_independent_energy_balances() {
    let cold = mixed(cold_parameters());
    let hot = mixed(hot_parameters());
    let loaded = cold
        .bonded_point()
        .trial([1.2, 0.8, -0.6], 1e-10)
        .unwrap()
        .point;
    let original = loaded.mixed_mode_history().unwrap();
    let rotated = loaded.trial([0.6, -0.8, 0.0], 1e-10).unwrap().point;
    let old_intact = 1.0 - rotated.damage();
    let transfer = rotated.transfer_temperature(&hot, 1e-10).unwrap();
    let beta = 4.0 / 5.44;
    let toughness = 0.5 + 0.5 * beta * beta;
    let failure = 4.0 * toughness;
    let radius = 2.92_f64.sqrt();
    let intact = 0.5 / radius * (failure - radius) / (failure - 0.5);
    let h = transfer.point.mixed_mode_history().unwrap();
    assert_eq!(h.onset_direction, original.onset_direction);
    assert_eq!(h.energy_shear_fraction, original.energy_shear_fraction);
    near(h.envelope_shear_fraction, beta);
    near(transfer.point.damage(), 1.0 - intact);
    near(transfer.parameter_work_j_m2, -2.14 * old_intact);
    near(
        transfer.damage_dissipated_j_m2,
        1.46 * (old_intact - intact),
    );
    near(
        transfer.point.dissipated_j_m2(),
        rotated.dissipated_j_m2() + transfer.damage_dissipated_j_m2,
    );
    let end = [0.0, 2.0, 0.0];
    let (whole, mean) = transfer.point.trial_with_mean_traction(end, 1e-10).unwrap();
    assert!(whole.point.failed());
    near(
        whole.energy.dissipated,
        toughness * (failure - radius) / (failure - 0.5),
    );
    let delta = core::array::from_fn::<_, 3, _>(|i| end[i] - transfer.point.jump_m()[i]);
    near(whole.energy.work, (0..3).map(|i| mean[i] * delta[i]).sum());
    // Independent midpoint integration across shear reversal and final rupture.
    let mut quadrature = 0.0;
    for k in 0..8192 {
        let s = (k as f64 + 0.5) / 8192.0;
        let gap = core::array::from_fn(|i| transfer.point.jump_m()[i] + s * delta[i]);
        let t = transfer.point.trial(gap, 1e-10).unwrap().traction_pa;
        quadrature += (0..3).map(|i| t[i] * delta[i]).sum::<f64>() / 8192.0;
    }
    assert!((whole.energy.work - quadrature).abs() < 2e-7);
    let middle = core::array::from_fn(|i| transfer.point.jump_m()[i] + 0.5 * delta[i]);
    let left = transfer.point.trial(middle, 1e-10).unwrap();
    let right = left.point.trial(end, 1e-10).unwrap();
    near(left.energy.work + right.energy.work, whole.energy.work);
    near(right.point.dissipated_j_m2(), whole.point.dissipated_j_m2());
    near(whole.energy.defect, 0.0);
}

#[test]
fn mixed_temperature_transfer_reduces_to_pure_opening_and_pure_shear() {
    let cold = mixed(cold_parameters());
    let hot = mixed(hot_parameters());
    for (component, k0, t0, gc0, k1, t1, gc1) in [
        (0, 4.0, 4.0, 10.0, 1.0, 0.5, 0.5),
        (1, 9.0, 6.0, 20.0, 4.0, 1.0, 1.0),
    ] {
        let common_cold = BilinearCohesive::new(k0, t0, gc0, 20.0, 50.0).unwrap();
        let common_hot = BilinearCohesive::new(k1, t1, gc1, 20.0, 50.0).unwrap();
        let mut gap = [0.0; 3];
        gap[component] = 1.2;
        let a = cold.bonded_point().trial(gap, 1e-10).unwrap().point;
        let b = common_cold.bonded_point().trial(gap, 1e-10).unwrap().point;
        let a = a.transfer_temperature(&hot, 1e-10).unwrap();
        let b = b.transfer_temperature(&common_hot, 1e-10).unwrap();
        near(a.parameter_work_j_m2, b.parameter_work_j_m2);
        near(a.damage_dissipated_j_m2, b.damage_dissipated_j_m2);
        near(a.point.damage(), b.point.damage());
        near(
            a.point.traction_pa()[component],
            b.point.traction_pa()[component],
        );
        near(a.point.dissipated_j_m2(), b.point.dissipated_j_m2());
        gap[component] = 4.0;
        let a = a.point.trial(gap, 1e-10).unwrap();
        let b = b.point.trial(gap, 1e-10).unwrap();
        assert!(a.point.failed() && b.point.failed());
        near(a.energy.work, b.energy.work);
        near(a.point.dissipated_j_m2(), b.point.dissipated_j_m2());
    }
}

#[test]
fn failed_mixed_bond_keeps_compression_and_drops_all_bond_viscosity() {
    use crate::cohesive::CohesiveViscosity;
    let cold =
        mixed(cold_parameters()).with_viscosity(CohesiveViscosity::new(2.0, 3.0, 50.0).unwrap());
    let hot = BilinearCohesive::new_mixed_mode(hot_parameters(), 40.0, 50.0)
        .unwrap()
        .with_viscosity(CohesiveViscosity::new(4.0, 5.0, 50.0).unwrap());
    let failed = cold
        .bonded_point()
        .trial_over_time([10.0, 0.0, 0.0], 1.0, 1e-10)
        .unwrap()
        .equilibrium
        .point;
    let closed = failed
        .trial_over_time([-0.4, 0.2, -0.1], 1.0, 1e-10)
        .unwrap()
        .equilibrium
        .point;
    let transfer = closed.transfer_temperature(&hot, 1e-10).unwrap();
    assert!(transfer.point.failed());
    near(transfer.parameter_work_j_m2, 1.6);
    near(transfer.damage_dissipated_j_m2, 0.0);
    near(transfer.point.dissipated_j_m2(), 10.0);
    assert_eq!(transfer.point.traction_pa(), [-16.0, 0.0, -0.0]);
    let opened = transfer
        .point
        .trial_over_time([0.2, -0.4, 0.3], 1.0, 1e-10)
        .unwrap();
    assert!(opened.equilibrium.point.failed());
    assert_eq!(opened.viscous_dissipated_j_m2, 0.0);
    near(opened.equilibrium.point.dissipated_j_m2(), 10.0);
    near(
        opened.total_work_j_m2()
            - opened.equilibrium.energy.stored_change
            - opened.equilibrium.energy.dissipated
            - opened.viscous_dissipated_j_m2,
        0.0,
    );
}
