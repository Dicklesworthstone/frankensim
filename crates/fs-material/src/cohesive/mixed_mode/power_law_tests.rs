use super::*;
use crate::cohesive::CohesiveViscosity;

fn coefficients() -> MixedModeCoefficients {
    MixedModeCoefficients {
        normal_stiffness_pa_m: 1e9,
        tangential_stiffness_pa_m: 4e8,
        normal_strength_pa: 2e5,
        tangential_strength_pa: 1.6e5,
        mode_i_fracture_energy_j_m2: 100.0,
        mode_ii_fracture_energy_j_m2: 400.0,
    }
}
fn law(alpha: f64) -> BilinearCohesive {
    BilinearCohesive::new_power_law_mixed_mode(coefficients(), alpha, 2e9, 0.02).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn expected(c: MixedModeCoefficients, g: [f64; 3], alpha: f64) -> f64 {
    let opening = c.normal_stiffness_pa_m * g[0].max(0.0).powi(2);
    let shear = c.tangential_stiffness_pa_m * (g[1] * g[1] + g[2] * g[2]);
    ((opening / (opening + shear) / c.mode_i_fracture_energy_j_m2).powf(alpha)
        + (shear / (opening + shear) / c.mode_ii_fracture_energy_j_m2).powf(alpha))
    .powf(-1.0 / alpha)
}

#[test]
fn pure_and_mixed_mode_work_follow_independent_power_law_envelopes() {
    for alpha in [0.5, 0.75, 1.0, 2.0, 4.0] {
        let law = law(alpha);
        assert!(law.mixed_mode_parameters().is_none());
        assert_eq!(law.mixed_mode_coefficients(), Some(coefficients()));
        assert_eq!(
            law.mixed_mode_criterion(),
            Some(MixedModeFractureCriterion::PowerLaw { exponent: alpha })
        );
        for direction in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.6, -0.8],
            [0.6, 0.64, -0.48],
        ] {
            let target = expected(coefficients(), direction, alpha);
            let end = direction.map(|g| g * 0.012);
            let one = law.bonded_point().trial(end, 1e-8).unwrap();
            assert!(one.point.failed());
            close(one.energy.work, target, 3e-12);
            let mut point = law.bonded_point();
            let mut work = 0.0;
            for i in 1..=40 {
                let step = point
                    .trial(end.map(|g| g * f64::from(i) / 40.0), 1e-8)
                    .unwrap();
                assert!(step.point.damage() >= point.damage());
                assert!(step.energy.dissipated >= 0.0);
                work += step.energy.work;
                point = step.point;
            }
            close(work, target, 3e-12);
            close(point.dissipated_j_m2(), target, 3e-12);
        }
    }
}

#[test]
fn weak_opening_survives_a_shear_fraction_rounded_to_one() {
    let law = law(0.5);
    let gap = [1e-12, 0.001, 0.0];
    let history = law
        .bonded_point()
        .trial(gap, 1e-8)
        .unwrap()
        .point
        .mixed_mode_history()
        .unwrap();
    assert_eq!(history.energy_shear_fraction, 1.0);
    assert!(history.fracture_energy_j_m2 < 400.0 - 1e-6);
    close(
        history.fracture_energy_j_m2,
        expected(coefficients(), gap, 0.5),
        2e-14,
    );
}

#[test]
fn interior_softening_and_complete_failure_bounds_are_admitted_for_all_modes() {
    let mut c = coefficients();
    c.normal_stiffness_pa_m = 2.0;
    c.tangential_stiffness_pa_m = 2.0;
    c.normal_strength_pa = 3.0_f64.sqrt();
    c.tangential_strength_pa = 3.0_f64.sqrt();
    c.mode_i_fracture_energy_j_m2 = 1.0;
    c.mode_ii_fracture_energy_j_m2 = 1.0;
    // Each pure onset consumes 0.75 Gc, but alpha=1/2 requires their sum < 1.
    assert!(BilinearCohesive::new_power_law_mixed_mode(c, 0.5, 2.0, 10.0).is_err());
    assert!(BilinearCohesive::new_power_law_mixed_mode(c, 1.0, 2.0, 10.0).is_ok());
    c.normal_strength_pa = 1.0;
    c.tangential_strength_pa = 1.0;
    // Pure failure is at 2, but equal-mode alpha=2 failure is at 2 sqrt(2).
    assert!(BilinearCohesive::new_power_law_mixed_mode(c, 2.0, 2.0, 2.5).is_err());
    assert!(BilinearCohesive::new_power_law_mixed_mode(c, 2.0, 2.0, 3.0).is_ok());
    for alpha in [f64::NAN, f64::INFINITY, 0.499, 0.0, -1.0] {
        assert!(BilinearCohesive::new_power_law_mixed_mode(c, alpha, 2.0, 10.0).is_err());
    }
}

#[test]
fn first_initiation_endpoint_and_path_tangents_include_envelope_sensitivity() {
    for alpha in [0.5, 0.75, 1.0, 2.0] {
        for law in [
            law(alpha),
            law(alpha).with_viscosity(CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap()),
        ] {
            let old = law
                .bonded_point()
                .trial_over_time([6e-5, 3e-5, -2e-5], 0.2, 1e-8)
                .unwrap()
                .equilibrium
                .point;
            for end in [[4e-4, 6e-4, -3e-4], [0.003, 0.008, -0.004]] {
                let (step, _, tangent) = old
                    .trial_over_time_with_mean_traction_tangent(end, 0.2, 1e-8)
                    .unwrap();
                let endpoint = step.equilibrium.tangent_pa_m.unwrap();
                let tangent = tangent.unwrap();
                for j in 0..3 {
                    let h = 2e-9;
                    let mut plus = end;
                    plus[j] += h;
                    let mut minus = end;
                    minus[j] -= h;
                    let (a, am, _) = old
                        .trial_over_time_with_mean_traction_tangent(plus, 0.2, 1e-8)
                        .unwrap();
                    let (b, bm, _) = old
                        .trial_over_time_with_mean_traction_tangent(minus, 0.2, 1e-8)
                        .unwrap();
                    for i in 0..3 {
                        let fd = (am[i] - bm[i]) / (2.0 * h);
                        assert!(
                            (fd - tangent[i][j]).abs() < 300.0,
                            "alpha={alpha} path[{i},{j}]: {fd} vs {}",
                            tangent[i][j]
                        );
                        let fd = (a.equilibrium.traction_pa[i] - b.equilibrium.traction_pa[i])
                            / (2.0 * h);
                        assert!(
                            (fd - endpoint[i][j]).abs() < 300.0,
                            "alpha={alpha} endpoint[{i},{j}]: {fd} vs {}",
                            endpoint[i][j]
                        );
                    }
                }
                assert!(old.mixed_mode_history().is_none());
            }
        }
    }
}

#[test]
fn half_power_pure_mode_cusp_declines_only_new_classical_derivatives() {
    let law = law(0.5);
    let (first, _, tangent) = law
        .bonded_point()
        .trial_with_mean_traction_tangent([0.0003, 0.0, 0.0], 1e-8)
        .unwrap();
    assert!(first.tangent_pa_m.is_none());
    assert!(tangent.is_none());
    assert!(
        first
            .branch_tangent_pa_m
            .iter()
            .flatten()
            .all(|v| v.is_finite())
    );
    let (next, _, tangent) = first
        .point
        .trial_with_mean_traction_tangent([0.0004, 0.0, 0.0], 1e-8)
        .unwrap();
    assert!(next.tangent_pa_m.is_some());
    assert!(tangent.is_some());
    assert_eq!(
        next.point.mixed_mode_history(),
        first.point.mixed_mode_history()
    );
    let smooth =
        BilinearCohesive::new_power_law_mixed_mode(coefficients(), 0.75, 2e9, 0.02).unwrap();
    assert!(
        smooth
            .bonded_point()
            .trial([0.0003, 0.0, 0.0], 1e-8)
            .unwrap()
            .tangent_pa_m
            .is_some()
    );
}

#[test]
fn thermal_transfer_retains_physical_onset_and_fracture_family_without_healing() {
    let cold = law(2.0);
    let first = cold
        .bonded_point()
        .trial([0.0004, 0.0003, 0.0002], 1e-8)
        .unwrap()
        .point;
    let point = first.trial([-0.0001, 0.0001, -0.0001], 1e-8).unwrap().point;
    let original = point.mixed_mode_history().unwrap();
    let mut c = coefficients();
    c.normal_stiffness_pa_m *= 0.8;
    c.tangential_stiffness_pa_m *= 1.2;
    let hot = BilinearCohesive::new_power_law_mixed_mode(c, 0.75, 2e9, 0.02).unwrap();
    let transfer = point.transfer_temperature(&hot, 1e-8).unwrap();
    let h = transfer.point.mixed_mode_history().unwrap();
    assert_eq!(h.onset_direction, original.onset_direction);
    assert_eq!(h.energy_shear_fraction, original.energy_shear_fraction);
    close(
        h.fracture_energy_j_m2,
        expected(c, original.onset_direction, 0.75),
        3e-12,
    );
    assert!(transfer.point.damage() >= point.damage());
    close(
        transfer.parameter_work_j_m2,
        transfer.stored_change_j_m2 + transfer.damage_dissipated_j_m2,
        1e-12,
    );
    let bk = BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: c.normal_stiffness_pa_m,
            tangential_stiffness_pa_m: c.tangential_stiffness_pa_m,
            normal_strength_pa: c.normal_strength_pa,
            tangential_strength_pa: c.tangential_strength_pa,
            mode_i_fracture_energy_j_m2: c.mode_i_fracture_energy_j_m2,
            mode_ii_fracture_energy_j_m2: c.mode_ii_fracture_energy_j_m2,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap();
    assert!(point.transfer_temperature(&bk, 1e-8).is_err());
    let common = BilinearCohesive::new(1e9, 2e5, 100.0, 2e9, 0.02).unwrap();
    assert!(point.transfer_temperature(&common, 1e-8).is_err());
}

#[test]
fn damaged_contact_pressure_gradient_tracks_power_law_initiation_and_reclosure() {
    for alpha in [0.5, 0.75, 2.0] {
        let law = law(alpha);
        let old = law
            .bonded_point()
            .trial([0.00018, 0.0001, 0.0], 1e-8)
            .unwrap()
            .point;
        let end = [-0.00015, 0.001, -0.0003];
        let (step, _, pressure, _, gradient) = old
            .trial_over_time_with_damage_weighted_compression_tangent(end, 0.2, 1e-8)
            .unwrap();
        assert!(pressure > 0.0);
        assert!(
            step.equilibrium
                .point
                .mixed_mode_history()
                .unwrap()
                .onset_direction[0]
                > 0.0
        );
        let gradient = gradient.unwrap();
        for j in 0..3 {
            let h = 2e-9;
            let mut plus = end;
            plus[j] += h;
            let mut minus = end;
            minus[j] -= h;
            let a = old
                .trial_over_time_with_damage_weighted_compression(plus, 0.2, 1e-8)
                .unwrap()
                .2;
            let b = old
                .trial_over_time_with_damage_weighted_compression(minus, 0.2, 1e-8)
                .unwrap()
                .2;
            assert!(((a - b) / (2.0 * h) - gradient[j]).abs() < 300.0);
        }
        assert!(old.mixed_mode_history().is_none());
    }
}
