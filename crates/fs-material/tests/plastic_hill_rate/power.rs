//! G0/G1/G3 nonlinear kinetics against independently posed scalar flow,
//! isotropic J2 correspondence, tensor derivatives and unit covariance.
//! These are authored model fixtures, not empirical material qualification.

use super::*;
use fs_material::plastic::hill::rate::PowerLawOverstressHill48;
use fs_material::plastic::rate::PowerLawOverstressJ2;

const S0: f64 = 7.0;

#[test]
fn prescribed_anisotropic_multiplier_has_the_power_law_endpoint_and_exact_slope() {
    let c = coefficients();
    for exponent in [1.25, 2.0, 4.0, 12.0] {
        for duration in [0.004, 0.04, 0.4] {
            let law =
                PowerLawOverstressHill48::new(native(c, Q), ETA, S0, exponent, duration).unwrap();
            for (axis, coefficient) in [(3, c.n), (4, c.l), (5, c.m)] {
                for alpha in [1e-7_f64, 0.002, 0.007] {
                    let overstress = S0 * (ETA * alpha / (duration * S0)).powf(1.0 / exponent);
                    let radius = Y + H * alpha + Q * -(-B * alpha).exp_m1() + overstress;
                    let plastic = alpha * (coefficient / 2.0).sqrt();
                    let stress = radius / (2.0 * coefficient).sqrt() + C * plastic;
                    let twice_mu = E / (1.0 + NU);
                    let mut strain = [0.0; 6];
                    strain[axis] = plastic + stress / twice_mu;
                    let step = law
                        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
                        .unwrap();
                    near(step.response.plastic_increment, alpha, 3e-14);
                    near(step.response.state.plastic_strain[axis], plastic, 3e-14);
                    near(step.response.stress_pa[axis], stress, 3e-10);
                    near(step.response.yield_residual_pa, overstress, 3e-10);
                    near(
                        law.endpoint_overstress_pa(alpha).unwrap(),
                        overstress,
                        3e-12,
                    );
                    let kinetic_slope = overstress / (exponent * alpha);
                    let endpoint_modulus =
                        C + (H + Q * B * (-B * alpha).exp() + kinetic_slope) / coefficient;
                    near(
                        step.response.tangent_pa[axis][axis],
                        twice_mu * endpoint_modulus / (twice_mu + endpoint_modulus),
                        3e-7,
                    );
                    near(
                        step.energy.physical_dissipation_j_m3,
                        (Y + overstress) * alpha,
                        3e-12,
                    );
                    assert_eq!(
                        step.energy.final_state,
                        law.local_law()
                            .stored_energy(&strain, &step.response.state)
                            .unwrap()
                    );
                    law.validate_history(&strain, &step.response.state).unwrap();
                }
            }
        }
    }
}

#[test]
fn variable_duration_shear_holds_and_reversals_match_independent_bisection() {
    let coefficient = coefficients().n;
    let twice_mu = E / (1.0 + NU);
    let r = (2.0 * coefficient).sqrt();
    let flow = (coefficient / 2.0).sqrt();
    for exponent in [1.4, 3.0, 8.0] {
        let mut old = native(coefficients(), Q).initial_state();
        let mut old_strain = [0.0; 6];
        for (shear, duration) in [
            (0.012, 0.01),
            (0.012, 0.03),
            (0.012, 0.1),
            (-0.010, 0.05),
            (-0.010, 0.2),
            (0.014, 0.015),
        ] {
            let law = PowerLawOverstressHill48::new(
                native(coefficients(), Q),
                ETA,
                S0,
                exponent,
                duration,
            )
            .unwrap();
            let relative_trial =
                twice_mu * (shear - old.plastic_strain[3]) - C * old.plastic_strain[3];
            let q_trial = r * relative_trial.abs();
            let old_radius = Y + H * old.alpha + Q * -(-B * old.alpha).exp_m1();
            // Direct physical alpha equation, without modal normalization,
            // an analytic tangent, or the production safeguarded Newton step.
            let (mut lower, mut upper) = (
                0.0,
                (q_trial - old_radius).max(0.0) / (H + (twice_mu + C) * coefficient),
            );
            for _ in 0..90 {
                let mid = 0.5 * lower + 0.5 * upper;
                let overstress = S0 * (ETA * mid / (duration * S0)).powf(1.0 / exponent);
                let residual = q_trial
                    - (twice_mu + C) * coefficient * mid
                    - Y
                    - H * (old.alpha + mid)
                    - Q * -(-B * (old.alpha + mid)).exp_m1()
                    - overstress;
                if residual > 0.0 {
                    lower = mid;
                } else {
                    upper = mid;
                }
            }
            let delta = 0.5 * lower + 0.5 * upper;
            let plastic = old.plastic_strain[3] + relative_trial.signum() * flow * delta;
            let stress = twice_mu * (shear - plastic);
            let mut strain = [0.0; 6];
            strain[3] = shear;
            let saved = old.clone();
            let step = law
                .advance_with_energy(&old_strain, &old, &strain, TOL)
                .unwrap();
            assert_eq!(old, saved);
            assert_eq!(
                step,
                law.advance_with_energy(&old_strain, &old, &strain, TOL)
                    .unwrap()
            );
            near(step.response.state.alpha, old.alpha + delta, 3e-14);
            near(step.response.state.plastic_strain[3], plastic, 3e-14);
            near(step.response.stress_pa[3], stress, 3e-10);
            near(
                step.energy.physical_dissipation_j_m3,
                (Y + S0 * (ETA * delta / (duration * S0)).powf(1.0 / exponent)) * delta,
                3e-12,
            );
            assert!(step.energy.physical_dissipation_j_m3 >= 0.0);
            assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
            near(trace(&step.response.state.plastic_strain), 0.0, 3e-17);
            if shear == old_strain[3] {
                assert!(step.response.state.alpha > old.alpha);
                near(step.energy.mechanical_work_j_m3, 0.0, 0.0);
                near(
                    -step.energy.stored_change_j_m3(),
                    step.energy.physical_dissipation_j_m3 + step.energy.numerical_dissipation_j_m3,
                    3e-12,
                );
            }
            old = step.response.state;
            old_strain = strain;
        }
    }
}

#[test]
fn rotated_nonproportional_power_flow_retains_the_full_tensor_tangent() {
    let rotation = [[0.8, -0.48, 0.36], [0.6, 0.64, -0.48], [0.0, 0.6, 0.8]];
    let weights = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
    for exponent in [1.3, 4.0, 12.0] {
        let backbone = native(coefficients(), Q);
        let local = PowerLawOverstressHill48::new(backbone.clone(), ETA, S0, exponent, DT).unwrap();
        let spatial = PowerLawOverstressHill48::new_oriented(
            OrientedHill48::new(backbone, rotation).unwrap(),
            ETA,
            S0,
            exponent,
            DT,
        )
        .unwrap();
        let old = [0.008, -0.003, -0.005, 0.002, 0.001, 0.0];
        let state = local.response(&old, &local.initial_state()).unwrap().state;
        let strain = [-0.006, 0.004, 0.002, -0.003, 0.002, 0.001];
        let expected = local
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let spatial_old = rotate(&old, &rotation);
        let spatial_strain = rotate(&strain, &rotation);
        let step = spatial
            .advance_with_energy(&spatial_old, &state, &spatial_strain, TOL)
            .unwrap();
        let stress = rotate(&expected.response.stress_pa, &rotation);
        for i in 0..6 {
            near(step.response.stress_pa[i], stress[i], 3e-10);
            near(
                step.response.state.plastic_strain[i],
                expected.response.state.plastic_strain[i],
                3e-14,
            );
        }
        near(
            step.energy.mechanical_work_j_m3,
            expected.energy.mechanical_work_j_m3,
            3e-12,
        );
        near(
            step.energy.physical_dissipation_j_m3,
            expected.energy.physical_dissipation_j_m3,
            3e-12,
        );
        let delta = 1e-7;
        for j in 0..6 {
            let (mut plus, mut minus) = (spatial_strain, spatial_strain);
            plus[j] += delta;
            minus[j] -= delta;
            let p = spatial.stress(&plus, &state);
            let m = spatial.stress(&minus, &state);
            for i in 0..6 {
                near(
                    step.response.tangent_pa[i][j],
                    (p[i] - m[i]) / (2.0 * delta),
                    5e-6,
                );
                near(
                    weights[i] * step.response.tangent_pa[i][j],
                    weights[j] * step.response.tangent_pa[j][i],
                    3e-9,
                );
            }
        }
    }
}

#[test]
fn linear_limit_is_exact_and_isotropic_power_flow_matches_the_native_j2_owner() {
    let backbone = native(coefficients(), Q);
    let power = PowerLawOverstressHill48::new(backbone.clone(), ETA, S0, 1.0, DT).unwrap();
    let linear = LinearOverstressHill48::new(backbone, ETA, DT).unwrap();
    let path = [
        [0.0001, -0.00005, -0.00005, 0.0002, 0.0, 0.0],
        [0.008, -0.003, -0.005, 0.002, 0.0, 0.001],
        [0.008, -0.003, -0.005, 0.002, 0.0, 0.001],
        [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
    ];
    let (mut old, mut state) = ([0.0; 6], power.initial_state());
    for strain in path {
        let p = power
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        assert_eq!(
            p,
            linear
                .advance_with_energy(&old, &state, &strain, TOL)
                .unwrap()
        );
        assert_eq!(
            power.stress(&strain, &state),
            linear.stress(&strain, &state)
        );
        assert_eq!(
            power.update_state(&strain, &state),
            linear.update_state(&strain, &state)
        );
        old = strain;
        state = p.response.state;
    }
    for exponent in [1.7, 4.0] {
        let backbone = native(Hill48Coefficients::isotropic(), 0.0);
        let h = PowerLawOverstressHill48::new(backbone.clone(), ETA, S0, exponent, DT).unwrap();
        let j = PowerLawOverstressJ2::new(
            J2Plasticity::new_mixed(backbone.elastic().clone(), Y, H, 1.5 * C).unwrap(),
            ETA,
            S0,
            exponent,
            DT,
        )
        .unwrap();
        let (mut hs, mut js, mut old) = (h.initial_state(), j.initial_state(), [0.0; 6]);
        for strain in path {
            let hp = h.advance_with_energy(&old, &hs, &strain, TOL).unwrap();
            let jp = j.advance_with_energy(&old, &js, &strain, TOL).unwrap();
            let jt = j.tangent(&strain, &js);
            near(hp.response.state.alpha, jp.state.alpha, 3e-14);
            for i in 0..6 {
                near(hp.response.stress_pa[i], jp.stress_pa[i], 3e-10);
                near(
                    hp.response.state.plastic_strain[i],
                    jp.state.plastic_strain[i],
                    3e-14,
                );
                for k in 0..6 {
                    near(hp.response.tangent_pa[i][k], jt[i][k], 3e-7);
                }
            }
            near(
                hp.energy.physical_dissipation_j_m3,
                jp.energy.physical_dissipation_j_m3,
                3e-12,
            );
            near(
                hp.energy.final_state.total_j_m3(),
                jp.energy.final_state.total_j_m3(),
                3e-12,
            );
            near(
                hp.energy.numerical_dissipation_j_m3,
                jp.energy.numerical_dissipation_j_m3,
                3e-12,
            );
            hs = hp.response.state;
            js = jp.state;
            old = strain;
        }
    }
}

#[test]
fn power_kinetics_preserves_stress_unit_covariance_and_refuses_invalid_parameters() {
    let strain = [0.008, -0.003, -0.005, 0.002, 0.0, 0.001];
    let reference =
        PowerLawOverstressHill48::new(native(coefficients(), Q), ETA, S0, 4.0, DT).unwrap();
    let state = reference.initial_state();
    let expected = reference.response(&strain, &state).unwrap();
    for scale in [1e-180, 1e180] {
        let backbone = Hill48::new(
            IsotropicElastic::new(E * scale, NU, 0.05).unwrap(),
            Y * scale,
            H * scale,
            coefficients(),
        )
        .unwrap()
        .with_linear_kinematic_hardening(C * scale)
        .unwrap()
        .with_voce_isotropic_hardening(Q * scale, B)
        .unwrap();
        let law =
            PowerLawOverstressHill48::new(backbone, ETA * scale, S0 * scale, 4.0, DT).unwrap();
        let actual = law.response(&strain, &state).unwrap();
        near(actual.state.alpha, expected.state.alpha, 3e-14);
        for i in 0..6 {
            near(actual.stress_pa[i] / scale, expected.stress_pa[i], 3e-10);
            for j in 0..6 {
                near(
                    actual.tangent_pa[i][j] / scale,
                    expected.tangent_pa[i][j],
                    3e-7,
                );
            }
        }
    }
    for (eta, s0, n, dt) in [
        (ETA, 0.0, 2.0, DT),
        (ETA, S0, 0.9, DT),
        (ETA, S0, f64::NAN, DT),
        (ETA, S0, f64::INFINITY, DT),
        (0.0, S0, 2.0, DT),
        (ETA, S0, 2.0, 0.0),
        (ETA, f64::INFINITY, 2.0, DT),
        (f64::MAX, S0, 2.0, f64::MIN_POSITIVE),
    ] {
        assert!(PowerLawOverstressHill48::new(native(coefficients(), Q), eta, s0, n, dt).is_err());
    }
    assert!(reference.endpoint_overstress_pa(-1.0).is_err());
    assert!(reference.endpoint_overstress_pa(f64::NAN).is_err());
    let mut invalid = state;
    invalid.plastic_strain = [0.001; 6];
    assert!(reference.validate_history(&strain, &invalid).is_err());
}
