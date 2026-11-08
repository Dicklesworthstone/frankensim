//! Synthetic anisotropic coupons, isotropic correspondence, consistent
//! derivatives and energetic cycles; no measured-alloy calibration is claimed.

use fs_material::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};
use fs_material::plastic::reduction::{PlaneStressReduction, ReductionConfig};
use fs_material::plastic::voce::VoceJ2;
use fs_material::tensor::{contract, trace};
use fs_material::{SmallStrainLaw, Voigt};

const E: f64 = 1200.0;
const NU: f64 = 0.2;
const Y: f64 = 3.0;
const H: f64 = 40.0;
const C: f64 = 100.0;
const Q: f64 = 5.0;
const B: f64 = 180.0;
const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};
const PATH: [Voigt; 6] = [
    [0.0001, -0.00005, -0.00005, 0.0002, 0.0, 0.0],
    [0.008, -0.003, -0.005, 0.002, 0.0, 0.001],
    [0.0078, -0.0029, -0.0049, 0.00195, 0.0, 0.000975],
    [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
    [0.01, -0.004, -0.006, 0.005, -0.002, 0.0],
    [0.0; 6],
];

fn coefficients() -> Hill48Coefficients {
    Hill48Coefficients {
        f: 0.25,
        g: 0.9,
        h: 0.45,
        l: 1.4,
        m: 2.1,
        n: 1.8,
    }
}

fn base(coefficients: Hill48Coefficients, kinematic: f64) -> Hill48 {
    Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        Y,
        H,
        coefficients,
    )
    .unwrap()
    .with_linear_kinematic_hardening(kinematic)
    .unwrap()
}

fn law() -> Hill48 {
    base(coefficients(), C)
        .with_voce_isotropic_hardening(Q, B)
        .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, allowance={tolerance:.3e}"
    );
}

#[test]
fn zero_saturation_keeps_native_cycles_and_kinematic_selection_order() {
    let original = base(coefficients(), C);
    let zero = original
        .clone()
        .with_voce_isotropic_hardening(0.0, B)
        .unwrap();
    let reordered = base(coefficients(), 0.0)
        .with_voce_isotropic_hardening(Q, B)
        .unwrap()
        .with_linear_kinematic_hardening(C)
        .unwrap();
    assert_eq!(reordered, law());
    let mut state = original.initial_state();
    let mut old = [0.0; 6];
    for strain in PATH {
        assert_eq!(
            zero.response(&strain, &state).unwrap(),
            original.response(&strain, &state).unwrap()
        );
        let expected = original
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        assert_eq!(
            zero.advance_with_energy(&old, &state, &strain, TOL)
                .unwrap(),
            expected
        );
        assert_eq!(zero.card(), original.card());
        state = expected.response.state;
        old = strain;
    }
}

#[test]
fn isotropic_modes_recover_voce_j2_stress_history_tangent_and_energy() {
    for kinematic in [0.0, C] {
        let hill = base(Hill48Coefficients::isotropic(), kinematic)
            .with_voce_isotropic_hardening(Q, B)
            .unwrap();
        let j2 = VoceJ2::new(
            J2Plasticity::new_mixed(hill.elastic().clone(), Y, H, 1.5 * kinematic).unwrap(),
            Q,
            B,
        )
        .unwrap();
        let mut hs = hill.initial_state();
        let mut js = j2.initial_state();
        let mut old = [0.0; 6];
        for strain in PATH {
            let h = hill.advance_with_energy(&old, &hs, &strain, TOL).unwrap();
            let j = j2.advance_with_energy(&old, &js, &strain, TOL).unwrap();
            let tangent = j2.tangent(&strain, &js);
            for i in 0..6 {
                near(h.response.stress_pa[i], j.stress_pa[i], 3e-12);
                near(h.response.state.plastic_strain[i], j.state.plastic_strain[i], 3e-15);
                for k in 0..6 {
                    near(h.response.tangent_pa[i][k], tangent[i][k], 3e-9);
                }
            }
            near(h.response.state.alpha, j.state.alpha, 3e-15);
            near(h.energy.final_state.total_j_m3(), j.energy.final_state.total_j_m3(), 3e-13);
            near(h.energy.physical_dissipation_j_m3, j.energy.physical_dissipation_j_m3, 3e-13);
            near(h.energy.numerical_dissipation_j_m3, j.energy.numerical_dissipation_j_m3, 3e-13);
            hs = h.response.state;
            js = j.state;
            old = strain;
        }
    }
}

#[test]
fn anisotropic_shear_matches_independent_saturation_curve_and_endpoint_slope() {
    let c = coefficients();
    let twice_mu = E / (1.0 + NU);
    for (q, b) in [(Q, B), (Q, 100_000.0), (2e12, 1e-12)] {
        let material = base(c, C).with_voce_isotropic_hardening(q, b).unwrap();
        for (axis, coefficient) in [(3, c.n), (4, c.l), (5, c.m)] {
            for alpha in [1e-8, 0.001, 0.007] {
                // Prescribe the plastic multiplier, then construct total strain
                // from the original Hill surface and associated shear flow.
                let saturation = q * -(-b * alpha).exp_m1();
                let radius = Y + H * alpha + saturation;
                let plastic = alpha * (coefficient / 2.0).sqrt();
                let stress = radius / (2.0 * coefficient).sqrt() + C * plastic;
                let mut strain = [0.0; 6];
                strain[axis] = plastic + stress / twice_mu;
                let step = material
                    .advance_with_energy(&[0.0; 6], &material.initial_state(), &strain, TOL)
                    .unwrap();
                near(step.response.state.alpha, alpha, 3e-16);
                near(step.response.state.plastic_strain[axis], plastic, 3e-16);
                near(step.response.stress_pa[axis], stress, 3e-12);
                let slope = H + q * b * (-b * alpha).exp();
                near(material.isotropic_hardening_pa(alpha), H * alpha + saturation, 1e-13);
                near(material.isotropic_hardening_derivative_pa(alpha), slope, 1e-9);
                near(
                    step.response.tangent_pa[axis][axis],
                    twice_mu * (slope + C * coefficient)
                        / (slope + (twice_mu + C) * coefficient),
                    3e-9,
                );
                if b < 1e-10 {
                    // b*alpha <= 7e-15: the leading Voce integral is
                    // Q*b*alpha^2/2, even though 1-exp(-b*alpha) would cancel.
                    near(step.energy.final_state.hardening_j_m3,
                        0.5 * (H + q * b) * alpha * alpha, 3e-17);
                }
            }
        }
    }
}

#[test]
fn anisotropic_reversal_tangent_differentiates_the_same_frozen_history() {
    for kinematic in [0.0, C] {
        let material = base(coefficients(), kinematic)
            .with_voce_isotropic_hardening(Q, B)
            .unwrap();
        let old = material.response(&PATH[1], &material.initial_state()).unwrap().state;
        let strain = PATH[3];
        let result = material.response(&strain, &old).unwrap();
        assert!(result.plastic_increment > 0.0);
        let h = 1e-7;
        for j in 0..6 {
            let mut plus = strain;
            let mut minus = strain;
            plus[j] += h;
            minus[j] -= h;
            let p = material.response(&plus, &old).unwrap();
            let m = material.response(&minus, &old).unwrap();
            for i in 0..6 {
                near(result.tangent_pa[i][j], (p.stress_pa[i] - m.stress_pa[i]) / (2.0 * h), 3e-6);
                let wi = if i < 3 { 1.0 } else { 2.0 };
                let wj = if j < 3 { 1.0 } else { 2.0 };
                near(wi * result.tangent_pa[i][j], wj * result.tangent_pa[j][i], 3e-12);
            }
        }
    }
}

#[test]
fn cyclic_work_separates_saturation_storage_from_physical_heat() {
    let material = law();
    let mut state = material.initial_state();
    let mut old = [0.0; 6];
    let storage = |alpha: f64| Q * (alpha + (-B * alpha).exp_m1() / B);
    for strain in PATH {
        let step = material.advance_with_energy(&old, &state, &strain, TOL).unwrap();
        let next = &step.response.state;
        let delta = next.alpha - state.alpha;
        let dp: Voigt = core::array::from_fn(|i| next.plastic_strain[i] - state.plastic_strain[i]);
        let de = core::array::from_fn(|i| strain[i] - old[i] - dp[i]);
        let saturation_loss = Q * -(-B * next.alpha).exp_m1() * delta
            - (storage(next.alpha) - storage(state.alpha));
        near(
            step.energy.final_state.hardening_j_m3,
            0.5 * H * next.alpha.powi(2) + storage(next.alpha),
            3e-14,
        );
        near(
            step.energy.final_state.kinematic_hardening_j_m3,
            0.5 * C * contract(&next.plastic_strain, &next.plastic_strain),
            3e-14,
        );
        near(step.energy.physical_dissipation_j_m3, Y * delta, 3e-14);
        near(step.energy.numerical_dissipation_j_m3,
            0.5 * contract(&de, &material.elastic().stress(&de, &()))
                + 0.5 * H * delta * delta + 0.5 * C * contract(&dp, &dp) + saturation_loss,
            3e-13);
        near(step.energy.residual_j_m3, 0.0, 3e-13);
        near(trace(&next.plastic_strain), 0.0, 3e-17);
        assert!(step.energy.physical_dissipation_j_m3 >= 0.0);
        assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
        let retained = material.advance_with_energy(&strain, next, &strain, TOL).unwrap();
        assert_eq!(retained.response.state, *next);
        assert_eq!(retained.energy.physical_dissipation_j_m3, 0.0);
        assert_eq!(retained.energy.numerical_dissipation_j_m3, 0.0);
        state = step.response.state;
        old = strain;
    }
}

#[test]
fn plane_stress_retains_the_full_hardening_law_and_condensed_tangent() {
    let material = law();
    let reduced = PlaneStressReduction::new(
        material.clone(),
        ReductionConfig { max_iterations: 32, stress_tolerance_pa: 1e-11, strain_limit: 0.05 },
    ).unwrap();
    let loaded = reduced.evaluate([0.009, -0.002, 0.003], &reduced.initial_state()).unwrap();
    let strain = [-0.008, 0.003, -0.004];
    let returned = reduced.evaluate(strain, &loaded.state).unwrap();
    assert!(returned.state.history.alpha > loaded.state.history.alpha);
    for component in [2, 4, 5] {
        near(returned.full_stress_pa[component], 0.0, 1e-11);
    }
    let history = &returned.state.history;
    let relative =
        core::array::from_fn(|i| returned.full_stress_pa[i] - C * history.plastic_strain[i]);
    near(material.equivalent_stress_pa(&relative).unwrap(),
        Y + H * history.alpha + Q * -(-B * history.alpha).exp_m1(), 3e-11);
    let h = 1e-7;
    for j in 0..3 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += h;
        minus[j] -= h;
        let p = reduced.evaluate(plus, &loaded.state).unwrap();
        let m = reduced.evaluate(minus, &loaded.state).unwrap();
        for i in 0..3 {
            near(returned.tangent_pa[i][j], (p.stress_pa[i] - m.stress_pa[i]) / (2.0 * h), 3e-6);
        }
    }
}

#[test]
fn saturation_return_is_covariant_under_stress_unit_changes() {
    let reference = law();
    let previous = reference.response(&PATH[1], &reference.initial_state()).unwrap().state;
    let expected = reference.response(&PATH[3], &previous).unwrap();
    for scale in [1e-180, 1e180] {
        let scaled = Hill48::new(
            IsotropicElastic::new(E * scale, NU, 0.05).unwrap(),
            Y * scale, H * scale, coefficients(),
        ).unwrap().with_voce_isotropic_hardening(Q * scale, B).unwrap()
            .with_linear_kinematic_hardening(C * scale).unwrap();
        let actual = scaled.response(&PATH[3], &previous).unwrap();
        near(actual.state.alpha, expected.state.alpha, 3e-15);
        for i in 0..6 {
            near(actual.stress_pa[i] / scale, expected.stress_pa[i], 3e-12);
            near(actual.state.plastic_strain[i], expected.state.plastic_strain[i], 3e-15);
            for j in 0..6 {
                near(actual.tangent_pa[i][j] / scale, expected.tangent_pa[i][j], 3e-9);
            }
        }
    }
}

#[test]
fn invalid_saturation_parameters_are_refused() {
    for (q, b) in [(-1.0, B), (f64::NAN, B), (f64::INFINITY, B),
        (Q, 0.0), (Q, -1.0), (Q, f64::NAN), (Q, f64::INFINITY),
        (f64::MAX, 2.0), (f64::MIN_POSITIVE, f64::MIN_POSITIVE)] {
        assert!(base(coefficients(), C).with_voce_isotropic_hardening(q, b).is_err());
    }
}
