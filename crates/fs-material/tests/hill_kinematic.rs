//! G0/G1 cyclic anisotropic shear for the original Hill return with Prager storage.
//! Synthetic coefficients test the model, not a measured cyclic metal calibration.
use fs_material::IsotropicElastic;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};
use fs_material::tensor::{contract, rotate};

const E: f64 = 200e9;
const NU: f64 = 0.3;
const M: f64 = E / (1.0 + NU);
const Y: f64 = 200e6;
const H: f64 = 2e9;
const C: f64 = 15e9;
const W: [f64; 6] = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];

fn base(hardening: f64) -> Hill48 {
    Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        Y,
        hardening,
        Hill48Coefficients::from_yield_ratios([1.0, 1.2, 1.5, 0.8, 1.2, 1.5]).unwrap(),
    )
    .unwrap()
}
fn law() -> Hill48 {
    base(H).with_linear_kinematic_hardening(C).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()).max(1.0),
        "{actual} != {expected}"
    );
}
fn tolerances() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-3,
        yield_pa: 2e-4,
    }
}

#[test]
fn reverse_shear_matches_all_directional_bauschinger_thresholds_and_hardening() {
    let law = law();
    let coefficients = law.coefficients();
    for (axis, a) in [coefficients.n, coefficients.l, coefficients.m]
        .into_iter()
        .enumerate()
    {
        let index = axis + 3;
        let radius_factor = (2.0 * a).sqrt();
        let flow_factor = (a / 2.0).sqrt();
        let mut strain = [0.0; 6];
        strain[index] = 0.004;
        let alpha = (radius_factor * M * strain[index] - Y) / (H + (M + C) * a);
        let plastic = alpha * flow_factor;
        let loaded = law.response(&strain, &law.initial_state()).unwrap();
        close(loaded.state.alpha, alpha, 2e-12);
        close(loaded.state.plastic_strain[index], plastic, 2e-12);
        close(
            loaded.stress_pa[index],
            M * (strain[index] - plastic),
            2e-12,
        );
        close(
            loaded.tangent_pa[index][index],
            M * (H + C * a) / (H + (M + C) * a),
            2e-12,
        );
        let radius = Y + H * alpha;
        let backstress = C * plastic;
        let reverse_yield = backstress - radius / radius_factor;
        // Translation advances reverse yield compared with the same expanded
        // surface centered at zero; this is the anisotropic Bauschinger effect.
        assert!(reverse_yield > -radius / radius_factor);
        strain[index] = plastic + (reverse_yield + 0.01 * Y / radius_factor) / M;
        let unloading = law.response(&strain, &loaded.state).unwrap();
        assert_eq!(unloading.state, loaded.state);
        close(
            unloading.stress_pa[index],
            reverse_yield + 0.01 * Y / radius_factor,
            2e-12,
        );
        strain[index] = plastic + (backstress - 1.3 * radius / radius_factor) / M;
        let reversed = law.response(&strain, &loaded.state).unwrap();
        let increment = 0.3 * radius / (H + (M + C) * a);
        let next_plastic = plastic - increment * flow_factor;
        close(reversed.plastic_increment, increment, 2e-12);
        close(reversed.state.alpha, alpha + increment, 2e-12);
        close(reversed.state.plastic_strain[index], next_plastic, 2e-12);
        close(
            reversed.stress_pa[index],
            M * (strain[index] - next_plastic),
            2e-12,
        );
        close(
            reversed.tangent_pa[index][index],
            M * (H + C * a) / (H + (M + C) * a),
            2e-12,
        );
        law.validate_committed_history(&strain, &reversed.state, 2e-4)
            .unwrap();
        assert_eq!(
            law.response(&strain, &reversed.state)
                .unwrap()
                .plastic_increment,
            0.0
        );
    }
}

#[test]
fn pure_kinematic_flow_retains_its_nonzero_exact_shear_tangent() {
    let law = base(0.0).with_linear_kinematic_hardening(C).unwrap();
    for index in 3..6 {
        let mut strain = [0.0; 6];
        strain[index] = 0.004;
        let result = law.response(&strain, &law.initial_state()).unwrap();
        assert!(result.plastic_increment > 0.0);
        close(result.tangent_pa[index][index], M * C / (M + C), 3e-12);
        assert!(result.tangent_pa[index][index] > 0.0);
    }
}

#[test]
fn mixed_loading_has_the_full_consistent_work_symmetric_tangent() {
    let law = law();
    let old = law
        .response(
            &[0.003, -0.0015, -0.0015, 0.002, 0.0, 0.0],
            &law.initial_state(),
        )
        .unwrap()
        .state;
    let strain = [-0.002, 0.0013, 0.0007, -0.002, 0.001, 0.0012];
    let result = law.response(&strain, &old).unwrap();
    assert!(result.plastic_increment > 0.0);
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        let step = 1e-8;
        plus[j] += step;
        minus[j] -= step;
        let p = law.response(&plus, &old).unwrap();
        let m = law.response(&minus, &old).unwrap();
        for i in 0..6 {
            let derivative = (p.stress_pa[i] - m.stress_pa[i]) / (2.0 * step);
            assert!(
                (result.tangent_pa[i][j] - derivative).abs() < 2e-7 * M,
                "({i},{j}): {} vs {derivative}",
                result.tangent_pa[i][j]
            );
            close(
                W[i] * result.tangent_pa[i][j],
                W[j] * result.tangent_pa[j][i],
                5e-14,
            );
        }
    }
}

#[test]
fn cyclic_work_separates_kinematic_storage_physical_heat_and_numerical_loss() {
    let law = law();
    let mut old = law.initial_state();
    let mut prior = [0.0; 6];
    for strain in [
        [0.0, 0.0, 0.0, 0.003, 0.001, 0.0],
        [0.001, -0.0005, -0.0005, 0.004, 0.0, 0.001],
        [0.001, -0.0005, -0.0005, 0.0038, 0.0, 0.001],
        [-0.001, 0.0005, 0.0005, -0.004, 0.001, -0.001],
        [0.0, 0.0, 0.0, 0.003, -0.001, 0.001],
    ] {
        let increment = law
            .advance_with_energy(&prior, &old, &strain, tolerances())
            .unwrap();
        let response = &increment.response;
        let energy = &increment.energy;
        close(
            energy.final_state.kinematic_hardening_j_m3,
            0.5 * C
                * contract(
                    &response.state.plastic_strain,
                    &response.state.plastic_strain,
                ),
            2e-13,
        );
        close(
            energy.physical_dissipation_j_m3,
            Y * response.plastic_increment,
            2e-13,
        );
        assert!(energy.physical_dissipation_j_m3 >= 0.0);
        assert!(energy.numerical_dissipation_j_m3 >= 0.0);
        assert!(energy.residual_j_m3.abs() < 1e-3);
        assert!(response.state.alpha >= old.alpha);
        old = response.state.clone();
        prior = strain;
    }
}

#[test]
fn zero_kinematic_modulus_preserves_original_responses_and_energy_exactly() {
    let original = base(H);
    let zero = original
        .clone()
        .with_linear_kinematic_hardening(0.0)
        .unwrap();
    assert_eq!(original, zero);
    let mut old = original.initial_state();
    let mut prior = [0.0; 6];
    for strain in [
        [0.0, 0.0, 0.0, 0.003, 0.001, 0.0],
        [-0.001, 0.0005, 0.0005, -0.004, 0.001, -0.001],
    ] {
        assert_eq!(
            original.response(&strain, &old).unwrap(),
            zero.response(&strain, &old).unwrap()
        );
        let a = original
            .advance_with_energy(&prior, &old, &strain, tolerances())
            .unwrap();
        let b = zero
            .advance_with_energy(&prior, &old, &strain, tolerances())
            .unwrap();
        assert_eq!(a, b);
        old = a.response.state;
        prior = strain;
    }
}

#[test]
fn fixed_orientation_preserves_the_same_backstress_history_and_work() {
    let local = law();
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let rotation = [[s, -s, 0.0], [s, s, 0.0], [0.0, 0.0, 1.0]];
    let oriented = OrientedHill48::new(local.clone(), rotation).unwrap();
    let mut old = local.initial_state();
    let mut prior = [0.0; 6];
    for strain in [
        [0.003, -0.0015, -0.0015, 0.002, 0.0, 0.001],
        [-0.002, 0.001, 0.001, -0.003, 0.001, 0.0],
    ] {
        let direct = local
            .advance_with_energy(&prior, &old, &strain, tolerances())
            .unwrap();
        let spatial = oriented
            .advance_with_energy(
                &rotate(&prior, &rotation),
                &old,
                &rotate(&strain, &rotation),
                tolerances(),
            )
            .unwrap();
        let expected = rotate(&direct.response.stress_pa, &rotation);
        for (a, b) in spatial.response.stress_pa.into_iter().zip(expected) {
            close(a, b, 3e-12);
        }
        for (a, b) in spatial
            .response
            .state
            .plastic_strain
            .into_iter()
            .zip(direct.response.state.plastic_strain)
        {
            close(a, b, 3e-12);
        }
        close(
            spatial.response.state.alpha,
            direct.response.state.alpha,
            3e-12,
        );
        close(
            spatial.energy.mechanical_work_j_m3,
            direct.energy.mechanical_work_j_m3,
            3e-12,
        );
        close(
            spatial.energy.final_state.total_j_m3(),
            direct.energy.final_state.total_j_m3(),
            3e-12,
        );
        old = direct.response.state;
        prior = strain;
    }
}

#[test]
fn invalid_kinematic_selection_refuses_without_changing_the_original_model() {
    let original = base(H);
    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            original
                .clone()
                .with_linear_kinematic_hardening(invalid)
                .is_err()
        );
        assert_eq!(original.kinematic_modulus_pa(), 0.0);
    }
}
