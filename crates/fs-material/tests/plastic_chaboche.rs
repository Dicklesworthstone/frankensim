//! Independent nonlinear backstress, tangent and thermodynamic checks.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2, ChabocheState};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::{PlaneStressReduction, ReductionConfig};
use fs_material::plastic::voce::VoceJ2;
use fs_material::tensor::{contract, trace};
use fs_material::{SmallStrainLaw, Voigt};

const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};

fn backbone(kinematic: f64) -> VoceJ2 {
    VoceJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(),
            3.0,
            40.0,
            kinematic,
        )
        .unwrap(),
        5.0,
        180.0,
    )
    .unwrap()
}
fn law() -> ChabocheJ2<2> {
    ChabocheJ2::new(
        backbone(60.0),
        [
            BackstressBranch::new(800.0, 150.0).unwrap(),
            BackstressBranch::new(250.0, 20.0).unwrap(),
        ],
    )
    .unwrap()
}
fn deviatoric(value: f64) -> Voigt {
    [value, -0.5 * value, -0.5 * value, 0.0, 0.0, 0.0]
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; allowance {tolerance}"
    );
}
fn path() -> [Voigt; 5] {
    [
        deviatoric(0.012),
        deviatoric(0.01),
        deviatoric(-0.014),
        [-0.006, 0.004, 0.002, -0.005, 0.002, 0.001],
        [0.01, -0.002, -0.008, 0.004, -0.002, -0.001],
    ]
}

#[test]
fn g1_prescribed_alpha_recovers_independent_backward_euler_branch_formula() {
    let material = law();
    let initial = material.initial_state();
    for alpha in [1e-7, 0.003, 0.018] {
        let a0 = 800.0 * alpha / (1.0 + 150.0 * alpha);
        let a1 = 250.0 * alpha / (1.0 + 20.0 * alpha);
        let radius = 3.0 + 40.0 * alpha + 5.0 * -(-180.0_f64 * alpha).exp_m1();
        let strain = deviatoric((radius + 1560.0 * alpha + a0 + a1) / 1500.0);
        let response = material.response(&strain, &initial).unwrap();
        near(response.state.plastic.alpha, alpha, 2e-16);
        near(response.state.plastic.plastic_strain[0], alpha, 2e-16);
        near(
            response.stress_pa[0],
            (2.0 / 3.0) * (radius + 60.0 * alpha + a0 + a1),
            2e-12,
        );
        for (index, expected) in [a0, a1].into_iter().enumerate() {
            near(
                response.state.backstress_pa[index][0],
                (2.0 / 3.0) * expected,
                2e-13,
            );
            near(
                response.state.backstress_pa[index][1],
                -expected / 3.0,
                2e-13,
            );
            near(trace(&response.state.backstress_pa[index]), 0.0, 2e-13);
        }
        let denominator = 1600.0
            + 900.0 * (-180.0_f64 * alpha).exp()
            + 800.0 / (1.0 + 150.0 * alpha).powi(2)
            + 250.0 / (1.0 + 20.0 * alpha).powi(2);
        near(
            response.tangent_pa[0][0]
                - 0.5 * response.tangent_pa[0][1]
                - 0.5 * response.tangent_pa[0][2],
            1000.0 * (1.0 - 1500.0 / denominator),
            2e-10,
        );
        near(
            material.yield_function_with_state(&response.stress_pa, &response.state),
            0.0,
            2e-12,
        );
        assert!(response.iterations <= 96);
    }
}

#[test]
fn g0_absent_branches_preserve_original_voce_response_and_energy_exactly() {
    let original = backbone(60.0);
    let empty = ChabocheJ2::<0>::new(original.clone(), []).unwrap();
    let disabled = ChabocheJ2::new(
        original.clone(),
        [BackstressBranch::new(0.0, 500.0).unwrap()],
    )
    .unwrap();
    let mut original_state = original.initial_state();
    let mut empty_state = empty.initial_state();
    let mut disabled_state = disabled.initial_state();
    let mut old = [0.0; 6];
    for strain in path() {
        let expected = original
            .advance_with_energy(&old, &original_state, &strain, TOL)
            .unwrap();
        for (stress, tangent) in [
            (
                empty.stress(&strain, &empty_state),
                empty.tangent(&strain, &empty_state),
            ),
            (
                disabled.stress(&strain, &disabled_state),
                disabled.tangent(&strain, &disabled_state),
            ),
        ] {
            assert_eq!(stress, original.stress(&strain, &original_state));
            assert_eq!(tangent, original.tangent(&strain, &original_state));
        }
        let a = empty
            .advance_with_energy(&old, &empty_state, &strain, TOL)
            .unwrap();
        let b = disabled
            .advance_with_energy(&old, &disabled_state, &strain, TOL)
            .unwrap();
        assert_eq!(a.energy, expected.energy);
        assert_eq!(b.energy, expected.energy);
        assert_eq!(a.state.plastic, expected.state);
        assert_eq!(b.state.plastic, expected.state);
        assert_eq!(a.recovery_dissipation_j_m3, 0.0);
        assert_eq!(b.recovery_dissipation_j_m3, 0.0);
        original_state = expected.state;
        empty_state = a.state;
        disabled_state = b.state;
        old = strain;
    }
}

#[test]
fn g1_zero_recovery_matches_independent_linear_prager_backbone_on_reversal() {
    let material = ChabocheJ2::new(
        backbone(60.0),
        [
            BackstressBranch::new(100.0, 0.0).unwrap(),
            BackstressBranch::new(200.0, 0.0).unwrap(),
        ],
    )
    .unwrap();
    let linear = backbone(360.0);
    let mut state = material.initial_state();
    let mut linear_state = linear.initial_state();
    let mut old = [0.0; 6];
    assert!(material.admissibility().tangent_symmetric);
    for strain in path() {
        let a = material
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let b = linear
            .advance_with_energy(&old, &linear_state, &strain, TOL)
            .unwrap();
        for i in 0..6 {
            near(a.stress_pa[i], b.stress_pa[i], 2e-12);
            near(
                a.state.plastic.plastic_strain[i],
                b.state.plastic_strain[i],
                2e-15,
            );
            near(
                a.state.backstress_pa[0][i],
                (2.0 / 3.0) * 100.0 * a.state.plastic.plastic_strain[i],
                2e-13,
            );
            near(
                a.state.backstress_pa[1][i],
                (2.0 / 3.0) * 200.0 * a.state.plastic.plastic_strain[i],
                2e-13,
            );
        }
        near(
            a.energy.final_state.total_j_m3(),
            b.energy.final_state.total_j_m3(),
            2e-13,
        );
        near(
            a.energy.numerical_dissipation_j_m3,
            b.energy.numerical_dissipation_j_m3,
            2e-13,
        );
        assert_eq!(a.recovery_dissipation_j_m3, 0.0);
        state = a.state;
        linear_state = b.state;
        old = strain;
    }
}

#[test]
fn g1_nonproportional_tangent_is_analytic_consistent_and_not_symmetrized() {
    let material = law();
    let mut state = material.initial_state();
    for strain in path() {
        let response = material.response(&strain, &state).unwrap();
        let h = 1e-7;
        for j in 0..6 {
            let mut plus = strain;
            let mut minus = strain;
            plus[j] += h;
            minus[j] -= h;
            let above = material.response(&plus, &state).unwrap();
            let below = material.response(&minus, &state).unwrap();
            for i in 0..6 {
                near(
                    response.tangent_pa[i][j],
                    (above.stress_pa[i] - below.stress_pa[i]) / (2.0 * h),
                    4e-6,
                );
            }
        }
        state = response.state;
    }
    let strain = [-0.008, 0.003, 0.005, -0.005, 0.002, 0.004];
    let tangent = material.tangent(&strain, &state);
    let mut asymmetry = 0.0_f64;
    for i in 0..6 {
        for j in 0..6 {
            let wi = if i < 3 { 1.0 } else { 2.0 };
            let wj = if j < 3 { 1.0 } else { 2.0 };
            asymmetry = asymmetry.max((wi * tangent[i][j] - wj * tangent[j][i]).abs());
        }
    }
    assert!(
        asymmetry > 0.1,
        "fixture must expose recovery-induced asymmetry: {asymmetry}"
    );
    assert!(!material.admissibility().tangent_symmetric);
}

#[test]
fn g1_recovery_is_physical_dissipation_while_storage_and_integration_loss_stay_separate() {
    let material = law();
    let mut state = material.initial_state();
    let mut old = [0.0; 6];
    let mut recovered = 0.0;
    for strain in path() {
        let next = material
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let alpha = next.state.plastic.alpha;
        let delta = alpha - state.plastic.alpha;
        let dp: Voigt = core::array::from_fn(|i| {
            next.state.plastic.plastic_strain[i] - state.plastic.plastic_strain[i]
        });
        let de: Voigt = core::array::from_fn(|i| strain[i] - old[i] - dp[i]);
        let d_stress = material.backbone().backbone().elastic.stress(&de, &());
        let sat_storage = |p: f64| 5.0 * (p + (-180.0 * p).exp_m1() / 180.0);
        let sat_loss = 5.0 * -(-180.0 * alpha).exp_m1() * delta
            - (sat_storage(alpha) - sat_storage(state.plastic.alpha));
        let mut branch_storage = 0.0;
        let mut recovery = 0.0;
        let mut branch_numerical = 0.0;
        for k in 0..2 {
            let c = material.branches()[k].modulus_pa();
            let gamma = material.branches()[k].recovery_per_strain();
            let a0 = state.backstress_pa[k];
            let a1 = next.state.backstress_pa[k];
            let da = core::array::from_fn(|i| a1[i] - a0[i]);
            for i in 0..6 {
                near(
                    a1[i],
                    (a0[i] + (2.0 / 3.0) * c * dp[i]) / (1.0 + gamma * delta),
                    2e-13,
                );
            }
            branch_storage += 0.75 / c * contract(&a1, &a1);
            recovery += 1.5 * gamma / c * delta * contract(&a1, &a1);
            branch_numerical += 0.75 / c * contract(&da, &da);
        }
        near(
            next.energy.final_state.hardening_j_m3,
            20.0 * alpha * alpha + sat_storage(alpha),
            2e-14,
        );
        near(
            next.energy.final_state.kinematic_hardening_j_m3,
            20.0 * contract(
                &next.state.plastic.plastic_strain,
                &next.state.plastic.plastic_strain,
            ) + branch_storage,
            2e-14,
        );
        near(next.recovery_dissipation_j_m3, recovery, 2e-14);
        near(
            next.energy.physical_dissipation_j_m3,
            3.0 * delta + recovery,
            2e-14,
        );
        near(
            next.energy.numerical_dissipation_j_m3,
            0.5 * contract(&de, &d_stress)
                + 20.0 * delta * delta
                + 20.0 * contract(&dp, &dp)
                + sat_loss
                + branch_numerical,
            2e-13,
        );
        near(next.energy.residual_j_m3, 0.0, 2e-13);
        assert!(next.energy.physical_dissipation_j_m3 >= 0.0);
        assert!(next.energy.numerical_dissipation_j_m3 >= 0.0);
        let repeated = material.response(&strain, &next.state).unwrap();
        for i in 0..6 {
            near(repeated.stress_pa[i], next.stress_pa[i], 2e-12);
        }
        recovered += recovery;
        state = next.state;
        old = strain;
    }
    assert!(recovered > 0.001);
}

#[test]
fn g1_branch_refinement_converges_to_continuous_monotonic_saturation() {
    let material = ChabocheJ2::new(
        backbone(60.0),
        [BackstressBranch::new(800.0, 150.0).unwrap()],
    )
    .unwrap();
    let target_alpha = 0.02;
    let exact = (800.0 / 150.0) * -(-150.0_f64 * target_alpha).exp_m1();
    let simulate = |steps: usize| {
        let mut state = material.initial_state();
        let delta = target_alpha / steps as f64;
        let mut expected_a = 0.0;
        for k in 1..=steps {
            expected_a = (expected_a + 800.0 * delta) / (1.0 + 150.0 * delta);
            let alpha = k as f64 * delta;
            let radius = 3.0 + 40.0 * alpha + 5.0 * -(-180.0 * alpha).exp_m1();
            let strain = deviatoric(alpha + (radius + 60.0 * alpha + expected_a) / 1500.0);
            state = material.response(&strain, &state).unwrap().state;
            near(1.5 * state.backstress_pa[0][0], expected_a, 3e-12);
        }
        (1.5 * state.backstress_pa[0][0] - exact).abs()
    };
    let coarse = simulate(20);
    let fine = simulate(40);
    assert!(fine < 0.55 * coarse && fine > 0.45 * coarse);
}

#[test]
fn g1_plane_stress_retains_recovery_history_and_nonsymmetric_condensed_tangent() {
    let material = law();
    let plane = PlaneStressReduction::new(
        material.clone(),
        ReductionConfig {
            max_iterations: 40,
            stress_tolerance_pa: 1e-11,
            strain_limit: 0.05,
        },
    )
    .unwrap();
    let loaded = plane
        .evaluate([0.012, -0.003, 0.003], &plane.initial_state())
        .unwrap();
    let strain = [-0.009, 0.003, -0.004];
    let next = plane.evaluate(strain, &loaded.state).unwrap();
    assert!(next.state.history.plastic.alpha > loaded.state.history.plastic.alpha);
    for i in [2, 4, 5] {
        near(next.full_stress_pa[i], 0.0, 1e-11);
    }
    near(
        material.yield_function_with_state(&next.full_stress_pa, &next.state.history),
        0.0,
        1e-10,
    );
    for j in 0..3 {
        let h = 1e-7;
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += h;
        minus[j] -= h;
        let above = plane.evaluate(plus, &loaded.state).unwrap();
        let below = plane.evaluate(minus, &loaded.state).unwrap();
        for i in 0..3 {
            near(
                next.tangent_pa[i][j],
                (above.stress_pa[i] - below.stress_pa[i]) / (2.0 * h),
                4e-6,
            );
        }
    }
}

#[test]
fn g3_legitimate_reversal_does_not_apply_elastic_domain_to_uncorrected_predictor() {
    let material = law();
    let load = deviatoric(0.038);
    let first = material
        .advance_with_energy(&[0.0; 6], &material.initial_state(), &load, TOL)
        .unwrap();
    let unload = deviatoric(-0.038);
    let trial: Voigt = core::array::from_fn(|i| unload[i] - first.state.plastic.plastic_strain[i]);
    assert!(contract(&trial, &trial).sqrt() > 0.05);
    let reversed = material
        .advance_with_energy(&load, &first.state, &unload, TOL)
        .unwrap();
    assert!(reversed.state.plastic.alpha > first.state.plastic.alpha);
    assert!(reversed.state.plastic.plastic_strain[0] < 0.0);
}

#[test]
fn g0_invalid_parameters_and_unreachable_independent_histories_refuse() {
    for c in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(BackstressBranch::new(c, 150.0).is_err());
    }
    for gamma in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(BackstressBranch::new(800.0, gamma).is_err());
    }
    assert!(
        ChabocheJ2::<9>::new(
            backbone(60.0),
            [BackstressBranch::new(1.0, 1.0).unwrap(); 9]
        )
        .is_err()
    );
    let material = law();
    let virgin = material.initial_state();
    let mut bad = virgin.clone();
    bad.plastic.alpha = 0.1;
    bad.backstress_pa[0] = deviatoric(100.0);
    let retained = bad.clone();
    assert!(material.response(&[0.0; 6], &bad).is_err());
    assert_eq!(bad, retained);
    bad.backstress_pa[0] = [0.01; 6];
    assert!(material.response(&[0.0; 6], &bad).is_err());
    assert!(material.response(&[f64::NAN; 6], &virgin).is_err());
    assert!(material.response(&[0.1; 6], &virgin).is_err());
    let disabled =
        ChabocheJ2::new(backbone(60.0), [BackstressBranch::new(0.0, 1.0).unwrap()]).unwrap();
    let disabled_bad = ChabocheState {
        plastic: bad.plastic,
        backstress_pa: [deviatoric(0.001)],
    };
    assert!(disabled.response(&[0.0; 6], &disabled_bad).is_err());
}
