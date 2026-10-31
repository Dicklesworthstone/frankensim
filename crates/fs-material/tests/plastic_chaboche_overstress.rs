//! Independent G0/G1 references for joint cyclic recovery and physical kinetics.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::rate::ChabocheOverstressJ2;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2, ChabocheState};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::VoceOverstressJ2;
use fs_material::plastic::voce::VoceJ2;
use fs_material::tensor::{contract, trace};
use fs_material::{SmallStrainLaw, Voigt};

const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};
const C: [f64; 2] = [300.0, 120.0];
const GAMMA: [f64; 2] = [150.0, 20.0];
fn voce(q: f64) -> VoceJ2 {
    VoceJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(),
            3.0,
            40.0,
            160.0,
        )
        .unwrap(),
        q,
        180.0,
    )
    .unwrap()
}
fn native() -> ChabocheJ2<2> {
    ChabocheJ2::new(
        voce(5.0),
        core::array::from_fn(|i| BackstressBranch::new(C[i], GAMMA[i]).unwrap()),
    )
    .unwrap()
}
fn law(n: f64, dt: f64) -> ChabocheOverstressJ2<2> {
    ChabocheOverstressJ2::new_power_law(native(), 20.0, 4.0, n, dt).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}
fn dev(value: f64) -> Voigt {
    [value, -0.5 * value, -0.5 * value, 0.0, 0.0, 0.0]
}

#[test]
fn g1_joint_return_matches_prescribed_plastic_increment_and_each_branch() {
    for n in [1.0, 1.7, 4.0, 12.0] {
        let law = law(n, 0.2);
        for alpha in [1e-8, 0.0008, 0.004, 0.012] {
            let saturation = 5.0 * -(-180.0_f64 * alpha).exp_m1();
            let viscous = 4.0 * (25.0_f64 * alpha).powf(1.0 / n);
            let branch_force: f64 = (0..2)
                .map(|i| C[i] * alpha / (1.0 + GAMMA[i] * alpha))
                .sum();
            let strain = dev((3.0 + 1700.0 * alpha + saturation + branch_force + viscous) / 1500.0);
            let step = law
                .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
                .unwrap();
            near(step.state.plastic.alpha, alpha, 1e-16);
            near(step.state.plastic.plastic_strain[0], alpha, 1e-16);
            near(trace(&step.state.plastic.plastic_strain), 0.0, 1e-16);
            for i in 0..2 {
                near(
                    step.state.backstress_pa[i][0],
                    (2.0 / 3.0) * C[i] * alpha / (1.0 + GAMMA[i] * alpha),
                    2e-14,
                );
            }
            near(
                step.stress_pa[0],
                (2.0 / 3.0) * (3.0 + 200.0 * alpha + saturation + branch_force + viscous),
                2e-12,
            );
            near(
                law.backbone()
                    .yield_function_with_state(&step.stress_pa, &step.state),
                viscous,
                2e-12,
            );
            let response = law.response(&strain, &law.initial_state()).unwrap();
            assert_eq!(response.state, step.state);
            assert!(response.iterations <= 96);
        }
    }
}

fn tangent_check(law: &ChabocheOverstressJ2<2>, strain: Voigt, state: &ChabocheState<2>) -> f64 {
    let tangent = law.tangent(&strain, state);
    let mut asymmetry = 0.0_f64;
    for j in 0..6 {
        let (mut plus, mut minus) = (strain, strain);
        plus[j] += 1e-7;
        minus[j] -= 1e-7;
        let above = law.stress(&plus, state);
        let below = law.stress(&minus, state);
        for i in 0..6 {
            near(tangent[i][j], (above[i] - below[i]) / 2e-7, 5e-6);
            let wi = if i < 3 { 1.0 } else { 2.0 };
            let wj = if j < 3 { 1.0 } else { 2.0 };
            asymmetry = asymmetry.max((wi * tangent[i][j] - wj * tangent[j][i]).abs());
        }
    }
    asymmetry
}

#[test]
fn g1_nonproportional_tangent_retains_direction_drift_and_rate_derivative() {
    for n in [1.0, 2.7, 8.0] {
        let law = law(n, 0.2);
        assert!(!law.admissibility().tangent_symmetric);
        let virgin = law.initial_state();
        tangent_check(
            &law,
            [0.0001, -0.00005, -0.00005, 0.0002, 0.0, 0.0],
            &virgin,
        );
        let loaded = [0.012, -0.004, -0.008, 0.002, 0.001, -0.001];
        tangent_check(&law, loaded, &virgin);
        let state = law.update_state(&loaded, &virgin);
        let reverse = [-0.008, 0.003, 0.005, -0.004, 0.002, 0.003];
        assert!(law.update_state(&reverse, &state).plastic.alpha > state.plastic.alpha);
        assert!(
            tangent_check(&law, reverse, &state) > 0.01,
            "must exercise a nonsymmetric tangent"
        );
    }
}

#[test]
fn g1_cyclic_energy_separates_viscosity_recovery_storage_and_numerical_loss() {
    let law = law(2.7, 0.2);
    let mut state = law.initial_state();
    let mut old = [0.0; 6];
    for strain in [
        dev(0.012),
        dev(0.012),
        [-0.01, 0.005, 0.005, 0.004, 0.001, 0.002],
        [0.01, -0.004, -0.006, -0.005, 0.003, 0.001],
    ] {
        let step = law.advance_with_energy(&old, &state, &strain, TOL).unwrap();
        let delta = step.state.plastic.alpha - state.plastic.alpha;
        let viscous = if delta == 0.0 {
            0.0
        } else {
            4.0 * (25.0 * delta).powf(1.0 / 2.7)
        };
        let qstore = |p: f64| 5.0 * (p + (-180.0 * p).exp_m1() / 180.0);
        let dp: Voigt = core::array::from_fn(|i| {
            step.state.plastic.plastic_strain[i] - state.plastic.plastic_strain[i]
        });
        let de: Voigt = core::array::from_fn(|i| strain[i] - old[i] - dp[i]);
        let alpha = step.state.plastic.alpha;
        let mut branch_storage = 160.0
            * contract(
                &step.state.plastic.plastic_strain,
                &step.state.plastic.plastic_strain,
            )
            / 3.0;
        let mut numerical = 0.5 * contract(&de, &voce(5.0).backbone().elastic.stress(&de, &()))
            + 20.0 * delta * delta
            + 160.0 * contract(&dp, &dp) / 3.0
            + 5.0 * -(-180.0 * alpha).exp_m1() * delta
            - (qstore(alpha) - qstore(state.plastic.alpha));
        let mut recovery = 0.0;
        for i in 0..2 {
            let a = &step.state.backstress_pa[i];
            let da: Voigt = core::array::from_fn(|j| a[j] - state.backstress_pa[i][j]);
            branch_storage += 0.75 * contract(a, a) / C[i];
            numerical += 0.75 * contract(&da, &da) / C[i];
            recovery += 1.5 * GAMMA[i] * delta * contract(a, a) / C[i];
        }
        near(
            step.energy.final_state.hardening_j_m3,
            20.0 * alpha * alpha + qstore(alpha),
            2e-15,
        );
        near(
            step.energy.final_state.kinematic_hardening_j_m3,
            branch_storage,
            2e-15,
        );
        near(step.recovery_dissipation_j_m3, recovery, 2e-15);
        near(
            step.energy.physical_dissipation_j_m3,
            (3.0 + viscous) * delta + recovery,
            2e-14,
        );
        near(step.energy.numerical_dissipation_j_m3, numerical, 2e-14);
        near(step.energy.residual_j_m3, 0.0, 2e-13);
        state = step.state;
        old = strain;
    }
}

fn absent_reduction<const N: usize>() {
    let native =
        ChabocheJ2::new(voce(5.0), [BackstressBranch::new(0.0, 150.0).unwrap(); N]).unwrap();
    let actual = ChabocheOverstressJ2::new_power_law(native, 20.0, 4.0, 2.7, 0.2).unwrap();
    let expected = VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 4.0, 2.7, 0.2).unwrap();
    let mut state = actual.initial_state();
    let mut old = [0.0; 6];
    for strain in [dev(0.012), [-0.009, 0.004, 0.005, -0.003, 0.002, 0.001]] {
        assert_eq!(
            actual.stress(&strain, &state),
            expected.stress(&strain, &state.plastic)
        );
        assert_eq!(
            actual.tangent(&strain, &state),
            expected.tangent(&strain, &state.plastic)
        );
        let a = actual
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let b = expected
            .advance_with_energy(&old, &state.plastic, &strain, TOL)
            .unwrap();
        assert_eq!(a.state.plastic, b.state);
        assert_eq!(a.stress_pa, b.stress_pa);
        assert_eq!(a.energy, b.energy);
        assert_eq!(a.recovery_dissipation_j_m3, 0.0);
        assert_eq!(actual.card(), expected.card());
        state = a.state;
        old = strain;
    }
}

#[test]
fn g3_absent_branches_and_n_one_preserve_their_exact_native_reductions() {
    absent_reduction::<0>();
    absent_reduction::<2>();
    let linear = ChabocheOverstressJ2::new_linear(native(), 20.0, 0.2).unwrap();
    let power = ChabocheOverstressJ2::new_power_law(native(), 20.0, 7.0, 1.0, 0.2).unwrap();
    let mut state = linear.initial_state();
    let mut old = [0.0; 6];
    for strain in [dev(0.012), [-0.009, 0.004, 0.005, -0.003, 0.002, 0.001]] {
        assert_eq!(
            linear.response(&strain, &state),
            power.response(&strain, &state)
        );
        let step = linear
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        assert_eq!(
            Ok(step.clone()),
            power.advance_with_energy(&old, &state, &strain, TOL)
        );
        state = step.state;
        old = strain;
    }
}

#[test]
fn g1_fixed_strain_relaxation_converges_to_the_analytical_power_law() {
    // gamma=0 and Q=0 are admitted physical limits with two retained linear
    // backstresses. f_dot=-A f^2/(eta s0), A=3mu+Hiso+Hk+C1+C2.
    let time = 0.05;
    let a = 2120.0;
    let initial = 12.0;
    let expected_f = initial / (1.0 + a * initial * time / (20.0 * 4.0));
    let expected_alpha = (initial - expected_f) / a;
    let solve = |steps: usize| {
        let native =
            ChabocheJ2::new(voce(0.0), C.map(|c| BackstressBranch::new(c, 0.0).unwrap())).unwrap();
        let law = ChabocheOverstressJ2::new_power_law(native, 20.0, 4.0, 2.0, time / steps as f64)
            .unwrap();
        let strain = dev(0.01);
        let mut state = law.initial_state();
        let mut losses = 0.0;
        let mut first = 0.0;
        let mut last = 0.0;
        for i in 0..steps {
            let step = law
                .advance_with_energy(&strain, &state, &strain, TOL)
                .unwrap();
            near(step.energy.mechanical_work_j_m3, 0.0, 0.0);
            if i == 0 {
                first = step.energy.initial.total_j_m3();
            }
            last = step.energy.final_state.total_j_m3();
            losses +=
                step.energy.physical_dissipation_j_m3 + step.energy.numerical_dissipation_j_m3;
            state = step.state;
        }
        near(first - last - losses, 0.0, 2e-12);
        state.plastic.alpha
    };
    let errors = [20, 40, 80].map(|steps| (solve(steps) - expected_alpha).abs());
    assert!(errors[1] < 0.57 * errors[0]);
    assert!(errors[2] < 0.57 * errors[1]);
    assert!(errors[2] < 2e-5);
}

#[test]
fn g0_tiny_rate_increments_remain_resolved_without_stress_power_overflow() {
    let law = ChabocheOverstressJ2::new_power_law(native(), 20.0, 1e6, 8.0, 0.2).unwrap();
    let strain = dev(0.012);
    let step = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
        .unwrap();
    assert!(step.state.plastic.alpha > 0.0 && step.state.plastic.alpha < 1e-25);
    let f = law
        .backbone()
        .yield_function_with_state(&step.stress_pa, &step.state);
    near(
        step.state.plastic.alpha / (0.2 * (1e6 / 20.0) * (f / 1e6).powf(8.0)),
        1.0,
        2e-12,
    );
    assert!(
        law.tangent(&strain, &law.initial_state())
            .iter()
            .flatten()
            .all(|x| x.is_finite())
    );
}

#[test]
fn g0_invalid_kinetics_and_unreachable_history_refuse_transactionally() {
    for eta in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(ChabocheOverstressJ2::new_linear(native(), eta, 0.2).is_err());
    }
    assert!(ChabocheOverstressJ2::new_linear(native(), 20.0, 0.0).is_err());
    assert!(ChabocheOverstressJ2::new_power_law(native(), 20.0, 0.0, 2.0, 0.2).is_err());
    assert!(ChabocheOverstressJ2::new_power_law(native(), 20.0, 4.0, 0.9, 0.2).is_err());
    let law = law(2.7, 0.2);
    let retained = law.initial_state();
    let mut invalid = retained.clone();
    invalid.backstress_pa[0][0] = 1.0;
    assert!(law.response(&dev(0.01), &invalid).is_err());
    assert!(
        law.advance_with_energy(&[0.0; 6], &retained, &dev(0.06), TOL)
            .is_err()
    );
    assert_eq!(retained, law.initial_state());
}
