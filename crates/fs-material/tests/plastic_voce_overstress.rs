//! Joint-rate references, retained cyclic history, physical heat and G0/G1 tangents.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{LinearOverstressJ2, PowerLawOverstressJ2, VoceOverstressJ2};
use fs_material::plastic::voce::VoceJ2;
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::tensor::contract;
use fs_material::{SmallStrainLaw, Voigt};

const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};
fn native() -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(),
        3.0,
        40.0,
        160.0,
    )
    .unwrap()
}
fn voce(q: f64) -> VoceJ2 {
    VoceJ2::new(native(), q, 180.0).unwrap()
}
fn dev(value: f64) -> Voigt {
    [value, -0.5 * value, -0.5 * value, 0.0, 0.0, 0.0]
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tol {tolerance:e}"
    );
}

#[test]
fn g1_joint_return_recovers_independently_prescribed_rate_and_saturation() {
    for n in [1.0, 1.7, 4.0, 12.0] {
        let law = VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 4.0, n, 0.2).unwrap();
        for alpha in [1e-8, 0.0008, 0.004, 0.015] {
            let hardening = 5.0 * -(-180.0_f64 * alpha).exp_m1();
            let overstress = 4.0 * (20.0_f64 * alpha / (0.2 * 4.0)).powf(1.0 / n);
            let strain = dev((3.0 + 1700.0 * alpha + hardening + overstress) / 1500.0);
            let step = law
                .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
                .unwrap();
            near(step.state.alpha, alpha, 4e-17);
            near(step.state.plastic_strain[0], alpha, 4e-17);
            near(
                step.stress_pa[0],
                (2.0 / 3.0) * (3.0 + 200.0 * alpha + hardening + overstress),
                2e-12,
            );
            near(
                law.backbone()
                    .yield_function_with_state(&step.stress_pa, &step.state),
                overstress,
                2e-12,
            );
            near(
                step.energy.physical_dissipation_j_m3,
                (3.0 + overstress) * alpha,
                2e-14,
            );
            assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
        }
    }
}

#[test]
fn g3_zero_saturation_is_exact_original_kinetics_and_n_one_is_linear() {
    let linear = LinearOverstressJ2::new(native(), 20.0, 0.2).unwrap();
    let combined_linear = VoceOverstressJ2::new_linear(voce(0.0), 20.0, 0.2).unwrap();
    let power = PowerLawOverstressJ2::new(native(), 20.0, 4.0, 2.7, 0.2).unwrap();
    let combined_power = VoceOverstressJ2::new_power_law(voce(0.0), 20.0, 4.0, 2.7, 0.2).unwrap();
    let saturation_linear = VoceOverstressJ2::new_linear(voce(5.0), 20.0, 0.2).unwrap();
    let saturation_n_one = VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 7.0, 1.0, 0.2).unwrap();
    assert_eq!(combined_linear.card(), linear.card());
    assert_eq!(combined_power.card(), power.card());
    assert_eq!(combined_linear.stress_scale_pa(), None);
    assert_eq!(saturation_n_one.stress_scale_pa(), Some(7.0));
    let mut state = power.initial_state();
    let mut old = [0.0; 6];
    for strain in [
        dev(0.01),
        [-0.008, 0.003, 0.005, -0.004, 0.002, 0.003],
        [0.004, 0.001, -0.005, 0.006, -0.002, -0.001],
    ] {
        assert_eq!(
            combined_linear.stress(&strain, &state),
            linear.stress(&strain, &state)
        );
        assert_eq!(
            combined_linear.tangent(&strain, &state),
            linear.tangent(&strain, &state)
        );
        assert_eq!(
            combined_linear.advance_with_energy(&old, &state, &strain, TOL),
            linear.advance_with_energy(&old, &state, &strain, TOL)
        );
        assert_eq!(
            combined_power.stress(&strain, &state),
            power.stress(&strain, &state)
        );
        assert_eq!(
            combined_power.tangent(&strain, &state),
            power.tangent(&strain, &state)
        );
        let expected = power
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        assert_eq!(
            combined_power
                .advance_with_energy(&old, &state, &strain, TOL)
                .unwrap(),
            expected
        );
        assert_eq!(
            saturation_linear.stress(&strain, &state),
            saturation_n_one.stress(&strain, &state)
        );
        assert_eq!(
            saturation_linear.tangent(&strain, &state),
            saturation_n_one.tangent(&strain, &state)
        );
        assert_eq!(
            saturation_linear.advance_with_energy(&old, &state, &strain, TOL),
            saturation_n_one.advance_with_energy(&old, &state, &strain, TOL)
        );
        state = expected.state;
        old = strain;
    }
}

fn check_tangent(law: &VoceOverstressJ2, strain: Voigt, state: &J2State) {
    let tangent = law.tangent(&strain, state);
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += 1e-7;
        minus[j] -= 1e-7;
        let above = law.stress(&plus, state);
        let below = law.stress(&minus, state);
        for i in 0..6 {
            near(tangent[i][j], (above[i] - below[i]) / 2e-7, 4e-6);
            let wi = if i < 3 { 1.0 } else { 2.0 };
            let wj = if j < 3 { 1.0 } else { 2.0 };
            near(wi * tangent[i][j], wj * tangent[j][i], 2e-12);
        }
    }
}

#[test]
fn g1_consistent_tangent_tracks_nonproportional_loading_and_reversal() {
    for n in [1.0, 2.7, 8.0] {
        let law = VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 4.0, n, 0.2).unwrap();
        let virgin = law.initial_state();
        check_tangent(
            &law,
            [0.0001, -0.00005, -0.00005, 0.0002, 0.0, 0.0],
            &virgin,
        );
        let strain = [0.009, -0.004, -0.005, 0.002, 0.001, -0.002];
        check_tangent(&law, strain, &virgin);
        let loaded = law.update_state(&strain, &virgin);
        let reversed = [-0.008, 0.003, 0.005, -0.004, 0.002, 0.003];
        assert!(law.update_state(&reversed, &loaded).alpha > loaded.alpha);
        check_tangent(&law, reversed, &loaded);
    }
}

#[test]
fn g1_cyclic_work_retains_voce_storage_and_separates_viscous_heat() {
    let law = VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 4.0, 2.7, 0.2).unwrap();
    let mut state = law.initial_state();
    let mut old = [0.0; 6];
    for strain in [
        dev(0.012),
        dev(0.012),
        [-0.01, 0.005, 0.005, 0.004, 0.001, 0.002],
        [0.01, -0.004, -0.006, -0.005, 0.003, 0.001],
    ] {
        let step = law.advance_with_energy(&old, &state, &strain, TOL).unwrap();
        let delta = step.state.alpha - state.alpha;
        let alpha = step.state.alpha;
        let viscous = if delta == 0.0 {
            0.0
        } else {
            4.0 * (25.0 * delta).powf(1.0 / 2.7)
        };
        let saturation_storage = |p: f64| 5.0 * (p + (-180.0 * p).exp_m1() / 180.0);
        near(
            step.energy.final_state.hardening_j_m3,
            20.0 * alpha * alpha + saturation_storage(alpha),
            1e-15,
        );
        near(
            step.energy.physical_dissipation_j_m3,
            (3.0 + viscous) * delta,
            1e-14,
        );
        let delta_p: Voigt =
            core::array::from_fn(|i| step.state.plastic_strain[i] - state.plastic_strain[i]);
        let delta_e: Voigt = core::array::from_fn(|i| strain[i] - old[i] - delta_p[i]);
        let elastic_loss = 0.5 * contract(&delta_e, &native().elastic.stress(&delta_e, &()));
        let hardening_loss = 20.0 * delta * delta + (160.0 / 3.0) * contract(&delta_p, &delta_p);
        let saturation_loss = 5.0 * -(-180.0 * alpha).exp_m1() * delta
            - (saturation_storage(alpha) - saturation_storage(state.alpha));
        near(
            step.energy.numerical_dissipation_j_m3,
            elastic_loss + hardening_loss + saturation_loss,
            1e-14,
        );
        near(step.energy.residual_j_m3, 0.0, 1e-13);
        state = step.state;
        old = strain;
    }
}

#[test]
fn g1_constant_stress_creep_converges_to_closed_form_voce_solution() {
    // With H_iso=H_k=0, q=12 Pa held, sigma_y=3 Pa, Q=5 Pa:
    // eta alpha_dot = 9-5+5 exp(-b alpha). For z=exp(b alpha),
    // z_dot=(b/eta)(4 z+5), giving z=(9 exp(4 b t/eta)-5)/4.
    let end_time = 0.02;
    let reference = ((9.0 * (4.0_f64 * 180.0 * end_time / 20.0).exp() - 5.0) / 4.0).ln() / 180.0;
    let solve = |steps: usize| {
        let backbone =
            J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(), 3.0, 0.0).unwrap();
        let law = VoceOverstressJ2::new_linear(
            VoceJ2::new(backbone, 5.0, 180.0).unwrap(),
            20.0,
            end_time / steps as f64,
        )
        .unwrap();
        let mut state = law.initial_state();
        let mut old = dev(12.0 / 1500.0);
        for _ in 0..steps {
            let mut lower = old[0];
            let mut upper = 0.03;
            for _ in 0..52 {
                let middle = 0.5 * (lower + upper);
                if law.stress(&dev(middle), &state)[0] < 8.0 {
                    lower = middle;
                } else {
                    upper = middle;
                }
            }
            let strain = dev(0.5 * (lower + upper));
            let step = law.advance_with_energy(&old, &state, &strain, TOL).unwrap();
            near(step.stress_pa[0], 8.0, 1e-11);
            state = step.state;
            old = strain;
        }
        state.alpha
    };
    let errors = [20, 40, 80].map(|steps| (solve(steps) - reference).abs());
    assert!(errors[1] < 0.53 * errors[0]);
    assert!(errors[2] < 0.53 * errors[1]);
    assert!(errors[2] < 2e-5);
}

#[test]
fn g0_small_rate_increment_and_actual_duration_remain_resolved() {
    let law = VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 1e6, 8.0, 0.2).unwrap();
    let strain = dev(0.012);
    let step = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
        .unwrap();
    assert!(step.state.alpha > 0.0 && step.state.alpha < 1e-25);
    let f = law
        .backbone()
        .yield_function_with_state(&step.stress_pa, &step.state);
    let rate_delta = 0.2 * (1e6 / 20.0) * (f / 1e6).powf(8.0);
    near(step.state.alpha / rate_delta, 1.0, 5e-13);
    assert!(
        law.tangent(&strain, &law.initial_state())
            .iter()
            .flatten()
            .all(|x| x.is_finite())
    );
    let fast = VoceOverstressJ2::new_linear(voce(5.0), 20.0, 0.002).unwrap();
    let slow = VoceOverstressJ2::new_linear(voce(5.0), 20.0, 2.0).unwrap();
    assert!(
        slow.update_state(&strain, &slow.initial_state()).alpha
            > fast.update_state(&strain, &fast.initial_state()).alpha
    );
    assert!(
        slow.stress(&strain, &slow.initial_state())[0]
            < fast.stress(&strain, &fast.initial_state())[0]
    );
}

#[test]
fn g0_invalid_kinetics_and_history_refuse_without_mutating_state() {
    for eta in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(VoceOverstressJ2::new_linear(voce(5.0), eta, 0.2).is_err());
    }
    assert!(VoceOverstressJ2::new_linear(voce(5.0), 20.0, 0.0).is_err());
    assert!(VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 0.0, 2.0, 0.2).is_err());
    assert!(VoceOverstressJ2::new_power_law(voce(5.0), 20.0, 4.0, 0.9, 0.2).is_err());
    let law = VoceOverstressJ2::new_linear(voce(5.0), 20.0, 0.2).unwrap();
    let state = law.initial_state();
    let retained = state.clone();
    assert!(
        law.advance_with_energy(&[0.0; 6], &state, &dev(0.06), TOL)
            .is_err()
    );
    let mut invalid = state.clone();
    invalid.alpha = -0.001;
    assert!(
        law.advance_with_energy(&[0.0; 6], &invalid, &dev(0.01), TOL)
            .is_err()
    );
    assert_eq!(state, retained);
}
