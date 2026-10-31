//! Independent saturation curves, tensor tangents and energetic cyclic history.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::{PlaneStressReduction, ReductionConfig};
use fs_material::plastic::voce::VoceJ2;
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::tensor::{contract, deviator, trace};
use fs_material::{SmallStrainLaw, Voigt};

const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};

fn backbone() -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(),
        3.0,
        40.0,
        160.0,
    )
    .unwrap()
}

fn law() -> VoceJ2 {
    VoceJ2::new(backbone(), 5.0, 180.0).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; allowance {tolerance}"
    );
}

fn uniaxial_deviator(value: f64) -> Voigt {
    [value, -0.5 * value, -0.5 * value, 0.0, 0.0, 0.0]
}

#[test]
fn g1_prescribed_plastic_strain_recovers_independent_nonlinear_yield_curve() {
    let material = law();
    let virgin = material.initial_state();
    // Prescribe alpha, then derive total strain independently from the yield
    // equation. Under this deviatoric loading epsilon_p_xx = alpha and
    // q_trial = 3 mu epsilon_xx, with mu = 500 Pa.
    for alpha in [1e-8, 0.0008, 0.004, 0.02] {
        let saturation = 5.0 * -(-180.0_f64 * alpha).exp_m1();
        let strain = uniaxial_deviator((3.0 + 1700.0 * alpha + saturation) / 1500.0);
        let step = material
            .advance_with_energy(&[0.0; 6], &virgin, &strain, TOL)
            .unwrap();
        near(step.state.alpha, alpha, 2e-17);
        near(step.state.plastic_strain[0], alpha, 2e-17);
        near(step.state.plastic_strain[1], -0.5 * alpha, 2e-17);
        near(step.state.plastic_strain[2], -0.5 * alpha, 2e-17);
        near(trace(&step.state.plastic_strain), 0.0, 2e-17);
        // The physical deviator includes both the expanded yield radius and
        // the independent translated backstress.
        near(
            step.stress_pa[0],
            (2.0 / 3.0) * (3.0 + 200.0 * alpha + saturation),
            3e-13,
        );
        near(
            material.yield_function_with_state(&step.stress_pa, &step.state),
            0.0,
            3e-13,
        );
        near(
            material.isotropic_hardening_derivative_pa(alpha),
            40.0 + 900.0 * (-180.0_f64 * alpha).exp(),
            2e-12,
        );
    }
}

#[test]
fn g0_zero_saturation_is_bitwise_original_through_nonproportional_cycles() {
    let original = backbone();
    let zero = VoceJ2::new(original.clone(), 0.0, 180.0).unwrap();
    let mut state = original.initial_state();
    let mut old_strain = [0.0; 6];
    for strain in [
        [0.0001, -0.00005, -0.00005, 0.0002, 0.0, 0.0],
        [0.008, -0.004, -0.004, 0.002, 0.0, 0.0],
        [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
        [0.01, -0.004, -0.006, 0.005, -0.002, 0.0],
    ] {
        assert_eq!(
            zero.stress(&strain, &state),
            original.stress(&strain, &state)
        );
        assert_eq!(
            zero.tangent(&strain, &state),
            original.tangent(&strain, &state)
        );
        assert_eq!(
            zero.update_state(&strain, &state),
            original.update_state(&strain, &state)
        );
        let expected = original
            .advance_with_energy(&old_strain, &state, &strain, TOL)
            .unwrap();
        let actual = zero
            .advance_with_energy(&old_strain, &state, &strain, TOL)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(zero.card(), original.card());
        state = expected.state;
        old_strain = strain;
    }
}

fn check_tangent(material: &VoceJ2, strain: Voigt, state: &J2State) {
    let tangent = material.tangent(&strain, state);
    let h = 1e-7;
    for direction in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[direction] += h;
        minus[direction] -= h;
        let above = material.stress(&plus, state);
        let below = material.stress(&minus, state);
        for component in 0..6 {
            near(
                tangent[component][direction],
                (above[component] - below[component]) / (2.0 * h),
                3e-6,
            );
            let weight_i = if component < 3 { 1.0 } else { 2.0 };
            let weight_j = if direction < 3 { 1.0 } else { 2.0 };
            near(
                weight_i * tangent[component][direction],
                weight_j * tangent[direction][component],
                1e-12,
            );
        }
    }
}

#[test]
fn g1_endpoint_tangent_tracks_reversed_and_nonproportional_saturation() {
    let material = law();
    let virgin = material.initial_state();
    check_tangent(
        &material,
        [0.0002, -0.0001, -0.0001, 0.0003, 0.0, 0.0],
        &virgin,
    );
    let strain = [0.009, -0.004, -0.005, 0.002, 0.001, -0.002];
    check_tangent(&material, strain, &virgin);
    let loaded = material.update_state(&strain, &virgin);
    let reversed = [-0.008, 0.003, 0.005, -0.004, 0.002, 0.003];
    assert!(material.update_state(&reversed, &loaded).alpha > loaded.alpha);
    check_tangent(&material, reversed, &loaded);
}

#[test]
fn g1_cyclic_storage_and_convex_loss_do_not_become_physical_plastic_heat() {
    let material = law();
    let mut old = [0.0; 6];
    let mut state = material.initial_state();
    for strain in [
        [0.008, -0.004, -0.004, 0.002, 0.0, 0.0],
        [0.007, -0.0035, -0.0035, 0.0015, 0.0, 0.0],
        [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
        [0.01, -0.004, -0.006, 0.005, -0.002, 0.0],
    ] {
        let step = material
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let alpha = step.state.alpha;
        let delta = alpha - state.alpha;
        let storage = |p: f64| 5.0 * (p + (-180.0 * p).exp_m1() / 180.0);
        let sat_work = 5.0 * -(-180.0 * alpha).exp_m1() * delta;
        let sat_loss = sat_work - (storage(alpha) - storage(state.alpha));
        let delta_p: Voigt =
            core::array::from_fn(|i| step.state.plastic_strain[i] - state.plastic_strain[i]);
        let delta_e = core::array::from_fn(|i| strain[i] - old[i] - delta_p[i]);
        let de_stress = material.backbone().elastic.stress(&delta_e, &());
        let p = deviator(&step.state.plastic_strain);
        near(
            step.energy.final_state.hardening_j_m3,
            20.0 * alpha.powi(2) + storage(alpha),
            1e-14,
        );
        near(
            step.energy.final_state.kinematic_hardening_j_m3,
            160.0 / 3.0 * contract(&p, &p),
            1e-14,
        );
        near(step.energy.physical_dissipation_j_m3, 3.0 * delta, 1e-14);
        near(
            step.energy.numerical_dissipation_j_m3,
            0.5 * contract(&delta_e, &de_stress)
                + 20.0 * delta.powi(2)
                + 160.0 / 3.0 * contract(&delta_p, &delta_p)
                + sat_loss,
            1e-13,
        );
        assert!(step.energy.physical_dissipation_j_m3 >= 0.0);
        assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
        near(step.energy.residual_j_m3, 0.0, 1e-13);
        assert!(alpha >= state.alpha);
        near(trace(&step.state.plastic_strain), 0.0, 2e-17);
        old = strain;
        state = step.state;
    }
}

#[test]
fn g1_saturated_limit_retains_linear_backbone_and_elastic_unloading() {
    let material = VoceJ2::new(backbone(), 5.0, 100_000.0).unwrap();
    let alpha = 0.01;
    let strain = uniaxial_deviator((8.0 + 1700.0 * alpha) / 1500.0);
    let step = material
        .advance_with_energy(&[0.0; 6], &material.initial_state(), &strain, TOL)
        .unwrap();
    near(step.state.alpha, alpha, 1e-17);
    near(
        material.isotropic_hardening_pa(alpha),
        40.0 * alpha + 5.0,
        0.0,
    );
    near(material.isotropic_hardening_derivative_pa(alpha), 40.0, 0.0);
    let unload = uniaxial_deviator(strain[0] - 0.0005);
    let next = material
        .advance_with_energy(&strain, &step.state, &unload, TOL)
        .unwrap();
    assert_eq!(next.state, step.state);
    assert_eq!(next.energy.physical_dissipation_j_m3, 0.0);
    assert_eq!(
        material.tangent(&unload, &step.state),
        backbone().elastic.stiffness()
    );
}

#[test]
fn g3_tiny_exponent_retains_representable_saturation_strength_and_storage() {
    let linear = J2Plasticity::new(backbone().elastic, 3.0, 0.0).unwrap();
    let material = VoceJ2::new(linear, 2e12, 1e-12).unwrap();
    let alpha = 0.002;
    let saturation = 2e12 * -(-1e-12_f64 * alpha).exp_m1();
    let strain = uniaxial_deviator((3.0 + 1500.0 * alpha + saturation) / 1500.0);
    let next = material
        .advance_with_energy(&[0.0; 6], &material.initial_state(), &strain, TOL)
        .unwrap();
    near(next.state.alpha, alpha, 2e-18);
    // Here b alpha = 2e-15; the independently known leading integral is
    // (Q b / 2) alpha^2 with a relative correction below 7e-16.
    near(next.energy.final_state.hardening_j_m3, alpha * alpha, 1e-20);
    near(material.isotropic_hardening_pa(alpha), 2.0 * alpha, 6e-18);
    near(next.energy.residual_j_m3, 0.0, 2e-17);
}

#[test]
fn g3_roundoff_at_retained_yield_surface_does_not_make_a_nonfinite_return() {
    let material = law();
    let alpha = 13.0 * 1e-6;
    let state = J2State {
        plastic_strain: uniaxial_deviator(alpha),
        alpha,
    };
    // This admissible boundary value differs by one ulp between the two
    // groupings sy+(H*alpha+R) and (sy+R)+H*alpha. Tensor and scalar decisions
    // must use the same grouping to avoid manufacturing a failed return.
    let strain = uniaxial_deviator(f64::from_bits(4566810045847161522));
    assert!(
        material
            .stress(&strain, &state)
            .iter()
            .all(|v| v.is_finite())
    );
    assert!(
        material
            .tangent(&strain, &state)
            .iter()
            .flatten()
            .all(|v| v.is_finite())
    );
    let returned = material
        .advance_with_energy(&strain, &state, &strain, TOL)
        .unwrap();
    near(returned.state.alpha, state.alpha, 1e-18);
}

#[test]
fn g1_plane_stress_reuses_nonlinear_history_and_condenses_actual_tangent() {
    let material = law();
    let reduced = PlaneStressReduction::new(
        material.clone(),
        ReductionConfig {
            max_iterations: 32,
            stress_tolerance_pa: 1e-11,
            strain_limit: 0.05,
        },
    )
    .unwrap();
    let loaded = reduced
        .evaluate([0.009, -0.002, 0.003], &reduced.initial_state())
        .unwrap();
    let strain = [-0.008, 0.003, -0.004];
    let returned = reduced.evaluate(strain, &loaded.state).unwrap();
    assert!(returned.state.history.alpha > loaded.state.history.alpha);
    for component in [2, 4, 5] {
        near(returned.full_stress_pa[component], 0.0, 1e-11);
    }
    near(
        material.yield_function_with_state(&returned.full_stress_pa, &returned.state.history),
        0.0,
        1e-10,
    );
    let h = 1e-7;
    for direction in 0..3 {
        let mut plus = strain;
        let mut minus = strain;
        plus[direction] += h;
        minus[direction] -= h;
        let above = reduced.evaluate(plus, &loaded.state).unwrap();
        let below = reduced.evaluate(minus, &loaded.state).unwrap();
        for component in 0..3 {
            near(
                returned.tangent_pa[component][direction],
                (above.stress_pa[component] - below.stress_pa[component]) / (2.0 * h),
                3e-6,
            );
        }
    }
}

#[test]
fn g0_invalid_parameters_strain_and_history_are_refused_without_mutation() {
    for saturation in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(VoceJ2::new(backbone(), saturation, 180.0).is_err());
    }
    for rate in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(VoceJ2::new(backbone(), 5.0, rate).is_err());
    }
    assert!(VoceJ2::new(backbone(), f64::MAX, 2.0).is_err());
    let material = law();
    let virgin = material.initial_state();
    for strain in [[f64::NAN; 6], [0.1; 6], [f64::MAX; 6]] {
        assert!(
            material
                .advance_with_energy(&[0.0; 6], &virgin, &strain, TOL)
                .is_err()
        );
    }
    for bad in [
        J2State {
            plastic_strain: [0.0; 6],
            alpha: -1.0,
        },
        J2State {
            plastic_strain: [0.001; 6],
            alpha: 0.01,
        },
        J2State {
            plastic_strain: uniaxial_deviator(0.01),
            alpha: 0.001,
        },
    ] {
        let retained = bad.clone();
        assert!(
            material
                .advance_with_energy(&[0.0; 6], &bad, &uniaxial_deviator(0.01), TOL)
                .is_err()
        );
        assert_eq!(bad, retained);
    }
    assert_eq!(virgin, material.initial_state());
}
