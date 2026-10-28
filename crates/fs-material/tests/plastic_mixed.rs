//! Independent cyclic limits, tangent actions and work of linear mixed J2.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::LinearOverstressJ2;
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::tensor::{contract, deviator};
use fs_material::{SmallStrainLaw, Voigt};

const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};
fn law(hi: f64, hk: f64) -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(1200.0, 0.2, 0.04).unwrap(),
        3.0,
        hi,
        hk,
    )
    .unwrap()
}
fn shear(value: f64) -> Voigt {
    [0.0, 0.0, 0.0, value, 0.0, 0.0]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn g1_identical_monotonic_tangent_has_distinct_source_selected_reverse_yield() {
    let mixed = law(40.0, 160.0);
    let isotropic = law(200.0, 0.0);
    let strain = shear(0.01);
    let loaded = mixed
        .advance_with_energy(&[0.0; 6], &mixed.initial_state(), &strain, TOL)
        .unwrap();
    let iso = isotropic
        .advance_with_energy(&[0.0; 6], &isotropic.initial_state(), &strain, TOL)
        .unwrap();
    assert_eq!(loaded.stress_pa, iso.stress_pa);
    assert_eq!(loaded.state, iso.state);
    // Pure shear has tau=2*mu*epsilon_xy and q=sqrt(3)*|tau-a_xy|.
    // Its reverse elastic limit is the translated center minus yield radius.
    let shift = (2.0 / 3.0) * 160.0 * loaded.state.plastic_strain[3];
    let reverse_stress = shift - (3.0 + 40.0 * loaded.state.alpha) / 3.0_f64.sqrt();
    let reverse_strain = 0.01 + (reverse_stress - loaded.stress_pa[3]) / 1000.0;
    near(mixed.backstress_pa(&loaded.state)[3], shift, 1e-14);
    let elastic = mixed
        .advance_with_energy(&strain, &loaded.state, &shear(reverse_strain + 1e-6), TOL)
        .unwrap();
    assert_eq!(elastic.state, loaded.state);
    assert_eq!(elastic.energy.physical_dissipation_j_m3, 0.0);
    let crossed = mixed
        .advance_with_energy(&strain, &loaded.state, &shear(reverse_strain - 1e-6), TOL)
        .unwrap();
    let still_elastic = isotropic
        .advance_with_energy(&strain, &iso.state, &shear(reverse_strain - 1e-6), TOL)
        .unwrap();
    assert!(crossed.state.alpha > loaded.state.alpha);
    assert!(crossed.energy.physical_dissipation_j_m3 > 0.0);
    assert_eq!(still_elastic.state, iso.state);
    assert_eq!(still_elastic.energy.physical_dissipation_j_m3, 0.0);
    near(
        mixed.yield_function_with_state(&crossed.stress_pa, &crossed.state),
        0.0,
        1e-12,
    );
}

#[test]
fn g1_cyclic_backstress_energy_stays_recoverable_and_rate_heat_uses_actual_overstress() {
    let mixed = law(40.0, 160.0);
    let mut state = mixed.initial_state();
    let mut old = [0.0; 6];
    for strain in [shear(0.01), shear(0.009), shear(-0.012), shear(0.008)] {
        let step = mixed
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let p = deviator(&step.state.plastic_strain);
        let delta_p =
            core::array::from_fn(|i| step.state.plastic_strain[i] - state.plastic_strain[i]);
        let delta_alpha = step.state.alpha - state.alpha;
        near(
            step.energy.final_state.kinematic_hardening_j_m3,
            160.0 / 3.0 * contract(&p, &p),
            1e-14,
        );
        near(
            step.energy.final_state.hardening_j_m3,
            20.0 * step.state.alpha.powi(2),
            1e-14,
        );
        near(
            step.energy.physical_dissipation_j_m3,
            3.0 * delta_alpha,
            1e-14,
        );
        let delta_e = core::array::from_fn(|i| strain[i] - old[i] - delta_p[i]);
        let de_stress = mixed.elastic.stress(&delta_e, &());
        let expected_numerical = 0.5 * contract(&delta_e, &de_stress)
            + 20.0 * delta_alpha.powi(2)
            + 160.0 / 3.0 * contract(&delta_p, &delta_p);
        near(
            step.energy.numerical_dissipation_j_m3,
            expected_numerical,
            1e-13,
        );
        near(
            step.energy.mechanical_work_j_m3,
            step.energy.stored_change_j_m3()
                + step.energy.physical_dissipation_j_m3
                + step.energy.numerical_dissipation_j_m3,
            1e-12,
        );
        state = step.state;
        old = strain;
    }
    let rate = LinearOverstressJ2::new(mixed.clone(), 30.0, 0.2).unwrap();
    let next = rate
        .advance_with_energy(&old, &state, &shear(-0.01), TOL)
        .unwrap();
    let delta_alpha = next.state.alpha - state.alpha;
    assert!(delta_alpha > 0.0);
    near(
        mixed.yield_function_with_state(&next.stress_pa, &next.state),
        150.0 * delta_alpha,
        1e-11,
    );
    near(
        next.energy.physical_dissipation_j_m3,
        3.0 * delta_alpha + 150.0 * delta_alpha.powi(2),
        1e-13,
    );
    assert!(next.energy.final_state.kinematic_hardening_j_m3 > 0.0);
}

fn tangent_check(law: &impl SmallStrainLaw<State = J2State>, state: &J2State, strain: Voigt) {
    let tangent = law.tangent(&strain, state);
    let h = 1e-7;
    for direction in 0..6 {
        let mut minus = strain;
        let mut plus = strain;
        minus[direction] -= h;
        plus[direction] += h;
        let a = law.stress(&minus, state);
        let b = law.stress(&plus, state);
        for component in 0..6 {
            near(
                tangent[component][direction],
                (b[component] - a[component]) / (2.0 * h),
                2e-5,
            );
        }
    }
}

#[test]
fn g1_returned_tangent_tracks_nonproportional_and_reversed_mixed_flow() {
    let mixed = law(40.0, 160.0);
    let state = mixed.update_state(
        &[0.008, -0.004, -0.004, 0.002, 0.0, 0.0],
        &mixed.initial_state(),
    );
    let strain = [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001];
    assert!(mixed.update_state(&strain, &state).alpha > state.alpha);
    tangent_check(&mixed, &state, strain);
    tangent_check(
        &LinearOverstressJ2::new(mixed, 30.0, 0.2).unwrap(),
        &state,
        strain,
    );
}

#[test]
fn g0_zero_kinematic_limit_preserves_original_updates_and_invalid_data_refuse() {
    let original =
        J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.04).unwrap(), 3.0, 40.0).unwrap();
    let zero = law(40.0, 0.0);
    let mut state = original.initial_state();
    let mut old = [0.0; 6];
    for strain in [shear(0.01), shear(0.008), shear(-0.012)] {
        assert_eq!(
            original.stress(&strain, &state),
            zero.stress(&strain, &state)
        );
        assert_eq!(
            original.tangent(&strain, &state),
            zero.tangent(&strain, &state)
        );
        let a = original
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        let b = zero
            .advance_with_energy(&old, &state, &strain, TOL)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(a.energy.final_state.kinematic_hardening_j_m3, 0.0);
        state = a.state;
        old = strain;
    }
    for hk in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(J2Plasticity::new_mixed(original.elastic.clone(), 3.0, 40.0, hk).is_err());
    }
    assert!(J2Plasticity::new_mixed(original.elastic.clone(), 3.0, f64::MAX, f64::MAX).is_err());
}
