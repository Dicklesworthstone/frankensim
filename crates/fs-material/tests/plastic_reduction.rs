//! G0/G1/G3 constraints, tangent and cyclic energy of native J2 reductions.
use fs_material::SmallStrainLaw;
use fs_material::elastic::IsotropicElastic;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::PowerLawOverstressJ2;
use fs_material::plastic::reduction::{
    PlaneStressReduction, PlaneStressWorkspace, ReducedState, ReductionConfig, UniaxialJ2,
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
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, allowance {tolerance}"
    );
}

#[test]
fn g1_free_axial_elastic_and_bilinear_limits_include_lateral_plastic_flow() {
    let law = UniaxialJ2::new(backbone()).unwrap();
    let virgin = law.initial_state();
    let elastic = law.response(0.001, &virgin);
    near(elastic.stress_pa, 1.2, 1e-14);
    near(elastic.tangent_pa, 1200.0, 0.0);
    near(elastic.strain[1], -0.0002, 1e-17);
    assert_eq!(elastic.state, virgin);
    let plastic = law.response(0.01, &virgin);
    let delta = (1200.0 * 0.01 - 3.0) / 1400.0;
    near(plastic.state.plastic_strain(), delta, 1e-17);
    near(plastic.state.equivalent_plastic_strain(), delta, 1e-17);
    near(plastic.stress_pa, 3.0 + 200.0 * delta, 1e-14);
    near(plastic.tangent_pa, 1200.0 * 200.0 / 1400.0, 1e-12);
    near(
        plastic.strain[1],
        -0.2 * plastic.stress_pa / 1200.0 - 0.5 * delta,
        1e-17,
    );
    let tensor = law
        .backbone()
        .stress(&plastic.strain, &virgin.tensor_history());
    near(tensor[0], plastic.stress_pa, 1e-13);
    for &stress in &tensor[1..] {
        near(stress, 0.0, 1e-13);
    }
    let history = law
        .backbone()
        .update_state(&plastic.strain, &virgin.tensor_history());
    near(history.alpha, delta, 1e-17);
    near(history.plastic_strain[0], delta, 1e-17);
    near(history.plastic_strain[1], -0.5 * delta, 1e-17);
}

#[test]
fn g1_axial_reversal_preserves_bauschinger_shift_and_consistent_tangent() {
    let law = UniaxialJ2::new(backbone()).unwrap();
    let loaded = law.response(0.01, &law.initial_state());
    let p = loaded.state.plastic_strain();
    let alpha = loaded.state.equivalent_plastic_strain();
    let reverse_yield_stress = 160.0 * p - (3.0 + 40.0 * alpha);
    let reverse_yield_strain = p + reverse_yield_stress / 1200.0;
    let inside = law.response(reverse_yield_strain + 1e-5, &loaded.state);
    assert_eq!(inside.state, loaded.state);
    let outside = law.response(reverse_yield_strain - 1e-5, &loaded.state);
    assert!(outside.state.equivalent_plastic_strain() > alpha);
    assert!(outside.state.plastic_strain() < p);
    for strain in [0.012, 0.009, -0.01] {
        let h = 1e-7;
        let difference = (law.stress(strain + h, &loaded.state)
            - law.stress(strain - h, &loaded.state))
            / (2.0 * h);
        near(law.tangent(strain, &loaded.state), difference, 1e-7);
        let response = law.response(strain, &loaded.state);
        let tensor = law
            .backbone()
            .stress(&response.strain, &loaded.state.tensor_history());
        near(tensor[0], response.stress_pa, 1e-13);
        for &stress in &tensor[1..] {
            near(stress, 0.0, 1e-13);
        }
    }
}

#[test]
fn g1_axial_cycles_use_native_energy_and_never_heat_with_recoverable_hardening() {
    let law = UniaxialJ2::new(backbone()).unwrap();
    let tolerances = EnergyTolerances {
        energy_j_m3: 1e-10,
        yield_pa: 1e-10,
    };
    let mut history = law.initial_state();
    let mut old_strain = 0.0;
    for strain in [0.01, 0.009, -0.012, 0.008] {
        let next = law.response(strain, &history);
        let energy = law
            .advance_with_energy(old_strain, &history, strain, tolerances)
            .unwrap();
        let alpha = next.state.equivalent_plastic_strain();
        let p = next.state.plastic_strain();
        near(
            energy.energy.mechanical_work_j_m3,
            next.stress_pa * (strain - old_strain),
            1e-13,
        );
        near(
            energy.energy.final_state.elastic_j_m3,
            600.0 * (strain - p).powi(2),
            1e-13,
        );
        near(
            energy.energy.final_state.hardening_j_m3,
            20.0 * alpha * alpha,
            1e-13,
        );
        near(
            energy.energy.final_state.kinematic_hardening_j_m3,
            80.0 * p * p,
            1e-13,
        );
        near(
            energy.energy.physical_dissipation_j_m3,
            3.0 * (alpha - history.equivalent_plastic_strain()),
            1e-13,
        );
        near(energy.state.alpha, alpha, 1e-17);
        near(energy.state.plastic_strain[0], p, 1e-17);
        near(energy.energy.residual_j_m3, 0.0, 1e-13);
        old_strain = strain;
        history = next.state;
    }
    assert!(
        law.advance_with_energy(old_strain, &history, 0.2, tolerances)
            .is_err()
    );
    assert!(
        law.advance_with_energy(old_strain, &history, f64::NAN, tolerances)
            .is_err()
    );
}

#[test]
fn g0_axial_perfect_plasticity_has_zero_loading_tangent_and_elastic_unloading() {
    let material = J2Plasticity::new(backbone().elastic, 3.0, 0.0).unwrap();
    let law = UniaxialJ2::new(material).unwrap();
    let loaded = law.response(0.01, &law.initial_state());
    near(loaded.stress_pa, 3.0, 1e-14);
    assert_eq!(loaded.tangent_pa, 0.0);
    let unloaded = law.response(0.009, &loaded.state);
    near(unloaded.stress_pa, 1.8, 1e-14);
    assert_eq!(unloaded.tangent_pa, 1200.0);
    assert_eq!(loaded.state, unloaded.state);
}

fn config() -> ReductionConfig {
    ReductionConfig {
        max_iterations: 32,
        stress_tolerance_pa: 1e-11,
        strain_limit: 0.05,
    }
}

#[test]
fn g1_plane_stress_elastic_limit_is_condensed_not_plane_strain() {
    let reduced = PlaneStressReduction::new(backbone(), config()).unwrap();
    let in_plane = [0.001, -0.0005, 0.0002];
    let next = reduced
        .evaluate(in_plane, &reduced.initial_state())
        .unwrap();
    let factor = 1200.0 / (1.0 - 0.2 * 0.2);
    near(
        next.stress_pa[0],
        factor * (in_plane[0] + 0.2 * in_plane[1]),
        1e-13,
    );
    near(
        next.stress_pa[1],
        factor * (in_plane[1] + 0.2 * in_plane[0]),
        1e-13,
    );
    near(next.stress_pa[2], 1000.0 * in_plane[2], 1e-13);
    near(
        next.state.strain[2],
        -0.2 / 0.8 * (in_plane[0] + in_plane[1]),
        1e-16,
    );
    near(next.tangent_pa[0][0], factor, 1e-12);
    near(next.tangent_pa[0][1], 0.2 * factor, 1e-12);
    near(next.tangent_pa[2][2], 1000.0, 1e-12);
    for i in [2, 4, 5] {
        near(next.full_stress_pa[i], 0.0, config().stress_tolerance_pa);
    }
    assert_eq!(next.state.history, reduced.initial_state().history);
}

#[test]
fn g1_plane_stress_agrees_with_independent_axial_return_and_free_contraction() {
    let axial = UniaxialJ2::new(backbone()).unwrap();
    let plane = PlaneStressReduction::new(backbone(), config()).unwrap();
    let mut axial_state = axial.initial_state();
    let mut plane_state = plane.initial_state();
    for strain in [0.01, 0.009, -0.012, 0.008] {
        let reference = axial.response(strain, &axial_state);
        let next = plane
            .evaluate([strain, reference.strain[1], 0.0], &plane_state)
            .unwrap();
        near(next.stress_pa[0], reference.stress_pa, 3e-11);
        near(next.stress_pa[1], 0.0, 3e-11);
        near(next.state.strain[2], reference.strain[2], 1e-13);
        near(
            next.state.history.alpha,
            reference.state.equivalent_plastic_strain(),
            1e-13,
        );
        near(
            next.state.history.plastic_strain[0],
            reference.state.plastic_strain(),
            1e-13,
        );
        let tangent = next.tangent_pa[0][0]
            - next.tangent_pa[0][1] * next.tangent_pa[1][0] / next.tangent_pa[1][1];
        near(tangent, reference.tangent_pa, 2e-7);
        axial_state = reference.state;
        plane_state = next.state;
    }
}

#[test]
fn g0_plane_stress_nonproportional_cycles_have_consistent_tensor_tangent() {
    let plane = PlaneStressReduction::new(backbone(), config()).unwrap();
    let mut state = plane.initial_state();
    let mut workspace = PlaneStressWorkspace::new().unwrap();
    for strain in [
        [0.008, -0.002, 0.004],
        [-0.006, 0.008, -0.003],
        [0.005, 0.004, 0.002],
    ] {
        let before = state.clone();
        let next = plane
            .evaluate_with_workspace(strain, &state, &mut workspace)
            .unwrap();
        assert_eq!(state, before);
        assert!(next.state.history.alpha >= state.history.alpha);
        for i in [2, 4, 5] {
            near(next.full_stress_pa[i], 0.0, config().stress_tolerance_pa);
        }
        for j in 0..3 {
            let h = 1e-7;
            let mut plus = strain;
            plus[j] += h;
            let mut minus = strain;
            minus[j] -= h;
            let positive = plane
                .evaluate_with_workspace(plus, &state, &mut workspace)
                .unwrap();
            let negative = plane
                .evaluate_with_workspace(minus, &state, &mut workspace)
                .unwrap();
            for i in 0..3 {
                near(
                    next.tangent_pa[i][j],
                    (positive.stress_pa[i] - negative.stress_pa[i]) / (2.0 * h),
                    2e-5,
                );
            }
        }
        // Matrix symmetry uses the tensor-shear work metric, not engineering shear.
        near(next.tangent_pa[0][2], 2.0 * next.tangent_pa[2][0], 1e-10);
        near(next.tangent_pa[1][2], 2.0 * next.tangent_pa[2][1], 1e-10);
        state = next.state;
    }
}

#[test]
fn g3_plane_stress_viscoplastic_iterations_advance_one_physical_interval() {
    let rate = PowerLawOverstressJ2::new(backbone(), 30.0, 4.0, 2.5, 0.2).unwrap();
    let plane = PlaneStressReduction::new(rate.clone(), config()).unwrap();
    let old = plane.initial_state();
    let strain = [0.008, 0.003, 0.004];
    let next = plane.evaluate(strain, &old).unwrap();
    assert!(next.iterations > 1);
    assert_eq!(
        next.state.history,
        rate.update_state(&next.state.strain, &old.history)
    );
    assert_eq!(
        next.full_stress_pa,
        rate.stress(&next.state.strain, &old.history)
    );
    assert_eq!(next, plane.evaluate(strain, &old).unwrap());
    // Repeating constitutive time at the accepted strain would relax again.
    assert!(
        rate.update_state(&next.state.strain, &next.state.history)
            .alpha
            > next.state.history.alpha
    );
    let hold = plane.evaluate(strain, &next.state).unwrap();
    assert!(hold.state.history.alpha > next.state.history.alpha);
    for i in [2, 4, 5] {
        near(hold.full_stress_pa[i], 0.0, config().stress_tolerance_pa);
    }
}

#[test]
fn g3_plane_stress_releases_transverse_shear_without_erasing_existing_history() {
    let native = backbone();
    let prior_strain = [0.005, -0.002, -0.003, 0.001, 0.003, -0.002];
    let history = native.update_state(&prior_strain, &native.initial_state());
    assert!(history.plastic_strain[4] != 0.0 && history.plastic_strain[5] != 0.0);
    native
        .validate_committed_history(&prior_strain, &history, 1e-10)
        .unwrap();
    let old = ReducedState {
        strain: prior_strain,
        history,
    };
    let plane = PlaneStressReduction::new(native, config()).unwrap();
    let next = plane.evaluate([0.004, -0.001, 0.002], &old).unwrap();
    for i in [2, 4, 5] {
        near(next.full_stress_pa[i], 0.0, config().stress_tolerance_pa);
    }
    assert!(next.state.strain[4] != 0.0 && next.state.strain[5] != 0.0);
    assert!(next.state.history.alpha >= old.history.alpha);
    assert_eq!(
        next.state.history,
        plane.law().update_state(&next.state.strain, &old.history)
    );
}

#[test]
fn g0_plane_stress_refusals_leave_accepted_history_unchanged() {
    let plane = PlaneStressReduction::new(backbone(), config()).unwrap();
    let old = plane
        .evaluate([0.008, -0.002, 0.003], &plane.initial_state())
        .unwrap()
        .state;
    let saved = old.clone();
    for strain in [[f64::NAN, 0.0, 0.0], [0.1, 0.0, 0.0]] {
        assert!(plane.evaluate(strain, &old).is_err());
        assert_eq!(old, saved);
    }
    let one = PlaneStressReduction::new(
        backbone(),
        ReductionConfig {
            max_iterations: 1,
            ..config()
        },
    )
    .unwrap();
    assert!(
        one.evaluate([0.008, 0.004, 0.004], &one.initial_state())
            .is_err()
    );
    assert_eq!(old, saved);
    for bad in [
        ReductionConfig {
            max_iterations: 0,
            ..config()
        },
        ReductionConfig {
            stress_tolerance_pa: 0.0,
            ..config()
        },
        ReductionConfig {
            strain_limit: 0.1,
            ..config()
        },
    ] {
        assert!(PlaneStressReduction::new(backbone(), bad).is_err());
    }
}
