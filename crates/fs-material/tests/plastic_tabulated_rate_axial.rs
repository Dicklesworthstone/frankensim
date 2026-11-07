//! Independent axial kinetics, free-side condensation, cyclic history and work.
//! Synthetic coefficients test the native law, not an alloy calibration.

use fs_material::elastic::IsotropicElastic;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{LinearOverstressJ2, PowerLawOverstressJ2};
use fs_material::plastic::reduction::{UniaxialJ2State, UniaxialTabulated};
use fs_material::plastic::tabulated::rate::axial::UniaxialTabulatedOverstress;
use fs_material::plastic::tabulated::{
    HardeningPoint, TabulatedError, TabulatedHardening, TabulatedJ2,
};
use fs_material::{SmallStrainLaw, Tangent6};

const E: f64 = 1200.0;
const NU: f64 = 0.2;
const HK: f64 = 160.0;
const ETA: f64 = 40.0;
const S0: f64 = 2.0;
const DT: f64 = 0.2;
const KNOTS: [(f64, f64); 5] = [
    (0.0, 3.0),
    (0.002, 4.0),
    (0.006, 4.0),
    (0.012, 7.0),
    (0.04, 8.4),
];
const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};

fn native(hardening: f64, kinematic: f64) -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        3.0,
        hardening,
        kinematic,
    )
    .unwrap()
}
fn table(knots: &[(f64, f64)], kinematic: f64) -> TabulatedJ2 {
    TabulatedJ2::new(
        native(0.0, kinematic),
        TabulatedHardening::new(
            knots
                .iter()
                .map(|&(plastic_strain, yield_stress_pa)| HardeningPoint {
                    plastic_strain,
                    yield_stress_pa,
                })
                .collect(),
        )
        .unwrap(),
    )
    .unwrap()
}
fn law(exponent: f64, duration: f64) -> UniaxialTabulatedOverstress {
    UniaxialTabulatedOverstress::new_power_law(table(&KNOTS, HK), ETA, S0, exponent, duration)
        .unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}; difference {:.8e}",
        (actual - expected).abs()
    );
}
fn radius(alpha: f64) -> f64 {
    for pair in KNOTS.windows(2) {
        let [(a, ya), (b, yb)] = [pair[0], pair[1]];
        if alpha <= b {
            return ya + (yb - ya) * ((alpha - a) / (b - a));
        }
    }
    panic!("independent curve oracle outside support");
}
fn storage(alpha: f64) -> f64 {
    let mut result = 0.0;
    for pair in KNOTS.windows(2) {
        let [(a, ya), (b, yb)] = [pair[0], pair[1]];
        if alpha <= a {
            break;
        }
        let width = alpha.min(b) - a;
        result += (ya - 3.0) * width + 0.5 * (yb - ya) / (b - a) * width * width;
    }
    result
}
fn viscous(delta: f64, exponent: f64, duration: f64) -> f64 {
    S0 * (ETA * delta / (duration * S0)).powf(1.0 / exponent)
}
fn condensed(tangent: &Tangent6) -> f64 {
    // Eliminate both independent transverse normal strains from the actual
    // native 3-D Jacobian. C_xxxx alone is not the axial-stress stiffness.
    let determinant = tangent[1][1] * tangent[2][2] - tangent[1][2] * tangent[2][1];
    let y = (tangent[2][2] * tangent[1][0] - tangent[1][2] * tangent[2][0]) / determinant;
    let z = (tangent[1][1] * tangent[2][0] - tangent[2][1] * tangent[1][0]) / determinant;
    tangent[0][0] - tangent[0][1] * y - tangent[0][2] * z
}

#[test]
fn axial_rate_matches_original_tensor_return_and_its_condensed_tangent() {
    for exponent in [1.0, 1.7, 4.0] {
        let material = law(exponent, DT);
        let virgin = material.initial_state();
        for target in [1e-6, 0.001, 0.004, 0.009, 0.018] {
            let excess = viscous(target, exponent, DT);
            let strain = target + (radius(target) + HK * target + excess) / E;
            let axial = material.response(strain, &virgin).unwrap();
            let full = material
                .backbone()
                .response(&axial.strain, &virgin.tensor_history())
                .unwrap();
            near(
                axial.state.equivalent_plastic_strain(),
                target,
                3e-16,
                3e-13,
            );
            near(axial.state.plastic_strain(), target, 3e-16, 3e-13);
            near(
                axial.stress_pa,
                radius(target) + HK * target + excess,
                3e-12,
                3e-13,
            );
            near(full.stress_pa[0], axial.stress_pa, 3e-12, 3e-13);
            for &stress in &full.stress_pa[1..] {
                near(stress, 0.0, 3e-12, 0.0);
            }
            near(full.state.alpha, target, 3e-16, 3e-13);
            let history = axial.state.tensor_history();
            for (&actual, expected) in full.state.plastic_strain.iter().zip(history.plastic_strain)
            {
                near(actual, expected, 3e-16, 3e-13);
            }
            near(axial.tangent_pa, condensed(&full.tangent_pa), 2e-9, 2e-11);
            near(
                axial.strain[1],
                -NU * axial.stress_pa / E - 0.5 * target,
                3e-16,
                3e-13,
            );
            assert_eq!(axial.strain[1], axial.strain[2]);
            assert_eq!(material.stress(strain, &virgin), axial.stress_pa);
            assert_eq!(material.tangent(strain, &virgin), axial.tangent_pa);
            assert_eq!(material.update_state(strain, &virgin), axial.state);
            assert_eq!(virgin.equivalent_plastic_strain(), 0.0);
        }
    }
}

fn check_derivatives(material: &UniaxialTabulatedOverstress, strain: f64, state: &UniaxialJ2State) {
    let response = material.response(strain, state).unwrap();
    let h = 1e-7;
    let plus = material.response(strain + h, state).unwrap();
    let minus = material.response(strain - h, state).unwrap();
    near(
        response.tangent_pa,
        (plus.stress_pa - minus.stress_pa) / (2.0 * h),
        2e-5,
        3e-8,
    );
    let (stress, _) = material.stress_and_potential(strain, state).unwrap();
    let wp = material.stress_and_potential(strain + h, state).unwrap().1;
    let wm = material.stress_and_potential(strain - h, state).unwrap().1;
    near(stress, response.stress_pa, 3e-12, 3e-13);
    near(response.stress_pa, (wp - wm) / (2.0 * h), 3e-8, 3e-8);
}

#[test]
fn cyclic_prager_history_and_rate_potential_keep_the_frozen_old_state() {
    let original = law(2.4, DT);
    let mut state = original.initial_state();
    let mut previous = 0.0;
    for (strain, duration) in [
        (0.0005, 0.05),
        (0.011, 0.2),
        (0.009, 0.07),
        (-0.009, 0.3),
        (-0.007, 0.1),
        (0.01, 0.15),
    ] {
        let material = original.for_duration(duration).unwrap();
        let retained = state.clone();
        check_derivatives(&material, strain, &state);
        let response = material.response(strain, &state).unwrap();
        let energy = material
            .advance_with_energy(previous, &state, strain, TOL)
            .unwrap();
        assert_eq!(state, retained);
        let delta = response.state.equivalent_plastic_strain() - state.equivalent_plastic_strain();
        assert!(delta >= 0.0);
        near(
            (response.state.plastic_strain() - state.plastic_strain()).abs(),
            delta,
            3e-16,
            3e-13,
        );
        if delta > 0.0 {
            let effective = response.stress_pa - HK * response.state.plastic_strain();
            near(
                effective.abs(),
                radius(response.state.equivalent_plastic_strain()) + viscous(delta, 2.4, duration),
                3e-12,
                3e-12,
            );
        }
        near(energy.energy.residual_j_m3, 0.0, 3e-13, 0.0);
        previous = strain;
        state = response.state;
    }
    assert!(state.equivalent_plastic_strain() > state.plastic_strain().abs());

    // Reverse yield translates by the retained Prager backstress, independent
    // of accumulated table hardening. A small elastic excursion crosses none.
    let loaded = original
        .response(0.011, &original.initial_state())
        .unwrap()
        .state;
    let reverse_yield = loaded.plastic_strain()
        + (HK * loaded.plastic_strain() - radius(loaded.equivalent_plastic_strain())) / E;
    let elastic = original.response(reverse_yield + 1e-5, &loaded).unwrap();
    assert_eq!(elastic.state, loaded);
    assert_eq!(elastic.tangent_pa, E);
    let reversed = original.response(reverse_yield - 0.001, &loaded).unwrap();
    assert!(reversed.state.plastic_strain() < loaded.plastic_strain());
    assert!(reversed.state.equivalent_plastic_strain() > loaded.equivalent_plastic_strain());
}

#[test]
fn axial_work_retains_exact_table_storage_prager_storage_and_only_physical_heat() {
    for exponent in [1.0, 3.2] {
        let material = law(exponent, DT);
        let mut state = material.initial_state();
        let mut old_strain = 0.0;
        for strain in [0.009, 0.014, 0.006, -0.008, -0.01, 0.008] {
            let next = material.response(strain, &state).unwrap();
            let checked = material
                .advance_with_energy(old_strain, &state, strain, TOL)
                .unwrap();
            let old_alpha = state.equivalent_plastic_strain();
            let alpha = next.state.equivalent_plastic_strain();
            let delta = alpha - old_alpha;
            let plastic = next.state.plastic_strain();
            let delta_plastic = plastic - state.plastic_strain();
            let elastic = strain - plastic;
            let delta_elastic = elastic - (old_strain - state.plastic_strain());
            let stored = material.stored_energy(strain, &next.state).unwrap();
            near(
                stored.elastic_j_m3,
                0.5 * E * elastic * elastic,
                3e-14,
                3e-12,
            );
            near(stored.hardening_j_m3, storage(alpha), 3e-14, 3e-12);
            near(
                stored.kinematic_hardening_j_m3,
                0.5 * HK * plastic * plastic,
                3e-14,
                3e-12,
            );
            near(
                checked.energy.final_state.total_j_m3(),
                stored.total_j_m3(),
                3e-14,
                3e-12,
            );
            near(
                checked.energy.physical_dissipation_j_m3,
                (3.0 + viscous(delta, exponent, DT)) * delta,
                3e-14,
                3e-12,
            );
            let hardening_loss =
                (radius(alpha) - 3.0) * delta - (storage(alpha) - storage(old_alpha));
            let numerical = 0.5 * E * delta_elastic * delta_elastic
                + 0.5 * HK * delta_plastic * delta_plastic
                + hardening_loss;
            near(
                checked.energy.numerical_dissipation_j_m3,
                numerical,
                3e-14,
                3e-12,
            );
            near(
                checked.energy.mechanical_work_j_m3,
                next.stress_pa * (strain - old_strain),
                3e-14,
                3e-12,
            );
            assert!(checked.energy.physical_dissipation_j_m3 >= 0.0);
            assert!(checked.energy.numerical_dissipation_j_m3 >= 0.0);
            near(checked.energy.residual_j_m3, 0.0, 3e-13, 0.0);
            old_strain = strain;
            state = next.state;
        }
    }
}

#[test]
fn variable_physical_intervals_relax_free_sided_stress_without_resetting_history() {
    let original =
        UniaxialTabulatedOverstress::new_linear(table(&[(0.0, 3.0), (0.05, 3.0)], 0.0), ETA, DT)
            .unwrap();
    let strain = 0.012;
    let first = original
        .response(strain, &original.initial_state())
        .unwrap();
    let mut state = first.state;
    let mut stress = first.stress_pa;
    for duration in [0.05, 0.4, 0.12, 0.03, 0.2] {
        let step_law = original.for_duration(duration).unwrap();
        assert_eq!(step_law.duration_s(), duration);
        assert_eq!(step_law.backbone().stress_scale_pa(), None);
        assert!(step_law.source().is_none());
        assert!(step_law.material().is_none());
        step_law.validate_history(strain, &state).unwrap();
        let next = step_law.response(strain, &state).unwrap();
        let expected = 3.0 + (stress - 3.0) / (1.0 + E * duration / ETA);
        near(next.stress_pa, expected, 3e-12, 3e-13);
        assert!(next.stress_pa < stress);
        assert!(next.state.equivalent_plastic_strain() > state.equivalent_plastic_strain());
        let energy = step_law
            .advance_with_energy(strain, &state, strain, TOL)
            .unwrap()
            .energy;
        near(energy.mechanical_work_j_m3, 0.0, 3e-14, 0.0);
        assert!(energy.physical_dissipation_j_m3 > 0.0);
        assert!(energy.stored_change_j_m3() < 0.0);
        near(
            -energy.stored_change_j_m3(),
            energy.physical_dissipation_j_m3 + energy.numerical_dissipation_j_m3,
            3e-14,
            3e-12,
        );
        stress = next.stress_pa;
        state = next.state;
    }
    assert_eq!(original.duration_s(), DT);
    for duration in [0.001, 10.0, 1000.0] {
        let material = original.for_duration(duration).unwrap();
        let virgin = material.initial_state();
        let subyield = material.response(0.001, &virgin).unwrap();
        assert_eq!(subyield.state, virgin);
        assert_eq!(subyield.tangent_pa, E);
        near(subyield.stress_pa, 1.2, 1e-15, 0.0);
    }
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(original.for_duration(duration).is_err());
    }
}

#[test]
fn affine_and_perfect_tables_match_closed_form_linear_and_quadratic_rate_returns() {
    for hardening in [0.0, 250.0] {
        for kinematic in [0.0, HK] {
            for exponent in [1.0, 2.0] {
                let material = UniaxialTabulatedOverstress::new_power_law(
                    table(&[(0.0, 3.0), (0.04, 3.0 + 0.04 * hardening)], kinematic),
                    ETA,
                    S0,
                    exponent,
                    DT,
                )
                .unwrap();
                let strain = 0.015;
                let trial_excess = E * strain - 3.0;
                let a = E + hardening + kinematic;
                let (delta, viscous_slope) = if exponent == 1.0 {
                    (trial_excess / (a + ETA / DT), ETA / DT)
                } else {
                    // x=sqrt(delta) solves a*x^2+b*x=f_trial. The stable
                    // quadratic formula is independent of the logarithmic solve.
                    let b = (ETA * S0 / DT).sqrt();
                    let x = 2.0 * trial_excess / (b + (b * b + 4.0 * a * trial_excess).sqrt());
                    (x * x, b / (2.0 * x))
                };
                let response = material
                    .response(strain, &material.initial_state())
                    .unwrap();
                near(
                    response.state.equivalent_plastic_strain(),
                    delta,
                    3e-16,
                    3e-13,
                );
                near(response.stress_pa, E * (strain - delta), 3e-12, 3e-13);
                let h = hardening + kinematic + viscous_slope;
                near(response.tangent_pa, E * h / (E + h), 3e-10, 3e-12);
                let full_linear = native(hardening, kinematic);
                let tensor = if exponent == 1.0 {
                    let law = LinearOverstressJ2::new(full_linear, ETA, DT).unwrap();
                    law.advance_with_energy(&[0.0; 6], &law.initial_state(), &response.strain, TOL)
                } else {
                    let law =
                        PowerLawOverstressJ2::new(full_linear, ETA, S0, exponent, DT).unwrap();
                    law.advance_with_energy(&[0.0; 6], &law.initial_state(), &response.strain, TOL)
                }
                .unwrap();
                let axial = material
                    .advance_with_energy(0.0, &material.initial_state(), strain, TOL)
                    .unwrap();
                near(axial.stress_pa[0], tensor.stress_pa[0], 3e-12, 3e-12);
                near(axial.state.alpha, tensor.state.alpha, 3e-16, 3e-12);
                near(
                    axial.energy.physical_dissipation_j_m3,
                    tensor.energy.physical_dissipation_j_m3,
                    3e-14,
                    3e-12,
                );
                near(
                    axial.energy.final_state.total_j_m3(),
                    tensor.energy.final_state.total_j_m3(),
                    3e-14,
                    3e-12,
                );
                near(
                    axial.energy.numerical_dissipation_j_m3,
                    tensor.energy.numerical_dissipation_j_m3,
                    3e-14,
                    3e-12,
                );
            }
        }
    }
}

#[test]
fn final_measured_endpoint_remains_bounded_through_native_energy_and_unloading() {
    let maximum = 0.012;
    for exponent in [1.0, 2.3] {
        let material = UniaxialTabulatedOverstress::new_power_law(
            table(&KNOTS[..4], HK),
            ETA,
            S0,
            exponent,
            DT,
        )
        .unwrap();
        let virgin = material.initial_state();
        let strain = maximum + (7.0 + HK * maximum + viscous(maximum, exponent, DT)) / E;
        let next = material.response(strain, &virgin).unwrap();
        assert!(next.state.equivalent_plastic_strain() <= maximum);
        near(next.state.equivalent_plastic_strain(), maximum, 3e-16, 0.0);
        let energy = material
            .advance_with_energy(0.0, &virgin, strain, TOL)
            .unwrap();
        assert!(energy.state.alpha <= maximum);
        material.validate_history(strain, &next.state).unwrap();
        material.stored_energy(strain, &next.state).unwrap();
        // Accepted positive overstress is valid, but another hold would require
        // new plastic history beyond the final measurement and must refuse.
        assert!(material.response(strain, &next.state).is_err());
        assert!(matches!(
            material.response(strain + 1e-5, &virgin),
            Err(TabulatedError::Domain { .. })
        ));
        let unload = next.state.plastic_strain() * (1.0 + HK / E);
        let elastic = material.response(unload, &next.state).unwrap();
        assert_eq!(elastic.state, next.state);
        assert_eq!(elastic.tangent_pa, E);
        material
            .advance_with_energy(strain, &next.state, unload, TOL)
            .unwrap();
        assert!(material.response(f64::NAN, &virgin).is_err());
        assert!(material.stress(f64::NAN, &virgin).is_nan());
        let poisoned = material.update_state(f64::NAN, &virgin);
        assert!(material.validate_history(0.0, &poisoned).is_err());
        assert!(material.response(0.1, &virgin).is_err());
        assert_eq!(virgin.equivalent_plastic_strain(), 0.0);
    }
}

#[test]
fn long_duration_limit_matches_rate_independent_axial_flow_without_changing_the_curve() {
    let backbone = table(&KNOTS, HK);
    let static_law = UniaxialTabulated::new(backbone.clone()).unwrap();
    let timed = UniaxialTabulatedOverstress::new_linear(backbone, ETA, DT).unwrap();
    let slow = timed.for_duration(1e12).unwrap();
    for strain in [0.001, 0.008, 0.014, -0.012] {
        let state = static_law.initial_state();
        let reference = static_law.response(strain, &state).unwrap();
        let response = slow.response(strain, &state).unwrap();
        near(response.stress_pa, reference.stress_pa, 3e-11, 3e-12);
        near(response.tangent_pa, reference.tangent_pa, 3e-10, 3e-12);
        near(
            response.state.equivalent_plastic_strain(),
            reference.state.equivalent_plastic_strain(),
            3e-14,
            3e-12,
        );
    }
}
