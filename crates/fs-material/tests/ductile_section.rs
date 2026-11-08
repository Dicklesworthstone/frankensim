//! Native progressive plane-stress sections; synthetic coefficients, no alloy fit.
#![deny(unsafe_code)]

use fs_material::IsotropicElastic;
use fs_material::ductile::evolution::section::DuctileSectionEnergyError;
use fs_material::ductile::evolution::{
    DuctileJ2, DuctileJ2Conditions, DuctileJ2Law, DuctileJ2Settings,
};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, ThicknessQuadrature,
};

const WIDTH: f64 = 0.1;
const GP: f64 = 0.00015;
const ONSET: f64 = 0.0002;

fn model(variable_locus: bool, damage_limit: f64) -> DuctileJ2<TabulatedDuctileLocus> {
    let values = if variable_locus {
        vec![(-2.0, 0.0004), (2.0, 0.0001)]
    } else {
        vec![(-2.0, ONSET), (2.0, ONSET)]
    };
    DuctileJ2::new(
        J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(), 1.0, 0.0).unwrap(),
        TabulatedDuctileLocus::new(
            values,
            DuctileDomain {
                triaxiality: [-2.0, 2.0],
                strain_rate_s_inv: [1e-10, 10.0],
                temperature_k: [300.0, 300.0],
            },
        )
        .unwrap(),
        DuctileJ2Settings {
            characteristic_length_m: WIDTH,
            post_initiation_plastic_work_j_m2: GP,
            maximum_damage_increment: damage_limit,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-10,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap()
}

fn section(
    model: &DuctileJ2<TabulatedDuctileLocus>,
    durations: [f64; 2],
) -> LayeredPlaneStressSection<DuctileJ2Law<'_, TabulatedDuctileLocus>> {
    LayeredPlaneStressSection::new(
        [(-0.05, 0.0), (0.0, 0.05)]
            .into_iter()
            .enumerate()
            .map(|(i, (bottom_z_m, top_z_m))| PlaneStressLayer {
                bottom_z_m,
                top_z_m,
                subdivisions: 1,
                law: model
                    .small_strain_law(DuctileJ2Conditions {
                        duration_s: durations[i],
                        temperature_k: 300.0,
                    })
                    .unwrap(),
            })
            .collect(),
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            max_iterations: 64,
            stress_tolerance_pa: 1e-11,
            strain_limit: 0.05,
        },
        SectionBudget {
            maximum_layers: 2,
            maximum_quadrature_points: 6,
        },
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}

#[test]
fn bending_retains_each_point_and_closes_separate_damage_and_plastic_accounts() {
    let model = model(false, 0.95);
    let section = section(&model, [1.0; 2]);
    let virgin = section.initial_state();
    let frozen = virgin.clone();
    let elastic = section
        .evaluate([0.0, 0.0, 0.0, 0.0, 0.0, 0.02], &virgin)
        .unwrap();
    let target = [0.0, 0.0, 0.0, 0.0, 0.0, 0.12];
    let step = section.evaluate(target, &elastic.state).unwrap();
    let account = section
        .audit_ductile_energy(&elastic.state, &step.state, 1e-11)
        .unwrap();
    assert_eq!(virgin, frozen);
    assert_eq!(account.resultants, step.resultants);
    for force in &step.resultants[..3] {
        near(*force, 0.0, 1e-14);
    }
    assert!(step.resultants[5] > 0.0);
    let mut damaged = 0;
    let mut equivalent_work = 0.0;
    let mut undamaged_moment = 0.0;
    for (point, value) in section.points().iter().zip(&step.state.points) {
        let shear = 0.5 * point.z_m * target[5];
        // Independent perfect-plastic pure-shear radial return: mu=500 Pa.
        let alpha = ((1000.0 * 3.0_f64.sqrt() * shear.abs() - 1.0) / 1500.0).max(0.0);
        let expected_damage = 1.0 - (-(WIDTH / GP) * (alpha - ONSET).max(0.0)).exp();
        near(value.history.native_state().alpha, alpha, 2e-14);
        near(value.history.damage(), expected_damage, 2e-12);
        near(value.strain[2], 0.0, 1e-14);
        near(value.strain[3], shear, 1e-14);
        damaged += usize::from(value.history.damage() > 0.0);
        equivalent_work += point.weight_m * GP * expected_damage / WIDTH;
        let native_shear = (1000.0 * shear.abs()).min(1.0 / 3.0_f64.sqrt()) * shear.signum();
        undamaged_moment += point.weight_m * point.z_m * native_shear;
    }
    assert_eq!(damaged, 4, "the two near-midplane points remain elastic");
    assert!(step.resultants[5] < 0.7 * undamaged_moment);
    assert!(account.plastic_dissipation_j_m2 > 0.0);
    assert!(account.damage_release_j_m2 > 0.0);
    assert!(account.numerical_dissipation_j_m2 > 0.0);
    near(
        account.degraded_equivalent_plastic_work_j_m2,
        equivalent_work,
        1e-15,
    );
    near(
        account.generalized_work_j_m2,
        step.resultants[5] * (target[5] - elastic.state.generalized_strain[5]),
        1e-15,
    );
    near(account.balance_residual_j_m2, 0.0, 1e-11);
    near(
        account.material_work_j_m2,
        account.stored_change_j_m2
            + account.plastic_dissipation_j_m2
            + account.damage_release_j_m2
            + account.numerical_dissipation_j_m2,
        1e-11,
    );
}

#[test]
fn unload_and_reverse_keep_damage_irreversible_without_a_second_return() {
    let model = model(false, 0.95);
    let section = section(&model, [1.0; 2]);
    let loaded = section
        .evaluate([0.0, 0.0, 0.0, 0.0, 0.0, 0.12], &section.initial_state())
        .unwrap();
    let unloaded = section
        .evaluate([0.0, 0.0, 0.0, 0.0, 0.0, 0.119], &loaded.state)
        .unwrap();
    let account = section
        .audit_ductile_energy(&loaded.state, &unloaded.state, 1e-11)
        .unwrap();
    assert_eq!(account.damage_release_j_m2, 0.0);
    assert_eq!(account.plastic_dissipation_j_m2, 0.0);
    assert_eq!(account.degraded_equivalent_plastic_work_j_m2, 0.0);
    for (before, after) in loaded.state.points.iter().zip(&unloaded.state.points) {
        assert_eq!(before.history.damage(), after.history.damage());
        assert_eq!(before.history.native_state(), after.history.native_state());
    }
    let reversed = section
        .evaluate([0.0, 0.0, 0.0, 0.0, 0.0, -0.12], &unloaded.state)
        .unwrap();
    let reverse_work = section
        .audit_ductile_energy(&unloaded.state, &reversed.state, 1e-11)
        .unwrap();
    assert!(reverse_work.damage_release_j_m2 > 0.0);
    assert!(reverse_work.plastic_dissipation_j_m2 > 0.0);
    for (before, after) in unloaded.state.points.iter().zip(&reversed.state.points) {
        assert!(after.history.damage() >= before.history.damage());
        assert!(after.history.native_state().alpha >= before.history.native_state().alpha);
    }
    let held = section
        .evaluate(reversed.state.generalized_strain, &reversed.state)
        .unwrap();
    assert_eq!(held.state, reversed.state);
    let hold = section
        .audit_ductile_energy(&reversed.state, &held.state, 1e-11)
        .unwrap();
    assert_eq!(hold.material_work_j_m2, 0.0);
    assert_eq!(hold.damage_release_j_m2, 0.0);
    assert_eq!(hold.numerical_dissipation_j_m2, 0.0);
}

#[test]
fn full_damaged_section_tangent_matches_frozen_history_finite_differences() {
    let model = model(true, 0.95);
    let section = section(&model, [1.0; 2]);
    let old = section.initial_state();
    // Every layer is in volumetric tension, away from the unilateral branch,
    // yield, onset and locus knots. Mixed membrane/bending/shear probes every
    // Schur-complement column and engineering-shear factor.
    let target = [0.003, 0.001, 0.0006, 0.01, 0.003, 0.02];
    let response = section.evaluate(target, &old).unwrap();
    assert!(
        response
            .state
            .points
            .iter()
            .all(|p| p.history.damage() > 0.0)
    );
    assert!(!section.tangent_symmetric());
    let asymmetry = (0..6)
        .flat_map(|i| (0..6).map(move |j| (i, j)))
        .map(|(i, j)| (response.tangent[i][j] - response.tangent[j][i]).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        asymmetry > 1e-6,
        "the onset chain rule remains nonsymmetric"
    );
    for column in 0..6 {
        let delta = if column < 3 { 1e-8 } else { 1e-7 };
        let mut plus = target;
        let mut minus = target;
        plus[column] += delta;
        minus[column] -= delta;
        let p = section.evaluate(plus, &old).unwrap();
        let m = section.evaluate(minus, &old).unwrap();
        for row in 0..6 {
            let finite_difference = (p.resultants[row] - m.resultants[row]) / (2.0 * delta);
            near(
                response.tangent[row][column],
                finite_difference,
                3e-5 * finite_difference.abs().max(1e-3),
            );
        }
    }
    assert_eq!(old, section.initial_state());
}

#[test]
fn owner_duration_and_failed_trial_are_transactional() {
    let model = model(false, 0.95);
    let material = section(&model, [1.0; 2]);
    let old = material.initial_state();
    let target = [0.0, 0.0, 0.0, 0.0, 0.0, 0.12];
    let next = material.evaluate(target, &old).unwrap();
    let saved = next.state.clone();
    let equivalent_clone = model.clone();
    let foreign = section(&equivalent_clone, [1.0; 2]);
    let other = foreign.evaluate(target, &foreign.initial_state()).unwrap();
    assert!(matches!(
        material.audit_ductile_energy(&old, &other.state, 1e-11),
        Err(DuctileSectionEnergyError::State { .. })
    ));
    let mixed_time = section(&model, [1.0, 0.5]);
    assert!(
        mixed_time
            .audit_ductile_energy(&old, &next.state, 1e-11)
            .is_err()
    );
    let mut forged = next.state.clone();
    forged.points[0].strain[2] += 1e-5;
    assert!(material.audit_ductile_energy(&old, &forged, 1e-11).is_err());
    assert!(
        material
            .audit_ductile_energy(&old, &next.state, 0.0)
            .is_err()
    );
    let restricted = self::model(false, 0.01);
    let limited = section(&restricted, [1.0; 2]);
    let untouched = limited.initial_state();
    assert!(limited.evaluate(target, &untouched).is_err());
    assert_eq!(untouched, limited.initial_state());
    assert_eq!(next.state, saved);
    assert_eq!(old, material.initial_state());
}
