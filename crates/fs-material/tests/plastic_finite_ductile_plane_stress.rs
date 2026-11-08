//! G1/G3: progressive finite-sheet softening, complete transverse equilibrium,
//! exact condensed derivatives, rate relaxation and transactional refusals.
//! Synthetic coefficients establish numerical behavior, not an alloy fit.
#![deny(unsafe_code)]

use fs_material::ductile::evolution::{DuctileJ2Conditions, DuctileJ2Settings};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::plane_stress::FiniteDuctilePlaneStressStep;
use fs_material::plastic::finite::ductile::{FiniteDuctileError, FiniteDuctileJ2};
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};

const MU: f64 = 500.0;
const K: f64 = 1000.0;
const Y: f64 = 3.0;
const H: f64 = 40.0;
const WIDTH: f64 = 0.02;
const WORK: f64 = 0.004;

fn material(onset: f64, slope: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    FiniteDuctileJ2::new(
        FiniteJ2::new(
            FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: K,
                yield_stress_pa: Y,
                hardening_modulus_pa: H,
            },
            FiniteJ2Limits {
                minimum_jacobian: 0.1,
                maximum_elastic_log_strain: 0.7,
                maximum_plastic_log_strain: 2.0,
                maximum_accumulated_plastic_strain: 3.0,
            },
        )
        .unwrap(),
        TabulatedDuctileLocus::new(
            vec![(-50.0, onset - 50.0 * slope), (50.0, onset + 50.0 * slope)],
            DuctileDomain {
                triaxiality: [-50.0, 50.0],
                strain_rate_s_inv: [1e-12, 1000.0],
                temperature_k: [300.0, 300.0],
            },
        )
        .unwrap(),
        DuctileJ2Settings {
            characteristic_length_m: WIDTH,
            post_initiation_plastic_work_j_m2: WORK,
            maximum_damage_increment: 1.0,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-9,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap()
}

fn conditions(duration_s: f64) -> DuctileJ2Conditions {
    DuctileJ2Conditions {
        duration_s,
        temperature_k: 300.0,
    }
}

fn local() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-12,
        ..PlaneStressSettings::default()
    }
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}, difference={:.8e}",
        (actual - expected).abs()
    );
}

fn traction_free<L>(response: &FiniteDuctilePlaneStressStep<'_, L>) {
    assert!(response.thickness_stretch() > 0.0);
    assert!(response.relative_transverse_residual <= local().relative_tolerance);
    for component in [2, 5, 8] {
        near(
            response.step.first_piola_pa[component] / response.step.state.integrity(),
            0.0,
            2e-8,
            0.0,
        );
    }
    let e = response.step.energy;
    near(
        e.mechanical_work_j_m3,
        e.final_stored_j_m3 - e.initial_stored_j_m3
            + e.plastic_dissipation_j_m3
            + e.damage_release_j_m3
            + e.endpoint_work_defect_j_m3,
        1e-9,
        1e-12,
    );
    assert!(e.plastic_dissipation_j_m3 >= 0.0 && e.damage_release_j_m3 >= 0.0);
}

#[test]
fn independent_biaxial_softening_recovers_thickness_native_heat_and_negative_loading_tangent() {
    let onset = 0.01;
    let dt = 0.4;
    for viscosity in [0.0, 30.0] {
        let base = material(onset, 0.0);
        let law = if viscosity > 0.0 {
            base.with_overstress(FiniteOverstressParameters {
                viscosity_pa_s: viscosity,
                stress_scale_pa: Y,
                exponent: 1.0,
            })
            .unwrap()
        } else {
            base
        };
        let old = law.initial_state();
        for stretch in [1.03_f64, 1.08, 1.13] {
            let response = law
                .plane_stress_trial(
                    &[stretch, 0.0, 0.0, stretch],
                    &old,
                    [0.0, 0.0, 1.0],
                    conditions(dt),
                    local(),
                )
                .unwrap();
            // Independent scalar logarithmic-strain equations for biaxial
            // tension: tau=(q,q,0), ep=(alpha/2,alpha/2,-alpha).
            let compliance = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);
            let resistance_slope = H + viscosity / dt;
            let alpha = (stretch.ln() - Y * compliance) / (0.5 + resistance_slope * compliance);
            let q = Y + resistance_slope * alpha;
            let thickness = (-2.0 * stretch.ln() + 2.0 * q / (3.0 * K)).exp();
            let work = Y * (alpha - onset) + 0.5 * H * (alpha * alpha - onset * onset);
            let integrity = (-WIDTH * work / WORK).exp();
            let stored = q * q / (6.0 * MU) + 2.0 * q * q / (9.0 * K) + 0.5 * H * alpha * alpha;
            traction_free(&response);
            near(
                response.step.state.native_state().alpha,
                alpha,
                2e-10,
                2e-10,
            );
            near(response.thickness_stretch(), thickness, 2e-10, 2e-10);
            near(response.step.state.integrity(), integrity, 2e-10, 2e-9);
            near(response.piola_pa[0], integrity * q / stretch, 2e-8, 2e-9);
            near(response.piola_pa[3], integrity * q / stretch, 2e-8, 2e-9);
            near(
                response.step.energy.final_stored_j_m3,
                integrity * stored,
                2e-9,
                2e-9,
            );
            near(
                response.step.energy.plastic_dissipation_j_m3,
                integrity * (Y + viscosity * alpha / dt) * alpha,
                2e-9,
                2e-9,
            );
            near(
                response.step.energy.equivalent_plastic_work_increment_j_m2,
                WIDTH * work,
                2e-11,
                2e-9,
            );
            assert_eq!(response.step.energy.damage_release_j_m3, 0.0);
            let alpha_derivative = 1.0 / (stretch * (0.5 + resistance_slope * compliance));
            let derivative = integrity / stretch
                * (resistance_slope * alpha_derivative
                    - WIDTH / WORK * (Y + H * alpha) * alpha_derivative * q
                    - q / stretch);
            assert!(derivative < 0.0);
            near(
                response.tangent_pa[0][0] + response.tangent_pa[0][3],
                derivative,
                2e-6,
                2e-7,
            );
            assert_eq!(old, law.initial_state());
        }
    }
}

#[test]
fn intact_reduction_matches_native_return_and_elastic_unloading_keeps_irreversible_damage() {
    let intact = material(10.0, 0.0);
    let in_plane = [1.05, 0.07, -0.02, 0.99];
    let native = intact
        .native()
        .plane_stress_response(
            &in_plane,
            &intact.native().initial_state(),
            [0.0, 0.0, 1.0],
            local(),
        )
        .unwrap();
    let response = intact
        .plane_stress_trial(
            &in_plane,
            &intact.initial_state(),
            [0.0, 0.0, 1.0],
            conditions(0.2),
            local(),
        )
        .unwrap();
    traction_free(&response);
    assert_eq!(response.step.state.damage(), 0.0);
    near(
        response.step.state.native_state().alpha,
        native.material.state.alpha,
        2e-10,
        2e-10,
    );
    for i in 0..4 {
        near(response.piola_pa[i], native.piola_pa[i], 2e-8, 2e-9);
        for j in 0..4 {
            near(
                response.tangent_pa[i][j],
                native.tangent_pa[i][j],
                2e-6,
                2e-8,
            );
        }
    }
    let law = material(0.01, 0.0);
    let loaded = law
        .plane_stress_trial(
            &[1.08, 0.0, 0.0, 1.08],
            &law.initial_state(),
            [0.0, 0.0, 1.0],
            conditions(0.2),
            local(),
        )
        .unwrap();
    let retained = loaded.step.state.clone();
    let alpha = retained.native_state().alpha;
    let residual_stretch = (0.5 * alpha).exp();
    let unloaded = law
        .plane_stress_trial(
            &[residual_stretch, 0.0, 0.0, residual_stretch],
            &retained,
            loaded.transverse_column(),
            conditions(0.2),
            local(),
        )
        .unwrap();
    traction_free(&unloaded);
    near(unloaded.thickness_stretch(), (-alpha).exp(), 2e-10, 2e-10);
    for p in unloaded.step.first_piola_pa {
        near(p, 0.0, 2e-8, 0.0);
    }
    assert_eq!(unloaded.step.state.native_state(), retained.native_state());
    assert_eq!(unloaded.step.state.damage(), retained.damage());
    assert_eq!(unloaded.step.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(unloaded.step.energy.damage_release_j_m3, 0.0);
    let held = law
        .plane_stress_trial(
            &[residual_stretch, 0.0, 0.0, residual_stretch],
            &unloaded.step.state,
            unloaded.transverse_column(),
            conditions(0.2),
            local(),
        )
        .unwrap();
    assert_eq!(held.step.state, unloaded.step.state);
    assert_eq!(held.step.energy.mechanical_work_j_m3, 0.0);
    assert_eq!(retained, loaded.step.state);
}

#[test]
fn exact_schur_tangent_retains_moving_onset_transverse_shear_and_compressive_closure() {
    let law = material(0.018, -0.0001);
    let virgin = law.initial_state();
    // This actual original 3-D loading creates both transverse plastic shear
    // components. No metric or damage state is forged for the sheet reduction.
    let retained = law
        .trial(
            &[1.03, 0.06, 0.10, 0.0, 0.99, -0.07, 0.0, 0.0, 1.0],
            &virgin,
            conditions(0.2),
        )
        .unwrap()
        .state;
    let mut maximum_asymmetry = 0.0_f64;
    for (in_plane, old, shear_expected) in [
        ([1.04, 0.09, -0.025, 1.01], &virgin, false),
        ([1.06, 0.11, -0.03, 0.99], &retained, true),
        ([0.96, 0.04, -0.02, 0.98], &retained, true),
    ] {
        let saved = old.clone();
        let response = law
            .plane_stress_trial(&in_plane, old, [0.0, 0.0, 1.0], conditions(0.2), local())
            .unwrap();
        traction_free(&response);
        assert!(response.step.state.damage() > 0.0);
        if shear_expected {
            assert!(response.transverse_column()[0].abs() > 1e-4);
            assert!(response.transverse_column()[1].abs() > 1e-4);
        }
        let mut error = 0.0;
        let mut reference = 0.0;
        for j in 0..4 {
            let h = 2e-7;
            let mut plus = in_plane;
            let mut minus = in_plane;
            plus[j] += h;
            minus[j] -= h;
            let positive = law
                .plane_stress_trial(
                    &plus,
                    old,
                    response.transverse_column(),
                    conditions(0.2),
                    local(),
                )
                .unwrap();
            let negative = law
                .plane_stress_trial(
                    &minus,
                    old,
                    response.transverse_column(),
                    conditions(0.2),
                    local(),
                )
                .unwrap();
            for i in 0..4 {
                let numerical = (positive.piola_pa[i] - negative.piola_pa[i]) / (2.0 * h);
                near(response.tangent_pa[i][j], numerical, 3e-4, 2e-5);
                error += (response.tangent_pa[i][j] - numerical).powi(2);
                reference += numerical * numerical;
                maximum_asymmetry = maximum_asymmetry
                    .max((response.tangent_pa[i][j] - response.tangent_pa[j][i]).abs());
            }
        }
        assert!((error / reference).sqrt() < 3e-6);
        assert_eq!(old, &saved);
    }
    assert!(
        maximum_asymmetry > 1e-3,
        "the nonsymmetric derivative must not be symmetrized"
    );
}

#[test]
fn held_sheet_relaxes_once_per_interval_and_kinetic_condensation_uses_fixed_old_history() {
    let law = material(0.008, 0.0)
        .with_overstress(FiniteOverstressParameters {
            viscosity_pa_s: 30.0,
            stress_scale_pa: Y,
            exponent: 3.0,
        })
        .unwrap();
    let in_plane = [1.08, 0.025, -0.01, 1.03];
    let loaded = law
        .plane_stress_trial(
            &in_plane,
            &law.initial_state(),
            [0.0, 0.0, 1.0],
            conditions(0.03),
            local(),
        )
        .unwrap();
    let old = loaded.step.state.clone();
    let short = law
        .plane_stress_trial(
            &in_plane,
            &old,
            loaded.transverse_column(),
            conditions(0.02),
            local(),
        )
        .unwrap();
    let relaxed = law
        .plane_stress_trial(
            &in_plane,
            &old,
            loaded.transverse_column(),
            conditions(0.2),
            local(),
        )
        .unwrap();
    traction_free(&relaxed);
    assert!(short.step.state.native_state().alpha > old.native_state().alpha);
    assert!(relaxed.step.state.native_state().alpha > short.step.state.native_state().alpha);
    assert!(relaxed.step.state.damage() > short.step.state.damage());
    assert!(relaxed.piola_pa[0] < loaded.piola_pa[0]);
    assert!((relaxed.thickness_stretch() - loaded.thickness_stretch()).abs() > 1e-6);
    assert!(relaxed.step.energy.plastic_dissipation_j_m3 > 0.0);
    assert!(relaxed.step.energy.damage_release_j_m3 > 0.0);
    near(relaxed.step.energy.mechanical_work_j_m3, 0.0, 2e-9, 0.0);
    let repeated = law
        .plane_stress_trial(
            &in_plane,
            &old,
            loaded.transverse_column(),
            conditions(0.2),
            local(),
        )
        .unwrap();
    assert_eq!(repeated, relaxed);
    law.validate_committed(&relaxed.step.state, 1e-8).unwrap();
    for j in 0..4 {
        let h = 2e-7;
        let mut plus = in_plane;
        let mut minus = in_plane;
        plus[j] += h;
        minus[j] -= h;
        let positive = law
            .plane_stress_trial(
                &plus,
                &old,
                relaxed.transverse_column(),
                conditions(0.2),
                local(),
            )
            .unwrap();
        let negative = law
            .plane_stress_trial(
                &minus,
                &old,
                relaxed.transverse_column(),
                conditions(0.2),
                local(),
            )
            .unwrap();
        for i in 0..4 {
            near(
                relaxed.tangent_pa[i][j],
                (positive.piola_pa[i] - negative.piola_pa[i]) / (2.0 * h),
                3e-4,
                2e-5,
            );
        }
    }
    assert_eq!(old, loaded.step.state);
}

#[test]
fn invalid_geometry_owner_time_local_budget_and_damage_limit_refuse_without_mutating_history() {
    let law = material(0.01, 0.0);
    let old = law.initial_state();
    let saved = old.clone();
    for (in_plane, transverse, dt, settings) in [
        ([1.0, 0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 0.2, local()),
        ([1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.2, local()),
        ([1.04, 0.0, 0.0, 1.04], [0.0, 0.0, 1.0], 0.0, local()),
        (
            [1.04, 0.0, 0.0, 1.04],
            [0.0, 0.0, 1.0],
            0.2,
            PlaneStressSettings {
                max_iterations: 0,
                ..local()
            },
        ),
        (
            [1.04, 0.0, 0.0, 1.04],
            [0.0, 0.0, 1.0],
            0.2,
            PlaneStressSettings {
                max_iterations: 1,
                max_backtracks: 1,
                relative_tolerance: 1e-15,
            },
        ),
    ] {
        assert!(
            law.plane_stress_trial(&in_plane, &old, transverse, conditions(dt), settings)
                .is_err()
        );
        assert_eq!(old, saved);
    }
    let duplicate = law.clone();
    assert_eq!(
        duplicate.plane_stress_trial(
            &[1.04, 0.0, 0.0, 1.04],
            &old,
            [0.0, 0.0, 1.0],
            conditions(0.2),
            local(),
        ),
        Err(FiniteDuctileError::StateMismatch)
    );
    let limited = FiniteDuctileJ2::new(
        law.native().clone(),
        law.initiation_law().clone(),
        DuctileJ2Settings {
            maximum_damage_increment: 1e-8,
            ..law.settings()
        },
    )
    .unwrap();
    let limited_old = limited.initial_state();
    assert!(
        limited
            .plane_stress_trial(
                &[1.04, 0.0, 0.0, 1.04],
                &limited_old,
                [0.0, 0.0, 1.0],
                conditions(0.2),
                local(),
            )
            .is_err()
    );
    assert_eq!(limited_old, limited.initial_state());
}
