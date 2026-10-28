//! G1/G3 independent SI references. The declared Hertz coefficients and
//! caloric sources are synthetic implementation fixtures, not measurements.
use super::common;
use fs_contact::normal_patch::*;
use fs_couple::tribothermal::normal::*;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{HeatPartition, InputAuthority, InterfaceMedium, InterfaceSystemRef};
use std::collections::BTreeMap;
use std::convert::Infallible;

const CELLS: [CellId; 2] = [CellId(30), CellId(10)];

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-7, "{a} != {b}");
}
fn request(indentation: f64, damping: f64) -> NormalPatchRequest {
    NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: "declared normal free energy".into(),
            source_id: "synthetic Hertz reference".into(),
            state_id: "retained normal history".into(),
        },
        interface: InterfaceSystemRef::new(
            "ordered A/B",
            "normal history",
            "authored",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        law: NormalPatchLaw::HuntCrossleySphere {
            effective_radius_m: 1.0,
            reduced_modulus_pa: 1e6,
            dissipation_s_per_m: damping,
        },
        geometry: NormalPatchGeometry::SpherePlane,
        indentation_m: indentation,
        indentation_rate_m_per_s: 0.0,
        step_s: 0.0,
        line_load_n_per_m: 0.0,
        applicability: ApplicabilityInput {
            half_space_depth_m: 10.0,
            layer_thickness_m: 10.0,
            yield_strength_pa: 1e12,
            characteristic_rate_m_per_s: 1.0,
            temperature_k: 300.0,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 0.2,
            max_strain: 0.05,
            max_patch_to_depth: 0.1,
            max_patch_to_layer: 0.1,
            max_pressure_to_yield: 0.1,
            max_rate_ratio: 1.0,
            min_temperature_k: 200.0,
            max_temperature_k: 600.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.0,
            modulus_relative: 0.0,
            load_relative: 0.0,
        },
    }
}
fn model(indentation: f64, damping: f64, slopes: [f64; 2]) -> AffineNormalModulus {
    AffineNormalModulus::new(
        request(indentation, damping),
        [300.0; 2],
        slopes,
        [[200.0, 600.0]; 2],
    )
    .unwrap()
}
fn bind<'a>(
    thermal: &DepositionSimulation<'_>,
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    law: AffineNormalModulus,
) -> CaloricNormalContact<'a> {
    CaloricNormalContact::new(
        thermal,
        law,
        [curves[0].curve(), curves[1].curve()],
        CELLS,
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
    )
    .unwrap()
}
fn controls() -> NormalThermalSolve {
    NormalThermalSolve {
        maximum_cells: 2,
        maximum_contact_trials: 800,
        maximum_newton_iterations: 40,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 2e-8,
        thermal_energy_tolerance_j: 4e-8,
        total_energy_tolerance_j: 1e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn other(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn guess(thermal: &DepositionSimulation<'_>) -> Vec<f64> {
    vec![0.0; thermal.cells().len()]
}
fn basis(d: f64) -> f64 {
    (8.0 / 15.0) * d.powf(2.5)
}

#[test]
fn loaded_elastic_contact_cools_both_finite_bodies_and_changes_patch() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let thermal = common::initial(&curves, 0.0, None);
    let contact = bind(&thermal, &curves, model(0.002, 0.0, [-1000.0, -500.0]));
    let candidate = contact
        .step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let dg = basis(0.01) - basis(0.002);
    // Independent entropy relation Ci*(Ti1-Ti0)=ai*Ti1*dG.
    let expected = [
        300.0 * 10.0 / (10.0 + 1000.0 * dg),
        300.0 * 20.0 / (20.0 + 500.0 * dg),
    ];
    let modulus = 1e6 - 1000.0 * (expected[0] - 300.0) - 500.0 * (expected[1] - 300.0);
    for i in 0..2 {
        near(candidate.contact().temperature_k()[i], expected[i]);
    }
    let report = candidate.report();
    near(report.work.elastic_work_j, modulus * dg);
    near(report.work.internal_energy_change_j, 1_450_000.0 * dg);
    near(
        report.endpoint.normal_force_n,
        (4.0 / 3.0) * modulus * 0.01_f64.powf(1.5),
    );
    near(report.endpoint.patch_radius_m, 0.1);
    near(report.work.mean_normal_force_n, modulus * dg / 0.008);
    assert!(report.work.mean_normal_force_n < report.endpoint.normal_force_n);
    assert_eq!(report.work.physical_dissipation_j, 0.0);
    assert_eq!(report.physical_heat_j, [0.0; 2]);
    assert!(report.reversible_heat_j.iter().all(|q| *q < 0.0));
    near(report.total_energy_residual_j, 0.0);
    near(contact.endpoint().approach_m, 0.002);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn heating_held_loaded_contact_changes_force_without_fake_heat() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let thermal = common::initial(&curves, 0.0, None);
    let contact = bind(&thermal, &curves, model(0.01, 2.0, [-1000.0, -500.0]));
    let initial_u = contact.thermodynamics().internal_energy_j;
    let candidate = contact
        .step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(CELLS[0], 100.0), (CELLS[1], 400.0)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || false,
        )
        .unwrap();
    let report = candidate.report();
    near(candidate.contact().temperature_k()[0], 310.0);
    near(candidate.contact().temperature_k()[1], 320.0);
    near(report.modulus.modulus_pa, 980000.0);
    near(report.thermodynamics.internal_energy_j, initial_u);
    near(report.work.work_into_contact_j, 0.0);
    near(
        report.work.coefficient_free_energy_change_j,
        -20000.0 * basis(0.01),
    );
    for q in report.deposited_heat_j {
        near(q, 0.0);
    }
    near(report.external_heat_j, 500.0);
    near(report.total_energy_residual_j, 0.0);
}

#[test]
fn hunt_crossley_heat_is_segment_integral_on_loading_and_release() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let mut thermal = common::initial(&curves, 0.0, None);
    let mut contact = bind(&thermal, &curves, model(0.0, 10.0, [0.0; 2]));
    let loaded = contact
        .step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let elastic = 1e6 * basis(0.01);
    let dissipated = 0.1 * elastic;
    near(loaded.report().work.physical_dissipation_j, dissipated);
    near(loaded.report().work.work_into_contact_j, 1.1 * elastic);
    // Native point receipt is an instantaneous sample, not endpoint quadrature.
    assert_eq!(loaded.report().endpoint.irreversible_work_j, 0.0);
    near(
        loaded.report().endpoint.dissipated_power_w,
        10.0 * (4.0 / 3.0) * 1e6 * 0.01_f64.powf(1.5) * 0.0001,
    );
    assert!(loaded.report().endpoint.dissipated_power_w > dissipated);
    near(loaded.report().physical_heat_j[0], 0.4 * dissipated);
    near(loaded.report().physical_heat_j[1], 0.6 * dissipated);
    (thermal, contact, _) = loaded.into_parts();
    let released = contact
        .step_candidate(
            &thermal,
            2.0,
            0.0,
            controls(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    near(released.report().work.physical_dissipation_j, dissipated);
    near(released.report().work.work_into_contact_j, -0.9 * elastic);
    near(released.contact().thermodynamics().internal_energy_j, 0.0);
    near(released.report().endpoint.normal_force_n, 0.0);
    near(
        released.contact().temperature_k()[0],
        300.0 + 0.8 * dissipated / 10.0,
    );
    near(
        released.contact().temperature_k()[1],
        300.0 + 1.2 * dissipated / 20.0,
    );
    near(released.report().total_energy_residual_j, 0.0);
}

#[test]
fn elliptic_circular_limit_and_tiny_increment_keep_native_force() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let thermal = common::initial(&curves, 0.0, None);
    let sphere = bind(&thermal, &curves, model(0.01, 0.0, [0.0; 2]));
    let mut input = request(0.01, 0.0);
    input.geometry = NormalPatchGeometry::EllipticParaboloid;
    input.law = NormalPatchLaw::HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: 1.0,
        minimum_principal_curvature_m_inverse: 1.0,
        reduced_modulus_pa: 1e6,
    };
    let ellipse = bind(
        &thermal,
        &curves,
        AffineNormalModulus::new(input, [300.0; 2], [0.0; 2], [[200.0, 600.0]; 2]).unwrap(),
    );
    let next_indent = 0.01 + 1e-12;
    let a = sphere
        .step_candidate(
            &thermal,
            1.0,
            next_indent,
            controls(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let b = ellipse
        .step_candidate(
            &thermal,
            1.0,
            next_indent,
            controls(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    near(
        a.report().work.mean_normal_force_n,
        (4.0 / 3.0) * 1e6 * 0.01_f64.powf(1.5),
    );
    near(
        a.report().work.mean_normal_force_n,
        b.report().work.mean_normal_force_n,
    );
    near(
        a.report().work.elastic_work_j,
        b.report().work.elastic_work_j,
    );
    assert!(a.report().work.elastic_work_j > 0.0);
}

#[test]
fn temperature_feedback_converges_to_independent_adiabatic_entropy_solution() {
    let curves = [
        common::material(1.0, 0.03, 0.0).1,
        common::material(1.0, 0.04, 0.0).1,
    ];
    let dg = basis(0.01) - basis(0.002);
    let exact = [
        300.0 * (-1000.0 * dg / 0.03).exp(),
        300.0 * (-500.0 * dg / 0.04).exp(),
    ];
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let mut thermal = common::initial(&curves, 0.0, None);
        let mut contact = bind(&thermal, &curves, model(0.002, 0.0, [-1000.0, -500.0]));
        for i in 1..=n {
            let t = i as f64 / n as f64;
            let candidate = contact
                .step_candidate(
                    &thermal,
                    t,
                    0.002 + 0.008 * t,
                    controls(),
                    &guess(&thermal),
                    other,
                    || false,
                )
                .unwrap();
            near(candidate.report().total_energy_residual_j, 0.0);
            (thermal, contact, _) = candidate.into_parts();
        }
        errors.push(
            (contact.temperature_k()[0] - exact[0]).abs()
                + (contact.temperature_k()[1] - exact[1]).abs(),
        );
    }
    assert!(
        errors[0] / errors[1] > 1.7 && errors[0] / errors[1] < 2.3,
        "{errors:?}"
    );
    assert!(
        errors[1] / errors[2] > 1.7 && errors[1] / errors[2] < 2.3,
        "{errors:?}"
    );
}

#[test]
fn source_refusal_cancellation_and_native_domain_failure_are_atomic() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let thermal = common::initial(&curves, 0.0, None);
    let contact = bind(&thermal, &curves, model(0.002, 0.0, [-1000.0, -500.0]));
    let run = || {
        contact
            .step_candidate(
                &thermal,
                1.0,
                0.01,
                controls(),
                &guess(&thermal),
                other,
                || false,
            )
            .unwrap()
    };
    let first = run();
    let mut budget = controls();
    budget.maximum_contact_trials = 1;
    assert!(matches!(
        contact.step_candidate(&thermal, 1.0, 0.01, budget, &guess(&thermal), other, || {
            false
        }),
        Err(NormalThermalError::TrialBudget)
    ));
    assert!(matches!(
        contact.step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            other,
            || true
        ),
        Err(NormalThermalError::Cancelled)
    ));
    let late = contact.step_candidate(
        &thermal,
        1.0,
        0.01,
        controls(),
        &guess(&thermal),
        |state| {
            if state.temperature_k(CELLS[0]).unwrap() < 299.99 {
                return Err("original source domain");
            }
            Ok(StepCoefficients {
                sources_w: BTreeMap::new(),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            })
        },
        || false,
    );
    assert!(matches!(
        late,
        Err(NormalThermalError::Source("original source domain"))
    ));
    assert!(matches!(
        contact.step_candidate(
            &thermal,
            1.0,
            0.1,
            controls(),
            &guess(&thermal),
            other,
            || false
        ),
        Err(NormalThermalError::Normal(_))
    ));
    let replay = run();
    assert_eq!(
        first.contact().temperature_k(),
        replay.contact().temperature_k()
    );
    assert_eq!(first.report().work, replay.report().work);
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    near(contact.endpoint().approach_m, 0.002);
}

#[derive(Clone)]
struct UnstablePotential(NormalPatchRequest);
impl NormalModulus for UnstablePotential {
    fn reference_request(&self) -> &NormalPatchRequest {
        &self.0
    }
    fn evaluate(&self, t: [f64; 2]) -> Result<NormalModulusState, NormalThermalError> {
        let dt = t[0] - 300.0;
        Ok(NormalModulusState {
            modulus_pa: 1e6 + 0.5e6 * dt * dt,
            gradient_pa_k: [1e6 * dt, 0.0],
            hessian_pa_k2: [[1e6, 0.0], [0.0, 0.0]],
        })
    }
}
#[test]
fn refuses_nonpositive_total_heat_capacity_and_nonoriginal_caloric_object() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let other_curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let thermal = common::initial(&curves, 0.0, None);
    assert!(matches!(
        CaloricNormalContact::new(
            &thermal,
            UnstablePotential(request(0.01, 0.0)),
            [curves[0].curve(), curves[1].curve()],
            CELLS,
            HeatPartition::new(0.4, 0.6, 0.0).unwrap()
        ),
        Err(NormalThermalError::Invalid(
            "nonpositive total normal-contact thermal capacity"
        ))
    ));
    assert!(
        CaloricNormalContact::new(
            &thermal,
            model(0.002, 0.0, [0.0; 2]),
            [other_curves[0].curve(), other_curves[1].curve()],
            CELLS,
            HeatPartition::new(0.4, 0.6, 0.0).unwrap()
        )
        .is_err()
    );
    assert!(
        AffineNormalModulus::new(
            request(0.002, 0.0),
            [300.0; 2],
            [-10000.0; 2],
            [[200.0, 600.0]; 2]
        )
        .is_err()
    );
}

#[path = "normal/material.rs"]
mod material;
