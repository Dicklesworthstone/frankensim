//! G0/G1/G3: original string motion -> source-bound tensile admission.
//! Manufactured elastic data, not measured wire fracture validation.

use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::supported::SupportedMaterialStrings;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
use fs_couple::render::plate::impact::supported::TranslatingSupport;
use fs_couple::string_specimen::with_uniform_circular_material_state;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, DynViscosity, Pressure};

const LIMIT: &str = "declared_tensile_limit";
const RADIUS: f64 = 0.0003;
const LENGTH: f64 = 0.5;
const YOUNG: f64 = 200e9;

#[path = "tensile/supported.rs"]
mod supported;

fn strength_material(
    young: f64,
    limit: f64,
    frequency_restriction: bool,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, young),
        (ETA, DynViscosity::DIMS, 10.0),
        (LIMIT, Pressure::DIMS, limit),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        if name == ETA || (name == LIMIT && frequency_restriction) {
            validity = validity.with("omega", 0.0, 200_000.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured elastic tensile string".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "tensile test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn tensile_source(limit: f64) -> MaterialStringSource {
    MaterialStringSource::Homogeneous(
        with_uniform_circular_material_state(
            template(3),
            RADIUS,
            &strength_material(YOUNG, limit, false),
        )
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap(),
    )
}

fn coordinates(amplitude: f64) -> Vec<f64> {
    initial(3, amplitude)
        .iter()
        .map(|state| state.displacement_m_sqrt_kg)
        .collect()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-6),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_tensile_envelope_uses_actual_stretching_and_outer_fibre_bending() {
    with_cx(|cx, gate| {
        let source = tensile_source(1e9);
        let MaterialStringSource::Homogeneous(specimen) = &source else {
            unreachable!()
        };
        let original = specimen.material().clone();
        let image = source
            .into_nonlinear(cx, gas(), 48_000, 0.2, initial(3, 0.0), gate)
            .unwrap()
            .with_tensile_limit(LIMIT)
            .unwrap();
        let amplitude = 0.001;
        let assessment = image
            .tensile_assessment(&coordinates(amplitude))
            .unwrap()
            .unwrap();
        let wave = core::f64::consts::PI / LENGTH;
        let strain = 0.25 * (wave * amplitude).powi(2);
        let tension = 20.0 + image.axial_rigidity_n() * strain;
        let axial_stress = tension / (core::f64::consts::PI * RADIUS * RADIUS);
        let curvature = wave * wave * amplitude;
        let stress = axial_stress + YOUNG * RADIUS * curvature;
        near(assessment.tension_n, tension);
        near(assessment.curvature_bound_per_m, curvature);
        near(assessment.governing_tensile_stress_bound_pa, stress);
        near(assessment.utilization, stress / 1e9);
        assert_eq!(
            image.tensile_capacity().unwrap().section().layers()[0].material,
            original
        );

        // Axial force alone passes; actual bending makes the selected source
        // limit inadmissible. The check must not reduce to T / force capacity.
        let limit = 0.5 * (axial_stress + stress);
        let unguarded = tensile_source(limit)
            .into_nonlinear(cx, gas(), 48_000, 0.2, initial(3, amplitude), gate)
            .unwrap();
        assert!(unguarded.with_tensile_limit(LIMIT).is_err());
    });
}

#[test]
fn g3_tensile_envelope_contains_the_whole_span_and_is_transverse_frame_invariant() {
    with_cx(|cx, gate| {
        let source = tensile_source(1e9)
            .with_transverse_polarizations(
                TransverseDirections::try_new([1.0, 0.0], [0.0, 1.0]).unwrap(),
            )
            .unwrap();
        let image = source
            .into_nonlinear(
                cx,
                gas(),
                48_000,
                0.2,
                vec![ModalAcousticState::default(); 6],
                gate,
            )
            .unwrap()
            .with_tensile_limit(LIMIT)
            .unwrap();
        let x = [0.002, -0.0005, 0.0003];
        let y = [0.0015, 0.0002, -0.0001];
        let norm = (0.5 * image.card().linear_density_kg_m * LENGTH).sqrt();
        let q: Vec<_> = x.into_iter().chain(y).map(|value| value * norm).collect();
        let mut rotated = q.clone();
        for i in 0..3 {
            rotated[i] = 0.6 * q[i] - 0.8 * q[3 + i];
            rotated[3 + i] = 0.8 * q[i] + 0.6 * q[3 + i];
        }
        let original = image.tensile_assessment(&q).unwrap().unwrap();
        let rotation = image.tensile_assessment(&rotated).unwrap().unwrap();
        near(
            rotation.curvature_bound_per_m,
            original.curvature_bound_per_m,
        );
        near(rotation.tension_n, original.tension_n);
        near(rotation.utilization, original.utilization);
        for station in 0..=128 {
            let mut curvature = [0.0; 2];
            for i in 0..3 {
                let wave = (i + 1) as f64 * core::f64::consts::PI / LENGTH;
                let shape = -wave * wave * (wave * LENGTH * station as f64 / 128.0).sin();
                curvature[0] += shape * x[i];
                curvature[1] += shape * y[i];
            }
            assert!(
                curvature[0].hypot(curvature[1]) <= original.curvature_bound_per_m * (1.0 + 1e-14)
            );
        }
    });
}

#[test]
fn g1_bonded_tensile_layers_keep_their_own_stiffness_geometry_and_limits() {
    with_cx(|cx, gate| {
        let section = ConcentricElasticSection::try_new(
            0.0,
            vec![
                ConcentricMaterialLayer {
                    outer_radius_m: 0.0002,
                    material: strength_material(200e9, 800e6, false),
                },
                ConcentricMaterialLayer {
                    outer_radius_m: RADIUS,
                    material: strength_material(50e9, 120e6, false),
                },
            ],
        )
        .unwrap();
        let expected = section.clone();
        let image = MaterialStringSource::Concentric {
            template: template(3),
            section: section.with_kelvin_voigt_bending(ETA).unwrap(),
            prestress: StringPrestress::FixedTension(20.0),
        }
        .into_nonlinear(cx, gas(), 48_000, 0.2, initial(3, 0.0), gate)
        .unwrap()
        .with_tensile_limit(LIMIT)
        .unwrap();
        let assessment = image
            .tensile_assessment(&coordinates(0.001))
            .unwrap()
            .unwrap();
        let strain = assessment.tension_n / expected.axial_stiffness_n();
        let expected_stress = 50e9 * (strain + RADIUS * assessment.curvature_bound_per_m);
        assert_eq!(assessment.governing_layer, 1);
        near(
            assessment.governing_tensile_stress_bound_pa,
            expected_stress,
        );
        near(assessment.utilization, expected_stress / 120e6);
        assert_eq!(image.tensile_capacity().unwrap().section(), &expected);
    });
}

#[test]
fn g0_tensile_crossing_rolls_back_reference_and_prepared_native_state() {
    with_cx(|cx, gate| {
        let static_stress = 20.0 / (core::f64::consts::PI * RADIUS * RADIUS);
        let mut state = initial(3, 0.0);
        state[0].velocity_m_sqrt_kg_per_s = 0.02;
        let image = tensile_source(static_stress + 1e4)
            .into_nonlinear(cx, gas(), 48_000, 0.2, state, gate)
            .unwrap()
            .with_tensile_limit(LIMIT)
            .unwrap();
        let mut native =
            ImpactSystem::new(image.bodies().to_vec(), vec![], vec![], vec![], config()).unwrap();
        let before = native.state().to_vec();
        let energy = native.stored_energy_j();
        let error = native.step(&[0.0; 3], gate).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(native.state(), before);
        assert_eq!(native.stored_energy_j(), energy);
        assert_eq!(native.samples(), 0);

        let mut prepared = image
            .prepare(vec![], vec![], vec![], vec![], config(), gate)
            .unwrap();
        let assessment = prepared.tensile_assessment().unwrap().unwrap();
        assert!(!assessment.exceeds_limit());
        let error = prepared.step(&[0.0; 3], gate).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(prepared.system().state(), before);
        assert_eq!(prepared.system().stored_energy_j(), energy);
        assert_eq!(prepared.system().samples(), 0);
        assert_eq!(prepared.tensile_assessment().unwrap().unwrap(), assessment);
    });
}

#[test]
fn g0_supported_conversion_retains_tensile_admission_and_atomic_refusal() {
    with_cx(|cx, gate| {
        let static_stress = 20.0 / (core::f64::consts::PI * RADIUS * RADIUS);
        let mut state = initial(3, 0.0);
        state[0].velocity_m_sqrt_kg_per_s = 0.02;
        let image = tensile_source(static_stress + 1e4)
            .into_nonlinear(cx, gas(), 48_000, 0.2, state, gate)
            .unwrap()
            .with_tensile_limit(LIMIT)
            .unwrap();
        let support = TranslatingSupport {
            mass_kg: 0.05,
            stiffness_n_m: 1000.0,
            damping_n_s_m: 0.0,
            initial_position_m: 0.0,
            initial_velocity_m_s: 0.0,
            maximum_travel_m: 0.02,
            maximum_slope: 0.2,
        };
        let supported =
            SupportedMaterialStrings::in_comoving_gas(vec![image], support, gate).unwrap();
        let mut run = supported.prepare(config(), gate).unwrap();
        let before = run.system().state().to_vec();
        let energy = run.system().stored_energy_j();
        let error = run.step(&[0.0; 4], gate).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(run.system().state(), before);
        assert_eq!(run.system().stored_energy_j(), energy);
        assert_eq!(run.system().samples(), 0);
    });
}

#[test]
fn g0_tensile_limits_are_explicit_static_and_do_not_change_admitted_mechanics() {
    with_cx(|cx, gate| {
        let missing = source(7800.0, YOUNG, 3)
            .into_nonlinear(cx, gas(), 48_000, 0.2, initial(3, 0.0), gate)
            .unwrap();
        assert!(missing.tensile_capacity().is_none());
        assert!(missing.with_tensile_limit(LIMIT).is_err());
        let restricted = MaterialStringSource::Homogeneous(
            with_uniform_circular_material_state(
                template(3),
                RADIUS,
                &strength_material(YOUNG, 1e9, true),
            )
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .unwrap(),
        );
        assert!(
            restricted
                .into_nonlinear(cx, gas(), 48_000, 0.2, initial(3, 0.0), gate)
                .unwrap()
                .with_tensile_limit(LIMIT)
                .is_err()
        );
        let original = tensile_source(1e9)
            .into_nonlinear(cx, gas(), 48_000, 0.2, initial(3, 0.001), gate)
            .unwrap();
        let guarded = original.clone().with_tensile_limit(LIMIT).unwrap();
        assert!(guarded.clone().with_tensile_limit(LIMIT).is_err());
        assert!(guarded.tensile_assessment(&[0.0; 2]).is_err());
        assert!(guarded.tensile_assessment(&[f64::NAN; 3]).is_err());
        let mut a = original
            .prepare(vec![], vec![], vec![], vec![], config(), gate)
            .unwrap();
        let mut b = guarded
            .prepare(vec![], vec![], vec![], vec![], config(), gate)
            .unwrap();
        for _ in 0..16 {
            assert_eq!(
                a.step(&[0.0; 3], gate).unwrap().stored_energy_j,
                b.step(&[0.0; 3], gate).unwrap().stored_energy_j
            );
            assert_eq!(a.system().state(), b.system().state());
        }
    });
}
