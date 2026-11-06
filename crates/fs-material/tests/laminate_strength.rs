//! G0/G1/G3 synthetic laminate strength oracles, not measured material validation.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::strength::{
    LaminaFace, LaminaFailureMode, LaminaMaximumStress, LaminaStrengthProperties,
    LaminaStrengthSelection, LaminateStrengthError, PreparedLaminateStrength,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

const NAMES: [&str; 5] = [
    "tensile_1",
    "compressive_1",
    "tensile_2",
    "compressive_2",
    "shear_12",
];
const LIMITS: [f64; 5] = [120e6, 300e6, 30e6, 90e6, 15e6];

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn at(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}
fn domain(lo: f64, hi: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("T", temperature_quantity(), lo, hi)
}
fn claim(name: &str, dims: Dims, knots: &[(f64, f64)]) -> PropertyClaim {
    let constant = knots.len() == 1;
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if constant {
            PropertyValue::Scalar {
                value: knots[0].1,
                dims,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.to_vec(),
                dims,
            }
        },
        interpolation: if constant {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        validity: domain(300.0, 400.0),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored laminate strength oracle".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn claims(limits: [f64; 5]) -> Vec<PropertyClaim> {
    let mut values = vec![
        claim(DENSITY_PROPERTY, Density::DIMS, &[(300.0, 1600.0)]),
        claim("young_modulus_1", Pressure::DIMS, &[(300.0, 100e9)]),
        claim("young_modulus_2", Pressure::DIMS, &[(300.0, 20e9)]),
        claim("poisson_ratio_12", Dims::NONE, &[(300.0, 0.0)]),
        claim("shear_modulus_12", Pressure::DIMS, &[(300.0, 10e9)]),
        claim("shear_modulus_31", Pressure::DIMS, &[(300.0, 12e9)]),
        claim("shear_modulus_23", Pressure::DIMS, &[(300.0, 8e9)]),
        claim(CTE, INVERSE_TEMPERATURE_DIMS, &[(300.0, 1e-5)]),
    ];
    values.extend(
        NAMES
            .into_iter()
            .zip(limits)
            .map(|(name, value)| claim(name, Pressure::DIMS, &[(300.0, value)])),
    );
    values
}
fn card(values: Vec<PropertyClaim>) -> Arc<MaterialCard> {
    let mut claims = ClaimSet::new();
    for value in values {
        claims.insert_claim(value).unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic directional laminate".into(),
                phase: "solid".into(),
                process: "numerical strength test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
fn properties() -> LaminaStrengthProperties {
    LaminaStrengthProperties {
        tensile: [NAMES[0].into(), NAMES[2].into()],
        compressive: [NAMES[1].into(), NAMES[3].into()],
        shear: NAMES[4].into(),
    }
}
fn selection() -> LaminaStrengthSelection {
    LaminaStrengthSelection {
        properties: properties(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn layer(card: &Arc<MaterialCard>, angle: f64, h: f64) -> LaminaSourceInput {
    LaminaSourceInput {
        card: Arc::clone(card),
        reference: at(300.0),
        thickness_m: h,
        orientation_rad: angle,
        elasticity: LaminaElasticSource::OrthotropicPlane12,
        expansion: LaminaExpansionSource::Isotropic,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn prepare(inputs: Vec<LaminaSourceInput>, temperatures: &[f64]) -> PreparedLaminateStrength {
    let selections = vec![selection(); inputs.len()];
    let source =
        LaminateThermalSource::try_new(inputs, ReferenceSurface::Midsurface, || false).unwrap();
    PreparedLaminateStrength::resolve(
        source.at_temperatures(temperatures, || false).unwrap(),
        &selections,
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * scale,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g0_signed_modes_use_independent_strengths_and_physical_shear() {
    let original = card(claims(LIMITS));
    let strength = LaminaMaximumStress::resolve(
        &original,
        &at(300.0),
        properties(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let result = strength.assess([60e6, -45e6, -3e6]).unwrap();
    assert_eq!(result.mode_utilizations(), [0.5, 0.0, 0.0, 0.5, 0.2]);
    assert_eq!(
        result.critical_mode(),
        Some(LaminaFailureMode::LongitudinalTension)
    );
    assert!(!result.reaches_limit());
    for (stress, mode) in [
        ([120e6, 0.0, 0.0], LaminaFailureMode::LongitudinalTension),
        (
            [-300e6, 0.0, 0.0],
            LaminaFailureMode::LongitudinalCompression,
        ),
        ([0.0, 30e6, 0.0], LaminaFailureMode::TransverseTension),
        ([0.0, -90e6, 0.0], LaminaFailureMode::TransverseCompression),
        ([0.0, 0.0, 15e6], LaminaFailureMode::InPlaneShear),
        ([0.0, 0.0, -15e6], LaminaFailureMode::InPlaneShear),
    ] {
        let r = strength.assess(stress).unwrap();
        assert_eq!(r.utilization(), 1.0);
        assert_eq!(r.critical_mode(), Some(mode));
        assert!(r.reaches_limit());
    }
    assert_eq!(strength.assess([0.0; 3]).unwrap().critical_mode(), None);
    assert!(strength.assess([f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn g3_rotated_material_axes_preserve_weak_direction_and_shear_failure() {
    let original = card(claims(LIMITS));
    for angle in [0.0, core::f64::consts::FRAC_PI_4, -0.31] {
        let prepared = prepare(vec![layer(&original, angle, 0.002)], &[300.0]);
        let (s, c) = angle.sin_cos();
        for (stress, mode, index) in [
            ([12e6, 45e6, 0.0], LaminaFailureMode::TransverseTension, 1.5),
            ([0.0, 0.0, -30e6], LaminaFailureMode::InPlaneShear, 2.0),
        ] {
            // Independent reciprocal compliance in material axes; tensor
            // rotation of strain then supplies native section kinematics.
            let [a, b, g] = [stress[0] / 100e9, stress[1] / 20e9, stress[2] / 10e9];
            let membrane = [
                c * c * a + s * s * b - c * s * g,
                s * s * a + c * c * b + c * s * g,
                2.0 * c * s * (a - b) + (c * c - s * s) * g,
            ];
            let result = prepared.assess(membrane, [0.0; 3], || false).unwrap();
            for (actual, expected) in result
                .assessment
                .material_stress_pa()
                .into_iter()
                .zip(stress)
            {
                near(actual, expected, 50e6);
            }
            near(result.assessment.utilization(), index, index);
            assert_eq!(result.assessment.critical_mode(), Some(mode));
            assert_eq!(result.face, LaminaFace::Bottom);
        }
    }
}

#[test]
fn g1_pure_bending_assesses_both_faces_with_distinct_compressive_capacity() {
    let original = card(claims([200e6, 1000e6, 30e6, 90e6, 15e6]));
    let prepared = prepare(vec![layer(&original, 0.0, 0.002)], &[300.0]);
    for (sign, critical_face) in [(1.0, LaminaFace::Top), (-1.0, LaminaFace::Bottom)] {
        let gradient = [3.0 * sign, 0.0, 0.0];
        let result = prepared.assess([0.0; 3], gradient, || false).unwrap();
        let faces = prepared.assess_layer(0, [0.0; 3], gradient).unwrap();
        assert_eq!(result.face, critical_face);
        assert_eq!(
            result.assessment.critical_mode(),
            Some(LaminaFailureMode::LongitudinalTension)
        );
        near(result.assessment.utilization(), 1.5, 1.5);
        let compression = if sign > 0.0 { faces.bottom } else { faces.top };
        near(compression.utilization(), 0.3, 1.0);
        assert_eq!(compression.mode_utilizations()[0], 0.0);
        let response = prepared
            .thermal()
            .prepared()
            .evaluate([0.0; 3], gradient, || false)
            .unwrap();
        near(response.membrane_force_n_m[0], 0.0, 1000.0);
    }
}

#[test]
fn g1_zero_resultants_do_not_hide_thermally_weakened_layer_self_stress() {
    let mut values = claims([60e6, 300e6, 200e6, 300e6, 15e6]);
    *values
        .iter_mut()
        .find(|p| p.key.name() == NAMES[0])
        .unwrap() = claim(NAMES[0], Pressure::DIMS, &[(300.0, 60e6), (400.0, 10e6)]);
    let original = card(values);
    let prepared = prepare(vec![layer(&original, 0.0, 0.001); 2], &[300.0, 400.0]);
    // The exact free extension/curvature of equal-stiffness bimetal layers.
    let membrane = [0.0005, 0.0005, 0.0];
    let gradient = [0.75, 0.75, 0.0];
    let mechanics = prepared
        .thermal()
        .prepared()
        .evaluate(membrane, gradient, || false)
        .unwrap();
    for force in mechanics.membrane_force_n_m {
        near(force, 0.0, 1000.0);
    }
    for moment in mechanics.bending_moment_n {
        near(moment, 0.0, 1.0);
    }
    assert!(mechanics.strain_energy_j_m2 > 0.0);
    let result = prepared.assess(membrane, gradient, || false).unwrap();
    assert_eq!((result.layer, result.face), (1, LaminaFace::Top));
    near(result.assessment.material_stress_pa()[0], 25e6, 25e6);
    near(result.assessment.utilization(), 2.5, 2.5);
    assert!(result.assessment.reaches_limit());
    // The two interface stresses belong to distinct layers; neither is lost
    // by evaluating only the laminate's outside surfaces or average force.
    let cold_faces = prepared.assess_layer(0, membrane, gradient).unwrap();
    let hot_faces = prepared.assess_layer(1, membrane, gradient).unwrap();
    near(cold_faces.top.material_stress_pa()[0], 50e6, 50e6);
    near(hot_faces.bottom.material_stress_pa()[0], -50e6, 50e6);
    for i in 0..2 {
        let strength = prepared.layer_strength(i).unwrap().resolved();
        let state = &prepared.thermal().layers()[i];
        assert_eq!(strength.card_identity(), original.content_hash());
        assert_eq!(strength.query_point(), state.elasticity().query_point());
        assert_eq!(
            strength.axis_quantities(),
            state.elasticity().axis_quantities()
        );
        assert_eq!(
            strength.query_point(),
            state.expansion().current().query_point()
        );
    }
}

#[test]
fn g0_temperature_refresh_pins_strength_claims_and_native_thermal_conditions() {
    let mut values = claims(LIMITS);
    *values
        .iter_mut()
        .find(|p| p.key.name() == NAMES[0])
        .unwrap() = claim(NAMES[0], Pressure::DIMS, &[(300.0, 120e6), (400.0, 60e6)]);
    let mut alternate = claim(NAMES[0], Pressure::DIMS, &[(350.0, 5e6)]);
    alternate.validity = domain(350.0, 400.0);
    values.push(alternate);
    let original = card(values);
    let cold = prepare(vec![layer(&original, 0.0, 0.002)], &[300.0]);
    let hot = cold.at_temperatures(&[400.0], || false).unwrap();
    let (cold_limit, hot_limit) = (
        cold.layer_strength(0).unwrap(),
        hot.layer_strength(0).unwrap(),
    );
    assert_eq!(cold_limit.limits_pa()[0], 120e6);
    assert_eq!(hot_limit.limits_pa()[0], 60e6);
    assert_eq!(
        cold_limit
            .resolved()
            .property(NAMES[0])
            .unwrap()
            .answer()
            .receipt
            .selected,
        hot_limit
            .resolved()
            .property(NAMES[0])
            .unwrap()
            .answer()
            .receipt
            .selected
    );
    assert_eq!(
        hot_limit.resolved().query_point(),
        hot.thermal().layers()[0].elasticity().query_point()
    );
    assert_eq!(hot_limit.resolved().state_coordinate("T"), Some(400.0));
    assert_eq!(cold_limit.resolved().state_coordinate("T"), Some(300.0));
    let cold_result = cold.assess([0.001, 0.0, 0.0], [0.0; 3], || false).unwrap();
    // Subtract the new 0.001 free strain to apply the same mechanical load.
    let hot_result = hot.assess([0.002, 0.001, 0.0], [0.0; 3], || false).unwrap();
    near(cold_result.assessment.material_stress_pa()[0], 100e6, 100e6);
    near(hot_result.assessment.material_stress_pa()[0], 100e6, 100e6);
    assert!(!cold_result.assessment.reaches_limit());
    assert!(hot_result.assessment.reaches_limit());
    // A fresh unpinned hot selection is ambiguous; the refresh never reselects.
    assert!(
        PreparedLaminateStrength::resolve(hot.thermal().clone(), &[selection()], || false).is_err()
    );
}

#[test]
fn g0_strength_domain_and_cancellation_refuse_without_changing_the_snapshot() {
    let mut values = claims(LIMITS);
    values
        .iter_mut()
        .find(|p| p.key.name() == NAMES[4])
        .unwrap()
        .validity = domain(300.0, 350.0);
    let original = card(values);
    let cold = prepare(vec![layer(&original, 0.0, 0.002)], &[300.0]);
    // Elasticity/CTE support 400 K, but the same original shear claim does not.
    assert!(matches!(
        cold.at_temperatures(&[400.0], || false),
        Err(LaminateStrengthError::Material { layer: 0, .. })
    ));
    assert_eq!(cold.layer_strength(0).unwrap().limits_pa(), LIMITS);
    assert_eq!(
        cold.layer_strength(0)
            .unwrap()
            .resolved()
            .state_coordinate("T"),
        Some(300.0)
    );
    assert!(matches!(
        cold.assess([0.0; 3], [0.0; 3], || true),
        Err(LaminateStrengthError::Cancelled)
    ));
    assert!(matches!(
        cold.at_temperatures(&[300.0], || true),
        Err(LaminateStrengthError::Cancelled)
    ));
    assert!(matches!(
        PreparedLaminateStrength::resolve(cold.thermal().clone(), &[], || false),
        Err(LaminateStrengthError::LayerCount {
            expected: 1,
            actual: 0
        })
    ));
    assert!(cold.assess_layer(1, [0.0; 3], [0.0; 3]).is_err());
    let mut missing = selection();
    missing.properties.compressive[1] = "not_measured".into();
    assert!(matches!(
        PreparedLaminateStrength::resolve(cold.thermal().clone(), &[missing], || false),
        Err(LaminateStrengthError::Material { layer: 0, .. })
    ));
}
