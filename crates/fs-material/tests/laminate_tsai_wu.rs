//! G0/G1/G3 coupled lamina failure oracles with synthetic source data.
//! These verify equations and native consumers, not a measured failure population.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::strength::tsai_wu::{
    LaminaTsaiWu, PreparedLaminateTsaiWu, TsaiWuInteractionSelection,
};
use fs_material::laminate::strength::{
    LaminaFace, LaminaMaximumStress, LaminaStrengthProperties, LaminaStrengthSelection,
    LaminateStrengthError, PreparedLaminateStrength,
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

const NAMES: [&str; 5] = ["Xt", "Xc", "Yt", "Yc", "S12"];
const INTERACTION: &str = "normalized_tsai_wu_interaction";
fn tq() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn at(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", tq(), t).unwrap()
}
fn claim(name: &str, dims: Dims, low: f64, high: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, low), (400.0, high)],
            dims,
        },
        interpolation: InterpolationPolicy::LinearInside,
        validity: ValidityDomain::unconstrained().with_quantity("T", tq(), 300.0, 400.0),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic Tsai-Wu oracle".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn claims(limits: [f64; 5], interaction: Option<(f64, f64)>) -> Vec<PropertyClaim> {
    let mut values = vec![
        claim(DENSITY_PROPERTY, Density::DIMS, 1600.0, 1600.0),
        claim("young_modulus_1", Pressure::DIMS, 100e9, 100e9),
        claim("young_modulus_2", Pressure::DIMS, 20e9, 20e9),
        claim("poisson_ratio_12", Dims::NONE, 0.0, 0.0),
        claim("shear_modulus_12", Pressure::DIMS, 10e9, 10e9),
        claim("shear_modulus_31", Pressure::DIMS, 12e9, 12e9),
        claim("shear_modulus_23", Pressure::DIMS, 8e9, 8e9),
        claim(CTE, INVERSE_TEMPERATURE_DIMS, 1e-5, 1e-5),
    ];
    values.extend(
        NAMES
            .into_iter()
            .zip(limits)
            .map(|(name, value)| claim(name, Pressure::DIMS, value, value)),
    );
    if let Some((lo, hi)) = interaction {
        values.push(claim(INTERACTION, Dims::NONE, lo, hi));
    }
    values
}
fn card(values: Vec<PropertyClaim>, tag: &str) -> Arc<MaterialCard> {
    let mut set = ClaimSet::new();
    for value in values {
        set.insert_claim(value).unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic lamina".into(),
                phase: "solid".into(),
                process: tag.into(),
                revision: 0,
            },
            set,
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
fn strength(card: &MaterialCard) -> LaminaMaximumStress {
    LaminaMaximumStress::resolve(
        card,
        &at(300.0),
        properties(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn interaction() -> TsaiWuInteractionSelection {
    TsaiWuInteractionSelection {
        property: INTERACTION.into(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn coupled(card: &MaterialCard) -> LaminaTsaiWu {
    LaminaTsaiWu::resolve(
        card,
        strength(card),
        INTERACTION,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn native(card: &Arc<MaterialCard>, angle: f64) -> PreparedLaminateStrength {
    let source = LaminateThermalSource::try_new(
        vec![LaminaSourceInput {
            card: Arc::clone(card),
            reference: at(300.0),
            thickness_m: 0.002,
            orientation_rad: angle,
            elasticity: LaminaElasticSource::OrthotropicPlane12,
            expansion: LaminaExpansionSource::Isotropic,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        ReferenceSurface::Midsurface,
        || false,
    )
    .unwrap();
    PreparedLaminateStrength::resolve(
        source.at_temperatures(&[300.0], || false).unwrap(),
        &[LaminaStrengthSelection {
            properties: properties(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        || false,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * b.abs().max(1e-300),
        "{a:e} != {b:e}"
    );
}
fn strain_for(stress: [f64; 3], angle: f64) -> [f64; 3] {
    let [a, b, g] = [stress[0] / 100e9, stress[1] / 20e9, stress[2] / 10e9];
    let (s, c) = angle.sin_cos();
    [
        c * c * a + s * s * b - c * s * g,
        s * s * a + c * c * b + c * s * g,
        2.0 * c * s * (a - b) + (c * c - s * s) * g,
    ]
}

#[test]
fn g0_coupled_envelope_retains_uniaxial_anchors_and_detects_interaction() {
    let limits = [120.0, 300.0, 30.0, 90.0, 15.0];
    let original = card(claims(limits, Some((0.35, 0.35))), "anchors");
    let law = coupled(&original);
    for stress in [
        [120.0, 0.0, 0.0],
        [-300.0, 0.0, 0.0],
        [0.0, 30.0, 0.0],
        [0.0, -90.0, 0.0],
        [0.0, 0.0, 15.0],
        [0.0, 0.0, -15.0],
    ] {
        near(law.assess(stress).unwrap().utilization(), 1.0);
        assert!(
            law.assess(stress.map(|v| 1.001 * v))
                .unwrap()
                .reaches_limit()
        );
        assert!(
            !law.assess(stress.map(|v| 0.999 * v))
                .unwrap()
                .reaches_limit()
        );
    }
    let symmetric = card(claims([100.0; 5], Some((0.5, 0.5))), "interaction");
    let law = coupled(&symmetric);
    let stress = [60.0; 3];
    assert!(!law.strength().assess(stress).unwrap().reaches_limit());
    // Independent polynomial: .36+.36+2*.5*.36+.36 = 1.44.
    near(law.assess(stress).unwrap().utilization(), 1.2);
    near(law.assess([60.0, 60.0, -60.0]).unwrap().utilization(), 1.2);
    assert!(law.assess(stress).unwrap().reaches_limit());
    assert_eq!(law.assess([0.0; 3]).unwrap().utilization(), 0.0);
}

#[test]
fn g0_stable_radial_index_avoids_squared_stress_overflow_and_underflow() {
    let original = card(
        claims([100.0, 200.0, 50.0, 125.0, 25.0], Some((-0.45, -0.45))),
        "scaling",
    );
    let law = coupled(&original);
    // Negative linear term exercises the cancellation-safe compression branch.
    for stress in [[-25.0, -10.0, 4.0], [25.0, 10.0, -4.0]] {
        let result = law.assess(stress).unwrap();
        let [s1, s2, t] = stress;
        let a = s1 * s1 / (100.0 * 200.0) + s2 * s2 / (50.0 * 125.0)
            - 0.9 * s1 * s2 / (100.0_f64 * 200.0 * 50.0 * 125.0).sqrt()
            + t * t / 25.0_f64.powi(2);
        let b = s1 * (1.0 / 100.0 - 1.0 / 200.0) + s2 * (1.0 / 50.0 - 1.0 / 125.0);
        let r = result.utilization();
        near(a / (r * r) + b / r, 1.0);
        for scale in [1e-200, 1e200] {
            near(
                law.assess(stress.map(|v| scale * v)).unwrap().utilization() / scale,
                r,
            );
        }
        near(
            law.assess(stress.map(|v| v / r)).unwrap().utilization(),
            1.0,
        );
    }
    assert!(law.assess([f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn g1_native_rotated_bending_faces_and_thermal_refresh_use_same_sources() {
    let original = card(claims([100e6; 5], Some((0.5, -0.5))), "native faces");
    for angle in [0.0, 0.37, core::f64::consts::FRAC_PI_4] {
        let prepared = native(&original, angle);
        let strain = strain_for([60e6; 3], angle);
        let gradient = strain.map(|v| v / 0.001);
        assert!(
            !prepared
                .assess([0.0; 3], gradient, || false)
                .unwrap()
                .assessment
                .reaches_limit()
        );
        let coupled =
            PreparedLaminateTsaiWu::from_maximum_stress(prepared, &[interaction()], || false)
                .unwrap();
        let result = coupled.assess([0.0; 3], gradient, || false).unwrap();
        near(result.assessment.utilization(), 1.2);
        assert!(result.assessment.reaches_limit());
        assert_eq!(result.face, LaminaFace::Bottom);
        let faces = coupled.assess_layer(0, [0.0; 3], gradient).unwrap();
        for (actual, expected) in faces.top.material_stress_pa().into_iter().zip([60e6; 3]) {
            near(actual, expected);
        }
        near(faces.bottom.utilization(), faces.top.utilization());
        let refreshed = coupled.at_temperatures(&[400.0], || false).unwrap();
        // Free thermal extension at alpha=1e-5/K offsets the same native
        // thermal eigenstrain, leaving the imposed bending stress unchanged.
        let hot = refreshed
            .assess([0.001, 0.001, 0.0], gradient, || false)
            .unwrap();
        near(hot.assessment.utilization(), 0.72_f64.sqrt());
        assert!(!hot.assessment.reaches_limit());
        let model = refreshed.layer_strength(0).unwrap();
        assert_eq!(
            model.strength().resolved().query(),
            model.interaction_source().query()
        );
        assert_eq!(
            model.interaction_source().query(),
            refreshed.thermal().layers()[0].elasticity().query()
        );
        assert_eq!(
            coupled.layer_strength(0).unwrap().normalized_interaction(),
            0.5
        );
    }
}

#[test]
fn g3_temperature_refresh_pins_interaction_claim_and_refuses_foreign_context() {
    let mut values = claims([100.0; 5], Some((0.25, -0.25)));
    let mut alternate = claim(INTERACTION, Dims::NONE, 0.75, 0.75);
    alternate.validity = ValidityDomain::unconstrained().with_quantity("T", tq(), 400.0, 400.0);
    alternate.provenance.source = "different synthetic interaction experiment".into();
    values.push(alternate);
    let original = card(values, "pinned sources");
    let law = coupled(&original);
    let refreshed = law.at_query_point(&original, &at(400.0)).unwrap();
    assert_eq!(refreshed.normalized_interaction(), -0.25);
    let hot_strength = law
        .strength()
        .at_query_point(&original, &at(400.0))
        .unwrap();
    assert!(
        LaminaTsaiWu::resolve(
            &original,
            hot_strength,
            INTERACTION,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    assert!(
        LaminaTsaiWu::try_from_resolved(
            law.strength().clone(),
            refreshed.interaction_source().clone(),
            INTERACTION
        )
        .is_err()
    );
    let foreign = card(
        claims([100.0; 5], Some((0.25, -0.25))),
        "different specimen",
    );
    assert!(
        LaminaTsaiWu::resolve(
            &foreign,
            law.strength().clone(),
            INTERACTION,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    assert!(law.at_query_point(&original, &at(401.0)).is_err());
    assert_eq!(law.normalized_interaction(), 0.25);
}

#[test]
fn g4_missing_nonconvex_and_cancelled_interactions_refuse_without_mutation() {
    for coefficient in [None, Some((1.0, 1.0)), Some((-1.0, -1.0)), Some((1.1, 1.1))] {
        let original = card(claims([100e6; 5], coefficient), "unsupported interaction");
        assert!(
            LaminaTsaiWu::resolve(
                &original,
                strength(&original),
                INTERACTION,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let original = card(claims([100e6; 5], Some((0.2, 0.2))), "refusals");
    assert!(matches!(
        PreparedLaminateTsaiWu::from_maximum_stress(native(&original, 0.0), &[], || false),
        Err(LaminateStrengthError::LayerCount {
            expected: 1,
            actual: 0
        })
    ));
    assert!(matches!(
        PreparedLaminateTsaiWu::from_maximum_stress(
            native(&original, 0.0),
            &[interaction()],
            || true
        ),
        Err(LaminateStrengthError::Cancelled)
    ));
    let prepared = PreparedLaminateTsaiWu::from_maximum_stress(
        native(&original, 0.0),
        &[interaction()],
        || false,
    )
    .unwrap();
    let before = prepared.assess([0.0001; 3], [0.0; 3], || false).unwrap();
    assert_eq!(
        prepared.assess([0.0001; 3], [0.0; 3], || true),
        Err(LaminateStrengthError::Cancelled)
    );
    assert!(prepared.assess_layer(1, [0.0; 3], [0.0; 3]).is_err());
    assert!(prepared.at_temperatures(&[401.0], || false).is_err());
    assert!(matches!(
        prepared.at_temperatures(&[300.0], || true),
        Err(LaminateStrengthError::Cancelled)
    ));
    assert_eq!(
        before,
        prepared.assess([0.0001; 3], [0.0; 3], || false).unwrap()
    );
}
