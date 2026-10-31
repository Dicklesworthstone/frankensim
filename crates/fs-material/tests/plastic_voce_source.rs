//! G0/G1 source admission and native response for explicitly fitted Voce data.
//! These synthetic coefficients verify integration, not empirical calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY, POST_YIELD_TANGENT_MODULUS_PROPERTY,
};
use fs_material::plastic::voce::source::{
    ResolvedVoceJ2, VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
    resolve_voce_j2, voce_requirements,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm};
use fs_qty::{Density, Dims, Frequency, Pressure, QuantitySpec, Temperature};

const ISO: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;
const MIXED: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const TANGENT: HardeningPropertyBasis = HardeningPropertyBasis::UniaxialTangentModulus;

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic independently specified Voce fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3.0),
        scalar(ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 40.0),
        scalar(Q, Pressure::DIMS, 2.0),
        scalar(B, Dims([0; 6]), 80.0),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 160.0),
    ]
}

fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic saturating solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(source: &MaterialCard, basis: HardeningPropertyBasis) -> ResolvedVoceJ2 {
    resolve_voce_j2(
        source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        basis,
        0.04,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1.0),
        "{actual} != {expected}",
    );
}

#[test]
fn independent_saturation_sources_drive_the_native_return_and_retain_the_original_split() {
    let source = card(properties());
    for (basis, expected_hk, count) in [(ISO, 0.0, 6), (MIXED, 160.0, 7)] {
        let resolved = resolve(&source, basis);
        assert_eq!(resolved.material().card_identity(), source.content_hash());
        assert_eq!(resolved.material().properties().len(), count);
        assert_eq!(resolved.hardening_basis(), basis);
        assert_eq!(resolved.law().saturation_pa(), 2.0);
        assert_eq!(resolved.law().rate_per_strain(), 80.0);
        assert_eq!(resolved.backbone().hardening_modulus_pa(), 40.0);
        assert_eq!(
            resolved.backbone().kinematic_hardening_modulus_pa(),
            expected_hk
        );
        assert_eq!(resolved.saturation_source().value_si(), 2.0);
        assert_eq!(resolved.rate_source().value_si(), 80.0);
        let strain = [0.0, 0.0, 0.0, 0.01, 0.0, 0.0];
        let initial = resolved.initial_state();
        let stress = resolved.stress(&strain, &initial);
        let next = resolved.update_state(&strain, &initial);
        assert!(next.alpha > 0.0);
        // Independent pure-shear yield check with the selected translation and
        // nonlinear saturation. No material adapter computes this equation.
        let shift = (2.0 / 3.0) * expected_hk * next.plastic_strain[3];
        let equivalent = 3.0_f64.sqrt() * (stress[3] - shift).abs();
        let yield_pa = 3.0 + 40.0 * next.alpha + 2.0 * (1.0 - (-80.0 * next.alpha).exp());
        near(equivalent, yield_pa);
        assert_eq!(
            resolved.stress(&strain, &initial),
            resolved.law().stress(&strain, &initial)
        );
        assert_eq!(
            resolved.tangent(&strain, &initial),
            resolved.law().tangent(&strain, &initial)
        );
        assert_eq!(next, resolved.law().update_state(&strain, &initial));
        for property in resolved.material().properties() {
            source
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn shared_bundles_preserve_extra_physics_and_independent_receipts_without_tangent_inference() {
    let mut all = properties();
    all.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    all.push(scalar(
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Pressure::DIMS,
        500.0,
    ));
    let source = card(all);
    let mut requirements = voce_requirements(MIXED).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    requirements.push(TANGENT.requirement().unwrap());
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let resolved = ResolvedVoceJ2::try_from_material(bundle.clone(), MIXED, 0.04).unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.backbone().material(), &bundle);
    assert_eq!(resolved.material().properties().len(), 9);
    assert_eq!(resolved.backbone().hardening_modulus_pa(), 40.0);
    assert_eq!(resolved.saturation_source(), bundle.property(Q).unwrap());
    assert_eq!(resolved.rate_source(), bundle.property(B).unwrap());
    assert_eq!(
        resolved
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        7800.0
    );
    assert!(ResolvedVoceJ2::try_from_material(bundle, TANGENT, 0.04).is_err());
    assert!(voce_requirements(TANGENT).is_err());
    assert!(
        resolve_voce_j2(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            TANGENT,
            0.04,
        )
        .is_err()
    );
}

#[test]
fn ambiguous_saturation_needs_explicit_pins_and_selected_values_keep_their_receipts() {
    let original = card(properties());
    let selected = resolve(&original, MIXED);
    let pins = selected
        .material()
        .properties()
        .iter()
        .map(|p| {
            (
                p.requirement().name().to_owned(),
                p.answer().receipt.selected,
            )
        })
        .collect();
    let mut ambiguous = properties();
    ambiguous.push(scalar(Q, Pressure::DIMS, 9.0));
    let source = card(ambiguous);
    assert!(
        resolve_voce_j2(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            MIXED,
            0.04,
        )
        .is_err()
    );
    let resolved = resolve_voce_j2(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::PinnedByProperty(pins),
        MIXED,
        0.04,
    )
    .unwrap();
    assert_eq!(resolved.law().saturation_pa(), 2.0);
    assert_eq!(
        resolved.saturation_source().answer().receipt.selected,
        selected.saturation_source().answer().receipt.selected
    );
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    for property in resolved.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
}

#[test]
fn fixed_temperature_sources_are_resolved_at_one_supported_state_without_extrapolation() {
    let mut all = properties();
    for (name, low, high) in [(Q, 1.0, 3.0), (B, 20.0, 60.0)] {
        let property = all.iter_mut().find(|p| p.key.name() == name).unwrap();
        property.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, low), (400.0, high)],
            dims: property.key.quantity().dims(),
        };
        property.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        property.interpolation = InterpolationPolicy::LinearInside;
    }
    let source = card(all);
    let point = QueryPoint::new().with("T", 350.0).unwrap();
    let resolved = resolve_voce_j2(
        &source,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        MIXED,
        0.04,
    )
    .unwrap();
    near(resolved.law().saturation_pa(), 2.0);
    near(resolved.law().rate_per_strain(), 40.0);
    for property in [resolved.saturation_source(), resolved.rate_source()] {
        assert_eq!(property.source_axes(), &[("T".into(), Temperature::DIMS)]);
        assert_eq!(
            property.answer().receipt.query_point,
            vec![("T".into(), 350.0)]
        );
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert!(
        resolve_voce_j2(
            &source,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            MIXED,
            0.04,
        )
        .is_err()
    );
}

#[test]
fn missing_wrong_quantity_and_invalid_saturation_data_refuse_but_explicit_zero_is_native_j2() {
    for missing in [Q, B] {
        let mut all = properties();
        all.retain(|p| p.key.name() != missing);
        assert!(
            resolve_voce_j2(
                &card(all),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                ISO,
                0.04,
            )
            .is_err()
        );
    }
    for (name, dims, value) in [
        (Q, Pressure::DIMS, -1.0),
        (Q, Density::DIMS, 2.0),
        (B, Dims([0; 6]), 0.0),
        (B, Dims([0; 6]), -1.0),
        (B, Frequency::DIMS, 80.0),
    ] {
        let mut all = properties();
        *all.iter_mut().find(|p| p.key.name() == name).unwrap() = scalar(name, dims, value);
        assert!(
            resolve_voce_j2(
                &card(all),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                ISO,
                0.04,
            )
            .is_err()
        );
    }
    let mut all = properties();
    *all.iter_mut().find(|p| p.key.name() == Q).unwrap() = scalar(Q, Pressure::DIMS, 0.0);
    let resolved = resolve(&card(all), MIXED);
    let backbone = resolved.backbone();
    let mut state = resolved.initial_state();
    for shear in [0.01, 0.008, -0.012] {
        let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
        assert_eq!(
            resolved.stress(&strain, &state),
            backbone.stress(&strain, &state)
        );
        assert_eq!(
            resolved.tangent(&strain, &state),
            backbone.tangent(&strain, &state)
        );
        let next = resolved.update_state(&strain, &state);
        assert_eq!(next, backbone.update_state(&strain, &state));
        state = next;
    }
}

#[test]
fn saturation_coefficients_cannot_freeze_evolving_mechanical_or_rate_conditions() {
    let static_quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for name in [Q, B] {
        for quantity in [
            static_quantity(QuantityKind::Stress),
            static_quantity(QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Normal,
            }),
            QuantitySpec::dimensional(Frequency::DIMS),
            QuantitySpec::dimensional(Dims([0; 6])),
        ] {
            let mut all = properties();
            let property = all.iter_mut().find(|p| p.key.name() == name).unwrap();
            property.validity =
                ValidityDomain::unconstrained().with_quantity("condition", quantity, 0.0, 1.0);
            let source = card(all);
            let point = QueryPoint::new()
                .with_quantity("condition", quantity, 0.1)
                .unwrap();
            // A legitimate point query does not authorize freezing its evolving
            // physical applicability while the constitutive history advances.
            let bundle = resolve_material_state_point(
                &source,
                &point,
                &voce_requirements(ISO).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            assert!(ResolvedVoceJ2::try_from_material(bundle, ISO, 0.04).is_err());
        }
    }
    let mut all = properties();
    let property = all.iter_mut().find(|p| p.key.name() == Q).unwrap();
    property.value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, 2.0), (0.02, 4.0)],
        dims: Pressure::DIMS,
    };
    property.validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.02);
    property.interpolation = InterpolationPolicy::LinearInside;
    assert!(
        resolve_voce_j2(
            &card(all),
            &QueryPoint::new().with("strain", 0.01).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ISO,
            0.04,
        )
        .is_err()
    );
}
