//! Synthetic coefficient fixtures exercising the source-to-native J2 boundary.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY,
    POST_YIELD_TANGENT_MODULUS_PROPERTY, ResolvedJ2Plasticity, j2_requirements,
    resolve_j2_plasticity,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

const H: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;
const ET: HardeningPropertyBasis = HardeningPropertyBasis::UniaxialTangentModulus;

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic J2 coefficient fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn claims(basis: HardeningPropertyBasis, modulus: f64) -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.3),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
        scalar(basis.property_name(), Pressure::DIMS, modulus),
    ]
}

fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic bilinear solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(card: &MaterialCard, basis: HardeningPropertyBasis) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(
        card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        basis,
        0.02,
    )
    .unwrap()
}

fn temperature_curve(property: &mut PropertyClaim, low: f64, high: f64) {
    property.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, low), (400.0, high)],
        dims: Pressure::DIMS,
    };
    property.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
    property.interpolation = InterpolationPolicy::LinearInside;
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn source_coefficients_feed_the_native_law_and_its_declared_strain_bound() {
    let source = card(claims(H, 100.0));
    let resolved = resolve(&source, H);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.material().properties().len(), 4);
    assert_eq!(resolved.hardening_basis(), H);
    near(resolved.law().elastic.youngs, 1000.0);
    near(resolved.law().elastic.poisson, 0.3);
    near(resolved.law().yield_stress, 1.0);
    near(resolved.hardening_modulus_pa(), 100.0);
    assert_eq!(
        resolved.law().card().validity.bound("strain-magnitude"),
        Some((0.0, 0.02)),
    );
    let native = fs_material::plastic::J2Plasticity::new(
        fs_material::elastic::IsotropicElastic::new(1000.0, 0.3, 0.02).unwrap(),
        1.0,
        100.0,
    )
    .unwrap();
    let strain = [0.006, -0.001, -0.001, 0.0002, 0.0, 0.0];
    let state = native.initial_state();
    assert_eq!(
        resolved.law().stress(&strain, &state),
        native.stress(&strain, &state)
    );
    assert_eq!(
        resolved.law().tangent(&strain, &state),
        native.tangent(&strain, &state)
    );
    assert_eq!(
        resolved.law().update_state(&strain, &state),
        native.update_state(&strain, &state)
    );
}

#[test]
fn mixed_hardening_requires_two_independent_sources_and_preserves_their_selection() {
    let mixed = HardeningPropertyBasis::MixedPlasticStrainModuli;
    let mut properties = claims(H, 40.0);
    properties.push(scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 160.0));
    let source = card(properties);
    let resolved = resolve(&source, mixed);
    assert_eq!(resolved.material().properties().len(), 5);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    near(resolved.hardening_modulus_pa(), 40.0);
    near(resolved.kinematic_hardening_modulus_pa(), 160.0);
    assert_eq!(resolved.kinematic_hardening_source().unwrap().requirement().name(),
        KINEMATIC_HARDENING_MODULUS_PROPERTY);
    let same_card_isotropic = resolve(&source, H);
    assert_eq!(same_card_isotropic.kinematic_hardening_modulus_pa(), 0.0);
    assert!(same_card_isotropic.kinematic_hardening_source().is_none());
    let strain = [0.006, -0.001, -0.001, 0.0002, 0.0, 0.0];
    let native = fs_material::plastic::J2Plasticity::new_mixed(
        resolved.law().elastic.clone(), 1.0, 40.0, 160.0).unwrap();
    let state = native.initial_state();
    assert_eq!(resolved.stress(&strain, &state), native.stress(&strain, &state));
    assert_eq!(resolved.update_state(&strain, &state), native.update_state(&strain, &state));
    for datum in [None,
        Some(scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, -1.0)),
        Some(scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Density::DIMS, 160.0))] {
        let mut properties = claims(H, 40.0);
        properties.extend(datum);
        assert!(resolve_j2_plasticity(&card(properties), &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly, mixed, 0.02).is_err());
    }
    // A monotonic Et source does not authorize inventing either part of a split.
    assert!(resolve_j2_plasticity(&card(claims(ET, 100.0)), &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly, mixed, 0.02).is_err());
}

#[test]
fn explicit_hardening_bases_match_but_keep_the_original_independent_receipts() {
    for hardening in [0.0, 100.0, 1e5] {
        let tangent = 1000.0 * hardening / (1000.0 + hardening);
        let by_h = resolve(&card(claims(H, hardening)), H);
        let by_et = resolve(&card(claims(ET, tangent)), ET);
        near(by_et.hardening_modulus_pa(), by_h.hardening_modulus_pa());
        near(by_et.hardening_source().value_si(), tangent);
        assert_eq!(
            by_et.hardening_source().requirement().name(),
            POST_YIELD_TANGENT_MODULUS_PROPERTY
        );
        assert!(
            by_et
                .material()
                .property(ISOTROPIC_HARDENING_MODULUS_PROPERTY)
                .is_none()
        );
    }
}

#[test]
fn adopting_a_shared_bundle_preserves_extra_properties_and_requires_a_basis_choice() {
    let mut properties = claims(H, 100.0);
    properties.push(scalar(
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Pressure::DIMS,
        50.0,
    ));
    properties.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    let source = card(properties);
    let mut requirements = j2_requirements(H).unwrap().to_vec();
    requirements.push(ET.requirement().unwrap());
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let by_h = ResolvedJ2Plasticity::try_from_material(bundle.clone(), H, 0.02).unwrap();
    let by_et = ResolvedJ2Plasticity::try_from_material(bundle.clone(), ET, 0.02).unwrap();
    assert_eq!(by_h.material(), &bundle);
    assert_eq!(by_et.material(), &bundle);
    near(by_h.hardening_modulus_pa(), 100.0);
    near(by_et.hardening_modulus_pa(), 50.0 * 1000.0 / 950.0);
    // Independent, inconsistent sources stay available; neither is overwritten.
    assert_ne!(by_h.hardening_modulus_pa(), by_et.hardening_modulus_pa());
}

#[test]
fn temperature_sources_are_evaluated_before_the_nonlinear_tangent_conversion() {
    let mut properties = claims(ET, 100.0);
    temperature_curve(&mut properties[0], 1000.0, 600.0);
    temperature_curve(&mut properties[2], 1.0, 0.5);
    temperature_curve(&mut properties[3], 100.0, 150.0);
    let source = card(properties);
    let resolved = resolve_j2_plasticity(
        &source,
        &QueryPoint::new().with("T", 350.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        ET,
        0.02,
    )
    .unwrap();
    near(resolved.law().elastic.youngs, 800.0);
    near(resolved.law().yield_stress, 0.75);
    near(resolved.hardening_source().value_si(), 125.0);
    near(resolved.hardening_modulus_pa(), 125.0 * 800.0 / 675.0);
    let wrong_order = 0.5 * (100.0 * 1000.0 / 900.0 + 150.0 * 600.0 / 450.0);
    assert!((resolved.hardening_modulus_pa() - wrong_order).abs() > 1.0);
    assert!(
        resolve_j2_plasticity(
            &source,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ET,
            0.02,
        )
        .is_err()
    );
}

#[test]
fn absent_invalid_or_unrepresentable_coefficients_refuse_without_fallbacks() {
    let mut missing = claims(H, 100.0);
    missing.pop();
    assert!(
        resolve_j2_plasticity(
            &card(missing),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            H,
            0.02,
        )
        .is_err()
    );
    for (index, dims, value) in [
        (0, Pressure::DIMS, 0.0),
        (1, Dims([0; 6]), 0.5),
        (2, Pressure::DIMS, -1.0),
        (3, Pressure::DIMS, -1.0),
        (3, Density::DIMS, 100.0),
        (0, Pressure::DIMS, f64::MAX),
    ] {
        let mut properties = claims(H, 100.0);
        let name = properties[index].key.name().to_owned();
        properties[index] = scalar(&name, dims, value);
        assert!(
            resolve_j2_plasticity(
                &card(properties),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                H,
                0.02,
            )
            .is_err()
        );
    }
    for tangent in [1000.0, 1001.0] {
        assert!(
            resolve_j2_plasticity(
                &card(claims(ET, tangent)),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                ET,
                0.02,
            )
            .is_err()
        );
    }
    let source = card(claims(H, 0.0));
    for limit in [0.0, -0.01, 0.051, f64::NAN, f64::INFINITY] {
        assert!(
            resolve_j2_plasticity(
                &source,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                H,
                limit,
            )
            .is_err()
        );
    }
}

#[test]
fn a_strain_curve_cannot_be_frozen_into_a_rate_independent_hardening_coefficient() {
    let mut properties = claims(H, 100.0);
    properties[3].value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, 100.0), (0.01, 200.0)],
        dims: Pressure::DIMS,
    };
    properties[3].validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.01);
    properties[3].interpolation = InterpolationPolicy::LinearInside;
    let source = card(properties);
    let point = QueryPoint::new().with("strain", 0.005).unwrap();
    // It is valid source data, but not a supported coefficient of this law.
    let bundle = resolve_material_state_point(
        &source,
        &point,
        &j2_requirements(H).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(ResolvedJ2Plasticity::try_from_material(bundle, H, 0.02).is_err());
}

fn static_quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

#[test]
fn pressure_surfaces_feed_the_same_native_j2_owner_and_hardening_conversion() {
    let mut properties = claims(ET, 100.0);
    for (property, values) in properties.iter_mut().zip([
        [1000.0, 600.0, 2000.0, 1200.0],
        [0.25; 4],
        [1.0, 0.5, 2.0, 1.0],
        [100.0, 150.0, 200.0, 300.0],
    ]) {
        property.value = PropertyValue::Surface2D {
            axes: ["T".into(), "pressure".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples: [[300.0, 1e5], [400.0, 1e5], [300.0, 3e5], [400.0, 3e5]]
                .into_iter()
                .zip(values)
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: property.key.quantity().dims(),
        };
        property.validity = ValidityDomain::unconstrained()
            .with("T", 300.0, 400.0)
            .with_quantity(
                "pressure",
                static_quantity(QuantityKind::Pressure),
                1e5,
                3e5,
            );
        property.interpolation = InterpolationPolicy::LinearInside;
    }
    let card = card(properties);
    for (pressure, youngs, yield_pa, tangent) in
        [(1e5, 800.0, 0.75, 125.0), (3e5, 1600.0, 1.5, 250.0)]
    {
        let point = QueryPoint::new()
            .with("T", 350.0)
            .unwrap()
            .with_quantity(
                "pressure",
                static_quantity(QuantityKind::Pressure),
                pressure,
            )
            .unwrap();
        let source = resolve_j2_plasticity(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            ET,
            0.02,
        )
        .unwrap();
        near(source.law().elastic.youngs, youngs);
        near(source.law().yield_stress, yield_pa);
        near(
            source.hardening_modulus_pa(),
            tangent * youngs / (youngs - tangent),
        );
        let native = fs_material::plastic::J2Plasticity::new(
            fs_material::elastic::IsotropicElastic::new(youngs, 0.25, 0.02).unwrap(),
            yield_pa,
            tangent * youngs / (youngs - tangent),
        )
        .unwrap();
        let strain = [0.006, -0.001, -0.001, 0.0002, 0.0, 0.0];
        let initial = native.initial_state();
        assert_eq!(
            source.update_state(&strain, &initial),
            native.update_state(&strain, &initial)
        );
        assert_eq!(
            source.stress(&strain, &initial),
            native.stress(&strain, &initial)
        );
        for property in source.material().properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn even_scalar_j2_coefficients_cannot_freeze_mechanical_applicability_conditions() {
    for quantity in [
        static_quantity(QuantityKind::Stress),
        static_quantity(QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Normal,
        }),
        QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
    ] {
        let mut properties = claims(H, 100.0);
        properties[0].validity =
            ValidityDomain::unconstrained().with_quantity("condition", quantity, 0.0, 1.0);
        let point = QueryPoint::new()
            .with_quantity("condition", quantity, 0.001)
            .unwrap();
        let source = card(properties);
        // The underlying property query is admissible, but this coefficient
        // cannot stay fixed while the constitutive owner changes strain/stress.
        let bundle = resolve_material_state_point(
            &source,
            &point,
            &j2_requirements(H).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(ResolvedJ2Plasticity::try_from_material(bundle, H, 0.02).is_err());
    }
}
