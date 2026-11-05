//! Source-to-native porous compression with explicitly authored coefficients.
//! These numerical fixtures do not claim an empirical foam calibration.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::porous::PorousPlasticity;
use fs_material::plastic::porous::source::{
    POROUS_HARDENING_MODULUS_PROPERTY as H, POROUS_PRESSURE_SENSITIVITY_PROPERTY as A,
    PorousSourceError, ResolvedPorousPlasticity, porous_requirements, resolve_porous_plasticity,
};
use fs_material::plastic::source::{
    ISOTROPIC_HARDENING_MODULUS_PROPERTY, POST_YIELD_TANGENT_MODULUS_PROPERTY,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic porous plasticity coefficients".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn properties() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 2.0),
        scalar(H, Pressure::DIMS, 40.0),
        scalar(A, Dims([0; 6]), 1.2),
    ]
}
fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored porous solid".into(),
            phase: "solid".into(),
            process: "numerical coefficient fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(card: &MaterialCard) -> Result<ResolvedPorousPlasticity, PorousSourceError> {
    resolve_porous_plasticity(
        card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.04,
    )
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn explicit_sources_drive_native_hydrostatic_and_mixed_compression() {
    let source = card(properties());
    let resolved = resolve(&source).unwrap();
    let native = PorousPlasticity::new(
        IsotropicElastic::new(1000.0, 0.2, 0.04).unwrap(),
        2.0,
        40.0,
        1.2,
    )
    .unwrap();
    assert_eq!(resolved.law(), &native);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.material().properties().len(), 5);
    assert_eq!(resolved.hardening_source().value_si(), 40.0);
    assert_eq!(resolved.pressure_sensitivity_source().value_si(), 1.2);
    for property in resolved.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let initial = resolved.initial_state();
    let hydro = [-0.01, -0.01, -0.01, 0.0, 0.0, 0.0];
    let trial = resolved.response(&hydro, &initial).unwrap();
    let bulk = 1000.0 / (3.0 * (1.0 - 2.0 * 0.2));
    let beta = 1.2 / (1.0_f64 + 1.2_f64.powi(2) / 9.0).sqrt();
    let expected_alpha = (beta * bulk * 0.03 - 2.0) / (40.0 + bulk * beta * beta);
    near(trial.state.alpha, expected_alpha);
    near(
        trial.state.plastic_strain[..3].iter().sum(),
        -beta * expected_alpha,
    );
    assert!(trial.state.plastic_strain[..3].iter().sum::<f64>() < 0.0);
    near(
        beta * trial.stress_pa[0].abs(),
        2.0 + 40.0 * trial.state.alpha,
    );
    assert_eq!(initial, native.initial_state());
    let mut state = initial;
    for strain in [
        hydro,
        [-0.004, -0.004, -0.004, 0.0, 0.0, 0.0],
        [-0.009, -0.007, -0.005, 0.001, 0.0, 0.0],
    ] {
        assert_eq!(
            resolved.response(&strain, &state),
            native.response(&strain, &state)
        );
        assert_eq!(
            resolved.stress(&strain, &state),
            native.stress(&strain, &state)
        );
        assert_eq!(
            resolved.tangent(&strain, &state),
            native.tangent(&strain, &state)
        );
        let next = resolved.update_state(&strain, &state);
        assert_eq!(next, native.update_state(&strain, &state));
        state = next;
    }
}

#[test]
fn fixed_temperature_bundle_keeps_extra_physics_query_and_original_source_receipts() {
    let mut all = properties();
    for (index, low, high) in [(0, 1000.0, 800.0), (2, 2.0, 1.0), (4, 1.0, 1.4)] {
        let dims = all[index].key.quantity().dims();
        all[index].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, low), (400.0, high)],
            dims,
        };
        all[index].validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        all[index].interpolation = InterpolationPolicy::LinearInside;
    }
    all.push(scalar(DENSITY_PROPERTY, Density::DIMS, 300.0));
    let source = card(all);
    let mut requirements = porous_requirements().unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let point = QueryPoint::new().with("T", 350.0).unwrap();
    let bundle = resolve_material_state_point(
        &source,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let resolved = ResolvedPorousPlasticity::try_from_material(bundle.clone(), 0.04).unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.material().query().axes().get("T"), Some(&350.0));
    assert_eq!(resolved.material().properties().len(), 6);
    near(resolved.law().elastic.youngs, 900.0);
    near(resolved.law().yield_stress, 1.5);
    near(resolved.law().pressure_sensitivity, 1.2);
    assert_eq!(resolved.hardening_source(), bundle.property(H).unwrap());
    assert_eq!(
        resolved.pressure_sensitivity_source(),
        bundle.property(A).unwrap()
    );
    for property in bundle.properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert!(
        resolve_porous_plasticity(
            &source,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.04
        )
        .is_err()
    );
}

#[test]
fn missing_wrong_or_other_law_calibrations_refuse_and_zero_limits_are_explicit() {
    for index in [3, 4] {
        let mut all = properties();
        all.remove(index);
        assert!(resolve(&card(all)).is_err());
    }
    for (index, name, dims, value) in [
        (3, H, Pressure::DIMS, -1.0),
        (3, H, Density::DIMS, 40.0),
        (4, A, Dims([0; 6]), -1.0),
        (4, A, Pressure::DIMS, 1.2),
    ] {
        let mut all = properties();
        all[index] = scalar(name, dims, value);
        assert!(resolve(&card(all)).is_err());
    }
    let mut other_law = properties();
    other_law[3] = scalar(ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 40.0);
    other_law.push(scalar(
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Pressure::DIMS,
        30.0,
    ));
    assert!(resolve(&card(other_law)).is_err());
    let mut zero = properties();
    zero[3] = scalar(H, Pressure::DIMS, 0.0);
    zero[4] = scalar(A, Dims([0; 6]), 0.0);
    let law = resolve(&card(zero)).unwrap();
    let initial = law.initial_state();
    assert_eq!(
        law.response(&[-0.01, -0.01, -0.01, 0.0, 0.0, 0.0], &initial)
            .unwrap()
            .state,
        initial
    );
    for limit in [0.0, 0.051, f64::NAN] {
        assert!(
            resolve_porous_plasticity(
                &card(properties()),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                limit
            )
            .is_err()
        );
    }
}

#[test]
fn ambiguity_requires_explicit_selection_and_evolving_sources_cannot_be_frozen() {
    let original = resolve(&card(properties())).unwrap();
    let pins = original
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
    ambiguous.push(scalar(A, Dims([0; 6]), 2.0));
    let source = card(ambiguous);
    assert!(resolve(&source).is_err());
    let selected = resolve_porous_plasticity(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::PinnedByProperty(pins),
        0.04,
    )
    .unwrap();
    assert_eq!(
        selected
            .pressure_sensitivity_source()
            .answer()
            .receipt
            .selected,
        original
            .pressure_sensitivity_source()
            .answer()
            .receipt
            .selected
    );
    assert_eq!(selected.law().pressure_sensitivity, 1.2);

    let mut evolving = properties();
    evolving[4].value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, 1.2), (0.02, 2.0)],
        dims: Dims([0; 6]),
    };
    evolving[4].validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.02);
    evolving[4].interpolation = InterpolationPolicy::LinearInside;
    let point = QueryPoint::new().with("strain", 0.01).unwrap();
    let bundle = resolve_material_state_point(
        &card(evolving),
        &point,
        &porous_requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        ResolvedPorousPlasticity::try_from_material(bundle, 0.04),
        Err(PorousSourceError::Invalid(_))
    ));

    // Even a separately resolved shared bundle must retain the exact units.
    let mut wrong = properties();
    wrong[4] = scalar(A, Pressure::DIMS, 1.2);
    let mut requirements = porous_requirements().unwrap();
    requirements[4] =
        ScalarPropertyRequirement::try_new(A, Pressure::DIMS, ScalarAdmissibility::NonNegative)
            .unwrap();
    let bundle = resolve_material_state_point(
        &card(wrong),
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        ResolvedPorousPlasticity::try_from_material(bundle, 0.04),
        Err(PorousSourceError::Invalid(_))
    ));
}
