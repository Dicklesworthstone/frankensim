//! G0/source binding: the finite-strain tensile seam retains the native law,
//! the exact material/model source and the complete typed external conditions.
//! Synthetic values test ownership and constitutive delegation, not calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::cards::{HyperelasticMaterial, ResolvedHyperelastic};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::strength::{IsotropicTensileLimit, ResolvedTensileHyperelastic};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const LIMIT: &str = "specified_cauchy_tensile_strength";
const ALTERNATIVE_LIMIT: &str = "alternative_cauchy_tensile_strength";

fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

fn model(mu: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            ("mu", mu, Pressure::DIMS),
            ("lambda", 3.0e6, Pressure::DIMS),
            ("stretch_limit", 4.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 280.0, 350.0),
        sources: vec![],
        provenance: provenance("synthetic hyperelastic calibration"),
    }
}

fn specimen(models: Vec<ConstitutiveModelCard>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    // Equal values distinguish property ownership from a copied scalar.
    for name in [LIMIT, ALTERNATIVE_LIMIT] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Pressure::DIMS),
                value: PropertyValue::Scalar {
                    value: 1.0e6,
                    dims: Pressure::DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 350.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance("synthetic Cauchy-strength observation"),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic polymer".into(),
            phase: "solid".into(),
            process: "binding test".into(),
            revision: 0,
        },
        claims,
        models,
    )
    .unwrap()
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}

fn limit(card: &MaterialCard, point: &QueryPoint) -> IsotropicTensileLimit {
    IsotropicTensileLimit::resolve(
        card,
        point,
        LIMIT,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn elastic(card: &MaterialCard, index: usize, point: &QueryPoint) -> ResolvedHyperelastic {
    ResolvedHyperelastic::try_from_model_card_at(&card.models()[index], point).unwrap()
}

fn refusal(result: Result<ResolvedTensileHyperelastic, MaterialStatePointError>) {
    assert!(matches!(
        result,
        Err(MaterialStatePointError::InvalidDerived { .. })
    ));
}

#[test]
fn g0_retains_original_sources_and_delegates_the_native_response() {
    let card = specimen(vec![model(2.0e6)]);
    let point = point(300.0);
    let elastic = elastic(&card, 0, &point);
    let limit = limit(&card, &point);
    let bound =
        ResolvedTensileHyperelastic::from_cauchy_limit(&card, elastic.clone(), limit.clone())
            .unwrap();
    assert_eq!(bound.elastic(), &elastic);
    assert_eq!(bound.limit(), &limit);
    assert_eq!(bound.material_card_identity(), card.content_hash());
    assert_eq!(bound.hyperelastic(), elastic.law());

    let f = [1.2, 0.1, 0.0, 0.0, 0.95, 0.05, 0.0, 0.0, 1.0];
    assert_eq!(
        bound.energy_density(&f).to_bits(),
        elastic.energy_density(&f).to_bits()
    );
    assert_eq!(
        bound.first_piola(&f).unwrap(),
        elastic.first_piola(&f).unwrap()
    );
    assert_eq!(
        bound.consistent_tangent(&f).unwrap(),
        elastic.consistent_tangent(&f).unwrap()
    );
    // A limit crossing does not damage the law or change its stored energy.
    assert!(
        bound
            .limit()
            .assess(&[2.0e6, 0.0, 0.0, 0.0, 0.0, 0.0])
            .unwrap()
            .reaches_limit()
    );
    assert_eq!(
        fs_material::hyper::cards::FiniteStrainPotential::first_piola(&bound, &f).unwrap(),
        elastic.first_piola(&f).unwrap()
    );
    assert_eq!(bound.limit(), &limit);
}

#[test]
fn g0_keeps_native_orientation_and_stretch_domain_refusals() {
    let card = specimen(vec![model(2.0e6)]);
    let point = point(300.0);
    let bound = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 0, &point),
        limit(&card, &point),
    )
    .unwrap();
    for stretch in [-1.0, 5.0, f64::NAN] {
        let f = [stretch, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert!(bound.first_piola(&f).is_err());
        assert!(bound.consistent_tangent(&f).is_err());
    }
}

#[test]
fn source_binding_refuses_strength_from_another_material_revision() {
    let card = specimen(vec![model(2.0e6)]);
    let successor =
        MaterialCard::supersede(&card, card.claims().clone(), card.models().to_vec()).unwrap();
    let point = point(300.0);
    refusal(ResolvedTensileHyperelastic::from_cauchy_limit(
        &successor,
        elastic(&successor, 0, &point),
        limit(&card, &point),
    ));
}

#[test]
fn source_binding_requires_the_exact_model_not_just_matching_parameters() {
    let original_model = model(2.0e6);
    let card = specimen(vec![original_model.clone()]);
    let point = point(300.0);
    let mut foreign_model = original_model;
    foreign_model.provenance.source = "different synthetic calibration".into();
    let foreign = ResolvedHyperelastic::try_from_model_card_at(&foreign_model, &point).unwrap();
    assert_eq!(foreign.law(), elastic(&card, 0, &point).law());
    refusal(ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        foreign,
        limit(&card, &point),
    ));

    // f64 equality alone would accept a different exact parameter encoding.
    // The model source identity must also match the actual card member.
    let mut zero_model = model(2.0e6);
    zero_model.parameters.get_mut("lambda").unwrap().value = 0.0;
    let card = specimen(vec![zero_model.clone()]);
    let mut different_bits = zero_model.clone();
    different_bits.parameters.get_mut("lambda").unwrap().value = -0.0;
    assert_eq!(zero_model, different_bits);
    assert_ne!(zero_model.content_hash(), different_bits.content_hash());
    refusal(ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        ResolvedHyperelastic::try_from_model_card_at(&different_bits, &point).unwrap(),
        limit(&card, &point),
    ));
}

#[test]
fn source_binding_compares_all_query_values_and_quantity_descriptors() {
    let card = specimen(vec![model(2.0e6)]);
    let original_point = point(300.0);
    let explicit_kelvin = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    for different_context in [
        point(310.0),
        explicit_kelvin.clone(),
        point(300.0).with("ambient_pressure", 101_325.0).unwrap(),
    ] {
        refusal(ResolvedTensileHyperelastic::from_cauchy_limit(
            &card,
            elastic(&card, 0, &different_context),
            limit(&card, &original_point),
        ));
    }
    let bound = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 0, &explicit_kelvin),
        limit(&card, &explicit_kelvin),
    )
    .unwrap();
    assert_eq!(bound.elastic().context(), &explicit_kelvin);
    assert_eq!(bound.limit().resolved().query(), explicit_kelvin);
}

#[test]
fn identity_distinguishes_model_property_and_conditions_inside_one_card() {
    let card = specimen(vec![model(2.0e6), model(3.0e6)]);
    let point = point(300.0);
    let requirements: Vec<_> = [LIMIT, ALTERNATIVE_LIMIT]
        .into_iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        })
        .collect();
    let bundle = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let first_limit = IsotropicTensileLimit::try_from_resolved(bundle.clone(), LIMIT).unwrap();
    let alternative = IsotropicTensileLimit::try_from_resolved(bundle, ALTERNATIVE_LIMIT).unwrap();
    assert_eq!(
        first_limit.resolved().identity(),
        alternative.resolved().identity()
    );
    assert_eq!(first_limit.limit_pa(), alternative.limit_pa());

    let first = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 0, &point),
        first_limit.clone(),
    )
    .unwrap();
    let same = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 0, &point),
        first_limit.clone(),
    )
    .unwrap();
    let second_model = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 1, &point),
        first_limit.clone(),
    )
    .unwrap();
    let second_property = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 0, &point),
        alternative,
    )
    .unwrap();
    let warmer_point = point.clone().with("T", 310.0).unwrap();
    let warmer = ResolvedTensileHyperelastic::from_cauchy_limit(
        &card,
        elastic(&card, 0, &warmer_point),
        first_limit.at_query_point(&card, &warmer_point).unwrap(),
    )
    .unwrap();
    assert_eq!(first.identity(), same.identity());
    assert_ne!(first.identity(), second_model.identity());
    assert_ne!(first.identity(), second_property.identity());
    assert_ne!(first.identity(), warmer.identity());
}
