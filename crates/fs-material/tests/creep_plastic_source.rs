//! Manufactured source-to-history regressions for simultaneous creep and yield.
//! These coefficients exercise the implementation, not measured alloy behavior.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::creep::plastic::source::{NortonJ2Source, norton_j2_requirements};
use fs_material::creep::source::*;
use fs_material::creep::{CreepLaw, CreepTolerances};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    POST_YIELD_TANGENT_MODULUS_PROPERTY,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use std::cell::Cell;

const H: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            500.0,
        ),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured joint creep/plasticity".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn curve(name: &str, knots: Vec<(f64, f64)>) -> PropertyClaim {
    let mut claim = scalar(name, Pressure::DIMS, knots[0].1);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots,
        dims: Pressure::DIMS,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim
}
fn coefficients() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, 2e6),
        scalar(CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 1e-5),
        scalar(CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
        scalar(
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            30_000.0,
        ),
        scalar(
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            400.0,
        ),
        curve(YIELD_STRESS_PROPERTY, vec![(300.0, 4e6), (500.0, 1e6)]),
        scalar(ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 4e8),
    ]
}
fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic combined solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn source(card: &MaterialCard) -> Result<NortonJ2Source<'_>, CreepSourceError> {
    NortonJ2Source::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        H,
        || false,
    )
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 1e-4,
        energy_j_m3: 1e-6,
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1e-9),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_temperature_softening_yields_while_retaining_permanent_creep() {
    let card = card(coefficients());
    let source = source(&card).unwrap();
    let cold = source.at_temperature(300.0, 0.1).unwrap();
    let strain = [0.002, -0.001, -0.001, 0.0, 0.0, 0.0];
    let initial = cold.initial_state();
    let loaded = cold
        .creep_increment(&[0.0; 6], &initial, &strain, tolerances())
        .unwrap();
    assert_eq!(loaded.state.plastic.alpha, 0.0);
    assert!(loaded.state.creep.equivalent_creep_strain > 0.0);
    let before = loaded.state.clone();
    let hot = source.at_temperature(500.0, 0.1).unwrap();
    let relaxed = hot
        .creep_increment(&strain, &before, &strain, tolerances())
        .unwrap();
    assert!(relaxed.state.plastic.alpha > 0.0);
    assert!(relaxed.state.creep.equivalent_creep_strain > before.creep.equivalent_creep_strain);
    assert!(relaxed.stress_pa[0] < loaded.stress_pa[0]);
    near(relaxed.energy.mechanical_work_j_m3, 0.0);
    assert!(relaxed.energy.creep_dissipation_j_m3 > 0.0);
    assert!(relaxed.energy.plastic_dissipation_j_m3 > 0.0);
    near(
        relaxed.energy.physical_dissipation_j_m3,
        relaxed.energy.creep_dissipation_j_m3 + relaxed.energy.plastic_dissipation_j_m3,
    );
    assert!(relaxed.energy.final_stored_j_m3 < relaxed.energy.initial_stored_j_m3);
    assert!(relaxed.energy.residual_j_m3.abs() < tolerances().energy_j_m3);
    assert_eq!(loaded.state, before);
    assert_eq!(initial, cold.initial_state());
    assert_eq!(source.reference().properties().len(), 9);
    for resolved in [&cold, &hot] {
        assert_eq!(resolved.material().card_identity(), card.content_hash());
        for property in resolved.material().properties() {
            let original = source
                .reference()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
            assert_eq!(
                property.answer().evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
    }
}

#[test]
fn g0_complete_pinned_selection_retains_one_material_and_both_constitutive_sources() {
    let card = card(coefficients());
    let reference = source(&card).unwrap();
    let selection = MaterialPropertySelection::PinnedByProperty(
        reference
            .reference()
            .properties()
            .iter()
            .map(|property| {
                (
                    property.requirement().name().into(),
                    property.answer().receipt.selected,
                )
            })
            .collect(),
    );
    let pinned = NortonJ2Source::new(
        &card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection,
            strain_limit: 0.02,
        },
        H,
        || false,
    )
    .unwrap();
    let resolved = pinned.at_temperature(400.0, 2.0).unwrap();
    assert_eq!(resolved.material().properties().len(), 9);
    assert_eq!(pinned.creep_source().reference().properties().len(), 7);
    assert_eq!(
        resolved.material().card_identity(),
        pinned.creep_source().reference().card_identity()
    );
    near(resolved.law().plastic().yield_stress, 2.5e6);
    near(resolved.law().plastic().hardening, 4e8);
    assert_eq!(resolved.hardening_basis(), H);
    assert!(pinned.at_temperature(500.001, 1.0).is_err());
    assert!(pinned.at_temperature(400.0, 0.0).is_err());
}

#[test]
fn g1_explicit_uniaxial_tangent_uses_the_original_hardening_conversion() {
    let mut properties = coefficients();
    properties.pop();
    properties.push(scalar(
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Pressure::DIMS,
        3e8,
    ));
    let card = card(properties);
    let source = NortonJ2Source::new(
        &card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        HardeningPropertyBasis::UniaxialTangentModulus,
        || false,
    )
    .unwrap();
    let law = source.at_temperature(450.0, 1.0).unwrap();
    near(law.law().plastic().hardening, 3e8 * 1.2e9 / (1.2e9 - 3e8));
    assert!(
        law.material()
            .property(ISOTROPIC_HARDENING_MODULUS_PROPERTY)
            .is_none()
    );
    assert!(
        law.material()
            .property(POST_YIELD_TANGENT_MODULUS_PROPERTY)
            .is_some()
    );
}

#[test]
fn g0_complete_path_rejects_hidden_hardening_storage_and_unsupported_yield() {
    let mut properties = coefficients();
    *properties.last_mut().unwrap() = curve(
        ISOTROPIC_HARDENING_MODULUS_PROPERTY,
        vec![(300.0, 4e8), (400.0, 5e8), (500.0, 4e8)],
    );
    assert!(source(&card(properties)).is_err());
    let mut properties = coefficients();
    properties[7] = curve(
        YIELD_STRESS_PROPERTY,
        vec![(300.0, 4e6), (400.0, 0.0), (500.0, 1e6)],
    );
    assert!(source(&card(properties)).is_err());
    let mut properties = coefficients();
    properties[7].validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        450.0,
    );
    assert!(source(&card(properties)).is_err());
}

#[test]
fn g0_missing_or_foreign_roles_never_supply_an_implicit_yield_law() {
    let mut properties = coefficients();
    properties.retain(|property| property.key.name() != YIELD_STRESS_PROPERTY);
    assert!(source(&card(properties)).is_err());
    assert!(norton_j2_requirements(HardeningPropertyBasis::MixedPlasticStrainModuli).is_err());
    let mut properties = coefficients();
    properties[8] = scalar(ISOTROPIC_HARDENING_MODULUS_PROPERTY, Dims([0; 6]), 4e8);
    assert!(source(&card(properties)).is_err());
}

#[test]
fn g0_source_preparation_cancellation_leaves_the_original_card_usable() {
    let card = card(coefficients());
    let identity = card.content_hash();
    let calls = Cell::new(0usize);
    let cancelled = NortonJ2Source::new(
        &card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        H,
        || {
            calls.set(calls.get() + 1);
            calls.get() == 4
        },
    );
    assert!(matches!(cancelled, Err(CreepSourceError::Cancelled)));
    assert_eq!(card.content_hash(), identity);
    assert!(source(&card).unwrap().at_temperature(400.0, 1.0).is_ok());
}
