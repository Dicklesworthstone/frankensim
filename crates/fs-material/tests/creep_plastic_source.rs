//! Manufactured source-to-history regressions for simultaneous creep and yield.
//! These coefficients exercise the implementation, not measured alloy behavior.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::source::{
    NortonJ2Source, norton_j2_requirements, norton_j2_strain_hardening_requirements,
};
use fs_material::creep::recovery::plastic::source::{
    RecoveringNortonJ2Source, recovering_norton_j2_requirements,
};
use fs_material::creep::recovery::source::{
    CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY, CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
    CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
};
use fs_material::creep::source::*;
use fs_material::creep::thermal::expansion::CreepHistoryStrain;
use fs_material::creep::{CreepLaw, CreepTolerances, NortonCreep, StrainHardeningCreep};
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
fn primary_coefficients() -> Vec<PropertyClaim> {
    let mut properties = coefficients();
    properties.extend([
        scalar(
            CREEP_HARDENING_REFERENCE_STRAIN_PROPERTY,
            Dims([0; 6]),
            1e-4,
        ),
        scalar(CREEP_HARDENING_EXPONENT_PROPERTY, Dims([0; 6]), 1.0),
        scalar(CREEP_RESIDUAL_RATE_FRACTION_PROPERTY, Dims([0; 6]), 0.2),
    ]);
    properties
}
fn recovery_coefficients() -> Vec<PropertyClaim> {
    let mut properties = primary_coefficients();
    properties.extend([
        scalar(CREEP_RECOVERY_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 0.5),
        scalar(
            CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            20_000.0,
        ),
        scalar(
            CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            350.0,
        ),
    ]);
    properties
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
fn primary_source(card: &MaterialCard) -> Result<NortonJ2Source<'_>, CreepSourceError> {
    NortonJ2Source::new_with_strain_hardening(
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
fn recovery_source(card: &MaterialCard) -> Result<RecoveringNortonJ2Source<'_>, CreepSourceError> {
    RecoveringNortonJ2Source::new(
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
fn g1_recovering_source_retains_three_histories_through_thermal_yield_and_dwell() {
    let card = card(recovery_coefficients());
    let source = recovery_source(&card).unwrap();
    assert_eq!(recovering_norton_j2_requirements(H).unwrap().len(), 15);
    assert_eq!(source.reference().properties().len(), 15);
    assert_eq!(source.recovery_source().reference().properties().len(), 13);
    let cold = source.at_temperature(300.0, 0.1).unwrap();
    let hot = source.at_temperature(500.0, 0.1).unwrap();
    for law in [&cold, &hot] {
        assert_eq!(law.material().card_identity(), card.content_hash());
        assert_eq!(law.material().properties().len(), 15);
        assert_eq!(law.hardening_basis(), H);
        assert_eq!(law.kinetics().duration_s(), 0.1);
        for property in law.material().properties() {
            assert_eq!(
                property.answer().receipt.selected,
                source
                    .reference()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
            );
        }
    }
    near(
        hot.law().recovery_rate_s() / cold.law().recovery_rate_s(),
        (20_000.0_f64 / fs_material::creep::GAS_CONSTANT_J_MOL_K * (1.0 / 300.0 - 1.0 / 500.0))
            .exp(),
    );
    let strain = [0.002, -0.001, -0.001, 0.0, 0.0, 0.0];
    let virgin = cold.initial_state();
    let loaded = cold
        .creep_increment(&[0.0; 6], &virgin, &strain, tolerances())
        .unwrap();
    assert_eq!(loaded.state.native.plastic.alpha, 0.0);
    assert!(loaded.state.native.creep.equivalent_creep_strain > 0.0);
    let frozen = loaded.state.clone();
    let yielded = hot
        .creep_increment(&strain, &loaded.state, &strain, tolerances())
        .unwrap();
    assert_eq!(loaded.state, frozen);
    assert!(yielded.state.native.plastic.alpha > 0.0);
    assert!(
        yielded.state.native.creep.equivalent_creep_strain
            > loaded.state.native.creep.equivalent_creep_strain
    );
    assert!(yielded.state.resistance_strain < yielded.state.native.creep.equivalent_creep_strain);
    assert!(yielded.energy.creep_dissipation_j_m3 > 0.0);
    assert!(yielded.energy.plastic_dissipation_j_m3 > 0.0);
    near(
        yielded.energy.physical_dissipation_j_m3,
        yielded.energy.creep_dissipation_j_m3 + yielded.energy.plastic_dissipation_j_m3,
    );
    assert_eq!(yielded.energy.mechanical_work_j_m3, 0.0);
    assert!(yielded.energy.residual_j_m3.abs() <= tolerances().energy_j_m3);

    let unloaded_strain = yielded.state.inelastic_strain();
    let unloaded = hot
        .creep_increment(&strain, &yielded.state, &unloaded_strain, tolerances())
        .unwrap();
    let dwell_law = source.at_temperature(500.0, 2.0).unwrap();
    let dwell = dwell_law
        .creep_increment(
            &unloaded_strain,
            &unloaded.state,
            &unloaded_strain,
            tolerances(),
        )
        .unwrap();
    for endpoint in [&unloaded, &dwell] {
        assert!(
            (endpoint.state.native.creep.equivalent_creep_strain
                - yielded.state.native.creep.equivalent_creep_strain)
                .abs()
                < 1e-15
        );
        assert_eq!(endpoint.state.native.plastic, yielded.state.native.plastic);
        assert!(endpoint.energy.physical_dissipation_j_m3.abs() < 1e-12);
    }
    near(
        dwell.state.resistance_strain,
        unloaded.state.resistance_strain / (1.0 + 2.0 * dwell_law.law().recovery_rate_s()),
    );
    assert!(dwell.state.resistance_strain < unloaded.state.resistance_strain);
    near(
        dwell.energy.final_stored_j_m3,
        dwell.energy.initial_stored_j_m3,
    );
    assert_eq!(virgin, cold.initial_state());
}

#[test]
fn g0_recovering_joint_source_keeps_original_yield_and_hardening_admission() {
    for absent in [7, 8] {
        let incomplete = recovery_coefficients()
            .into_iter()
            .enumerate()
            .filter(|(index, _)| *index != absent)
            .map(|(_, property)| property)
            .collect();
        assert!(recovery_source(&card(incomplete)).is_err());
    }
    for (index, knots) in [
        (7, vec![(300.0, 4e6), (400.0, -1.0), (500.0, 1e6)]),
        (8, vec![(300.0, 4e8), (400.0, 5e8), (500.0, 4e8)]),
    ] {
        let mut changed = recovery_coefficients();
        changed[index] = curve(changed[index].key.name(), knots);
        assert!(recovery_source(&card(changed)).is_err());
    }
    assert!(
        recovering_norton_j2_requirements(HardeningPropertyBasis::MixedPlasticStrainModuli)
            .is_err()
    );
    let mut tangent_basis = recovery_coefficients();
    tangent_basis[8] = scalar(POST_YIELD_TANGENT_MODULUS_PROPERTY, Pressure::DIMS, 3e8);
    let card = card(tangent_basis);
    let source = RecoveringNortonJ2Source::new(
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
    let resolved = source.at_temperature(450.0, 0.5).unwrap();
    near(resolved.law().norton_j2().plastic().hardening, 4e8);
    assert_eq!(resolved.material().properties().len(), 15);
    assert!(source.at_temperature(501.0, 0.5).is_err());
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

#[test]
fn g1_primary_creep_and_heated_yield_keep_twelve_pinned_claims_and_native_history() {
    let card = card(primary_coefficients());
    let original = primary_source(&card).unwrap();
    let secondary = source(&card).unwrap();
    assert_eq!(secondary.reference().properties().len(), 9);
    assert!(secondary.creep_source().strain_hardening().is_none());
    assert_eq!(
        norton_j2_strain_hardening_requirements(H).unwrap().len(),
        12
    );
    assert_eq!(original.reference().properties().len(), 12);
    assert_eq!(original.creep_source().reference().properties().len(), 10);
    let parameters = StrainHardeningCreep {
        reference_strain: 1e-4,
        exponent: 1.0,
        residual_rate_fraction: 0.2,
    };
    let selection = MaterialPropertySelection::PinnedByProperty(
        original
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
    let pinned = NortonJ2Source::new_with_strain_hardening(
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
    let cold = pinned.at_temperature(300.0, 1.0).unwrap();
    let strain = [0.002, -0.001, -0.001, 0.0, 0.0, 0.0];
    let loaded = cold
        .creep_increment(&[0.0; 6], &cold.initial_state(), &strain, tolerances())
        .unwrap();
    assert_eq!(loaded.state.plastic.alpha, 0.0);
    assert!(loaded.state.creep.equivalent_creep_strain > 0.0);
    let hot = pinned.at_temperature(500.0, 0.5).unwrap();
    let native = NortonJ2::new(
        NortonCreep::new(
            hot.law().creep().elastic().clone(),
            pinned.creep_source().parameters().clone(),
            500.0,
            0.5,
        )
        .unwrap()
        .with_strain_hardening(parameters)
        .unwrap(),
        hot.law().plastic().clone(),
    )
    .unwrap();
    let relaxed = hot
        .creep_increment(&strain, &loaded.state, &strain, tolerances())
        .unwrap();
    assert_eq!(
        relaxed,
        native
            .advance_with_energy(&strain, &loaded.state, &strain, tolerances())
            .unwrap(),
    );
    assert!(relaxed.state.plastic.alpha > 0.0);
    assert!(
        relaxed.state.creep.equivalent_creep_strain > loaded.state.creep.equivalent_creep_strain
    );
    assert!(relaxed.energy.creep_dissipation_j_m3 > 0.0);
    assert!(relaxed.energy.plastic_dissipation_j_m3 > 0.0);
    near(relaxed.energy.mechanical_work_j_m3, 0.0);
    for resolved in [&cold, &hot] {
        assert_eq!(resolved.law().creep().strain_hardening(), Some(&parameters));
        assert_eq!(resolved.material().card_identity(), card.content_hash());
        assert_eq!(resolved.material().properties().len(), 12);
        for property in resolved.material().properties() {
            assert_eq!(
                property.answer().receipt.selected,
                original
                    .reference()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
            );
            assert_eq!(
                property.answer().evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
    }
    assert_eq!(hot.material().state_coordinate("T"), Some(500.0));
}

#[test]
fn g0_combined_primary_source_rejects_missing_invalid_or_changing_primary_claims() {
    for index in 9..12 {
        let incomplete = primary_coefficients()
            .into_iter()
            .enumerate()
            .filter(|(position, _)| *position != index)
            .map(|(_, property)| property)
            .collect();
        assert!(
            primary_source(&card(incomplete)).is_err(),
            "missing primary claim {index}"
        );
        let mut changed = primary_coefficients();
        let (value, dims) = match changed[index].value {
            PropertyValue::Scalar { value, dims } => (value, dims),
            _ => unreachable!(),
        };
        changed[index].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, value), (400.0, value * 1.1), (500.0, value)],
            dims,
        };
        changed[index].interpolation = InterpolationPolicy::LinearInside;
        assert!(matches!(
            primary_source(&card(changed)),
            Err(CreepSourceError::Invalid(_))
        ));
    }
    for (index, value) in [(9, 0.0), (10, -1.0), (11, 1.01)] {
        let mut invalid = primary_coefficients();
        if let PropertyValue::Scalar { value: actual, .. } = &mut invalid[index].value {
            *actual = value;
        }
        assert!(primary_source(&card(invalid)).is_err());
    }
    assert!(
        norton_j2_strain_hardening_requirements(HardeningPropertyBasis::MixedPlasticStrainModuli)
            .is_err()
    );
}
