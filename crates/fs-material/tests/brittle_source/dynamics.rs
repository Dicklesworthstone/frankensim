//! G0/G1/G3: the original material bundle supplies native brittle reference mass.
use super::*;
use fs_material::brittle::source::dynamics::{
    brittle_dynamic_requirements, resolve_brittle_dynamic_material,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const RHO: f64 = 2300.0;

fn dynamic_properties() -> Vec<PropertyClaim> {
    let mut claims = properties();
    claims.push(scalar(DENSITY_PROPERTY, Density::DIMS, RHO));
    claims
}

fn dynamic_at(
    source: &MaterialCard,
    query: &QueryPoint,
) -> Result<ResolvedBrittleCrackBand, BrittleSourceError> {
    resolve_brittle_dynamic_material(
        source,
        query,
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
    )
}

#[test]
fn original_density_and_fracture_response_retain_one_complete_bundle() {
    let source = card(dynamic_properties());
    let requirements = brittle_dynamic_requirements().unwrap();
    assert_eq!(requirements.len(), 8);
    assert_eq!(
        &requirements[..7],
        brittle_crack_band_requirements().unwrap().as_slice()
    );
    let dynamic = dynamic_at(&source, &QueryPoint::new()).unwrap();
    let original_bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(dynamic.material(), &original_bundle);
    assert_eq!(dynamic.material().card_identity(), source.content_hash());
    assert_eq!(dynamic.reference_density_kg_m3().unwrap(), RHO);
    let quasistatic = resolve(&source, 0.01);
    assert_eq!(dynamic.law(), quasistatic.law());
    // A fracture-only query does not acquire unselected density by coincidence.
    assert!(quasistatic.reference_density_kg_m3().is_err());

    let damaged = dynamic
        .law()
        .trial(
            &[0.004, 0.0, 0.0, 0.001, 0.0, 0.0],
            &dynamic.law().initial_state(),
        )
        .unwrap();
    assert!(damaged.state.damage() > 0.0);
    let direct = quasistatic
        .law()
        .trial(
            &[0.004, 0.0, 0.0, 0.001, 0.0, 0.0],
            &quasistatic.law().initial_state(),
        )
        .unwrap();
    assert_eq!(damaged.stress_pa, direct.stress_pa);
    assert_eq!(damaged.energy, direct.energy);
    assert_eq!(dynamic.reference_density_kg_m3().unwrap(), RHO);

    // Adopting a larger original bundle retains every extra receipt.
    let mut claims = dynamic_properties();
    claims.push(scalar("numerical_material_marker", Dims::NONE, 17.0));
    let source = card(claims);
    let mut larger_requirements = requirements;
    larger_requirements.push(
        ScalarPropertyRequirement::try_new(
            "numerical_material_marker",
            Dims::NONE,
            ScalarAdmissibility::Finite,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &larger_requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let admitted =
        ResolvedBrittleCrackBand::try_from_material(bundle.clone(), settings(0.01)).unwrap();
    assert_eq!(admitted.material(), &bundle);
    assert_eq!(admitted.material().properties().len(), 9);
    assert_eq!(admitted.reference_density_kg_m3().unwrap(), RHO);
}

#[test]
fn absent_wrong_unit_nonpositive_and_equal_dimensional_other_quantities_refuse() {
    assert!(dynamic_at(&card(properties()), &QueryPoint::new()).is_err());
    let mut wrong_units = properties();
    wrong_units.push(scalar(DENSITY_PROPERTY, Pressure::DIMS, RHO));
    assert!(dynamic_at(&card(wrong_units), &QueryPoint::new()).is_err());

    for value in [0.0, -RHO] {
        let mut claims = properties();
        claims.push(scalar(DENSITY_PROPERTY, Density::DIMS, value));
        let source = card(claims);
        assert!(dynamic_at(&source, &QueryPoint::new()).is_err());
        // Deliberately adopt a weaker multiphysics requirement. The reference
        // mass getter must enforce its own positive-inventory contract.
        let mut requirements = brittle_crack_band_requirements().unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::Finite,
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
        let adopted = ResolvedBrittleCrackBand::try_from_material(bundle, settings(0.01)).unwrap();
        assert!(matches!(
            adopted.reference_density_kg_m3(),
            Err(BrittleSourceError::Invalid(_))
        ));
    }

    let concentration = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::MassConcentration,
        ValueForm::Static,
    ));
    assert_eq!(concentration.dims(), Density::DIMS);
    let mut claims = properties();
    let mut claim = scalar(DENSITY_PROPERTY, Density::DIMS, RHO);
    claim.key = PropertyKey::with_quantity(DENSITY_PROPERTY, concentration);
    claims.push(claim);
    let source = card(claims);
    assert!(dynamic_at(&source, &QueryPoint::new()).is_err());
    let mut requirements = brittle_crack_band_requirements().unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_with_quantity(
            DENSITY_PROPERTY,
            concentration,
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
    let adopted = ResolvedBrittleCrackBand::try_from_material(bundle, settings(0.01)).unwrap();
    assert!(matches!(
        adopted.reference_density_kg_m3(),
        Err(BrittleSourceError::Invalid(_))
    ));
}

#[test]
fn density_keeps_its_original_temperature_context_without_mutating_the_reference_snapshot() {
    let mut claims = dynamic_properties();
    claims[7].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, RHO), (600.0, 2000.0)],
        dims: Density::DIMS,
    };
    claims[7].validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        600.0,
    );
    claims[7].interpolation = InterpolationPolicy::LinearInside;
    let source = card(claims);
    let query = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let cold = dynamic_at(&source, &query).unwrap();
    let hot = cold.at_temperature(&source, "T", 600.0).unwrap();
    assert_eq!(cold.reference_density_kg_m3().unwrap(), RHO);
    assert_eq!(hot.reference_density_kg_m3().unwrap(), 2000.0);
    assert_eq!(cold.material().state_coordinate("T"), Some(300.0));
    assert_eq!(hot.material().state_coordinate("T"), Some(600.0));
    assert_eq!(
        cold.material().card_identity(),
        hot.material().card_identity()
    );
    assert_ne!(cold.material().identity(), hot.material().identity());
    assert_eq!(cold.law().parameters(), hot.law().parameters());
    assert_eq!(
        cold.material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected,
        hot.material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected
    );
    assert!(cold.at_temperature(&source, "T", 601.0).is_err());
    assert_eq!(cold.reference_density_kg_m3().unwrap(), RHO);
    // The dynamics owner retains the density it bound to reference geometry;
    // this material endpoint query itself supplies no mass-transfer operation.
}

#[test]
fn mechanical_or_rate_dependent_density_cannot_be_frozen_into_reference_mass() {
    for (axis, dims, value) in [
        ("strain", Dims::NONE, 0.001),
        ("stress", Pressure::DIMS, 1.0),
        ("strain_rate", Dims([0, 0, -1, 0, 0, 0]), 1.0),
        ("frequency", Dims([0, 0, -1, 0, 0, 0]), 1.0),
    ] {
        let mut claims = dynamic_properties();
        claims[7].validity = ValidityDomain::unconstrained().with_quantity(
            axis,
            QuantitySpec::dimensional(dims),
            0.0,
            2.0,
        );
        let source = card(claims);
        let query = QueryPoint::new()
            .with_quantity(axis, QuantitySpec::dimensional(dims), value)
            .unwrap();
        assert!(dynamic_at(&source, &query).is_err());
        // The unrelated quasistatic fracture coefficients remain usable.
        assert!(
            resolve_brittle_crack_band(
                &source,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                settings(0.01),
            )
            .is_ok()
        );
    }
}
