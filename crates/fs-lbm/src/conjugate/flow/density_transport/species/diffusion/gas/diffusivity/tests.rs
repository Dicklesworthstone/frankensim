//! Numerical material fixtures; these coefficients make no measured-gas claim.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Temperature;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 1e-5), (1000.0, 4e-5)],
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    250.0,
                    1000.0,
                )
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7)
                .with("conditioning", 0.4, 0.6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured temperature-dependent binary diffusivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "explicit A/B pair".into(),
            phase: "gas".into(),
            process: "analytic diffusion fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn binary_source_refreshes_temperature_and_preserves_other_material_axes() {
    let material = card();
    let point = QueryPoint::new().with("conditioning", 0.5).unwrap();
    let binding = BinaryDiffusionMaterial {
        species: ["A", "B"].map(|id| SpeciesId::new(id).unwrap()),
        material: &material,
        context: &point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let source = GasSpeciesDiffusivity::Binary(&binding);
    source.validate_species(&binding.species).unwrap();
    source
        .validate_species(&[binding.species[1].clone(), binding.species[0].clone()])
        .unwrap();
    let gate = CancelGate::new_clock_free();
    for (temperature, expected) in [(300.0, 1.2e-5), (600.0, 2.4e-5), (900.0, 3.6e-5)] {
        let actual = source.sample(0, temperature, 100_000.0, &gate).unwrap();
        assert!((actual - expected).abs() < 1e-18);
    }
    // Each query changes only physical T/p, retaining the original specimen axis.
    assert_eq!(point.axes().get("conditioning"), Some(&0.5));
    assert!(!point.axes().contains_key("temperature"));
    let outside = point.clone().with("conditioning", 0.9).unwrap();
    let wrong_condition = BinaryDiffusionMaterial {
        context: &outside,
        ..binding.clone()
    };
    assert!(
        GasSpeciesDiffusivity::Binary(&wrong_condition)
            .sample(0, 300.0, 100_000.0, &gate)
            .is_err()
    );
}

#[test]
fn binary_source_refuses_wrong_species_and_source_domain() {
    let material = card();
    let point = QueryPoint::new().with("conditioning", 0.5).unwrap();
    let binding = BinaryDiffusionMaterial {
        species: ["A", "B"].map(|id| SpeciesId::new(id).unwrap()),
        material: &material,
        context: &point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let source = GasSpeciesDiffusivity::Binary(&binding);
    assert!(
        source
            .validate_species(&["A", "C"].map(|id| SpeciesId::new(id).unwrap()))
            .is_err()
    );
    assert!(
        source
            .validate_species(&["A", "B", "C"].map(|id| SpeciesId::new(id).unwrap()))
            .is_err()
    );
    let gate = CancelGate::new_clock_free();
    for (temperature, pressure) in [
        (200.0, 1e5),
        (1100.0, 1e5),
        (300.0, 1e8),
        (0.0, 1e5),
        (300.0, f64::NAN),
    ] {
        assert!(source.sample(0, temperature, pressure, &gate).is_err());
    }
    gate.request();
    assert!(matches!(
        source.sample(0, 300.0, 1e5, &gate),
        Err(ChtError::Cancelled)
    ));
}

#[test]
fn authored_common_field_retains_zero_barriers_without_inventing_molecular_data() {
    let values = [0.0, 0.125, f64::NAN];
    let source = GasSpeciesDiffusivity::Frozen(&values);
    let species = ["A", "B", "C"].map(|id| SpeciesId::new(id).unwrap());
    source.validate_species(&species).unwrap();
    let gate = CancelGate::new_clock_free();
    assert_eq!(source.sample(0, 400.0, 2e5, &gate).unwrap(), 0.0);
    assert_eq!(source.sample(1, 800.0, 1e5, &gate).unwrap(), 0.125);
    assert!(source.sample(2, 400.0, 1e5, &gate).is_err());
    assert!(source.sample(3, 400.0, 1e5, &gate).is_err());
}
