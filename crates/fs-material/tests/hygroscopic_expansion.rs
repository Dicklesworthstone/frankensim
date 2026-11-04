//! Analytic source-path fixtures, not a measured wood swelling data pack.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};

fn qty(k: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(k, ValueForm::Static))
}
fn point(m: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            m,
        )
        .unwrap()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
}
fn card(reference: f64, count: usize, maximum: f64) -> MaterialCard {
    let mut set = ClaimSet::new();
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            0.,
            maximum,
        )
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
            293.15,
        );
    let make = |key, value| PropertyClaim {
        key,
        value,
        validity: validity.clone(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "authored analytic hygroscopic derivative/datum".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    for (name, beta) in ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES
        .into_iter()
        .zip([0.1, 0.2, -0.05])
        .take(count)
    {
        set.insert_claim(make(
            PropertyKey::new(name, Dims::NONE),
            PropertyValue::Curve {
                abscissa: "moisture_content".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0., beta), (0.1, beta), (0.2, 3. * beta), (0.3, 2. * beta)],
                dims: Dims::NONE,
            },
        ))
        .unwrap();
    }
    set.insert_claim(make(
        PropertyKey::with_quantity(
            HYGROSCOPIC_REFERENCE_MOISTURE_PROPERTY,
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
        ),
        PropertyValue::Scalar {
            value: reference,
            dims: Dims::NONE,
        },
    ))
    .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
#[test]
fn exact_dry_basis_integral_uses_three_signed_axes_and_fixed_datum() {
    let c = card(0.1, 3, 0.3);
    let path = integrate_orthotropic_hygroscopic_expansion(
        &c,
        &point(0.1),
        &point(0.3),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    // Integral: .1*(beta+3beta)/2 + .1*(3beta+2beta)/2 = .45 beta.
    for (a, b) in path
        .free_principal_strains()
        .iter()
        .zip([0.045, 0.09, -0.0225])
    {
        assert!((a - b).abs() < 1e-14);
    }
    assert_eq!(
        path,
        integrate_orthotropic_hygroscopic_expansion(
            &c,
            &point(0.1),
            &point(0.3),
            MaterialPropertySelection::SingleClaimOnly
        )
        .unwrap()
    );
    let dry = integrate_orthotropic_hygroscopic_expansion(
        &c,
        &point(0.1),
        &point(0.),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    for (a, b) in dry
        .free_principal_strains()
        .iter()
        .zip([-0.01, -0.02, 0.005])
    {
        assert!((a - b).abs() < 1e-14);
    }
    for key in ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES {
        assert_eq!(
            path.reference()
                .property(key)
                .unwrap()
                .answer()
                .receipt
                .selected,
            path.current()
                .property(key)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
}
#[test]
fn missing_direction_wrong_datum_changed_temperature_or_wet_basis_refuses() {
    let run = |c: &MaterialCard, p: &QueryPoint| {
        integrate_orthotropic_hygroscopic_expansion(
            c,
            &point(0.1),
            p,
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    assert!(run(&card(0.1, 2, 0.3), &point(0.2)).is_err());
    assert!(run(&card(0.12, 3, 0.3), &point(0.2)).is_err());
    assert!(run(&card(0.1, 3, 0.15), &point(0.2)).is_err());
    assert!(
        run(
            &card(0.1, 3, 0.3),
            &point(0.2).with("temperature", 300.).unwrap()
        )
        .is_err()
    );
    let wet = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::WetMass)),
            0.2,
        )
        .unwrap();
    assert!(run(&card(0.1, 3, 0.3), &wet).is_err());
}
