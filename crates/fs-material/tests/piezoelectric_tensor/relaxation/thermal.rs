//! Source support of held coefficients, independently of the thermal integrator.
use super::*;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

fn at_query(
    card: &MaterialCard,
    mechanical: TensorBranchProperties,
    dielectric: DirectedDebyeProperties,
    query: &QueryPoint,
) -> RelaxingTensorPiezoelectric {
    let base = TensorPiezoelectric::resolve(
        card,
        query,
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    RelaxingTensorPiezoelectric::resolve(
        base,
        card,
        vec![mechanical],
        vec![dielectric],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g0_held_temperature_requires_complete_source_support_and_absolute_axis() {
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let identity = material.identity();
    let selected = material.material().identity();
    let coefficients = material.frequency_coefficients(170.0).unwrap();
    for temperature in [250.0, 300.0, 350.0, 400.0] {
        material.admit_held_temperature("T", temperature).unwrap();
    }
    for temperature in [249.0, 401.0, 0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(material.admit_held_temperature("T", temperature).is_err());
    }
    for axis in ["omega", "missing_temperature"] {
        assert!(material.admit_held_temperature(axis, 300.0).is_err());
    }
    assert_eq!(material.identity(), identity);
    assert_eq!(material.material().identity(), selected);
    assert_eq!(
        material.base().material().state_coordinate("T"),
        Some(300.0)
    );
    assert_eq!(
        material
            .frequency_coefficients(170.0)
            .unwrap()
            .stiffness_mandel_pa,
        coefficients.stiffness_mandel_pa
    );
    assert_eq!(
        material
            .frequency_coefficients(170.0)
            .unwrap()
            .permittivity_f_m,
        coefficients.permittivity_f_m
    );

    // Restrict each original bundle independently; a broad companion bundle
    // cannot extend either the base material or one auxiliary branch source.
    for (property, upper) in [("density", 310.0), ("piezo_polarization_tau", 305.0)] {
        let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
            if claim.key.name() == property {
                claim.validity = ValidityDomain::unconstrained().with("T", 295.0, upper);
            }
        });
        let material = bound(&card, &mechanical, &dielectric);
        assert!(material.admit_held_temperature("T", upper).is_ok());
        assert!(material.admit_held_temperature("T", upper + 0.01).is_err());
    }

    // These scalar declarations have no temperature restriction. Their actual
    // query still must name an absolute-temperature coordinate, not a delta T.
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
        claim.validity = ValidityDomain::unconstrained();
    });
    let difference = QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap();
    let wrong = at_query(&card, mechanical.clone(), dielectric.clone(), &difference);
    assert!(wrong.admit_held_temperature("T", 300.0).is_err());
    let absolute = QueryPoint::new()
        .with_quantity(
            "ceramic_temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap();
    let named = at_query(&card, mechanical, dielectric, &absolute);
    assert!(
        named
            .admit_held_temperature("ceramic_temperature", 330.0)
            .is_ok()
    );
    assert!(named.admit_held_temperature("T", 330.0).is_err());
}

#[test]
fn g0_frozen_temperature_samples_cannot_authorize_a_changed_ceramic_temperature() {
    for (property, flat) in [("density", false), ("piezo_memory_tau", true)] {
        let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
            if claim.key.name() == property {
                let PropertyValue::Scalar { value, dims } = &claim.value else {
                    panic!("manufactured coefficient is scalar before the explicit curve");
                };
                let (value, dims) = (*value, *dims);
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (250.0, value),
                        (300.0, value),
                        (400.0, if flat { value } else { 1.1 * value }),
                    ],
                    dims,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
        });
        let material = bound(&card, &mechanical, &dielectric);
        assert!(material.admit_held_temperature("T", 300.0).is_ok());
        assert!(material.admit_held_temperature("T", 300.01).is_err());
        assert!(material.admit_held_temperature("T", 299.99).is_err());
        // A retained sample does not contain the full curve needed to prove
        // even an apparently flat neighbouring interval; identity stays pinned.
        assert_eq!(
            material.base().material().state_coordinate("T"),
            Some(300.0)
        );
    }
}
