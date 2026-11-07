//! Actual next-source activation prepares geometry without inventing history.

use super::*;

#[test]
fn temperature_strength_activation_keeps_the_original_complete_source_and_plane_preparation() {
    let mut claims = properties();
    claims[2].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, FT), (600.0, 0.5 * FT)],
        dims: Pressure::DIMS,
    };
    claims[2].validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        600.0,
    );
    claims[2].interpolation = InterpolationPolicy::LinearInside;
    claims.push(scalar(DENSITY_PROPERTY, Density::DIMS, 2300.0));
    let source = card(claims);
    let mut requirements = brittle_crack_band_requirements().unwrap();
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
        &QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = ResolvedBrittleCrackBand::try_from_material(bundle, settings(0.01)).unwrap();
    let next = original.at_temperature(&source, "T", 600.0).unwrap();
    let held = uniaxial(1, 0.75 * FT);
    assert!(
        original
            .tensile_initiator(resolution())
            .unwrap()
            .activation_at(&held, None)
            .unwrap()
            .is_none()
    );
    let detector = next.tensile_initiator(resolution()).unwrap();
    let observed = detector.activation_at(&held, None).unwrap().unwrap();
    assert!(core::ptr::eq(observed.source(), &next));
    assert_eq!(observed.native().strain(), &held);
    near(observed.native().stress_pa()[1], 0.75 * FT);
    near(observed.native().tensile_strength_pa(), 0.5 * FT);
    assert!(observed.native().normal()[1].abs() > 1.0 - 1e-14);
    assert_eq!(
        observed.source().material().state_coordinate("T"),
        Some(600.0)
    );

    let receiving = observed.resolve_band(0.008).unwrap();
    assert_eq!(receiving.material(), next.material());
    assert_eq!(receiving.material().properties().len(), 8);
    assert_eq!(receiving.law().parameters(), next.law().parameters());
    assert_eq!(receiving.law().elastic(), next.law().elastic());
    assert_eq!(
        receiving.law().settings().maximum_strain,
        next.law().settings().maximum_strain
    );
    assert_eq!(
        receiving.law().settings().energy_tolerance_j_m3,
        next.law().settings().energy_tolerance_j_m3
    );
    assert_eq!(receiving.law().band_width_m(), 0.008);
    assert_eq!(receiving.law().normal(), observed.native().normal());
    assert!(observed.resolve_band(0.0).is_err());
    assert_eq!(original.law().normal(), [1.0, 0.0, 0.0]);
    assert_eq!(next.law().normal(), [1.0, 0.0, 0.0]);
    assert!(matches!(
        detector.between(&held, &held, None),
        Err(ResolvedTensileInitiationError::Native(
            TensileInitiationError::StartNotElastic { .. }
        ))
    ));
    assert!(matches!(
        detector.activation_at_controlled(&held, None, &mut || true),
        Err(ResolvedTensileInitiationError::Native(
            TensileInitiationError::Cancelled
        ))
    ));
}
