//! G0/G1/G4 uniform-field reductions from synthetic data, not measurements.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{
        CompositionBasis, FrequencyConvention, MoistureBasis, QuantityKind, SemanticType,
        StrainBasis, StrainComponent, ValueForm,
    },
};

fn card(key: &str, dims: Dims, value: PropertyValue, validity: ValidityDomain) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic uniform-field storage fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic linear material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn constant(key: &str, dims: Dims, value: f64) -> MaterialCard {
    card(
        key,
        dims,
        PropertyValue::Scalar { value, dims },
        ValidityDomain::unconstrained(),
    )
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1e-300), "{a} != {b}");
}

#[test]
fn geometry_and_material_set_capacitance_and_the_same_field_energy() {
    let source = constant(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        4e-11,
    );
    let capacitor = ResolvedCapacitor::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
        0.001,
    )
    .unwrap();
    near(capacitor.capacitance_f(), 8e-10);
    let voltage = 12.0;
    let q = capacitor.capacitance_f() * voltage;
    let field = voltage / capacitor.gap_m();
    // Independent field integral: integral(epsilon E^2/2 dV).
    near(
        capacitor.stored_energy_j(q).unwrap(),
        0.5 * 4e-11 * field * field * capacitor.area_m2() * capacitor.gap_m(),
    );
    near(
        capacitor.stored_energy_j(-q).unwrap(),
        capacitor.stored_energy_j(q).unwrap(),
    );
    assert_eq!(capacitor.material().card_identity(), source.content_hash());
    let larger =
        ResolvedCapacitor::try_from_material(capacitor.material().clone(), 0.04, 0.001).unwrap();
    near(larger.capacitance_f(), 2.0 * capacitor.capacitance_f());
    assert_eq!(larger.material(), capacitor.material());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ResolvedCapacitor::try_from_material(capacitor.material().clone(), bad, 1.0).is_err()
        );
        assert!(
            ResolvedCapacitor::try_from_material(capacitor.material().clone(), 1.0, bad).is_err()
        );
        if !bad.is_finite() {
            assert!(capacitor.stored_energy_j(bad).is_err());
        }
    }
    assert_eq!(capacitor.stored_energy_j(0.0).unwrap(), 0.0);
}

#[test]
fn magnetic_material_and_turns_set_flux_linkage_energy() {
    let source = constant(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        2e-4,
    );
    let inductor = ResolvedInductor::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.001,
        0.1,
        50,
    )
    .unwrap();
    near(inductor.inductance_h(), 0.005);
    let current = 3.0;
    let linkage = inductor.inductance_h() * current;
    let h = f64::from(inductor.turns()) * current / inductor.path_length_m();
    near(
        inductor.stored_energy_j(linkage).unwrap(),
        0.5 * 2e-4 * h * h * inductor.area_m2() * inductor.path_length_m(),
    );
    let doubled =
        ResolvedInductor::try_from_material(inductor.material().clone(), 0.001, 0.1, 100).unwrap();
    near(doubled.inductance_h(), 4.0 * inductor.inductance_h());
    assert_eq!(doubled.material(), inductor.material());
    assert!(
        ResolvedInductor::try_from_material(inductor.material().clone(), 0.001, 0.1, 0).is_err()
    );
}

#[test]
fn temperature_curve_refresh_keeps_source_and_refuses_extrapolation() {
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let source = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 4e-11), (400.0, 8e-11)],
            dims,
        },
        ValidityDomain::unconstrained().with("T", 300.0, 400.0),
    );
    let resolve = |t| {
        ResolvedCapacitor::resolve(
            &source,
            &QueryPoint::new().with("T", t).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.02,
            0.001,
        )
    };
    let cold = resolve(300.0).unwrap();
    let warm = resolve(350.0).unwrap();
    near(warm.capacitance_f() / cold.capacitance_f(), 1.5);
    assert_eq!(
        cold.material().properties()[0].answer().receipt.selected,
        warm.material().properties()[0].answer().receipt.selected
    );
    assert_ne!(cold.material(), warm.material());
    assert!(resolve(401.0).is_err());
}

#[test]
fn incompatible_quantities_dispersive_curves_and_field_support_refuse() {
    let relative = constant(ELECTRIC_PERMITTIVITY_PROPERTY, Dims([0; 6]), 4.0);
    assert!(
        ResolvedCapacitor::resolve(
            &relative,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0
        )
        .is_err()
    );
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let dispersive = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Curve {
            abscissa: "omega".into(),
            abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
            knots: vec![(1.0, 4e-11), (100.0, 8e-11)],
            dims,
        },
        ValidityDomain::unconstrained().with("omega", 1.0, 100.0),
    );
    assert!(
        ResolvedCapacitor::resolve(
            &dispersive,
            &QueryPoint::new().with("omega", 10.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0
        )
        .is_err()
    );
    let restricted = card(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        PropertyValue::Scalar {
            value: 2e-4,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained().with("B", 0.0, 1.0),
    );
    assert!(
        ResolvedInductor::resolve(
            &restricted,
            &QueryPoint::new().with("B", 0.5).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
            1
        )
        .is_err()
    );
    let frequency_limited = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Scalar { value: 4e-11, dims },
        ValidityDomain::unconstrained().with("frequency", 100.0, 1000.0),
    );
    assert!(
        ResolvedCapacitor::resolve(
            &frequency_limited,
            &QueryPoint::new().with("frequency", 500.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0
        )
        .is_err()
    );
}

fn static_quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

#[test]
fn g1_typed_fixed_external_curves_set_magnetic_storage_without_losing_sources() {
    for kind in [
        QuantityKind::Pressure,
        QuantityKind::Composition(CompositionBasis::MassFraction),
        QuantityKind::Moisture(MoistureBasis::RelativeHumidity),
        QuantityKind::MassConcentration,
        QuantityKind::AmountConcentration,
    ] {
        let quantity = static_quantity(kind);
        let source = card(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            PropertyValue::Curve {
                abscissa: "environment".into(),
                abscissa_dims: quantity.dims(),
                knots: vec![(0.25, 2e-4), (0.75, 4e-4)],
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
            ValidityDomain::unconstrained()
                .with_quantity("environment", quantity, 0.25, 0.75),
        );
        let point = QueryPoint::new()
            .with_quantity("environment", quantity, 0.5)
            .unwrap();
        let first = ResolvedInductor::resolve(
            &source, &point, MaterialPropertySelection::SingleClaimOnly, 0.001, 0.1, 50,
        )
        .unwrap();
        near(first.inductance_h(), 0.0075);
        near(linear_magnetic_permeability(first.material()).unwrap(), 3e-4);
        assert_eq!(first.material().card_identity(), source.content_hash());
        assert_eq!(first.material().axis_quantities().get("environment"), Some(&quantity));
        assert_eq!(first.material().properties()[0].source_axes().len(), 1);
        assert!(first.material().properties()[0].is_static_at_fixed_external_conditions());

        // Another fixed physical context is a new snapshot, not a live zero-work update.
        let next_point = point.clone().with("environment", 0.625).unwrap();
        let refreshed = first.material().at_query_point(&source, &next_point).unwrap();
        let second = ResolvedInductor::try_from_material(refreshed, 0.001, 0.1, 50).unwrap();
        near(second.inductance_h(), 0.00875);
        near(first.inductance_h(), 0.0075);
        assert_eq!(
            first.material().properties()[0].answer().receipt.selected,
            second.material().properties()[0].answer().receipt.selected
        );
        assert!(ResolvedInductor::resolve(
            &source,
            &point.clone().with("environment", 0.8).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.001,
            0.1,
            50,
        )
        .is_err());
    }
}

#[test]
fn g1_fixed_pressure_composition_and_temperature_surfaces_set_electric_storage() {
    let pressure = static_quantity(QuantityKind::Pressure);
    for (kind, low, high) in [
        (QuantityKind::Composition(CompositionBasis::MoleFraction), 0.0, 1.0),
        (QuantityKind::AbsoluteTemperature, 300.0, 400.0),
    ] {
        let context = static_quantity(kind);
        let source = card(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            PropertyValue::Surface2D {
                axes: ["pressure".into(), "context".into()],
                axis_dims: [pressure.dims(), context.dims()],
                samples: vec![
                    ([1e5, low], 4e-11),
                    ([3e5, low], 8e-11),
                    ([1e5, high], 6e-11),
                    ([3e5, high], 10e-11),
                ],
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            },
            ValidityDomain::unconstrained()
                .with_quantity("pressure", pressure, 1e5, 3e5)
                .with_quantity("context", context, low, high),
        );
        let point = QueryPoint::new()
            .with_quantity("pressure", pressure, 2e5)
            .unwrap()
            .with_quantity("context", context, f64::midpoint(low, high))
            .unwrap();
        let first = ResolvedCapacitor::resolve(
            &source, &point, MaterialPropertySelection::SingleClaimOnly, 0.01, 0.001,
        )
        .unwrap();
        // Independent affine field coefficient: epsilon=4e-11+4e-11*p_hat+2e-11*c_hat.
        near(first.capacitance_f(), 7e-10);
        near(first.stored_energy_j(2e-9).unwrap(), 0.5 * 4e-18 / 7e-10);
        assert_eq!(first.material().properties()[0].source_axes().len(), 2);
        let changed = first.material().at_query_point(
            &source, &point.clone().with("pressure", 2.5e5).unwrap(),
        )
        .unwrap();
        let second = ResolvedCapacitor::try_from_material(changed, 0.01, 0.001).unwrap();
        near(second.capacitance_f(), 8e-10);
        assert_eq!(first.material().card_identity(), second.material().card_identity());
        assert_eq!(
            first.material().properties()[0].answer().receipt.selected,
            second.material().properties()[0].answer().receipt.selected
        );
    }
}

#[test]
fn g4_dynamic_or_untyped_dependencies_and_applicability_bounds_are_not_fixed_context() {
    let cases = [
        ("frequency", static_quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic))),
        ("stress", static_quantity(QuantityKind::Stress)),
        ("strain", static_quantity(QuantityKind::Strain {
            basis: StrainBasis::Tensor,
            component: StrainComponent::Normal,
        })),
        ("B", QuantitySpec::dimensional(Dims([0, 1, -2, 0, -1, 0]))),
        ("E", QuantitySpec::dimensional(Dims([1, 1, -3, 0, -1, 0]))),
        ("strain_rate", QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0]))),
        ("pressure", QuantitySpec::dimensional(Dims([-1, 1, -2, 0, 0, 0]))),
        ("composition", QuantitySpec::dimensional(Dims([0; 6]))),
    ];
    for (axis, quantity) in cases {
        for varying in [false, true] {
            let dims = ELECTRIC_PERMITTIVITY_DIMS;
            let value = if varying {
                PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: quantity.dims(),
                    knots: vec![(0.25, 4e-11), (0.75, 8e-11)],
                    dims,
                }
            } else {
                PropertyValue::Scalar { value: 4e-11, dims }
            };
            let source = card(
                ELECTRIC_PERMITTIVITY_PROPERTY,
                dims,
                value,
                ValidityDomain::unconstrained().with_quantity(axis, quantity, 0.25, 0.75),
            );
            let result = ResolvedCapacitor::resolve(
                &source,
                &QueryPoint::new().with_quantity(axis, quantity, 0.5).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                1.0,
                1.0,
            );
            assert!(
                matches!(result, Err(FieldStorageError::Invalid(_))),
                "unsupported {axis} must be refused by field-storage admission; varying={varying}"
            );
        }
    }
}
