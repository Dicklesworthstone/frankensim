//! G0/G1: original electrical paths on sensible and latent caloric intervals.
use super::*;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};
use fs_qty::Dims;

fn phase_node(
    basis: ElectricalPropertyBasis,
    electrical: PropertyValue,
    density: PropertyValue,
    maximum_fraction: f64,
) -> (MaterialCard, MaterialOhmicNode) {
    let key = PropertyKey::new("original-electrical-measurement", basis.dimensions());
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (key.clone(), electrical),
        (PropertyKey::new("density", Density::DIMS), density),
    ] {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute(), 300.0, 400.0)
                    .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured original phase transport".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductor".into(),
            phase: "confined solid-liquid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", absolute(), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &[
            basis.requirement_for(&key).unwrap(),
            ScalarPropertyRequirement::try_new(
                "density",
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = ResolvedConductorState::try_from_material_with_property(material, basis, &key)
        .unwrap()
        .with_liquid_mass_fraction_axis("f")
        .unwrap();
    let law = MaterialOhmicNode::try_new_with_liquid_mass_fraction(
        &card,
        original,
        "T",
        [300.0, 400.0],
        [0.0, maximum_fraction],
        2.0,
        0.5,
    )
    .unwrap();
    (card, law)
}

fn density() -> PropertyValue {
    PropertyValue::Scalar {
        value: 8500.0,
        dims: Density::DIMS,
    }
}

fn caloric(card: &MaterialCard, points: &[(f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        points
            .iter()
            .map(|&(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 8500.0,
            })
            .collect(),
    )
    .unwrap()
}

fn expected_bounds(basis: ElectricalPropertyBasis) -> [(OhmicDrive, f64); 2] {
    let (current, voltage) = match basis {
        ElectricalPropertyBasis::Resistivity => (0.72, 0.36),
        ElectricalPropertyBasis::Conductivity => (0.36, 0.72),
    };
    [
        (OhmicDrive::Current(-3.0), current),
        (OhmicDrive::Voltage(-12.0), voltage),
    ]
}

#[test]
fn g1_latent_feedback_retains_internal_source_knots_and_both_drive_bases() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let (card, law) = phase_node(
            basis,
            PropertyValue::Curve {
                abscissa: "f".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 4.0), (0.5, 2.0), (1.0, 3.0)],
                dims: basis.dimensions(),
            },
            density(),
            1.0,
        );
        let original = law.source().original().clone();
        let curve = caloric(
            &card,
            &[
                (-100.0, 300.0, 0.0),
                (0.0, 350.0, 0.0),
                (100.0, 350.0, 1.0),
                (200.0, 400.0, 1.0),
            ],
        );
        assert_eq!(
            curve
                .temperature_derivative_at_specific_enthalpy(50.0)
                .unwrap(),
            0.0
        );
        assert_eq!(law.liquid_mass_fraction_axis(), Some("f"));
        assert_eq!(law.liquid_mass_fraction_range(), Some([0.0, 1.0]));
        for (drive, expected) in expected_bounds(basis) {
            let bound = law
                .maximum_positive_joule_specific_enthalpy_slope(drive, &curve)
                .unwrap();
            assert!(bound >= expected);
            near(bound, expected);
            assert!(
                law.maximum_positive_joule_temperature_slope_w_k(drive)
                    .is_err()
            );
        }
        assert_eq!(
            law.maximum_positive_joule_specific_enthalpy_slope(OhmicDrive::Voltage(0.0), &curve,)
                .unwrap(),
            0.0
        );
        assert!(
            law.maximum_positive_joule_specific_enthalpy_slope(
                OhmicDrive::Current(f64::NAN),
                &curve,
            )
            .is_err()
        );
        assert_eq!(law.source().original(), &original);
    }
}

#[test]
fn g1_temperature_only_source_uses_original_cp_minimum_and_ignores_unowned_fraction() {
    let (card, original) = super::curve(
        ElectricalPropertyBasis::Resistivity,
        vec![(300.0, 2.0), (400.0, 4.0)],
    );
    let law = node(&card, original);
    let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 2.0,
                bulk_density_kg_m3: 8500.0,
            },
            HeatCapacityKnot {
                temperature_k: 400.0,
                specific_heat_capacity_j_kg_k: 8.0,
                bulk_density_kg_m3: 8500.0,
            },
        ],
    )
    .unwrap();
    assert_eq!(curve.knots()[1].specific_enthalpy_j_kg, 500.0);
    let bound = law
        .maximum_positive_joule_specific_enthalpy_slope(OhmicDrive::Current(3.0), &curve)
        .unwrap();
    // P'=0.72 W/K; the real derivative reaches 0.36 W/(J/kg) at Cp=2.
    // Dividing the endpoint power difference by 500 J/kg would give only 0.144.
    assert!(bound >= 0.36);
    near(bound, 0.36);
    let latent = caloric(&card, &[(0.0, 350.0, 0.0), (100.0, 350.0, 1.0)]);
    assert_eq!(
        law.maximum_positive_joule_specific_enthalpy_slope(OhmicDrive::Current(3.0), &latent,)
            .unwrap(),
        0.0
    );
    assert_eq!(law.liquid_mass_fraction_axis(), None);
}

fn surface(dims: Dims, reverse: bool, hole: bool) -> PropertyValue {
    let coordinates = if hole {
        vec![
            [300.0, 0.0],
            [400.0, 0.0],
            [300.0, 0.4],
            [400.0, 0.4],
            [300.0, 0.6],
            [400.0, 0.6],
            [300.0, 1.0],
            [400.0, 1.0],
        ]
    } else {
        vec![[300.0, 0.0], [400.0, 0.0], [300.0, 1.0], [400.0, 1.0]]
    };
    PropertyValue::Surface2D {
        axes: if reverse {
            ["f".into(), "T".into()]
        } else {
            ["T".into(), "f".into()]
        },
        axis_dims: if reverse {
            [Dims::NONE, Temperature::DIMS]
        } else {
            [Temperature::DIMS, Dims::NONE]
        },
        samples: coordinates
            .into_iter()
            .enumerate()
            .map(|(i, p)| {
                (
                    if reverse { [p[1], p[0]] } else { p },
                    if hole { 2.0 } else { [4.0, 2.0, 2.0, 3.0][i] },
                )
            })
            .collect(),
        triangles: if hole {
            vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]]
        } else {
            vec![[0, 1, 2], [1, 3, 2]]
        },
        dims,
    }
}

#[test]
fn g1_simultaneous_temperature_fraction_path_keeps_surface_crossings_and_axis_order() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for reverse in [false, true] {
            let (card, law) = phase_node(
                basis,
                surface(basis.dimensions(), reverse, false),
                density(),
                1.0,
            );
            let curve = caloric(&card, &[(0.0, 300.0, 0.0), (100.0, 400.0, 1.0)]);
            // The original diagonal crosses the declared triangle edge: q4→2→3.
            for (drive, expected) in expected_bounds(basis) {
                let bound = law
                    .maximum_positive_joule_specific_enthalpy_slope(drive, &curve)
                    .unwrap();
                assert!(bound >= expected);
                near(bound, expected);
            }
        }
    }
}

#[test]
fn g0_electrical_and_auxiliary_source_holes_refuse_even_at_zero_drive() {
    let basis = ElectricalPropertyBasis::Resistivity;
    for electrical_hole in [false, true] {
        let electrical = if electrical_hole {
            surface(basis.dimensions(), false, true)
        } else {
            PropertyValue::Scalar {
                value: 2.0,
                dims: basis.dimensions(),
            }
        };
        let auxiliary = if electrical_hole {
            density()
        } else {
            surface(Density::DIMS, false, true)
        };
        let (card, law) = phase_node(basis, electrical, auxiliary, 1.0);
        let curve = caloric(&card, &[(0.0, 300.0, 0.0), (100.0, 400.0, 1.0)]);
        assert!(
            law.maximum_positive_joule_specific_enthalpy_slope(OhmicDrive::Current(0.0), &curve,)
                .is_err()
        );
    }
}

#[test]
fn g0_bound_retains_original_chart_card_and_declared_phase_limits() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (card, law) = phase_node(
        basis,
        PropertyValue::Scalar {
            value: 2.0,
            dims: basis.dimensions(),
        },
        density(),
        0.75,
    );
    for points in [
        [(0.0, 300.0, 0.0), (100.0, 400.0, 1.0)],
        [(0.0, 300.0, 0.0), (100.0, 401.0, 1.0)],
    ] {
        assert!(
            law.maximum_positive_joule_specific_enthalpy_slope(
                OhmicDrive::Current(1.0),
                &caloric(&card, &points),
            )
            .is_err()
        );
    }
    let (foreign, _) = phase_node(
        basis,
        PropertyValue::Scalar {
            value: 3.0,
            dims: basis.dimensions(),
        },
        density(),
        1.0,
    );
    assert!(
        law.maximum_positive_joule_specific_enthalpy_slope(
            OhmicDrive::Current(1.0),
            &caloric(&foreign, &[(0.0, 300.0, 0.0), (100.0, 400.0, 1.0)]),
        )
        .is_err()
    );
}
