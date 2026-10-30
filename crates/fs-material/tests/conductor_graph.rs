//! Native source resolution and conductor physics through the public law graph.
//! Manufactured property tables verify integration; no measured-material claim.
#[path = "conductor_graph/expansion.rs"]
mod expansion;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{
    ElectricalPropertyBasis, ResolvedConductorState, graph::MaterialOhmicNode,
};
use fs_material::graph::{ConstitutiveGraph, GraphError, LawNode, check_consistent_tangent};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};
use std::collections::BTreeMap;

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn fixture(
    basis: ElectricalPropertyBasis,
    value: PropertyValue,
    interpolation: InterpolationPolicy,
    mut validity: ValidityDomain,
    pressure: Option<f64>,
    auxiliary_high: f64,
) -> (MaterialCard, ResolvedConductorState) {
    let key = PropertyKey::new("original-electrical-measurement", basis.dimensions());
    if pressure.is_some() {
        validity = validity.with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
            2e5,
        );
    }
    let mut claims = ClaimSet::new();
    for (key, value, interpolation, validity) in [
        (key.clone(), value, interpolation, validity),
        (
            PropertyKey::new("density", Density::DIMS),
            PropertyValue::Scalar {
                value: 8500.0,
                dims: Density::DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, auxiliary_high),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                interpolation,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured electrical property table".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductor".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut query = QueryPoint::new()
        .with_quantity("T", absolute(), 300.0)
        .unwrap();
    if let Some(p) = pressure {
        query = query
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                p,
            )
            .unwrap();
    }
    let material = resolve_material_state_point(
        &card,
        &query,
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
    let original =
        ResolvedConductorState::try_from_material_with_property(material, basis, &key).unwrap();
    (card, original)
}
fn curve(
    basis: ElectricalPropertyBasis,
    knots: Vec<(f64, f64)>,
) -> (MaterialCard, ResolvedConductorState) {
    let end = knots.last().unwrap().0;
    fixture(
        basis,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: basis.dimensions(),
        },
        InterpolationPolicy::LinearInside,
        ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, end),
        None,
        end,
    )
}
fn node(card: &MaterialCard, state: ResolvedConductorState) -> MaterialOhmicNode {
    MaterialOhmicNode::try_new(card, state, "T", [300.0, 400.0], 2.0, 0.5).unwrap()
}

#[test]
fn exact_native_ohmic_tangents_preserve_each_independent_source_basis() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let (card, original) = curve(basis, vec![(300.0, 2.0), (400.0, 3.0)]);
        let law = node(&card, original.clone());
        let material = law.resolve_at_temperature_k(350.0).unwrap();
        let resistance = if basis == ElectricalPropertyBasis::Resistivity {
            10.0
        } else {
            1.6
        };
        let slope = if basis == ElectricalPropertyBasis::Resistivity {
            0.04
        } else {
            -0.0064
        };
        near(material.resistance_ohm(), resistance);
        assert_eq!(material.source_basis(), basis);
        assert_eq!(material.source_key(), original.source_key());
        assert_eq!(material.material().card_identity(), card.content_hash());
        assert_eq!(material.material().query().axes()["T"], 350.0);
        assert_eq!(material.material().properties().len(), 2);
        for property in material.material().properties() {
            let previous = original
                .material()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                previous.answer().receipt.selected
            );
            assert_eq!(property.requirement(), previous.requirement());
        }
        let value = law.evaluate(&[], &[50.0, -3.0]).unwrap();
        near(value.outputs[0], -3.0 * resistance);
        near(value.outputs[1], 9.0 * resistance);
        assert_eq!(value.dissipation_rate, Some(value.outputs[1]));
        assert!(value.next_state.is_empty());
        let tangent = law.tangent(&[], &[50.0, -3.0]).unwrap();
        for (actual, expected) in
            tangent
                .iter()
                .zip([-3.0 * slope, resistance, 9.0 * slope, -6.0 * resistance])
        {
            near(*actual, expected);
        }
        check_consistent_tangent("source", &law, &[], &[50.0, -3.0], 1e-7).unwrap();
        // Evaluate sigma(350) before reciprocal; interpolating endpoint rho is wrong.
        if basis == ElectricalPropertyBasis::Conductivity {
            let reciprocal_interpolation = 4.0 * (0.5 / 2.0 + 0.5 / 3.0);
            assert!((material.resistance_ohm() - reciprocal_interpolation).abs() > 0.05);
        }
    }
}

#[test]
fn graph_accepts_zero_reference_and_material_data_alone_changes_voltage_and_heat() {
    let mut results = Vec::new();
    for scale in [1.0, 2.0] {
        let (card, original) = curve(
            ElectricalPropertyBasis::Resistivity,
            vec![(300.0, 2.0 * scale), (400.0, 3.0 * scale)],
        );
        let mut graph = ConstitutiveGraph::new();
        graph
            .add_node("wire", Box::new(node(&card, original)))
            .unwrap();
        let external: BTreeMap<_, _> = [
            (("wire".into(), "temperature_change".into()), 50.0),
            (("wire".into(), "current".into()), 3.0),
        ]
        .into();
        let value = graph
            .execute(graph.state_schema().version(), &[], &external)
            .unwrap();
        near(
            value.outputs[&("wire".into(), "voltage".into())],
            30.0 * scale,
        );
        near(value.outputs[&("wire".into(), "heat".into())], 90.0 * scale);
        near(value.total_dissipation, 90.0 * scale);
        assert!(value.next_state.is_empty());
        let direction = [
            (("wire".into(), "temperature_change".into()), 2.0),
            (("wire".into(), "current".into()), -0.5),
        ]
        .into();
        let linear = graph
            .linearize(
                graph.state_schema().version(),
                &[],
                &external,
                &direction,
                || false,
            )
            .unwrap();
        near(
            linear.output_direction[&("wire".into(), "voltage".into())],
            -4.76 * scale,
        );
        near(
            linear.output_direction[&("wire".into(), "heat".into())],
            -29.28 * scale,
        );
        results.push(value.total_dissipation);
    }
    near(results[1], 2.0 * results[0]);
}

#[test]
fn source_knots_keep_the_declared_one_sided_piecewise_tangent() {
    let (card, original) = curve(
        ElectricalPropertyBasis::Resistivity,
        vec![(300.0, 2.0), (350.0, 3.0), (400.0, 3.25)],
    );
    let law = node(&card, original);
    for (change, source_slope) in [(0.0, 0.02), (49.0, 0.02), (50.0, 0.005), (100.0, 0.005)] {
        near(
            law.tangent(&[], &[change, 1.0]).unwrap()[0],
            4.0 * source_slope,
        );
    }
    assert!(law.evaluate(&[], &[-1.0, 1.0]).is_err());
    assert!(law.evaluate(&[], &[101.0, 1.0]).is_err());
    assert!(matches!(
        law.evaluate(&[0.0], &[1.0, 1.0]),
        Err(GraphError::ArityMismatch { .. })
    ));
    assert!(law.tangent(&[], &[f64::NAN, 1.0]).is_none());
}

fn surface(with_hole: bool) -> (MaterialCard, ResolvedConductorState) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let temperatures = if with_hole {
        vec![300.0, 340.0, 360.0, 400.0]
    } else {
        vec![300.0, 400.0]
    };
    let samples = temperatures
        .iter()
        .flat_map(|&t| {
            [1e5, 2e5]
                .into_iter()
                .map(move |p| ([t, p], 2.0 + 0.01 * (t - 300.0) + (p - 1e5) / 1e5))
        })
        .collect();
    let triangles = if with_hole {
        vec![[0, 2, 1], [1, 2, 3], [4, 6, 5], [5, 6, 7]]
    } else {
        vec![[0, 2, 1], [1, 2, 3]]
    };
    fixture(
        basis,
        PropertyValue::Surface2D {
            axes: ["T".into(), "pressure".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples,
            triangles,
            dims: basis.dimensions(),
        },
        InterpolationPolicy::LinearInside,
        ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, 400.0),
        Some(1.5e5),
        400.0,
    )
}

#[test]
fn native_fixed_pressure_surface_slice_keeps_exact_tangent_and_refuses_holes() {
    let (card, original) = surface(false);
    let law = node(&card, original);
    let conductor = law.resolve_at_temperature_k(350.0).unwrap();
    near(conductor.resistance_ohm(), 4.0 / 3.0);
    near(conductor.material().query().axes()["pressure"], 1.5e5);
    near(
        law.tangent(&[], &[50.0, 2.0]).unwrap()[0],
        -8.0 * 0.01 / 9.0,
    );
    check_consistent_tangent("source", &law, &[], &[47.0, 2.0], 1e-7).unwrap();
    let (hole_card, hole_state) = surface(true);
    assert!(
        MaterialOhmicNode::try_new(&hole_card, hole_state, "T", [300.0, 400.0], 2.0, 0.5).is_err()
    );
}

#[test]
fn complete_source_bundle_discrete_data_and_geometry_refusals_remain_native() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (card, original) = curve(basis, vec![(300.0, 2.0), (400.0, 3.0)]);
    for (range, length, area) in [
        ([300.0, 300.0], 1.0, 1.0),
        ([301.0, 400.0], 1.0, 1.0),
        ([300.0, 401.0], 1.0, 1.0),
        ([300.0, 400.0], 0.0, 1.0),
        ([300.0, 400.0], 1.0, f64::NAN),
    ] {
        assert!(
            MaterialOhmicNode::try_new(&card, original.clone(), "T", range, length, area).is_err()
        );
    }
    let value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 2.0), (400.0, 3.0)],
        dims: basis.dimensions(),
    };
    let (discrete, sampled) = fixture(
        basis,
        value.clone(),
        InterpolationPolicy::TabulatedOnly,
        ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, 400.0),
        None,
        400.0,
    );
    assert!(MaterialOhmicNode::try_new(&discrete, sampled, "T", [300.0, 400.0], 1.0, 1.0).is_err());
    let (narrow_card, narrow_state) = fixture(
        basis,
        value,
        InterpolationPolicy::LinearInside,
        ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, 400.0),
        None,
        350.0,
    );
    assert!(
        MaterialOhmicNode::try_new(&narrow_card, narrow_state, "T", [300.0, 400.0], 1.0, 1.0)
            .is_err()
    );
    assert!(
        MaterialOhmicNode::try_new(
            &narrow_card,
            original.clone(),
            "T",
            [300.0, 350.0],
            1.0,
            1.0
        )
        .is_err()
    );
    // A permissive source requirement cannot weaken the native conductor's
    // positive-resistivity domain anywhere inside the advertised interval.
    let (signed_card, signed_state) = curve(basis, vec![(300.0, 2.0), (400.0, -1.0)]);
    let permissive = resolve_material_state_point(
        &signed_card,
        &signed_state.material().query(),
        &[ScalarPropertyRequirement::try_new(
            signed_state.source_key().name(),
            basis.dimensions(),
            ScalarAdmissibility::Finite,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let signed_state = ResolvedConductorState::try_from_material_with_property(
        permissive,
        basis,
        signed_state.source_key(),
    )
    .unwrap();
    assert!(
        MaterialOhmicNode::try_new(&signed_card, signed_state, "T", [300.0, 400.0], 1.0, 1.0)
            .is_err()
    );
    let (tiny_card, tiny_state) = curve(
        ElectricalPropertyBasis::Conductivity,
        vec![(300.0, 2.0), (400.0, f64::from_bits(1))],
    );
    assert!(
        MaterialOhmicNode::try_new(&tiny_card, tiny_state, "T", [300.0, 400.0], 1.0, 1.0).is_err()
    );
    assert_eq!(original.material().query().axes()["T"], 300.0);
}
