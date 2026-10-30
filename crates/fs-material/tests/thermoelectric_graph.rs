//! Kelvin source data and exact tangents through the ordinary law graph.
//! All values are manufactured numerical fixtures, not measured material data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::graph::{
    ConstitutiveGraph, EnergyBehavior, GraphError, LawNode, NodeDeclaration, NodeOutput,
    TimeParity, check_consistent_tangent,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermoelectric::{
    ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS, graph::MaterialSeebeckNode,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use std::collections::BTreeMap;

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), 350.0)
        .unwrap()
}
fn source(
    value: PropertyValue,
    policy: InterpolationPolicy,
    validity: ValidityDomain,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS),
            value,
            validity,
            interpolation: policy,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "manufactured absolute thermopower".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermoelectric".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}
fn domain() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("T", absolute(), 250.0, 500.0)
}
fn constant(s: f64) -> MaterialCard {
    source(
        PropertyValue::Scalar {
            value: s,
            dims: SEEBECK_DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
        domain(),
    )
}
fn curve(knots: Vec<(f64, f64)>) -> MaterialCard {
    source(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: SEEBECK_DIMS,
        },
        InterpolationPolicy::LinearInside,
        domain(),
    )
}
fn node(card: &MaterialCard) -> MaterialSeebeckNode {
    MaterialSeebeckNode::try_new(
        card,
        &context(),
        "T",
        [250.0, 500.0],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn inputs(t: [f64; 3], current: f64) -> [f64; 4] {
    [t[0] - 350.0, t[1] - 350.0, t[2] - 350.0, current]
}
fn check_balances(values: &[f64], t: [f64; 3], current: f64) {
    near(values[1] + values[2] + values[3] + values[4], 0.0);
    near(values[1] / t[0] + values[2] / t[1] + values[3] / t[2], 0.0);
    near(values[4], -current * values[0]);
}

#[test]
fn source_data_controls_native_kelvin_power_ports_and_exact_derivatives() {
    for s in [-0.002, 0.0, 0.002] {
        let card = constant(s);
        let law = node(&card);
        assert_eq!(law.declaration().energy, EnergyBehavior::ReversibleTransfer);
        assert_eq!(
            law.declaration()
                .inputs
                .iter()
                .map(|p| p.parity)
                .collect::<Vec<_>>(),
            [
                TimeParity::Even,
                TimeParity::Even,
                TimeParity::Even,
                TimeParity::Odd
            ],
        );
        assert_eq!(
            law.declaration()
                .outputs
                .iter()
                .map(|p| p.parity)
                .collect::<Vec<_>>(),
            [
                TimeParity::Even,
                TimeParity::Odd,
                TimeParity::Odd,
                TimeParity::Odd,
                TimeParity::Odd
            ],
        );
        assert_eq!(law.free_energy(&[], &[0.0; 4]), None);
        for current in [-2.0, 0.0, 2.0] {
            let t = [300.0, 450.0, 400.0];
            let x = inputs(t, current);
            let value = law.evaluate(&[], &x).unwrap();
            for (actual, expected) in value.outputs.iter().zip([
                100.0 * s,
                -300.0 * s * current,
                0.0,
                400.0 * s * current,
                -100.0 * s * current,
            ]) {
                near(*actual, expected);
            }
            assert_eq!(value.dissipation_rate, Some(0.0));
            assert!(value.next_state.is_empty());
            check_balances(&value.outputs, t, current);
            check_consistent_tangent("kelvin", &law, &[], &x, 2e-7).unwrap();
            let span = law.resolve_at_temperatures_k(t).unwrap();
            for material in [
                &span.material()[0],
                span.bulk_material(),
                &span.material()[1],
            ] {
                assert_eq!(material.card_identity(), card.content_hash());
                assert_eq!(
                    material.properties()[0].answer().receipt.selected,
                    law.source().bulk_material().properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
                assert_eq!(material.axis_quantities(), context().axis_quantities());
            }
        }
    }
    // S(T) = 1e-6 T - 200e-6. The independent integral means are 125e-6
    // and 225e-6 for [300,350,500], with nonzero original bulk exchange.
    let law = node(&curve(vec![(250.0, 50e-6), (500.0, 300e-6)]));
    let t = [300.0, 350.0, 500.0];
    let value = law.evaluate(&[], &inputs(t, 2.0)).unwrap();
    for (actual, expected) in value
        .outputs
        .iter()
        .zip([0.04, -0.075, -0.07, 0.225, -0.08])
    {
        near(*actual, expected);
    }
    check_balances(&value.outputs, t, 2.0);
    let reversed = law.evaluate(&[], &inputs(t, -2.0)).unwrap();
    near(reversed.outputs[0], value.outputs[0]);
    for i in 1..5 {
        near(reversed.outputs[i], -value.outputs[i]);
    }
    let tangent = law.tangent(&[], &inputs(t, 2.0)).unwrap();
    near(tangent[0], -100e-6);
    near(tangent[2], 300e-6);
    near(tangent[4], -2.0 * (125e-6 + 300.0 * 0.5e-6));
    near(tangent[9], 2.0 * (125e-6 - 225e-6));
    near(tangent[11], -0.035);
}

#[test]
fn piecewise_sources_keep_reversed_hot_bulk_and_nearly_isothermal_tangents() {
    let law = node(&curve(vec![
        (250.0, -200e-6),
        (350.0, 0.0),
        (500.0, 150e-6),
    ]));
    for t in [
        [280.0, 370.0, 470.0],
        [470.0, 370.0, 280.0],
        [400.0, 490.0, 420.0],
        [320.0, 320.0, 320.0],
        [321.0, 321.0 + 1e-10, 410.0],
    ] {
        for current in [-3.0, 0.0, 3.0] {
            let x = inputs(t, current);
            let value = law.evaluate(&[], &x).unwrap();
            check_balances(&value.outputs, t, current);
            check_consistent_tangent("kelvin", &law, &[], &x, 3e-7).unwrap();
        }
    }
    // Exact isothermal knot: chosen right derivative is half 1e-6 per input.
    // Recover it from d(heat_a)/d(T_bulk), so no production helper is reused.
    let at_knot = law.tangent(&[], &inputs([350.0; 3], 2.0)).unwrap();
    near(at_knot[5], -2.0 * 350.0 * 0.5e-6);
    let at_end = law.tangent(&[], &inputs([500.0; 3], 2.0)).unwrap();
    near(at_end[5], -2.0 * 500.0 * 0.5e-6);
}

#[test]
fn public_graph_jvp_preserves_signed_heat_and_zero_dissipation() {
    let law = node(&curve(vec![(250.0, 50e-6), (500.0, 300e-6)]));
    let x = inputs([310.0, 360.0, 420.0], -2.0);
    let direction = [0.5, -0.75, 1.0, 0.1];
    let mut graph = ConstitutiveGraph::new();
    let ports: Vec<_> = law
        .declaration()
        .inputs
        .iter()
        .map(|p| p.name.clone())
        .collect();
    graph.add_node("kelvin", Box::new(law)).unwrap();
    let bind = |values: &[f64]| {
        ports
            .iter()
            .zip(values)
            .map(|(p, v)| (("kelvin".into(), p.clone()), *v))
            .collect::<BTreeMap<_, _>>()
    };
    let version = graph.state_schema().version();
    let jvp = graph
        .linearize(version, &[], &bind(&x), &bind(&direction), || false)
        .unwrap();
    assert_eq!(jvp.value.total_dissipation, 0.0);
    assert!(jvp.value.next_state.is_empty());
    let h = 1e-4;
    let plus = std::array::from_fn::<_, 4, _>(|i| x[i] + h * direction[i]);
    let minus = std::array::from_fn::<_, 4, _>(|i| x[i] - h * direction[i]);
    let a = graph.execute(version, &[], &bind(&plus)).unwrap();
    let b = graph.execute(version, &[], &bind(&minus)).unwrap();
    for (key, action) in &jvp.output_direction {
        assert!((action - (a.outputs[key] - b.outputs[key]) / (2.0 * h)).abs() < 1e-9);
    }
    assert!(matches!(
        graph.execute_cancellable(version, &[], &bind(&x), || true),
        Err(GraphError::Cancelled)
    ));
}

fn surface(hole: bool) -> (MaterialCard, QueryPoint) {
    let temperatures = if hole {
        vec![250.0, 330.0, 370.0, 500.0]
    } else {
        vec![250.0, 500.0]
    };
    let samples = temperatures
        .iter()
        .flat_map(|&t| {
            [1e5, 2e5]
                .into_iter()
                .map(move |p| ([t, p], 1e-6 * t + 1e-9 * p - 200e-6))
        })
        .collect();
    let triangles = if hole {
        vec![[0, 2, 1], [1, 2, 3], [4, 6, 5], [5, 6, 7]]
    } else {
        vec![[0, 2, 1], [1, 2, 3]]
    };
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let card = source(
        PropertyValue::Surface2D {
            axes: ["T".into(), "p".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples,
            triangles,
            dims: SEEBECK_DIMS,
        },
        InterpolationPolicy::LinearInside,
        domain().with_quantity("p", pressure, 1e5, 2e5),
    );
    (
        card,
        context()
            .with("T", 300.0)
            .unwrap()
            .with_quantity("p", pressure, 1.5e5)
            .unwrap(),
    )
}

#[test]
fn source_support_context_and_runtime_refusals_remain_native() {
    let (card, p) = surface(false);
    let law = MaterialSeebeckNode::try_new(
        &card,
        &p,
        "T",
        [250.0, 500.0],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let value = law.evaluate(&[], &[10.0, 60.0, 120.0, 2.0]).unwrap();
    near(value.outputs[0], 110.0 * (1e-6 * 365.0 - 50e-6));
    check_consistent_tangent("surface", &law, &[], &[10.0, 60.0, 120.0, 2.0], 3e-7).unwrap();
    let span = law
        .resolve_at_temperatures_k([310.0, 360.0, 420.0])
        .unwrap();
    assert_eq!(span.bulk_material().query().axes()["p"], 1.5e5);
    let (hole, p) = surface(true);
    assert!(
        MaterialSeebeckNode::try_new(
            &hole,
            &p,
            "T",
            [250.0, 500.0],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let discrete = source(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(250.0, 1e-4), (350.0, 2e-4), (500.0, 3e-4)],
            dims: SEEBECK_DIMS,
        },
        InterpolationPolicy::TabulatedOnly,
        domain(),
    );
    assert!(
        MaterialSeebeckNode::try_new(
            &discrete,
            &context(),
            "T",
            [250.0, 500.0],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let phase = source(
        PropertyValue::Scalar {
            value: 1e-4,
            dims: SEEBECK_DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
        domain().with("liquid_mass_fraction", 0.0, 1.0),
    );
    assert!(
        MaterialSeebeckNode::try_new(
            &phase,
            &context().with("liquid_mass_fraction", 0.0).unwrap(),
            "T",
            [250.0, 500.0],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    for x in [
        [-51.0, 0.0, 0.0, 1.0],
        [0.0, 201.0, 0.0, 1.0],
        [0.0, 0.0, f64::NAN, 1.0],
        [0.0, 0.0, 0.0, f64::INFINITY],
    ] {
        assert!(law.evaluate(&[], &x).is_err());
        assert!(law.tangent(&[], &x).is_none());
    }
    assert!(matches!(
        law.evaluate(&[0.0], &[0.0; 4]),
        Err(GraphError::ArityMismatch { .. })
    ));
    assert!(matches!(
        law.evaluate(&[], &[0.0; 3]),
        Err(GraphError::ArityMismatch { .. })
    ));
}

#[test]
fn reversible_graph_claim_requires_explicit_zero_dissipation() {
    struct Report(MaterialSeebeckNode, Option<f64>);
    impl LawNode for Report {
        fn declaration(&self) -> &NodeDeclaration {
            self.0.declaration()
        }
        fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
            let mut value = self.0.evaluate(state, inputs)?;
            value.dissipation_rate = self.1;
            Ok(value)
        }
        fn tangent(&self, state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
            self.0.tangent(state, inputs)
        }
    }
    for rate in [None, Some(-1.0), Some(1.0), Some(f64::NAN)] {
        let law = Report(node(&constant(1e-4)), rate);
        let external = law
            .declaration()
            .inputs
            .iter()
            .map(|p| (("kelvin".into(), p.name.clone()), 0.0))
            .collect();
        let mut graph = ConstitutiveGraph::new();
        graph.add_node("kelvin", Box::new(law)).unwrap();
        assert!(matches!(
            graph.execute(graph.state_schema().version(), &[], &external),
            Err(GraphError::CardMismatch { .. })
        ));
    }
}
