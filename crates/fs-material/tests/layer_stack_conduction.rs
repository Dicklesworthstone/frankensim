//! Local source domains on the original ordered conductivity-integral owner.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, LayerStackConduction, TwoLayerConduction,
};
use fs_material::thermal::transport::{
    ThermalTransportError,
    stack::{
        LayerStackSourceQuery, PreparedLayerStackConduction,
        transient::{PreparedLayeredThermalBody, TransientLayerSourceQuery},
    },
    storage::{HeatStorageAmount, HeatStorageComponentQuery},
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::DepositionError;

fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance}"
    );
}

fn constant(lo: f64, hi: f64, k: f64) -> ConductivityCurve {
    ConductivityCurve::new(vec![[lo, k], [hi, k]]).unwrap()
}

#[test]
fn disjoint_three_layer_sources_close_series_heat_without_source_extrapolation() {
    let curves = [
        constant(690.0, 810.0, 2.0),
        constant(490.0, 610.0, 4.0),
        constant(290.0, 410.0, 1.0),
    ];
    let lengths = [0.02, 0.04, 0.01];
    let contacts = [0.01; 2];
    let stack = LayerStackConduction::new(&curves, &lengths, &contacts, 0.003).unwrap();
    let got = stack.evaluate([800.0, 300.0], config(), || Ok(())).unwrap();
    // Three 0.01 m² K/W layers and two 0.01 m² K/W contacts carry 10000 W/m².
    close(got.power_a_to_b_w, 30.0, 1e-9);
    for (actual, expected) in got.layer_temperatures_k().iter().flatten().zip(
        [[800.0, 700.0], [600.0, 500.0], [400.0, 300.0]]
            .iter()
            .flatten(),
    ) {
        close(*actual, *expected, 1e-9);
    }
    close(got.endpoint_derivative_w_k[0], 0.06, 1e-12);
    close(got.endpoint_derivative_w_k[1], -0.06, 1e-12);
    let reversed_curves = [curves[2].clone(), curves[1].clone(), curves[0].clone()];
    let reversed =
        LayerStackConduction::new(&reversed_curves, &[0.01, 0.04, 0.02], &contacts, 0.003)
            .unwrap()
            .evaluate([300.0, 800.0], config(), || Ok(()))
            .unwrap();
    close(reversed.power_a_to_b_w, -got.power_a_to_b_w, 1e-9);
    for (actual, expected) in reversed
        .layer_temperatures_k()
        .iter()
        .zip(got.layer_temperatures_k().iter().rev())
    {
        close(actual[0], expected[1], 1e-9);
        close(actual[1], expected[0], 1e-9);
    }
}

#[test]
fn touching_sources_admit_an_exact_perfect_contact_interface() {
    let curves = [
        constant(600.0, 800.0, 1.0),
        constant(400.0, 600.0, 1.0),
        constant(300.0, 400.0, 1.0),
    ];
    let stack = LayerStackConduction::new(&curves, &[0.02, 0.02, 0.01], &[0.0; 2], 0.01).unwrap();
    let got = stack.evaluate([800.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 100.0, 1e-10);
    assert_eq!(
        got.layer_temperatures_k(),
        &[[800.0, 600.0], [600.0, 400.0], [400.0, 300.0]]
    );
}

#[test]
fn nonlinear_four_layer_sources_retain_integrated_flux_and_endpoint_tangents() {
    let faces = [
        [720.0, 640.0],
        [620.0, 530.0],
        [500.0, 410.0],
        [390.0, 300.0],
    ];
    let bounds = [
        [610.0, 750.0],
        [510.0, 640.0],
        [390.0, 520.0],
        [280.0, 410.0],
    ];
    let laws = [[2.0, 0.01], [12.0, -0.01], [1.0, 0.02], [8.0, -0.01]];
    let curves: Vec<_> = bounds
        .iter()
        .zip(laws)
        .map(|(span, [a, b])| {
            ConductivityCurve::new(span.iter().map(|&t| [t, a + b * t]).collect()).unwrap()
        })
        .collect();
    let lengths: Vec<_> = faces
        .iter()
        .zip(laws)
        .map(|([left, right], [a, b])| (left - right) * (a + b * 0.5 * (left + right)) / 1000.0)
        .collect();
    let contacts: Vec<_> = faces
        .windows(2)
        .map(|pair| (pair[0][1] - pair[1][0]) / 1000.0)
        .collect();
    let stack = LayerStackConduction::new(&curves, &lengths, &contacts, 0.005).unwrap();
    let got = stack.evaluate([720.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 5.0, 1e-9);
    for (actual, expected) in got
        .layer_temperatures_k()
        .iter()
        .flatten()
        .zip(faces.iter().flatten())
    {
        close(*actual, *expected, 1e-7);
    }
    for endpoint in 0..2 {
        let mut plus = [720.0, 300.0];
        let mut minus = plus;
        plus[endpoint] += 0.001;
        minus[endpoint] -= 0.001;
        let difference = (stack
            .evaluate(plus, config(), || Ok(()))
            .unwrap()
            .power_a_to_b_w
            - stack
                .evaluate(minus, config(), || Ok(()))
                .unwrap()
                .power_a_to_b_w)
            / 0.002;
        close(got.endpoint_derivative_w_k[endpoint], difference, 1e-7);
    }
}

#[test]
fn source_misses_budget_and_cancellation_do_not_publish_a_stack_response() {
    let curves = [
        constant(690.0, 810.0, 2.0),
        constant(490.0, 610.0, 4.0),
        constant(290.0, 410.0, 1.0),
    ];
    let lengths = [0.02, 0.04, 0.01];
    let unsupported = LayerStackConduction::new(&curves, &lengths, &[0.0; 2], 0.003).unwrap();
    assert_eq!(
        unsupported.evaluate([800.0, 300.0], config(), || Ok(())),
        Err(DepositionError::MaterialDomain)
    );
    let supported = LayerStackConduction::new(&curves, &lengths, &[0.01; 2], 0.003).unwrap();
    assert!(matches!(
        supported.evaluate(
            [800.0, 300.0],
            LayerConductionConfig {
                max_iterations: 1,
                ..config()
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    let mut polls = 0;
    let expected = supported
        .evaluate([800.0, 300.0], config(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut again = 0;
    assert_eq!(
        supported.evaluate([800.0, 300.0], config(), || {
            again += 1;
            if again == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    );
    assert_eq!(
        supported
            .evaluate([800.0, 300.0], config(), || Ok(()))
            .unwrap(),
        expected
    );
}

#[test]
fn two_layer_stacks_use_the_same_physical_owner_and_diagnostics() {
    let curves = [constant(500.0, 800.0, 8.0), constant(250.0, 450.0, 2.0)];
    let pair = TwoLayerConduction::new([&curves[0], &curves[1]], [0.02, 0.03], 0.004, 0.04)
        .unwrap()
        .evaluate([600.0, 300.0], config(), || Ok(()))
        .unwrap();
    let stack = LayerStackConduction::new(&curves, &[0.02, 0.03], &[0.04], 0.004)
        .unwrap()
        .evaluate([600.0, 300.0], config(), || Ok(()))
        .unwrap();
    assert_eq!(stack.power_a_to_b_w, pair.power_a_to_b_w);
    assert_eq!(stack.endpoint_derivative_w_k, pair.endpoint_derivative_w_k);
    assert_eq!(stack.iterations, pair.iterations);
    assert_eq!(
        stack.layer_temperatures_k(),
        &[
            [600.0, pair.interface_temperature_k[0]],
            [pair.interface_temperature_k[1], 300.0]
        ]
    );
}

fn temperature_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}

fn source(name: &str, bounds: [f64; 2], k: f64) -> MaterialCard {
    let k_dims = Dims([1, 1, -3, -1, 0, 0]);
    let cp_dims = Dims([2, 0, -2, -1, 0, 0]);
    let mut claims = ClaimSet::new();
    for (name, dims, value, interpolation) in [
        (
            "k",
            k_dims,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(bounds[0], k), (bounds[1], k)],
                dims: k_dims,
            },
            InterpolationPolicy::LinearInside,
        ),
        (
            "cp",
            cp_dims,
            PropertyValue::Scalar {
                value: 100.0,
                dims: cp_dims,
            },
            InterpolationPolicy::ConstantWithinValidity,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    bounds[0],
                    bounds[1],
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "authored series-conduction oracle".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "source-local thermal fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn original_source_stack_admits_disjoint_domains_without_inventing_a_common_storage_span() {
    let bounds = [[690.0, 810.0], [490.0, 610.0], [290.0, 410.0]];
    let cards = [
        source("hot", bounds[0], 2.0),
        source("middle", bounds[1], 4.0),
        source("cold", bounds[2], 1.0),
    ];
    let lower = bounds.map(|span| temperature_point(span[0]));
    let upper = bounds.map(|span| temperature_point(span[1]));
    let lengths = [0.02, 0.04, 0.01];
    let key = PropertyKey::new("k", Dims([1, 1, -3, -1, 0, 0]));
    let queries: [LayerStackSourceQuery<'_>; 3] = std::array::from_fn(|i| LayerStackSourceQuery {
        id: i as u64,
        card: &cards[i],
        conductivity: &key,
        thickness_m: lengths[i],
        temperature_axis: "T",
        lower: &lower[i],
        upper: &upper[i],
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    });
    let stack = PreparedLayerStackConduction::try_new(&queries, &[0.01; 2], 0.003).unwrap();
    assert_eq!(stack.common_temperature_span_k(), None);
    let got = stack.evaluate([800.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 30.0, 1e-9);
    for (i, layer) in stack.layers().iter().enumerate() {
        assert_eq!(layer.material_card_identity(), cards[i].content_hash());
        assert_eq!(
            layer.query_endpoints(),
            &[lower[i].clone(), upper[i].clone()]
        );
        assert_eq!(layer.source_line().knots.first().unwrap().0, bounds[i][0]);
        assert_eq!(layer.source_line().knots.last().unwrap().0, bounds[i][1]);
    }
    // The finite-cell adapter still declares one common evolving caloric span.
    let cp = PropertyKey::new("cp", Dims([2, 0, -2, -1, 0, 0]));
    let transient: Vec<_> = queries
        .iter()
        .enumerate()
        .map(|(i, &conduction)| TransientLayerSourceQuery {
            conduction,
            capacity: HeatStorageComponentQuery {
                id: i as u64,
                card: &cards[i],
                capacity: &cp,
                amount: HeatStorageAmount::MassKg(1.0),
                temperature_axis: "T",
                lower: &lower[i],
                upper: &upper[i],
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            mass_kg: 1.0,
            cell_count: 1,
        })
        .collect();
    assert!(matches!(
        PreparedLayeredThermalBody::try_new(&transient, &[0.01; 2], 0.003, 300.0),
        Err(ThermalTransportError::Invalid(
            "transient layered body needs a positive-width common conductivity interval"
        ))
    ));
}

#[test]
fn source_adapter_reports_a_touching_common_interval_without_dropping_the_boundary() {
    let cards = [
        source("hot", [400.0, 700.0], 2.0),
        source("cold", [250.0, 400.0], 3.0),
    ];
    let lower = [temperature_point(400.0), temperature_point(250.0)];
    let upper = [temperature_point(700.0), temperature_point(400.0)];
    let key = PropertyKey::new("k", Dims([1, 1, -3, -1, 0, 0]));
    let queries: [LayerStackSourceQuery<'_>; 2] = std::array::from_fn(|i| LayerStackSourceQuery {
        id: i as u64,
        card: &cards[i],
        conductivity: &key,
        thickness_m: [0.02, 0.015][i],
        temperature_axis: "T",
        lower: &lower[i],
        upper: &upper[i],
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    });
    let stack = PreparedLayerStackConduction::try_new(&queries, &[0.0], 0.004).unwrap();
    assert_eq!(stack.common_temperature_span_k(), Some([400.0, 400.0]));
    let got = stack.evaluate([600.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 80.0, 1e-10);
    assert_eq!(
        got.layer_temperatures_k(),
        &[[600.0, 400.0], [400.0, 300.0]]
    );
}
