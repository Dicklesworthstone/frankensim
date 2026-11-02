//! Physical series laws, source-resolved spatial refinement and native heat transport.
use super::super::material_response::MaterialExchangeLaw;
use super::*;
use crate::conjugate::AirPath;
use crate::graph::tests::{boundary, config, edge, with_cx};
use crate::graph::thermal::medium::FlowHeatGeometry;
use crate::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_convection::{CorrelationId, ThermalDirection};
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Area, Length};
use std::cell::Cell;
use std::rc::Rc;

#[derive(Debug)]
struct Source {
    label: &'static str,
    state: FluidTransportState,
    queries: Rc<Cell<usize>>,
}
impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}
fn source(label: &'static str, rho: f64, mu: f64, k: f64) -> Source {
    Source {
        label,
        state: FluidTransportState::try_new(rho, 2000.0, k, mu).unwrap(),
        queries: Rc::new(Cell::new(0)),
    }
}
fn law(a: f64, b: f64, region: Option<&str>) -> MaterialBranchLaw {
    MaterialBranchLaw {
        resistance: MaterialFlowResistance::try_new(
            a,
            b,
            Area::new(1e-4),
            Length::new(0.01),
            (0.0, 1e7),
        )
        .unwrap(),
        source: SourceProvenance::new("authored analytic section geometry", "section-geometry-v1"),
        uncertainty_rel: 0.0,
        tolerance_basis: ToleranceBasis::Analytic,
        exchange: region
            .into_iter()
            .map(|name| MaterialExchangeLaw {
                region: name.into(),
                geometry: FlowHeatGeometry {
                    flow_area: Area::new(1e-4),
                    characteristic_length: Length::new(0.01),
                    // Fully developed CWT length. Subdivision does not reset it.
                    heated_length: Length::new(1.0),
                    wetted_area: Area::new(0.02),
                    aspect_ratio: None,
                },
                correlation: CorrelationId::CircularDuctLaminarCwt,
                direction: ThermalDirection::HeatingFluid,
            })
            .collect(),
    }
}
fn density() -> BoussinesqDensity {
    BoussinesqDensity::try_new(1000.0, 0.25).unwrap()
}
fn template() -> LossGraph {
    LossGraph::new(2, vec![edge("channel", 0, 1, 1.0)]).unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-8),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_unequal_section_viscosities_control_native_parallel_flow_in_both_directions() {
    with_cx(|cx| {
        let graph = LossGraph::new(
            2,
            vec![edge("channel", 0, 1, 1.0), edge("bypass", 0, 1, 2.0)],
        )
        .unwrap();
        let mut hot = law(9e9, 1e8, Some("hot"));
        hot.uncertainty_rel = 0.05;
        hot.source.identifier = "long-hot-section".into();
        let mut cold = law(1e9, 2e8, Some("cold"));
        cold.uncertainty_rel = 0.2;
        cold.tolerance_basis = ToleranceBasis::ManufacturerDeclared;
        cold.source.identifier = "short-cold-section".into();
        for q in [-1e-5_f64, 1e-5] {
            let result = evaluate_distributed_material_graph(
                cx,
                &graph,
                density(),
                &[vec![hot.clone(), cold.clone()], vec![law(1e10, 0.0, None)]],
                &[q, 2.0 * q],
                vec![
                    vec![
                        source("hot bulk", 900.0, 0.001, 0.5),
                        source("cold bulk", 1100.0, 0.01, 1.0),
                    ],
                    vec![source("bypass bulk", 1000.0, 0.0011, 0.5)],
                ],
                vec![
                    vec![
                        vec![source("hot film", 900.0, 0.001, 0.5)],
                        vec![source("cold film", 1100.0, 0.01, 1.0)],
                    ],
                    vec![vec![]],
                ],
            )
            .unwrap();
            let a = 9e9 * 0.001 + 1e9 * 0.01;
            let b = (1e8 + 2e8) * 1000.0;
            let dp = a * q + b * q * q.abs();
            let solved = result
                .graph
                .solve(&[boundary(0, dp), boundary(1, 0.0)], config(), cx)
                .unwrap();
            near(solved.branches[0].flow.value(), q);
            near(solved.branches[1].flow.value(), 2.0 * q);
            near(solved.node_outflows[0].value(), 3.0 * q);
            near(
                result.sections[0]
                    .iter()
                    .map(|s| s.momentum.evaluation.pressure_drop_pa)
                    .sum(),
                dp,
            );
            near(
                result.sections[0]
                    .iter()
                    .map(|s| s.momentum.evaluation.dissipated_power_w)
                    .sum(),
                dp * q,
            );
            assert_eq!(result.graph.node_count(), graph.node_count());
            for (output, original) in result.graph.branches().iter().zip(graph.branches()) {
                assert_eq!(
                    (output.from, output.to, &output.loss.name),
                    (original.from, original.to, &original.loss.name)
                );
            }
            assert!(result.graph.gravity().is_none());
            assert_eq!(result.sections[0][0].momentum.medium.label, "hot bulk");
            assert_eq!(result.sections[0][1].momentum.medium.label, "cold bulk");
            assert_eq!(result.sections[1][0].momentum.medium.label, "bypass bulk");
            assert_eq!(result.sections[0][0].exchange_range, 0..1);
            assert_eq!(result.sections[0][1].exchange_range, 1..2);
            assert_eq!(result.sections[1][0].exchange_range, 0..0);
            for section in result.sections.iter().flatten() {
                assert_eq!(section.momentum.medium.queries.get(), 1);
                for heat in &section.heat {
                    assert_eq!(heat.medium.queries.get(), 1);
                }
            }
            let linear = result.graph.linear_losses()[0].as_ref().unwrap();
            near(linear.resistance.value(), a);
            near(linear.uncertainty_rel, 0.2);
            near(result.graph.branches()[0].loss.resistance.value(), b);
            near(result.graph.branches()[0].loss.uncertainty_rel, 0.2);
            assert_eq!(linear.tolerance_basis, ToleranceBasis::EngineeringAllowance);
            assert_eq!(
                result.sections[0][1].law.tolerance_basis,
                ToleranceBasis::ManufacturerDeclared
            );
            assert!(linear.source.identifier.contains("long-hot-section"));
            assert!(linear.source.identifier.contains("short-cold-section"));
            let BranchThermalModel::Exchange(segments) = &result.models[0] else {
                panic!("missing exchange")
            };
            assert_eq!(
                segments.iter().map(|s| s.region()).collect::<Vec<_>>(),
                vec!["hot", "cold"]
            );
            near(segments[0].htc_w_per_m2_k(), 3.66 * 0.5 / 0.01);
            near(segments[1].htc_w_per_m2_k(), 3.66 * 1.0 / 0.01);
            assert_eq!(result.models[1], BranchThermalModel::Adiabatic);
        }
    });
}

#[test]
fn g1_physical_subdivision_preserves_drag_and_native_exponential_wall_heat() {
    with_cx(|cx| {
        let graph = template();
        let h = 3.66 * 0.5 / 0.01;
        let expected_outlet = 350.0 - 50.0 * (-h * 0.04 / 20.0_f64).exp();
        for count in [1_usize, 2, 4] {
            let laws = (0..count)
                .map(|i| {
                    let mut local = law(1e10 / count as f64, 0.0, Some(&format!("wall-{i}")));
                    local.exchange[0].geometry.wetted_area = Area::new(0.04 / count as f64);
                    local
                })
                .collect::<Vec<_>>();
            let response = evaluate_distributed_material_graph(
                cx,
                &graph,
                density(),
                &[laws],
                &[1e-5],
                vec![
                    (0..count)
                        .map(|_| source("bulk", 1000.0, 0.002, 0.5))
                        .collect(),
                ],
                vec![
                    (0..count)
                        .map(|_| vec![source("film", 1000.0, 0.002, 0.5)])
                        .collect(),
                ],
            )
            .unwrap();
            let solved = response
                .graph
                .solve(&[boundary(0, 200.0), boundary(1, 0.0)], config(), cx)
                .unwrap();
            near(solved.branches[0].flow.value(), 1e-5);
            let BranchThermalModel::Exchange(segments) = &response.models[0] else {
                panic!("missing exchange")
            };
            let path = AirPath::new(300.0, 0.01, 2000.0, segments.clone()).unwrap();
            let heat = path.march(&vec![350.0; count]).unwrap();
            near(heat.outlet_temperature_k, expected_outlet);
            near(heat.total_heat_rate_w, 20.0 * (expected_outlet - 300.0));
            for section in &response.sections[0] {
                near(section.heat[0].geometry.heated_length.value(), 1.0);
            }
        }
    });
}

#[test]
fn g3_spatial_material_sampling_converges_to_independent_integrated_viscous_drag() {
    with_cx(|cx| {
        // Exercise a pure-linear input template as well as the mixed output
        // case above. No artificial quadratic placeholder is needed.
        let graph = LossGraph::new_mixed(
            2,
            vec![PassiveBranch::Linear {
                from: 0,
                to: 1,
                loss: LinearLossElement::new(
                    "channel",
                    LinearLossResistance::new(1.0),
                    0.0,
                    SourceProvenance::new("analytic template", "template"),
                    ToleranceBasis::Analytic,
                )
                .unwrap(),
            }],
        )
        .unwrap();
        let mut previous_error = None;
        let mut previous_flow = f64::INFINITY;
        for count in [1_usize, 2, 4, 8] {
            let materials = (0..count)
                .map(|i| {
                    let x = (i as f64 + 0.5) / count as f64;
                    source("spatial bulk", 1000.0, 0.001 + 0.003 * x * x, 0.5)
                })
                .collect();
            let response = evaluate_distributed_material_graph(
                cx,
                &graph,
                density(),
                &[vec![law(1e10 / count as f64, 0.0, None); count]],
                &[1e-5],
                vec![materials],
                vec![(0..count).map(|_| vec![]).collect()],
            )
            .unwrap();
            let a = response.graph.linear_losses()[0]
                .as_ref()
                .unwrap()
                .resistance
                .value();
            // Integral_0^1 (0.001 + 0.003*x^2) dx = 0.002. The independent
            // midpoint error is 0.00025/N^2, so resistance error quarters.
            near(a, 2e7 - 2.5e6 / (count * count) as f64);
            let error = 2e7 - a;
            if let Some(previous) = previous_error {
                near(previous / error, 4.0);
            }
            previous_error = Some(error);
            let solved = response
                .graph
                .solve(&[boundary(0, 200.0), boundary(1, 0.0)], config(), cx)
                .unwrap();
            let flow = solved.branches[0].flow.value();
            near(flow, 200.0 / a);
            assert!(flow > 1e-5 && flow < previous_flow);
            previous_flow = flow;
            assert_eq!(response.models[0], BranchThermalModel::Adiabatic);
        }
    });
}

#[test]
fn g0_local_support_counts_and_overflow_refuse_without_partial_response() {
    with_cx(|cx| {
        let graph = template();
        let mut unsupported = law(1e10, 0.0, Some("wall"));
        unsupported.exchange[0].correlation = CorrelationId::CircularDuctLaminarChf;
        let result = evaluate_distributed_material_graph(
            cx,
            &graph,
            density(),
            &[vec![law(1e10, 0.0, None), unsupported]],
            &[1e-5],
            vec![vec![
                source("first", 1000.0, 0.002, 0.5),
                source("second", 1000.0, 0.002, 0.5),
            ]],
            vec![vec![vec![], vec![source("film", 1000.0, 0.002, 0.5)]]],
        );
        assert!(matches!(
            result,
            Err(DistributedResponseError::Section {
                branch: 0,
                section: 1,
                source: MaterialResponseError::Branch { .. },
            })
        ));
        let result = evaluate_distributed_material_graph(
            cx,
            &graph,
            density(),
            &[vec![law(1e10, 0.0, None); 2]],
            &[1e-5],
            vec![vec![
                source("valid", 1000.0, 0.002, 0.5),
                source("invalid density", 500.0, 0.002, 0.5),
            ]],
            vec![vec![vec![], vec![]]],
        );
        assert!(matches!(
            result,
            Err(DistributedResponseError::Section {
                branch: 0,
                section: 1,
                ..
            })
        ));
        let original = source("unqueried", 1000.0, 0.002, 0.5);
        let queries = original.queries.clone();
        assert!(
            evaluate_distributed_material_graph(
                cx,
                &graph,
                density(),
                &[vec![law(1e10, 0.0, None)]],
                &[1e-5],
                vec![vec![original]],
                vec![vec![]],
            )
            .is_err()
        );
        assert_eq!(queries.get(), 0);

        let count = MAX_MATERIAL_FLOW_SECTIONS + 1;
        let reads = Rc::new(Cell::new(0));
        let mut bulk = Vec::new();
        for _ in 0..count {
            let mut material = source("over budget", 1000.0, 0.002, 0.5);
            material.queries = reads.clone();
            bulk.push(material);
        }
        assert!(matches!(
            evaluate_distributed_material_graph(
                cx,
                &graph,
                density(),
                &[vec![law(1e10, 0.0, None); count]],
                &[1e-5],
                vec![bulk],
                vec![(0..count).map(|_| vec![]).collect()],
            ),
            Err(DistributedResponseError::Invalid(_))
        ));
        assert_eq!(reads.get(), 0);

        assert!(matches!(
            evaluate_distributed_material_graph(
                cx,
                &graph,
                density(),
                &[vec![law(1e308, 0.0, None); 2]],
                &[1e-100],
                vec![vec![
                    source("huge a", 1000.0, 1.0, 0.5),
                    source("huge b", 1000.0, 1.0, 0.5)
                ]],
                vec![vec![vec![], vec![]]],
            ),
            Err(DistributedResponseError::Branch { branch: 0, .. })
        ));
    });
}
