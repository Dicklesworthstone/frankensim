//! Original material, pressure graph and correlation owners against independent laws.
use super::*;
use crate::graph::tests::{boundary, config, edge, with_cx};
use fs_material::gas::{GasSpec, GasState};
use fs_qty::{Area, Length};

#[derive(Debug)]
struct Source {
    label: &'static str,
    state: FluidTransportState,
    queries: std::cell::Cell<usize>,
}
impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}
fn source(label: &'static str, rho: f64, mu: f64, conductivity: f64) -> Source {
    Source {
        label,
        state: FluidTransportState::try_new(rho, 2000.0, conductivity, mu).unwrap(),
        queries: std::cell::Cell::new(0),
    }
}
fn geometry() -> FlowHeatGeometry {
    FlowHeatGeometry {
        flow_area: Area::new(1e-4),
        characteristic_length: Length::new(0.01),
        heated_length: Length::new(1.0),
        wetted_area: Area::new(0.02),
        aspect_ratio: None,
    }
}
fn law(a: f64, b: f64) -> MaterialBranchLaw {
    MaterialBranchLaw {
        resistance: MaterialFlowResistance::try_new(
            a,
            b,
            Area::new(1e-4),
            Length::new(0.01),
            (0.0, 1e7),
        )
        .unwrap(),
        source: SourceProvenance::new(
            "explicit analytic geometry fixture",
            "material-response-test-v1",
        ),
        uncertainty_rel: 0.0,
        tolerance_basis: ToleranceBasis::Analytic,
        exchange: vec![MaterialExchangeLaw {
            region: "wall".into(),
            geometry: geometry(),
            correlation: CorrelationId::CircularDuctLaminarCwt,
            direction: ThermalDirection::HeatingFluid,
        }],
    }
}
fn template() -> LossGraph {
    LossGraph::new(2, vec![edge("native", 0, 1, 1.0)]).unwrap()
}
fn density() -> BoussinesqDensity {
    BoussinesqDensity::try_new(1000.0, 0.25).unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_rebuilt_native_drag_matches_linear_quadratic_and_signed_material_laws() {
    with_cx(|cx| {
        let graph = template();
        for (a, b) in [(1e10, 0.0), (0.0, 2e8), (1e10, 2e8)] {
            for q in [-1e-5_f64, 1e-5] {
                let result = evaluate_material_graph(
                    cx,
                    &graph,
                    density(),
                    &[law(a, b)],
                    &[q],
                    vec![source("bulk", 900.0, 0.002, 0.5)],
                    vec![vec![source("film", 900.0, 0.002, 0.5)]],
                )
                .unwrap();
                let expected_dp = a * 0.002 * q + b * 1000.0 * q * q.abs();
                let physical = result
                    .graph
                    .solve(&[boundary(0, expected_dp), boundary(1, 0.0)], config(), cx)
                    .unwrap();
                near(physical.branches[0].flow.value(), q);
                near(result.momentum[0].evaluation.pressure_drop_pa, expected_dp);
                near(
                    result.momentum[0].evaluation.dissipated_power_w,
                    expected_dp * q,
                );
                assert_eq!(
                    result.momentum[0].evaluation.transport.density_kg_m3(),
                    900.0
                );
                assert_eq!(result.momentum[0].evaluation.momentum_density_kg_m3, 1000.0);
                assert_eq!(result.momentum[0].medium.label, "bulk");
                assert_eq!(result.heat[0][0].medium.label, "film");
                assert_eq!(result.momentum[0].medium.queries.get(), 1);
                assert_eq!(result.heat[0][0].medium.queries.get(), 1);
                assert_eq!(result.graph.branches()[0].from, 0);
                assert_eq!(result.graph.branches()[0].to, 1);
                assert_eq!(result.graph.branches()[0].loss.name, "native");
                assert!(result.graph.gravity().is_none());
                near(
                    result.heat[0][0].coefficient.value.value(),
                    3.66 * 0.5 / 0.01,
                );
            }
        }
    });
}

#[test]
fn g3_temperature_updated_viscosity_and_conductivity_change_real_flow_and_heat() {
    with_cx(|cx| {
        let graph = template();
        let run = |mu, k, q| {
            evaluate_material_graph(
                cx,
                &graph,
                density(),
                &[law(1e10, 0.0)],
                &[q],
                vec![source("bulk", 1000.0, mu, k)],
                vec![vec![source("film", 1000.0, mu, k)]],
            )
            .unwrap()
        };
        let cold = run(0.002, 0.5, 1e-5);
        let warm = run(0.001, 1.0, 2e-5);
        let solve = |response: &MaterialGraphResponse<Source>| {
            response
                .graph
                .solve(&[boundary(0, 200.0), boundary(1, 0.0)], config(), cx)
                .unwrap()
        };
        near(
            solve(&warm).branches[0].flow.value(),
            2.0 * solve(&cold).branches[0].flow.value(),
        );
        near(
            warm.heat[0][0].coefficient.value.value(),
            2.0 * cold.heat[0][0].coefficient.value.value(),
        );
        assert_eq!(cold.models.len(), 1);
        if let BranchThermalModel::Exchange(segments) = &warm.models[0] {
            near(
                segments[0].htc_w_per_m2_k(),
                warm.heat[0][0].coefficient.value.value(),
            );
            near(segments[0].area_m2(), geometry().wetted_area.value());
        } else {
            panic!("explicit exchange disappeared");
        }
    });
}

#[test]
fn g1_actual_gas_and_solved_speed_drive_hausen_without_detaching_materials() {
    with_cx(|cx| {
        let gas = || GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, 101325.0).unwrap();
        let rho = gas().transport_state().unwrap().density_kg_m3();
        let density = BoussinesqDensity::try_new(rho, 0.1).unwrap();
        let mut branch = law(1e10, 0.0);
        branch.exchange[0].correlation = CorrelationId::CircularDuctHausen;
        branch.exchange[0].geometry.heated_length = Length::new(0.1);
        let run = |q| {
            evaluate_material_graph(
                cx,
                &template(),
                density,
                &[branch.clone()],
                &[q],
                vec![gas()],
                vec![vec![gas()]],
            )
            .unwrap()
        };
        let low = run(1e-4);
        let high = run(2e-4);
        let reverse = run(-2e-4);
        assert!(
            high.heat[0][0].coefficient.value.value() > low.heat[0][0].coefficient.value.value()
        );
        near(
            high.heat[0][0].coefficient.value.value(),
            reverse.heat[0][0].coefficient.value.value(),
        );
        let state = high.heat[0][0].medium.transport_state().unwrap();
        let re = state.density_kg_m3() * (2e-4 / 1e-4) * 0.01 / state.dynamic_viscosity_pa_s();
        let pr = state.specific_heat_j_kg_k() * state.dynamic_viscosity_pa_s()
            / state.thermal_conductivity_w_m_k();
        let gz = re * pr / 10.0;
        let nu = 3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0));
        near(
            high.heat[0][0].coefficient.value.value(),
            nu * state.thermal_conductivity_w_m_k() / 0.01,
        );
        assert_eq!(high.heat[0][0].region, reverse.heat[0][0].region);
        assert_eq!(high.heat[0][0].geometry.wetted_area, geometry().wetted_area);
    });
}

#[test]
fn g0_source_density_reynolds_and_boundary_condition_refusals_remain_terminal() {
    with_cx(|cx| {
        let graph = template();
        let run = |branch: MaterialBranchLaw, q, rho| {
            evaluate_material_graph(
                cx,
                &graph,
                density(),
                &[branch],
                &[q],
                vec![source("bulk", rho, 0.002, 0.5)],
                vec![vec![source("film", rho, 0.002, 0.5)]],
            )
        };
        let mut wrong_wall = law(1e10, 0.0);
        wrong_wall.exchange[0].correlation = CorrelationId::CircularDuctLaminarChf;
        assert!(run(wrong_wall, 1e-5, 1000.0).is_err());
        assert!(matches!(
            run(law(1e10, 0.0), 1e-5, 500.0),
            Err(MaterialResponseError::Branch { .. })
        ));
        assert!(matches!(
            run(law(1e10, 0.0), 1e-3, 1000.0),
            Err(MaterialResponseError::Convection { .. })
        ));
        let mut unsupported = law(1e10, 0.0);
        unsupported.resistance = MaterialFlowResistance::try_new(
            1e10,
            0.0,
            Area::new(1e-4),
            Length::new(0.01),
            (0.0, 100.0),
        )
        .unwrap();
        assert!(matches!(
            run(unsupported, 1e-5, 1000.0),
            Err(MaterialResponseError::Resistance { .. })
        ));
        assert!(run(law(1e10, 0.0), 0.0, 1000.0).is_err());
        let mut bad_source = law(1e10, 0.0);
        bad_source.source.identifier.clear();
        assert!(run(bad_source, 1e-5, 1000.0).is_err());
        assert!(
            evaluate_material_graph::<Source>(cx, &graph, density(), &[], &[], vec![], vec![])
                .is_err()
        );
    });
}


fn turbulent_law(length_m: f64, roughness_m: f64) -> MaterialBranchLaw {
    let mut result = law(0.0, 1.0);
    result.resistance = MaterialFlowResistance::circular_turbulent(
        Length::new(length_m), Length::new(0.01), Length::new(roughness_m), (4000.0, 1e6),
    ).unwrap();
    result.exchange.clear();
    result
}

// Independent algebraic inversion of the original Colebrook equation using
// platform sqrt/log10, with explicit geometry rather than production helpers.
fn colebrook_flow(dp: f64, length: f64, roughness: f64, mu: f64) -> f64 {
    let diameter = 0.01;
    let area = core::f64::consts::FRAC_PI_4 * diameter * diameter;
    let t = area * (2.0 * diameter * dp.abs() / (length * 1000.0)).sqrt();
    let q = -2.0 * t * (roughness / (3.7 * diameter)
        + 2.51 * mu * area / (1000.0 * diameter * t)).log10();
    q.copysign(dp)
}

fn turbulent_response(cx: &Cx<'_>, length: f64, mu: f64) -> MaterialGraphResponse<Source> {
    evaluate_material_graph(cx, &template(), density(), &[turbulent_law(length, 1e-5)],
        &[1e-4], vec![source("current bulk", 900.0, mu, 0.5)], vec![vec![]]).unwrap()
}

#[test]
fn turbulent_junction_preserves_material_laws_reverse_flow_and_pressure_work() {
    with_cx(|cx| {
        let topology = LossGraph::new(3, vec![
            edge("feed", 0, 1, 1.0), edge("smooth outlet", 1, 2, 1.0),
            edge("rough reverse outlet", 2, 1, 1.0),
        ]).unwrap();
        let pressure = 5000.0;
        let forward = colebrook_flow(pressure, 1.0, 1e-5, 0.001);
        let reverse = colebrook_flow(pressure, 2.5, 5e-5, 0.001);
        let feed_resistance = 3.7e7;
        let mut feed = law(feed_resistance / 0.001, 0.0);
        feed.exchange.clear();
        let laws = [feed, turbulent_law(1.0, 1e-5), turbulent_law(2.5, 5e-5)];
        let expected = [forward + reverse, forward, -reverse];
        let response = evaluate_material_graph(cx, &topology, density(), &laws, &expected,
            (0..3).map(|_| source("bulk", 900.0, 0.001, 0.5)).collect(),
            vec![vec![], vec![], vec![]]).unwrap();
        let inlet = pressure + feed_resistance * expected[0];
        let mut controls = config();
        controls.absolute_flow_tolerance = fs_qty::VolumetricFlowRate::new(1e-16);
        controls.relative_flow_tolerance = 1e-12;
        let solved = response.graph.solve(&[boundary(0, inlet), boundary(2, 0.0)], controls, cx).unwrap();
        assert!((solved.pressures[1].value() / pressure - 1.0).abs() < 2e-10);
        for (index, (branch, &q)) in solved.branches.iter().zip(&expected).enumerate() {
            assert!((branch.flow.value() / q - 1.0).abs() < 2e-10);
            assert_eq!(response.momentum[index].medium.queries.get(), 1);
            assert_eq!(response.momentum[index].evaluation.transport.density_kg_m3(), 900.0);
            if index > 0 {
                let loss = branch.turbulent_loss.as_ref().unwrap();
                assert_eq!(loss.law(), laws[index].resistance);
                assert_eq!(loss.source(), &laws[index].source);
                let actual = loss.evaluate(branch.flow).unwrap();
                near(actual.pressure_drop_pa, branch.dissipative_pressure_drop.value());
                near(actual.dissipated_power_w,
                    branch.dissipative_pressure_drop.value() * branch.flow.value());
                assert!(actual.dissipated_power_w > 0.0);
                assert!(actual.dpressure_dmass_flow > 0.0);
                assert_eq!(branch.loss.resistance.value(), 0.0);
                assert!(branch.linear_loss.is_none());
            }
        }
        let dissipated: f64 = solved.branches.iter().map(|b|
            b.dissipative_pressure_drop.value() * b.flow.value()).sum();
        assert!((dissipated / (inlet * solved.node_outflows[0].value()) - 1.0).abs() < 2e-10);
        assert!(solved.node_outflows[1].value().abs() < 1e-15);
    });
}

#[test]
fn turbulent_exact_graph_tangent_and_adjoint_follow_length_and_pressure_perturbations() {
    use crate::graph::sensitivity::{FlowWeight, GraphSensitivityConfig, PressureWeight};
    use fs_qty::{Pressure, VolumetricFlowRate};
    use fs_solver::CheckedCgConfig;
    with_cx(|cx| {
        let length = 1.25;
        let pressure = 4000.0;
        let response = turbulent_response(cx, length, 0.001);
        let linear = response.graph.solve_linearized(
            &[boundary(0, pressure), boundary(1, 0.0)], config(),
            GraphSensitivityConfig { minimum_absolute_flow: VolumetricFlowRate::new(0.0),
                linear: CheckedCgConfig { tolerance: 1e-12, max_iterations: 64, max_corrections: 2 } },
            cx).unwrap();
        let tangent = linear.tangent(&[0.4], &[Pressure::new(120.0), Pressure::new(-30.0)], cx).unwrap();
        let gradient = linear.adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], cx).unwrap();
        let step = 1e-5;
        let expected = (colebrook_flow(pressure + step * 150.0, length * (1.0 + step * 0.4), 1e-5, 0.001)
            - colebrook_flow(pressure - step * 150.0, length * (1.0 - step * 0.4), 1e-5, 0.001))
            / (2.0 * step);
        assert!((tangent.flows[0].value() / expected - 1.0).abs() < 1e-8);
        let length_derivative = (colebrook_flow(pressure, length * (1.0 + step), 1e-5, 0.001)
            - colebrook_flow(pressure, length * (1.0 - step), 1e-5, 0.001)) / (2.0 * step);
        assert!((gradient.log_resistances[0] / length_derivative - 1.0).abs() < 1e-8);
        let dual = 0.4 * gradient.log_resistances[0]
            + 120.0 * gradient.boundary_pressures[0].value()
            - 30.0 * gradient.boundary_pressures[1].value();
        assert!((dual / tangent.flows[0].value() - 1.0).abs() < 1e-12);
        let q = linear.solution().branches[0].flow.value();
        // Freezing f would incorrectly recover the quadratic q/2 response.
        assert!((gradient.log_resistances[0] + 0.5 * q).abs() > 0.01 * q);
        let viscous = turbulent_response(cx, length, 0.002).graph
            .solve(&[boundary(0, pressure), boundary(1, 0.0)], config(), cx).unwrap();
        let expected_viscous = colebrook_flow(pressure, length, 1e-5, 0.002);
        assert!((viscous.branches[0].flow.value() / expected_viscous - 1.0).abs() < 1e-12);
        assert!(viscous.branches[0].flow.value() < q);
    });
}

#[test]
fn turbulent_fan_intersection_uses_actual_nonlinear_curve_and_private_zero_flow_brackets() {
    use crate::graph::GraphFanConfig;
    use crate::{FanArrangement, FanBank, FanCurve, FanPoint};
    use fs_qty::{Pressure, VolumetricFlowRate};
    with_cx(|cx| {
        let pressure = 6000.0;
        let q = colebrook_flow(pressure, 1.0, 1e-5, 0.001);
        let fan = FanBank::new(FanCurve::new("analytic crossing",
            vec![FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(2.0 * pressure)),
                FanPoint::new(VolumetricFlowRate::new(2.0 * q), Pressure::new(0.0))],
            SourceProvenance::new("authored fan fixture", "colebrook-crossing"),
            0.0, ToleranceBasis::Analytic, VolumetricFlowRate::new(0.0), (0.5, 2.0),
        ).unwrap(), 1, FanArrangement::Series, 1.0).unwrap();
        let response = turbulent_response(cx, 1.0, 0.001);
        let point = response.graph.solve_with_fan(0, 1, &fan, config(),
            GraphFanConfig { max_iterations: 128,
                absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
                absolute_pressure_tolerance: Pressure::new(1e-7) }, cx).unwrap();
        assert!((point.network.branches[0].flow.value() / q - 1.0).abs() < 1e-9);
        assert!((point.network.pressures[0].value() / pressure - 1.0).abs() < 1e-9);
        assert!(point.network.branches[0].turbulent_loss.is_some());
    });
}

#[test]
fn turbulent_continuation_refuses_unsupported_physics_and_native_inertia_retains_the_law() {
    use crate::graph::inertial::{InertialBranch, InertialGraph};
    use crate::graph::thermal::natural_circulation::distributed_response::{
        DistributedResponseError, evaluate_distributed_material_graph,
    };
    use fs_qty::{Density, VolumetricFlowRate};
    with_cx(|cx| {
        let response = turbulent_response(cx, 1.0, 0.001);
        let pipe = response.graph.turbulent_losses()[0].as_ref().unwrap();
        for re in [0.0, 2000.0, 2e6] {
            // Endpoint losses bound the original supported pressure range.
            let endpoint = if re < 4000.0 { 5000.0 } else { 9e5 };
            let mass = endpoint * 0.001 * pipe.law().flow_area().value() / 0.01;
            let admitted = pipe.law().evaluate(&source("reference", 1000.0, 0.001, 0.5),
                MassFlowRate::new(mass)).unwrap();
            let dp = if re == 0.0 { 0.0 } else if re < 4000.0 {
                admitted.pressure_drop_pa * 0.01
            } else {
                admitted.pressure_drop_pa * 100.0
            };
            assert!(matches!(response.graph.solve(&[boundary(0, dp), boundary(1, 0.0)], config(), cx),
                Err(GraphError::Resistance { branch: 0, source: ResistanceError::Reynolds { .. } })));
        }
        for (re, lower) in [(4000.0, true), (1e6, false)] {
            let mut q = re * 0.001 * pipe.law().flow_area().value() / (1000.0 * 0.01);
            let mut supported = None;
            for _ in 0..16 {
                if let Ok(value) = pipe.evaluate(VolumetricFlowRate::new(q)) {
                    supported = Some(value);
                    break;
                }
                q = if lower { q.next_up() } else { q.next_down() };
            }
            let dp = supported.unwrap().pressure_drop_pa;
            for sign in [-1.0, 1.0] {
                let value = response.graph.solve(
                    &[boundary(0, sign * dp), boundary(1, 0.0)], config(), cx).unwrap();
                assert!((value.branches[0].flow.value() / (sign * q) - 1.0).abs() < 1e-12);
            }
        }
        assert!(pipe.evaluate(VolumetricFlowRate::new(0.0)).is_err());
        let inertial = InertialGraph::try_new(cx, &response.graph, Density::new(1000.0),
            &[InertialBranch { length_m: 1.0, flow_area_m2: pipe.law().flow_area().value() }],
            &[boundary(0, 6000.0), boundary(1, 0.0)], &[1e-4], 0.0, 8).unwrap();
        assert_eq!(inertial.volume_flows_m3_s(), vec![1e-4]);
        assert!(matches!(evaluate_distributed_material_graph(cx, &template(), density(),
            &[vec![turbulent_law(1.0, 1e-5)]], &[1e-4],
            vec![vec![source("bulk", 900.0, 0.001, 0.5)]], vec![vec![vec![]]]),
            Err(DistributedResponseError::Branch { branch: 0, .. })));
        assert!(response.graph.clone().with_gravity(
            &source("wrong reference", 900.0, 0.001, 0.5), &[[0.0; 3], [0.0, 0.0, 1.0]],
            [0.0, 0.0, -9.81], cx).is_err());
    });
}

#[test]
fn churchill_graph_crosses_transition_and_reverses_with_original_materials() {
    with_cx(|cx| {
        let mut selected = law(1.0, 0.0);
        selected.resistance = MaterialFlowResistance::circular_churchill(
            Length::new(2.0), Length::new(0.01), Length::new(1e-6), (0.0, 1e6),
        ).unwrap();
        selected.exchange.clear();
        let area = core::f64::consts::PI * 0.01_f64.powi(2) / 4.0;
        let q = 3000.0 * 0.001 * area / (1000.0 * 0.01);
        for (mu, re) in [(0.003, 1000.0_f64), (0.001, 3000.0), (0.0003, 10000.0)] {
            // Direct Churchill equation with the bounded Boussinesq momentum
            // density, while the original 900 kg/m^3 snapshot stays attached.
            let a = (-2.457 * ((7.0 / re).powf(0.9) + 0.27e-4).ln()).powi(16);
            let b = (37530.0 / re).powi(16);
            let f = 8.0 * ((8.0 / re).powi(12) + (a + b).powf(-1.5)).powf(1.0 / 12.0);
            let dp = f * 200.0 * 500.0 * (q / area).powi(2);
            for sign in [-1.0, 1.0] {
                let response = evaluate_material_graph(cx, &template(), density(),
                    &[selected.clone()], &[sign * q],
                    vec![source("current bulk", 900.0, mu, 0.5)], vec![vec![]]).unwrap();
                let mut controls = config();
                controls.absolute_flow_tolerance = fs_qty::VolumetricFlowRate::new(1e-16);
                controls.relative_flow_tolerance = 1e-12;
                let solved = response.graph.solve(
                    &[boundary(0, sign * dp), boundary(1, 0.0)], controls, cx,
                ).unwrap();
                let flow = solved.branches[0].flow.value();
                assert!((flow / (sign * q) - 1.0).abs() < 3e-12);
                let original = solved.branches[0].turbulent_loss.as_ref().unwrap();
                let value = original.evaluate(solved.branches[0].flow).unwrap();
                assert_eq!(original.law(), selected.resistance);
                assert_eq!(original.source(), &selected.source);
                assert_eq!(value.transport.density_kg_m3(), 900.0);
                assert_eq!(value.transport.dynamic_viscosity_pa_s(), mu);
                assert_eq!(value.momentum_density_kg_m3, 1000.0);
                assert!((value.reynolds / re - 1.0).abs() < 3e-12);
                assert!((value.dissipated_power_w / (dp * q) - 1.0).abs() < 6e-12);
                assert_eq!(response.momentum[0].medium.queries.get(), 1);
            }
        }
    });
}

#[test]
fn churchill_graph_rest_has_poiseuille_tangent_and_continuous_signed_pressure_response() {
    use crate::graph::sensitivity::GraphSensitivityConfig;
    use fs_qty::{Pressure, VolumetricFlowRate};
    use fs_solver::CheckedCgConfig;
    with_cx(|cx| {
        let mut selected = law(1.0, 0.0);
        selected.resistance = MaterialFlowResistance::circular_churchill(
            Length::new(2.0), Length::new(0.01), Length::new(1e-6), (0.0, 1e6),
        ).unwrap();
        selected.exchange.clear();
        let response = evaluate_material_graph(cx, &template(), density(), &[selected],
            &[0.0], vec![source("rest", 900.0, 0.001, 0.5)], vec![vec![]]).unwrap();
        let linear = response.graph.solve_linearized(
            &[boundary(0, 0.0), boundary(1, 0.0)], config(),
            GraphSensitivityConfig { minimum_absolute_flow: VolumetricFlowRate::new(0.0),
                linear: CheckedCgConfig { tolerance: 1e-12, max_iterations: 64, max_corrections: 2 } },
            cx,
        ).unwrap();
        assert_eq!(linear.solution().branches[0].flow.value(), 0.0);
        let original = linear.solution().branches[0].turbulent_loss.as_ref().unwrap();
        assert_eq!(original.evaluate(VolumetricFlowRate::new(0.0)).unwrap().dissipated_power_w, 0.0);
        let tangent = linear.tangent(&[0.0], &[Pressure::new(1.0), Pressure::new(0.0)], cx).unwrap();
        let conductance = core::f64::consts::PI * 0.01_f64.powi(4) / (128.0 * 0.001 * 2.0);
        assert!((tangent.flows[0].value() / conductance - 1.0).abs() < 1e-14);
        for dp in [-1e-8, 1e-8] {
            let solved = response.graph.solve(
                &[boundary(0, dp), boundary(1, 0.0)], config(), cx,
            ).unwrap();
            let q = solved.branches[0].flow.value();
            assert!((q / (conductance * dp) - 1.0).abs() < 1e-13);
            assert!(dp * q > 0.0);
        }
        assert_eq!(response.momentum[0].medium.queries.get(), 1);
    });
}
