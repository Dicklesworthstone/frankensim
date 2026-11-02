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
