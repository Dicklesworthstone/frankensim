//! Actual native graph/material/heat owners against an independent ideal-gas
//! thermosiphon root and exact density integrals along exponential profiles.

use super::*;
use crate::conjugate::AirSegment;
use crate::graph::tests::{boundary, config as hydraulic_config, edge, with_cx};
use fs_material::gas::{GasSpec, GasState};

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101325.0).unwrap()
}

fn controls() -> NaturalCirculationConfig {
    let mut hydraulics = hydraulic_config();
    hydraulics.absolute_flow_tolerance = VolumetricFlowRate::new(1e-14);
    hydraulics.relative_flow_tolerance = 1e-13;
    NaturalCirculationConfig {
        hydraulics,
        transport: TransportConfig {
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
            relative_flow_tolerance: 1e-11,
            absolute_heat_tolerance_w: 1e-6,
            relative_heat_tolerance: 1e-9,
        },
        max_iterations: 100,
        max_material_evaluations: 200_000,
        max_density_quadrature_panels: 256,
        absolute_density_tolerance_kg_m3: 1e-9,
        relative_density_tolerance: 1e-9,
        absolute_flow_change_tolerance: VolumetricFlowRate::new(1e-11),
        relative_flow_change_tolerance: 1e-8,
        relaxation: 0.5,
        acceleration: Some(IqnIlsConfig::default()),
        maximum_neglected_dissipation_w: 1.0,
    }
}

fn exchange(name: &str, conductance: f64, wall: f64) -> NaturalCirculationBranch {
    NaturalCirculationBranch {
        model: BranchThermalModel::Exchange(vec![AirSegment::new(name, 1.0, conductance).unwrap()]),
        axial_fractions: vec![1.0],
        wall_temperatures_k: vec![wall],
    }
}

struct Fixture {
    graph: LossGraph,
    branches: Vec<NaturalCirculationBranch>,
    positions: [[f64; 3]; 2],
    datums: [FixedPressure; 1],
}
impl Fixture {
    fn new(reverse_return: bool, hot_wall: f64) -> Self {
        Self {
            graph: LossGraph::new(2, vec![edge("rise", 0, 1, 10_000.0),
                if reverse_return { edge("return", 0, 1, 10_000.0) }
                else { edge("return", 1, 0, 10_000.0) }]).unwrap(),
            branches: vec![exchange("hot", 50.0, hot_wall), exchange("cold", 80.0, 300.0)],
            positions: [[0.0; 3], [0.0, 0.0, 1.0]],
            datums: [boundary(0, 0.0)],
        }
    }
    fn problem(&self) -> NaturalCirculationProblem<'_> {
        NaturalCirculationProblem {
            graph: &self.graph,
            reference: gas(325.0).transport_state().unwrap(),
            node_positions_m: &self.positions,
            acceleration_m_s2: [0.0, 0.0, -9.81],
            maximum_relative_density_departure: 0.25,
            pressure_datums: &self.datums,
            branches: &self.branches,
        }
    }
    fn solve(&self, cx: &Cx<'_>, config: NaturalCirculationConfig)
        -> Result<NaturalCirculationSolution, NaturalCirculationError> {
        solve_natural_circulation(cx, &self.problem(), config, &[335.0, 310.0],
            |_, _, temperature| GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature.value(), 101325.0))
    }
}

// For T(x)=Tw+(Tin-Tw)exp(-N x), integral_0^1 dx/T(x)
// = (1 + log(Tout/Tin)/N)/Tw. This oracle never calls production quadrature.
fn mean_gas_density(inlet: f64, outlet: f64, wall: f64, ntu: f64) -> f64 {
    let rho_temperature = gas(325.0).density * 325.0;
    rho_temperature/wall * (1.0 + (outlet/inlet).ln()/ntu)
}

fn analytic_loop_flow(hot_wall: f64) -> f64 {
    let reference = gas(325.0);
    let residual = |q: f64| {
        let capacity = reference.density*reference.specific_heat_cp*q;
        let hot_ntu = 50.0/capacity;
        let cold_ntu = 80.0/capacity;
        let (a, b) = ((-hot_ntu).exp(), (-cold_ntu).exp());
        let cold = ((1.0-b)*300.0 + b*(1.0-a)*hot_wall)/(1.0-a*b);
        let hot = hot_wall+(cold-hot_wall)*a;
        let rho_hot = mean_gas_density(cold, hot, hot_wall, hot_ntu);
        let rho_cold = mean_gas_density(hot, cold, 300.0, cold_ntu);
        20_000.0*q*q - 9.81*(rho_cold-rho_hot)
    };
    let (mut lo, mut hi) = (1e-5, 0.05);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let middle = 0.5*(lo+hi);
        if residual(middle) > 0.0 { hi = middle; } else { lo = middle; }
    }
    0.5*(lo+hi)
}

#[test]
fn g1_temperature_feedback_selects_the_analytic_ideal_gas_circulation() {
    with_cx(|cx| {
        let fixture = Fixture::new(false, 350.0);
        let result = fixture.solve(cx, controls()).unwrap();
        let expected = analytic_loop_flow(350.0);
        for branch in &result.hydraulics.branches {
            assert!((branch.flow.value()-expected).abs() < 2e-9,
                "flow {}, analytic {expected}", branch.flow.value());
        }
        assert!(result.iterations > 1 && result.hydraulic_solves > result.iterations);
        assert!(result.max_density_residual_kg_m3 < 3e-9);
        assert!(result.max_confirmed_flow_change_m3_s < 2e-10);
        assert!(result.transport.heat_imbalance_w.abs() < 1e-6);
        assert_eq!(result.transport.external_heat_gain_w, 0.0);
        assert!(result.transport.branches[0].march.as_ref().unwrap().total_heat_rate_w > 0.0);
        assert!(result.neglected_dissipation_w > 0.0 && result.neglected_dissipation_w < 1.0);
        // The same source law and real wall change must change the operating point.
        let hotter = Fixture::new(false, 370.0).solve(cx, controls()).unwrap();
        assert!(hotter.hydraulics.branches[0].flow.value() > expected);
        assert!((hotter.hydraulics.branches[0].flow.value()-analytic_loop_flow(370.0)).abs() < 2e-9);
        let mut unaccelerated = controls();
        unaccelerated.acceleration = None;
        let plain = fixture.solve(cx, unaccelerated).unwrap();
        assert!((plain.hydraulics.branches[0].flow.value()-expected).abs() < 2e-9);
    });
}

#[test]
fn g3_segment_lengths_and_reverse_declarations_preserve_integrated_gravity() {
    with_cx(|cx| {
        let mut forward = Fixture::new(false, 350.0);
        forward.branches[1] = NaturalCirculationBranch {
            model: BranchThermalModel::Exchange(vec![
                AirSegment::new("mild return", 1.0, 10.0).unwrap(),
                AirSegment::new("cold return", 1.0, 70.0).unwrap(),
            ]),
            axial_fractions: vec![0.1, 0.9],
            wall_temperatures_k: vec![315.0, 295.0],
        };
        let mut reversed = Fixture::new(true, 350.0);
        reversed.branches[1] = forward.branches[1].clone();
        if let BranchThermalModel::Exchange(rows) = &mut reversed.branches[1].model { rows.reverse(); }
        reversed.branches[1].axial_fractions.reverse();
        reversed.branches[1].wall_temperatures_k.reverse();
        let a = forward.solve(cx, controls()).unwrap();
        let b = reversed.solve(cx, controls()).unwrap();
        assert!((a.hydraulics.branches[0].flow.value()-b.hydraulics.branches[0].flow.value()).abs() < 1e-11);
        assert!((a.hydraulics.branches[1].flow.value()+b.hydraulics.branches[1].flow.value()).abs() < 1e-11);
        for i in 0..2 {
            assert!((a.branch_mean_temperatures_k[i]-b.branch_mean_temperatures_k[i]).abs() < 1e-8);
            assert!((a.thermal_mean_densities_kg_m3[i]-b.thermal_mean_densities_kg_m3[i]).abs() < 1e-10);
        }
        let segments = &a.transport.branches[1].march.as_ref().unwrap().segments;
        let means = [segments[0].reference_temperature_k, segments[1].reference_temperature_k];
        let expected_mean = 0.1*means[0] + 0.9*means[1];
        assert!((a.branch_mean_temperatures_k[1]-expected_mean).abs() < 1e-11);
        assert!((expected_mean-0.5*(means[0]+means[1])).abs() > 1.0);
        let expected_density: f64 = segments.iter().zip([315.0, 295.0]).zip([0.1, 0.9])
            .map(|((segment, wall), weight)| weight*mean_gas_density(segment.inlet_temperature_k,
                segment.outlet_temperature_k, wall, segment.ntu)).sum();
        assert!((a.thermal_mean_densities_kg_m3[1]-expected_density).abs() < 1e-9);
    });
}

#[test]
fn g0_tiny_relaxation_does_not_pass_unclosed_feedback_or_source_budgets() {
    with_cx(|cx| {
        let fixture = Fixture::new(false, 350.0);
        let mut config = controls();
        config.max_iterations = 2;
        config.relaxation = 1e-12;
        config.acceleration = None;
        assert!(matches!(fixture.solve(cx, config), Err(NaturalCirculationError::NotConverged { .. })));
        config = controls();
        config.max_material_evaluations = 2;
        assert!(matches!(fixture.solve(cx, config), Err(NaturalCirculationError::MaterialBudget { evaluations: 2 })));
        config = controls();
        config.max_density_quadrature_panels = 2;
        assert!(matches!(fixture.solve(cx, config), Err(NaturalCirculationError::Quadrature { .. })));
        config = controls();
        config.maximum_neglected_dissipation_w = 0.0;
        assert!(matches!(fixture.solve(cx, config), Err(NaturalCirculationError::Dissipation { .. })));
        let mut problem = fixture.problem();
        problem.maximum_relative_density_departure = 0.001;
        assert!(matches!(solve_natural_circulation(cx, &problem, controls(), &[335.0, 310.0],
            |_, _, temperature| Ok(gas(temperature.value()))), Err(NaturalCirculationError::Branch { .. })));
    });
}

#[test]
fn g0_undetermined_stagnant_heat_and_malformed_geometry_are_not_reservoirs() {
    with_cx(|cx| {
        let mut fixture = Fixture::new(false, 350.0);
        assert!(matches!(solve_natural_circulation(cx, &fixture.problem(), controls(), &[325.0; 2],
            |_, _, temperature| Ok(gas(temperature.value()))), Err(NaturalCirculationError::Transport(_))));
        fixture.branches[0].axial_fractions[0] = 0.8;
        let mut calls = 0;
        assert!(matches!(solve_natural_circulation(cx, &fixture.problem(), controls(), &[335.0, 310.0],
            |_, _, temperature| { calls += 1; Ok(gas(temperature.value())) }), Err(NaturalCirculationError::Branch { .. })));
        assert_eq!(calls, 0);
    });
}

#[test]
fn g0_density_quadrature_admits_actual_endpoints_and_queries_constant_profiles_once() {
    with_cx(|cx| {
        let reference = gas(325.0).transport_state().unwrap();
        let outlet = 300.0 + 50.0*(-det::expm1(-0.5));
        for reject_inlet in [true, false] {
            let endpoint = if reject_inlet { 300.0 } else { outlet };
            let mut sampler = Sampler {
                material: |_: &Cx<'_>, _, temperature: Temperature| {
                    if temperature.value() == endpoint {
                        Err(MaterialError::Parameters { what: "endpoint outside source domain".into() })
                    } else { Ok(reference) }
                },
                reference, maximum_departure: 0.25, evaluations: 0, limit: 100,
            };
            assert!(matches!(sampler.integrate(cx, 0, 0, 300.0, 350.0, 0.5, 16, 1e-9),
                Err(NaturalCirculationError::Material(_))));
            assert_eq!(sampler.evaluations, if reject_inlet { 1 } else { 2 });
        }
        let mut constant = Sampler {
            material: |_: &Cx<'_>, _, _: Temperature| Ok(reference),
            reference, maximum_departure: 0.25, evaluations: 0, limit: 1,
        };
        assert_eq!(constant.integrate(cx, 0, 0, 300.0, 350.0, 0.0, 16, 1e-9).unwrap(),
            (reference.density_kg_m3(), 0.0));
        assert_eq!(constant.evaluations, 1);
    });
}

#[test]
fn g4_cancelled_circulation_does_not_query_materials() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 1, kernel_id: 1713, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let fixture = Fixture::new(false, 350.0);
        let mut calls = 0;
        assert!(matches!(solve_natural_circulation(&cx, &fixture.problem(), controls(), &[335.0, 310.0],
            |_, _, temperature| { calls += 1; Ok(gas(temperature.value())) }), Err(NaturalCirculationError::Interrupted)));
        assert_eq!(calls, 0);
    });
}
