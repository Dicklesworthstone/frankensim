//! Size two heat exchangers around a warm bypass to meet an outlet limit.
//!
//! Run: cargo run -p fs-airflow --example size_cooling -- 310
//! Test: cargo test -p fs-airflow --example size_cooling
//!
//! The argument is a kelvin limit. Both external streams enter at 340 K; all
//! exchanger walls are prescribed at 290 K. A real hydraulic graph determines
//! the flows before adjoint-guided thermal sizing. Geometry, fluid properties
//! and fixed wall temperatures are illustrative declarations, not validation.
//! This is NOT a shared-solid FEM design: the wall thermostat is an input.

use std::error::Error;
use fs_airflow::conjugate::AirSegment;
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::graph::thermal::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork};
use fs_airflow::graph::thermal::transport::sensitivity::TransportGradient;
use fs_airflow::graph::thermal::transport::sensitivity::design::{UniformCoolingDesign, UniformCoolingRequest};
use fs_airflow::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};

type Failure = Box<dyn Error>;
const WALL: f64 = 290.0;
const SUPPLY: f64 = 340.0;
const DENSITY: f64 = 1.2;
const CP: f64 = 1007.0;

struct CoolingRun {
    design: UniformCoolingDesign,
    gradient: TransportGradient,
    baseline_outlet_k: f64,
    flows_m3_s: Vec<f64>,
}
fn branch(name: &str, from: usize, to: usize, resistance: f64) -> Result<GraphBranch, Failure> {
    Ok(GraphBranch { from, to, loss: LossElement::new(name, LossResistance::new(resistance), 0.0,
        SourceProvenance::new("declared quadratic sizing example", "size-cooling-example-v1"), ToleranceBasis::Analytic)? })
}
fn run(cx: &Cx<'_>, limit: f64) -> Result<CoolingRun, Failure> {
    // Nominal p_mix = 5 Pa, with 3 L/s through the upstream exchanger,
    // 1 L/s through the bypass, and 4 L/s through the downstream exchanger.
    let graph = LossGraph::new(4, vec![
        branch("upstream-exchanger", 0, 2, 15.0 / 0.003_f64.powi(2))?,
        branch("unheated-bypass", 1, 2, 15.0 / 0.001_f64.powi(2))?,
        branch("downstream-exchanger", 2, 3, 5.0 / 0.004_f64.powi(2))?,
    ])?;
    let flow = graph.solve(&[
        FixedPressure { node: 0, pressure: Pressure::new(20.0) },
        FixedPressure { node: 1, pressure: Pressure::new(20.0) },
        FixedPressure { node: 3, pressure: Pressure::new(0.0) },
    ], GraphSolveConfig { max_sweeps: 4096, max_node_iterations: 80,
        absolute_flow_tolerance: VolumetricFlowRate::new(1.0e-14), relative_flow_tolerance: 1.0e-12 }, cx)?;
    let network = TransportNetwork::new(cx, &flow,
        TransportAir { density: Density::new(DENSITY), specific_heat_j_kg_k: CP }, vec![
            BranchThermalModel::Exchange(vec![AirSegment::new("upstream-wall", 0.01, 50.0)?]),
            BranchThermalModel::Adiabatic,
            BranchThermalModel::Exchange(vec![AirSegment::new("downstream-wall", 0.01, 80.0)?]),
        ], &[
            TransportInlet { node: 0, temperature: Temperature::new(SUPPLY) },
            TransportInlet { node: 1, temperature: Temperature::new(SUPPLY) },
        ], TransportConfig { absolute_flow_tolerance: VolumetricFlowRate::new(1.0e-12),
            relative_flow_tolerance: 1.0e-10, absolute_heat_tolerance_w: 1.0e-7, relative_heat_tolerance: 1.0e-9 })?;
    let baseline_outlet_k = network.march(cx, &[WALL, WALL])?.node_temperatures_k[3]
        .ok_or("the outlet has no supplied air")?;
    let design = network.size_uniform_cooling(cx, UniformCoolingRequest { node: 3,
        wall_temperature: Temperature::new(WALL), outlet_limit: Temperature::new(limit),
        minimum_scale: 0.05, maximum_scale: 50.0, temperature_tolerance_k: 1.0e-7,
        log_scale_tolerance: 1.0e-8, max_evaluations: 100 })?;
    let chosen = network.scaled_conductance(cx, design.scale)?;
    let linearization = chosen.linearize(cx, &[WALL, WALL])?;
    let mut objective = linearization.zero_objective(); objective.node_temperatures[3] = 1.0;
    let gradient = linearization.pullback(cx, &objective)?;
    Ok(CoolingRun { design, gradient, baseline_outlet_k,
        flows_m3_s: flow.branches.iter().map(|edge| edge.flow.value()).collect() })
}
fn main() -> Result<(), Failure> {
    let mut args = std::env::args().skip(1);
    let limit = args.next().map_or(Ok(310.0), |s| s.parse::<f64>())?;
    if args.next().is_some() { return Err("usage: size_cooling [outlet-limit-kelvin]".into()); }
    let gate = CancelGate::new();
    let result = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 7, kernel_id: 713, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        run(&cx, limit)
    })?;
    println!("Nominal fixed-flow/fixed-wall design; no FEM wall feedback or uncertainty certificate.");
    println!("Branch flows (m^3/s): {:?}", result.flows_m3_s);
    println!("Baseline outlet: {:.9} K; requested limit: {:.9} K", result.baseline_outlet_k, limit);
    println!("Conductance scale: {:.12}; passing outlet: {:.9} K", result.design.scale, result.design.outlet_temperature.value());
    println!("Chosen conductances: {:.9}, {:.9} W/K", 0.5 * result.design.scale, 0.8 * result.design.scale);
    println!("Lower scale: {:.12}; lower outlet: {:.9} K; already feasible at minimum: {}",
        result.design.lower_scale, result.design.lower_outlet_temperature.value(), result.design.at_lower_bound);
    println!("Sizing evaluations: {}; dT_out/dln(scale): {:.9} K", result.design.evaluations, result.design.slope_per_log_scale_k);
    println!("dT_out/dln(U_i) (K): {:?}", result.gradient.log_conductances);
    println!("dT_out/dT_wall_i: {:?}; dT_out/dT_supply_node: {:?}", result.gradient.walls, result.gradient.inlets);
    println!("Heat leaving walls: {:.9} W; raw heat imbalance: {:.3e} W",
        result.design.march.wall_heat_rate_w, result.design.march.heat_imbalance_w);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn reference(scale: f64, flows: &[f64]) -> f64 {
        let first = WALL + (SUPPLY-WALL)*(-0.5*scale/(DENSITY*CP*flows[0])).exp();
        let mixed = (flows[0]*first+flows[1]*SUPPLY)/(flows[0]+flows[1]);
        WALL + (mixed-WALL)*(-0.8*scale/(DENSITY*CP*flows[2])).exp()
    }
    #[test]
    fn sizing_example_matches_an_independent_series_mixing_formula() {
        let gate=CancelGate::new();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx=Cx::new(&gate,arena,StreamKey {seed:7,kernel_id:713,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
            for limit in [300.0,310.0,320.0] {
                let result=run(&cx,limit).unwrap(); let mut low=0.05; let mut high=50.0;
                for _ in 0..100 { let mid=0.5*(low+high);
                    if reference(mid,&result.flows_m3_s)>limit {low=mid;} else {high=mid;} }
                assert!((result.design.scale-high).abs()<1e-6);
                assert!(result.design.outlet_temperature.value()<=limit);
                assert!(result.design.lower_outlet_temperature.value()>limit);
                assert!((reference(result.design.scale,&result.flows_m3_s)-result.design.outlet_temperature.value()).abs()<1e-9);
                let slope:f64=result.gradient.log_conductances.iter().sum();
                assert!((slope-result.design.slope_per_log_scale_k).abs()<1e-12);
            }
            assert!(run(&cx,280.0).is_err());
        });
    }
}
