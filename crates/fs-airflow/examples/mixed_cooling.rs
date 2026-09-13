//! Pressure-driven cooling with a heated stream, an explicit fresh-air bypass,
//! a mixing junction, and a downstream surface of the SAME finite-element slab.
//!
//! Run: cargo run -p fs-airflow --example mixed_cooling -- 330 290
//! Test: cargo test -p fs-airflow --example mixed_cooling
//!
//! Inputs are the two reservoir temperatures in kelvin. Constant density,
//! specific heat, wall coefficients, geometry and conductivity are illustrative
//! declarations, not experimentally validated data. This example uses the real
//! hydraulic, transport, coupling and FEM implementations without extra crates.

use std::error::Error;

use fs_airflow::conjugate::{AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::thermal::coupled_transport::{CoupledTransportSolution, solve_coupled_transport};
use fs_airflow::graph::thermal::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::{AirflowError, LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, ConductionSolution, InitialGuess, SolveConfig, solve};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};

type Failure = Box<dyn Error>;
const LENGTH: f64 = 0.05;
const WIDTH: f64 = 0.1;
const HEIGHT: f64 = 0.1;
const AREA: f64 = WIDTH * HEIGHT;
const CONDUCTIVITY: f64 = 10.0;
const DENSITY: f64 = 1.2;
const CP: f64 = 1007.0;
const H_FIRST: f64 = 50.0;
const H_LAST: f64 = 80.0;
const FIRST_FLOW: f64 = 0.003;
const BYPASS_FLOW: f64 = 0.001;

#[derive(Debug)]
struct CoolingRun {
    coupled: CoupledTransportSolution,
    flows_m3_s: Vec<f64>,
    analytic_heat_w: f64,
    analytic_faces_k: [f64; 2],
    max_nodal_error_k: f64,
    robin_total_w: f64,
}

fn branch(name: &str, from: usize, to: usize, flow_at_ten_pa: f64) -> Result<GraphBranch, Failure> {
    Ok(GraphBranch { from, to, loss: LossElement::new(name,
        LossResistance::new(10.0 / flow_at_ten_pa.powi(2)), 0.0,
        SourceProvenance::new("declared quadratic mixed-cooling example", "mixed-cooling-example-v1"),
        ToleranceBasis::Analytic)? })
}

fn slab(cx: &Cx<'_>, mesh: &ConductionMesh, references: &[f64]) -> Result<ConductionSolution, Failure> {
    let material = ConductivityModel::isotropic_declared(CONDUCTIVITY)?;
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(mesh)
        .region("first-face", |face| on_box_face(face.centroid[0], 0.0), ThermalBc::robin(H_FIRST, references[0])?)?
        .region("last-face", |face| on_box_face(face.centroid[0], LENGTH), ThermalBc::robin(H_LAST, references[1])?)?
        .adiabatic_remainder().finish()?;
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(0.5 * references[0] + 0.5 * references[1]);
    config.linear.tolerance = 1e-12;
    config.stop.residual_rtol = 1e-10;
    Ok(solve(cx, ConductionProblem { mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source }, config)?)
}

fn simulate(first_inlet: f64, bypass_inlet: f64) -> Result<CoolingRun, Failure> {
    for temperature in [first_inlet, bypass_inlet] {
        if !(temperature.is_finite() && temperature > 0.0) {
            return Err("reservoir temperatures must be finite positive kelvin values".into());
        }
    }
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<CoolingRun, Failure> {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 41, kernel_id: 714, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        // Both supplies are at 20 Pa, the merge solves to 10 Pa, the sink is
        // at 0 Pa. The declared resistances give exactly 3/1/4 litres per second.
        let graph = LossGraph::new(4, vec![branch("first-stream", 0, 2, FIRST_FLOW)?,
            branch("fresh-bypass", 1, 2, BYPASS_FLOW)?, branch("mixed-stream", 2, 3, FIRST_FLOW + BYPASS_FLOW)?])?;
        let flow = graph.solve(&[
            FixedPressure { node: 0, pressure: Pressure::new(20.0) },
            FixedPressure { node: 1, pressure: Pressure::new(20.0) },
            FixedPressure { node: 3, pressure: Pressure::new(0.0) },
        ], GraphSolveConfig { max_sweeps: 4096, max_node_iterations: 80,
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-13), relative_flow_tolerance: 1e-11 }, &cx)?;
        let network = TransportNetwork::new(&cx, &flow,
            TransportAir { density: Density::new(DENSITY), specific_heat_j_kg_k: CP },
            vec![BranchThermalModel::Exchange(vec![AirSegment::new("first-face", AREA, H_FIRST)?]),
                BranchThermalModel::Adiabatic,
                BranchThermalModel::Exchange(vec![AirSegment::new("last-face", AREA, H_LAST)?])],
            &[TransportInlet { node: 0, temperature: Temperature::new(first_inlet) },
                TransportInlet { node: 1, temperature: Temperature::new(bypass_inlet) }],
            TransportConfig { absolute_flow_tolerance: VolumetricFlowRate::new(1e-12), relative_flow_tolerance: 1e-10,
                absolute_heat_tolerance_w: 1e-7, relative_heat_tolerance: 1e-8 })?;
        let (complex, positions) = box_grid([4, 2, 2], [LENGTH, WIDTH, HEIGHT]);
        let mesh = ConductionMesh::new(complex, positions)?;
        let mut last_solid = None;
        let mut solid_failure: Option<Failure> = None;
        let config = ConjugateConfig { balance_tolerance_w: 1e-8, max_iterations: 200, ..ConjugateConfig::default() };
        let coupled = solve_coupled_transport(&cx, &network, &config, |cx, references| {
            let evaluated = (|| -> Result<Vec<SolidRegionState>, Failure> {
                let solution = slab(cx, &mesh, references)?;
                let states = network.regions().into_iter().map(|name| {
                    solution.report.robin_fluxes.iter().find(|flux| flux.region == name)
                        .map(SolidRegionState::from_robin_flux)
                        .ok_or_else(|| -> Failure { format!("FEM report has no Robin region {name}").into() })
                }).collect::<Result<Vec<_>, Failure>>()?;
                last_solid = Some(solution);
                Ok(states)
            })();
            evaluated.map_err(|error| {
                // Preserve the original FEM error outside the channel callback;
                // this sentinel merely stops its enclosing coupling iteration.
                solid_failure = Some(error);
                AirflowError::Cancelled { iteration: 0, references_k: references.to_vec() }
            })
        });
        if let Some(error) = solid_failure { return Err(error); }
        let coupled = coupled?;
        let solid = last_solid.ok_or("coupling returned without a solid solve")?;

        // Independent closed form. If q is heat conducted first-face -> last-
        // face, the first air stream loses q, the mixed stream gains q, and
        // T_mix = T_bar - q/C_total. Each film has effective conductance
        // E = C(1-exp(-hA/C)). Thus the series denominator contains -1/C_total:
        // dropping that term would incorrectly restart the downstream inlet.
        let c_first = DENSITY * FIRST_FLOW * CP;
        let c_bypass = DENSITY * BYPASS_FLOW * CP;
        let c_total = c_first + c_bypass;
        let mean_inlet = (c_first * first_inlet + c_bypass * bypass_inlet) / c_total;
        let e_first = c_first * (-(-H_FIRST * AREA / c_first).exp_m1());
        let e_last = c_total * (-(-H_LAST * AREA / c_total).exp_m1());
        let resistance = LENGTH / (CONDUCTIVITY * AREA);
        let heat = (first_inlet - mean_inlet) / (resistance + 1.0 / e_first + 1.0 / e_last - 1.0 / c_total);
        let faces = [first_inlet - heat / e_first, mean_inlet - heat / c_total + heat / e_last];
        let max_nodal_error_k = mesh.positions().iter().zip(&solid.temperature)
            .map(|(position, &temperature)| (temperature - (faces[0] + (faces[1] - faces[0]) * position[0] / LENGTH)).abs())
            .fold(0.0_f64, f64::max);
        Ok(CoolingRun { coupled, flows_m3_s: flow.branches.iter().map(|branch| branch.flow.value()).collect(),
            analytic_heat_w: heat, analytic_faces_k: faces, max_nodal_error_k,
            robin_total_w: solid.report.energy.robin_out_w })
    })
}

fn main() -> Result<(), Failure> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!("Usage: mixed_cooling [FIRST_INLET_K BYPASS_INLET_K]\nDefaults: 330 290. Frozen-property, nominal hydraulic/thermal/FEM example.");
        return Ok(());
    }
    let (first, bypass) = match args.as_slice() {
        [] => (330.0, 290.0),
        [first, bypass] => (first.parse::<f64>()?, bypass.parse::<f64>()?),
        _ => return Err("usage: mixed_cooling [FIRST_INLET_K BYPASS_INLET_K]".into()),
    };
    let run = simulate(first, bypass)?;
    println!("Nominal mixed-stream / shared-FEM cooling; frozen illustrative inputs, no physical validation.");
    println!("Shared solid solves: {}", run.coupled.iterations);
    for (index, (branch, flow)) in run.coupled.transport.branches.iter().zip(&run.flows_m3_s).enumerate() {
        println!("Branch {index}: {flow:.9} m^3/s; inlet {:?} K -> outlet {:?} K", branch.inlet_temperature_k, branch.outlet_temperature_k);
    }
    for (index, solid) in run.coupled.solid.iter().enumerate() {
        println!("{}: wall {:.9} K (analytic {:.9} K), outward heat {:.9} W", solid.region,
            solid.mean_wall_temperature_k, run.analytic_faces_k[index], solid.heat_rate_w);
    }
    println!("Analytic through-slab heat: {:.9} W", run.analytic_heat_w);
    println!("Max FEM/analytic nodal difference: {:.3e} K", run.max_nodal_error_k);
    println!("Whole-domain Robin total: {:.3e} W; air heat imbalance: {:.3e} W", run.robin_total_w, run.coupled.transport.heat_imbalance_w);
    Ok(())
}

#[test]
fn real_fem_and_mixed_transport_match_the_linear_slab_closed_form() {
    for (first, bypass) in [(330.0, 290.0), (290.0, 330.0), (310.0, 310.0)] {
        let run = simulate(first, bypass).expect("real hydraulic, transport and FEM solve");
        assert!(run.max_nodal_error_k < 1e-5, "FEM error {} K", run.max_nodal_error_k);
        assert!((run.coupled.solid[1].heat_rate_w - run.analytic_heat_w).abs() < 1e-5);
        assert!((run.coupled.solid[0].heat_rate_w + run.analytic_heat_w).abs() < 1e-5);
        assert!(run.robin_total_w.abs() < 1e-6);
        let expected_outlet = (FIRST_FLOW * first + BYPASS_FLOW * bypass) / (FIRST_FLOW + BYPASS_FLOW);
        assert!((run.coupled.transport.node_temperatures_k[3].unwrap() - expected_outlet).abs() < 1e-6);
        for (&actual, expected) in run.flows_m3_s.iter().zip([FIRST_FLOW, BYPASS_FLOW, FIRST_FLOW + BYPASS_FLOW]) {
            assert!((actual - expected).abs() < 1e-10);
        }
    }
}
