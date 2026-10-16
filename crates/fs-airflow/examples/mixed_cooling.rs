//! Pressure-driven cooling with a heated stream, an explicit fresh-air bypass,
//! a mixing junction, and a downstream surface of the SAME finite-element slab.
//!
//! Run: cargo run -p fs-airflow --example mixed_cooling -- 330 290
//! Design: cargo run -p fs-airflow --example mixed_cooling -- 330 290 323
//! Test: cargo test -p fs-airflow --example mixed_cooling
//!
//! The optional third argument is an upper limit on the first wall temperature,
//! in kelvin. Bounded adjoint-guided sizing varies the last face's h in [10,1000]
//! W/(m^2 K); BOTH the solid Robin operator and the transported air respond.
//! This sizes an effective coefficient at fixed geometry/hydraulics, not fins,
//! fan speed or a physically validated hardware design. Other properties are
//! illustrative frozen declarations. No outside crates or solver are used.

use std::error::Error;

use fs_airflow::conjugate::{AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::thermal::coupled_transport::{CoupledTransportSolution, solve_coupled_transport};
use fs_airflow::graph::thermal::coupled_transport::sensitivity::{CoupledGradient, CoupledLinearization, InterfaceSolveConfig};
use fs_airflow::graph::thermal::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::{AirflowError, LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::adjoint::robin::RobinLinearization;
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, InitialGuess, SolveConfig};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::det;
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
const DESIGN_H_BOUNDS: [f64; 2] = [10.0, 1000.0];
const DESIGN_EVALUATIONS: usize = 48;
const DESIGN_TEMPERATURE_TOLERANCE_K: f64 = 1e-5;
const DESIGN_LOG_H_TOLERANCE: f64 = 1e-5;

#[derive(Debug)]
struct CoolingRun {
    coupled: CoupledTransportSolution,
    first_wall_gradient: CoupledGradient,
    last_htc: f64,
    flows_m3_s: Vec<f64>,
    analytic_heat_w: f64,
    analytic_faces_k: [f64; 2],
    max_nodal_error_k: f64,
    robin_total_w: f64,
}

#[derive(Debug)]
struct SizedCooling {
    passing: CoolingRun,
    /// None means the declared lower coefficient already satisfies the limit.
    failing_htc: Option<f64>,
    evaluations: usize,
}

fn branch(name: &str, from: usize, to: usize, flow_at_ten_pa: f64) -> Result<GraphBranch, Failure> {
    Ok(GraphBranch { from, to, loss: LossElement::new(name,
        LossResistance::new(10.0 / flow_at_ten_pa.powi(2)), 0.0,
        SourceProvenance::new("declared quadratic mixed-cooling example", "mixed-cooling-example-v1"),
        ToleranceBasis::Analytic)? })
}

fn slab(cx: &Cx<'_>, mesh: &ConductionMesh, references: &[f64], last_htc: f64) -> Result<RobinLinearization, Failure> {
    let material = ConductivityModel::isotropic_declared(CONDUCTIVITY)?;
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(mesh)
        .region("first-face", |face| on_box_face(face.centroid[0], 0.0), ThermalBc::robin(H_FIRST, references[0])?)?
        .region("last-face", |face| on_box_face(face.centroid[0], LENGTH), ThermalBc::robin(last_htc, references[1])?)?
        .adiabatic_remainder().finish()?;
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(0.5 * references[0] + 0.5 * references[1]);
    config.linear.tolerance = 1e-12;
    config.stop.residual_rtol = 1e-10;
    Ok(RobinLinearization::new(cx, ConductionProblem { mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source }, config, &["first-face", "last-face"])?)
}

fn simulate(first_inlet: f64, bypass_inlet: f64) -> Result<CoolingRun, Failure> {
    simulate_at_htc(first_inlet, bypass_inlet, H_LAST)
}

fn simulate_at_htc(first_inlet: f64, bypass_inlet: f64, last_htc: f64) -> Result<CoolingRun, Failure> {
    for value in [first_inlet, bypass_inlet, last_htc] {
        if !(value.is_finite() && value > 0.0) {
            return Err("reservoir temperatures and transfer coefficients must be finite and positive".into());
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
                BranchThermalModel::Exchange(vec![AirSegment::new("last-face", AREA, last_htc)?])],
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
                let linear = slab(cx, &mesh, references, last_htc)?;
                let states = network.regions().into_iter().map(|name| {
                    linear.primal().report.robin_fluxes.iter().find(|flux| flux.region == name)
                        .map(SolidRegionState::from_robin_flux)
                        .ok_or_else(|| -> Failure { format!("FEM report has no Robin region {name}").into() })
                }).collect::<Result<Vec<_>, Failure>>()?;
                last_solid = Some(linear);
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
        // No re-solve at a proposed air reference: this is the very FEM field
        // accepted by the primal driver, including its retained Robin operator.
        let linear = CoupledLinearization::new(&cx, &network, &solid, &config)?;
        let mut objective = linear.zero_objective();
        objective.wall_temperatures[0] = 1.0;
        let first_wall_gradient = linear.pullback(&cx, &objective, InterfaceSolveConfig {
            max_iterations: 400, absolute_tolerance: 1e-10, relative_tolerance: 1e-10, relaxation: 1.0,
        })?;

        // Independent continuous slab/film oracle, not used to pick a design.
        let (heat, faces) = closed_form(first_inlet, bypass_inlet, H_FIRST, last_htc);
        let max_nodal_error_k = mesh.positions().iter().zip(&solid.primal().temperature)
            .map(|(position, &temperature)| (temperature - (faces[0] + (faces[1] - faces[0]) * position[0] / LENGTH)).abs())
            .fold(0.0_f64, f64::max);
        Ok(CoolingRun { coupled, first_wall_gradient, last_htc,
            flows_m3_s: flow.branches.iter().map(|branch| branch.flow.value()).collect(),
            analytic_heat_w: heat, analytic_faces_k: faces, max_nodal_error_k,
            robin_total_w: solid.primal().report.energy.robin_out_w })
    })
}

fn closed_form(first: f64, bypass: f64, first_htc: f64, last_htc: f64) -> (f64, [f64; 2]) {
    let c_first = DENSITY * FIRST_FLOW * CP;
    let c_total = DENSITY * (FIRST_FLOW + BYPASS_FLOW) * CP;
    let mean_inlet = (FIRST_FLOW * first + BYPASS_FLOW * bypass) / (FIRST_FLOW + BYPASS_FLOW);
    let e_first = c_first * (-(-first_htc * AREA / c_first).exp_m1());
    let e_last = c_total * (-(-last_htc * AREA / c_total).exp_m1());
    // T_mix = mean_inlet - q/C_total. Omitting -1/C_total would incorrectly
    // restart the downstream inlet and remove part of the thermal feedback.
    let resistance = LENGTH / (CONDUCTIVITY * AREA);
    let heat = (first - mean_inlet) / (resistance + 1.0 / e_first + 1.0 / e_last - 1.0 / c_total);
    (heat, [first - heat / e_first, mean_inlet - heat / c_total + heat / e_last])
}

/// Specialized to this fixed-flow, unpowered slab with first inlet hotter than
/// bypass. Its first-face temperature decreases as last h increases. No generic
/// monotonicity assumption is imposed on arbitrary conjugate cooling systems.
fn size_last_exchanger(first: f64, bypass: f64, target: f64) -> Result<SizedCooling, Failure> {
    if !(first.is_finite() && bypass.is_finite() && target.is_finite()
        && first > bypass && bypass > 0.0 && target > 0.0)
    { return Err("sizing requires finite positive kelvin inputs and FIRST_INLET_K > BYPASS_INLET_K".into()); }
    let low_run = simulate_at_htc(first, bypass, DESIGN_H_BOUNDS[0])?;
    if low_run.coupled.solid[0].mean_wall_temperature_k <= target {
        return Ok(SizedCooling { passing: low_run, failing_htc: None, evaluations: 1 });
    }
    let mut passing = simulate_at_htc(first, bypass, DESIGN_H_BOUNDS[1])?;
    if passing.coupled.solid[0].mean_wall_temperature_k > target {
        return Err(format!("target {target} K is unattainable within h = {:?} W/(m^2 K); upper-bound wall is {} K",
            DESIGN_H_BOUNDS, passing.coupled.solid[0].mean_wall_temperature_k).into());
    }
    let mut low = det::ln(DESIGN_H_BOUNDS[0]);
    let mut high = det::ln(DESIGN_H_BOUNDS[1]);
    let mut evaluations = 2;
    loop {
        let wall = passing.coupled.solid[0].mean_wall_temperature_k;
        if high - low <= DESIGN_LOG_H_TOLERANCE && target - wall <= DESIGN_TEMPERATURE_TOLERANCE_K {
            return Ok(SizedCooling { passing, failing_htc: Some(det::exp(low)), evaluations });
        }
        if evaluations == DESIGN_EVALUATIONS {
            return Err(format!("sizing exhausted {evaluations} full coupled evaluations before both kelvin and log(h) tolerances held").into());
        }
        // The derivative includes solid feedback. The independently printed
        // closed form never enters this Newton step or its acceptance decision.
        let slope = passing.first_wall_gradient.log_htc[1];
        let proposal = high - (wall - target) / slope;
        let padding = 0.05 * (high - low);
        let candidate = if slope < 0.0 && proposal.is_finite()
            && proposal > low + padding && proposal < high - padding { proposal }
            else { 0.5 * low + 0.5 * high };
        let trial = simulate_at_htc(first, bypass, det::exp(candidate))?;
        evaluations += 1;
        if trial.coupled.solid[0].mean_wall_temperature_k <= target {
            high = candidate;
            passing = trial;
        } else { low = candidate; }
    }
}

fn main() -> Result<(), Failure> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!("Usage: mixed_cooling [FIRST_INLET_K BYPASS_INLET_K [TARGET_FIRST_WALL_K]]\nDefaults: 330 290. Optional target sizes last-face h in [10,1000] W/(m^2 K).\nFrozen-property, nominal hydraulic/thermal/FEM example; no hardware validation.");
        return Ok(());
    }
    let (first, bypass, target) = match args.as_slice() {
        [] => (330.0, 290.0, None),
        [first, bypass] => (first.parse::<f64>()?, bypass.parse::<f64>()?, None),
        [first, bypass, target] => (first.parse::<f64>()?, bypass.parse::<f64>()?, Some(target.parse::<f64>()?)),
        _ => return Err("usage: mixed_cooling [FIRST_INLET_K BYPASS_INLET_K [TARGET_FIRST_WALL_K]]".into()),
    };
    let run = if let Some(target) = target {
        let sized = size_last_exchanger(first, bypass, target)?;
        println!("First-wall limit {target} K; {} full coupled design evaluations; failing h {:?}, passing h {:.9} W/(m^2 K)",
            sized.evaluations, sized.failing_htc, sized.passing.last_htc);
        sized.passing
    } else { simulate(first, bypass)? };
    println!("Nominal mixed-stream / shared-FEM cooling; frozen illustrative inputs, no physical validation.");
    println!("Shared solid solves at returned design: {}", run.coupled.iterations);
    for (index, (branch, flow)) in run.coupled.transport.branches.iter().zip(&run.flows_m3_s).enumerate() {
        println!("Branch {index}: {flow:.9} m^3/s; inlet {:?} K -> outlet {:?} K", branch.inlet_temperature_k, branch.outlet_temperature_k);
    }
    for (index, solid) in run.coupled.solid.iter().enumerate() {
        println!("{}: wall {:.9} K (analytic {:.9} K), outward heat {:.9} W", solid.region,
            solid.mean_wall_temperature_k, run.analytic_faces_k[index], solid.heat_rate_w);
    }
    println!("First-wall dT/dln(h): {:?} K; dT/dinlet: {:?}; interface adjoint residual {:.3e} in {} sweeps",
        run.first_wall_gradient.log_htc, run.first_wall_gradient.inlets,
        run.first_wall_gradient.interface_residual, run.first_wall_gradient.iterations);
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

#[test]
fn example_adjoint_matches_the_independent_continuous_slab_derivative() {
    let run = simulate(330.0, 290.0).unwrap();
    let eps = 1e-4_f64;
    for control in 0..2 {
        let mut plus = [H_FIRST, H_LAST]; let mut minus = plus;
        plus[control] *= eps.exp(); minus[control] *= (-eps).exp();
        let fd = (closed_form(330.0, 290.0, plus[0], plus[1]).1[0]
            - closed_form(330.0, 290.0, minus[0], minus[1]).1[0]) / (2.0 * eps);
        assert!((run.first_wall_gradient.log_htc[control] - fd).abs() < 1e-5);
    }
    assert!((run.first_wall_gradient.inlets[0] + run.first_wall_gradient.inlets[1] - 1.0).abs() < 1e-7);
}

#[test]
fn sizing_returns_the_passing_fem_design_and_matches_an_independent_inverse() {
    let target = 323.0;
    let sized = size_last_exchanger(330.0, 290.0, target).unwrap();
    let c_first = DENSITY * FIRST_FLOW * CP;
    let c_total = DENSITY * (FIRST_FLOW + BYPASS_FLOW) * CP;
    let e_first = c_first * (-(-H_FIRST * AREA / c_first).exp_m1());
    let heat = (330.0 - target) * e_first;
    let required_r_last = 10.0 / heat - LENGTH / (CONDUCTIVITY * AREA) - 1.0 / e_first + 1.0 / c_total;
    let exact_h = -c_total / AREA * (-1.0 / (required_r_last * c_total)).ln_1p();
    assert!(sized.passing.coupled.solid[0].mean_wall_temperature_k <= target);
    assert!(target - sized.passing.coupled.solid[0].mean_wall_temperature_k <= DESIGN_TEMPERATURE_TOLERANCE_K);
    assert!((sized.passing.last_htc / exact_h - 1.0).abs() < 2e-5);
    let lower = sized.failing_htc.expect("nontrivial sizing bracket");
    assert!(closed_form(330.0, 290.0, H_FIRST, lower).1[0] > target - 1e-7);
    assert!((sized.passing.last_htc / lower).ln() <= DESIGN_LOG_H_TOLERANCE * 1.01);
    assert!(sized.evaluations <= DESIGN_EVALUATIONS);
    let minimum = size_last_exchanger(330.0, 290.0, 329.0).unwrap();
    assert_eq!(minimum.evaluations, 1);
    assert!(minimum.failing_htc.is_none());
    assert_eq!(minimum.passing.last_htc, DESIGN_H_BOUNDS[0]);
    assert!(size_last_exchanger(330.0, 290.0, 321.0).unwrap_err().to_string().contains("unattainable"));
    assert!(size_last_exchanger(290.0, 330.0, 310.0).is_err());
}
