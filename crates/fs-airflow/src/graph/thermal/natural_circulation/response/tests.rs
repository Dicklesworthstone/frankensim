//! Real material viscosity and source correlations in the native buoyant graph.
use super::*;
use crate::conjugate::AirSegment;
use crate::graph::{FixedPressure, GraphSolveConfig, LinearLossElement, LinearLossResistance,
    PassiveBranch};
use crate::graph::tests::{boundary, config as hydraulic_config, edge, with_cx};
use crate::graph::thermal::medium::convection::{ForcedConvectionSettings, evaluate_forced_from_material};
use crate::graph::thermal::transport::{TransportConfig, BranchThermalModel};
use fs_convection::CorrelationId;
use fs_material::gas::{GasSpec, GasState};
use fs_qty::Length;

const DIAMETER: f64 = 0.02;
const FLOW_AREA: f64 = std::f64::consts::PI*DIAMETER*DIAMETER/4.0;
const WALL_AREA: f64 = std::f64::consts::PI*DIAMETER;

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101325.0).unwrap()
}
fn config() -> CirculationClosureConfig {
    let hydraulics = GraphSolveConfig {
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
        relative_flow_tolerance: 1e-13,
        ..hydraulic_config()
    };
    CirculationClosureConfig {
        inner: NaturalCirculationConfig {
            hydraulics,
            transport: TransportConfig {
                absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
                relative_flow_tolerance: 1e-10,
                absolute_heat_tolerance_w: 1e-8,
                relative_heat_tolerance: 1e-9,
            },
            max_iterations: 100, max_material_evaluations: 200_000,
            max_density_quadrature_panels: 256,
            absolute_density_tolerance_kg_m3: 1e-9, relative_density_tolerance: 1e-9,
            absolute_flow_change_tolerance: VolumetricFlowRate::new(1e-12),
            relative_flow_change_tolerance: 1e-8, relaxation: 0.5,
            acceleration: Some(IqnIlsConfig::default()),
            maximum_neglected_dissipation_w: 0.01,
        },
        max_iterations: 60, max_response_evaluations: 120,
        max_hydraulic_solves: 4000, max_buoyancy_material_evaluations: 2_000_000,
        temperature_tolerance_k: 2e-6,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-11),
        relative_flow_tolerance: 1e-7,
        relative_coefficient_tolerance: 1e-7,
        relaxation: 0.6, temperature_scale_k: 325.0,
        flow_scale: VolumetricFlowRate::new(1e-4), acceleration: Some(IqnIlsConfig::default()),
    }
}
struct Fixture {
    graph: LossGraph, branches: Vec<NaturalCirculationBranch>,
    positions: [[f64;3];2], datums: [FixedPressure;1],
}
impl Fixture {
    fn new(hot: f64) -> Self {
        Self {
            graph: LossGraph::new(2, vec![edge("rise",0,1,1.0),edge("return",1,0,1.0)]).unwrap(),
            branches: [("hot",hot),("cold",300.0)].into_iter().map(|(name,wall)|
                NaturalCirculationBranch {
                    model: BranchThermalModel::Exchange(vec![
                        AirSegment::new(name,WALL_AREA,5.0).unwrap()]),
                    axial_fractions: vec![1.0], wall_temperatures_k: vec![wall],
                }).collect(),
            positions: [[0.0;3],[0.0,0.0,1.0]], datums: [boundary(0,0.0)],
        }
    }
    fn problem(&self) -> NaturalCirculationProblem<'_> {
        NaturalCirculationProblem {
            graph: &self.graph, reference: gas(325.0).transport_state().unwrap(),
            node_positions_m: &self.positions, acceleration_m_s2: [0.0,0.0,-9.81],
            maximum_relative_density_departure: 0.3, pressure_datums: &self.datums,
            branches: &self.branches,
        }
    }
    fn initial(&self) -> CirculationClosureState {
        CirculationClosureState {
            branch_mean_temperatures_k: vec![335.0,315.0],
            segment_mean_temperatures_k: vec![vec![335.0],vec![315.0]],
            volume_flows_m3_s: vec![1e-4,1e-4],
        }
    }
}

// Independent Hagen-Poiseuille geometry and actual material viscosity. Heat
// uses the original Hausen correlation, not a guessed temperature-to-h fit.
fn response(fixture: &Fixture, state: &CirculationClosureState)
    -> CirculationResponse<Vec<(GasState,GasState,f64)>> {
    let radius = DIAMETER/2.0;
    let geometry = 8.0/(std::f64::consts::PI*radius.powi(4));
    let mut paths = Vec::new();
    let mut models = Vec::new();
    let mut owners = Vec::new();
    for (i, branch) in fixture.graph.branches().iter().enumerate() {
        let bulk = gas(state.branch_mean_temperatures_k[i]);
        let properties = gas(state.segment_mean_temperatures_k[i][0]);
        let loss = LinearLossElement::new(branch.loss.name.clone(),
            LinearLossResistance::new(geometry*bulk.transport_state().unwrap().dynamic_viscosity_pa_s()),
            branch.loss.uncertainty_rel, branch.loss.source.clone(), branch.loss.tolerance_basis).unwrap();
        paths.push(PassiveBranch::Linear { from:branch.from, to:branch.to, loss });
        let heat = evaluate_forced_from_material(CorrelationId::CircularDuctHausen,
            &properties, state.volume_flows_m3_s[i].abs()/FLOW_AREA, Length::new(DIAMETER),
            ForcedConvectionSettings { length_over_hydraulic_diameter: Some(1.0/DIAMETER),
                ..ForcedConvectionSettings::default() }).unwrap();
        let h = heat.coefficient().value.value();
        let name = if i == 0 { "hot" } else { "cold" };
        models.push(BranchThermalModel::Exchange(vec![AirSegment::new(name,WALL_AREA,h).unwrap()]));
        owners.push((bulk,properties,h));
    }
    CirculationResponse { graph:LossGraph::new_mixed(2,paths).unwrap(), models, payload:owners }
}

fn run(fixture: &Fixture, cx: &Cx<'_>, controls: CirculationClosureConfig,
    seed: &CirculationClosureState) -> Result<CoupledCirculationSolution<Vec<(GasState,GasState,f64)>>,
        CirculationClosureError<&'static str>> {
    solve_circulation_with_response(cx,&fixture.problem(),controls,seed,
        |_,_,temperature| GasState::try_new(&GasSpec::dry_air_ussa1976(),temperature.value(),101325.0),
        |_,state| Ok(response(fixture,state)))
}

#[test]
fn material_viscosity_and_flow_dependent_h_reach_one_native_circulation() {
    with_cx(|cx| {
        let fixture=Fixture::new(350.0);
        let controls=config();
        let result=run(&fixture,cx,controls,&fixture.initial()).unwrap();
        assert!(result.max_temperature_residual_k <= controls.temperature_tolerance_k);
        assert!(result.max_relative_coefficient_residual <= controls.relative_coefficient_tolerance);
        let actual=solved_state(&result.solution);
        let confirmed=response(&fixture,&actual);
        assert!(coefficient_residual(&result.response,&confirmed) <= controls.relative_coefficient_tolerance);
        let q=result.solution.hydraulics.branches[0].flow.value();
        assert!(q > 0.0 && q < 1e-3);
        assert!((result.solution.hydraulics.branches[1].flow.value()-q).abs() < 1e-12);
        assert!((result.response.payload[0].2-result.response.payload[1].2).abs() > 1e-4);
        let old=gas(325.0).transport_state().unwrap().dynamic_viscosity_pa_s();
        let actual_mu=result.response.payload[0].0.transport_state().unwrap().dynamic_viscosity_pa_s();
        assert!((actual_mu/old-1.0).abs() > 1e-3);
        assert!(result.hydraulic_solves <= controls.max_hydraulic_solves);
        assert!(result.buoyancy_material_evaluations <= controls.max_buoyancy_material_evaluations);
        let hotter=Fixture::new(370.0);
        let more=run(&hotter,cx,controls,&hotter.initial()).unwrap();
        assert!(more.solution.hydraulics.branches[0].flow.value() > q);
    });
}

#[test]
fn opposite_property_flow_seeds_do_not_invent_stagnation() {
    with_cx(|cx| {
        let fixture=Fixture::new(350.0);
        let base=run(&fixture,cx,config(),&fixture.initial()).unwrap();
        let mut seed=solved_state(&base.solution);
        for q in &mut seed.volume_flows_m3_s { *q = -*q; }
        let mut controls=config();
        controls.acceleration=None;
        controls.relaxation=0.5;
        let reversed=run(&fixture,cx,controls,&seed).unwrap();
        for (a,b) in reversed.solution.hydraulics.branches.iter().zip(&base.solution.hydraulics.branches) {
            assert!(a.flow.value() > 0.0);
            assert!((a.flow.value()-b.flow.value()).abs() < 5e-10);
        }
    });
}

#[test]
fn tiny_relaxation_and_changed_coefficients_cannot_bypass_physical_residuals() {
    with_cx(|cx| {
        let fixture=Fixture::new(350.0);
        let mut controls=config();
        controls.max_iterations=1;
        controls.acceleration=None;
        controls.relaxation=1e-14;
        assert!(matches!(run(&fixture,cx,controls,&fixture.initial()),
            Err(CirculationClosureError::NotConverged { .. })));

        // Deliberately loose state tolerances isolate the fresh coefficient
        // gate. This declared bounded coefficient jump must still refuse.
        controls.temperature_tolerance_k=100.0;
        controls.absolute_flow_tolerance=VolumetricFlowRate::new(0.1);
        let mut seed=fixture.initial();
        seed.branch_mean_temperatures_k[0]=350.0;
        seed.segment_mean_temperatures_k[0][0]=350.0;
        let result=solve_circulation_with_response(cx,&fixture.problem(),controls,&seed,
            |_,_,t| GasState::try_new(&GasSpec::dry_air_ussa1976(),t.value(),101325.0),
            |_,state| {
                let mut value=response(&fixture,state);
                let h=if state.branch_mean_temperatures_k[0] >= 349.0 { 10.0 } else { 5.0 };
                value.models[0]=BranchThermalModel::Exchange(vec![
                    AirSegment::new("hot",WALL_AREA,h).unwrap()]);
                Ok::<_,&'static str>(value)
            });
        assert!(matches!(result,Err(CirculationClosureError::NotConverged { coefficient_residual, .. })
            if coefficient_residual > 0.1));
    });
}

#[test]
fn response_cannot_change_physical_regions_and_global_work_budgets_are_enforced() {
    with_cx(|cx| {
        let fixture=Fixture::new(350.0);
        let malformed=solve_circulation_with_response(cx,&fixture.problem(),config(),&fixture.initial(),
            |_,_,t| GasState::try_new(&GasSpec::dry_air_ussa1976(),t.value(),101325.0),
            |_,state| {
                let mut value=response(&fixture,state);
                value.models[0]=BranchThermalModel::Exchange(vec![
                    AirSegment::new("hot",2.0*WALL_AREA,5.0).unwrap()]);
                Ok::<_,&'static str>(value)
            });
        assert!(matches!(malformed,Err(CirculationClosureError::Invalid(_))));
        let mut limits=config();
        limits.max_buoyancy_material_evaluations=1;
        assert!(matches!(run(&fixture,cx,limits,&fixture.initial()),
            Err(CirculationClosureError::Natural(NaturalCirculationError::MaterialBudget { .. }))));
    });
}
