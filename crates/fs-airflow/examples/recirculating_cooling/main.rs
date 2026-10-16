//! Material/pump tables -> finite recirculating inventory <-> actual solid FEM.
//! Source readers and every physics/iteration algorithm are existing owners.
#[path = "../pressure_cooling/source.rs"]
mod source;
#[path = "../pressure_cooling/pump_source.rs"]
mod pump_source;
mod surface;
use surface::SurfaceTable;

use fs_airflow::AirflowError;
use fs_airflow::conjugate::{ConjugateConfig, Relaxation, SolidRegionState};
use fs_airflow::graph::thermal::medium::FlowHeatGeometry;
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{LiquidResistanceSection, PressureDrivenBudget};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::{ParallelLiquidBranch, ParallelLiquidNetwork, ParallelPressureBudget};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::solid::ParallelLiquidCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::PumpOperatingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::solid::PumpCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::recirculation::ReservoirStepControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::{
    StoredCoolingControls, StoredCoolingState, StoredCoolingStep, StoredSolidProblem,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::radiation::RadiatingStoredCoolingControls;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, ThermalBoundary,
    InitialGuess, LinearConfig, Nonlinearity, ScalarField, SolveConfig, StopRule, ThermalBoundaryBuilder, solve};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_convection::{CorrelationId, CorrelationRobinBoundary};
use fs_conduction::transient::enthalpy::{EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig, EnthalpyStepSolution};
use fs_solver::NewtonKrylovConfig;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use source::{LiquidTable, ThermalTable};
use pump_source::PumpTable;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const REGIONS: [&str;2] = ["first","last"];
const INITIAL_K: f64 = 300.0;

#[derive(Clone, Copy, Debug, PartialEq)]
struct Settings {
    mass_kg: f64, dt_s: f64, steps: usize, body_w: f64, reservoir_w: f64,
    inlet_bracket_k: (f64,f64), pressure_bracket_pa: (f64,f64),
}
#[derive(Debug)]
struct Sample {
    time_s: f64, reservoir_k: f64, return_k: f64, mass_flow_kg_s: f64,
    walls_k: [f64;2], pressure_pa: f64, stored_j: f64, return_j: f64,
    external_j: f64, residual_j: f64, ignored_hydraulic_j: f64,
    solid_heat_w: f64, solid_calls: usize,
    storage: Option<SolidStorageSample>,
    radiation: Option<RadiationSample>,
}
#[derive(Debug)]
struct SolidStorageSample {
    stored_j: f64, external_j: f64, interface_residual_j: f64, combined_residual_j: f64,
}

#[derive(Debug)]
struct RadiationSample {
    outward_j: f64, external_convection_j: f64, physical_residual_j: f64,
    enthalpy_solves: usize,
}

fn stored_sample<R>(result: &StoredCoolingStep<'_, '_, R>, solid: &EnthalpyStepSolution,
    coolant_out_w: f64, radiation: Option<RadiationSample>) -> Sample {
    let coupled = &result.liquid.circuit.state;
    Sample {
        time_s: result.liquid.end_time_s, reservoir_k: result.liquid.final_state.temperature_k(),
        return_k: coupled.fluid.mixed_outlet.temperature_k(), mass_flow_kg_s: coupled.fluid.total_mass_flow_kg_per_s,
        walls_k: [coupled.solid.regions[0].mean_wall_temperature_k, coupled.solid.regions[1].mean_wall_temperature_k],
        pressure_pa: result.liquid.circuit.network_pressure_pa,
        stored_j: result.liquid.stored_enthalpy_change_j, return_j: result.liquid.returned_heat_j,
        external_j: result.liquid.external_heat_j, residual_j: result.liquid.energy_residual_j,
        ignored_hydraulic_j: result.liquid.neglected_hydraulic_energy_j,
        solid_heat_w: coolant_out_w, solid_calls: result.liquid.total_solid_solves,
        storage: Some(SolidStorageSample { stored_j: solid.stored_energy_change_j,
            external_j: result.solid_external_heat_j, interface_residual_j: result.interface_energy_residual_j,
            combined_residual_j: result.combined_energy_residual_j }), radiation,
    }
}

fn mapped_boundary(mesh: &ConductionMesh, rows: &[CorrelationRobinBoundary])
    -> std::result::Result<ThermalBoundary, ConductionError>
{
    mapped_boundary_with_surface(mesh, rows, None)
}

fn mapped_boundary_with_surface(mesh: &ConductionMesh, rows: &[CorrelationRobinBoundary], surface: Option<&SurfaceTable>)
    -> std::result::Result<ThermalBoundary, ConductionError>
{
    let mut builder = ThermalBoundaryBuilder::new(mesh);
    for (index, name) in REGIONS.iter().enumerate() {
        builder = builder.region(name, |face| face.outward_normal[2]>0.99
            && on_box_face(face.centroid[2],0.01)
            && if index==0 {face.centroid[0]<0.03} else {face.centroid[0]>=0.03},
            rows[index].boundary_condition().clone())?;
    }
    if let Some(surface) = surface {
        builder = builder.region(surface::EXPOSED_REGION,
            |face| face.outward_normal[2] < -0.99 && on_box_face(face.centroid[2],0.0),
            fs_conduction::ThermalBc::robin(surface.convection_w_m2_k, surface.air_temperature_k)?)?;
    }
    builder.adiabatic_remainder().finish()
}

fn run(table: &LiquidTable, pump: &PumpTable, settings: Settings) -> Result<Vec<Sample>> {
    run_with_solid(table, pump, None, settings)
}

fn run_with_solid(table: &LiquidTable, pump: &PumpTable, solid_table: Option<&ThermalTable>,
    settings: Settings) -> Result<Vec<Sample>> {
    run_with_surface(table, pump, solid_table, None, settings)
}

fn run_with_surface(table: &LiquidTable, pump: &PumpTable, solid_table: Option<&ThermalTable>,
    surface: Option<&SurfaceTable>, settings: Settings) -> Result<Vec<Sample>> {
    if surface.is_some() && solid_table.is_none() {
        return Err("--radiation requires an explicit --solid storage material".into());
    }
    if !(1..=128).contains(&settings.steps) || !settings.dt_s.is_finite() || settings.dt_s <= 0.0
        || !(settings.dt_s*settings.steps as f64).is_finite()
        || !settings.body_w.is_finite() || !settings.reservoir_w.is_finite()
    { return Err("finite heat rates, positive dt and 1..=128 steps required".into()); }
    let law = LiquidTemperatureLaw::try_new(&table.card,&table.point(INITIAL_K)?,"T",
        LiquidViscosityBasis::Dynamic,MaterialPropertySelection::SingleClaimOnly)?;
    let cp = PropertyKey::with_quantity("specific_heat_capacity",QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density",QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    let caloric = ResolvedHeatCapacityCurve::try_new(&table.card,HeatCapacitySourceQuery {
        heat_capacity:&cp,density:&rho,temperature_axis:"T",lower:&table.point(table.lower_k)?,
        upper:&table.point(table.upper_k)?,selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Liquid,0.0,0.5)?;
    let geometry = FlowHeatGeometry { flow_area:Area::new(core::f64::consts::PI*0.01*0.01/4.0),
        characteristic_length:Length::new(0.01),heated_length:Length::new(2.0),wetted_area:Area::new(0.003),aspect_ratio:None };
    let sections = REGIONS.map(|region|LiquidHeatSection {region,geometry,correlation:CorrelationId::CircularDuctLaminarCwt});
    let resistance = MaterialFlowResistance::try_new(5e9,0.0,geometry.flow_area,geometry.characteristic_length,(0.0,2300.0))?;
    let losses = REGIONS.map(|region|[LiquidResistanceSection {region,resistance}]);
    let (topology,positions) = box_grid([2,1,1],[0.06,0.1,0.01]);
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(),topology.tets),positions)?;
    // Optional solid data supplies BOTH storage and conductivity at its own
    // fixed source context. The old no-storage path keeps its original 5 W/m/K.
    let solid_source = solid_table.map(|solid| -> Result<_> {
        let chart = ResolvedHeatCapacityCurve::try_new(&solid.card, HeatCapacitySourceQuery {
            heat_capacity: &cp, density: &rho, temperature_axis: "T",
            lower: &solid.point(solid.lower_k)?, upper: &solid.point(solid.upper_k)?,
            selection: SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, 0.0, 0.01)?;
        let conductivity = PropertyKey::with_quantity("thermal_conductivity", QuantitySpec::semantic(
            SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static)));
        let model = ConductivityModel::isotropic(ConductivityTable::from_claims_at_query_points(
            solid.card.claims(), &conductivity, "T", &[solid.point(solid.lower_k)?, solid.point(solid.upper_k)?],
            SelectionPolicy::SingleClaimOnly)?);
        let initial_density = solid.card.claims().query_typed(&rho, &solid.point(INITIAL_K)?,
            SelectionPolicy::SingleClaimOnly)?.evidence.value.value;
        // Query the SAME chart for initial h; do not integrate another Cp law.
        let initial_h = chart.state_at_temperature(INITIAL_K)?.specific_enthalpy_j_kg();
        Ok((chart, model, initial_density, initial_h))
    }).transpose()?;
    let declared_material = ConductivityModel::isotropic_declared(5.0)?;
    let solid_material = solid_source.as_ref().map_or(&declared_material, |s| &s.1);
    let power = ScalarField::nodal("declared recirculating solid heat",mesh.vertex_count(),
        vec![settings.body_w/(0.06*0.1*0.01);mesh.vertex_count()])?;
    let branch = PressureDrivenBudget {mass_flow_bracket_kg_per_s:(0.0002,0.012),
        pressure_absolute_tolerance_pa:1e-9,pressure_relative_tolerance:1e-9,max_flow_evaluations:64,
        maximum_neglected_hydraulic_power_w:0.01};
    let controls = PumpCouplingControls {
        network:ParallelLiquidCouplingControls {
            coupling:ConjugateConfig {temperature_tolerance_k:1e-7,max_iterations:64,
                balance_tolerance_w:1e-6,balance_relative_tolerance:1e-8,relaxation:Relaxation::Fixed {omega:1.0}},
            local:LiquidFeedbackBudget {temperature_tolerance_k:1e-9,coefficient_relative_tolerance:1e-9,max_region_iterations:64},
            acceleration:Default::default(),mass_flow_relative_tolerance:1e-8,
            work:ParallelPressureBudget {max_flow_evaluations:1_000_000,maximum_neglected_hydraulic_power_w:0.01},
        },
        operating_point:PumpOperatingControls {pressure_bracket_pa:settings.pressure_bracket_pa,
            absolute_pressure_tolerance_pa:1e-7,relative_pressure_tolerance:1e-9,max_network_evaluations:64},
        max_total_solid_solves:8192,
    };
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<_> {
        let cx = Cx::new(&gate,arena,StreamKey {seed:947,kernel_id:755,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let paths = sections.iter().map(|s|LiquidHeatPath::try_new(&cx,&law,&caloric,
            MassFlowRate::new(0.001),Temperature::new(INITIAL_K),core::slice::from_ref(s),LiquidPropertyTemperature::BulkMean))
            .collect::<std::result::Result<Vec<_>,_>>()?;
        let descriptions:Vec<_> = (0..2).map(|i|ParallelLiquidBranch {
            name:REGIONS[i],path:&paths[i],resistances:&losses[i],pressure_budget:branch,
        }).collect();
        let network = ParallelLiquidNetwork::try_new(&cx,&descriptions)?;
        let mut tank = network.reservoir(&cx,settings.mass_kg,Temperature::new(INITIAL_K),0.0)?;
        let spatial_storage = solid_source.as_ref().map(|(chart, _, density, _)| {
            EnthalpyBackwardEuler::uniform(&cx, &mesh, chart.curve(), *density,
                EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() })
        }).transpose()?;
        let mut stored_state = match (&spatial_storage, &solid_source) {
            (Some(storage), Some((_, _, _, h))) => Some(StoredCoolingState::try_new(
                &cx, tank.clone(), storage, vec![*h; mesh.vertex_count()])?),
            _ => None,
        };
        let mut initial_walls = [INITIAL_K;2];
        let mut samples = Vec::with_capacity(settings.steps);
        for i in 0..settings.steps {
            let end = settings.dt_s*(i+1) as f64;
            let duration = end-stored_state.as_ref().map_or(tank.time_s(), |s| s.reservoir().time_s());
            let step = ReservoirStepControls {temperature_bracket_k:settings.inlet_bracket_k,
                energy_absolute_tolerance_j:1e-4,max_inlet_evaluations:64,
                maximum_neglected_hydraulic_energy_j:0.01*duration};
            if let Some(state) = stored_state.as_mut() {
                let problem = StoredSolidProblem { mesh: &mesh, material: solid_material, source: &power,
                    element_materials: None, interfaces: None };
                let thermal = StoredCoolingControls {
                    pump: controls, reservoir: step,
                    solid: EnthalpyStepConfig { newton: NewtonKrylovConfig {
                        absolute_tolerance: 1e-10, relative_tolerance: 1e-11,
                        linear_restart: 24, max_linear_cycles: 16, forcing_minimum: 1e-12,
                        forcing_maximum: 1e-3, ..Default::default()
                    }, max_newton_iterations: 32, energy_tolerance_j: 1e-8 },
                    joint_energy_tolerance_j: 5e-4,
                };
                let sample = if let Some(surface) = surface {
                    let answer = network.advance_stored_solid_with_pump_and_radiation(&cx, state, end,
                        settings.reservoir_w*duration, &initial_walls, &pump.curve, problem,
                        core::slice::from_ref(&surface.patch), RadiatingStoredCoolingControls {
                            thermal, radiation: fs_conduction::AmbientRadiationConfig {
                                max_iterations: 48, temperature_tolerance_k: 1e-9,
                                balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-10, relaxation: 1.0,
                            }, max_total_enthalpy_solves: 65536,
                        }, |_, rows| mapped_boundary_with_surface(&mesh, rows, Some(surface)))?;
                    let response = &answer.step.liquid.circuit.state.solid.state;
                    stored_sample(&answer.step, &response.thermal.conduction, response.coolant_out_w,
                        Some(RadiationSample {
                            outward_j: duration*response.thermal.radiation.nonlinear_radiation_out_w,
                            external_convection_j: duration*response.external_convection_out_w,
                            physical_residual_j: response.thermal.physical_residual_norm_j,
                            enthalpy_solves: answer.total_enthalpy_solves,
                        }))
                } else {
                    let result = network.advance_stored_solid_with_pump(&cx, state, end,
                        settings.reservoir_w*duration, &initial_walls, &pump.curve, problem,
                        thermal, |_, rows| mapped_boundary(&mesh, rows))?;
                    let solid = &result.liquid.circuit.state.solid.state;
                    stored_sample(&result, solid, solid.robin_out_w, None)
                };
                initial_walls = sample.walls_k;
                samples.push(sample);
                continue;
            }
            let mut original_error = None;
            let result = network.advance_reservoir_with_pump_coupled_iqn(&cx,&mut tank,end,
                settings.reservoir_w*duration,&initial_walls,&pump.curve,&controls,step,|cx,rows| {
                let attempted = (|| {
                    let boundary = mapped_boundary(&mesh, rows)?;
                    solve(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,material:solid_material,
                        source:&power,element_materials:None},SolveConfig {
                        nonlinearity:Nonlinearity::FixedPoint {relaxation:1.0,max_backtracks:8},
                        stop:StopRule {residual_rtol:1e-10,residual_atol:1e-18,step_atol:0.0,max_iterations:12},
                        linear:LinearConfig {tolerance:1e-10,max_iterations:10000,restart:64},
                        initial:InitialGuess::Uniform(INITIAL_K),
                    })
                })();
                let field = attempted.map_err(|error:ConductionError| {
                    original_error=Some(error.to_string());
                    AirflowError::InvalidConjugateInput {field:"recirculation example solid; original error retained",value_bits:0}
                })?;
                Ok(LiquidSolidResponse {regions:field.report.robin_fluxes.iter().map(SolidRegionState::from_robin_flux).collect(),state:field})
            }).map_err(|error| match original_error {
                Some(original)=>format!("{error}; original solid error: {original}"),None=>error.to_string(),
            })?;
            let coupled = &result.circuit.state;
            initial_walls = [coupled.solid.regions[0].mean_wall_temperature_k,coupled.solid.regions[1].mean_wall_temperature_k];
            samples.push(Sample {time_s:tank.time_s(),reservoir_k:tank.state().temperature_k(),
                return_k:coupled.fluid.mixed_outlet.temperature_k(),mass_flow_kg_s:coupled.fluid.total_mass_flow_kg_per_s,
                walls_k:initial_walls,pressure_pa:result.circuit.network_pressure_pa,
                stored_j:result.stored_enthalpy_change_j,return_j:result.returned_heat_j,external_j:result.external_heat_j,
                residual_j:result.energy_residual_j,ignored_hydraulic_j:result.neglected_hydraulic_energy_j,
                solid_heat_w:coupled.solid.state.report.energy.robin_out_w,solid_calls:result.total_solid_solves,storage:None,radiation:None});
        }
        Ok(samples)
    })
}

fn parse_args(a: &[String]) -> Result<(Option<&str>, &str, &str, Settings)> {
    let (solid, a) = if a.first().is_some_and(|v| v == "--solid") {
        if a.len() < 2 { return Err("--solid requires a solid material table".into()); }
        (Some(a[1].as_str()), &a[2..])
    } else { (None, a) };
    if a.len()!=11 {return Err("usage: recirculating_cooling [--solid SOLID.txt] MATERIAL.txt PUMP.txt MASS_KG DT_S STEPS BODY_W RESERVOIR_W MIN_T_K MAX_T_K MIN_DP_PA MAX_DP_PA".into());}
    Ok((solid, &a[0], &a[1], Settings {mass_kg:a[2].parse()?,dt_s:a[3].parse()?,steps:a[4].parse()?,
        body_w:a[5].parse()?,reservoir_w:a[6].parse()?,inlet_bracket_k:(a[7].parse()?,a[8].parse()?),
        pressure_bracket_pa:(a[9].parse()?,a[10].parse()?)}))
}

fn parse_radiation_args(a: &[String]) -> Result<(Option<&str>, Option<&str>, &str, &str, Settings)> {
    let (surface, args) = if a.first().is_some_and(|v| v == "--radiation") {
        if a.len()<2 { return Err("--radiation needs a surface file followed by --solid SOLID.txt".into()); }
        (Some(a[1].as_str()), &a[2..])
    } else { (None, a) };
    let (solid, fluid, pump, settings) = parse_args(args)?;
    if surface.is_some() && solid.is_none() {
        return Err("--radiation requires --solid SOLID.txt before the existing arguments".into());
    }
    Ok((surface, solid, fluid, pump, settings))
}

fn main() -> Result<()> {
    let a:Vec<_> = std::env::args().skip(1).collect();
    let (surface_path, solid_path, fluid_path, pump_path, settings) = parse_radiation_args(&a)?;
    let table = LiquidTable::parse(&std::fs::read_to_string(fluid_path)?)?;
    let pump = PumpTable::parse(&std::fs::read_to_string(pump_path)?)?;
    let solid = solid_path.map(|path| -> Result<_> {
        ThermalTable::parse_solid(&std::fs::read_to_string(path)?)
    }).transpose()?;
    let surface = surface_path.map(|path| -> Result<_> {
        SurfaceTable::parse(&std::fs::read_to_string(path)?)
    }).transpose()?;
    let samples = match &surface {
        Some(surface) => run_with_surface(&table, &pump, solid.as_ref(), Some(surface), settings)?,
        None => match &solid {
            Some(solid) => run_with_solid(&table, &pump, Some(solid), settings)?,
            None => run(&table, &pump, settings)?,
        },
    };
    eprintln!("model=well_mixed_fixed_pressure_reservoir; solid={}; hydraulic_heat=bounded_omission; calibration=not_established",
        if solid.is_some() { "transient_reference_enthalpy_FEM" } else { "quasi_steady_FEM" });
    eprintln!("material={:?}; pump={:?}; pump_license={}",table.card.content_hash(),pump.curve.source(),pump.license);
    if let Some(solid) = &solid { eprintln!("solid_material={:?}", solid.card.content_hash()); }
    if let Some(surface) = &surface {
        eprintln!("surface_material={:?}; surface_license={}; radiation=implicit_area_mean_gray_patch; emissivity=constant_source_snapshot",
            surface.patch.emissivity().card_identity(), surface.license);
    }
    print!("time_s,reservoir_k,return_k,total_mass_kg_s,wall_first_k,wall_last_k,pressure_pa,stored_change_j,return_heat_j,external_heat_j,balance_j,ignored_hydraulic_j,solid_heat_w,solid_calls");
    if solid.is_some() { print!(",solid_stored_change_j,solid_external_heat_j,interface_balance_j,combined_balance_j"); }
    if surface.is_some() { print!(",radiation_out_j,external_convection_out_j,physical_solid_residual_j,enthalpy_solves"); }
    println!();
    for s in samples {
        print!("{:.9},{:.9},{:.9},{:.12e},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.9},{}",
            s.time_s,s.reservoir_k,s.return_k,s.mass_flow_kg_s,s.walls_k[0],s.walls_k[1],s.pressure_pa,
            s.stored_j,s.return_j,s.external_j,s.residual_j,s.ignored_hydraulic_j,s.solid_heat_w,s.solid_calls);
        if let Some(storage) = s.storage { print!(",{:.9},{:.9},{:.12e},{:.12e}",
            storage.stored_j,storage.external_j,storage.interface_residual_j,storage.combined_residual_j); }
        if let Some(radiation) = s.radiation { print!(",{:.9},{:.9},{:.12e},{}",
            radiation.outward_j,radiation.external_convection_j,radiation.physical_residual_j,radiation.enthalpy_solves); }
        println!();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const FIXED:&str=include_str!("../pressure_cooling/constant.txt");
    const VARIABLE:&str=include_str!("../pressure_cooling/temperature_dependent.txt");
    const PUMP:&str=include_str!("../pressure_cooling/pump.txt");
    fn settings()->Settings {Settings {mass_kg:0.1,dt_s:5.0,steps:2,body_w:20.0,reservoir_w:-5.0,
        inlet_bracket_k:(290.0,315.0),pressure_bracket_pa:(20.0,80.0)}}
    #[test]
    fn real_fem_recirculation_closes_step_energy_and_accumulates_the_actual_inventory() {
        let table=LiquidTable::parse(FIXED).unwrap();let pump=PumpTable::parse(PUMP).unwrap();
        let rows=run(&table,&pump,settings()).unwrap();
        for s in rows {
            assert!((s.reservoir_k-(300.0+15.0*s.time_s/200.0)).abs()<2e-5);
            assert!((s.stored_j-75.0).abs()<2e-4);
            assert!((s.return_j-100.0).abs()<2e-4);
            assert!((s.solid_heat_w-20.0).abs()<2e-5);
            assert!(s.residual_j.abs()<=1e-4);
            assert!(s.walls_k.iter().all(|&t|t>s.reservoir_k));
        }
        let cold=run(&table,&pump,Settings {body_w:-5.0,reservoir_w:0.0,steps:1,..settings()}).unwrap();
        assert!((cold[0].reservoir_k-299.875).abs()<2e-5);
    }
    #[test]
    fn data_only_swap_changes_thermal_history_pump_flow_and_solid_field() {
        let pump=PumpTable::parse(PUMP).unwrap();
        let fixed=run(&LiquidTable::parse(FIXED).unwrap(),&pump,settings()).unwrap();
        let varied=run(&LiquidTable::parse(VARIABLE).unwrap(),&pump,settings()).unwrap();
        assert!((fixed[1].reservoir_k-varied[1].reservoir_k).abs()>1e-4);
        assert!((fixed[1].mass_flow_kg_s-varied[1].mass_flow_kg_s).abs()>1e-6);
        assert_ne!(fixed[1].walls_k,varied[1].walls_k);
        assert!(run(&LiquidTable::parse(FIXED).unwrap(),&pump,Settings {steps:0,..settings()}).is_err());
    }

    #[test]
    fn stored_solid_option_preserves_the_legacy_argument_values() {
        let original: Vec<_> = ["fluid.txt", "pump.txt", "0.1", "5", "2", "20", "-5", "290", "315", "20", "80"]
            .into_iter().map(str::to_owned).collect();
        let (solid, fluid, pump, settings) = parse_args(&original).unwrap();
        assert!(solid.is_none()); assert_eq!(fluid, "fluid.txt"); assert_eq!(pump, "pump.txt");
        let mut stored = vec!["--solid".into(), "solid.txt".into()]; stored.extend(original);
        let (solid, _, _, with_storage) = parse_args(&stored).unwrap();
        assert_eq!(solid, Some("solid.txt")); assert_eq!(with_storage, settings);
        stored.push("unexpected".into()); assert!(parse_args(&stored).is_err());
        assert!(parse_args(&["--solid".into()]).is_err());
    }

    #[test]
    fn solid_material_data_changes_storage_and_coolant_history_through_real_fem() {
        let fluid = LiquidTable::parse(FIXED).unwrap();
        let pump = PumpTable::parse(PUMP).unwrap();
        let low = ThermalTable::parse_solid(include_str!("solid_low_capacity.txt")).unwrap();
        let high = ThermalTable::parse_solid(include_str!("solid_high_capacity.txt")).unwrap();
        let inputs = Settings { reservoir_w: 0.0, dt_s: 1.0, ..settings() };
        let a = run_with_solid(&fluid, &pump, Some(&low), inputs).unwrap();
        let b = run_with_solid(&fluid, &pump, Some(&high), inputs).unwrap();
        assert!(a[1].reservoir_k > b[1].reservoir_k);
        assert!(a[1].walls_k[0] > b[1].walls_k[0]);
        // A 20 W quasi-steady body would deliver 20 J per one-second interval.
        assert!(a[0].return_j > 0.0 && a[0].return_j < 5.0);
        for samples in [&a, &b] {
            let mut stored = 0.0;
            for sample in samples {
                let body = sample.storage.as_ref().unwrap();
                stored += body.stored_j + sample.stored_j;
                assert!((body.external_j - 20.0).abs() < 1e-8);
                assert!((stored - 20.0*sample.time_s).abs() < 1e-3);
                assert!(body.combined_residual_j.abs() <= 5e-4);
                assert!(body.interface_residual_j.abs() <= 5e-4);
            }
        }
        let cooled = run_with_solid(&fluid, &pump, Some(&low), Settings { body_w: -5.0, steps: 1, ..inputs }).unwrap();
        assert!(cooled[0].storage.as_ref().unwrap().stored_j < 0.0);
        assert!(cooled[0].reservoir_k < INITIAL_K);
        assert!((cooled[0].stored_j + cooled[0].storage.as_ref().unwrap().stored_j + 5.0).abs() < 5e-4);
    }
}
