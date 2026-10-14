//! Material table -> pressure-selected flow -> coolant enthalpy <-> actual FEM.
//! This example owns no material law, heat equation, or coupling algorithm.
mod source;
mod pump_source;

use fs_airflow::AirflowError;
use fs_airflow::conjugate::{ConjugateConfig, Relaxation, SolidRegionState};
use fs_airflow::graph::thermal::medium::FlowHeatGeometry;
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{
    LiquidResistanceSection, PressureDrivenBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::solid::PressureSolidControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::{
    ParallelLiquidBranch, ParallelLiquidNetwork, ParallelPressureBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::solid::ParallelLiquidCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::ParallelLiquidSolution;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{PumpHeadCurve, PumpOperatingControls};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::solid::PumpCouplingControls;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    InitialGuess, LinearConfig, Nonlinearity, ScalarField, SolveConfig, StopRule,
    ThermalBoundaryBuilder, solve};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_convection::{CorrelationId, CorrelationRobinBoundary};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, Pressure, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use source::LiquidTable;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const INLET_K: f64 = 300.0;
const REGIONS: [&str; 2] = ["first", "last"];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Topology { Series, Parallel }

#[derive(Clone, Copy, Debug, PartialEq)]
struct Settings {
    drive_pa: f64, heat_w: f64, minimum_mass: f64, maximum_mass: f64,
    topology: Topology,
}

#[derive(Debug)]
struct CoolingResult {
    mass_flow_kg_s: f64,
    outlet_k: f64,
    walls_k: Vec<f64>,
    heat_w: f64,
    solid_heat_w: f64,
    // Accepted pressure along the first complete inlet-to-outlet path.
    // All parallel path drops are reported separately, never summed.
    pressure_drop_pa: f64,
    branch_mass_flows_kg_s: Vec<f64>,
    branch_outlets_k: Vec<f64>,
    branch_pressure_drops_pa: Vec<f64>,
    neglected_hydraulic_power_w: f64,
    maximum_heat_imbalance_w: f64,
    solid_solves: usize,
    pump: Option<PumpReport>,
}

#[derive(Debug)]
struct PumpReport {
    head_m: f64, available_pressure_pa: f64, residual_pa: f64,
    inlet_volume_flow_m3_s: f64, hydraulic_power_w: f64,
    network_evaluations: usize, total_flow_evaluations: usize,
}

fn run(table: &LiquidTable, settings: Settings) -> Result<CoolingResult> {
    run_with_drive(table, settings, None)
}

fn run_with_drive(table: &LiquidTable, settings: Settings,
    pump: Option<(&PumpHeadCurve, PumpOperatingControls)>) -> Result<CoolingResult> {
    if !settings.heat_w.is_finite() { return Err("heat power must be finite (negative means cooling)".into()); }
    let law = LiquidTemperatureLaw::try_new(&table.card, &table.point(INLET_K)?, "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly)?;
    let cp = PropertyKey::with_quantity("specific_heat_capacity", QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])));
    let rho = PropertyKey::with_quantity("density", QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])));
    let caloric = ResolvedHeatCapacityCurve::try_new(&table.card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T",
        lower: &table.point(table.lower_k)?, upper: &table.point(table.upper_k)?,
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Liquid, 0.0, 0.5)?;
    let geometry = FlowHeatGeometry {
        flow_area: Area::new(core::f64::consts::PI * 0.01 * 0.01 / 4.0),
        characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
        wetted_area: Area::new(0.003), aspect_ratio: None,
    };
    let sections = REGIONS.map(|region| LiquidHeatSection {
        region, geometry, correlation: CorrelationId::CircularDuctLaminarCwt,
    });
    // Authored lumped hydraulic geometry; not a friction law inferred from the heat-transfer diameter.
    let loss = MaterialFlowResistance::try_new(5e9, 0.0, geometry.flow_area,
        geometry.characteristic_length, (0.0, 2300.0))?;
    let resistances = REGIONS.map(|region| LiquidResistanceSection { region, resistance: loss });
    let (topology, positions) = box_grid([2, 1, 1], [0.06, 0.1, 0.01]);
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), topology.tets), positions)?;
    let material = ConductivityModel::isotropic_declared(5.0)?;
    let heat = ScalarField::nodal("declared body power", mesh.vertex_count(),
        vec![settings.heat_w / (0.06 * 0.1 * 0.01); mesh.vertex_count()])?;
    let controls = PressureSolidControls {
        pressure: PressureDrivenBudget {
            mass_flow_bracket_kg_per_s: (settings.minimum_mass, settings.maximum_mass),
            pressure_absolute_tolerance_pa: 1e-7, pressure_relative_tolerance: 1e-7,
            max_flow_evaluations: 64, maximum_neglected_hydraulic_power_w: 0.01,
        },
        coupling: ConjugateConfig {
            temperature_tolerance_k: 1e-7, max_iterations: 64,
            balance_tolerance_w: 1e-5, balance_relative_tolerance: 1e-7,
            relaxation: Relaxation::Fixed { omega: 1.0 },
        },
        local: LiquidFeedbackBudget { temperature_tolerance_k: 1e-9,
            coefficient_relative_tolerance: 1e-9, max_region_iterations: 64 },
        acceleration: Default::default(), max_total_solid_solves: 2048,
    };
    let cancel = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<_> {
        let cx = Cx::new(&cancel, arena, StreamKey { seed: 91, kernel_id: 741, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let mut solid_error = None;
        // Both topologies use this SAME production FEM callback and physical
        // solid. Only the fluid connectivity changes; no surrogate solid or
        // topology-specific boundary coefficient is introduced.
        let mut solid = |cx: &Cx<'_>, rows: &[CorrelationRobinBoundary]| {
            let attempted = (|| {
                let mut builder = ThermalBoundaryBuilder::new(&mesh);
                for (index, name) in REGIONS.iter().enumerate() {
                    builder = builder.region(name, |face| {
                        face.outward_normal[2] > 0.99 && on_box_face(face.centroid[2], 0.01)
                            && if index == 0 { face.centroid[0] < 0.03 } else { face.centroid[0] >= 0.03 }
                    }, rows[index].boundary_condition().clone())?;
                }
                let boundary = builder.adiabatic_remainder().finish()?;
                solve(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
                    material: &material, source: &heat, element_materials: None }, SolveConfig {
                    nonlinearity: Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 8 },
                    stop: StopRule { residual_rtol: 1e-10, residual_atol: 1e-18, step_atol: 0.0, max_iterations: 12 },
                    linear: LinearConfig { tolerance: 1e-10, max_iterations: 10000, restart: 64 },
                    initial: InitialGuess::Uniform(INLET_K),
                })
            })();
            let solution = attempted.map_err(|error: ConductionError| {
                solid_error = Some(error.to_string());
                AirflowError::InvalidConjugateInput { field: "example solid solve; original error retained", value_bits: 0 }
            })?;
            let regions = solution.report.robin_fluxes.iter().map(SolidRegionState::from_robin_flux).collect();
            Ok(LiquidSolidResponse { regions, state: solution })
        };
        let attempted: std::result::Result<CoolingResult, String> = (|| {
            if let Some((curve, operating_point)) = pump {
                // One branch for the original series circuit, two for parallel.
                // Both use the exact same source/material and FEM callback above.
                let ranges = match settings.topology {
                    Topology::Series => vec![0..2], Topology::Parallel => vec![0..1, 1..2],
                };
                let paths = ranges.iter().map(|range| LiquidHeatPath::try_new(&cx, &law, &caloric,
                    MassFlowRate::new(settings.minimum_mass), Temperature::new(INLET_K),
                    &sections[range.clone()], LiquidPropertyTemperature::BulkMean))
                    .collect::<std::result::Result<Vec<_>, _>>().map_err(|e| e.to_string())?;
                let pressure_budget = PressureDrivenBudget {
                    pressure_absolute_tolerance_pa: 1e-10, pressure_relative_tolerance: 1e-11,
                    ..controls.pressure
                };
                let branches: Vec<_> = paths.iter().zip(&ranges).enumerate().map(|(i, (path, range))|
                    ParallelLiquidBranch { name: REGIONS[i], path,
                        resistances: &resistances[range.clone()], pressure_budget }).collect();
                let network = ParallelLiquidNetwork::try_new(&cx, &branches).map_err(|e| e.to_string())?;
                let pumped = PumpCouplingControls {
                    network: ParallelLiquidCouplingControls {
                        coupling: controls.coupling, local: controls.local, acceleration: controls.acceleration,
                        mass_flow_relative_tolerance: 1e-8,
                        work: ParallelPressureBudget { max_flow_evaluations: 250000,
                            maximum_neglected_hydraulic_power_w: 0.01 },
                    }, operating_point, max_total_solid_solves: 4096,
                };
                let point = network.solve_with_pump_coupled_iqn(&cx, &[INLET_K; 2], curve,
                    &pumped, &mut solid).map_err(|e| e.to_string())?;
                let report = PumpReport { head_m: point.supply.head_m,
                    available_pressure_pa: point.supply.available_pressure_pa, residual_pa: point.pressure_residual_pa,
                    inlet_volume_flow_m3_s: point.supply.inlet_volume_flow_m3_s,
                    hydraulic_power_w: point.supply.hydraulic_power_w,
                    network_evaluations: point.network_evaluations, total_flow_evaluations: point.total_flow_evaluations };
                let answer = point.state;
                return network_result(&answer.fluid, &answer.solid.regions,
                    answer.solid.state.report.energy.robin_out_w, answer.residuals.heat_imbalance_w,
                    point.total_solid_solves, Some(report));
            }
            match settings.topology {
                Topology::Series => {
                    let path = LiquidHeatPath::try_new(&cx, &law, &caloric,
                        MassFlowRate::new(settings.minimum_mass), Temperature::new(INLET_K),
                        &sections, LiquidPropertyTemperature::BulkMean).map_err(|e| e.to_string())?;
                    let answer = path.solve_pressure_coupled_iqn(&cx, &[INLET_K; 2], &resistances,
                        Pressure::new(settings.drive_pa), &controls, &mut solid).map_err(|e| e.to_string())?;
                    let point = answer.operating_point;
                    let coupled = point.thermal;
                    let outlet_k = coupled.fluid.regions.last()
                        .ok_or("accepted path lacks an outlet")?.outlet.temperature_k();
                    Ok(CoolingResult {
                        mass_flow_kg_s: point.mass_flow_kg_per_s, outlet_k,
                        walls_k: coupled.solid.regions.iter().map(|r| r.mean_wall_temperature_k).collect(),
                        heat_w: coupled.fluid.total_heat_rate_w,
                        solid_heat_w: coupled.solid.state.report.energy.robin_out_w,
                        pressure_drop_pa: point.hydraulics.total_pressure_drop_pa,
                        branch_mass_flows_kg_s: vec![point.mass_flow_kg_per_s],
                        branch_outlets_k: vec![outlet_k],
                        branch_pressure_drops_pa: vec![point.hydraulics.total_pressure_drop_pa],
                        neglected_hydraulic_power_w: point.hydraulics.total_dissipated_power_w,
                        maximum_heat_imbalance_w: coupled.maximum_heat_imbalance_w,
                        solid_solves: answer.total_solid_solves, pump: None,
                    })
                }
                Topology::Parallel => {
                    // Each existing surface now has its own stream from the
                    // common inlet; both streams exchange with the same solid.
                    let first = LiquidHeatPath::try_new(&cx, &law, &caloric,
                        MassFlowRate::new(settings.minimum_mass), Temperature::new(INLET_K),
                        &sections[..1], LiquidPropertyTemperature::BulkMean).map_err(|e| e.to_string())?;
                    let last = LiquidHeatPath::try_new(&cx, &law, &caloric,
                        MassFlowRate::new(settings.minimum_mass), Temperature::new(INLET_K),
                        &sections[1..], LiquidPropertyTemperature::BulkMean).map_err(|e| e.to_string())?;
                    // Resolve pressure more tightly than the outer wall and
                    // coefficient gates; branch-flow noise must not set their floor.
                    let branch_pressure = PressureDrivenBudget {
                        pressure_absolute_tolerance_pa: 1e-9, pressure_relative_tolerance: 1e-9,
                        ..controls.pressure
                    };
                    let network = ParallelLiquidNetwork::try_new(&cx, &[
                        ParallelLiquidBranch { name: REGIONS[0], path: &first,
                            resistances: &resistances[..1], pressure_budget: branch_pressure },
                        ParallelLiquidBranch { name: REGIONS[1], path: &last,
                            resistances: &resistances[1..], pressure_budget: branch_pressure },
                    ]).map_err(|e| e.to_string())?;
                    let parallel = ParallelLiquidCouplingControls {
                        coupling: controls.coupling, local: controls.local,
                        acceleration: controls.acceleration, mass_flow_relative_tolerance: 1e-8,
                        work: ParallelPressureBudget { max_flow_evaluations: 32768,
                            maximum_neglected_hydraulic_power_w: 0.01 },
                    };
                    let answer = network.solve_coupled_iqn(&cx, &[INLET_K; 2],
                        Pressure::new(settings.drive_pa), &parallel, &mut solid).map_err(|e| e.to_string())?;
                    network_result(&answer.fluid, &answer.solid.regions,
                        answer.solid.state.report.energy.robin_out_w, answer.residuals.heat_imbalance_w,
                        answer.iterations, None)
                }
            }
        })();
        // End the callback's mutable borrow before exposing its original error.
        drop(solid);
        attempted.map_err(|error| match solid_error {
            Some(original) => format!("{error}; solid: {original}").into(),
            None => error.into(),
        })
    })
}

fn parse_args(args: &[String]) -> Result<(&str, Settings)> {
    let (topology, args) = if args.first().is_some_and(|arg| arg == "--parallel") {
        (Topology::Parallel, &args[1..])
    } else { (Topology::Series, args) };
    if args.len() != 5 {
        return Err("usage: pressure_cooling [--parallel] MATERIAL.txt DRIVE_PA BODY_HEAT_W MIN_MASS_KG_S MAX_MASS_KG_S".into());
    }
    Ok((&args[0], Settings { drive_pa: args[1].parse()?, heat_w: args[2].parse()?,
        minimum_mass: args[3].parse()?, maximum_mass: args[4].parse()?, topology }))
}

fn network_result(fluid: &ParallelLiquidSolution<'_>, solid: &[SolidRegionState],
    solid_heat_w: f64, maximum_heat_imbalance_w: f64, solid_solves: usize, pump: Option<PumpReport>)
    -> std::result::Result<CoolingResult, String> {
    let branch_mass_flows_kg_s = fluid.branches.iter().map(|b| b.solution.mass_flow_kg_per_s).collect();
    let branch_outlets_k = fluid.branches.iter().map(|b| b.solution.thermal.regions.last()
        .map(|r| r.outlet.temperature_k()).ok_or("accepted branch lacks an outlet"))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let branch_pressure_drops_pa: Vec<_> = fluid.branches.iter()
        .map(|b| b.solution.hydraulics.total_pressure_drop_pa).collect();
    Ok(CoolingResult {
        mass_flow_kg_s: fluid.total_mass_flow_kg_per_s, outlet_k: fluid.mixed_outlet.temperature_k(),
        walls_k: solid.iter().map(|r| r.mean_wall_temperature_k).collect(),
        heat_w: fluid.total_heat_rate_w, solid_heat_w,
        pressure_drop_pa: branch_pressure_drops_pa[0],
        branch_mass_flows_kg_s, branch_outlets_k, branch_pressure_drops_pa,
        neglected_hydraulic_power_w: fluid.total_dissipated_power_w,
        maximum_heat_imbalance_w, solid_solves, pump,
    })
}

// Legacy invocations still use parse_args unchanged. A pump invocation supplies
// a pressure bracket instead of pretending its operating pressure is known.
fn parse_command(args: &[String]) -> Result<(&str, Settings, Option<(&str, PumpOperatingControls)>)> {
    let (topology, rest) = if args.first().is_some_and(|a| a == "--parallel") {
        (Topology::Parallel, &args[1..])
    } else { (Topology::Series, args) };
    if rest.first().is_some_and(|a| a == "--pump") {
        if rest.len() != 8 {
            return Err("usage: pressure_cooling [--parallel] --pump PUMP.txt MIN_DP_PA MAX_DP_PA MATERIAL.txt HEAT_W MIN_MASS_KG_S MAX_MASS_KG_S".into());
        }
        let bracket = (rest[2].parse()?, rest[3].parse()?);
        let controls = PumpOperatingControls { pressure_bracket_pa: bracket,
            absolute_pressure_tolerance_pa: 1e-6, relative_pressure_tolerance: 1e-8,
            max_network_evaluations: 48 };
        // drive_pa is unused by run_with_drive on the pump path; use the initial
        // bracket endpoint rather than an invented final pressure or NaN sentinel.
        Ok((&rest[4], Settings { drive_pa: bracket.0, heat_w: rest[5].parse()?,
            minimum_mass: rest[6].parse()?, maximum_mass: rest[7].parse()?, topology }, Some((&rest[1], controls))))
    } else {
        let (path, settings) = parse_args(args)?;
        Ok((path, settings, None))
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let (path, settings, pump_request) = parse_command(&args)?;
    let table = LiquidTable::parse(&std::fs::read_to_string(path)?)?;
    let pump_table = pump_request.map(|(path, _)| -> Result<_> {
        pump_source::PumpTable::parse(&std::fs::read_to_string(path)?)
    }).transpose()?;
    let result = match (&pump_table, pump_request) {
        (Some(table_source), Some((_, controls))) => run_with_drive(&table, settings, Some((&table_source.curve, controls)))?,
        _ => run(&table, settings)?,
    };
    println!("model=regional_liquid_pressure_cooling; phase_stability=assumed; calibration=not_established");
    println!("topology={}", match settings.topology { Topology::Series => "series", Topology::Parallel => "parallel" });
    println!("material_card={:?}", table.card.content_hash());
    println!("mass_flow_kg_s={:.12e}\noutlet_k={:.9}\nwall_temperatures_k={:?}", result.mass_flow_kg_s, result.outlet_k, result.walls_k);
    println!("coolant_heat_w={:.9}\nsolid_robin_heat_w={:.9}\npressure_drop_pa={:.9}", result.heat_w, result.solid_heat_w, result.pressure_drop_pa);
    println!("branch_mass_flows_kg_s={:?}\nbranch_outlets_k={:?}\nbranch_pressure_drops_pa={:?}",
        result.branch_mass_flows_kg_s, result.branch_outlets_k, result.branch_pressure_drops_pa);
    println!("neglected_hydraulic_power_w={:.12e}\nmaximum_regional_heat_imbalance_w={:.12e}\nsolid_solves={}",
        result.neglected_hydraulic_power_w, result.maximum_heat_imbalance_w, result.solid_solves);
    if let Some(pump) = &result.pump {
        println!("supply=pump; head_m={:.9}; available_pressure_pa={:.9}; pressure_residual_pa={:.12e}",
            pump.head_m, pump.available_pressure_pa, pump.residual_pa);
        println!("pump_inlet_volume_flow_m3_s={:.12e}; pump_hydraulic_power_w={:.12e}; network_evaluations={}; total_flow_evaluations={}",
            pump.inlet_volume_flow_m3_s, pump.hydraulic_power_w, pump.network_evaluations, pump.total_flow_evaluations);
        if let Some(table_source) = &pump_table {
            println!("pump_source={:?}; pump_license={:?}", table_source.curve.source(), table_source.license);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const FIXED: &str = include_str!("constant.txt");
    const VARIABLE: &str = include_str!("temperature_dependent.txt");
    fn settings(heat_w: f64) -> Settings {
        Settings { drive_pa: 50.0, heat_w, minimum_mass: 0.0002, maximum_mass: 0.008, topology: Topology::Series }
    }
    #[test]
    fn material_ingress_rejects_missing_context_bad_columns_and_nonphysical_rows() {
        assert!(LiquidTable::parse(FIXED).is_ok());
        for text in [FIXED.replace("pressure_pa=100000", ""),
            FIXED.replace("450 1000", "250 1000"), FIXED.replace("0.002", "NaN"),
            format!("{FIXED}\nsource=duplicate"), FIXED.replace("250 1000", "250 0 1000")]
        { assert!(LiquidTable::parse(&text).is_err()); }
    }
    #[test]
    fn data_only_material_swap_reaches_actual_pressure_and_fem_temperature() {
        let a = LiquidTable::parse(FIXED).unwrap();
        let b = LiquidTable::parse(VARIABLE).unwrap();
        assert_ne!(a.card.content_hash(), b.card.content_hash());
        let fixed = run(&a, settings(20.0)).unwrap();
        let variable = run(&b, settings(20.0)).unwrap();
        assert!((fixed.mass_flow_kg_s - 0.0025).abs() < 1e-9);
        assert!(variable.mass_flow_kg_s > fixed.mass_flow_kg_s);
        assert!(variable.outlet_k < fixed.outlet_k);
        assert_ne!(variable.walls_k, fixed.walls_k);
        for result in [&fixed, &variable] {
            assert!((result.heat_w - 20.0).abs() < 3e-5);
            assert!((result.solid_heat_w - 20.0).abs() < 3e-5);
            assert!((result.pressure_drop_pa - 50.0).abs() < 1e-5);
            assert!(result.maximum_heat_imbalance_w < 1e-5);
            assert!(result.neglected_hydraulic_power_w <= 0.01);
        }
    }
    #[test]
    fn cooling_preserves_signed_heat_and_rejects_an_unbracketed_drive() {
        let table = LiquidTable::parse(FIXED).unwrap();
        let result = run(&table, settings(-5.0)).unwrap();
        assert!(result.outlet_k < INLET_K);
        assert!((result.heat_w + 5.0).abs() < 3e-5);
        assert!(run(&table, Settings { drive_pa: 10000.0, ..settings(0.0) }).is_err());
    }

    #[test]
    fn parallel_selection_preserves_the_original_series_command() {
        let series: Vec<_> = ["constant.txt", "50", "20", "0.0002", "0.008"]
            .into_iter().map(str::to_owned).collect();
        let (path, original) = parse_args(&series).unwrap();
        assert_eq!(path, "constant.txt");
        assert_eq!(original, settings(20.0));
        let mut parallel = vec!["--parallel".to_owned()];
        parallel.extend(series.clone());
        assert_eq!(parse_args(&parallel).unwrap().1, Settings { topology: Topology::Parallel, ..original });
        parallel.push("unexpected".to_owned());
        assert!(parse_args(&parallel).is_err());
        assert!(parse_args(&["--parallel".to_owned()]).is_err());
    }
    #[test]
    fn parallel_topology_and_material_tables_reach_actual_fem_and_enthalpy_mixing() {
        let fixed_table = LiquidTable::parse(FIXED).unwrap();
        let variable_table = LiquidTable::parse(VARIABLE).unwrap();
        let settings = Settings { topology: Topology::Parallel, ..settings(20.0) };
        let fixed = run(&fixed_table, settings).unwrap();
        let variable = run(&variable_table, settings).unwrap();
        // Two 5e9/m3 resistances in parallel, not the old series equivalent:
        // each branch carries 0.005 kg/s; total 20 W raises the mixed outlet 1 K.
        assert!((fixed.mass_flow_kg_s - 0.01).abs() < 2e-9);
        assert!((fixed.outlet_k - 301.0).abs() < 1e-5);
        assert!(fixed.branch_mass_flows_kg_s.iter().all(|m| (m - 0.005).abs() < 1e-9));
        assert!(variable.mass_flow_kg_s > fixed.mass_flow_kg_s);
        assert!(variable.outlet_k < fixed.outlet_k);
        assert_ne!(variable.walls_k, fixed.walls_k);
        for result in [&fixed, &variable] {
            assert_eq!(result.branch_mass_flows_kg_s.len(), 2);
            assert_eq!(result.branch_outlets_k.len(), 2);
            assert!((result.mass_flow_kg_s - result.branch_mass_flows_kg_s.iter().sum::<f64>()).abs() < 1e-12);
            assert!((result.heat_w - 20.0).abs() < 3e-5);
            assert!((result.solid_heat_w - 20.0).abs() < 3e-5);
            assert!(result.branch_pressure_drops_pa.iter().all(|p| (p - 50.0).abs() < 1e-5));
            assert!(result.maximum_heat_imbalance_w < 1e-5);
            assert!(result.neglected_hydraulic_power_w <= 0.01);
            assert!(result.solid_solves <= 64);
        }
        let cooled = run(&fixed_table, Settings { heat_w: -5.0, ..settings }).unwrap();
        assert!((cooled.outlet_k - 299.75).abs() < 1e-5);
        assert!((cooled.heat_w + 5.0).abs() < 3e-5);
        assert!((cooled.solid_heat_w + 5.0).abs() < 3e-5);
    }

    #[test]
    fn pump_command_preserves_legacy_forms_and_selects_an_explicit_pressure_bracket() {
        let args: Vec<_> = ["--parallel", "--pump", "pump.txt", "20", "90", "constant.txt", "20", "0.0002", "0.012"]
            .into_iter().map(str::to_owned).collect();
        let (material, settings, pump) = parse_command(&args).unwrap();
        assert_eq!(material, "constant.txt");
        assert_eq!(settings.topology, Topology::Parallel);
        let (name, controls) = pump.unwrap();
        assert_eq!(name, "pump.txt");
        assert_eq!(controls.pressure_bracket_pa, (20.0, 90.0));
        assert_eq!(parse_command(&args[1..]).unwrap().1.topology, Topology::Series);
        assert!(parse_command(&args[..args.len()-1]).is_err());
        let old: Vec<_> = ["constant.txt", "50", "20", "0.0002", "0.008"].into_iter().map(str::to_owned).collect();
        assert!(parse_command(&old).unwrap().2.is_none());
        assert_eq!(parse_command(&old).unwrap().1, parse_args(&old).unwrap().1);
    }

    #[test]
    fn pump_and_material_tables_drive_the_same_series_and_parallel_fem() {
        let curve = pump_source::PumpTable::parse(include_str!("pump.txt")).unwrap();
        let root = PumpOperatingControls { pressure_bracket_pa: (20.0, 90.0),
            absolute_pressure_tolerance_pa: 1e-6, relative_pressure_tolerance: 1e-8,
            max_network_evaluations: 48 };
        let fixed = LiquidTable::parse(FIXED).unwrap();
        let variable = LiquidTable::parse(VARIABLE).unwrap();
        for (topology, expected_pressure, expected_mass) in [
            (Topology::Series, 80.0, 0.004), (Topology::Parallel, 50.0, 0.01),
        ] {
            let setup = Settings { topology, maximum_mass: 0.012, ..settings(20.0) };
            let a = run_with_drive(&fixed, setup, Some((&curve.curve, root))).unwrap();
            let b = run_with_drive(&variable, setup, Some((&curve.curve, root))).unwrap();
            assert!((a.pressure_drop_pa-expected_pressure).abs() < 1e-5);
            assert!((a.mass_flow_kg_s-expected_mass).abs() < 1e-8);
            assert!(b.pressure_drop_pa < a.pressure_drop_pa);
            assert!(b.mass_flow_kg_s > a.mass_flow_kg_s);
            assert_ne!(a.walls_k, b.walls_k);
            for answer in [&a, &b] {
                assert!((answer.heat_w-20.0).abs() < 3e-5);
                assert!((answer.solid_heat_w-20.0).abs() < 3e-5);
                assert!(answer.pump.as_ref().unwrap().residual_pa.abs() < 2e-6);
                assert!(answer.maximum_heat_imbalance_w < 1e-5);
                assert!(answer.neglected_hydraulic_power_w <= 0.01);
                assert!(answer.solid_solves <= 4096);
            }
        }
    }

    #[test]
    fn pump_table_substitution_and_signed_heat_are_not_fixed_pressure_shortcuts() {
        let curve = pump_source::PumpTable::parse(include_str!("pump.txt")).unwrap();
        let stronger = pump_source::PumpTable::parse(&include_str!("pump.txt").replace("0 0.01", "0 0.02")).unwrap();
        let root = PumpOperatingControls { pressure_bracket_pa: (20.0, 75.0),
            absolute_pressure_tolerance_pa: 1e-6, relative_pressure_tolerance: 1e-8,
            max_network_evaluations: 48 };
        let table = LiquidTable::parse(FIXED).unwrap();
        let setup = Settings { topology: Topology::Parallel, maximum_mass: 0.012, ..settings(20.0) };
        let a = run_with_drive(&table, setup, Some((&curve.curve, root))).unwrap();
        let b = run_with_drive(&table, setup, Some((&stronger.curve, root))).unwrap();
        assert!(b.mass_flow_kg_s > a.mass_flow_kg_s);
        assert!(b.pressure_drop_pa > a.pressure_drop_pa);
        assert!(b.outlet_k < a.outlet_k);
        let cooled = run_with_drive(&table, Settings { heat_w: -5.0, ..setup }, Some((&curve.curve, root))).unwrap();
        assert!(cooled.outlet_k < INLET_K);
        assert!((cooled.heat_w+5.0).abs() < 3e-5);
    }

}
