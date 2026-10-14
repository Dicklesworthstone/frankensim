//! Material table -> pressure-selected flow -> coolant enthalpy <-> actual FEM.
//! This example owns no material law, heat equation, or coupling algorithm.
mod source;

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
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    InitialGuess, LinearConfig, Nonlinearity, ScalarField, SolveConfig, StopRule,
    ThermalBoundaryBuilder, solve};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_convection::CorrelationId;
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

#[derive(Clone, Copy)]
struct Settings { drive_pa: f64, heat_w: f64, minimum_mass: f64, maximum_mass: f64 }

#[derive(Debug)]
struct CoolingResult {
    mass_flow_kg_s: f64,
    outlet_k: f64,
    walls_k: Vec<f64>,
    heat_w: f64,
    solid_heat_w: f64,
    pressure_drop_pa: f64,
    neglected_hydraulic_power_w: f64,
    maximum_heat_imbalance_w: f64,
    solid_solves: usize,
}

fn run(table: &LiquidTable, settings: Settings) -> Result<CoolingResult> {
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
        let path = LiquidHeatPath::try_new(&cx, &law, &caloric, MassFlowRate::new(settings.minimum_mass),
            Temperature::new(INLET_K), &sections, LiquidPropertyTemperature::BulkMean)?;
        let mut solid_error = None;
        let answer = path.solve_pressure_coupled_iqn(&cx, &[INLET_K; 2], &resistances,
            Pressure::new(settings.drive_pa), &controls, |cx, rows| {
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
            }).map_err(|error| match solid_error {
                Some(original) => format!("{error}; solid: {original}"),
                None => error.to_string(),
            })?;
        let point = answer.operating_point;
        let coupled = point.thermal;
        let outlet = coupled.fluid.regions.last().ok_or("accepted path lacks an outlet")?;
        Ok(CoolingResult {
            mass_flow_kg_s: point.mass_flow_kg_per_s, outlet_k: outlet.outlet.temperature_k(),
            walls_k: coupled.solid.regions.iter().map(|r| r.mean_wall_temperature_k).collect(),
            heat_w: coupled.fluid.total_heat_rate_w, solid_heat_w: coupled.solid.state.report.energy.robin_out_w,
            pressure_drop_pa: point.hydraulics.total_pressure_drop_pa,
            neglected_hydraulic_power_w: point.hydraulics.total_dissipated_power_w,
            maximum_heat_imbalance_w: coupled.maximum_heat_imbalance_w, solid_solves: answer.total_solid_solves,
        })
    })
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 5 {
        return Err("usage: pressure_cooling MATERIAL.txt DRIVE_PA BODY_HEAT_W MIN_MASS_KG_S MAX_MASS_KG_S".into());
    }
    let table = LiquidTable::parse(&std::fs::read_to_string(&args[0])?)?;
    let result = run(&table, Settings { drive_pa: args[1].parse()?, heat_w: args[2].parse()?,
        minimum_mass: args[3].parse()?, maximum_mass: args[4].parse()? })?;
    println!("model=regional_liquid_pressure_cooling; phase_stability=assumed; calibration=not_established");
    println!("material_card={:?}", table.card.content_hash());
    println!("mass_flow_kg_s={:.12e}\noutlet_k={:.9}\nwall_temperatures_k={:?}", result.mass_flow_kg_s, result.outlet_k, result.walls_k);
    println!("coolant_heat_w={:.9}\nsolid_robin_heat_w={:.9}\npressure_drop_pa={:.9}", result.heat_w, result.solid_heat_w, result.pressure_drop_pa);
    println!("neglected_hydraulic_power_w={:.12e}\nmaximum_regional_heat_imbalance_w={:.12e}\nsolid_solves={}",
        result.neglected_hydraulic_power_w, result.maximum_heat_imbalance_w, result.solid_solves);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const FIXED: &str = include_str!("constant.txt");
    const VARIABLE: &str = include_str!("temperature_dependent.txt");
    fn settings(heat_w: f64) -> Settings {
        Settings { drive_pa: 50.0, heat_w, minimum_mass: 0.0002, maximum_mass: 0.008 }
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
}
