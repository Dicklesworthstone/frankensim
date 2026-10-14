//! Synthetic pump/material data through the production cooling-network owners.
use fs_airflow::SourceProvenance;
use fs_airflow::graph::thermal::medium::FlowHeatGeometry;
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{
    LiquidResistanceSection, PressureDrivenBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::{
    ParallelLiquidBranch, ParallelLiquidNetwork, ParallelPressureBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpHeadCurve, PumpHeadPoint, PumpOperatingControls, PumpOperatingError,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::CorrelationId;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, Pressure, QuantitySpec, Temperature, VolumetricFlowRate};

// Exercise the same material ingress as the runnable cooling example.
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../examples/pressure_cooling/source.rs"]
mod source;
const FIXED: &str = include_str!("../examples/pressure_cooling/constant.txt");
const VARIABLE: &str = include_str!("../examples/pressure_cooling/temperature_dependent.txt");
fn pump(flat: bool) -> PumpHeadCurve {
    PumpHeadCurve::try_new(vec![
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.0), head: Length::new(if flat { 0.005 } else { 0.01 }) },
        PumpHeadPoint { flow: VolumetricFlowRate::new(1e-5), head: Length::new(if flat { 0.005 } else { 0.0 }) },
    ], 10.0, (1e-6, 4e-6), SourceProvenance::new("synthetic fixed-speed pump, not measurements", "pump-test-v1")).unwrap()
}
fn root() -> PumpOperatingControls {
    PumpOperatingControls { pressure_bracket_pa: (20.0, 70.0), absolute_pressure_tolerance_pa: 1e-7,
        relative_pressure_tolerance: 1e-9, max_network_evaluations: 40 }
}
fn work() -> ParallelPressureBudget {
    ParallelPressureBudget { max_flow_evaluations: 10000, maximum_neglected_hydraulic_power_w: 0.01 }
}
fn local() -> LiquidFeedbackBudget {
    LiquidFeedbackBudget { temperature_tolerance_k: 1e-9, coefficient_relative_tolerance: 1e-9, max_region_iterations: 64 }
}
fn with_network<R>(text: &str, run: impl FnOnce(&Cx<'_>, &ParallelLiquidNetwork<'_, '_, '_>) -> R) -> R {
    let table = source::LiquidTable::parse(text).unwrap();
    let law = LiquidTemperatureLaw::try_new(&table.card, &table.point(300.0).unwrap(), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let cp = PropertyKey::with_quantity("specific_heat_capacity", QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density", QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    let caloric = ResolvedHeatCapacityCurve::try_new(&table.card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T",
        lower: &table.point(table.lower_k).unwrap(), upper: &table.point(table.upper_k).unwrap(),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Liquid, 0.0, 0.5).unwrap();
    let geometry = FlowHeatGeometry { flow_area: Area::new(core::f64::consts::PI*0.01*0.01/4.0),
        characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
        wetted_area: Area::new(0.003), aspect_ratio: None };
    let resistance = MaterialFlowResistance::try_new(1e10, 0.0, geometry.flow_area,
        geometry.characteristic_length, (0.0,2300.0)).unwrap();
    let branch_budget = PressureDrivenBudget { mass_flow_bracket_kg_per_s: (0.0002,0.008),
        pressure_absolute_tolerance_pa: 1e-10, pressure_relative_tolerance: 1e-11,
        max_flow_evaluations: 64, maximum_neglected_hydraulic_power_w: 0.01 };
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate,arena,StreamKey { seed: 91,kernel_id: 744,tile: 0,iteration: 0 },
            Budget::INFINITE,ExecMode::Deterministic);
        let make = |name| LiquidHeatPath::try_new(&cx,&law,&caloric,MassFlowRate::new(0.001),
            Temperature::new(300.0), &[LiquidHeatSection { region: name,geometry,
                correlation: CorrelationId::CircularDuctLaminarCwt }],LiquidPropertyTemperature::BulkMean).unwrap();
        let a = make("first"); let b = make("last");
        let ar = [LiquidResistanceSection { region: "first",resistance }];
        let br = [LiquidResistanceSection { region: "last",resistance }];
        let network = ParallelLiquidNetwork::try_new(&cx,&[
            ParallelLiquidBranch { name: "first",path: &a,resistances: &ar,pressure_budget: branch_budget },
            ParallelLiquidBranch { name: "last",path: &b,resistances: &br,pressure_budget: branch_budget },
        ]).unwrap();
        run(&cx,&network)
    })
}
fn close(a: f64,b: f64,tol: f64) { assert!((a-b).abs() <= tol,"{a} != {b}"); }

#[test]
fn pump_intersection_matches_analytic_pressure_flow_and_actual_heat_state() {
    with_network(FIXED,|cx,network| {
        let curve = pump(false);
        let result = network.solve_with_pump(cx,&[400.0,280.0],&curve,local(),work(),root()).unwrap();
        // Q_total = dp/1e7; pump dp = 100 - 1e7*Q_total => dp = 50 Pa.
        close(result.network_pressure_pa,50.0,1e-6);
        close(result.state.total_mass_flow_kg_per_s,0.005,1e-10);
        close(result.supply.available_pressure_pa,50.0,1e-6);
        assert!(result.pressure_residual_pa.abs() < 1e-6);
        let fixed = network.solve_pressure_driven(cx,&[400.0,280.0],Pressure::new(result.network_pressure_pa),local(),work()).unwrap();
        assert_eq!(result.state.mixed_outlet,fixed.mixed_outlet);
        close(result.state.total_heat_rate_w,fixed.total_heat_rate_w,1e-10);
        close(result.supply.hydraulic_power_w,result.state.total_dissipated_power_w,1e-10);
        assert!(result.total_flow_evaluations > result.state.flow_evaluations);
        assert_eq!(result.total_solid_solves,0);
        assert_eq!(result.curve,&curve);
    });
}

#[test]
fn material_temperature_changes_move_the_pump_pressure_and_the_branch_split() {
    with_network(VARIABLE,|cx,network| {
        let curve = pump(false);
        let cold = network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),work(),root()).unwrap();
        let hot = network.solve_with_pump(cx,&[400.0,300.0],&curve,local(),work(),root()).unwrap();
        assert!(hot.network_pressure_pa < cold.network_pressure_pa);
        assert!(hot.state.total_mass_flow_kg_per_s > cold.state.total_mass_flow_kg_per_s);
        assert!(hot.state.branches[0].solution.mass_flow_kg_per_s > hot.state.branches[1].solution.mass_flow_kg_per_s);
        close(hot.supply.inlet_volume_flow_m3_s,hot.state.total_mass_flow_kg_per_s/hot.state.inlet_material.density_kg_m3(),1e-18);
        let reference_h = network.branches()[0].path.inlet().specific_enthalpy_j_kg();
        close(hot.state.total_heat_rate_w, hot.state.total_mass_flow_kg_per_s*
            (hot.state.mixed_outlet.specific_enthalpy_j_kg()-reference_h),1e-7);
    });
}

#[test]
fn head_uses_changed_inlet_density_and_a_flat_curve_matches_constant_pressure() {
    let denser = FIXED.replace("250 1000", "250 1200").replace("450 1000", "450 1200");
    with_network(&denser,|cx,network| {
        let curve = pump(true);
        let result = network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),work(),root()).unwrap();
        close(result.supply.available_pressure_pa,60.0,1e-12);
        close(result.network_pressure_pa,60.0,1e-6);
        close(result.state.total_mass_flow_kg_per_s,0.0072,1e-9);
        close(result.state.inlet_volume_flow_m3_per_s,6e-6,1e-12);
    });
}

#[test]
fn bad_curves_unsupported_viscosity_and_flow_never_extrapolate() {
    let curve = pump(false);
    let mut reversed = curve.points().to_vec(); reversed.reverse();
    assert!(PumpHeadCurve::try_new(reversed,10.0,(1e-6,4e-6),curve.source().clone()).is_err());
    assert!(PumpHeadCurve::try_new(curve.points().to_vec(),0.0,(1e-6,4e-6),curve.source().clone()).is_err());
    with_network(FIXED,|cx,network| {
        let state = network.solve_pressure_driven(cx,&[300.0,300.0],Pressure::new(50.0),local(),work()).unwrap();
        assert!(matches!(curve.evaluate(&state.inlet_material,VolumetricFlowRate::new(2e-5)),Err(PumpOperatingError::OutsideCurve { .. })));
        let unsupported = PumpHeadCurve::try_new(curve.points().to_vec(),10.0,(3e-6,4e-6),curve.source().clone()).unwrap();
        assert!(matches!(network.solve_with_pump(cx,&[300.0,300.0],&unsupported,local(),work(),root()),Err(PumpOperatingError::OutsideCurve { .. })));
        assert_eq!(curve.evaluate(&state.inlet_material,VolumetricFlowRate::new(1e-5)).unwrap().available_pressure_pa,0.0);
    });
}

#[test]
fn pump_matching_preserves_whole_network_budgets_and_explicit_brackets() {
    with_network(FIXED,|cx,network| {
        let curve = pump(false);
        assert!(matches!(network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),work(),
            PumpOperatingControls { pressure_bracket_pa: (20.0,30.0),..root() }),Err(PumpOperatingError::NotBracketed { .. })));
        assert!(matches!(network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),work(),
            PumpOperatingControls { max_network_evaluations: 2,..root() }),Err(PumpOperatingError::NotConverged { evaluations: 2,.. })));
        assert!(network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),
            ParallelPressureBudget { max_flow_evaluations: 2,..work() },root()).is_err());
        assert!(network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),
            ParallelPressureBudget { maximum_neglected_hydraulic_power_w: 0.0,..work() },root()).is_err());
        let a = network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),work(),root()).unwrap();
        let b = network.solve_with_pump(cx,&[300.0,300.0],&curve,local(),work(),root()).unwrap();
        assert_eq!(a.supply,b.supply);
        assert_eq!(a.state.mixed_outlet,b.state.mixed_outlet);
    });
}
