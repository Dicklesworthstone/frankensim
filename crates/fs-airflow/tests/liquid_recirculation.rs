//! Finite-inventory recirculation through the real material/pump/flow owners.
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
    PumpHeadCurve, PumpHeadPoint, PumpOperatingControls,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::recirculation::{
    RecirculationError, ReservoirStepControls,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::CorrelationId;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, QuantitySpec, Temperature, VolumetricFlowRate};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../examples/pressure_cooling/source.rs"]
mod source;
const FIXED: &str = include_str!("../examples/pressure_cooling/constant.txt");
const VARIABLE: &str = include_str!("../examples/pressure_cooling/temperature_dependent.txt");

fn with_network<R>(text: &str, run: impl FnOnce(&Cx<'_>, &ParallelLiquidNetwork<'_, '_, '_>) -> R) -> R {
    with_network_gate(text, |cx, network, _| run(cx, network))
}
fn with_network_gate<R>(text: &str, run: impl FnOnce(&Cx<'_>, &ParallelLiquidNetwork<'_, '_, '_>, &CancelGate) -> R) -> R {
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
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 947, kernel_id: 754, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let geometry = FlowHeatGeometry { flow_area: Area::new(core::f64::consts::PI*0.01*0.01/4.0),
            characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
            wetted_area: Area::new(0.003), aspect_ratio: None };
        let names = ["first", "last"];
        let paths: Vec<_> = names.iter().map(|name| LiquidHeatPath::try_new(&cx, &law, &caloric,
            MassFlowRate::new(0.001), Temperature::new(300.0),
            &[LiquidHeatSection { region: name, geometry, correlation: CorrelationId::CircularDuctLaminarCwt }],
            LiquidPropertyTemperature::BulkMean).unwrap()).collect();
        let losses: Vec<_> = names.iter().map(|name| [LiquidResistanceSection { region: name,
            resistance: MaterialFlowResistance::try_new(5e9, 0.0, geometry.flow_area,
                geometry.characteristic_length, (0.0,2300.0)).unwrap() }]).collect();
        let branches: Vec<_> = (0..2).map(|i| ParallelLiquidBranch { name: names[i], path: &paths[i],
            resistances: &losses[i], pressure_budget: PressureDrivenBudget {
                mass_flow_bracket_kg_per_s: (0.0002,0.012), pressure_absolute_tolerance_pa: 1e-10,
                pressure_relative_tolerance: 1e-10, max_flow_evaluations: 64,
                maximum_neglected_hydraulic_power_w: 0.01,
            } }).collect();
        let network = ParallelLiquidNetwork::try_new(&cx, &branches).unwrap();
        run(&cx, &network, &gate)
    })
}
fn pump() -> PumpHeadCurve {
    PumpHeadCurve::try_new(vec![PumpHeadPoint { flow: VolumetricFlowRate::new(0.0), head: Length::new(0.005) },
        PumpHeadPoint { flow: VolumetricFlowRate::new(4e-5), head: Length::new(0.005) }],
        10.0, (1e-7,1e-4), fs_airflow::SourceProvenance::new("synthetic constant head", "reservoir-test")).unwrap()
}
fn local() -> LiquidFeedbackBudget {
    LiquidFeedbackBudget { temperature_tolerance_k: 1e-10, coefficient_relative_tolerance: 1e-10, max_region_iterations: 64 }
}
fn work() -> ParallelPressureBudget {
    ParallelPressureBudget { max_flow_evaluations: 1_000_000, maximum_neglected_hydraulic_power_w: 0.01 }
}
fn operating() -> PumpOperatingControls {
    PumpOperatingControls { pressure_bracket_pa: (40.0,60.0), absolute_pressure_tolerance_pa: 1e-9,
        relative_pressure_tolerance: 1e-10, max_network_evaluations: 64 }
}
fn step() -> ReservoirStepControls {
    ReservoirStepControls { temperature_bracket_k: (290.0,340.0), energy_absolute_tolerance_j: 1e-6,
        max_inlet_evaluations: 64, maximum_neglected_hydraulic_energy_j: 0.1 }
}

#[test]
fn finite_inventory_matches_implicit_heat_balance_and_preserves_exact_inlet() {
    with_network(FIXED, |cx, network| {
        let pump = pump();
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let result = network.advance_reservoir_with_pump(cx,&mut tank,10.0,0.0,&[380.0;2],
            &pump,local(),work(),operating(),step()).unwrap();
        let conductance = 2.0*10.0*(-(-183.0_f64*0.003/10.0).exp_m1());
        let expected = (200.0*300.0 + 10.0*conductance*380.0)/(200.0+10.0*conductance);
        assert!((tank.state().temperature_k()-expected).abs() < 1e-7);
        assert_eq!(tank.time_s(),10.0);
        assert_eq!(result.final_state,tank.state());
        assert!(result.energy_residual_j.abs() <= step().energy_absolute_tolerance_j);
        assert!(result.returned_heat_j > 0.0);
        assert!(result.inlet_evaluations <= 6, "linear residual should use bracketed secant");
        for branch in &result.circuit.state.branches {
            assert_eq!(branch.solution.thermal.regions[0].inlet,tank.state());
        }
        // The source description remains reusable; it has not become mutable history.
        assert_eq!(network.branches()[0].path.inlet().temperature_k(),300.0);
    });
}

#[test]
fn source_heat_capacity_and_transport_changes_reach_reservoir_and_pump() {
    let run = |text| with_network(text, |cx, network| {
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let pump = pump();
        let answer = network.advance_reservoir_with_pump(cx,&mut tank,10.0,100.0,&[380.0;2],
            &pump,local(),work(),operating(),step()).unwrap();
        assert!((answer.stored_enthalpy_change_j-answer.returned_heat_j-100.0).abs() < 1e-6);
        assert!((answer.circuit.supply.inlet_density_kg_m3 - answer.circuit.state.inlet_material.density_kg_m3()).abs() < 1e-12);
        (tank.state().temperature_k(),answer.circuit.state.total_mass_flow_kg_per_s)
    });
    let fixed = run(FIXED);
    let variable = run(VARIABLE);
    assert!((fixed.0-variable.0).abs() > 1e-3);
    assert!((fixed.1-variable.1).abs() > 1e-6);
}

#[test]
fn cooling_and_external_reservoir_heat_are_signed_and_counted_once() {
    with_network(FIXED, |cx, network| {
        let pump = pump();
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),2.0).unwrap();
        let initial = tank.clone();
        let answer = network.advance_reservoir_with_pump(cx,&mut tank,12.0,-100.0,&[280.0;2],
            &pump,local(),work(),operating(),step()).unwrap();
        assert!(tank.state().temperature_k() < initial.state().temperature_k());
        assert!(answer.returned_heat_j < 0.0);
        assert!((answer.stored_enthalpy_change_j-answer.returned_heat_j+100.0).abs() < 1e-6);
        assert_eq!(answer.start_time_s,2.0);
        assert_eq!(answer.end_time_s,12.0);
        let mut replay = initial;
        network.advance_reservoir_with_pump(cx,&mut replay,12.0,-100.0,&[280.0;2],
            &pump,local(),work(),operating(),step()).unwrap();
        assert_eq!(tank,replay);
    });
}

#[test]
fn bad_source_bracket_energy_or_total_work_never_advances_inventory() {
    with_network(FIXED, |cx, network| {
        let pump = pump();
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let original = tank.clone();
        for control in [ReservoirStepControls { temperature_bracket_k: (290.0,291.0), ..step() },
            ReservoirStepControls { temperature_bracket_k: (100.0,340.0), ..step() },
            ReservoirStepControls { maximum_neglected_hydraulic_energy_j: 0.0, ..step() },
            ReservoirStepControls { max_inlet_evaluations: 2, ..step() }] {
            assert!(network.advance_reservoir_with_pump(cx,&mut tank,10.0,0.0,&[380.0;2],
                &pump,local(),work(),operating(),control).is_err());
            assert_eq!(tank,original);
        }
        assert!(network.advance_reservoir_with_pump(cx,&mut tank,10.0,0.0,&[380.0;2],
            &pump,local(),ParallelPressureBudget { max_flow_evaluations: 2, ..work() },operating(),step()).is_err());
        assert!(matches!(network.advance_reservoir_with_pump(cx,&mut tank,0.0,0.0,&[380.0;2],
            &pump,local(),work(),operating(),step()),Err(RecirculationError::Invalid(_))));
        assert_eq!(tank,original);
    });
}

use fs_airflow::conjugate::{ConjugateConfig, Relaxation, SolidRegionState};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::solid::ParallelLiquidCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::solid::PumpCouplingControls;
use fs_convection::CorrelationRobinBoundary;

fn coupling() -> PumpCouplingControls {
    PumpCouplingControls {
        network: ParallelLiquidCouplingControls {
            coupling: ConjugateConfig { temperature_tolerance_k: 1e-8, max_iterations: 64,
                balance_tolerance_w: 1e-7, balance_relative_tolerance: 1e-9,
                relaxation: Relaxation::Fixed { omega: 1.0 } },
            local: local(), acceleration: Default::default(), mass_flow_relative_tolerance: 1e-8,
            work: work(),
        }, operating_point: operating(), max_total_solid_solves: 4096,
    }
}
// Independent two-node steady solid with a 0.3 W/K internal thermal link.
fn body(rows: &[CorrelationRobinBoundary], powers: [f64;2]) -> LiquidSolidResponse<[f64;2]> {
    let k: Vec<_> = rows.iter().map(|r|r.coefficient().value.value()*0.003).collect();
    let t: Vec<_> = rows.iter().map(|r|r.reference_temperature().value()).collect();
    let g = 0.3;
    let determinant = k[0]*k[1]+g*(k[0]+k[1]);
    let rhs = [powers[0]+k[0]*t[0],powers[1]+k[1]*t[1]];
    let walls = [(rhs[0]*(k[1]+g)+g*rhs[1])/determinant,
        (rhs[1]*(k[0]+g)+g*rhs[0])/determinant];
    LiquidSolidResponse { regions: (0..2).map(|i|SolidRegionState {
        region: ["first","last"][i].to_owned(), area_m2: 0.003,
        mean_wall_temperature_k: walls[i], mean_reference_temperature_k: Some(t[i]),
        heat_rate_w: k[i]*(walls[i]-t[i]),
    }).collect(), state: walls }
}

#[test]
fn common_solid_power_drives_inventory_with_exact_returned_field_and_call_counts() {
    with_network(FIXED, |cx, network| {
        let pump = pump();
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let mut calls = 0;
        let result = network.advance_reservoir_with_pump_coupled_iqn(cx,&mut tank,10.0,-50.0,&[325.0;2],
            &pump,&coupling(),step(),|_,rows| { calls+=1; Ok(body(rows,[8.0,12.0])) }).unwrap();
        // (20 W * 10 s - 50 J) / (0.1 kg * 2000 J/kg/K) = 0.75 K.
        assert!((tank.state().temperature_k()-300.75).abs()<1e-7);
        assert!((result.returned_heat_j-200.0).abs()<1e-6);
        assert!((result.stored_enthalpy_change_j-150.0).abs()<1e-6);
        assert_eq!(calls,result.total_solid_solves);
        let state = &result.circuit.state;
        for (i,branch) in state.boundary_state.branches.iter().enumerate() {
            let row = &branch.solution.thermal.regions[0].boundary;
            assert_eq!(state.solid.regions[i].mean_reference_temperature_k,Some(row.reference_temperature().value()));
            assert_eq!(state.solid.state[i],state.solid.regions[i].mean_wall_temperature_k);
        }
        assert!(result.energy_residual_j.abs()<=1e-6);
    });
}

#[test]
fn repeated_variable_material_steps_keep_enthalpy_history_not_a_reset_inlet() {
    with_network(VARIABLE, |cx, network| {
        let pump = pump();
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let h0 = tank.state().specific_enthalpy_j_kg();
        let mut flows = Vec::new();
        for end in [5.0,10.0] {
            let result = network.advance_reservoir_with_pump_coupled_iqn(cx,&mut tank,end,0.0,&[325.0;2],
                &pump,&coupling(),step(),|_,rows|Ok(body(rows,[8.0,12.0]))).unwrap();
            assert!((0.1*(tank.state().specific_enthalpy_j_kg()-h0)-20.0*end).abs()<3e-6);
            flows.push(result.circuit.state.fluid.total_mass_flow_kg_per_s);
        }
        assert!(flows[1]>flows[0]);
        assert_eq!(tank.time_s(),10.0);
    });
}

#[test]
fn regional_power_failure_cancellation_and_shared_solid_budget_roll_back_inventory() {
    with_network_gate(FIXED, |cx, network, gate| {
        let pump = pump();
        let mut tank = network.reservoir(cx,0.1,Temperature::new(300.0),0.0).unwrap();
        let original=tank.clone();
        assert!(network.advance_reservoir_with_pump_coupled_iqn(cx,&mut tank,10.0,0.0,&[325.0;2],
            &pump,&coupling(),step(),|_,rows| {
                let mut body=body(rows,[8.0,12.0]);
                body.regions[0].heat_rate_w+=0.1; body.regions[1].heat_rate_w-=0.1;
                Ok(body)
            }).is_err());
        assert_eq!(tank,original);
        let mut calls=0;
        let tiny=PumpCouplingControls { max_total_solid_solves:1,..coupling() };
        assert!(network.advance_reservoir_with_pump_coupled_iqn(cx,&mut tank,10.0,0.0,&[325.0;2],
            &pump,&tiny,step(),|_,rows|{calls+=1;Ok(body(rows,[8.0,12.0]))}).is_err());
        assert_eq!(calls,1);
        assert_eq!(tank,original);
        assert!(network.advance_reservoir_with_pump_coupled_iqn(cx,&mut tank,10.0,0.0,&[325.0;2],
            &pump,&coupling(),step(),|_,rows|{gate.request();Ok(body(rows,[8.0,12.0]))}).is_err());
        assert_eq!(tank,original);
    });
}
