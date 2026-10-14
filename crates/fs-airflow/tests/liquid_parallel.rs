//! Synthetic material data through the actual branch pressure/thermal owners.
use fs_airflow::graph::thermal::medium::FlowHeatGeometry;
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{
    LiquidResistanceSection, PressureDrivenBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::{
    ParallelLiquidBranch, ParallelLiquidNetwork, ParallelPressureBudget, ParallelPressureError,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::CorrelationId;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn geometry() -> FlowHeatGeometry {
    FlowHeatGeometry { flow_area: Area::new(core::f64::consts::PI*0.01*0.01/4.0),
        characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
        wetted_area: Area::new(0.003), aspect_ratio: None }
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T",t).unwrap().with("pressure",100000.0).unwrap()
}
fn material(variable_mu: bool, variable_cp: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density",1000.0,Dims([-3,1,0,0,0,0])),
        ("specific_heat_capacity",2000.0,Dims([2,0,-2,-1,0,0])),
        ("thermal_conductivity",0.5,Dims([1,1,-3,-1,0,0])),
        ("dynamic_viscosity",0.002,Dims([-1,1,-1,0,0,0])),
    ] {
        let knots = match name {
            "dynamic_viscosity" if variable_mu => Some(vec![(250.0,0.0025),(450.0,0.0005)]),
            "specific_heat_capacity" if variable_cp => Some(vec![(250.0,1000.0),(450.0,5000.0)]),
            _ => None,
        };
        let varying = knots.is_some();
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity,ValueForm::Static))
        } else { QuantitySpec::dimensional(dims) };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name,quantity),
            value: match knots {
                Some(knots) => PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS, knots, dims },
                None => PropertyValue::Scalar { value,dims },
            },
            validity: ValidityDomain::unconstrained().with("T",250.0,450.0).with("pressure",100000.0,100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if varying { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "synthetic parallel-liquid regression, not measurements".into(),
                license: "test-only".into(),artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic-liquid".into(), phase: "liquid".into(),
        process: "numerical fixture".into(),revision: 0 },claims,vec![]).unwrap()
}
fn with_source<R>(variable_mu: bool, variable_cp: bool,
    run: impl FnOnce(&Cx<'_>, &LiquidTemperatureLaw<'_>, &ResolvedHeatCapacityCurve) -> R) -> R {
    let card = material(variable_mu,variable_cp);
    let law = LiquidTemperatureLaw::try_new(&card,&point(300.0),"T",LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let cp = PropertyKey::with_quantity("specific_heat_capacity",QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density",QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    let source = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity: &cp,density: &rho,temperature_axis: "T",lower: &point(250.0),upper: &point(450.0),
        selection: SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Liquid,0.0,0.5).unwrap();
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(&gate,arena,StreamKey { seed: 93,kernel_id: 742,tile: 0,iteration: 0 },
            Budget::INFINITE,ExecMode::Deterministic),&law,&source)
    })
}
fn path<'a,'c>(cx: &Cx<'_>, law: &'a LiquidTemperatureLaw<'c>, source: &'a ResolvedHeatCapacityCurve,
    name: &str) -> LiquidHeatPath<'a,'c> {
    LiquidHeatPath::try_new(cx,law,source,MassFlowRate::new(0.001),Temperature::new(300.0),
        &[LiquidHeatSection { region: name,geometry: geometry(),correlation: CorrelationId::CircularDuctLaminarCwt }],
        LiquidPropertyTemperature::BulkMean).unwrap()
}
fn loss(name: &str, scale: f64) -> [LiquidResistanceSection<'_>;1] {
    [LiquidResistanceSection { region: name, resistance: MaterialFlowResistance::try_new(
        scale*1e10,0.0,geometry().flow_area,Length::new(0.01),(0.0,2300.0)).unwrap() }]
}
fn local() -> LiquidFeedbackBudget {
    LiquidFeedbackBudget { temperature_tolerance_k: 1e-9,coefficient_relative_tolerance: 1e-9,max_region_iterations: 64 }
}
fn branch_budget() -> PressureDrivenBudget {
    PressureDrivenBudget { mass_flow_bracket_kg_per_s: (0.0002,0.008),pressure_absolute_tolerance_pa: 1e-9,
        pressure_relative_tolerance: 1e-9,max_flow_evaluations: 64,maximum_neglected_hydraulic_power_w: 0.01 }
}
fn budget() -> ParallelPressureBudget {
    ParallelPressureBudget { max_flow_evaluations: 10000, maximum_neglected_hydraulic_power_w: 0.02 }
}
fn near(a: f64,b: f64) { assert!((a-b).abs() < 3e-8*(1.0+a.abs().max(b.abs())),"{a} != {b}"); }
fn with_network<R>(variable_mu: bool, variable_cp: bool, second_scale: f64,
    run: impl FnOnce(&Cx<'_>, &ParallelLiquidNetwork<'_, '_, '_>) -> R) -> R {
    with_source(variable_mu,variable_cp,|cx,law,source| {
        let a = path(cx,law,source,"a-wall");
        let b = path(cx,law,source,"b-wall");
        let ar = loss("a-wall",1.0);
        let br = loss("b-wall",second_scale);
        let network = ParallelLiquidNetwork::try_new(cx,&[
            ParallelLiquidBranch { name: "a",path: &a,resistances: &ar,pressure_budget: branch_budget() },
            ParallelLiquidBranch { name: "b",path: &b,resistances: &br,pressure_budget: branch_budget() },
        ]).unwrap();
        run(cx,&network)
    })
}

#[test]
fn parallel_constant_limit_matches_independent_pressure_and_heat_equations() {
    with_network(false,false,2.0,|cx,network| {
        let walls = [400.0,280.0];
        let solved = network.solve_pressure_driven(cx,&walls,Pressure::new(50.0),local(),budget()).unwrap();
        let expected_mass = [0.0025,0.00125];
        let mut expected_heat = 0.0;
        for ((branch,&mass),&wall) in solved.branches.iter().zip(&expected_mass).zip(&walls) {
            near(branch.solution.mass_flow_kg_per_s,mass);
            near(branch.solution.hydraulics.total_pressure_drop_pa,50.0);
            let capacity = mass*2000.0;
            let outlet = 300.0+(wall-300.0)*(-(-183.0*0.003/capacity).exp_m1());
            near(branch.solution.thermal.regions[0].outlet.temperature_k(),outlet);
            expected_heat += capacity*(outlet-300.0);
        }
        near(solved.total_mass_flow_kg_per_s,0.00375);
        near(solved.inlet_volume_flow_m3_per_s,0.00375/1000.0);
        near(solved.total_heat_rate_w,expected_heat);
        near(solved.mixed_outlet.temperature_k(),300.0+expected_heat/(0.00375*2000.0));
        assert_eq!(solved.flow_evaluations,solved.branches.iter().map(|b| b.solution.flow_evaluations).sum::<usize>());
    });
}

#[test]
fn heating_redistributes_flow_without_changing_geometry_or_the_common_source() {
    with_network(true,false,1.0,|cx,network| {
        let cold = network.solve_pressure_driven(cx,&[300.0,300.0],Pressure::new(50.0),local(),budget()).unwrap();
        let hot = network.solve_pressure_driven(cx,&[400.0,300.0],Pressure::new(50.0),local(),budget()).unwrap();
        near(cold.branches[0].solution.mass_flow_kg_per_s,cold.branches[1].solution.mass_flow_kg_per_s);
        assert!(hot.branches[0].solution.mass_flow_kg_per_s > hot.branches[1].solution.mass_flow_kg_per_s);
        near(hot.branches[1].solution.mass_flow_kg_per_s,cold.branches[1].solution.mass_flow_kg_per_s);
        assert!(hot.total_heat_rate_w > 0.0);
        assert_eq!(cold.inlet_material,hot.inlet_material);
        let unchanged = network.solve_pressure_driven(cx,&[300.0,300.0],Pressure::new(50.0),local(),budget()).unwrap();
        assert_eq!(cold.mixed_outlet,unchanged.mixed_outlet);
    });
}

#[test]
fn mixing_variable_heat_capacity_conserves_enthalpy_not_average_temperature() {
    with_network(false,true,2.0,|cx,network| {
        let solved = network.solve_pressure_driven(cx,&[440.0,270.0],Pressure::new(50.0),local(),budget()).unwrap();
        let inlet_h = network.branches()[0].path.inlet().specific_enthalpy_j_kg();
        let expected_h = inlet_h+solved.branches.iter().map(|b| {
            b.solution.mass_flow_kg_per_s*(b.solution.thermal.regions[0].outlet.specific_enthalpy_j_kg()-inlet_h)
        }).sum::<f64>()/solved.total_mass_flow_kg_per_s;
        near(solved.mixed_outlet.specific_enthalpy_j_kg(),expected_h);
        let mean_temperature = solved.branches.iter().map(|b| b.solution.mass_flow_kg_per_s
            * b.solution.thermal.regions[0].outlet.temperature_k()).sum::<f64>()/solved.total_mass_flow_kg_per_s;
        assert!((solved.mixed_outlet.temperature_k()-mean_temperature).abs() > 1e-3);
    });
}

#[test]
fn branch_order_changes_addresses_not_physics() {
    with_network(true,true,2.0,|cx,network| {
        let original = network.solve_pressure_driven(cx,&[400.0,280.0],Pressure::new(50.0),local(),budget()).unwrap();
        let branches = network.branches();
        let swapped = ParallelLiquidNetwork::try_new(cx,&[branches[1],branches[0]]).unwrap();
        let actual = swapped.solve_pressure_driven(cx,&[280.0,400.0],Pressure::new(50.0),local(),budget()).unwrap();
        near(actual.mixed_outlet.temperature_k(),original.mixed_outlet.temperature_k());
        near(actual.total_mass_flow_kg_per_s,original.total_mass_flow_kg_per_s);
        assert_eq!(actual.branches[0].name,"b");
        near(actual.branches[0].solution.mass_flow_kg_per_s,original.branches[1].solution.mass_flow_kg_per_s);
    });
}

#[test]
fn whole_network_budgets_do_not_reset_for_each_branch() {
    with_network(false,false,1.0,|cx,network| {
        assert!(network.solve_pressure_driven(cx,&[400.0,400.0],Pressure::new(50.0),local(),
            ParallelPressureBudget { max_flow_evaluations: 2,..budget() }).is_err());
        // Each branch dissipates 0.000125 W, below 0.0002 W. Their sum is not.
        assert!(matches!(network.solve_pressure_driven(cx,&[300.0,300.0],Pressure::new(50.0),local(),
            ParallelPressureBudget { maximum_neglected_hydraulic_power_w: 0.0002,..budget() }),
            Err(ParallelPressureError::NeglectedDissipation { .. })));
    });
}

#[test]
fn duplicate_regions_and_malformed_late_branch_refuse() {
    with_network(false,false,1.0,|cx,network| {
        let first = network.branches()[0];
        assert!(ParallelLiquidNetwork::try_new(cx,&[first,ParallelLiquidBranch { name: "other",..first }]).is_err());
        let mut bad = network.branches().to_vec();
        bad[1].pressure_budget.mass_flow_bracket_kg_per_s = (0.001,0.0002);
        let bad = ParallelLiquidNetwork::try_new(cx,&bad).unwrap();
        assert!(bad.solve_pressure_driven(cx,&[300.0,300.0],Pressure::new(50.0),local(),budget()).is_err());
        assert!(network.solve_pressure_driven(cx,&[300.0,f64::NAN],Pressure::new(50.0),local(),budget()).is_err());
        assert!(network.solve_pressure_driven(cx,&[300.0],Pressure::new(50.0),local(),budget()).is_err());
    });
}

#[test]
fn mixing_different_material_charts_is_not_silently_permitted() {
    with_source(false,false,|cx,law,source| {
        with_source(true,false,|_,other_law,other_source| {
            let a = path(cx,law,source,"a-wall");
            let b = path(cx,other_law,other_source,"b-wall");
            let ar = loss("a-wall",1.0);
            let br = loss("b-wall",1.0);
            assert!(ParallelLiquidNetwork::try_new(cx,&[
                ParallelLiquidBranch { name: "a",path: &a,resistances: &ar,pressure_budget: branch_budget() },
                ParallelLiquidBranch { name: "b",path: &b,resistances: &br,pressure_budget: branch_budget() },
            ]).is_err());
        });
    });
}
