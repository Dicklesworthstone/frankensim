//! Actual radiative FEM / material coolant / pump / storage integration.
use fs_airflow::conjugate::{ConjugateConfig, Relaxation};
use fs_airflow::graph::thermal::medium::FlowHeatGeometry;
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{
    LiquidResistanceSection, PressureDrivenBudget};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::{
    ParallelLiquidBranch, ParallelLiquidNetwork, ParallelPressureBudget};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::solid::ParallelLiquidCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpHeadCurve, PumpHeadPoint, PumpOperatingControls};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::solid::PumpCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::recirculation::ReservoirStepControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::{
    StoredCoolingState, StoredCoolingControls, StoredCoolingError, StoredSolidProblem};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::radiation::RadiatingStoredCoolingControls;
use fs_airflow::SourceProvenance;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, SurfaceEmissivity,
    ConductionError, ConductionMesh, ConductivityModel, ScalarField, ThermalBc,
    ThermalBoundary, ThermalBoundaryBuilder, SURFACE_EMISSIVITY_PROPERTY};
use fs_conduction::transient::enthalpy::{EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig};
use fs_convection::{CorrelationId, CorrelationRobinBoundary};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, QuantitySpec, Temperature, VolumetricFlowRate};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../examples/pressure_cooling/source.rs"]
mod source;

fn chart(table: &source::ThermalTable, phase: SolidLiquidPhase) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::with_quantity("specific_heat_capacity", QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0])));
    let rho = PropertyKey::with_quantity("density", QuantitySpec::dimensional(Dims([-3,1,0,0,0,0])));
    ResolvedHeatCapacityCurve::try_new(&table.card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T",
        lower: &table.point(table.lower_k).unwrap(), upper: &table.point(table.upper_k).unwrap(),
        selection: SelectionPolicy::SingleClaimOnly,
    }, phase, 0.0, 0.01).unwrap()
}
fn patch(epsilon: f64, ambient: f64, span: (f64,f64)) -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::with_quantity(SURFACE_EMISSIVITY_PROPERTY, QuantitySpec::dimensional(Dims::NONE)),
        value: PropertyValue::Scalar { value: epsilon, dims: Dims::NONE },
        validity: ValidityDomain::unconstrained().with("T", span.0, span.1),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic finish, not measurements".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic surface".into(),
        phase: "solid".into(), process: "test coating".into(), revision: 0 }, claims, vec![]).unwrap();
    AmbientRadiationPatch::new("exposed", SurfaceEmissivity::from_card("exposed", &card,
        300.0, SelectionPolicy::SingleClaimOnly).unwrap(), ambient).unwrap()
}
fn controls() -> RadiatingStoredCoolingControls {
    RadiatingStoredCoolingControls {
        thermal: StoredCoolingControls {
            pump: PumpCouplingControls {
                network: ParallelLiquidCouplingControls {
                    coupling: ConjugateConfig { temperature_tolerance_k: 1e-8, max_iterations: 48,
                        balance_tolerance_w: 1e-7, balance_relative_tolerance: 1e-8,
                        relaxation: Relaxation::Fixed { omega: 1.0 } },
                    local: LiquidFeedbackBudget { temperature_tolerance_k: 1e-10,
                        coefficient_relative_tolerance: 1e-9, max_region_iterations: 64 },
                    acceleration: Default::default(), mass_flow_relative_tolerance: 1e-8,
                    work: ParallelPressureBudget { max_flow_evaluations: 1_000_000,
                        maximum_neglected_hydraulic_power_w: 0.01 },
                },
                operating_point: PumpOperatingControls { pressure_bracket_pa: (20.0,80.0),
                    absolute_pressure_tolerance_pa: 1e-8, relative_pressure_tolerance: 1e-9,
                    max_network_evaluations: 48 }, max_total_solid_solves: 8192,
            },
            reservoir: ReservoirStepControls { temperature_bracket_k: (299.0,302.0),
                energy_absolute_tolerance_j: 1e-7, max_inlet_evaluations: 48,
                maximum_neglected_hydraulic_energy_j: 0.1 },
            solid: EnthalpyStepConfig { newton: fs_solver::NewtonKrylovConfig {
                absolute_tolerance: 1e-10, relative_tolerance: 1e-11,
                linear_restart: 24, max_linear_cycles: 16, forcing_minimum: 1e-12,
                forcing_maximum: 1e-3, ..Default::default()
            }, max_newton_iterations: 32, energy_tolerance_j: 1e-8 },
            joint_energy_tolerance_j: 2e-6,
        },
        radiation: AmbientRadiationConfig { max_iterations: 48, relaxation: 1.0,
            temperature_tolerance_k: 1e-9, balance_tolerance_w: 1e-8,
            balance_relative_tolerance: 1e-10 }, max_total_enthalpy_solves: 65536,
    }
}
fn boundary(mesh: &ConductionMesh, rows: &[CorrelationRobinBoundary], air: f64)
    -> std::result::Result<ThermalBoundary, ConductionError> {
    // Deliberately declare the external trace FIRST: it must not become coolant.
    ThermalBoundaryBuilder::new(mesh)
        .region("exposed", |f| f.outward_normal[2] < -0.99, ThermalBc::robin(5.0, air)?)?
        .region("body", |f| f.outward_normal[2] > 0.99, rows[0].boundary_condition().clone())?
        .adiabatic_remainder().finish()
}
// Preserve the shared fluid chart lifetime across the callback's mutable state.
fn fixture<R>(run: impl for<'fluid> FnOnce(&Cx<'_>, &CancelGate, &ParallelLiquidNetwork<'_, 'fluid, '_>,
    &PumpHeadCurve, &mut StoredCoolingState<'fluid, '_, '_, '_>, &ConductionMesh, &ConductivityModel) -> R) -> R {
    let liquid = source::ThermalTable::parse(include_str!("../examples/pressure_cooling/constant.txt")).unwrap();
    let solid = source::ThermalTable::parse_solid(include_str!("../examples/recirculating_cooling/solid_low_capacity.txt")).unwrap();
    let law = LiquidTemperatureLaw::try_new(&liquid.card, &liquid.point(300.0).unwrap(), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let fluid_chart = chart(&liquid, SolidLiquidPhase::Liquid);
    let solid_chart = chart(&solid, SolidLiquidPhase::Solid);
    let (complex, positions) = fs_conduction::fixtures::box_grid([2,1,1], [0.06,0.1,0.01]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let material = ConductivityModel::isotropic_declared(5.0).unwrap();
    let pump = PumpHeadCurve::try_new(vec![
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.0), head: Length::new(0.005) },
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.001), head: Length::new(0.005) },
    ], 10.0, (1e-6,3e-6), SourceProvenance::new("synthetic head", "radiation-test")).unwrap();
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 95, kernel_id: 747, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let geometry = FlowHeatGeometry { flow_area: Area::new(core::f64::consts::PI*0.01*0.01/4.0),
            characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
            wetted_area: Area::new(0.006), aspect_ratio: None };
        let path = LiquidHeatPath::try_new(&cx, &law, &fluid_chart, MassFlowRate::new(0.001),
            Temperature::new(300.0), &[LiquidHeatSection { region: "body", geometry,
                correlation: CorrelationId::CircularDuctLaminarCwt }], LiquidPropertyTemperature::BulkMean).unwrap();
        let resistance = [LiquidResistanceSection { region: "body", resistance: MaterialFlowResistance::try_new(
            1e10, 0.0, geometry.flow_area, geometry.characteristic_length, (0.0,2300.0)).unwrap() }];
        let network = ParallelLiquidNetwork::try_new(&cx, &[ParallelLiquidBranch { name: "loop", path: &path,
            resistances: &resistance, pressure_budget: PressureDrivenBudget {
                mass_flow_bracket_kg_per_s: (0.0002,0.008), pressure_absolute_tolerance_pa: 1e-10,
                pressure_relative_tolerance: 1e-10, max_flow_evaluations: 64,
                maximum_neglected_hydraulic_power_w: 0.01 } }]).unwrap();
        let storage = EnthalpyBackwardEuler::uniform(&cx, &mesh, solid_chart.curve(), 1000.0,
            EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() }).unwrap();
        let h = solid_chart.state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg();
        let tank = network.reservoir(&cx, 0.1, Temperature::new(300.0), 0.0).unwrap();
        let mut state = StoredCoolingState::try_new(&cx, tank, &storage, vec![h; mesh.vertex_count()]).unwrap();
        run(&cx, &gate, &network, &pump, &mut state, &mesh, &material)
    })
}

#[test]
fn radiation_and_external_convection_do_not_enter_the_coolant_heat_account() {
    fixture(|cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("body", 20.0/0.00006).unwrap();
        let patches = [patch(0.9,250.0,(250.0,450.0))];
        let result = network.advance_stored_solid_with_pump_and_radiation(cx, state, 1.0, 0.0,
            &[300.0], pump, StoredSolidProblem { mesh, material, source: &heat,
                element_materials: None, interfaces: None }, &patches, controls(),
            |_, rows| boundary(mesh, rows, 310.0)).unwrap();
        let response = &result.step.liquid.circuit.state.solid.state;
        let solid = &response.thermal.conduction;
        let radiation = response.thermal.radiation.nonlinear_radiation_out_w;
        assert!(radiation > 1.0);
        assert!(response.external_convection_out_w < 0.0); // separate HOT air reservoir
        assert_eq!(response.coolant_robin_fluxes.len(),1);
        assert_eq!(response.coolant_robin_fluxes[0].region,"body");
        assert_eq!(response.thermal.convective_robin_fluxes.len(),2);
        assert!((response.coolant_out_w-result.step.liquid.returned_heat_j).abs()<2e-6);
        assert!((solid.robin_out_w-result.step.liquid.returned_heat_j).abs()>0.5);
        assert!((result.step.solid_external_heat_j-(20.0-radiation-response.external_convection_out_w)).abs()<1e-8);
        assert!(result.step.combined_energy_residual_j.abs()<2e-6);
        assert!(result.step.interface_energy_residual_j.abs()<2e-6);
        assert!(response.thermal.physical_residual_norm_j<=response.thermal.physical_residual_tolerance_j);
        assert!(response.thermal.physical_energy_residual_j.abs()<=controls().thermal.solid.energy_tolerance_j);
        assert_eq!(response.thermal.radiation.patches[0].patch,patches[0]);
        assert_eq!(state.solid_specific_enthalpy_j_kg(),solid.specific_enthalpy_j_kg.as_slice());
        assert!(result.total_enthalpy_solves>=result.step.liquid.total_solid_solves);
        assert!(result.total_enthalpy_solves<=controls().max_total_enthalpy_solves);
        assert!(result.total_krylov_iterations>0);
    });
}

#[test]
fn changing_only_surface_material_changes_solid_and_coolant_history() {
    let run = |epsilon| fixture(|cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("body",20.0/0.00006).unwrap();
        let patches = [patch(epsilon,250.0,(250.0,450.0))];
        let result = network.advance_stored_solid_with_pump_and_radiation(cx,state,1.0,0.0,
            &[300.0],pump,StoredSolidProblem {mesh,material,source:&heat,element_materials:None,interfaces:None},
            &patches,controls(),|_,rows|boundary(mesh,rows,300.0)).unwrap();
        (state.reservoir().state().temperature_k(),result.step.liquid.circuit.state.solid.state.thermal.conduction.temperature)
    });
    let low = run(0.1); let high = run(0.9);
    assert!(high.0<low.0);
    assert!(high.1.iter().zip(&low.1).all(|(a,b)|a<b));
}

#[test]
fn a_hot_radiative_reservoir_heats_both_histories_without_a_body_source() {
    fixture(|cx,_,network,pump,state,mesh,material| {
        let heat=ScalarField::uniform("no body heat",0.0).unwrap();
        let result=network.advance_stored_solid_with_pump_and_radiation(cx,state,1.0,0.0,&[300.0],pump,
            StoredSolidProblem {mesh,material,source:&heat,element_materials:None,interfaces:None},
            &[patch(0.9,400.0,(250.0,450.0))],controls(),|_,rows|boundary(mesh,rows,300.0)).unwrap();
        let response=&result.step.liquid.circuit.state.solid.state;
        assert!(response.thermal.radiation.nonlinear_radiation_out_w<0.0);
        assert!(state.reservoir().state().temperature_k()>300.0);
        assert!(result.step.solid_external_heat_j>0.0);
        assert!(response.thermal.conduction.stored_energy_change_j>0.0);
        assert!(result.step.combined_energy_residual_j.abs()<2e-6);
    });
}

#[test]
fn source_failure_and_inner_work_exhaustion_leave_both_histories_unchanged() {
    for (span, maximum) in [((250.0,450.0),1), ((300.0,300.0),65536)] {
        fixture(|cx,_,network,pump,state,mesh,material| {
            let old=state.solid_specific_enthalpy_j_kg().to_vec();
            let tank=state.reservoir().clone();
            let heat=ScalarField::uniform("body",20.0/0.00006).unwrap();
            let result=network.advance_stored_solid_with_pump_and_radiation(cx,state,1.0,0.0,&[300.0],pump,
                StoredSolidProblem {mesh,material,source:&heat,element_materials:None,interfaces:None},
                &[patch(0.9,250.0,span)], RadiatingStoredCoolingControls {max_total_enthalpy_solves:maximum,..controls()},
                |_,rows|boundary(mesh,rows,300.0));
            assert!(result.is_err());
            assert_eq!(state.reservoir(),&tank);
            assert_eq!(state.solid_specific_enthalpy_j_kg(),old.as_slice());
        });
    }
}

#[test]
fn undeclared_heat_ownership_and_cancellation_refuse_without_publication() {
    for cancel in [false,true] {
        fixture(|cx,gate,network,pump,state,mesh,material| {
            let old=state.solid_specific_enthalpy_j_kg().to_vec(); let tank=state.reservoir().clone();
            let heat=ScalarField::uniform("body",20.0/0.00006).unwrap();
            let result=network.advance_stored_solid_with_pump_and_radiation(cx,state,1.0,0.0,&[300.0],pump,
                StoredSolidProblem {mesh,material,source:&heat,element_materials:None,interfaces:None},
                &[patch(0.9,250.0,(250.0,450.0))],controls(),|_,rows| {
                    if cancel {gate.request();return boundary(mesh,rows,300.0);}
                    ThermalBoundaryBuilder::new(mesh)
                        .region("body",|f|f.outward_normal[2]>0.99,rows[0].boundary_condition().clone())?
                        .region("undeclared",|f|f.outward_normal[2] < -0.99,ThermalBc::robin(5.0,300.0)?)?
                        .adiabatic_remainder().finish()
                });
            assert!(result.is_err());
            if !cancel {assert!(matches!(result,Err(StoredCoolingError::Invalid(_))));}
            assert_eq!(state.reservoir(),&tank);assert_eq!(state.solid_specific_enthalpy_j_kg(),old.as_slice());
        });
    }
}
