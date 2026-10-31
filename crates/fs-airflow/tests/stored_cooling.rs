//! Actual spatial enthalpy / pump / reservoir integration, not a surrogate callback.
use fs_airflow::conjugate::{ConjugateConfig, Relaxation};
use fs_airflow::graph::thermal::medium::FlowHeatGeometry;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::{
    LiquidFeedbackBudget, LiquidHeatPath, LiquidHeatSection, LiquidPropertyTemperature,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::{
    LiquidResistanceSection, PressureDrivenBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::{
    ParallelLiquidBranch, ParallelLiquidNetwork, ParallelPressureBudget,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::solid::ParallelLiquidCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpHeadCurve, PumpHeadPoint, PumpOperatingControls,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::solid::PumpCouplingControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::recirculation::ReservoirStepControls;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::{
    StoredCoolingState, StoredCoolingControls, StoredCoolingError, StoredSolidProblem,
};
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::SourceProvenance;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductivityModel, ConductivityTable,
    ScalarField, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::transient::enthalpy::{EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig};
use fs_convection::{CorrelationId, CorrelationRobinBoundary};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Area, Dims, Length, MassFlowRate, QuantitySpec, Temperature, VolumetricFlowRate};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../examples/pressure_cooling/source.rs"]
mod source;
#[path = "stored_cooling/heterogeneous.rs"]
mod heterogeneous;
#[path = "stored_cooling/hydraulic_heat.rs"]
mod hydraulic_heat;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("pressure", 100000.0).unwrap()
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    PropertyKey::with_quantity(name, if name == "thermal_conductivity" {
        QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
    } else { QuantitySpec::dimensional(dims) })
}
fn caloric(card: &MaterialCard, phase: SolidLiquidPhase) -> ResolvedHeatCapacityCurve {
    let cp = key("specific_heat_capacity", Dims([2,0,-2,-1,0,0]));
    let rho = key("density", Dims([-3,1,0,0,0,0]));
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T",
        lower: &point(250.0), upper: &point(450.0), selection: SelectionPolicy::SingleClaimOnly,
    }, phase, 0.0, 0.01).unwrap()
}
fn solid_card(cp_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3,1,0,0,0,0]), [1000.0, 1000.0]),
        ("specific_heat_capacity", Dims([2,0,-2,-1,0,0]), [400.0*cp_scale, 800.0*cp_scale]),
        ("thermal_conductivity", Dims([1,1,-3,-1,0,0]), [4.0, 8.0]),
    ] {
        claims.insert_claim(PropertyClaim {
            key: key(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, values[0]), (450.0, values[1])], dims },
            validity: ValidityDomain::unconstrained().with("T",250.0,450.0).with("pressure",100000.0,100000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic solid storage regression, not measurements".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic solid".into(), phase: "solid".into(),
        process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn controls() -> StoredCoolingControls {
    StoredCoolingControls {
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
                max_network_evaluations: 48 },
            max_total_solid_solves: 8192,
        },
        reservoir: ReservoirStepControls { temperature_bracket_k: (299.0,302.0),
            energy_absolute_tolerance_j: 1e-7, max_inlet_evaluations: 48,
            maximum_neglected_hydraulic_energy_j: 0.1 },
        solid: EnthalpyStepConfig { newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10, relative_tolerance: 1e-11,
            linear_restart: 24, max_linear_cycles: 16, forcing_minimum: 1e-12,
            forcing_maximum: 1e-3, ..Default::default()
        }, max_newton_iterations: 32, energy_tolerance_j: 1e-8 },
        joint_energy_tolerance_j: 2e-6,
    }
}
fn boundary(mesh: &ConductionMesh, rows: &[CorrelationRobinBoundary]) -> std::result::Result<ThermalBoundary, ConductionError> {
    ThermalBoundaryBuilder::new(mesh).region("body", |f| f.outward_normal[2] > 0.99,
        rows[0].boundary_condition().clone())?.adiabatic_remainder().finish()
}
// The network and mutable reservoir retain the SAME fluid chart lifetime.
// Independent elided lifetimes lose that relation at the callback boundary.
fn fixture<R>(cp_scale: f64, run: impl for<'fluid> FnOnce(&Cx<'_>, &CancelGate,
    &ParallelLiquidNetwork<'_, 'fluid, '_>, &PumpHeadCurve,
    &mut StoredCoolingState<'fluid, '_, '_, '_>, &ConductionMesh, &ConductivityModel) -> R) -> R {
    let fluid = source::LiquidTable::parse(include_str!("../examples/pressure_cooling/constant.txt")).unwrap();
    let law = LiquidTemperatureLaw::try_new(&fluid.card, &point(300.0), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let fluid_chart = caloric(&fluid.card, SolidLiquidPhase::Liquid);
    let body = solid_card(cp_scale);
    let solid_chart = caloric(&body, SolidLiquidPhase::Solid);
    let conductivity = ConductivityModel::isotropic(ConductivityTable::from_claims_at_query_points(
        body.claims(), &key("thermal_conductivity", Dims([1,1,-3,-1,0,0])), "T",
        &[point(250.0), point(450.0)], SelectionPolicy::SingleClaimOnly).unwrap());
    let (complex, positions) = fs_conduction::fixtures::box_grid([2,1,1], [0.06,0.1,0.01]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let pump = PumpHeadCurve::try_new(vec![
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.0), head: Length::new(0.005) },
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.001), head: Length::new(0.005) },
    ], 10.0, (1e-6, 3e-6), SourceProvenance::new("synthetic constant head", "stored-cooling-test")).unwrap();
    let cancel = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&cancel, arena, StreamKey { seed: 94, kernel_id: 746, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let geometry = FlowHeatGeometry { flow_area: Area::new(core::f64::consts::PI*0.01*0.01/4.0),
            characteristic_length: Length::new(0.01), heated_length: Length::new(2.0),
            wetted_area: Area::new(0.006), aspect_ratio: None };
        let path = LiquidHeatPath::try_new(&cx, &law, &fluid_chart, MassFlowRate::new(0.001),
            Temperature::new(300.0), &[LiquidHeatSection { region: "body", geometry,
                correlation: CorrelationId::CircularDuctLaminarCwt }], LiquidPropertyTemperature::BulkMean).unwrap();
        let resistance = [LiquidResistanceSection { region: "body", resistance: MaterialFlowResistance::try_new(
            1e10, 0.0, geometry.flow_area, geometry.characteristic_length, (0.0,2300.0)).unwrap() }];
        let network = ParallelLiquidNetwork::try_new(&cx, &[ParallelLiquidBranch {
            name: "loop", path: &path, resistances: &resistance,
            pressure_budget: PressureDrivenBudget { mass_flow_bracket_kg_per_s: (0.0002,0.008),
                pressure_absolute_tolerance_pa: 1e-10, pressure_relative_tolerance: 1e-10,
                max_flow_evaluations: 64, maximum_neglected_hydraulic_power_w: 0.01 },
        }]).unwrap();
        let storage = EnthalpyBackwardEuler::uniform(&cx, &mesh, solid_chart.curve(), 1000.0,
            EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() }).unwrap();
        let initial_h = solid_chart.state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg();
        let tank = network.reservoir(&cx, 0.1, Temperature::new(300.0), 0.0).unwrap();
        let mut state = StoredCoolingState::try_new(&cx, tank, &storage, vec![initial_h; mesh.vertex_count()]).unwrap();
        run(&cx, &cancel, &network, &pump, &mut state, &mesh, &conductivity)
    })
}

#[test]
fn actual_solid_storage_delays_tank_heating_and_both_histories_are_published() {
    fixture(1.0, |cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("20 W body source", 20.0/0.00006).unwrap();
        let old = state.solid_specific_enthalpy_j_kg().to_vec();
        let result = network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
            controls(), |_, rows| boundary(mesh, rows)).unwrap();
        let solid = &result.liquid.circuit.state.solid.state;
        assert!(solid.stored_energy_change_j > 15.0);
        assert!(result.liquid.returned_heat_j > 0.0 && result.liquid.returned_heat_j < 5.0);
        assert!(state.reservoir().state().temperature_k() < 300.025); // steady body would give 300.1 K
        assert!((solid.stored_energy_change_j+result.liquid.stored_enthalpy_change_j-20.0).abs() < 2e-6);
        assert_eq!(state.reservoir().time_s(), 1.0);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), solid.specific_enthalpy_j_kg.as_slice());
        assert_ne!(old, state.solid_specific_enthalpy_j_kg());
        assert!(result.combined_energy_residual_j.abs() <= controls().joint_energy_tolerance_j);
        assert!(result.interface_energy_residual_j.abs() <= controls().joint_energy_tolerance_j);
    });
}

#[test]
fn pulse_off_preserves_stored_heat_instead_of_resetting_the_solid() {
    fixture(1.0, |cx, _, network, pump, state, mesh, material| {
        let initial_h = state.solid_specific_enthalpy_j_kg().to_vec();
        let initial_tank_h = state.reservoir().state().specific_enthalpy_j_kg();
        let mut walls = vec![300.0];
        let mut tank_after_pulse = 0.0;
        for (end, watts) in [(1.0,20.0),(2.0,0.0)] {
            let heat = ScalarField::uniform("pulse", watts/0.00006).unwrap();
            let result = network.advance_stored_solid_with_pump(cx, state, end, 0.0, &walls, pump,
                StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
                controls(), |_, rows| boundary(mesh, rows)).unwrap();
            walls = result.liquid.circuit.state.solid.regions.iter().map(|r| r.mean_wall_temperature_k).collect();
            if end == 1.0 { tank_after_pulse = state.reservoir().state().temperature_k(); }
            else {
                assert!(state.reservoir().state().temperature_k() > tank_after_pulse);
                assert!(result.liquid.circuit.state.solid.state.stored_energy_change_j < 0.0);
            }
        }
        let solid_j: f64 = state.storage().reference_nodal_masses_kg().iter()
            .zip(state.solid_specific_enthalpy_j_kg()).zip(&initial_h).map(|((&m,&h),&old)| m*(h-old)).sum();
        let tank_j = state.reservoir().mass_kg()*(state.reservoir().state().specific_enthalpy_j_kg()-initial_tank_h);
        assert!((solid_j+tank_j-20.0).abs() < 4e-6);
    });
}

#[test]
fn changing_only_solid_caloric_data_changes_the_coupled_response() {
    let response = |scale| fixture(scale, |cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("20 W", 20.0/0.00006).unwrap();
        let result = network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
            controls(), |_, rows| boundary(mesh, rows)).unwrap();
        (state.reservoir().state().temperature_k(), result.liquid.circuit.state.solid.state.temperature)
    });
    let a = response(1.0);
    let b = response(2.0);
    assert!(b.0 < a.0);
    assert!(b.1.iter().sum::<f64>() < a.1.iter().sum::<f64>());
}

#[test]
fn original_boundary_failure_and_cancel_leave_both_histories_unchanged() {
    fixture(1.0, |cx, cancel, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("20 W", 20.0/0.00006).unwrap();
        let problem = StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None };
        let old_tank = state.reservoir().clone();
        let old_h = state.solid_specific_enthalpy_j_kg().to_vec();
        let failed = network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump, problem,
            controls(), |_, _| Err(ConductionError::Mesh { what: "test boundary refusal".into(), fix: "none".into() }));
        assert!(matches!(failed, Err(StoredCoolingError::Boundary(ConductionError::Mesh { .. }))));
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_h.as_slice());
        let failed = network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump, problem,
            controls(), |_, rows| { cancel.request(); boundary(mesh, rows) });
        assert!(failed.is_err());
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_h.as_slice());
    });
}

#[test]
fn prescribed_bath_is_refused_before_advancing_either_storage_history() {
    fixture(1.0, |cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("20 W", 20.0/0.00006).unwrap();
        let old_tank = state.reservoir().clone();
        let old_h = state.solid_specific_enthalpy_j_kg().to_vec();
        let result = network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
            controls(), |_, rows| {
                ThermalBoundaryBuilder::new(mesh)
                    .region("body", |f| f.outward_normal[2] > 0.99,
                        rows[0].boundary_condition().clone())?
                    .region("bath", |f| f.outward_normal[2] < -0.99,
                        fs_conduction::bc::ThermalBc::dirichlet(320.0)?)?
                    .adiabatic_remainder().finish()
            });
        assert!(matches!(result, Err(StoredCoolingError::Solid(
            fs_conduction::transient::enthalpy::EnthalpyError::UnsupportedDirichlet))));
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_h.as_slice());
    });
}

#[test]
fn budget_or_invalid_joint_control_cannot_publish_either_storage_state() {
    fixture(1.0, |cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("20 W", 20.0/0.00006).unwrap();
        let old_tank = state.reservoir().clone();
        let old_h = state.solid_specific_enthalpy_j_kg().to_vec();
        let mut limits = controls();
        limits.pump.max_total_solid_solves = 1;
        assert!(network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
            limits, |_, rows| boundary(mesh, rows)).is_err());
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_h.as_slice());
        let called = std::cell::Cell::new(false);
        let mut invalid = controls(); invalid.joint_energy_tolerance_j = 0.0;
        assert!(network.advance_stored_solid_with_pump(cx, state, 1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
            invalid, |_, rows| { called.set(true); boundary(mesh, rows) }).is_err());
        assert!(!called.get());
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_h.as_slice());
    });
}
