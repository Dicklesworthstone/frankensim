//! Manufactured sources, independent thermoelastic U and physical nodal work.
use super::*;
use crate::air_path::cylinder::{CylinderFlowResolution, frame::{FrameFluidBinding, FrameFluidBudget}};
use crate::frame_vibration::FrameVibrationBudget;
use crate::circuit_enthalpy::material::{MaterialCell, CellResistor};
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::circuit_enthalpy::ResistorCell;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::phase::{ResolvedHeatCapacityCurve, HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point, SPECIFIC_HEAT_CAPACITY_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS};
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, Geometry};
use fs_thermal::enthalpy::*;
use fs_phs::circuit::{Branch, CircuitGraph};
use crate::material_circuit::MaterialCircuitSolve;

const CELL: CellId = CellId(17);
const RHO: f64 = 1000.0;
const RADIUS: f64 = 0.02;
const YOUNG: f64 = 2e6;
const CP: f64 = 0.2;
fn area() -> f64 { core::f64::consts::PI*RADIUS*RADIUS }
fn moment() -> f64 { area()*RADIUS*RADIUS/4.0 }
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
        .with("omega", 1.0).unwrap()
}
fn sources(slope: f64, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, RHO), ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, CP),
        ("young_modulus", Pressure::DIMS, YOUNG),
        ("axial_eta", DynViscosity::DIMS, 1000.0), ("bending_eta", DynViscosity::DIMS, 1000.0),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, alpha),
        ("electrical_resistivity", fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS, 0.02),
    ] {
        let varying = name == "young_modulus" && slope != 0.0;
        let value = if varying { PropertyValue::Curve { abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, value-100.0*slope),(600.0,value+300.0*slope)], dims }
        } else { PropertyValue::Scalar { value, dims } };
        let mut validity = ValidityDomain::unconstrained().with("T",200.0,600.0);
        if name.ends_with("_eta") { validity = validity.with("omega",0.0,1e7); }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims), value, validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if varying { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "synthetic wet thermoelastic fixture".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "authored solid".into(),
        phase: "solid".into(), process: "analytic reference only".into(), revision: 0 }, claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density",Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(600.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, curve)
}
fn frame(card: &MaterialCard) -> Frame2d {
    let requirements = [("density",Density::DIMS),("young_modulus",Pressure::DIMS),
        ("axial_eta",DynViscosity::DIMS),("bending_eta",DynViscosity::DIMS),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS)]
        .into_iter().map(|(n,d)| ScalarPropertyRequirement::try_new(n,d,ScalarAdmissibility::Finite).unwrap()).collect::<Vec<_>>();
    let material = resolve_material_state_point(card,&point(300.0),&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let section = ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {
        outer_radius_m: RADIUS, material }]).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0,0.0); frame.add_node(1.0,0.0);
    frame.fix(0,[true;3]); frame.fix(1,[false,false,true]);
    frame.add_material_viscous_element(0,1,3,
        FrameViscousSection::try_new(section,"axial_eta","bending_eta").unwrap()).unwrap();
    frame
}
fn thermal(curve: &ResolvedHeatCapacityCurve, power: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(vec![TimedEvent { time_s: 0.0,
        event: Event::Birth(DepositedCell { id: CELL, segment_id: 0, body_id: 0,
            material_id: MaterialId(0), mass_kg: RHO*area(), specific_h_j_kg: 100.0*CP,
            geometry: CellGeometry { start_m: [0.0;3], end_m: [1.0,0.0,0.0], volume_m3: area(),
                surface_area_m2: 2.0*core::f64::consts::PI*RADIUS },
            boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 0.0,reservoir: None,source_w: power },
        }) }], &[(MaterialId(0),curve.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2,max_links: 1,max_events: 2,max_substeps: 1000 }).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-10},||Ok(())).unwrap().0
}
fn model(frame: &Frame2d, wet: bool) -> MaterialFluidFrame<'_,FluidTransportState> {
    MaterialFluidFrame::new(frame,2, if wet { vec![FrameFluidBinding { element:0,
        fluid: FluidTransportState::try_new(400.0,4000.0,0.5,0.02).unwrap() }] } else { vec![] },
        CylinderFlowResolution {minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:2.0,intervals:2},
        FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36},
        FrameFluidBudget {max_states:12,max_dense_entries:144,max_projection_terms:4096,mass_inverse_tolerance:1e-10},
        &CancelGate::new_clock_free()).unwrap()
}
fn config() -> FrameFluidThermoelasticConfig {
    FrameFluidThermoelasticConfig { sample_rate_hz:1000,modal:ModalAcousticTimeBudget::audible_reference(),
        coupling: ModalCouplingConfig { max_modes:16,max_connections:32,max_setup_terms:32768,
            nyquist_guard_fraction:0.9,maximum_total_energy_j:10.0,maximum_abs_pressure_pa:100.0,
            maximum_abs_connection_force_n:1e6,solve_relative_tolerance:1e-10,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-9 },
        prepared:None,maximum_projection_terms:4096,maximum_absolute_strain:0.05,initial_energy_tolerance_j:1e-9 }
}
fn controls() -> FrameThermalSolve {
    FrameThermalSolve { maximum_cells:2,maximum_material_trials:200,maximum_newton_iterations:16,
        difference_step_j:1e-6,cell_energy_tolerance_j:1e-9,thermal_energy_tolerance_j:1e-9,
        total_energy_tolerance_j:2e-8,newton:NewtonKrylovConfig {absolute_tolerance:1e-11,
            relative_tolerance:1e-12,..NewtonKrylovConfig::default()} }
}
fn near(a: f64,b: f64,tol: f64) { assert!((a-b).abs()<=tol,"{a:e} != {b:e}"); }
fn fixture<R>(slope:f64,alpha:f64,power:f64,wet:bool,
    run:impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,DepositionSimulation<'_>,&MaterialNetwork<'_>)->R) -> R {
    let (card,curve)=sources(slope,alpha);
    let frame=frame(&card); let model=model(&frame,wet); let thermal=thermal(&curve,power);
    let context=point(300.0);
    let cells=[MaterialCell {cell:CELL,material_id:MaterialId(0),curve:curve.curve(),
        point:ThermalMaterialPoint {material:&card,context:&context,temperature_axis:"T",liquid_fraction_axis:None},conductivity:None}];
    let source=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    run(&model,thermal,&source)
}
fn initial(model:&MaterialFluidFrame<'_,FluidTransportState>) -> Vec<f64> {
    model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:1e-3,
        velocity_m_sqrt_kg_per_s:0.02},ModalAcousticState::default()],
        &vec![1e-4;model.memory_masses_kg().len()]).unwrap()
}
const DESTINATIONS:[FrameLayerCell;1]=[FrameLayerCell {element:0,layer:0,cell:CELL}];

#[test]
fn binding_keeps_original_motion_energy_and_every_fluid_momentum() {
    fixture(-2000.0,1e-4,10.0,true,|model,thermal,source| {
        let state=initial(model); let gate=CancelGate::new_clock_free();
        let owner=FrameFluidThermoelasticSystem::new(model,&state,&thermal,source,&DESTINATIONS,config(),&gate).unwrap();
        assert!(core::ptr::eq(owner.source(),model));
        for (a,b) in state.iter().zip(owner.physical_state().unwrap()) {near(*a,b,2e-12);}
        near(owner.mechanical().total_energy_j().unwrap(),model.system().hamiltonian(&state),1e-10);
        assert_eq!(owner.mechanical().components()[1].states().len(),model.memory_masses_kg().len());
        let mut too_small=config();too_small.coupling.max_modes=9;
        assert!(FrameFluidThermoelasticSystem::new(model,&state,&thermal,source,&DESTINATIONS,too_small,&gate).is_err());
    });
}

fn independent_total(owner:&FrameFluidThermoelasticSystem<'_, '_, '_, '_, FluidTransportState>,
    thermal:&DepositionSimulation<'_>,slope:f64,alpha:f64)->f64 {
    let state=owner.physical_state().unwrap();let motion=owner.modal_states().unwrap();
    let (mut u,mut v)=([0.0;6],[0.0;6]);
    owner.reconstruct_into(&mut u,&mut v,&CancelGate::new_clock_free()).unwrap();
    let m=motion.len();
    let kinetic=0.5*(0..m).map(|i| state[m+i]*motion[i].velocity_m_sqrt_kg_per_s).sum::<f64>()
        +0.5*state[2*m..].iter().zip(owner.source().memory_masses_kg()).map(|(p,m)| p*p/m).sum::<f64>();
    let t=thermal.temperature_k(CELL).unwrap();let e=YOUNG+slope*(t-300.0);
    let elastic=u[3]-alpha*(t-300.0);
    // U=Psi-T*dPsi/dT. This is NOT the native spring Hamiltonian.
    let axial=area()*(0.5*(e-t*slope)*elastic*elastic+t*e*alpha*elastic);
    let bending=0.5*(e-t*slope)*12.0*moment()*u[4]*u[4];
    kinetic+axial+bending+RHO*area()*CP*(t-200.0)
}

#[test]
fn heated_affine_modulus_and_expansion_close_independent_total_energy_in_fluid() {
    let mut endpoints=Vec::new();
    for (slope,alpha) in [(0.0,0.0),(-2000.0,1e-4)] {
        fixture(slope,alpha,10.0,true,|model,mut thermal,source| {
            let gate=CancelGate::new_clock_free();let state=initial(model);
            let mut owner=FrameFluidThermoelasticSystem::new(model,&state,&thermal,source,&DESTINATIONS,config(),&gate).unwrap();
            let initial=independent_total(&owner,&thermal,slope,alpha);
            let (mut work,mut heat,mut outside)=(0.0,0.0,0.0);
            for _ in 0..24 {
                let (mut before,mut velocity)=([0.0;6],[0.0;6]);
                owner.reconstruct_into(&mut before,&mut velocity,&gate).unwrap();
                let force=[0.0,0.0,0.0,0.0,0.001,0.0];
                let report=owner.step(&mut thermal,&force,controls(),&[0.0],&gate).unwrap();
                let mut after=[0.0;6];owner.reconstruct_into(&mut after,&mut velocity,&gate).unwrap();
                let nodal:f64=force.iter().zip(after.iter().zip(before)).map(|(f,(b,a))|f*(b-a)).sum();
                near(nodal,report.mechanical.external_work_j,1e-11);
                let t=thermal.temperature_k(CELL).unwrap();
                near(owner.members()[0].source.layers[0].young_modulus_pa,YOUNG+slope*(t-300.0),1e-6);
                work+=nodal;heat+=report.external_thermal_heat_j;outside+=report.external_mechanical_loss_j;
                near(independent_total(&owner,&thermal,slope,alpha)-initial+outside-work-heat,0.0,1e-7);
                near(report.total_energy_residual_j,0.0,2e-8);
            }
            let (mut u,mut v)=([0.0;6],[0.0;6]);owner.reconstruct_into(&mut u,&mut v,&gate).unwrap();
            endpoints.push(u);
            assert!(outside>0.0);
            assert!(owner.physical_state().unwrap()[4..].iter().any(|p| *p!=0.0));
        });
    }
    assert!(endpoints[1][3]>1e-8);
    assert!((endpoints[1][4]-endpoints[0][4]).abs()>1e-10);
    near(endpoints[0][3],0.0,1e-12);
}

#[test]
fn late_thermal_refusal_and_cancellation_leave_complete_histories_for_retry() {
    fixture(-2000.0,1e-4,10.0,true,|model,mut thermal,source| {
        let state=initial(model);let gate=CancelGate::new_clock_free();
        let mut owner=FrameFluidThermoelasticSystem::new(model,&state,&thermal,source,&DESTINATIONS,config(),&gate).unwrap();
        let accepted=owner.physical_state().unwrap();let cells=thermal.cells().clone();
        let mut exhausted=controls();exhausted.maximum_material_trials=1;
        assert!(owner.step(&mut thermal,&[0.0;6],exhausted,&[0.0],&gate).is_err());
        assert_eq!(owner.physical_state().unwrap(),accepted);assert_eq!(thermal.cells(),&cells);
        assert_eq!(thermal.time_s(),0.0);assert_eq!(owner.mechanical().samples_rendered(),0);
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(owner.step(&mut thermal,&[0.0;6],controls(),&[0.0],&cancelled).is_err());
        assert_eq!(owner.physical_state().unwrap(),accepted);
        owner.step(&mut thermal,&[0.0;6],controls(),&[0.0],&gate).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);near(thermal.time_s(),0.001,0.0);
    });
}

#[test]
fn empty_fluid_binding_retains_original_dry_thermoelastic_trajectory() {
    fixture(-2000.0,1e-4,10.0,false,|model,mut thermal,source| {
        let gate=CancelGate::new_clock_free();let state=initial(model);let mut direct_heat=thermal.clone();
        let mut component=model.basis().coupled_component(1000,&[0.0;2],config().modal,4096,||false).unwrap();
        component.restore_states(&model.modal_states(&state).unwrap()).unwrap();
        let mut direct=component.into_thermoelastic_system(&direct_heat,source,&DESTINATIONS,0.05,
            config().coupling,None,&gate).unwrap();
        let mut wet=FrameFluidThermoelasticSystem::new(model,&state,&thermal,source,&DESTINATIONS,config(),&gate).unwrap();
        for _ in 0..8 {
            wet.step(&mut thermal,&[0.0;6],controls(),&[0.0],&gate).unwrap();
            direct.step(&mut direct_heat,&[0.0;2],controls(),&[0.0],&gate).unwrap();
            for (a,b) in wet.modal_states().unwrap().iter().zip(direct.mechanical().components()[0].states()) {
                near(a.displacement_m_sqrt_kg,b.displacement_m_sqrt_kg,1e-12);
                near(a.velocity_m_sqrt_kg_per_s,b.velocity_m_sqrt_kg_per_s,1e-12);
            }
            near(thermal.temperature_k(CELL).unwrap(),direct_heat.temperature_k(CELL).unwrap(),1e-10);
        }
    });
}

#[test]
fn original_resistor_heat_drives_expansion_on_same_wet_frame_clock() {
    fixture(-2000.0,1e-4,0.0,true,|model,mut thermal,thermal_source| {
        let gate=CancelGate::new_clock_free();
        let state=model.state_from_motion(&[ModalAcousticState::default();2],&[0.0;8]).unwrap();
        let resistors=[CellResistor { destination:ResistorCell {branch:0,cell:CELL},
            selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area() }];
        let source=MaterialNetwork { cells:thermal_source.cells,links:thermal_source.links,resistors:&resistors };
        let resolved=source.resolve(&thermal,||false).unwrap();
        let circuit=MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
            (1,0,Branch::Resistor {ohms:999.0}), (0,1,Branch::CurrentSource {port:0})],transformers:vec![]},
            resolved.resistors,vec![],vec![],&[1.0]).unwrap();
        let mut owner=FrameFluidThermoelasticSystem::new(model,&state,&thermal,thermal_source,&DESTINATIONS,config(),&gate).unwrap()
            .with_material_circuit(&thermal,&source,circuit,&[(0,ElectricalPropertyBasis::Resistivity)],&gate).unwrap();
        assert!(owner.step(&mut thermal,&[0.0;6],controls(),&[0.0],&gate).is_err());
        for _ in 0..8 {
            let report=owner.step_electrical(&mut thermal,&[0.0;6],FrameElectricalInput {
                ports:&[1.0],solve:MaterialCircuitSolve {energy_tolerance_j:1e-9,power_tolerance_w:1e-9}},
                controls(),&[0.0],&gate).unwrap();
            near(report.deposited_joule_heat_j,0.001*0.02/area(),1e-9);
            assert!(report.electrical.is_some());near(report.total_energy_residual_j,0.0,2e-8);
        }
        let (mut u,mut v)=([0.0;6],[0.0;6]);owner.reconstruct_into(&mut u,&mut v,&gate).unwrap();
        assert!(u[3]>0.0 && thermal.temperature_k(CELL).unwrap()>300.0);
        assert_eq!(owner.mechanical().samples_rendered(),8);
    });
}
