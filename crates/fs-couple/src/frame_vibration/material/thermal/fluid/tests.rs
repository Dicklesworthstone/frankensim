//! Original source/enthalpy/prepared-mechanical owners, with manufactured data.
use super::*;
use crate::air_path::cylinder::frame::viscosity::tests::{
    frame, model, near, point, solid_card,
};
use crate::circuit_enthalpy::material::MaterialCell;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::modal_acoustic_time::ModalAcousticState;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, SPECIFIC_HEAT_CAPACITY_DIMS, resolve_material_state_point};
use fs_qty::{Density, DynViscosity, Pressure, Temperature};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, Geometry};
use fs_thermal::enthalpy::*;

const CELL: CellId = CellId(41);
const RADIUS: f64 = 0.02;
fn exterior() -> FluidTransportState {
    // A declared numerical Newtonian transport fixture, not a measured liquid.
    FluidTransportState::try_new(400.0, 4000.0, 0.5, 0.1).unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &heat, density: &density, temperature_axis: "T",
        lower: &point(290.0), upper: &point(500.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap()
}
fn controls() -> FrameThermalSolve {
    FrameThermalSolve {
        maximum_cells: 4, maximum_material_trials: 200, maximum_newton_iterations: 16,
        difference_step_j: 1e-5, cell_energy_tolerance_j: 2e-9,
        thermal_energy_tolerance_j: 2e-9, total_energy_tolerance_j: 4e-9,
        newton: NewtonKrylovConfig {absolute_tolerance: 1e-10, relative_tolerance: 1e-13,
            ..NewtonKrylovConfig::default()},
    }
}
fn initial<'a>(curves: &[&'a ResolvedHeatCapacityCurve], radii: &[f64],
    powers: &[f64], conductance: f64) -> DepositionSimulation<'a> {
    assert_eq!(curves.len(), radii.len());
    assert_eq!(curves.len(), powers.len());
    let mut inner = 0.0;
    let mut events = Vec::new();
    for (i, (&radius, &power)) in radii.iter().zip(powers).enumerate() {
        let volume = core::f64::consts::PI*(radius-inner)*(radius+inner);
        // Unit length, explicit annular surface; no hidden finite fluid volume.
        let surface = 2.0*core::f64::consts::PI*(radius+inner)+2.0*volume;
        events.push(TimedEvent {time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(41+i as u64), segment_id: i as u64, body_id: 0,
            material_id: MaterialId(i as u64), mass_kg: 1000.0*volume,
            specific_h_j_kg: 100.0,
            geometry: CellGeometry {start_m: [0.0; 3], end_m: [1.0,0.0,0.0],
                volume_m3: volume, surface_area_m2: surface},
            boundary: Boundary {ambient_k: 300.0, htc_w_m2_k: conductance/surface,
                reservoir: None, source_w: power},
        })});
        inner=radius;
    }
    let laws: Vec<_> = curves.iter().enumerate().map(|(i,c)|
        (MaterialId(i as u64), c.curve() as &dyn ThermalLaw)).collect();
    DepositionSimulation::new(events,&laws,DepositionBudget {max_cells:4,max_links:4,
        max_events:8,max_substeps:1000}).unwrap().advance_to(0.0,
        StepConfig {max_dt_s:0.01,energy_tolerance_j:1e-10},||Ok(())).unwrap().0
}
fn binding<'a>(i: usize, card: &'a MaterialCard, curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint) -> MaterialCell<'a> {
    MaterialCell {cell:CellId(41+i as u64),material_id:MaterialId(i as u64),curve:curve.curve(),
        point:ThermalMaterialPoint {material:card,context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:None}
}
fn routes(count:usize) -> Vec<FrameLayerCell> {
    (0..count).map(|i| FrameLayerCell {element:0,layer:i,cell:CellId(41+i as u64)}).collect()
}
fn motion<F>(model:&MaterialFluidFrame<'_,F>) -> Vec<f64> {
    let phi=model.basis().modes()[0].displacement[4];
    model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:0.002/phi,
        velocity_m_sqrt_kg_per_s:0.03/phi}],&[0.001;12]).unwrap()
}

#[test]
fn actual_solid_heat_and_exterior_sink_close_the_complete_native_history() {
    let card=solid_card(0.0); let f=frame(&card); let curve=caloric(&card); let context=point(300.0);
    let mut thermal=initial(&[&curve],&[RADIUS],&[0.0],0.0);
    let cells=[binding(0,&card,&curve,&context)];
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let model=model(&f,exterior()); let state=motion(&model); let gate=CancelGate::new_clock_free();
    let saved_frame=f.displacements().to_vec();
    let initial_energy=model.system().hamiltonian(&state);
    let initial_h=thermal.cells()[&CELL].specific_h_j_kg;
    let mass=thermal.cells()[&CELL].mass_kg;
    let mut coupled=FrameFluidThermalSystem::new(&model,state,&thermal,&network,&routes(1),&gate).unwrap();
    let (mut solid,mut fluid,mut work)=(0.0,0.0,0.0);
    for step in 0..30 {
        let force=if step<10 {0.03} else {-0.01};
        let load=[0.0,0.0,0.0,0.0,force,0.0];
        let frozen=model.step_candidate(coupled.state(),&load,0.005,1e-9,&gate).unwrap();
        let target=coupled.time_s()+0.005;
        let report=coupled.step(&mut thermal,target,&load,controls(),&[0.0],&gate).unwrap();
        for (a,b) in coupled.state().iter().zip(frozen.native.x) {near(*a,b,3e-9);}
        assert_eq!(coupled.time_s().to_bits(),thermal.time_s().to_bits());
        near(report.deposited_solid_heat_j,report.layer_heat[0].total_j(),1e-16);
        near(report.exported_fluid_heat_j,report.mechanical.fluid_member_dissipation_j[0],1e-16);
        assert!(report.entropy.material_loss_generation_j_k>0.0);
        solid+=report.deposited_solid_heat_j;fluid+=report.exported_fluid_heat_j;
        work+=report.mechanical.nodal_work_j;
    }
    let thermal_change=mass*(thermal.cells()[&CELL].specific_h_j_kg-initial_h);
    near(thermal_change,solid,1e-8);
    near(model.system().hamiltonian(coupled.state())-initial_energy+thermal_change+fluid,work,1e-8);
    assert!(solid>0.0 && fluid>100.0*solid);
    assert!(thermal.temperature_k(CELL).unwrap()>300.0);
    assert!(coupled.state()[2..].iter().all(|p| *p!=0.0));
    assert_eq!(f.displacements(),saved_frame);
}

#[test]
fn source_only_temperature_dependence_reaches_motion_in_the_same_accepted_interval() {
    let gate=CancelGate::new_clock_free(); let mut outputs=Vec::new();
    for alpha in [0.0,0.05] {
        let card=solid_card(alpha); let f=frame(&card); let curve=caloric(&card); let context=point(300.0);
        let mut thermal=initial(&[&curve],&[RADIUS],&[1000.0],2.0);
        let cells=[binding(0,&card,&curve,&context)];
        let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
        let model=model(&f,exterior());
        let original=model.frame().material_viscous_section(0).unwrap().bending().viscous_stiffness_n_m2_s();
        let mut coupled=FrameFluidThermalSystem::new(&model,motion(&model),&thermal,&network,&routes(1),&gate).unwrap();
        for step in 0..30 {
            let force=if step<10 {0.03} else {-0.01};
            let target=coupled.time_s()+0.005;
            let report=coupled.step(&mut thermal,target,&[0.0,0.0,0.0,0.0,force,0.0],
                controls(),&[5.0],&gate).unwrap();
            let temperature=thermal.temperature_k(CELL).unwrap();
            let viscosity=report.viscosity[0].1.bending().viscous_stiffness_n_m2_s();
            near(viscosity/original,1.0+alpha*(temperature-300.0),2e-12);
            let layer=&report.viscosity[0].1.elastic().layers()[0].material;
            assert!(layer.query_point().contains(&("T".into(),temperature)));
            if alpha>0.0 {assert!(viscosity>original);}
        }
        let u=coupled.state()[0]*model.basis().modes()[0].displacement[4];
        outputs.push(u);
        // Independent physical-coordinate midpoint/enthalpy equations, not the
        // modal projection or native solver used above. Source Bessel rounding
        // and nonlinear tolerances are allowed explicitly in these comparisons.
        near(thermal.temperature_k(CELL).unwrap(),311.79062998,2e-7);
        let expected=if alpha==0.0 {0.005645025955270403} else {0.005644968621268697};
        near(u,expected,2e-9);
    }
    assert!(outputs[0]-outputs[1]>1e-9);
}

fn layered_frame(cards:&[MaterialCard],radii:&[f64]) -> Frame2d {
    let requirements:Vec<_>=[("density",Density::DIMS),("young_modulus",Pressure::DIMS),
        ("axial_eta",DynViscosity::DIMS),("bending_eta",DynViscosity::DIMS)].into_iter()
        .map(|(n,d)|ScalarPropertyRequirement::try_new(n,d,ScalarAdmissibility::Finite).unwrap()).collect();
    let layers=cards.iter().zip(radii).map(|(card,&outer_radius_m)|ConcentricMaterialLayer {
        outer_radius_m,material:resolve_material_state_point(card,&point(300.0),&requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()}).collect();
    let section=ConcentricElasticSection::try_new(0.0,layers).unwrap();
    let section=FrameViscousSection::try_new(section,"axial_eta","bending_eta").unwrap();
    let mut frame=Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0,0.0);frame.add_node(1.0,0.0);frame.fix(0,[true;3]);frame.fix(1,[true,false,true]);
    frame.add_material_viscous_element(0,1,3,section).unwrap();frame
}

#[test]
fn separate_layer_temperatures_use_fresh_source_heat_partitions() {
    let cards=[solid_card(0.05),solid_card(0.0)];
    let radii=[RADIUS/2.0_f64.sqrt(),RADIUS];let f=layered_frame(&cards,&radii);
    let curves=[caloric(&cards[0]),caloric(&cards[1])];let context=point(300.0);
    let mut thermal=initial(&[&curves[0],&curves[1]],&radii,&[1000.0,0.0],0.0);
    let cells=[binding(0,&cards[0],&curves[0],&context),binding(1,&cards[1],&curves[1],&context)];
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};let model=model(&f,exterior());
    let gate=CancelGate::new_clock_free();
    let mut coupled=FrameFluidThermalSystem::new(&model,motion(&model),&thermal,&network,&routes(2),&gate).unwrap();
    let mut last_share=0.0;
    for _ in 0..20 {
        let target=coupled.time_s()+0.005;
        let report=coupled.step(&mut thermal,target,&[0.0;6],controls(),&[5.0,0.0],&gate).unwrap();
        let section=&report.viscosity[0].1;
        last_share=section.bending().layer_viscous_stiffness_n_m2_s()[0]/section.bending().viscous_stiffness_n_m2_s();
        let total=report.layer_heat.iter().map(|h|h.bending_heat_j).sum::<f64>();
        assert!(total>0.0);
        near(report.layer_heat[0].bending_heat_j/total,last_share,2e-12);
        near(report.layer_heat.iter().map(|h|h.total_j()).sum(),report.deposited_solid_heat_j,1e-15);
    }
    assert!(last_share>0.35); // Cold core fraction is exactly 0.25, not the hot fraction.
    assert!(thermal.temperature_k(CELL).unwrap()>thermal.temperature_k(CellId(42)).unwrap()+10.0);
}

#[test]
fn exhausted_joint_trials_and_late_source_refusal_publish_neither_history() {
    let card=solid_card(0.05);let f=frame(&card);let curve=caloric(&card);let context=point(300.0);
    let mut thermal=initial(&[&curve],&[RADIUS],&[1e6],0.0);
    let cells=[binding(0,&card,&curve,&context)];let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let model=model(&f,exterior());let gate=CancelGate::new_clock_free();
    let mut coupled=FrameFluidThermalSystem::new(&model,motion(&model),&thermal,&network,&routes(1),&gate).unwrap();
    let saved=coupled.state().to_vec();let before=thermal.cells()[&CELL].specific_h_j_kg.to_bits();
    let limited=FrameThermalSolve {maximum_material_trials:1,..controls()};
    assert!(matches!(coupled.step(&mut thermal,0.0001,&[0.0;6],limited,&[0.0],&gate),
        Err(CircuitEnthalpyError::CircuitBudget)));
    // Initial mechanical trial is valid, but the joint target would require
    // heating beyond the ORIGINAL caloric/material support. No fallback source.
    assert!(coupled.step(&mut thermal,0.01,&[0.0;6],controls(),&[0.0],&gate).is_err());
    assert_eq!(coupled.state(),saved);assert_eq!(coupled.time_s(),0.0);
    assert_eq!(thermal.time_s(),0.0);assert_eq!(thermal.cells()[&CELL].specific_h_j_kg.to_bits(),before);
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(coupled.step(&mut thermal,0.0001,&[0.0;6],controls(),&[0.0],&cancelled),
        Err(CircuitEnthalpyError::Cancelled)));
    let report=coupled.step(&mut thermal,0.0001,&[0.0;6],controls(),&[100.0],&gate).unwrap();
    assert!(report.deposited_solid_heat_j>0.0);
    assert!(thermal.temperature_k(CELL).unwrap()>307.0);
    assert_eq!(coupled.time_s(),thermal.time_s());
}

#[test]
fn thermal_inventory_must_be_the_original_whole_material_layer() {
    let card=solid_card(0.0);let f=frame(&card);let curve=caloric(&card);let context=point(300.0);
    let thermal=initial(&[&curve],&[RADIUS*1.1_f64.sqrt()],&[0.0],0.0);
    let cells=[binding(0,&card,&curve,&context)];let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let model=model(&f,exterior());let state=motion(&model);let saved=state.clone();
    assert!(FrameFluidThermalSystem::new(&model,state,&thermal,&network,&routes(1),
        &CancelGate::new_clock_free()).is_err());
    assert_eq!(saved,motion(&model));assert_eq!(thermal.time_s(),0.0);
}

fn zero_reference_card() -> MaterialCard {
    let mut claims=ClaimSet::new();
    for (name,dims,value) in [("density",Density::DIMS,1000.0),("young_modulus",Pressure::DIMS,2e6),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,10.0),
        ("axial_eta",DynViscosity::DIMS,0.0),("bending_eta",DynViscosity::DIMS,0.0)] {
        let viscosity=name.ends_with("_eta");
        let mut validity=ValidityDomain::unconstrained().with("T",290.0,500.0);
        if viscosity {validity=validity.with("omega",0.0,1e7);}
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:if viscosity {PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(290.0,0.0),(300.0,0.0),(500.0,5000.0)],dims}}
                else {PropertyValue::Scalar {value,dims}},
            validity,interpolation:if viscosity {InterpolationPolicy::LinearInside}
                else {InterpolationPolicy::ConstantWithinValidity},
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"synthetic zero-to-positive viscosity".into(),license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic switching viscosity".into(),phase:"solid".into(),
        process:"continuous source curve".into(),revision:0},claims,vec![]).unwrap()
}

#[test]
fn zero_reference_viscosity_can_become_positive_without_a_divide_by_zero_override() {
    let card=zero_reference_card();let f=frame(&card);let curve=caloric(&card);let context=point(300.0);
    let mut thermal=initial(&[&curve],&[RADIUS],&[10000.0],0.0);
    let cells=[binding(0,&card,&curve,&context)];let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let model=model(&f,exterior());let gate=CancelGate::new_clock_free();let state=motion(&model);
    let frozen=model.step_candidate(&state,&[0.0;6],0.005,1e-9,&gate).unwrap();
    assert_eq!(frozen.material_port_dissipation_j.iter().sum::<f64>(),0.0);
    let mut coupled=FrameFluidThermalSystem::new(&model,state,&thermal,&network,&routes(1),&gate).unwrap();
    let report=coupled.step(&mut thermal,0.005,&[0.0;6],controls(),&[50.0],&gate).unwrap();
    assert!(report.deposited_solid_heat_j>0.0);
    assert!(report.viscosity[0].1.bending().viscous_stiffness_n_m2_s()>0.0);
    assert!(report.mechanical.state[2..].iter().any(|p| *p!=0.0));
}
