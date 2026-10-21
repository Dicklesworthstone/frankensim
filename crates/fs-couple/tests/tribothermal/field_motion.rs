//! Native field work, source contact traction and finite heat share one step.
use super::*;
#[path = "../electrical_contacts/source.rs"]
mod electrical_source;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitThermalCoefficients};
use fs_couple::circuit_enthalpy::contacts::{ContactHeatCells, ContactThermalCoefficients};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialResistor};
use fs_couple::material_circuit::contacts::SeriesElectricalContact;
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_couple::tribothermal::field_motion::{Controls, Input};
use fs_ga::Vec3;
use fs_material::field_storage::{ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS,ResolvedCapacitor};
use fs_material::electrical_interface::resolve_electrical_interface_state;
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch,CircuitGraph};
use fs_thermal::enthalpy::{DepositedCell,DepositionBudget,Event,MaterialId,StepConfig,ThermalLaw,TimedEvent};
use fs_tribo::FrictionRegime;

fn field_controls() -> Controls {
    let mut sliding=controls();
    // Compliant return maps also resolve the outer displacement Jacobian's
    // perturbed trajectories; prescribed-motion iteration limits are too low.
    sliding.tangential.max_iterations=64;
    sliding.tangential.force_absolute_tolerance_n=1e-12;
    sliding.tangential.force_relative_tolerance=1e-12;
    Controls { sliding,circuit:MaterialCircuitSolve {
        energy_tolerance_j:1e-8,power_tolerance_w:1e-8 },maximum_motion_trials:100,
        maximum_field_trials:10000,maximum_newton_iterations:20,difference_step_m:1e-6,
        displacement_tolerance_m:1e-8,energy_tolerance_j:2e-6,
        newton:NewtonKrylovConfig { absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:2,..NewtonKrylovConfig::default() } }
}
fn empty(_: &DepositionSimulation<'_>) -> Result<ContactThermalCoefficients<()>,Infallible> {
    Ok(ContactThermalCoefficients { bulk:CircuitThermalCoefficients {
        resistors:vec![],sources_w:BTreeMap::new(),conductances_w_k:BTreeMap::new(),source:() },contacts:vec![] })
}
fn actuator(material: &MaterialCard,charge: f64,velocity: f64,
    resistors: Vec<MaterialResistor>,contacts: Vec<SeriesElectricalContact>) -> MaterialFieldMotion {
    let dielectric = electrical_source::scalar_card(ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS,1.0);
    let storage = MaterialStorage { branch:0,element:MaterialStorageElement::Capacitor(
        ResolvedCapacitor::resolve(&dielectric,&QueryPoint::new(),MaterialPropertySelection::SingleClaimOnly,
            1.0,1.0).unwrap()) };
    let mut branches=vec![(1,0,Branch::Capacitor { farads:f64::NAN })];
    let nodes=if resistors.is_empty() { 2 } else {
        branches.extend([(1,2,Branch::Resistor { ohms:f64::NAN }),(2,0,Branch::Resistor { ohms:f64::NAN })]); 3 };
    let circuit=MaterialCircuit::new_with_storage(CircuitGraph { node_count:nodes,branches,
        transformers:vec![] },resistors,vec![storage],vec![],vec![charge],&[]).unwrap()
        .rebind_contacts_candidate(contacts,&[]).unwrap();
    let body=resolve_body(&[RegionBinding { key:"moving surface A",material,point:&common::point(300.0),
        selection:MaterialPropertySelection::SingleClaimOnly,
        geometry:VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap() }],1,||false).unwrap();
    MaterialFieldMotion::new(circuit,0,body,velocity,[0.5,1.5]).unwrap()
}
fn input(motion:&MaterialFieldMotion,normal:f64) -> Input<'_> {
    Input { motion,electrical_ports:&[],external_force_n:0.0,normal_load_n:normal,
        resistor_cells:&[],contact_cells:&[] }
}
fn initial_with<'a>(curves:&'a [ResolvedHeatCapacityCurve;2],edit:impl Fn(&mut DepositedCell))
    -> DepositionSimulation<'a> {
    let original=common::initial(curves,0.0,None);
    let mut schedule=original.cells().values().cloned().map(|mut cell| {
        edit(&mut cell); TimedEvent {time_s:0.0,event:Event::Birth(cell)}
    }).collect::<Vec<_>>();
    schedule.extend(original.links().values().cloned().map(|link|
        TimedEvent {time_s:0.0,event:Event::Contact(link)}));
    DepositionSimulation::new(schedule,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],DepositionBudget {
        max_cells:2,max_links:1,max_events:3,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-7},||Ok(())).unwrap().0
}

#[test]
fn both_sliding_directions_match_maxwell_force_and_original_friction_work() {
    let (a,ca)=common::material(1.0,2.0,0.0); let (b,cb)=common::material(1.0,3.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=context();
    let card=interface_with_static(&cards,0.3,"speed",SPEED,[(0.0,0.2),(2.0,0.2)],2000.0);
    let thermal=common::initial(&curves,0.0,None);
    for velocity in [-1.0_f64,1.0] {
        let contact=MaterialSlidingContact::new(&thermal,binding(&cards,&curves,&card,&point,
            HeatPartition::new(0.4,0.6,0.0).unwrap()),1.0,1.0,InputAuthority::SyntheticFixture).unwrap();
        let motion=actuator(&cards[0],0.2,velocity,vec![],vec![]); let before=motion.clone();
        let next=contact.field_motion_candidate(&thermal,0.1,input(&motion,1.0),field_controls(),
            [0.1*velocity,0.1*velocity],&[0.0;2],empty,||false).unwrap();
        let traction=0.2*velocity.signum(); let acceleration=-0.5*0.2*0.2-traction;
        let v1=velocity+0.1*acceleration; let displacement=0.1*(velocity+v1)/2.0;
        let heat=traction*displacement; let report=next.report();
        close(next.motion().gap_m(),1.0+displacement); close(next.motion().velocity_m_s(),v1);
        close(report.sliding.friction.traction_n[0],traction);
        close(report.sliding.friction.slip_increment_m[0],displacement);
        close(report.sliding.material.speed_mps,displacement.abs()/0.1);
        close(report.sliding.deposited_heat_j,heat);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(),300.0+0.4*heat/2.0);
        close(next.thermal().temperature_k(CELLS[1]).unwrap(),300.0+0.6*heat/3.0);
        close(report.field.mechanical_supply_j,-heat);
        close(report.field.kinetic_energy_j[1]-report.field.kinetic_energy_j[0]
            +report.field.electrical_energy_j[1]-report.field.electrical_energy_j[0]+heat,0.0);
        close(report.contact_transfer_residual_j,0.0); close(report.total_energy_residual_j,0.0);
        assert_eq!(report.field.heat_j,0.0); assert_eq!(report.deposited_joule_j,0.0);
        assert_eq!(report.external_heat_j,0.0); assert_eq!(motion,before);
        assert_eq!(contact.state().slip_m(),[0.0;2]);
    }
}

#[test]
fn static_capacity_holds_field_force_without_creating_motion_or_heat() {
    let (a,ca)=common::material(1.0,2.0,0.0); let (b,cb)=common::material(1.0,3.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=context();
    let card=interface_with_static(&cards,0.8,"speed",SPEED,[(0.0,0.2),(2.0,0.2)],2000.0);
    let thermal=common::initial(&curves,0.0,None);
    let contact=MaterialSlidingContact::new(&thermal,binding(&cards,&curves,&card,&point,
        HeatPartition::new(0.5,0.5,0.0).unwrap()),1.0,0.0,InputAuthority::SyntheticFixture).unwrap();
    let motion=actuator(&cards[0],1.0,0.0,vec![],vec![]);
    for force in [0.0,1.0] {
        // The tiny query seed is below Newton's tolerance: accepted sticking
        // must still re-query EXACT zero speed/distance, in either direction.
        let next=contact.field_motion_candidate(&thermal,0.1,Input {
            external_force_n:force,..input(&motion,1.0) },field_controls(),
            [0.005*(force-0.5),1e-14],&[0.0;2],empty,||false).unwrap();
        assert_eq!(next.report().sliding.friction.regime,FrictionRegime::Sticking);
        close(next.report().sliding.friction.traction_n[0],force-0.5);
        close(next.motion().gap_m(),1.0); close(next.motion().velocity_m_s(),0.0);
        assert_eq!(next.report().sliding.material.speed_mps,0.0);
        assert_eq!(next.report().sliding.material.accumulated_distance_m,0.0);
        assert_eq!(next.contact().state().slip_m(),[0.0;2]);
        assert_eq!(next.report().sliding.deposited_heat_j,0.0);
        assert_eq!(next.thermal().cells(),thermal.cells());
    }
}

#[test]
fn body_source_mass_and_reference_volume_must_match_the_actual_thermal_surface() {
    let (a,ca)=common::material(1.0,2.0,0.0); let (b,cb)=common::material(1.0,3.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=context();
    let card=interface_with_static(&cards,0.3,"speed",SPEED,[(0.0,0.2),(2.0,0.2)],2000.0);
    for (mass,volume,source) in [(2.0,1.0,0),(1.0,2.0,0),(1.0,1.0,1)] {
        let thermal=initial_with(&curves,|cell| if cell.id==CELLS[0] {
            cell.mass_kg=mass; cell.geometry.volume_m3=volume;
        });
        let contact=MaterialSlidingContact::new(&thermal,binding(&cards,&curves,&card,&point,
            HeatPartition::new(0.5,0.5,0.0).unwrap()),1.0,1.0,InputAuthority::SyntheticFixture).unwrap();
        let motion=actuator(&cards[source],0.2,1.0,vec![],vec![]);
        assert_eq!(motion.body().regions()[0].inertia().mass(),1.0);
        assert_eq!(motion.body().regions()[0].geometry().volume_m3(),1.0);
        let calls=Cell::new(0usize);
        assert!(matches!(contact.field_motion_candidate(&thermal,0.1,input(&motion,1.0),
            field_controls(),[0.1;2],&[0.0;2],|state| {calls.set(calls.get()+1);empty(state)},||false),
            Err(SlidingError::Invalid("surface A must be the same complete single-region field body, source, mass and reference volume"))));
        assert_eq!(calls.get(),0); assert_eq!(contact.accumulated_distance_m(),0.0);
    }
}

#[test]
fn authored_boundary_and_callback_heat_are_separate_from_internal_losses_and_force_work() {
    let (a,ca)=common::material(1.0,2.0,0.0); let (b,cb)=common::material(1.0,3.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=context();
    let card=interface_with_static(&cards,0.3,"speed",SPEED,[(0.0,0.2),(2.0,0.2)],2000.0);
    let thermal=initial_with(&curves,|cell| if cell.id==CELLS[0] {cell.boundary.source_w=3.0;});
    let contact=MaterialSlidingContact::new(&thermal,binding(&cards,&curves,&card,&point,
        HeatPartition::new(0.4,0.6,0.0).unwrap()),1.0,1.0,InputAuthority::SyntheticFixture).unwrap();
    let motion=actuator(&cards[0],0.2,1.0,vec![],vec![]);
    let next=contact.field_motion_candidate(&thermal,0.1,Input {external_force_n:0.4,
        ..input(&motion,1.0)},field_controls(),[0.1;2],&[0.0;2],|state| {
            let mut source=empty(state)?; source.bulk.sources_w.insert(CELLS[1],2.0); Ok::<_,Infallible>(source)
        },||false).unwrap();
    let displacement=0.1*(1.0+(1.0+0.1*(0.4-0.02-0.2)))/2.0;
    let heat=0.2*displacement; let report=next.report();
    close(next.motion().gap_m(),1.0+displacement);
    close(report.external_heat_j,0.5); close(report.external_mechanical_work_j,0.4*displacement);
    close(report.sliding.thermal.stored_energy_change_j,0.5+heat);
    close(next.thermal().temperature_k(CELLS[0]).unwrap(),300.0+(0.3+0.4*heat)/2.0);
    close(next.thermal().temperature_k(CELLS[1]).unwrap(),300.0+(0.2+0.6*heat)/3.0);
    close(report.total_energy_residual_j,0.0);
    close(report.field.electrical_energy_j[1]-report.field.electrical_energy_j[0]
        +report.field.kinetic_energy_j[1]-report.field.kinetic_energy_j[0]
        +report.sliding.thermal.stored_energy_change_j,
        report.external_heat_j+report.external_mechanical_work_j);
}

#[test]
fn caloric_data_changes_field_motion_through_same_trial_friction_and_electrical_contacts() {
    let run=|cp| {
        let (a,ca)=common::material(1.0,cp,0.005); let (b,cb)=common::material(2.0,2.0,0.0);
        let cards=[a,b]; let curves=[ca,cb]; let point=context();
        let card=interface_with_static(&cards,1.5,"Ta",Temperature::DIMS,[(299.0,0.1),(310.0,1.2)],2000.0);
        let thermal=common::initial(&curves,0.0,None);
        let contact=MaterialSlidingContact::new(&thermal,binding(&cards,&curves,&card,&point,
            HeatPartition::new(0.5,0.5,0.0).unwrap()),1.0,1.0,InputAuthority::SyntheticFixture).unwrap();
        let electrical=electrical_source::card(&cards[0],&cards[1],0.5,0.01);
        let resolve=|state:&DepositionSimulation<'_>| {
            let query=electrical_source::point().with("Ta",state.temperature_k(CELLS[0]).unwrap()).unwrap()
                .with("Tb",state.temperature_k(CELLS[1]).unwrap()).unwrap();
            Ok::<_,Infallible>(ContactThermalCoefficients {
                bulk:common::coefficients(state,&cards,&curves).unwrap(),
                contacts:vec![SeriesElectricalContact::new(7,1,resolve_electrical_interface_state(
                    &electrical,&query,MaterialPropertySelection::SingleClaimOnly).unwrap(),1.0).unwrap()] })
        };
        let initial=resolve(&thermal).unwrap();
        let motion=actuator(&cards[0],2.0,1.0,initial.bulk.resistors,initial.contacts);
        let routes=[ContactHeatCells { id:7,cells:CELLS,heat_fraction_to_first:0.25 }];
        let next=contact.field_motion_candidate(&thermal,0.05,Input { resistor_cells:&common::DESTINATIONS,
            contact_cells:&routes,..input(&motion,1.0) },field_controls(),[0.05,0.047],&[0.0;2],resolve,||false).unwrap();
        let r=next.report(); let ta=next.thermal().temperature_k(CELLS[0]).unwrap();
        close(r.sliding.friction.traction_n[0],0.2+0.1*(ta-300.0));
        close(r.field.contacts[0].contact.resistance_ohm(),0.5*(1.0+0.01*(ta-300.0)));
        close(r.field.resistors[0].material.conductor.resistance_ohm(),1.0+0.005*(ta-300.0));
        let bulk=r.field.resistors.iter().map(|r|r.heat_j).sum::<f64>();
        close(r.deposited_joule_j,bulk+r.field.contacts[0].heat_j);
        close(r.sliding.thermal.stored_energy_change_j,r.deposited_joule_j+r.sliding.deposited_heat_j);
        close(r.total_energy_residual_j,0.0); close(r.external_heat_j,0.0);
        close(next.motion().circuit().charges_c()[0].powi(2)*next.motion().gap_m()/2.0
            +next.motion().velocity_m_s().powi(2)/2.0
            +next.thermal().energy_j()-thermal.energy_j(),2.5);
        (ta,next.motion().velocity_m_s(),next.motion().circuit().charges_c()[0])
    };
    let hotter=run(0.01); let cooler=run(0.2);
    assert!(hotter.0>cooler.0+1.0);
    assert!(hotter.1<cooler.1-1e-3,"actual source friction must alter body momentum");
    assert!((hotter.2-cooler.2).abs()>1e-5,"same source temperatures must alter circuit history");
}

#[test]
fn cancelled_exhausted_or_reversing_field_trajectory_publishes_no_partial_history() {
    let (a,ca)=common::material(1.0,2.0,0.0); let (b,cb)=common::material(1.0,3.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=context();
    let card=interface_with_static(&cards,0.3,"speed",SPEED,[(0.0,0.2),(2.0,0.2)],2000.0);
    let thermal=common::initial(&curves,0.0,None); let before=thermal.cells().clone();
    let contact=MaterialSlidingContact::new(&thermal,binding(&cards,&curves,&card,&point,
        HeatPartition::new(0.5,0.5,0.0).unwrap()),1.0,1.0,InputAuthority::SyntheticFixture).unwrap();
    let motion=actuator(&cards[0],0.2,1.0,vec![],vec![]); let original=motion.clone();
    let polls=Cell::new(0usize);
    let good=contact.field_motion_candidate(&thermal,0.1,input(&motion,1.0),field_controls(),
        [0.1;2],&[0.0;2],empty,|| {polls.set(polls.get()+1);false}).unwrap();
    let last=polls.get(); polls.set(0);
    assert!(matches!(contact.field_motion_candidate(&thermal,0.1,input(&motion,1.0),field_controls(),
        [0.1;2],&[0.0;2],empty,|| {polls.set(polls.get()+1);polls.get()==last}),Err(SlidingError::Cancelled)));
    let mut controls=field_controls(); controls.maximum_field_trials=good.report().field_trials-1;
    assert!(matches!(contact.field_motion_candidate(&thermal,0.1,input(&motion,1.0),controls,
        [0.1;2],&[0.0;2],empty,||false),Err(SlidingError::Source(CircuitEnthalpyError::CircuitBudget))));
    controls=field_controls(); controls.maximum_motion_trials=good.report().motion_trials-1;
    assert!(matches!(contact.field_motion_candidate(&thermal,0.1,input(&motion,1.0),controls,
        [0.1;2],&[0.0;2],empty,||false),Err(SlidingError::MotionBudget)));
    let slow=actuator(&cards[0],0.2,0.01,vec![],vec![]);
    assert!(contact.field_motion_candidate(&thermal,0.5,input(&slow,1.0),field_controls(),
        [0.005;2],&[0.0;2],empty,||false).is_err());
    assert_eq!(motion,original); assert_eq!(thermal.cells(),&before);
    assert_eq!(contact.state().slip_m(),[0.0;2]); assert_eq!(contact.work().dissipated_work_j(),0.0);
}
