//! Actual thermal histories, material substitutions and existing dynamics.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_ga::{Se3,Vec3,Wrench};
use fs_material::phase::SolidLiquidPhase;
use fs_mbd::articulated::{ArticulatedModel,BaseState,JointModel,Link,SpatialInertia,forward_dynamics};
use fs_mbd_material::thermal::heating::{HeatedBody,HeatingError};
use fs_thermal::enthalpy::*;
use fs_thermal::observers::ObserverCouplingError;
use std::collections::BTreeMap;

fn acceleration(inertia:SpatialInertia)->f64 {
    let model=ArticulatedModel::new(vec![
        Link::new("base",None,Se3::identity(),JointModel::FIXED,inertia),
        Link::new("rotor",Some(0),Se3::identity(),JointModel::revolute(Vec3::new(0.,0.,1.),None).unwrap(),inertia),
    ]).unwrap();
    forward_dynamics(&model,BaseState::stationary(Se3::identity()),&[0.],&[0.],&[1.],
        Vec3::default(),&[Wrench::default();2]).unwrap().generalized_acceleration[0]
}

#[test]
fn cp_only_material_change_propagates_heat_to_geometry_and_actual_rotor_response() {
    let a=card(10.,1e-3,500.); let b=card(20.,1e-3,500.);
    let ra=reference(&a,1); let rb=reference(&b,1);
    let ca=caloric(&a,1e5,SolidLiquidPhase::Solid); let cb=caloric(&b,1e5,SolidLiquidPhase::Solid);
    let off=vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),boundary(0.))}];
    let old_a=heated(&ra,&a,&ca,40.,off.clone()); let old_b=heated(&rb,&b,&cb,40.,off);
    let (hot_a,aa)=old_a.advance_to(2.,config(),||Ok(())).unwrap();
    let (hot_b,ab)=old_b.advance_to(2.,config(),||Ok(())).unwrap();
    assert_eq!(ra.inertia(),rb.inertia());
    close(hot_a.thermal().temperature_k(CellId(0)).unwrap(),302.);
    close(hot_b.thermal().temperature_k(CellId(0)).unwrap(),301.);
    close(aa.external_heat_j,40.); close(ab.external_heat_j,40.);
    close(hot_a.configuration().regions()[0].geometry().volume_m3(),1.002_f64.powi(3));
    close(acceleration(hot_a.configuration().inertia())/acceleration(hot_b.configuration().inertia()),
        (1.001_f64/1.002).powi(2));
    for body in [&hot_a,&hot_b] {
        assert_eq!(body.thermal().mass_kg().to_bits(),ra.inertia().mass().to_bits());
        assert_eq!(body.configuration().inertia().mass().to_bits(),ra.inertia().mass().to_bits());
        let (later,audit)=body.advance_to(3.,config(),||Ok(())).unwrap();
        assert_eq!(later.configuration(),body.configuration());
        close(audit.external_heat_j,0.);
    }
    close(old_a.thermal().time_s(),0.); close(old_a.thermal().energy_j(),0.);
}

#[test]
fn shared_contact_cfl_and_cooling_match_the_original_thermal_path() {
    let card=card(10.,1e-3,500.); let reference=reference(&card,2);
    let chart=caloric(&card,1e5,SolidLiquidPhase::Solid);
    let link=ThermalLink {id:77,a:CellId(0),b:CellId(1),area_m2:1.,conductance_w_k:20.};
    let events=vec![TimedEvent {time_s:0.,event:Event::Contact(link)}];
    let body=HeatedBody::new(&reference,vec![region(&card,&chart,"b",1,0.),region(&card,&chart,"a",0,40.)],
        events,budget(),config(),||Ok(())).unwrap();
    let cfg=StepConfig {max_dt_s:10.,energy_tolerance_j:1e-9};
    let (plain,expected)=body.thermal().advance_to(1.3,cfg,||Ok(())).unwrap();
    let (coupled,actual)=body.advance_to(1.3,cfg,||Ok(())).unwrap();
    assert!(actual.substeps>1); assert_eq!(actual,expected);
    assert_eq!(plain.cells(),coupled.thermal().cells());
    close(coupled.thermal().energy_j(),52.);
    assert!(coupled.thermal().temperature_k(CellId(1)).unwrap()>300.);
    for (i,region) in coupled.configuration().regions().iter().enumerate() {
        let t=coupled.thermal().temperature_k(CellId(i as u64)).unwrap();
        close(region.expansion().free_linear_strain(),1e-3*(t-300.));
        assert_eq!(region.inertia().mass(),reference.regions()[i].inertia().mass());
    }
    let single_reference=common::reference(&card,1);
    let single=heated(&single_reference,&card,&chart,40.,cooling_events());
    let (hot,_)=single.advance_to(1.,config(),||Ok(())).unwrap();
    let (cool,audit)=hot.advance_to(3.,config(),||Ok(())).unwrap();
    assert!(cool.thermal().temperature_k(CellId(0)).unwrap()<hot.thermal().temperature_k(CellId(0)).unwrap());
    assert!(cool.configuration().regions()[0].geometry().volume_m3()<hot.configuration().regions()[0].geometry().volume_m3());
    assert!(audit.external_heat_j<0.);
}

// Cooling helper takes a reference owned by its caller in the test below.
fn cooling_events()->Vec<TimedEvent> {
    let mut cool=boundary(0.); cool.htc_w_m2_k=1.;
    vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),cool)}]
}

#[derive(Clone,Debug,PartialEq)]
struct GeneratedHeat { cell:CellId,watts:f64,calls:usize,work:f64,fail:bool }
impl DepositionObserver for GeneratedHeat {
    type Error=DepositionError;
    fn sources(&mut self,_:&DepositionSimulation<'_>,check:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {
        check()?;self.calls+=1;Ok(BTreeMap::from([(self.cell,self.watts)]))
    }
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error>{
        self.work+=dt*self.watts;
        if self.fail {return Err(DepositionError::InvalidInput("source refuses after a physical step"));} Ok(())
    }
}

#[test]
fn generated_heat_is_counted_once_and_both_histories_roll_back_together() {
    let card=card(10.,1e-3,500.);let reference=reference(&card,1);
    let chart=caloric(&card,1e5,SolidLiquidPhase::Solid);
    let body=heated(&reference,&card,&chart,5.,vec![]);
    let source=GeneratedHeat {cell:CellId(0),watts:35.,calls:0,work:0.,fail:false};
    let (next,accepted,audit)=body.advance_with(&source,1.3,config(),||Ok(())).unwrap();
    close(accepted.work,35.*1.3);close(audit.external_heat_j,40.*1.3);
    close(next.thermal().temperature_k(CellId(0)).unwrap(),302.6);
    assert_eq!(accepted.calls,audit.substeps);assert_eq!(source.calls,0);
    let mut bad=source.clone();bad.fail=true;
    assert!(matches!(body.advance_with(&bad,1.,config(),||Ok(())),Err(ObserverCouplingError::First(_))));
    assert_eq!(body.thermal().time_s(),0.);assert_eq!(bad.calls,0);assert_eq!(bad.work,0.);
    let mut checkpoints=0;
    body.advance_with(&source,1.,config(),||{checkpoints+=1;Ok(())}).unwrap();
    let mut n=0;
    let failed=body.advance_with(&source,1.,config(),||{n+=1;if n==checkpoints {Err(DepositionError::Cancelled)}else{Ok(())}});
    assert!(matches!(failed,Err(ObserverCouplingError::Thermal(DepositionError::Cancelled))));
    assert_eq!(body.thermal().energy_j(),0.);assert_eq!(source.calls,0);
}

#[test]
fn invalid_intermediate_expansion_cannot_hide_behind_a_cool_final_state() {
    let card=card(10.,1e-3,301.);let reference=reference(&card,1);
    let chart=caloric(&card,1e5,SolidLiquidPhase::Solid);
    let events=vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),boundary(-40.))}];
    let body=heated(&reference,&card,&chart,40.,events);
    // Endpoint-only coupling would accept this trajectory at 300 K, although
    // it leaves the expansion source domain during heating.
    let (plain,_)=body.thermal().advance_to(2.,config(),||Ok(())).unwrap();
    close(plain.temperature_k(CellId(0)).unwrap(),300.);
    assert!(matches!(body.advance_to(2.,config(),||Ok(())),Err(HeatingError::Mechanical(_))));
    assert_eq!(body.thermal().time_s(),0.);assert_eq!(body.thermal().energy_j(),0.);
    assert_eq!(body.configuration().inertia(),reference.inertia());
}

#[test]
fn source_context_initial_state_and_inventory_mismatches_are_refused() {
    let card=card(10.,1e-3,500.);let reference=reference(&card,1);
    let chart=caloric(&card,1e5,SolidLiquidPhase::Solid);
    let wrong_context=caloric(&card,2e5,SolidLiquidPhase::Solid);
    let liquid=caloric(&card,1e5,SolidLiquidPhase::Liquid);
    for curve in [&wrong_context,&liquid] {
        assert!(HeatedBody::new(&reference,vec![region(&card,curve,"a",0,0.)],vec![],budget(),config(),||Ok(())).is_err());
    }
    let mut bad=region(&card,&chart,"a",0,0.);bad.geometry.volume_m3=2.;
    assert!(HeatedBody::new(&reference,vec![bad],vec![],budget(),config(),||Ok(())).is_err());
    let mut bad=region(&card,&chart,"a",0,0.);bad.initial_specific_enthalpy_j_kg=10.;
    assert!(HeatedBody::new(&reference,vec![bad],vec![],budget(),config(),||Ok(())).is_err());
    let too_small=DepositionBudget {max_events:0,..budget()};
    assert!(HeatedBody::new(&reference,vec![region(&card,&chart,"a",0,0.)],vec![],too_small,config(),||Ok(())).is_err());
}

#[test]
fn aligned_resume_preserves_enthalpy_configuration_and_source_off_history() {
    let card=card(10.,1e-3,500.);let reference=reference(&card,1);
    let chart=caloric(&card,1e5,SolidLiquidPhase::Solid);
    let body=heated(&reference,&card,&chart,40.,cooling_events());
    let (whole,a)=body.advance_to(3.,config(),||Ok(())).unwrap();
    let (half,b)=body.advance_to(1.,config(),||Ok(())).unwrap();
    let (rest,c)=half.advance_to(3.,config(),||Ok(())).unwrap();
    assert_eq!(whole.thermal().cells(),rest.thermal().cells());
    assert_eq!(whole.configuration(),rest.configuration());
    close(a.external_heat_j,b.external_heat_j+c.external_heat_j);
    let limited=DepositionBudget {max_substeps:1,..budget()};
    let small=HeatedBody::new(&reference,vec![region(&card,&chart,"a",0,40.)],vec![],limited,config(),||Ok(())).unwrap();
    assert!(matches!(small.advance_to(3.,config(),||Ok(())),Err(HeatingError::Thermal(DepositionError::Budget(_)))));
    assert_eq!(small.thermal().time_s(),0.);
}
