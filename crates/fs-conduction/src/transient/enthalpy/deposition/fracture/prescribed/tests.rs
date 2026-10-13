use super::*;
use fs_material::cohesive::BilinearCohesive;
use crate::transient::deposition::{Boundary, CellGeometry, DepositedCell,
    DepositionBudget, MaterialId, ThermalLaw, TimedEvent};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static HEAT: Sensible = Sensible;
fn law() -> BilinearCohesive { BilinearCohesive::new(1e9,2e5,100.0,2e9,0.01).unwrap() }
fn config() -> Config {
    Config { thermal: StepConfig { max_dt_s: 0.001,energy_tolerance_j: 1e-7 },
        event_time_tolerance_s: 1e-9,local_energy_tolerance_j_m2: 1e-8,
        max_event_trials: 1024,max_crossings: 16,max_thermal_substeps: 20_000 }
}
fn endpoint(id: u64,jump_m: [f64; 3]) -> Endpoint { Endpoint { interface: id,jump_m } }
fn binding(point: CohesivePoint<'_>) -> Interface<'_> {
    Interface { reference: ThermalLink { id: 7,a: CellId(1),b: CellId(2),area_m2: 0.0002,
        conductance_w_k: 0.2 },point,closed: true,
        policy: ContactPolicy { open_gap_m: 0.0001,reclose_gap_m: 0.0,
            reclosed_conductance_w_k: 0.05,maximum_slip_m: 0.005 } }
}
fn coupled<'a>(interfaces: Vec<Interface<'a>>,extra: Vec<TimedEvent>) -> FractureCooling<'static,'a> {
    let mut events = Vec::new();
    for (id,t) in [(1,500.0),(2,300.0)] {
        events.push(TimedEvent { time_s: 0.0,event: Event::Birth(DepositedCell {
            id: CellId(id),segment_id: id,body_id: id,material_id: MaterialId(1),mass_kg: 0.001,
            specific_h_j_kg: 2000.0*(t-300.0),
            geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],end_m: [0.01,0.0,id as f64*0.001],
                volume_m3: 1e-6,surface_area_m2: 0.001 },
            boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 0.0,reservoir: None,source_w: 0.0 },
        }) });
    }
    for i in &interfaces {
        if let Some(link) = i.thermal_link() { events.push(TimedEvent { time_s: 0.0,event: Event::Contact(link) }); }
    }
    events.extend(extra);
    let thermal = DepositionSimulation::new(events,&[(MaterialId(1),&HEAT)],DepositionBudget {
        max_cells: 4,max_links: 8,max_events: 32,max_substeps: 20_000,
    }).unwrap().advance_to(0.0,config().thermal,|| Ok(())).unwrap().0;
    FractureCooling::new(thermal,interfaces,8,|| Ok(())).unwrap()
}
fn near(a: f64,b: f64,tolerance: f64) { assert!((a-b).abs() < tolerance,"{a:e} != {b:e}"); }
fn failure_time(r: &Report) -> f64 {
    r.crossings.iter().find(|e| e.kind == CrossingKind::Failure).unwrap().after_s
}

// G1: real thermal integration must stop exchanging heat at t=5, not at t=10.
#[test]
fn failure_inside_motion_changes_heat_before_the_final_endpoint() {
    let law = law();
    let initial = coupled(vec![binding(law.bonded_point())],vec![]);
    let (result,r) = initial.advance_prescribed_to(10.0,&[endpoint(7,[0.002,0.0,0.0])],config(),|| Ok(())).unwrap();
    near(failure_time(&r),5.0,1e-9);
    near(result.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-1.0_f64).exp(),0.005);
    assert_eq!(r.contacts.opened,1);
    assert!(result.thermal().links().is_empty());
    near(r.contacts.fracture_dissipation_j,0.02,1e-12);
    near(r.thermal.residual_j,0.0,1e-7);
    assert!(result.interfaces()[&7].point.failed());
    for event in &r.crossings { assert!(event.after_s-event.before_s <= config().event_time_tolerance_s); }
    let (lagged,_) = initial.advance_to(10.0,config().thermal,|| Ok(())).unwrap();
    assert!(result.thermal().temperature_k(CellId(1)).unwrap()-lagged.thermal().temperature_k(CellId(1)).unwrap() > 20.0);
}

#[test]
fn shear_failure_changes_conductance_without_inventing_open_surface() {
    let law = law();
    let initial = coupled(vec![binding(law.bonded_point())],vec![]);
    let (result,r) = initial.advance_prescribed_to(10.0,&[endpoint(7,[0.0,0.002,0.0])],config(),|| Ok(())).unwrap();
    near(failure_time(&r),5.0,1e-9);
    near(result.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-1.25_f64).exp(),0.005);
    assert_eq!(r.contacts.changed_conductances,1);
    assert_eq!(r.contacts.opened,0);
    assert_eq!(result.thermal().links()[&7].conductance_w_k,0.05);
    near(result.thermal().exposed_area_m2(CellId(1)).unwrap(),0.0008,1e-15);
}

// G0: final closed state cannot reveal this intervening open interval.
#[test]
fn one_nonproportional_segment_can_fail_open_and_reclose() {
    let law = law();
    let point = law.bonded_point().trial([0.0008,0.0,0.0],1e-8).unwrap().point;
    let initial = coupled(vec![binding(point)],vec![]);
    let (result,r) = initial.advance_prescribed_to(10.0,&[endpoint(7,[-0.0002,0.0018,0.0])],config(),|| Ok(())).unwrap();
    // Solve (0.8-f)^2+(1.8f)^2=1 in millimetre coordinates.
    let expected_failure = 10.0*(1.6+(1.6_f64*1.6+4.0*4.24*0.36).sqrt())/(2.0*4.24);
    near(failure_time(&r),expected_failure,1e-9);
    let close = r.crossings.iter().find(|e| e.kind == CrossingKind::ReclosingGap).unwrap();
    near(close.after_s,8.0,1e-9);
    let expected_hot = 400.0+100.0*(-0.2*expected_failure-0.05*2.0).exp();
    near(result.thermal().temperature_k(CellId(1)).unwrap(),expected_hot,0.005);
    assert_eq!(r.contacts.opened,1); assert_eq!(r.contacts.reclosed,1);
    assert!(result.interfaces()[&7].closed && result.interfaces()[&7].point.failed());
    assert_eq!(result.interfaces()[&7].point.traction_pa()[1],0.0);
    assert_eq!(result.thermal().links()[&7].conductance_w_k,0.05);
}

#[test]
fn reclosure_is_localized_without_healing_a_failed_bond() {
    let law = law();
    let failed = law.bonded_point().trial([0.002,0.0,0.0],1e-8).unwrap().point;
    let initial = coupled(vec![Interface { closed: false,..binding(failed) }],vec![]);
    let (result,r) = initial.advance_prescribed_to(4.0,&[endpoint(7,[-0.002,0.0,0.0])],config(),|| Ok(())).unwrap();
    assert_eq!(r.contacts.reclosed,1); assert_eq!(r.contacts.fracture_dissipation_j,0.0);
    assert!(!r.crossings.iter().any(|e| e.kind == CrossingKind::Failure));
    near(result.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-0.1_f64).exp(),0.002);
    assert_eq!(result.interfaces()[&7].point.damage(),1.0);
}

// G3: same physical area and simultaneous motion, independent of input order.
#[test]
fn simultaneous_partitioned_contacts_and_endpoint_permutations_agree() {
    let law = law();
    let whole = coupled(vec![binding(law.bonded_point())],vec![]);
    let mut a = binding(law.bonded_point());
    a.reference.area_m2 *= 0.5; a.reference.conductance_w_k *= 0.5;
    a.policy.reclosed_conductance_w_k *= 0.5;
    let b = Interface { reference: ThermalLink { id: 8,..a.reference },..a };
    let split = coupled(vec![a,b],vec![]);
    let ends = [endpoint(8,[0.002,0.0,0.0]),endpoint(7,[0.002,0.0,0.0])];
    let (one,or) = whole.advance_prescribed_to(2.0,&ends[1..],config(),|| Ok(())).unwrap();
    let (two,tr) = split.advance_prescribed_to(2.0,&ends,config(),|| Ok(())).unwrap();
    let (permuted,pr) = split.advance_prescribed_to(2.0,&[ends[1],ends[0]],config(),|| Ok(())).unwrap();
    assert_eq!(two.thermal().cells(),permuted.thermal().cells());
    assert_eq!(tr.crossings,pr.crossings);
    assert_eq!(or.contacts.fracture_dissipation_j,tr.contacts.fracture_dissipation_j);
    near(one.thermal().temperature_k(CellId(1)).unwrap(),two.thermal().temperature_k(CellId(1)).unwrap(),1e-10);
}

#[test]
fn shared_work_caps_and_invalid_paths_do_not_publish_any_prefix() {
    let law = law();
    let initial = coupled(vec![binding(law.bonded_point())],vec![]);
    let end = [endpoint(7,[0.002,0.0,0.0])];
    for cfg in [Config { max_event_trials: 0,..config() },Config { max_crossings: 0,..config() },
        Config { max_thermal_substeps: 600,..config() }] {
        assert!(initial.advance_prescribed_to(10.0,&end,cfg,|| Ok(())).is_err());
    }
    for ends in [vec![],vec![endpoint(8,[0.0;3])],vec![endpoint(7,[f64::NAN,0.0,0.0])],
        vec![endpoint(7,[0.0,0.006,0.0])]] {
        assert!(initial.advance_prescribed_to(1.0,&ends,config(),|| Ok(())).is_err());
    }
    assert!(initial.advance_prescribed_to(0.0,&end,config(),|| Ok(())).is_err());
    for allowed in [0,20,200,5000] {
        let mut calls = 0;
        assert!(initial.advance_prescribed_to(10.0,&end,config(),|| {
            calls += 1;
            if calls > allowed { Err(DepositionError::Cancelled) } else { Ok(()) }
        }).is_err());
    }
    assert_eq!(initial.thermal().time_s(),0.0);
    assert_eq!(initial.thermal().temperature_k(CellId(1)).unwrap(),500.0);
    assert_eq!(initial.interfaces()[&7].point.damage(),0.0);
    assert_eq!(initial.thermal().links()[&7].conductance_w_k,0.2);
}

#[test]
fn unrelated_scheduled_heat_is_counted_and_owned_contact_conflicts_refuse() {
    let law = law();
    let extra = vec![TimedEvent { time_s: 2.0,event: Event::Boundary(CellId(1),Boundary {
        ambient_k: 300.0,htc_w_m2_k: 0.0,reservoir: None,source_w: 1.0 }) }];
    let initial = coupled(vec![binding(law.bonded_point())],extra);
    let (result,r) = initial.advance_prescribed_to(10.0,&[endpoint(7,[0.002,0.0,0.0])],config(),|| Ok(())).unwrap();
    near(r.thermal.external_heat_j,8.0,1e-9);
    near(result.thermal().energy_j(),408.0,1e-7);
    assert_eq!(r.thermal.events,1);
    let conflict = coupled(vec![binding(law.bonded_point())],vec![TimedEvent {
        time_s: 5.0,event: Event::Separation(7) }]);
    assert!(conflict.advance_prescribed_to(10.0,&[endpoint(7,[0.002,0.0,0.0])],config(),|| Ok(())).is_err());
    assert_eq!(conflict.thermal().time_s(),0.0);
}

#[test]
fn motion_without_crossings_preserves_the_ordinary_thermal_path_bitwise() {
    let law = law();
    let initial = coupled(vec![binding(law.bonded_point())],vec![]);
    let (plain,a) = initial.advance_to(1.0,config().thermal,|| Ok(())).unwrap();
    let (path,r) = initial.advance_prescribed_to(1.0,&[endpoint(7,[0.00001,0.0,0.0])],config(),|| Ok(())).unwrap();
    assert!(r.crossings.is_empty());
    assert_eq!(plain.thermal().cells(),path.thermal().cells());
    assert_eq!(a,r.thermal);
    assert_eq!(r.thermal_intervals,1);
}
