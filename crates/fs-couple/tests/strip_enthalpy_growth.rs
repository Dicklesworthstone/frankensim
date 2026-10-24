//! Real coupled intervals followed by actual thermal events, not an event mock.
#[path = "strip_enthalpy/common.rs"]
mod common;
use common::*;
use fs_couple::strip_enthalpy::{self, LayerThermalBinding, StripEnthalpyError};
use fs_couple::strip_enthalpy::growth::{GrowingStrip, GrowthError, LayerBirth};
use fs_solid::deposition::{self, AxialControl, BendingControl, Control, LayeredStrip};
use fs_thermal::enthalpy::*;

fn stretched(e: f64) -> Control {
    Control { axial: AxialControl::Strain(e), bending: BendingControl::CurvaturePerM(0.0) }
}
fn birth(materials: &[Material; 2], id: u64, index: usize) -> LayerBirth<'_> {
    LayerBirth { binding: LayerThermalBinding { layer: id, cell: CellId(id),
        material: &materials[0].mechanical, caloric: &materials[0].caloric }, geometry: geometry(index) }
}
fn event(id: u64, index: usize) -> TimedEvent {
    let mut incoming = cell(CellId(id), index);
    incoming.material_id = MaterialId(0);
    TimedEvent { time_s: 0.4, event: Event::Birth(incoming) }
}

fn growth_materials(slope: f64) -> [Material; 2] {
    [5.0,10.0].map(|cp|material_with_affine_scale(cp,3e8,0.05,0.0,
        vec![(280.0,0.002),(320.0,0.002)],None,
        vec![(280.0,0.1),(320.0,0.1)],slope))
}

// Independent fixture law U=(1-300*b)*Xi, integrated over the original two
// Gauss points. No source coefficient or accepted branch history is reset.
fn internal_energy(strip: &LayeredStrip<'_>, slope: f64) -> f64 {
    strip.layers().values().map(|layer| {
        0.5*strip.length_m()*layer.geometry.area_m2*layer.material_points().iter().map(|p| {
            let e=p.total_strain()-p.birth_strain();
            let q=p.branch_strains()[0];
            (1.0-300.0*slope)*0.5*(1e8*e*e+3e8*q*q)
        }).sum::<f64>()
    }).sum()
}

#[test]
fn joint_birth_matches_original_owners_and_keeps_existing_histories() {
    joint_birth_case(0.0);
}

fn joint_birth_case(slope: f64) {
    let materials = growth_materials(slope);
    let (heat, strip) = initial(&materials,0.002,vec![event(99,2)]);
    let map = bindings(&materials);
    let (old, initialization) = GrowingStrip::try_new(heat.clone(),strip.clone(),1,
        map.to_vec(),vec![birth(&materials,99,2)],3,1e-9,||false).unwrap();
    assert_eq!(initialization.events,0);
    assert!(old.step_candidate(0.5,stretched(0.01),config(),&[0.0;2],coefficients,||false).is_err());
    let left = strip_enthalpy::step_candidate(&heat,&strip,&map,0.4,stretched(0.01),config(),
        &[0.0;2],coefficients,||false).unwrap();
    let (manual_heat, _) = left.thermal().advance_to(0.4,
        StepConfig { max_dt_s:0.1,energy_tolerance_j:1e-9 },||Ok(())).unwrap();
    let manual_strip = left.strip().deposit(99,geometry(2),materials[0].mechanical.law(),
        manual_heat.temperature_k(CellId(99)).unwrap()).unwrap();
    let (next, report) = old.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],coefficients,||false).unwrap();
    assert_eq!(next.thermal().cells(),manual_heat.cells());
    assert_eq!(next.strip(),&manual_strip);
    assert_eq!(report.events.events,1);
    assert_eq!(report.events.substeps,0);
    assert_eq!(report.born_layers,1);
    assert!(report.birth_equilibrium_pending);
    assert_eq!(next.remaining_births(),0);
    near(report.events.incoming_mass_kg,MASS,1e-16);
    for id in LAYERS { assert_eq!(next.strip().layers()[&id],left.strip().layers()[&id]); }
    for point in next.strip().layers()[&99].material_points() {
        assert_eq!(point.stress_pa(),0.0);
        assert_eq!(point.branch_strains(),&[0.0]);
        assert_eq!(point.total_strain(),point.birth_strain());
    }
    near(internal_energy(next.strip(),slope)-internal_energy(old.strip(),slope),
        report.interval.mechanical_internal_energy_change_j,1e-10);
    near(internal_energy(next.strip(),slope)-internal_energy(old.strip(),slope)
        +next.thermal().energy_j()-old.thermal().energy_j()
        +report.interval.mechanical.energy_j.numerical_dissipation
        -report.interval.external_mechanical_work_j-report.interval.non_mechanical_heat_j
        -report.events.incoming_energy_j-report.events.external_heat_j,
        0.0,config().total_energy_tolerance_j);
    if slope!=0.0 {
        assert!(report.interval.reversible_heat_j!=0.0);
        assert!(report.interval.mechanical_entropy_energy_change_j.abs()>1e-8);
    }
    // Birth has no prior physical heat or aging. The subsequent real interval does.
    assert!(!report.interval.layer_energy_j.contains_key(&99));
    let (after, r) = next.step_candidate(0.6,stretched(0.015),config(),&[0.0;3],coefficients,||false).unwrap();
    let manual = strip_enthalpy::step_candidate(next.thermal(),next.strip(),next.bindings(),0.6,
        stretched(0.015),config(),&[0.0;3],coefficients,||false).unwrap();
    assert_eq!(after.strip(),manual.strip());
    assert_eq!(after.thermal().cells(),manual.thermal().cells());
    assert!(r.interval.layer_energy_j[&99].dissipated>0.0);
    assert!(!r.birth_equilibrium_pending);
    assert_eq!(after.strip().layers()[&99].birth_time_s,0.4);
    near(internal_energy(after.strip(),slope)-internal_energy(next.strip(),slope),
        r.interval.mechanical_internal_energy_change_j,1e-10);
    assert!(report.total_energy_residual_j.abs()<=config().total_energy_tolerance_j);
    assert!(r.total_energy_residual_j.abs()<=config().total_energy_tolerance_j);
    assert_eq!(old.thermal().time_s(),0.0);
    assert_eq!(old.remaining_births(),1);
}

#[test]
fn last_bad_birth_discards_the_solved_interval_and_all_earlier_births() {
    bad_birth_case(0.0);
}

fn bad_birth_case(slope: f64) {
    let materials = growth_materials(slope);
    for wrong_body in [false,true] {
        let mut bad = event(100,3);
        if let Event::Birth(c) = &mut bad.event {
            if wrong_body { c.body_id=2; } else { c.mass_kg *= 2.0; }
        }
        let (heat,strip) = initial(&materials,0.002,vec![event(99,2),bad]);
        let (old,_) = GrowingStrip::try_new(heat,strip,1,bindings(&materials).to_vec(),
            vec![birth(&materials,99,2),birth(&materials,100,3)],4,1e-9,||false).unwrap();
        let before = old.clone();
        let calls = std::cell::Cell::new(0);
        let result = old.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],|s| {
            calls.set(calls.get()+1); coefficients(s)
        },||false);
        assert!(calls.get()>0,"the structural/thermal interval really ran before birth admission");
        assert!(matches!(result,Err(GrowthError::Events(_))));
        assert_eq!(old.thermal().cells(),before.thermal().cells());
        assert_eq!(old.strip(),before.strip());
        assert_eq!(old.remaining_births(),2);
        assert_eq!(old.thermal().time_s(),0.0);
        assert_eq!(internal_energy(old.strip(),slope),internal_energy(before.strip(),slope));
    }
}

#[test]
fn affine_hereditary_internal_energy_survives_actual_birth_and_late_refusal() {
    joint_birth_case(-0.005);
    bad_birth_case(-0.005);
}

#[test]
fn final_cancellation_and_checkpoint_replay_keep_all_state_atomic() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let (heat,strip) = initial(&materials,0.002,vec![event(99,2)]);
    let (old,_) = GrowingStrip::try_new(heat,strip,1,bindings(&materials).to_vec(),
        vec![birth(&materials,99,2)],3,1e-9,||false).unwrap();
    let saved = old.clone();
    let mut polls = 0;
    let (next,_) = old.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],coefficients,
        ||{polls+=1;false}).unwrap();
    let mut actual = 0;
    assert!(matches!(old.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],coefficients,
        ||{actual+=1;actual==polls}),Err(GrowthError::Events(StripEnthalpyError::Cancelled))));
    assert_eq!(old.strip(),saved.strip());
    assert_eq!(old.thermal().cells(),saved.thermal().cells());
    let (replay,_) = saved.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],coefficients,||false).unwrap();
    assert_eq!(next.strip(),replay.strip());
    assert_eq!(next.thermal().cells(),replay.thermal().cells());
    assert_eq!(next.remaining_births(),replay.remaining_births());
}

#[test]
fn scheduled_initial_births_use_real_enthalpy_and_create_no_prebirth_history() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let events = CELLS.into_iter().enumerate().map(|(i,id)|TimedEvent {
        time_s:0.0,event:Event::Birth(cell(id,i)) }).collect();
    let thermal = DepositionSimulation::new(events,&[
        (MaterialId(0),materials[0].caloric.curve() as &dyn ThermalLaw),
        (MaterialId(1),materials[1].caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:1000}).unwrap();
    let strip = LayeredStrip::new(LENGTH,deposition::Budget {max_layers:8,max_branch_states:32}).unwrap();
    let map = bindings(&materials);
    let births = [0,1].map(|i|LayerBirth {binding:map[i],geometry:geometry(i)}).to_vec();
    let (state,audit) = GrowingStrip::try_new(thermal,strip,1,vec![],births,2,1e-9,||false).unwrap();
    near(audit.incoming_mass_kg,2.0*MASS,1e-16);
    assert_eq!(audit.events,2);
    assert_eq!(audit.substeps,0);
    assert_eq!(state.remaining_births(),0);
    assert_eq!(state.strip().stored_energy_j(),0.0);
    assert!(state.birth_equilibrium_pending());
    let (_,r) = state.step_candidate(0.2,stretched(0.01),config(),&[0.0;2],coefficients,||false).unwrap();
    assert!(r.interval.deposited_heat_j>0.0);
}

#[test]
fn passive_body_births_retain_enthalpy_without_invented_mechanical_layers() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let mut incoming = event(99,2);
    if let Event::Birth(c) = &mut incoming.event { c.body_id=2; c.specific_h_j_kg=50.0; }
    let (heat,strip) = initial(&materials,0.002,vec![incoming]);
    let (old,_) = GrowingStrip::try_new(heat,strip,1,bindings(&materials).to_vec(),vec![],2,1e-9,||false).unwrap();
    let (next,r) = old.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],coefficients,||false).unwrap();
    assert_eq!(next.thermal().cells().len(),3);
    assert_eq!(next.strip().layers().len(),2);
    assert_eq!(r.born_layers,0);
    assert!(!r.birth_equilibrium_pending);
    near(r.events.incoming_energy_j,MASS*50.0,1e-14);
    near(next.thermal().temperature_k(CellId(99)).unwrap(),310.0,1e-10);
    assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
}

#[test]
fn registry_and_selected_body_completeness_are_checked_before_advancement() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let (heat,strip) = initial(&materials,0.002,vec![event(99,2)]);
    let map = bindings(&materials);
    assert!(GrowingStrip::try_new(heat.clone(),strip.clone(),1,map.to_vec(),
        vec![birth(&materials,99,2)],2,1e-9,||false).is_err());
    assert!(GrowingStrip::try_new(heat.clone(),strip.clone(),1,map.to_vec(),
        vec![birth(&materials,99,2),birth(&materials,99,3)],4,1e-9,||false).is_err());
    assert!(GrowingStrip::try_new(heat.clone(),strip.clone(),2,map.to_vec(),vec![],2,1e-9,||false).is_err());
    assert!(GrowingStrip::try_new(heat.clone(),strip.clone(),1,map[..1].to_vec(),vec![],2,1e-9,||false).is_err());
    let (unregistered,_) = GrowingStrip::try_new(heat,strip,1,map.to_vec(),vec![],2,1e-9,||false).unwrap();
    assert!(matches!(unregistered.step_candidate(0.4,stretched(0.01),config(),&[0.0;2],coefficients,||false),
        Err(GrowthError::Events(StripEnthalpyError::Invalid(_)))));
}
