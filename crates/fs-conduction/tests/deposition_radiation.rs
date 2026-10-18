//! Actual shared thermal evolution under source-resolved enclosure radiation.
mod support;
#[path="radiation_material/common.rs"] mod common;
use common::*;
use fs_conduction::transient::deposition::*;
use fs_conduction::transient::deposition::radiation::*;
use fs_conduction::STEFAN_BOLTZMANN_W_M2_K4;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_solver::NewtonKrylovConfig;
use fs_matdb::{MaterialCard, QueryPoint};
use std::collections::BTreeMap;

fn id(i:usize)->CellId { CellId(if i==0 {30} else {10}) }
fn initial<'a>(curves:&'a [EquilibriumEnthalpyPhaseCurve;2],h:[f64;2],power:f64)->DepositionSimulation<'a> {
    let events=[0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:id(i),segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),
        mass_kg:1.0,specific_h_j_kg:h[i],geometry:CellGeometry {
            start_m:[i as f64,0.0,0.0],end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:1.5},
        boundary:Boundary {ambient_k:350.0,htc_w_m2_k:0.0,reservoir:None,source_w:if i==0 {power} else {0.0}},
    })}).to_vec();
    DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),
        (MaterialId(1),&curves[1] as &dyn ThermalLaw)],DepositionBudget {
        max_cells:2,max_links:0,max_events:2,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-8},||Ok(())).unwrap().0
}
fn bindings<'a>(cards:&'a [MaterialCard;2],curves:&'a [EquilibriumEnthalpyPhaseCurve;2],context:&'a QueryPoint)
    ->[RadiatingCell<'a>;2] {
    // Deliberately reversed input order; matrix order is not CellId order either.
    [1,0].map(|i|RadiatingCell {surface:i,cell:id(i),material_id:MaterialId(i as u64),curve:&curves[i],
        material:binding(&cards[i],context),selection:MaterialPropertySelection::SingleClaimOnly})
}
fn controls()->RadiationStepSolve {
    RadiationStepSolve {maximum_cells:2,maximum_evaluations:300,maximum_iterations:32,difference_step_j:1e-4,
        cell_tolerance_j:1e-7,total_tolerance_j:1e-7,newton:NewtonKrylovConfig {
            absolute_tolerance:1e-9,relative_tolerance:0.0,forcing_minimum:1e-4,linear_restart:4,
            ..NewtonKrylovConfig::default()}}
}
#[test]
fn actual_implicit_cell_heat_matches_independent_nonlinear_parallel_plate_root() {
    let cards=[card(0.4,0.4),card(0.6,0.6)]; let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
    let old=initial(&curves,[150.0,-50.0],0.0);
    let (next,report)=support::with_cx(|cx|model.step_candidate(cx,&old,0.01,controls(),&[0.0;2]).unwrap());
    let (mut lo,mut hi)=(0.0,50.0);
    for _ in 0..60 {
        let q=0.5*(lo+hi);
        let radiation=STEFAN_BOLTZMANN_W_M2_K4*((450.0_f64-q).powi(4)-(350.0_f64+q).powi(4))/(1.0/0.4+1.0/0.6-1.0);
        if q-0.01*radiation<0.0 {lo=q;} else {hi=q;}
    }
    let q=0.5*(lo+hi);
    near(next.cells()[&id(0)].specific_h_j_kg,150.0-q);
    near(next.cells()[&id(1)].specific_h_j_kg,-50.0+q);
    near(next.energy_j(),old.energy_j());near(report.external_heat_j,0.0);
    assert!(report.maximum_cell_residual_j<=controls().cell_tolerance_j);
    assert_eq!(old.time_s(),0.0);assert_eq!(old.cells()[&id(0)].specific_h_j_kg,150.0);
}
#[test]
fn emissivity_data_changes_phase_evolution_at_fixed_kelvin_and_preserves_independent_heat() {
    let run=|scale:f64| {
        let cards=[card(0.2*scale,0.8*scale),card(0.6,0.6)];
        let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
        let context=context();let surfaces=bindings(&cards,&curves,&context);
        let enclosure=enclosure([emissivity(&cards[0],&curves[0],50.0),emissivity(&cards[1],&curves[1],-50.0)]);
        let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
        let old=initial(&curves,[50.0,-50.0],5.0);let mut state=old.clone();
        support::with_cx(|cx| {
            for i in 1..=20 {
                let (next,r)=model.step_candidate(cx,&state,i as f64*0.002,controls(),&[0.0;2]).unwrap();
                near(r.external_heat_j,0.002*5.0);
                let hot=&r.radiation.surfaces()[0];
                assert_eq!(hot.phase.temperature_k(),400.0);
                near(r.radiation.enclosure().surfaces()[0].emissivity().value(),scale*(0.2+0.6*hot.phase.liquid_mass_fraction()));
                near(next.energy_j()-old.energy_j(),i as f64*0.002*5.0);
                state=next;
            }
        });
        (curves[0].state_at_specific_enthalpy(state.cells()[&id(0)].specific_h_j_kg).unwrap().liquid_mass_fraction(),
            state.temperature_k(id(1)).unwrap())
    };
    let low=run(0.5);let high=run(1.0);
    assert!(high.0<low.0 && high.1>low.1);
}
#[test]
fn complete_binding_summed_coverage_and_transactional_source_addition_are_enforced() {
    let cards=[card(0.4,0.4),card(0.6,0.6)];let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let mut surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let old=initial(&curves,[150.0,-50.0],0.0);
    fn make<'a>(enclosure: &'a fs_conduction::GrayDiffuseEnclosure, cells: &'a [RadiatingCell<'a>]) -> MaterialRadiation<'a> {
        MaterialRadiation {enclosure,cells,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8}
    }
    support::with_cx(|cx| {
        let source=make(&enclosure,&surfaces).resolve(cx,&old).unwrap();
        let mut map=BTreeMap::from([(id(0),5.0),(id(1),-2.0)]);
        source.add_to(&mut map).unwrap();near(map.values().copied().sum(),3.0);
        let mut bad=BTreeMap::from([(id(0),f64::INFINITY),(id(1),7.0)]);let saved=bad.clone();
        assert!(source.add_to(&mut bad).is_err());assert_eq!(bad,saved);
        assert!(make(&enclosure,&surfaces[..1]).resolve(cx,&old).is_err());
    });
    // Each trace is individually 1 m2; jointly 2 m2 exceeds this cell's 1.5 m2.
    surfaces[0].cell=id(0);surfaces[0].material_id=MaterialId(0);surfaces[0].curve=&curves[0];
    surfaces[0].material=binding(&cards[0],&context);
    support::with_cx(|cx|assert!(make(&enclosure,&surfaces).resolve(cx,&old).is_err()));
}
#[test]
fn late_work_budget_and_cancellation_leave_the_original_timeline_unchanged() {
    let cards=[card(0.4,0.4),card(0.6,0.6)];let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
    let old=initial(&curves,[150.0,-50.0],0.0);let saved=old.cells().clone();
    support::with_cx(|cx| {
        let (_,report)=model.step_candidate(cx,&old,0.01,controls(),&[0.0;2]).unwrap();
        assert!(matches!(model.step_candidate(cx,&old,0.01,
            RadiationStepSolve {maximum_evaluations:report.evaluations-1,..controls()},&[0.0;2]),
            Err(RadiationCellError::EvaluationBudget)));
    });
    support::with_cancelled_cx(|cx|assert!(matches!(model.step_candidate(cx,&old,0.01,controls(),&[0.0;2]),
        Err(RadiationCellError::Cancelled))));
    assert_eq!(old.cells(),&saved);assert_eq!(old.time_s(),0.0);
}
