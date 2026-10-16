use super::*;
use crate::fixture;

fn near(a: f64,b: f64,tolerance: f64) { assert!((a-b).abs() <= tolerance,"{a:e} != {b:e}"); }

#[test]
fn scheduled_births_do_not_preload_or_preage_future_material() {
    let law = fixture::material();
    let state = fixture::build(&law).unwrap();
    let (early,r) = state.advance_to(0.25,fixture::stepping(0.03125),|| Ok(())).unwrap();
    assert_eq!(early.beam().cells().len(),4);
    assert!(!early.thermal().cells().contains_key(&CellId(4)));
    assert_eq!(r.born_cells,0);
    near(r.last_mechanics.unwrap().vertical_reactions_n[0],early.beam().mass_kg()*9.81,1e-12);
    let (born,r) = early.advance_to(0.5,fixture::stepping(0.03125),|| Ok(())).unwrap();
    assert_eq!(born.beam().cells().len(),6); assert_eq!(r.born_cells,2);
    assert!(r.birth_equilibrium_pending); assert!(r.last_mechanics.is_none());
    near(r.thermal.incoming_mass_kg,2.0*fixture::source_cell(4,500.0).mass_kg,1e-16);
    for cell in &born.beam().cells()[4..] {
        let layer = &cell.section().layers()[&9];
        assert_eq!(layer.birth_time_s,0.5);
        for p in layer.material_points() {
            assert_eq!(p.temperature_k(),500.0); assert_eq!(p.stress_pa(),0.0);
            assert!(p.branch_strains().iter().all(|&q| q==0.0));
        }
    }
    let (later,r) = born.advance_to(0.75,fixture::stepping(0.03125),|| Ok(())).unwrap();
    assert!(!r.birth_equilibrium_pending);
    near(r.last_mechanics.unwrap().vertical_reactions_n[0],later.beam().mass_kg()*9.81,1e-12);
    assert!(later.beam().cells()[4].section().layers()[&9].material_points()[0].branch_strains()[0] != 0.0);
}

#[test]
fn every_actual_thermal_step_updates_the_matching_material_temperature() {
    let law = fixture::material(); let state = fixture::build(&law).unwrap();
    let (next,r) = state.advance_to(0.25,fixture::stepping(0.03125),|| Ok(())).unwrap();
    assert_eq!(r.material_point_updates,8*r.thermal.substeps);
    for cell in next.beam().cells() {
        for p in cell.section().layers()[&9].material_points() {
            assert_eq!(p.temperature_k(),next.thermal().temperature_k(CellId(cell.id())).unwrap());
        }
    }
    near(next.beam().mass_kg(),next.thermal().mass_kg(),1e-14);
    assert!(r.thermal.external_heat_j < 0.0); near(r.thermal.residual_j,0.0,1e-9);
    assert!(next.beam().deformed_position_m(next.beam().length_m(),0.0).unwrap()[1] < 0.0);
}

#[test]
fn same_grid_restart_keeps_birth_clocks_heat_and_relaxation_history() {
    let law=fixture::material(); let state=fixture::build(&law).unwrap();
    let (whole,_) = state.advance_to(2.0,fixture::stepping(0.03125),|| Ok(())).unwrap();
    let (prefix,_) = state.advance_to(0.5,fixture::stepping(0.03125),|| Ok(())).unwrap();
    let (resumed,_) = prefix.clone().advance_to(2.0,fixture::stepping(0.03125),|| Ok(())).unwrap();
    assert_eq!(whole.beam(),resumed.beam());
    assert_eq!(whole.thermal().cells(),resumed.thermal().cells());
    assert_eq!(whole.thermal().links(),resumed.thermal().links());
}

#[test]
fn late_mechanical_refusal_discards_already_staged_thermal_births_and_heat() {
    let law=fixture::material(); let (thermal,beam)=fixture::initial(&law,fixture::schedule()).unwrap();
    let mut bindings=fixture::bindings(&law); bindings[2].thermal_material=MaterialId(99);
    let state=GrowingRoad::new(thermal,beam,fixture::BODY,bindings,fixture::config(),|| Ok(())).unwrap();
    assert!(matches!(state.advance_to(2.0,fixture::stepping(0.03125),|| Ok(())),Err(Error::Invalid(_))));
    assert_eq!(state.thermal().time_s(),0.0); assert_eq!(state.beam().time_s(),0.0);
    assert_eq!(state.thermal().cells().len(),4); assert_eq!(state.beam().cells().len(),4);
    for cell in state.thermal().cells().values() { assert_eq!(cell.specific_h_j_kg,300_000.0); }
    assert!(state.beam().cells().iter().all(|c| c.section().layers()[&9].material_points()[0].stress_pa()==0.0));
}

#[test]
fn cumulative_material_budget_and_late_cancellation_do_not_renew_at_births() {
    let law=fixture::material(); let (thermal,beam)=fixture::initial(&law,fixture::schedule()).unwrap();
    let cfg=Config { max_material_point_updates: 130,..fixture::config() };
    let state=GrowingRoad::new(thermal,beam,fixture::BODY,fixture::bindings(&law),cfg,|| Ok(())).unwrap();
    // 16 steps * 8 points consumes 128; the second same-time birth exceeds 130.
    assert!(matches!(state.advance_to(0.5,fixture::stepping(0.03125),|| Ok(())),Err(Error::Budget)));
    assert_eq!(state.beam().cells().len(),4); assert_eq!(state.thermal().cells().len(),4);
    let state=fixture::build(&law).unwrap(); let mut full_polls=0;
    state.advance_to(1.0,fixture::stepping(0.03125),|| {full_polls+=1;Ok(())}).unwrap();
    let mut calls=0;
    let result=state.advance_to(1.0,fixture::stepping(0.03125),|| {
        calls+=1;
        if calls>=full_polls-1 { Err(DepositionError::Cancelled) } else { Ok(()) }
    });
    assert!(matches!(result,Err(Error::Thermal(DepositionError::Cancelled))));
    assert_eq!(state.beam().time_s(),0.0); assert_eq!(state.thermal().cells().len(),4);
}

#[test]
fn reference_mismatch_and_unbound_road_births_refuse_without_fake_coupling() {
    let law=fixture::material();
    let mut events=fixture::schedule();
    if let Event::Birth(c)=&mut events[0].event { c.geometry.start_m[1]=0.001; }
    let (thermal,beam)=fixture::initial(&law,events).unwrap();
    assert!(GrowingRoad::new(thermal,beam,fixture::BODY,fixture::bindings(&law),fixture::config(),|| Ok(())).is_err());
    let (thermal,beam)=fixture::initial(&law,fixture::schedule()).unwrap();
    let state=GrowingRoad::new(thermal,beam,fixture::BODY,vec![],fixture::config(),|| Ok(())).unwrap();
    assert!(matches!(state.advance_to(0.5,fixture::stepping(0.03125),|| Ok(())),Err(Error::Invalid(_))));
    let (thermal,beam)=fixture::initial(&law,fixture::schedule()).unwrap();
    let mut duplicate=fixture::bindings(&law); duplicate.push(duplicate[0]);
    assert!(GrowingRoad::new(thermal,beam,fixture::BODY,duplicate,fixture::config(),|| Ok(())).is_err());
}


#[test]
fn declared_drawing_tension_changes_growth_sag_without_changing_birth_heat() {
    use fs_solid::bridge::drawing::Drawing;
    let law=fixture::material();
    let slack=fixture::build(&law).unwrap();
    let drawn=fixture::build_with_drawing(&law,Drawing::TensionN(0.002)).unwrap();
    let (slack,rs)=slack.advance_to(2.0,fixture::stepping(0.03125),|| Ok(())).unwrap();
    let (drawn,rd)=drawn.advance_to(2.0,fixture::stepping(0.03125),|| Ok(())).unwrap();
    assert_eq!(slack.thermal().cells(),drawn.thermal().cells());
    assert_eq!(slack.thermal().links(),drawn.thermal().links());
    assert_eq!(rs.thermal,rd.thermal);
    assert_eq!(slack.beam().cells().len(),drawn.beam().cells().len());
    let zs=slack.beam().deformed_position_m(slack.beam().length_m(),0.0).unwrap()[1];
    let zd=drawn.beam().deformed_position_m(drawn.beam().length_m(),0.0).unwrap()[1];
    assert!(zd>zs,"drawing tension should reduce downward sag: {zd:e} <= {zs:e}");
    let report=rd.last_mechanics.unwrap();
    assert_eq!(report.drawing,Drawing::TensionN(0.002));
    near(report.drawing_force_n,0.002,1e-12);
    for cell in &drawn.beam().cells()[4..] {
        assert!(cell.section().layers()[&9].birth_time_s>=0.5);
    }
}
