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
fn policy() -> ContactPolicy {
    ContactPolicy { open_gap_m: 0.0001, reclose_gap_m: 0.0,
        reclosed_conductance_w_k: 0.05, maximum_slip_m: 0.005 }
}
fn link(id: u64) -> ThermalLink {
    ThermalLink { id, a: CellId(1), b: CellId(2), area_m2: 0.0002, conductance_w_k: 0.2 }
}
fn birth(id: u64, temperature: f64, htc: f64) -> TimedEvent {
    TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CellId(id), segment_id: id, body_id: id, material_id: MaterialId(1),
        mass_kg: 0.001, specific_h_j_kg: 2000.0*(temperature-300.0),
        geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],
            end_m: [0.01,0.0,id as f64*0.001], volume_m3: 1e-6, surface_area_m2: 0.001 },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: htc, reservoir: None, source_w: 0.0 },
    }) }
}
fn thermal(links: &[ThermalLink], htc: f64, extra: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    let mut events = vec![birth(1,500.0,htc),birth(2,300.0,htc)];
    events.extend(links.iter().map(|l| TimedEvent { time_s: 0.0,event: Event::Contact(*l) }));
    events.extend(extra);
    DepositionSimulation::new(events,&[(MaterialId(1),&HEAT)],DepositionBudget {
        max_cells: 2,max_links: 4,max_events: 16,max_substeps: 10000,
    }).unwrap().advance_to(0.0,config(),|| Ok(())).unwrap().0
}
fn config() -> StepConfig { StepConfig { max_dt_s: 0.01,energy_tolerance_j: 1e-8 } }
fn binding(law: &BilinearCohesive, reference: ThermalLink) -> Interface<'_> {
    Interface { reference,point: law.bonded_point(),closed: true,policy: policy() }
}
fn coupled(law: &BilinearCohesive, htc: f64) -> FractureCooling<'static,'_> {
    FractureCooling::new(thermal(&[link(7)],htc,vec![]),vec![binding(law,link(7))],1,|| Ok(())).unwrap()
}
fn near(a: f64,b: f64,tol: f64) { assert!((a-b).abs() < tol,"{a:e} != {b:e}"); }

// G1: without ambient loss, contact exchange has an independent two-body limit.
#[test]
fn failure_opens_actual_contact_and_changes_subsequent_temperatures() {
    let law = law();
    let state = coupled(&law,0.0);
    let (held,audit) = state.advance_to(1.0,config(),|| Ok(())).unwrap();
    near(held.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-0.2_f64).exp(),0.02);
    near(audit.residual_j,0.0,1e-8);
    let broken = law.bonded_point().trial([law.failure_separation_m(),0.0,0.0],1e-8).unwrap().point;
    let (open,report) = state.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    assert_eq!(report.opened,1);
    near(report.fracture_dissipation_j,0.02,1e-12);
    assert!(open.thermal().links().is_empty());
    assert_eq!(open.thermal().cells(),state.thermal().cells());
    near(open.thermal().exposed_area_m2(CellId(1)).unwrap(),0.001,1e-15);
    let (cooled,audit) = open.advance_to(1.0,config(),|| Ok(())).unwrap();
    assert_eq!(cooled.thermal().temperature_k(CellId(1)).unwrap(),500.0);
    assert_eq!(cooled.thermal().temperature_k(CellId(2)).unwrap(),300.0);
    assert_eq!(audit.external_heat_j,0.0);
    assert!(cooled.interfaces()[&7].point.failed());
    assert_eq!(state.thermal().links()[&7],link(7));
}

#[test]
fn reclosure_restores_separate_conductance_but_never_tensile_strength() {
    let law = law();
    let state = coupled(&law,0.0);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (open,_) = state.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    let compressed = broken.trial([-0.00001,0.0,0.0],1e-8).unwrap().point;
    let (closed,r) = open.synchronize(0.0,&[(7,compressed)],1e-8,|| Ok(())).unwrap();
    assert_eq!(r.reclosed,1);
    assert_eq!(r.fracture_dissipation_j,0.0);
    assert_eq!(closed.thermal().links()[&7].conductance_w_k,0.05);
    assert_eq!(closed.interfaces()[&7].point.damage(),1.0);
    let inside_band = compressed.trial([0.00005,0.0,0.0],1e-8).unwrap().point;
    assert_eq!(inside_band.traction_pa(),[0.0;3]);
    let (still_closed,_) = closed.synchronize(0.0,&[(7,inside_band)],1e-8,|| Ok(())).unwrap();
    assert!(still_closed.interfaces()[&7].closed);
    let (heated,_) = still_closed.advance_to(1.0,config(),|| Ok(())).unwrap();
    near(heated.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-0.05_f64).exp(),0.002);
    // The same gap stays open when approached from the open branch.
    let approach = broken.trial([0.00005,0.0,0.0],1e-8).unwrap().point;
    let (still_open,_) = open.synchronize(0.0,&[(7,approach)],1e-8,|| Ok(())).unwrap();
    assert!(!still_open.interfaces()[&7].closed);
}

#[test]
fn shear_failure_and_partial_damage_do_not_invent_missing_contact_area() {
    let law = law();
    let state = coupled(&law,0.0);
    let partial = law.bonded_point().trial([0.0005,0.0,0.0],1e-8).unwrap().point;
    let (damaged,_) = state.synchronize(0.0,&[(7,partial)],1e-8,|| Ok(())).unwrap();
    assert!(partial.damage() > 0.0 && !partial.failed());
    assert_eq!(damaged.thermal().links()[&7],link(7));
    let shear = law.bonded_point().trial([0.0,0.001,0.0],1e-8).unwrap().point;
    let (slipped,r) = state.synchronize(0.0,&[(7,shear)],1e-8,|| Ok(())).unwrap();
    assert_eq!(r.changed_conductances,1);
    assert_eq!(r.opened,0);
    assert_eq!(slipped.thermal().links()[&7].area_m2,link(7).area_m2);
    assert_eq!(slipped.thermal().links()[&7].conductance_w_k,0.05);
}

#[test]
fn inconsistent_history_law_ids_time_and_slip_refuse_without_mutation() {
    let law = law();
    let other = law.clone();
    let state = coupled(&law,0.0);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (open,_) = state.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    let excessive_slip = broken.trial([0.0,0.006,0.0],1e-8).unwrap().point;
    for endpoints in [vec![],vec![(8,broken)],vec![(7,law.bonded_point())],
        vec![(7,other.bonded_point())],vec![(7,excessive_slip)]] {
        assert!(open.synchronize(0.0,&endpoints,1e-8,|| Ok(())).is_err());
    }
    assert!(open.synchronize(1.0,&[(7,broken)],1e-8,|| Ok(())).is_err());
    // A larger excursion unseen by this owner cannot be smuggled in at unloading.
    let unload = broken.trial([0.0003,0.0,0.0],1e-8).unwrap().point;
    assert!(state.synchronize(0.0,&[(7,unload)],1e-8,|| Ok(())).is_err());
    assert!(open.interfaces()[&7].point.failed());
    assert!(open.thermal().links().is_empty());
    assert_eq!(open.thermal().energy_j(),state.thermal().energy_j());
}

#[test]
fn cancelled_contact_replacement_discards_the_entire_mechanical_proposal() {
    let law = law();
    let state = coupled(&law,0.0);
    let point = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    for allowed in [0,2,4,6,8] {
        let mut calls = 0;
        assert!(state.synchronize(0.0,&[(7,point)],1e-8,|| {
            calls += 1;
            if calls > allowed { Err(DepositionError::Cancelled) } else { Ok(()) }
        }).is_err());
        assert_eq!(state.interfaces()[&7].point.damage(),0.0);
        assert_eq!(state.thermal().links()[&7],link(7));
    }
}

#[test]
fn same_time_external_changes_cannot_silently_override_mechanical_ownership() {
    let law = law();
    let events = vec![TimedEvent { time_s: 0.5,event: Event::Separation(7) },
        TimedEvent { time_s: 0.5,event: Event::Contact(link(7)) }];
    let state = FractureCooling::new(thermal(&[link(7)],0.0,events),
        vec![binding(&law,link(7))],1,|| Ok(())).unwrap();
    assert!(state.advance_to(1.0,config(),|| Ok(())).is_err());
    assert_eq!(state.thermal().time_s(),0.0);
    assert_eq!(state.thermal().temperature_k(CellId(1)).unwrap(),500.0);
    assert_eq!(state.interfaces()[&7].point.damage(),0.0);
}

#[test]
fn checked_contact_updates_refuse_stale_duplicate_and_invalid_area_batches() {
    let state = thermal(&[link(7)],0.0,vec![]);
    let open = ContactUpdate { id: 7,expected: Some(link(7)),replacement: None };
    assert!(state.replace_contacts(&[open,open],|| Ok(())).is_err());
    assert!(state.replace_contacts(&[ContactUpdate { expected: None,..open }],|| Ok(())).is_err());
    let huge = ThermalLink { area_m2: 0.1,..link(7) };
    assert!(state.replace_contacts(&[ContactUpdate { replacement: Some(huge),..open }],|| Ok(())).is_err());
    assert_eq!(state.links()[&7],link(7));
    let opened = state.replace_contacts(&[open],|| Ok(())).unwrap();
    assert_eq!(state.cells(),opened.cells());
    assert_eq!(state.time_s(),opened.time_s());
    assert!(opened.replace_contacts(&[open],|| Ok(())).is_err());
}

#[test]
fn area_partition_and_same_grid_restart_preserve_heat_and_fracture_work() {
    let law = law();
    let whole = coupled(&law,100.0);
    let a = ThermalLink { area_m2: 0.0001,conductance_w_k: 0.1,..link(7) };
    let b = ThermalLink { id: 8,..a };
    let mut ia = binding(&law,a);
    ia.policy.reclosed_conductance_w_k = 0.025;
    let split = FractureCooling::new(thermal(&[a,b],100.0,vec![]),
        vec![ia,Interface { reference: b,..ia }],2,|| Ok(())).unwrap();
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (whole,wr) = whole.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    let (split,sr) = split.synchronize(0.0,&[(8,broken),(7,broken)],1e-8,|| Ok(())).unwrap();
    assert_eq!(wr.fracture_dissipation_j,sr.fracture_dissipation_j);
    let cfg = StepConfig { max_dt_s: 0.125,..config() };
    let (whole,_) = whole.advance_to(2.0,cfg,|| Ok(())).unwrap();
    let (prefix,_) = split.advance_to(1.0,cfg,|| Ok(())).unwrap();
    let (resumed,_) = prefix.clone().advance_to(2.0,cfg,|| Ok(())).unwrap();
    assert_eq!(whole.thermal().cells(),resumed.thermal().cells());
    assert_eq!(resumed.interfaces()[&7].point,broken);
    assert_eq!(resumed.interfaces()[&8].point,broken);
}

fn heat_partition(id: u64, fraction: f64, first: f64) -> FractureHeatPartition {
    FractureHeatPartition { interface: id, fraction_to_heat: fraction, fraction_to_first: first }
}
fn heat_fixture<'t, 'm>(law: &'m BilinearCohesive, caloric: &'t dyn ThermalLaw,
    masses: [f64; 2], enthalpies: [f64; 2],
) -> FractureCooling<'t, 'm> {
    let mut events = vec![birth(1,300.0,0.0),birth(2,300.0,0.0)];
    for (i, event) in events.iter_mut().enumerate() {
        if let Event::Birth(cell) = &mut event.event {
            cell.mass_kg = masses[i];
            cell.specific_h_j_kg = enthalpies[i];
        }
    }
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(link(7)) });
    let thermal = DepositionSimulation::new(events, &[(MaterialId(1),caloric)],
        DepositionBudget { max_cells: 2,max_links: 4,max_events: 16,max_substeps: 100 },
    ).unwrap().advance_to(0.0,config(),|| Ok(())).unwrap().0;
    FractureCooling::new(thermal,vec![binding(law,link(7))],1,|| Ok(())).unwrap()
}
struct LatentHeat;
impl ThermalLaw for LatentHeat {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=400.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(if h < 100.0 { 300.0+h/10.0 } else if h <= 200.0 { 310.0 }
            else { 310.0+(h-200.0)/10.0 })
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.1 }
}
// A bounded synthetic latent plateau used only for floating-point stress tests.
struct WidePlateau;
impl ThermalLaw for WidePlateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=36_028_797_018_963_968.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
}

#[test]
fn fracture_heat_partitions_new_work_into_original_unequal_mass_cells() {
    let law = law();
    let state = heat_fixture(&law,&HEAT,[0.001,0.002],[400_000.0,0.0]);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (next,r) = state.synchronize_with_heat(0.0,&[(7,broken)],
        &[heat_partition(7,0.75,0.25)],1e-8,1e-8,|| Ok(())).unwrap();
    // Independent triangular traction law: 100 J/m2 * 0.0002 m2 = 0.02 J.
    near(r.contact.fracture_dissipation_j,0.02,1e-12);
    near(r.heat_by_cell_j[&CellId(1)],0.00375,1e-12);
    near(r.heat_by_cell_j[&CellId(2)],0.01125,1e-12);
    near(r.nonthermal_fracture_energy_j,0.005,1e-12);
    near(next.thermal().temperature_k(CellId(1)).unwrap(),500.001875,1e-10);
    near(next.thermal().temperature_k(CellId(2)).unwrap(),300.0028125,1e-10);
    near(r.thermal.stored_energy_change_j,0.015,1e-10);
    near(r.transfer_residual_j,0.0,1e-10);
    assert_eq!(r.contact.opened,1);
    assert!(next.thermal().links().is_empty());
    assert_eq!(next.thermal().time_s(),0.0);
    assert_eq!((r.thermal.substeps,r.thermal.events),(0,0));
    assert_eq!(state.thermal().cells()[&CellId(1)].specific_h_j_kg,400_000.0);
}

#[test]
fn fracture_heat_uses_latent_enthalpy_instead_of_a_frozen_heat_capacity() {
    let law = law();
    let state = heat_fixture(&law,&LatentHeat,[0.001,0.002],[98.0,150.0]);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (next,r) = state.synchronize_with_heat(0.0,&[(7,broken)],
        &[heat_partition(7,1.0,0.25)],1e-8,1e-10,|| Ok(())).unwrap();
    near(next.thermal().cells()[&CellId(1)].specific_h_j_kg,103.0,1e-12);
    near(next.thermal().cells()[&CellId(2)].specific_h_j_kg,157.5,1e-12);
    assert_eq!(next.thermal().temperature_k(CellId(1)).unwrap(),310.0);
    assert_eq!(next.thermal().temperature_k(CellId(2)).unwrap(),310.0);
    near(r.thermal.stored_energy_change_j,0.02,1e-12);
    assert_eq!(r.nonthermal_fracture_energy_j,0.0);
}

#[test]
fn fracture_heat_does_not_repeat_on_replay_unloading_reclosure_or_imported_damage() {
    let law = law();
    let state = coupled(&law,0.0);
    let partition = [heat_partition(7,1.0,0.5)];
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (heated,_) = state.synchronize_with_heat(0.0,&[(7,broken)],&partition,
        1e-8,1e-8,|| Ok(())).unwrap();
    for point in [broken,broken.trial([0.0002,0.0,0.0],1e-8).unwrap().point,
        broken.trial([-0.00001,0.0,0.0],1e-8).unwrap().point] {
        let (next,r) = heated.synchronize_with_heat(0.0,&[(7,point)],&partition,
            1e-8,0.0,|| Ok(())).unwrap();
        assert_eq!(r.contact.fracture_dissipation_j,0.0);
        assert_eq!(r.thermal.external_heat_j,0.0);
        assert_eq!(next.thermal().cells(),heated.thermal().cells());
        assert_eq!(next.interfaces()[&7].point.damage(),1.0);
    }
    let partial = law.bonded_point().trial([0.0005,0.0,0.0],1e-8).unwrap().point;
    let imported = FractureCooling::new(thermal(&[link(7)],0.0,vec![]),
        vec![Interface { point: partial,..binding(&law,link(7)) }],1,|| Ok(())).unwrap();
    let (next,r) = imported.synchronize_with_heat(0.0,&[(7,partial)],&partition,
        1e-8,0.0,|| Ok(())).unwrap();
    assert_eq!(r.thermal.external_heat_j,0.0);
    assert_eq!(next.thermal().cells(),imported.thermal().cells());
}

#[test]
fn zero_fracture_heat_fraction_keeps_all_work_in_the_nonthermal_account() {
    let law = law();
    let state = coupled(&law,0.0);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (next,r) = state.synchronize_with_heat(0.0,&[(7,broken)],
        &[heat_partition(7,0.0,0.3)],1e-8,0.0,|| Ok(())).unwrap();
    assert_eq!(next.thermal().cells(),state.thermal().cells());
    assert!(next.thermal().links().is_empty());
    near(r.nonthermal_fracture_energy_j,0.02,1e-12);
    assert_eq!(r.thermal.external_heat_j,0.0);
    assert_eq!(r.transfer_residual_j,0.0);
}

#[test]
fn fracture_heat_shared_cell_sums_are_order_invariant_and_preserve_area_partition() {
    let law = law();
    let whole = coupled(&law,0.0);
    let a = ThermalLink { area_m2: 0.0001,conductance_w_k: 0.1,..link(7) };
    let b = ThermalLink { id: 8,..a };
    let split = FractureCooling::new(thermal(&[a,b],0.0,vec![]),
        vec![binding(&law,a),binding(&law,b)],2,|| Ok(())).unwrap();
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (one,wr) = whole.synchronize_with_heat(0.0,&[(7,broken)],
        &[heat_partition(7,0.5,0.25)],1e-8,1e-8,|| Ok(())).unwrap();
    let parts = [heat_partition(7,0.5,0.25),heat_partition(8,0.5,0.25)];
    let (two,sr) = split.synchronize_with_heat(0.0,&[(7,broken),(8,broken)],
        &parts,1e-8,1e-8,|| Ok(())).unwrap();
    let (reversed,rr) = split.synchronize_with_heat(0.0,&[(8,broken),(7,broken)],
        &[parts[1],parts[0]],1e-8,1e-8,|| Ok(())).unwrap();
    assert_eq!(one.thermal().cells(),two.thermal().cells());
    assert_eq!(two.thermal().cells(),reversed.thermal().cells());
    assert_eq!(wr.heat_by_cell_j,sr.heat_by_cell_j);
    assert_eq!(sr.heat_by_cell_j,rr.heat_by_cell_j);
    assert_eq!(wr.contact.fracture_dissipation_j,sr.contact.fracture_dissipation_j);
    assert_eq!(sr.transfer_residual_j,rr.transfer_residual_j);
}

#[test]
fn fracture_heat_invalid_partitions_domains_energy_and_late_cancel_are_atomic() {
    let law = law();
    let state = coupled(&law,0.0);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let good = heat_partition(7,1.0,1.0);
    for parts in [vec![],vec![good,good],vec![heat_partition(8,1.0,0.5)],
        vec![heat_partition(7,-0.1,0.5)],vec![heat_partition(7,1.1,0.5)],
        vec![heat_partition(7,f64::NAN,0.5)],vec![heat_partition(7,1.0,f64::INFINITY)]] {
        assert!(state.synchronize_with_heat(0.0,&[(7,broken)],&parts,
            1e-8,1e-8,|| Ok(())).is_err());
    }
    let hot = heat_fixture(&law,&HEAT,[0.001,0.001],[600_000.0,0.0]);
    assert!(matches!(hot.synchronize_with_heat(0.0,&[(7,broken)],&[good],
        1e-8,1e-8,|| Ok(())),Err(Error::Thermal(DepositionError::MaterialDomain))));
    let coarse = heat_fixture(&law,&WidePlateau,[1.0,1.0],[9_007_199_254_740_992.0,0.0]);
    assert!(matches!(coarse.synchronize_with_heat(0.0,&[(7,broken)],&[good],
        1e-8,0.0,|| Ok(())),Err(Error::Thermal(DepositionError::EnergyBalance(_)))));
    let mut total_calls = 0;
    state.synchronize_with_heat(0.0,&[(7,broken)],&[good],1e-8,1e-8,|| {
        total_calls += 1; Ok(())
    }).unwrap();
    for allowed in [0,total_calls-1] {
        let mut calls = 0;
        assert!(matches!(state.synchronize_with_heat(0.0,&[(7,broken)],&[good],
            1e-8,1e-8,|| { calls += 1;
                if calls > allowed { Err(DepositionError::Cancelled) } else { Ok(()) }
            }),Err(Error::Thermal(DepositionError::Cancelled))));
    }
    assert_eq!(state.interfaces()[&7].point.damage(),0.0);
    assert_eq!(state.thermal().links()[&7],link(7));
    assert_eq!(state.thermal().temperature_k(CellId(1)).unwrap(),500.0);
    assert_eq!(coarse.thermal().cells()[&CellId(1)].specific_h_j_kg,9_007_199_254_740_992.0);
}

#[test]
fn fracture_heat_reports_and_budgets_disparate_work_and_nonthermal_rounding() {
    let large = 9_007_199_254_740_992.0; // 2^53; adding 1 J rounds away.
    let law = BilinearCohesive::new(4.0,2.0,1.0,4.0,2.0).unwrap();
    let links = [ThermalLink { area_m2: large,conductance_w_k: 0.0,..link(7) },
        ThermalLink { id: 8,area_m2: 1.0,conductance_w_k: 0.0,..link(7) },
        ThermalLink { id: 9,area_m2: 1.0,conductance_w_k: 0.0,..link(7) }];
    let mut events = vec![birth(1,300.0,0.0),birth(2,300.0,0.0)];
    for event in &mut events {
        if let Event::Birth(cell) = &mut event.event {
            cell.mass_kg = 1.0;
            cell.geometry.surface_area_m2 = 2.0*large;
        }
    }
    events.extend(links.iter().map(|&l| TimedEvent { time_s: 0.0,event: Event::Contact(l) }));
    let heat = DepositionSimulation::new(events,&[(MaterialId(1),&WidePlateau)],
        DepositionBudget { max_cells: 2,max_links: 4,max_events: 16,max_substeps: 100 },
    ).unwrap().advance_to(0.0,config(),|| Ok(())).unwrap().0;
    let state = FractureCooling::new(heat,links.iter().map(|&l| binding(&law,l)).collect(),
        3,|| Ok(())).unwrap();
    let broken = law.bonded_point().trial([1.0,0.0,0.0],1e-8).unwrap().point;
    let points = [(7,broken),(8,broken),(9,broken)];
    let parts = [heat_partition(7,0.0,1.0),heat_partition(8,0.0,1.0),heat_partition(9,1.0,1.0)];
    assert!(matches!(state.synchronize_with_heat(0.0,&points,&parts,1e-8,0.5,|| Ok(())),
        Err(Error::Thermal(DepositionError::EnergyBalance(_)))));
    let (next,r) = state.synchronize_with_heat(0.0,&points,&parts,1e-8,1.0,|| Ok(())).unwrap();
    assert_eq!(r.contact.fracture_dissipation_j,large+2.0);
    assert_eq!(r.nonthermal_fracture_energy_j,large);
    assert_eq!(r.thermal.stored_energy_change_j,1.0);
    assert_eq!(r.thermal.external_heat_j,1.0);
    assert_eq!(r.transfer_residual_j,-1.0);
    assert_eq!(next.thermal().cells()[&CellId(1)].specific_h_j_kg,1.0);
    // The old contact-only API retains its original arithmetic and no heating.
    let (legacy,old_report) = state.synchronize(0.0,&points,1e-8,|| Ok(())).unwrap();
    assert_eq!(old_report.fracture_dissipation_j,large);
    assert_eq!(legacy.thermal().cells(),state.thermal().cells());
    assert_eq!(state.thermal().cells()[&CellId(1)].specific_h_j_kg,0.0);
}
