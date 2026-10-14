use super::*;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> { Ok(300.0+h/2000.0) }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static LAW: Sensible = Sensible;
fn config(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-8 } }
fn link() -> ThermalLink {
    ThermalLink { id: 7, a: CellId(1), b: CellId(2), area_m2: 0.0003, conductance_w_k: 0.2 }
}
fn patch() -> ContactPatch { ContactPatch { link: link(), covered_area_m2: [0.0002,0.0008] } }
fn initial(htc: f64, extra: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for (id,t) in [(1,500.0),(2,300.0)] {
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(id), segment_id: id, body_id: id, material_id: MaterialId(1),
            mass_kg: 0.001, specific_h_j_kg: 2000.0*(t-300.0),
            geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],
                end_m: [0.01,0.0,id as f64*0.001], volume_m3: 1e-6, surface_area_m2: 0.001 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: htc, reservoir: None, source_w: 0.0 },
        }) });
    }
    events.extend(extra);
    DepositionSimulation::new(events,&[(MaterialId(1),&LAW)],DepositionBudget {
        max_cells: 3,max_links: 8,max_events: 16,max_substeps: 10000,
    }).unwrap().advance_to(0.0,config(0.1),|| Ok(())).unwrap().0
}
fn change(s: &DepositionSimulation<'_>, replacement: Option<ContactPatch>) -> ContactPatchUpdate {
    ContactPatchUpdate { id: 7, expected: s.contact_patch(7), replacement }
}
fn near(a: f64,b: f64) { assert!((a-b).abs()<1e-9,"{a:e} != {b:e}"); }

// G1: independent one-step balance, with different exposed reference areas.
#[test]
fn side_coverage_controls_actual_ambient_loss_but_not_exchange_reciprocity() {
    let s = initial(100.0,vec![]);
    let s = s.replace_contact_patches(&[change(&s,Some(patch()))],|| Ok(())).unwrap();
    near(s.exposed_area_m2(CellId(1)).unwrap(),0.0008);
    near(s.exposed_area_m2(CellId(2)).unwrap(),0.0002);
    let (next,a) = s.advance_to(0.1,config(0.1),|| Ok(())).unwrap();
    near(next.temperature_k(CellId(1)).unwrap(),497.2);
    near(next.temperature_k(CellId(2)).unwrap(),302.0);
    near(a.external_heat_j,-1.6);
    near(a.residual_j,0.0);
    assert_eq!(next.contact_patch(7),Some(patch()));
    let s = initial(0.0,vec![]);
    let s = s.replace_contact_patches(&[change(&s,Some(patch()))],|| Ok(())).unwrap();
    let (next,a) = s.advance_to(1.0,config(0.01),|| Ok(())).unwrap();
    near(next.energy_j(),s.energy_j());
    assert_eq!(a.external_heat_j,0.0);
}

#[test]
fn zero_conductance_still_covers_each_side_and_updates_the_cfl_bound() {
    let s = initial(100.0,vec![]);
    let p = ContactPatch { link: ThermalLink { conductance_w_k: 0.0,..link() },..patch() };
    let s = s.replace_contact_patches(&[change(&s,Some(p))],|| Ok(())).unwrap();
    let (next,a) = s.advance_to(20.0,config(100.0),|| Ok(())).unwrap();
    assert_eq!(a.substeps,1); // dt = 0.8*C/(h*A_exposed) = 20 s, not nominal-area 22.86.
    near(next.temperature_k(CellId(1)).unwrap(),340.0);
    near(next.temperature_k(CellId(2)).unwrap(),300.0);
    let (open,_) = next.advance_to(20.0,config(1.0),|| Ok(())).unwrap();
    let open = open.replace_contact_patches(&[change(&open,None)],|| Ok(())).unwrap();
    near(open.exposed_area_m2(CellId(1)).unwrap(),0.001);
    near(open.exposed_area_m2(CellId(2)).unwrap(),0.001);
}

// G3: symmetric API equivalence and exchange-orientation invariance.
#[test]
fn symmetric_existing_api_and_reversed_pair_reproduce_the_same_heat() {
    let s = initial(100.0,vec![]);
    let old = s.replace_contacts(&[ContactUpdate { id: 7,expected: None,replacement: Some(link()) }],|| Ok(())).unwrap();
    let new = s.replace_contact_patches(&[change(&s,Some(ContactPatch::symmetric(link())))],|| Ok(())).unwrap();
    let (a,aa) = old.advance_to(2.0,config(0.125),|| Ok(())).unwrap();
    let (b,ba) = new.advance_to(2.0,config(0.125),|| Ok(())).unwrap();
    assert_eq!(a.cells(),b.cells()); assert_eq!(aa,ba);
    let p = patch();
    let reverse = ContactPatch { link: ThermalLink { a: p.link.b,b: p.link.a,..p.link },
        covered_area_m2: [p.covered_area_m2[1],p.covered_area_m2[0]] };
    let a = s.replace_contact_patches(&[change(&s,Some(p))],|| Ok(())).unwrap()
        .advance_to(2.0,config(0.125),|| Ok(())).unwrap().0;
    let b = s.replace_contact_patches(&[change(&s,Some(reverse))],|| Ok(())).unwrap()
        .advance_to(2.0,config(0.125),|| Ok(())).unwrap().0;
    assert_eq!(a.cells(),b.cells());
}

#[test]
fn stale_side_areas_and_symmetric_proposals_cannot_erase_an_asymmetric_patch() {
    let s = initial(100.0,vec![]);
    let s = s.replace_contact_patches(&[change(&s,Some(patch()))],|| Ok(())).unwrap();
    let stale = change(&s,None);
    let changed = s.replace_contact_patches(&[change(&s,Some(ContactPatch {
        covered_area_m2: [0.0004,0.0006],..patch() }))],|| Ok(())).unwrap();
    assert_eq!(changed.links(),s.links());
    assert!(changed.replace_contact_patches(&[stale],|| Ok(())).is_err());
    assert!(s.replace_contacts(&[ContactUpdate { id: 7,expected: Some(link()),
        replacement: Some(link()) }],|| Ok(())).is_err());
    assert_eq!(s.contact_patch(7),Some(patch()));
}

#[test]
fn separation_and_rebirth_do_not_reuse_old_side_areas() {
    let s = initial(100.0,vec![
        TimedEvent { time_s: 1.0,event: Event::Separation(7) },
        TimedEvent { time_s: 1.0,event: Event::Contact(link()) },
    ]);
    let s = s.replace_contact_patches(&[change(&s,Some(patch()))],|| Ok(())).unwrap();
    let (next,_) = s.advance_to(1.0,config(0.125),|| Ok(())).unwrap();
    assert_eq!(next.contact_patch(7),Some(ContactPatch::symmetric(link())));
    near(next.exposed_area_m2(CellId(1)).unwrap(),0.0007);
    near(next.exposed_area_m2(CellId(2)).unwrap(),0.0007);
    assert!(next.covered_areas.is_empty());
}

// G0/G4: reject invalid second-side area even for insulation; no partial batch.
#[test]
fn invalid_area_topology_budget_or_cancelled_batch_keeps_the_original_state() {
    let s = initial(100.0,vec![]);
    for p in [ContactPatch { covered_area_m2: [0.0001,0.002],..patch() },
        ContactPatch { covered_area_m2: [0.0,0.0002],..patch() },
        ContactPatch { covered_area_m2: [0.0001,f64::NAN],..patch() },
        ContactPatch { link: ThermalLink { b: CellId(99),..link() },..patch() }] {
        assert!(s.replace_contact_patches(&[change(&s,Some(p))],|| Ok(())).is_err());
        assert!(s.links().is_empty()); assert!(s.covered_areas.is_empty());
    }
    let good = change(&s,Some(patch()));
    let bad = ContactPatchUpdate { id: 8,expected: None,replacement: Some(ContactPatch {
        link: ThermalLink { id: 8,conductance_w_k: 0.0,..link() },
        covered_area_m2: [0.0009,0.0002] }) };
    assert!(s.replace_contact_patches(&[good,bad],|| Ok(())).is_err());
    assert!(s.replace_contact_patches(&[good,good],|| Ok(())).is_err());
    let mut limited = s.clone(); limited.budget.max_events = 0;
    assert!(limited.replace_contact_patches(&[good],|| Ok(())).is_err());
    let mut count = 0;
    s.replace_contact_patches(&[good],|| { count+=1; Ok(()) }).unwrap();
    let mut calls = 0;
    assert!(s.replace_contact_patches(&[good],|| { calls+=1;
        if calls == count { Err(DepositionError::Cancelled) } else { Ok(()) } }).is_err());
    assert!(s.links().is_empty()); near(s.energy_j(),400.0);
}

// G3/G5: area partition and existing-grid restart preserve heat/coverage.
#[test]
fn patch_partition_and_checkpoint_replay_preserve_areas_and_energy() {
    let s = initial(100.0,vec![]);
    let p = patch();
    let whole = s.replace_contact_patches(&[change(&s,Some(p))],|| Ok(())).unwrap();
    let half = ContactPatch { link: ThermalLink { area_m2: 0.00015,conductance_w_k: 0.1,..link() },
        covered_area_m2: [0.0001,0.0004] };
    let split = s.replace_contact_patches(&[change(&s,Some(half)),ContactPatchUpdate {
        id: 8,expected: None,replacement: Some(ContactPatch {
            link: ThermalLink { id: 8,..half.link },..half }) }],|| Ok(())).unwrap();
    let (a,_) = whole.advance_to(2.0,config(0.125),|| Ok(())).unwrap();
    let (b,_) = split.advance_to(2.0,config(0.125),|| Ok(())).unwrap();
    for id in [CellId(1),CellId(2)] { near(a.temperature_k(id).unwrap(),b.temperature_k(id).unwrap()); }
    let mid = whole.advance_to(1.0,config(0.125),|| Ok(())).unwrap().0;
    let replay = mid.clone().advance_to(2.0,config(0.125),|| Ok(())).unwrap().0;
    assert_eq!(a.cells(),replay.cells()); assert_eq!(a.contact_patch(7),replay.contact_patch(7));
}
