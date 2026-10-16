//! Actual enthalpy kernel with state-resolved links, not a replacement stencil.
use std::collections::BTreeMap;
use fs_conduction::transient::deposition::*;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(-100.0..=100.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0 + h)
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
}
static LAW: Caloric = Caloric;
fn config(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-10 } }
fn cell(id: u64, h: f64) -> DepositedCell {
    DepositedCell { id:CellId(id), segment_id:id, body_id:1, material_id:MaterialId(1),
        mass_kg:1.0, specific_h_j_kg:h,
        geometry:CellGeometry {start_m:[id as f64,0.0,0.0],end_m:[id as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:6.0},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0} }
}
fn link(g: f64) -> ThermalLink { ThermalLink {id:7,a:CellId(1),b:CellId(2),area_m2:1.0,conductance_w_k:g} }
fn simulation(tail: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    let mut events=vec![
        TimedEvent {time_s:0.0,event:Event::Birth(cell(1,0.0))},
        TimedEvent {time_s:0.0,event:Event::Birth(cell(2,10.0))},
        TimedEvent {time_s:0.0,event:Event::Contact(link(1.0))},
    ];
    events.extend(tail);
    DepositionSimulation::new(events,&[(MaterialId(1),&LAW)],DepositionBudget {
        max_cells:2,max_links:1,max_events:10,max_substeps:100_000,
    }).unwrap().advance_to(0.0,config(0.01),||Ok(())).unwrap().0
}
#[derive(Clone)]
struct Transport {
    g: Option<f64>, nonlinear:bool, calls:usize, fail_after:Option<usize>,
    times:Vec<f64>, pending:bool, observed:usize,
}
impl Transport {
    fn fixed(g:Option<f64>) -> Self {
        Self {g,nonlinear:false,calls:0,fail_after:None,times:vec![],pending:false,observed:0}
    }
}
impl DepositionObserver for Transport {
    type Error=DepositionError;
    fn sources(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {
        assert!(!self.pending); self.pending=true; Ok(BTreeMap::new())
    }
    fn conductances(&mut self,s:&DepositionSimulation<'_>,check:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<u64,f64>,Self::Error> {
        check()?; assert!(self.pending);
        if self.fail_after==Some(self.calls) {return Err(DepositionError::MaterialDomain);}
        self.calls+=1; self.times.push(s.time_s());
        let mut out=BTreeMap::new();
        for (&id,l) in s.links() {
            let g=if self.nonlinear {
                1.0+(s.temperature_k(l.b)?-s.temperature_k(l.a)?).abs()/10.0
            } else {self.g.unwrap_or(l.conductance_w_k)};
            out.insert(id,g);
        }
        Ok(out)
    }
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,s:&DepositionSimulation<'_>,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        assert!(self.pending); self.pending=false; self.observed+=1;
        // No dynamic conductance leaks into authored contact templates.
        for l in s.links().values() {assert!(l.conductance_w_k==1.0 || l.conductance_w_k==3.0);}
        Ok(())
    }
}

#[test]
fn constant_transport_matches_original_and_does_not_persist() {
    let old=simulation(vec![]);
    let (plain,a)=old.advance_to(0.3,config(0.01),||Ok(())).unwrap();
    let observer=Transport::fixed(Some(1.0));
    let (local,o,b)=old.advance_observed(&observer,0.3,config(0.01),||Ok(())).unwrap();
    assert_eq!(plain.cells(),local.cells()); assert_eq!(a,b);
    assert_eq!(o.calls,b.substeps);assert_eq!(o.observed,o.calls);assert_eq!(observer.calls,0);
    assert_eq!(local.contact_patch(7),old.contact_patch(7));
    let (next,_,_)=old.advance_observed(&Transport::fixed(Some(5.0)),0.1,config(0.01),||Ok(())).unwrap();
    assert!(next.temperature_k(CellId(1)).unwrap()>plain.temperature_k(CellId(1)).unwrap());
    assert_eq!(next.contact_patch(7),old.contact_patch(7));
}

#[test]
fn applied_conductance_controls_the_actual_cfl_and_conserves_heat() {
    let old=simulation(vec![]);
    let (next,o,a)=old.advance_observed(&Transport::fixed(Some(1000.0)),0.1,config(1.0),||Ok(())).unwrap();
    assert!(a.substeps>100);assert_eq!(a.substeps,o.calls);
    for c in next.cells().keys() {assert!((300.0..=310.0).contains(&next.temperature_k(*c).unwrap()));}
    assert!(a.residual_j.abs()<1e-10);assert_eq!(a.external_heat_j,0.0);
    assert_eq!(next.links()[&7],link(1.0));
}

#[test]
fn nonlinear_transport_is_resampled_with_first_order_time_refinement() {
    let old=simulation(vec![]);let mut observer=Transport::fixed(None);observer.nonlinear=true;
    // dD/dt=-2 D(1+D/10), with D(0)=10; mean temperature remains 305 K.
    let e=(-1.0_f64).exp(); let exact=305.0-5.0*e/(2.0-e);
    let mut errors=Vec::new();
    for dt in [0.01,0.005,0.0025] {
        let (next,o,a)=old.advance_observed(&observer,0.5,config(dt),||Ok(())).unwrap();
        errors.push((next.temperature_k(CellId(1)).unwrap()-exact).abs());
        assert_eq!(o.calls,a.substeps);assert!(a.residual_j.abs()<1e-10);
        assert!(o.times.windows(2).all(|t|t[0]<t[1]));
    }
    for p in errors.windows(2) {assert!((1.9..2.1).contains(&(p[0]/p[1])));}
}

#[test]
fn insulation_preserves_asymmetric_coverage_and_scheduled_contact_replacement() {
    let old=simulation(vec![
        TimedEvent {time_s:0.1,event:Event::Separation(7)},
        TimedEvent {time_s:0.2,event:Event::Contact(link(3.0))},
    ]);
    let original=old.contact_patch(7).unwrap();let mut patch=original;patch.covered_area_m2=[0.5,1.5];
    let old=old.replace_contact_patches(&[ContactPatchUpdate {id:7,expected:Some(original),replacement:Some(patch)}],||Ok(())).unwrap();
    let (held,_,_)=old.advance_observed(&Transport::fixed(Some(0.0)),0.05,config(1.0),||Ok(())).unwrap();
    assert_eq!(held.cells(),old.cells());assert_eq!(held.contact_patch(7),Some(patch));
    assert_eq!(held.exposed_area_m2(CellId(1)).unwrap(),5.5);
    assert_eq!(held.exposed_area_m2(CellId(2)).unwrap(),4.5);
    let (next,o,a)=old.advance_observed(&Transport::fixed(None),0.3,config(0.08),||Ok(())).unwrap();
    assert!(o.times.contains(&0.1) && o.times.contains(&0.2));assert_eq!(a.events,2);
    assert_eq!(next.contact_patch(7),Some(ContactPatch::symmetric(link(3.0))));
}

#[test]
fn bad_conductance_late_material_failure_and_cancellation_are_atomic() {
    let old=simulation(vec![]);let cells=old.cells().clone();let patch=old.contact_patch(7);
    for g in [f64::NAN,f64::INFINITY,-1.0] {
        assert!(old.advance_observed(&Transport::fixed(Some(g)),0.1,config(0.01),||Ok(())).is_err());
    }
    let mut observer=Transport::fixed(Some(3.0));observer.fail_after=Some(2);
    assert!(matches!(old.advance_observed(&observer,0.1,config(0.01),||Ok(())),Err(DepositionError::MaterialDomain)));
    let mut polls=0;
    let cancelled=old.advance_observed(&Transport::fixed(Some(3.0)),0.1,config(0.01),|| {
        polls+=1;if polls>25 {Err(DepositionError::Cancelled)} else {Ok(())}
    });
    assert!(matches!(cancelled,Err(DepositionError::Cancelled)));
    assert_eq!(old.cells(),&cells);assert_eq!(old.contact_patch(7),patch);assert_eq!(old.time_s(),0.0);
    assert_eq!(observer.calls,0);
    assert!(old.advance_to(0.1,config(0.01),||Ok(())).is_ok());
}
