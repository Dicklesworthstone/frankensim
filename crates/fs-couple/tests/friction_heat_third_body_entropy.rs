//! Finite-layer entropy with independent three-capacity balances.
//! Synthetic caloric/transport fixtures are not measured-material validation.
use fs_couple::friction_heat::{FrictionHeatEntropy, third_body::*};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 { Ok(300.0 + h) }
        else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
    fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64)
        -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b-a)/t).ln_1p()))
    }
}
static SENSIBLE: Sensible = Sensible;
const CAPACITY: [f64; 3] = [2.0, 3.0, 0.25];
const K: [f64; 3] = [4.0, 1.0, 2.0];
const RESISTANCE: [f64; 2] = [0.2, 0.8];
const AMBIENT: [f64; 3] = [285.0, 320.0, 290.0];
const LINKS: [u64; 2] = [7, 8];
const LAYER: CellId = CellId(2);
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol, "{a:.16e} != {b:.16e}; tolerance {tol:e}");
}
fn controls() -> ThirdBodyHeatTolerances {
    ThirdBodyHeatTolerances { cell_j: 1e-10, thermal_j: 1e-10,
        layer_j: 1e-10, total_j: 1e-10 }
}
fn fixture(layer_law: &'static dyn ThermalLaw, h: [f64; 3], k: [f64; 3], htc: [f64; 3])
    -> (DepositionSimulation<'static>, ThirdBodyHeatRoute) {
    let mut events = Vec::new();
    for i in 0..3 {
        let x = [-1.0, 1.0, 0.0][i];
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: i as u64, body_id: i as u64,
            material_id: MaterialId(if i == 2 { 1 } else { 0 }),
            mass_kg: CAPACITY[i], specific_h_j_kg: h[i],
            geometry: CellGeometry { start_m: [x,0.0,0.0], end_m: [x+1.0,0.0,0.0],
                volume_m3: 1.0, surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: AMBIENT[i], htc_w_m2_k: htc[i],
                reservoir: None, source_w: 0.0 },
        }) });
    }
    for i in 0..2 {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: LINKS[i], a: CellId(i as u64), b: LAYER, area_m2: 1.0,
            conductance_w_k: series_contact_conductance_w_k(1.0, [0.5;2],
                [k[i],k[2]], RESISTANCE[i]).unwrap(),
        }) });
    }
    let state = DepositionSimulation::new(events,
        &[(MaterialId(0), &SENSIBLE as &dyn ThermalLaw), (MaterialId(1), layer_law)],
        DepositionBudget { max_cells: 3, max_links: 2, max_events: 5, max_substeps: 1000 })
        .unwrap().advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-10 },
            || Ok(())).unwrap().0;
    let route = ThirdBodyHeatRoute::from_links(&state, LAYER, LINKS, k, RESISTANCE).unwrap();
    (state,route)
}
fn source(work: f64, k: [f64;3], extra: [f64;3]) -> ThirdBodyHeatInput<()> {
    ThirdBodyHeatInput { dissipated_j: work, conductivity_w_m_k: k,
        area_resistance_m2_k_w: RESISTANCE, transport: StepCoefficients {
            sources_w: extra.into_iter().enumerate().filter(|(_,q)| *q != 0.0)
                .map(|(i,q)| (CellId(i as u64),q)).collect(),
            conductances_w_k: BTreeMap::new(), payload: (),
        } }
}
// Eliminate body temperatures analytically; this does not call the heat owner.
fn endpoint(old: &DepositionSimulation<'_>, dt: f64, work: f64, k: [f64;3], extra: [f64;3]) -> [f64;3] {
    let g: [f64;2] = std::array::from_fn(|i| if k[i] == 0.0 || k[2] == 0.0 { 0.0 }
        else { 1.0/(0.5/k[i]+0.5/k[2]+RESISTANCE[i]) });
    let bath: [f64;3] = std::array::from_fn(|i| [5.0,5.0,4.0][i]
        * old.cells()[&CellId(i as u64)].boundary.htc_w_m2_k);
    let b: [f64;3] = std::array::from_fn(|i| CAPACITY[i]/dt
        * old.temperature_k(CellId(i as u64)).unwrap()+bath[i]*AMBIENT[i]+extra[i]);
    let a = [CAPACITY[0]/dt+bath[0]+g[0], CAPACITY[1]/dt+bath[1]+g[1]];
    let layer = (b[2]+work/dt+g[0]*b[0]/a[0]+g[1]*b[1]/a[1])
        /(CAPACITY[2]/dt+bath[2]+g[0]+g[1]-g[0]*g[0]/a[0]-g[1]*g[1]/a[1]);
    let t = [(b[0]+g[0]*layer)/a[0], (b[1]+g[1]*layer)/a[1], layer];
    std::array::from_fn(|i| CAPACITY[i]*(t[i]-old.temperature_k(CellId(i as u64)).unwrap()))
}
fn balance(e: FrictionHeatEntropy) {
    close(e.thermal_storage_change_j_k, e.reservoir_entropy_input_j_k
        +e.network_conduction_generation_j_k+e.interface_conduction_generation_j_k
        +e.friction_generation_j_k+e.unclassified_source_at_cell_j_k
        +e.temporal_defect_j_k+e.interface_residual_entropy_j_k, 2e-12);
    assert_eq!(e.interface_conduction_generation_j_k,0.0);
    assert_eq!(e.interface_residual_entropy_j_k,0.0);
}

#[test]
fn finite_layer_entropy_keeps_storage_and_counts_each_contact_once() {
    let (old,route)=fixture(&SENSIBLE,[0.0;3],K,[0.0;3]);
    let energy=endpoint(&old,0.1,0.4,K,[0.0;3]);
    let step=ThirdBodyHeatStep::new(&old,0.1,&route,||Ok(())).unwrap();
    let (plain,(),ordinary)=step.evaluate(&energy,
        |_,_|Ok::<_,DepositionError>(source(0.4,K,[0.0;3])),||Ok(()))
        .unwrap().accept(controls(),||Ok(())).unwrap();
    let calls=Cell::new(0);
    let (next,(),report)=step.evaluate_with_entropy(&energy,|_,_| {
        calls.set(calls.get()+1); Ok::<_,DepositionError>(source(0.4,K,[0.0;3]))
    },||Ok(())).unwrap().accept(controls(),||Ok(())).unwrap();
    assert_eq!(calls.get(),1);
    assert_eq!(next.cells(),plain.cells()); assert_eq!(report.thermal,ordinary.thermal);
    assert!(ordinary.entropy.is_none());
    assert!(report.layer_stored_change_j>0.0);
    assert!(report.body_heat_j.iter().sum::<f64>()<report.dissipated_j);
    let e=report.entropy.unwrap();
    let exact:f64=(0..3).map(|i|CAPACITY[i]*(next.temperature_k(CellId(i as u64)).unwrap()/300.0).ln()).sum();
    let conduction:f64=(0..2).map(|i|report.body_heat_j[i]
        *(1.0/next.temperature_k(CellId(i as u64)).unwrap()-1.0/report.layer_temperature_k)).sum();
    close(e.thermal_storage_change_j_k,exact,2e-12);
    close(e.network_conduction_generation_j_k,conduction,2e-15);
    close(e.friction_generation_j_k,0.4/report.layer_temperature_k,1e-15);
    assert!(!e.has_unclassified_sources); assert!(e.temporal_defect_j_k>0.0);
    balance(e);
}

#[test]
fn insulated_finite_layer_stores_all_friction_without_fictitious_conduction() {
    let (old,route)=fixture(&SENSIBLE,[0.0;3],[0.0;3],[0.0;3]);
    let step=ThirdBodyHeatStep::new(&old,0.5,&route,||Ok(())).unwrap();
    let (next,(),r)=step.evaluate_with_entropy(&[0.0,0.0,0.8],
        |_,_|Ok::<_,DepositionError>(source(0.8,[0.0;3],[0.0;3])),||Ok(()))
        .unwrap().accept(controls(),||Ok(())).unwrap();
    let e=r.entropy.unwrap();
    assert_eq!(r.body_heat_j,[0.0;2]); close(r.layer_stored_change_j,0.8,1e-15);
    assert_eq!(e.network_conduction_generation_j_k,0.0);
    close(e.thermal_storage_change_j_k,0.25*(303.2_f64/300.0).ln(),2e-15);
    assert_eq!(next.temperature_k(CellId(0)).unwrap(),300.0);
    balance(e);
}

#[test]
fn cooling_layer_transfers_entropy_without_inventing_new_friction() {
    let (old,route)=fixture(&SENSIBLE,[10.0,-10.0,30.0],K,[0.0;3]);
    let energy=endpoint(&old,0.2,0.0,K,[0.0;3]);
    let step=ThirdBodyHeatStep::new(&old,0.2,&route,||Ok(())).unwrap();
    let (next,(),r)=step.evaluate_with_entropy(&energy,
        |_,_|Ok::<_,DepositionError>(source(0.0,K,[0.0;3])),||Ok(()))
        .unwrap().accept(controls(),||Ok(())).unwrap();
    let e=r.entropy.unwrap();
    assert!(r.layer_stored_change_j<0.0); assert!(e.network_conduction_generation_j_k>0.0);
    assert_eq!(e.friction_generation_j_k,0.0); assert!(e.thermal_storage_change_j_k>0.0);
    close(next.energy_j(),old.energy_j(),1e-10); balance(e);
}

#[test]
fn reservoirs_and_independent_sources_keep_their_own_entropy_semantics() {
    let (old,route)=fixture(&SENSIBLE,[5.0,-3.0,8.0],K,[0.2,0.3,0.4]);
    let extra=[0.1,-0.2,0.3]; let dt=0.2;
    let energy=endpoint(&old,dt,0.4,K,extra);
    let step=ThirdBodyHeatStep::new(&old,dt,&route,||Ok(())).unwrap();
    let (next,(),r)=step.evaluate_with_entropy(&energy,
        |_,_|Ok::<_,DepositionError>(source(0.4,K,extra)),||Ok(()))
        .unwrap().accept(controls(),||Ok(())).unwrap();
    let e=r.entropy.unwrap();
    let reservoir:f64=(0..3).map(|i|dt*[5.0,5.0,4.0][i]*[0.2,0.3,0.4][i]
        *(AMBIENT[i]-next.temperature_k(CellId(i as u64)).unwrap())/AMBIENT[i]).sum();
    let other:f64=(0..3).map(|i|dt*extra[i]/next.temperature_k(CellId(i as u64)).unwrap()).sum();
    close(e.reservoir_entropy_input_j_k,reservoir,2e-14);
    close(e.unclassified_source_at_cell_j_k,other,2e-15);
    assert!(e.has_unclassified_sources); balance(e);
}

#[test]
fn latent_layer_storage_is_integrated_once_on_its_original_chart() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self,h:f64)->Result<f64,DepositionError> {
            if (0.0..=100.0).contains(&h) {Ok(300.0)} else {Err(DepositionError::MaterialDomain)}
        }
        fn maximum_temperature_slope(&self)->f64 {0.0}
        fn specific_entropy_change_j_kg_k(&self,a:f64,b:f64)->Result<Option<f64>,DepositionError> {
            self.temperature_k(a)?;self.temperature_k(b)?;Ok(Some((b-a)/300.0))
        }
    }
    static PLATEAU:Plateau=Plateau;
    let (old,route)=fixture(&PLATEAU,[0.0,0.0,10.0],K,[0.0;3]);
    let step=ThirdBodyHeatStep::new(&old,0.5,&route,||Ok(())).unwrap();
    let (next,(),r)=step.evaluate_with_entropy(&[0.0,0.0,0.8],
        |_,_|Ok::<_,DepositionError>(source(0.8,K,[0.0;3])),||Ok(()))
        .unwrap().accept(controls(),||Ok(())).unwrap();
    let e=r.entropy.unwrap();
    assert_eq!(next.temperature_k(LAYER).unwrap(),300.0);
    assert_eq!(r.body_heat_j,[0.0;2]); assert_eq!(e.network_conduction_generation_j_k,0.0);
    close(e.thermal_storage_change_j_k,0.8/300.0,1e-15);
    close(e.friction_generation_j_k,0.8/300.0,1e-15);
    close(e.temporal_defect_j_k,0.0,1e-15); balance(e);
}

#[test]
fn entropy_refusal_and_energy_gates_preserve_both_physical_histories() {
    struct Missing;
    impl ThermalLaw for Missing {
        fn temperature_k(&self,h:f64)->Result<f64,DepositionError>{SENSIBLE.temperature_k(h)}
        fn maximum_temperature_slope(&self)->f64{1.0}
    }
    static MISSING:Missing=Missing;
    let (old,route)=fixture(&MISSING,[0.0;3],K,[0.0;3]);
    let saved=old.cells().clone(); let energy=endpoint(&old,0.1,0.4,K,[0.0;3]);
    let step=ThirdBodyHeatStep::new(&old,0.1,&route,||Ok(())).unwrap();
    assert!(step.evaluate(&energy,|_,_|Ok::<_,DepositionError>(source(0.4,K,[0.0;3])),||Ok(())).is_ok());
    assert!(matches!(step.evaluate_with_entropy(&energy,
        |_,_|Ok::<_,DepositionError>(source(0.4,K,[0.0;3])),||Ok(())),
        Err(DepositionError::InvalidInput(_))));
    let (supported,route)=fixture(&SENSIBLE,[0.0;3],K,[0.0;3]);
    let step=ThirdBodyHeatStep::new(&supported,0.1,&route,||Ok(())).unwrap();
    let candidate=step.evaluate_with_entropy(&energy,
        |_,_|Ok::<_,DepositionError>(source(0.4,K,[0.0;3])),||Ok(())).unwrap();
    assert!(matches!(candidate.accept(controls(),||Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)));
    let defective=step.evaluate_with_entropy(&[0.0;3],
        |_,_|Ok::<_,DepositionError>(source(0.4,K,[0.0;3])),||Ok(())).unwrap();
    assert!(defective.accept(controls(),||Ok(())).is_err());
    assert_eq!(old.cells(),&saved);assert_eq!(supported.cells(),&saved);
    assert_eq!(old.time_s(),0.0);assert_eq!(supported.time_s(),0.0);
}

#[test]
fn time_refinement_removes_numerical_entropy_without_changing_exact_storage() {
    let defects:Vec<f64>=[10,20,40].into_iter().map(|n| {
        let (mut state,route)=fixture(&SENSIBLE,[0.0;3],[0.0;3],[0.0;3]);
        let dt=1.0/f64::from(n);let mut defect=0.0;let mut storage=0.0;
        for _ in 0..n {
            let step=ThirdBodyHeatStep::new(&state,state.time_s()+dt,&route,||Ok(())).unwrap();
            let (next,(),r)=step.evaluate_with_entropy(&[0.0,0.0,0.8*dt],
                |_,_|Ok::<_,DepositionError>(source(0.8*dt,[0.0;3],[0.0;3])),||Ok(()))
                .unwrap().accept(controls(),||Ok(())).unwrap();
            let e=r.entropy.unwrap();defect+=e.temporal_defect_j_k;
            storage+=e.thermal_storage_change_j_k;state=next;
        }
        close(storage,0.25*(303.2_f64/300.0).ln(),2e-13);defect
    }).collect();
    for w in defects.windows(2) {assert!((1.9..2.1).contains(&(w[0]/w[1])),"{defects:?}");}
}
