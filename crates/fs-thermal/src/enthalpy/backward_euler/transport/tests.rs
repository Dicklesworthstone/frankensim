//! G1/G3: independent upwind balances, caloric datum and native energy/entropy.
use super::*;
use super::super::super::advection::{MassBoundary, MassFlow};

struct Law { datum: f64 }
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h - self.datum;
        if !t.is_finite() || !(100.0..=1000.0).contains(&t) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(t)
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
    fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64)
        -> Result<Option<f64>, DepositionError> {
        Ok(Some((self.temperature_k(b)? / self.temperature_k(a)?).ln()))
    }
}
fn thermal<'a>(law: &'a Law, values: &[(f64, f64)], max_links: usize)
    -> DepositionSimulation<'a> {
    let events = values.iter().enumerate().map(|(i, &(mass, h))| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: 0, body_id: 0, material_id: MaterialId(1),
            mass_kg: mass, specific_h_j_kg: law.datum + h,
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0],
                volume_m3: mass, surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0,
                reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    DepositionSimulation::new(events, &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 4, max_links, max_events: 4, max_substeps: 100 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-9 }, || Ok(())).unwrap().0
}
fn flow(a: u64, b: u64, id: u64, rate: f64) -> MassFlow {
    MassFlow { a: CellId(a), b: CellId(b), id, kg_s: rate }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 2e-10, "{a:e} != {b:e}"); }

#[test]
fn closed_cycle_uses_exact_internal_transfers_above_explicit_cfl() {
    let law = Law { datum: 0.0 };
    let old = thermal(&law, &[(1.0,0.0),(1.0,40.0),(1.0,100.0)], 3);
    let graph = Advected::new((), vec![flow(0,1,8,2.0),flow(1,2,2,2.0),flow(2,0,4,2.0)],3).unwrap();
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    // (I + dt*A) h_new = h_old for a directed 3-cycle at rate 2 kg/s.
    // Independent exact solution is [40,40,60]; no Fourier symmetry assumed.
    let trial = step.evaluate_with_transport_and_entropy(&[40.0,0.0,-40.0],
        |state| graph.endpoint_transport(state, || Ok(())), || Ok(())).unwrap();
    near(trial.maximum_cell_residual_j(),0.0);
    let entropy = *trial.entropy_audit().unwrap();
    near(entropy.unclassified_source_at_cell_j_k,40.0/340.0-40.0/360.0);
    assert!(entropy.has_unclassified_sources);
    assert_eq!(entropy.conduction_generation_j_k,0.0);
    near(entropy.thermal_storage_change_j_k,
        (340.0_f64/300.0).ln()+(360.0_f64/400.0).ln());
    let (accepted,(),audit) = trial.accept(1e-9,1e-9,||Ok(())).unwrap();
    assert_eq!(audit.external_heat_j,0.0);
    assert_eq!(audit.incoming_mass_kg,0.0);
    for (i,t) in [340.0,340.0,360.0].into_iter().enumerate() {
        near(accepted.temperature_k(CellId(i as u64)).unwrap(),t);
        assert_eq!(accepted.cells()[&CellId(i as u64)].mass_kg,1.0);
    }
    assert_eq!(graph.audit(),Default::default());
    assert_eq!(old.time_s(),0.0);
    assert_eq!(accepted.time_s(),1.0);
}

#[test]
fn inlet_outlet_use_current_donor_enthalpy_and_preserve_large_common_datum() {
    for datum in [0.0,1099511627776.0] {
        let law = Law { datum };
        let old = thermal(&law,&[(1.0,0.0),(1.0,0.0)],1);
        let graph = Advected::new((),vec![flow(0,1,2,1.0)],3).unwrap().with_boundaries(vec![
            MassBoundary::Inlet {id:1,cell:CellId(0),kg_s:1.0,material_id:MaterialId(1),specific_h_j_kg:datum+120.0},
            MassBoundary::Outlet {id:3,cell:CellId(1),kg_s:1.0},
        ],3).unwrap();
        let step = BackwardEulerStep::new(&old,1.0,||Ok(())).unwrap();
        // h0=120/2=60; h1=h0/2=30; exported heat is 120-h1=90 J.
        let trial = step.evaluate_with_transport_and_entropy(&[60.0,30.0],
            |state| graph.endpoint_transport(state,||Ok(())),||Ok(())).unwrap();
        near(trial.maximum_cell_residual_j(),0.0);
        let (next,(),audit) = trial.accept(1e-9,1e-9,||Ok(())).unwrap();
        near(audit.external_heat_j,90.0);
        near(audit.stored_energy_change_j,90.0);
        near(next.temperature_k(CellId(0)).unwrap(),360.0);
        near(next.temperature_k(CellId(1)).unwrap(),330.0);
        assert_eq!(next.cells()[&CellId(1)].boundary.source_w,0.0);
    }
}

#[test]
fn empty_transport_is_the_original_implicit_boundary_and_entropy_path() {
    let law = Law { datum: 0.0 };
    let old = thermal(&law,&[(2.0,0.0)],1);
    let step = BackwardEulerStep::new(&old,0.5,||Ok(())).unwrap();
    let coefficients = || BoundaryStepCoefficients {
        bulk: StepCoefficients {sources_w:[(CellId(0),6.0)].into(),
            conductances_w_k:BTreeMap::new(),payload:17usize},
        boundaries:BoundaryConductances::default(),
    };
    let a = step.evaluate_with_boundary_conductances_and_entropy(&[3.0],
        |_| Ok::<_,DepositionError>(coefficients()),||Ok(())).unwrap();
    let b = step.evaluate_with_transport_and_entropy(&[3.0],
        |_| Ok::<_,DepositionError>(TransportStepCoefficients::from(coefficients())),||Ok(())).unwrap();
    assert_eq!(a.residual_j(),b.residual_j());
    assert_eq!(a.entropy_audit(),b.entropy_audit());
    let (a,pa,aa)=a.accept(0.0,0.0,||Ok(())).unwrap();
    let (b,pb,ab)=b.accept(0.0,0.0,||Ok(())).unwrap();
    assert_eq!(a.cells(),b.cells());assert_eq!(pa,pb);
    assert_eq!(aa.external_heat_j,ab.external_heat_j);
    assert_eq!(a.cells()[&CellId(0)].boundary.source_w,0.0);
}

#[test]
fn invalid_mass_coverage_transfer_budget_and_late_cancellation_never_publish() {
    let law = Law { datum: 0.0 };
    let old = thermal(&law,&[(1.0,30.0),(1.0,0.0)],1);
    let saved=old.cells().clone();
    let unbalanced=Advected::new((),vec![flow(0,1,1,1.0)],1).unwrap();
    assert!(unbalanced.endpoint_transport(&old,||Ok(())).is_err());
    let bad_material=unbalanced.clone().with_boundaries(vec![
        MassBoundary::Inlet {id:2,cell:CellId(0),kg_s:1.0,material_id:MaterialId(2),specific_h_j_kg:30.0},
        MassBoundary::Outlet {id:3,cell:CellId(1),kg_s:1.0}],3).unwrap();
    assert!(bad_material.endpoint_transport(&old,||Ok(())).is_err());
    let cycle=Advected::new((),vec![flow(0,1,1,1.0),flow(1,0,2,1.0)],2).unwrap();
    let step=BackwardEulerStep::new(&old,1.0,||Ok(())).unwrap();
    assert!(step.evaluate_with_transport(&[-10.0,10.0],
        |state| cycle.endpoint_transport(state,||Ok(())),||Ok(())).is_err()); // original max_links=1
    let transfer=|a,b,power_w| TransportStepCoefficients {
        heat:StepCoefficients {sources_w:BTreeMap::new(),conductances_w_k:BTreeMap::new(),payload:()}.into(),
        transfers:vec![HeatTransfer {a:CellId(a),b:CellId(b),power_w}],
    };
    for (a,b,power) in [(0,0,1.0),(0,9,1.0),(0,1,f64::NAN)] {
        assert!(step.evaluate_with_transport(&[0.0;2],
            |_|Ok::<_,DepositionError>(transfer(a,b,power)),||Ok(())).is_err());
    }
    let mut checks=0;
    step.evaluate_with_transport(&[-10.0,10.0],
        |_|Ok::<_,DepositionError>(transfer(0,1,10.0)),||{checks+=1;Ok(())}).unwrap();
    let mut seen=0;
    assert!(matches!(step.evaluate_with_transport(&[-10.0,10.0],
        |_|Ok::<_,DepositionError>(transfer(0,1,10.0)),||{
            seen+=1;if seen==checks {Err(DepositionError::Cancelled)} else {Ok(())}
        }),Err(DepositionError::Cancelled)));
    assert_eq!(old.cells(),&saved);assert_eq!(old.time_s(),0.0);
    let retry=step.evaluate_with_transport(&[-10.0,10.0],
        |_|Ok::<_,DepositionError>(transfer(0,1,10.0)),||Ok(())).unwrap();
    retry.accept(0.0,0.0,||Ok(())).unwrap();
    // A zero total defect cannot hide a missing pair of local debits/credits.
    let wrong=step.evaluate_with_transport(&[0.0;2],
        |_|Ok::<_,DepositionError>(transfer(0,1,10.0)),||Ok(())).unwrap();
    assert_eq!(wrong.audit().residual_j,0.0);
    assert!(wrong.accept(1e-10,1e-10,||Ok(())).is_err());
}
