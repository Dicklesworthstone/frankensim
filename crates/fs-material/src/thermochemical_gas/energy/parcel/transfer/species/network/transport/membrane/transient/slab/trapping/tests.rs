//! Independent kinetic/conservation limits and physical breakthrough history.
use super::*;

const AREA: f64 = 0.01;
const L: f64 = 2e-4;
const D: f64 = 2e-10;
const S: f64 = 1e-5;
const P: f64 = 80000.0;

fn layer(cells: usize) -> MembraneDiffusionLayer {
    MembraneDiffusionLayer { thickness_m: L, cells,
        diffusivity: MembraneDiffusivity::try_new(D).unwrap(),
        solubility: MolarHenrySolubility::try_new(S).unwrap() }
}
fn mobile(cells: usize, pressure: f64, sealed: bool) -> MembraneDiffusionState {
    let seal = sealed.then(|| MolarPermeance::try_new(0.0).unwrap());
    MembraneDiffusion::new_with_surfaces(AREA,&[layer(cells)],&[],[seal;2],||false).unwrap()
        .initial_state(&vec![pressure;cells],||false).unwrap()
}
fn controls() -> TrappingMembraneSolve {
    TrappingMembraneSolve { balance: MembraneDiffusionSolve {
        absolute_amount_tolerance_mol: 1e-19, relative_amount_tolerance: 1e-11 },
        ..Default::default() }
}
fn trapped(cells: usize, pressure: f64, capacity: f64, affinity: f64, rate: f64,
    loading: f64, sealed: bool) -> TrappingMembraneDiffusionState {
    TrappingMembraneDiffusionState::new(mobile(cells,pressure,sealed),
        &[Some(LangmuirMembraneTrap::try_new(capacity,affinity,rate).unwrap())],
        &vec![loading;cells],||false).unwrap()
}
fn near(a: f64, b: f64, absolute: f64) { assert!((a-b).abs() <= absolute,"{a} != {b}"); }

#[test]
fn sealed_cell_coupled_uptake_matches_independent_quadratic_and_retains_all_amount() {
    let state = trapped(1,P,2.0,2e-5,0.5,0.2e-6,true);
    let mut work = state.workspace().unwrap(); let dt = 10.0;
    let result = state.step_candidate([0.0;2],dt,controls(),&mut work,||false).unwrap();
    let cm = AREA*L*S; let ct = AREA*L*2.0;
    let total = cm*P+0.2e-6;
    let a = dt*0.5*2e-5; let d = 1.0+dt*0.5;
    let aa = cm*a; let bb = cm*d+a*(ct-total); let cc = 0.2e-6-total*d;
    let exact = 2.0*(-cc)/(bb+(bb*bb-4.0*aa*cc).sqrt());
    near(result.state().mobile().partial_pressures_pa()[0],exact,1e-7);
    near(result.state().stored_amount_mol().unwrap(),total,1e-18);
    assert_eq!(result.step().boundary_transfers_mol,[0.0;2]);
    assert!(result.step().trapped_change_mol>0.0);
    assert!(result.step().mobile_change_mol<0.0);
    assert!(result.step().nonlinear_iterations>1);
    near(result.step().gross_adsorption_mol[0]-result.step().gross_desorption_mol[0],
        result.step().trapped_change_mol,1e-19);
}

#[test]
fn no_traps_and_inert_traps_recover_original_henry_transport() {
    let base = mobile(12,P*0.25,false);
    let mut original_work = MembraneDiffusionWorkspace::new(12).unwrap();
    let reference = base.step_candidate([P,0.0],5.0,controls().balance,&mut original_work,||false).unwrap();
    for trap in [None,Some(LangmuirMembraneTrap::try_new(4.0,2e-5,0.0).unwrap())] {
        let initial = if trap.is_some() { 1e-8 } else { 0.0 };
        let state = TrappingMembraneDiffusionState::new(base.clone(),&[trap],&[initial;12],||false).unwrap();
        let mut work = state.workspace().unwrap();
        let result = state.step_candidate([P,0.0],5.0,controls(),&mut work,||false).unwrap();
        for (&a,&b) in result.state().mobile().partial_pressures_pa().iter().zip(reference.state().partial_pressures_pa()) {
            near(a,b,1e-8);
        }
        assert_eq!(result.state().trapped_moles(),&[initial;12]);
        assert_eq!(result.step().gross_adsorption_mol,vec![0.0;12]);
        assert_eq!(result.step().gross_desorption_mol,vec![0.0;12]);
    }
}

#[test]
fn a_sub_ulp_trap_change_cannot_silently_erase_physical_uptake() {
    // Only 1e-8 mol is mobile beside a 1e12 mol preloaded population.
    // Its initial kinetic imbalance is 1e-6 mol/s, but the bound-state ULP
    // exceeds all available mobile material. Unchanged stored amounts would
    // appear conservative while violating the original finite-rate site law.
    let state=trapped(1,500.0,1e18,0.004,1e-18,1e12,true);
    let saved=state.clone();
    let mut work=state.workspace().unwrap();
    assert!(matches!(state.step_candidate([0.0;2],1.0,controls(),&mut work,||false),
        Err(TrappingMembraneError::Diffusion(MembraneDiffusionError::AmountBalance {
            residual_mol,tolerance_mol,
        })) if residual_mol > 1e-7 && tolerance_mol < 1e-18));
    assert_eq!(state,saved);
}

#[test]
fn finite_reservoirs_relax_to_the_nonlinear_equilibrium_with_traps_in_the_inventory() {
    let state = trapped(1,0.0,2.0,2e-5,0.5,0.0,false);
    let capacity = MolarPressureCapacity::try_new(AREA*L*S).unwrap();
    let mut state = FiniteTrappingMembraneDiffusionState::new(state,[capacity;2],[P,0.0]).unwrap();
    let initial_total = state.total_amount_mol().unwrap();
    let cm = 3.0*capacity.value_si(); let ct = AREA*L*2.0; let b = 2e-5;
    let bb = cm+b*(ct-initial_total);
    let exact = 2.0*initial_total/(bb+(bb*bb+4.0*cm*b*initial_total).sqrt());
    let mut work = state.workspace().unwrap();
    for _ in 0..256 {
        let next = state.step_candidate(20.0,controls(),&mut work,||false).unwrap();
        let balance = next.step().finite_reservoirs.as_ref().unwrap();
        assert!(balance.closed_amount_residual_mol.abs()<1e-18);
        state = next.into_parts().0;
    }
    for q in state.reservoir_pressures_pa().into_iter().chain(state.membrane().mobile().partial_pressures_pa().iter().copied()) {
        near(q,exact,0.01);
    }
    near(state.total_amount_mol().unwrap(),initial_total,1e-16);
    assert!(state.membrane().trapped_moles()[0]>0.0);
}

#[test]
fn saturable_traps_delay_breakthrough_and_preloading_changes_the_transient() {
    let n = 24; let dt = 0.4;
    let mut bare = mobile(n,0.0,false);
    let mut bare_work = MembraneDiffusionWorkspace::new(n).unwrap();
    let mut empty = trapped(n,0.0,4.0,2e-5,0.2,0.0,false);
    let cell_capacity = 4.0*AREA*L/n as f64;
    let loaded = cell_capacity*(2e-5*P)/(1.0+2e-5*P);
    let mut preloaded = trapped(n,P,4.0,2e-5,0.2,loaded,false);
    let mut empty_work = empty.workspace().unwrap(); let mut loaded_work = preloaded.workspace().unwrap();
    let mut transmitted = [0.0;3];
    let mut unloaded = 0.0;
    for _ in 0..100 {
        let next = bare.step_candidate([P,0.0],dt,controls().balance,&mut bare_work,||false).unwrap();
        transmitted[0] += next.step().membrane_to_right_reservoir_mol;
        bare = next.into_parts().0;
        let next = empty.step_candidate([P,0.0],dt,controls(),&mut empty_work,||false).unwrap();
        transmitted[1] += next.step().boundary_transfers_mol[1];
        empty = next.into_parts().0;
        let next = preloaded.step_candidate([P,0.0],dt,controls(),&mut loaded_work,||false).unwrap();
        transmitted[2] += next.step().boundary_transfers_mol[1];
        unloaded += next.step().trapped_change_mol;
        preloaded = next.into_parts().0;
    }
    assert!(transmitted[1]<0.25*transmitted[0],"{transmitted:?}");
    assert!(transmitted[2]>transmitted[0],"{transmitted:?}");
    assert!(unloaded<0.0);
    assert!(empty.trapped_moles().iter().all(|&amount|amount>=0.0&&amount<=cell_capacity));
}

#[test]
fn separate_sealed_kinetic_equilibria_preserve_both_stores_exactly() {
    let seal = Some(MolarPermeance::try_new(0.0).unwrap());
    let model = MembraneDiffusion::new_with_surfaces(AREA,&[layer(3),layer(4)],&[seal],[seal;2],||false).unwrap();
    let q = [P,P,P,0.5*P,0.5*P,0.5*P,0.5*P];
    let mobile = model.initial_state(&q,||false).unwrap();
    let mut initial = vec![0.5*(2.0*(L/3.0*AREA));3];
    initial.extend([0.5*(7.0*(L/4.0*AREA));4]);
    let state = TrappingMembraneDiffusionState::new(mobile,&[
        Some(LangmuirMembraneTrap::try_new(2.0,1.0/P,0.5).unwrap()),
        Some(LangmuirMembraneTrap::try_new(7.0,1.0/(0.5*P),0.5).unwrap())],&initial,||false).unwrap();
    let mut work = state.workspace().unwrap();
    let result = state.step_candidate([2.0*P,0.0],10.0,controls(),&mut work,||false).unwrap();
    assert_eq!(result.state().mobile().partial_pressures_pa(),&q);
    assert_eq!(result.state().trapped_moles(),initial);
    assert_eq!(result.step().transport.face_molar_fluxes_mol_m2_s,vec![0.0;8]);
    assert_eq!(result.step().nonlinear_iterations,0);
    assert_eq!(result.step().trapped_change_mol,0.0);
    assert!(result.step().gross_adsorption_mol.iter().all(|&x|x>0.0));
}

#[test]
fn nonlinear_work_refusal_and_late_cancellation_leave_both_inventories_unchanged() {
    let state = trapped(8,0.0,4.0,2e-5,0.2,0.0,false); let saved = state.clone();
    let mut work = state.workspace().unwrap();
    let limited = TrappingMembraneSolve { maximum_iterations:1,..controls() };
    assert!(matches!(state.step_candidate([P,0.0],20.0,limited,&mut work,||false),
        Err(TrappingMembraneError::NotConverged { .. })));
    assert_eq!(state,saved);
    let mut polls = 0;
    let reference = state.step_candidate([P,0.0],20.0,controls(),&mut work,
        ||{polls+=1;false}).unwrap();
    for stop in [1,polls/2,polls] {
        let mut calls = 0;
        assert!(matches!(state.step_candidate([P,0.0],20.0,controls(),&mut work,
            ||{calls+=1;calls==stop}),Err(TrappingMembraneError::Cancelled)));
        assert_eq!(state,saved);
    }
    assert_eq!(reference,state.step_candidate([P,0.0],20.0,controls(),&mut work,||false).unwrap());
    let mut fresh = state.workspace().unwrap();
    assert_eq!(reference,state.step_candidate([P,0.0],20.0,controls(),&mut fresh,||false).unwrap());
}
