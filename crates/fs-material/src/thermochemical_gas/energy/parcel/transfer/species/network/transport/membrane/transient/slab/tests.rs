//! G0/G1/G3 analytical diffusion, partition, time-lag and inventory oracles.
use super::*;
use super::super::super::coefficients::{MEMBRANE_DIFFUSIVITY_DIMS,MOLAR_HENRY_SOLUBILITY_DIMS,
    MOLAR_PERMEABILITY_DIMS};
use fs_math::det;

const D:f64=2e-10;
const S:f64=1e-5;
const L:f64=2e-4;
const AREA:f64=0.01;
const PRESSURE:f64=80000.0;

fn layer(thickness_m:f64,cells:usize,d:f64,s:f64)->MembraneDiffusionLayer {
    MembraneDiffusionLayer {thickness_m,cells,diffusivity:MembraneDiffusivity::try_new(d).unwrap(),
        solubility:MolarHenrySolubility::try_new(s).unwrap()}
}
fn controls()->MembraneDiffusionSolve {
    MembraneDiffusionSolve {absolute_amount_tolerance_mol:1e-14,relative_amount_tolerance:1e-9}
}
fn initial(n:usize)->MembraneDiffusionState {
    MembraneDiffusion::new(AREA,&[layer(L,n,D,S)],&[],||false).unwrap()
        .initial_state(&vec![0.0;n],||false).unwrap()
}
fn relative(a:f64,b:f64,tolerance:f64) {
    assert!((a-b).abs()<=tolerance*b.abs(),"{a} versus {b}, relative tolerance {tolerance}");
}
fn analytic_outlet_fraction(time_s:f64)->f64 {
    let mut result=1.0;
    for n in 1..=100 {
        let sign=if n%2==0 {1.0}else{-1.0};
        result+=2.0*sign*det::exp(-(n*n) as f64*core::f64::consts::PI.powi(2)*D*time_s/(L*L));
    }
    result
}
fn analytic_stored_fraction(time_s:f64)->f64 {
    let mut result=0.5;
    for n in (1..=199).step_by(2) {
        result-=4.0/((n*n) as f64*core::f64::consts::PI.powi(2))
            *det::exp(-(n*n) as f64*core::f64::consts::PI.powi(2)*D*time_s/(L*L));
    }
    result
}
fn run(n:usize,dt:f64,steps:usize)->(MembraneDiffusionState,MembraneDiffusionStep,f64) {
    let mut state=initial(n);
    let mut work=MembraneDiffusionWorkspace::new(n).unwrap();
    let mut out=Sum::default();let mut last=None;
    for _ in 0..steps {
        let candidate=state.step_candidate([PRESSURE,0.0],dt,controls(),&mut work,||false).unwrap();
        let (next,report)=candidate.into_parts();
        out.add(report.membrane_to_right_reservoir_mol).unwrap();
        assert!(report.amount_residual_mol.abs()<1e-14);
        state=next;last=Some(report);
    }
    (state,last.unwrap(),out.value)
}

#[test]
fn typed_diffusivity_solubility_product_has_the_original_permeability_basis() {
    assert_eq!(MEMBRANE_DIFFUSIVITY_DIMS.checked_plus(MOLAR_HENRY_SOLUBILITY_DIMS),Some(MOLAR_PERMEABILITY_DIMS));
    let p=MolarPermeability::from_diffusion(MembraneDiffusivity::try_new(D).unwrap(),
        MolarHenrySolubility::try_new(S).unwrap()).unwrap();
    relative(p.value_si(),D*S,1e-15);
    relative(p.permeance(L).unwrap().value_si(),D*S/L,1e-15);
    for bad in [-1.0,f64::NAN,f64::INFINITY] {
        assert!(MembraneDiffusivity::try_new(bad).is_err());
        assert!(MolarHenrySolubility::try_new(bad).is_err());
    }
}

#[test]
fn degassed_slab_has_real_breakthrough_delay_and_matches_fourier_solution() {
    let flux=initial(96).sample_fluxes([PRESSURE,0.0],||false).unwrap();
    assert!(flux.face_molar_fluxes_mol_m2_s[0]>0.0);
    assert_eq!(flux.face_molar_fluxes_mol_m2_s[96],0.0);
    let tau=L*L/D;
    let (state,report,_)=run(96,tau/4000.0,1200);
    let steady=D*S*PRESSURE/L;
    relative(report.face_molar_fluxes_mol_m2_s[96]/steady,analytic_outlet_fraction(0.3*tau),0.001);
    relative(state.stored_amount_mol().unwrap()/(AREA*L*S*PRESSURE),analytic_stored_fraction(0.3*tau),0.0003);
    let (_,early,_)=run(96,tau/4000.0,40);
    assert!(early.face_molar_fluxes_mol_m2_s[96]<steady*1e-6);
    assert!(early.left_reservoir_to_membrane_mol>early.membrane_to_right_reservoir_mol);
    assert!(early.stored_after_mol>0.0);
}

#[test]
fn backward_euler_transient_error_converges_at_first_order_in_time() {
    let tau=L*L/D;
    let exact=analytic_outlet_fraction(0.2*tau);
    let mut errors=Vec::new();
    for divisions in [80,160,320] {
        let (_,step,_)=run(256,tau/f64::from(divisions),divisions as usize/5);
        errors.push((step.face_molar_fluxes_mol_m2_s[256]/(D*S*PRESSURE/L)-exact).abs());
    }
    for pair in errors.windows(2) {assert!((1.7..2.3).contains(&(pair[0]/pair[1])),"{errors:?}");}
}

#[test]
fn cumulative_permeation_time_lag_converges_to_l_squared_over_six_d() {
    let tau=L*L/D;let continuum=tau/6.0;let mut errors=Vec::new();
    for n in [8,16,32] {
        let (state,_,out)=run(n,tau/100.0,300);
        let lag=state.elapsed_s()-out/(AREA*D*S*PRESSURE/L);
        let finite_volume=continuum+(L/n as f64).powi(2)/(12.0*D);
        relative(lag,finite_volume,2e-9);
        errors.push((lag-continuum).abs());
    }
    for pair in errors.windows(2) {relative(pair[0]/pair[1],4.0,1e-5);}
}

#[test]
fn layered_storage_and_interface_resistance_approach_the_same_steady_series_flux() {
    let layers=[layer(1e-4,8,1e-10,2e-5),layer(2e-4,16,4e-10,1e-5)];
    let interface=MolarPermeance::try_new(2e-11).unwrap();
    let model=MembraneDiffusion::new(AREA,&layers,&[Some(interface)],||false).unwrap();
    let mut state=model.initial_state(&[0.0;24],||false).unwrap();
    let mut work=MembraneDiffusionWorkspace::new(24).unwrap();let mut last=None;
    for _ in 0..40 {
        let candidate=state.step_candidate([100000.0,0.0],1000.0,controls(),&mut work,||false).unwrap();
        let (next,report)=candidate.into_parts();state=next;last=Some(report);
    }
    let report=last.unwrap();
    let expected=100000.0/(1e-4/(1e-10*2e-5)+1.0/2e-11+2e-4/(4e-10*1e-5));
    for &flux in &report.face_molar_fluxes_mol_m2_s {relative(flux,expected,2e-11);}
    let [Some(left),Some(right)]=report.face_partial_pressures_pa[8] else {panic!("active interface missing traces")};
    relative((left-right)*interface.value_si(),expected,2e-11);
    relative(left,100000.0-expected*1e-4/(1e-10*2e-5),2e-11);
    relative(right,expected*2e-4/(4e-10*1e-5),2e-11);
}

#[test]
fn partition_equilibrium_preserves_common_potential_and_different_concentrations() {
    let layers=[layer(L,4,D,S),layer(L,4,D,3.0*S)];
    let state=MembraneDiffusion::new(AREA,&layers,&[None],||false).unwrap()
        .initial_state(&[PRESSURE;8],||false).unwrap();
    let mut work=MembraneDiffusionWorkspace::new(8).unwrap();
    let candidate=state.step_candidate([PRESSURE;2],100.0,controls(),&mut work,||false).unwrap();
    for &q in candidate.state().partial_pressures_pa() {relative(q,PRESSURE,2e-14);}
    let c=candidate.state().concentrations_mol_m3().collect::<Vec<_>>();
    relative(c[4]/c[3],3.0,2e-14);
    relative(candidate.state().stored_amount_mol().unwrap(),state.stored_amount_mol().unwrap(),2e-14);
}

#[test]
fn zero_diffusion_and_multiple_seals_preserve_authored_trapped_inventory() {
    let state=MembraneDiffusion::new(AREA,&[layer(L,4,0.0,S)],&[],||false).unwrap()
        .initial_state(&[1000.0,2000.0,5000.0,7000.0],||false).unwrap();
    let mut work=MembraneDiffusionWorkspace::new(4).unwrap();
    let candidate=state.step_candidate([0.0,PRESSURE],1e8,controls(),&mut work,||false).unwrap();
    assert_eq!(candidate.state().partial_pressures_pa(),state.partial_pressures_pa());
    assert_eq!(candidate.step().face_molar_fluxes_mol_m2_s,[0.0;5]);
    assert_eq!(candidate.step().stored_change_mol,0.0);
    let zero=MolarPermeance::try_new(0.0).unwrap();
    let layers=[layer(L,1,D,S);3];
    let state=MembraneDiffusion::new(AREA,&layers,&[Some(zero);2],||false).unwrap()
        .initial_state(&[0.0,12345.0,0.0],||false).unwrap();
    let mut work=MembraneDiffusionWorkspace::new(3).unwrap();
    let candidate=state.step_candidate([PRESSURE,0.0],100.0,controls(),&mut work,||false).unwrap();
    assert_eq!(candidate.state().partial_pressures_pa()[1],12345.0);
    assert_eq!(&candidate.step().face_molar_fluxes_mol_m2_s[1..3],&[0.0;2]);
    assert!(candidate.state().partial_pressures_pa()[0]>0.0);
}

#[test]
fn stiff_steps_remain_nonnegative_and_unloading_accounts_both_reservoirs() {
    let mut work=MembraneDiffusionWorkspace::new(32).unwrap();
    let saved=initial(32);
    let loaded=saved.step_candidate([PRESSURE,0.0],1e5,controls(),&mut work,||false).unwrap();
    assert!(loaded.state().partial_pressures_pa().iter().all(|&q|q>=0.0&&q<=PRESSURE));
    let unloaded=loaded.state().step_candidate([0.0;2],10.0,controls(),&mut work,||false).unwrap();
    assert!(unloaded.step().left_reservoir_to_membrane_mol<0.0);
    assert!(unloaded.step().membrane_to_right_reservoir_mol>0.0);
    assert!(unloaded.step().stored_change_mol<0.0);
    assert!(unloaded.step().amount_residual_mol.abs()<1e-14);
    assert_eq!(saved.stored_amount_mol().unwrap(),0.0);
}

#[test]
fn geometry_storage_budget_and_cancellation_refuse_without_a_partial_state() {
    assert!(MembraneDiffusion::new(AREA,&[],&[],||false).is_err());
    assert!(MembraneDiffusion::new(AREA,&[layer(L,0,D,S)],&[],||false).is_err());
    assert!(MembraneDiffusion::new(AREA,&[layer(L,1,D,0.0)],&[],||false).is_err());
    assert!(MembraneDiffusion::new(AREA,&[layer(L,MAX_TRANSIENT_MEMBRANE_CELLS+1,D,S)],&[],||false).is_err());
    let tiny=MembraneDiffusion::new(1e-200,&[layer(1e-200,1,0.0,1e200)],&[],||false).unwrap();
    assert!(tiny.initial_state(&[1e200],||false).is_err());
    let state=initial(256);let saved=state.clone();let mut work=MembraneDiffusionWorkspace::new(256).unwrap();
    let mut polls=0;
    assert!(matches!(state.step_candidate([PRESSURE,0.0],1.0,controls(),&mut work,||{polls+=1;polls==8}),Err(Error::Cancelled)));
    assert_eq!(state,saved);
    let a=state.step_candidate([PRESSURE,0.0],1.0,controls(),&mut work,||false).unwrap();
    let mut fresh=MembraneDiffusionWorkspace::new(256).unwrap();
    let b=state.step_candidate([PRESSURE,0.0],1.0,controls(),&mut fresh,||false).unwrap();
    assert_eq!(a,b);
    assert!(state.step_candidate([-1.0,0.0],1.0,controls(),&mut work,||false).is_err());
}
