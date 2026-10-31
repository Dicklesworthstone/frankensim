//! Independent capacity/eigenmode oracles for the closed amount balance.
use super::*;
use fs_math::det;

const AREA:f64=0.01;
const L:f64=2e-4;
const D:f64=2e-10;
const S:f64=1e-5;
const P:f64=80000.0;

fn layer(length:f64,cells:usize,d:f64,s:f64)->MembraneDiffusionLayer {
    MembraneDiffusionLayer {thickness_m:length,cells,diffusivity:MembraneDiffusivity::try_new(d).unwrap(),
        solubility:MolarHenrySolubility::try_new(s).unwrap()}
}
fn capacity(value:f64)->MolarPressureCapacity {MolarPressureCapacity::try_new(value).unwrap()}
fn controls()->MembraneDiffusionSolve {
    MembraneDiffusionSolve {absolute_amount_tolerance_mol:1e-14,relative_amount_tolerance:1e-10}
}
fn membrane(n:usize)->MembraneDiffusionState {
    MembraneDiffusion::new(AREA,&[layer(L,n,D,S)],&[],||false).unwrap()
        .initial_state(&vec![0.0;n],||false).unwrap()
}
fn symmetric()->FiniteMembraneDiffusionState {
    FiniteMembraneDiffusionState::new(membrane(1),[capacity(AREA*S*L);2],[P,0.0],||false).unwrap()
}
fn near(a:f64,b:f64,tolerance:f64) {
    assert!((a-b).abs()<=tolerance*b.abs().max(1e-30),"{a} != {b}, tolerance {tolerance}");
}

#[test]
fn symmetric_three_inventory_transient_matches_independent_eigenmodes() {
    // With equal capacities and one cell, the chain has rates 0, lambda and
    // 3*lambda, with eigenvectors [1,1,1], [1,0,-1], [1,-2,1].
    let lambda=2.0*D/(L*L);let dt=2.0;
    let mut state=symmetric();let total=state.total_amount_mol().unwrap();
    let mut work=FiniteMembraneDiffusionWorkspace::new(1).unwrap();
    for k in 1..=80 {
        let candidate=state.step_candidate(dt,controls(),&mut work,||false).unwrap();
        assert!(candidate.step().reservoir_amount_changes_mol[0]<0.0);
        assert!(candidate.step().reservoir_amount_changes_mol[1]>0.0);
        assert!(candidate.step().amount_residual_mol.abs()<1e-14);
        let (next,_)=candidate.into_parts();state=next;
        let slow=(1.0+lambda*dt).powi(-k);let fast=(1.0+3.0*lambda*dt).powi(-k);
        let [left,right]=state.reservoir_pressures_pa();
        near(left,P*(1.0/3.0+0.5*slow+fast/6.0),3e-12);
        near(state.membrane().partial_pressures_pa()[0],P*(1.0-fast)/3.0,3e-12);
        near(right,P*(1.0/3.0-0.5*slow+fast/6.0),3e-11);
        near(state.total_amount_mol().unwrap(),total,3e-13);
    }
    // Finite initial storage delays transmission while the inlet depletes.
    let first=symmetric().step_candidate(0.1,controls(),&mut work,||false).unwrap();
    assert!(first.step().membrane.left_reservoir_to_membrane_mol
        >first.step().membrane.membrane_to_right_reservoir_mol*100.0);
    assert!(first.state().reservoir_pressures_pa()[0]<P);
}

#[test]
fn finite_reservoir_transient_converges_in_time_to_the_continuous_capacity_model() {
    let lambda=2.0*D/(L*L);let time=0.8/lambda;
    let exact=P*(1.0/3.0-0.5*det::exp(-lambda*time)+det::exp(-3.0*lambda*time)/6.0);
    let mut errors=Vec::new();
    for steps in [20,40,80] {
        let mut state=symmetric();let mut work=FiniteMembraneDiffusionWorkspace::new(1).unwrap();
        for _ in 0..steps {
            state=state.step_candidate(time/f64::from(steps),controls(),&mut work,||false).unwrap().into_parts().0;
        }
        errors.push((state.reservoir_pressures_pa()[1]-exact).abs());
    }
    for pair in errors.windows(2) {assert!((1.8..2.2).contains(&(pair[0]/pair[1])),"{errors:?}");}
}

#[test]
fn unequal_reservoirs_and_layers_equilibrate_with_sorbed_capacity_in_the_denominator() {
    let layers=[layer(L,8,D,S),layer(0.5*L,4,0.5*D,4.0*S)];
    let model=MembraneDiffusion::new(AREA,&layers,&[Some(MolarPermeance::try_new(1e-11).unwrap())],||false).unwrap();
    let film=model.initial_state(&[0.0;12],||false).unwrap();
    let capacities=[capacity(5e-12),capacity(2e-10)];
    let mut state=FiniteMembraneDiffusionState::new(film,capacities,[P,2000.0],||false).unwrap();
    let total=state.total_amount_mol().unwrap();
    let equilibrium=total/(5e-12+2e-10+AREA*(S*L+4.0*S*0.5*L));
    let mut work=FiniteMembraneDiffusionWorkspace::new(12).unwrap();
    for _ in 0..60 {
        let candidate=state.step_candidate(2000.0,controls(),&mut work,||false).unwrap();
        assert!(candidate.state().reservoir_pressures_pa().into_iter().all(|p|p>=0.0&&p<=P));
        assert!(candidate.state().membrane().partial_pressures_pa().iter().all(|&p|p>=0.0&&p<=P));
        near(candidate.state().total_amount_mol().unwrap(),total,2e-10);
        state=candidate.into_parts().0;
    }
    for p in state.reservoir_pressures_pa().into_iter().chain(state.membrane().partial_pressures_pa().iter().copied()) {
        near(p,equilibrium,2e-9);
    }
    // A six-order capacity ratio remains nonnegative without inventory clips.
    let state=FiniteMembraneDiffusionState::new(membrane(16),[capacity(1e-14),capacity(1e-8)],[P,0.0],||false).unwrap();
    let mut work=FiniteMembraneDiffusionWorkspace::new(16).unwrap();
    let next=state.step_candidate(1e4,controls(),&mut work,||false).unwrap();
    assert!(next.state().reservoir_pressures_pa().into_iter().all(|p|p>=0.0&&p<=P));
    near(next.state().total_amount_mol().unwrap(),state.total_amount_mol().unwrap(),2e-8);
}

#[test]
fn seals_preserve_separate_connected_inventory_and_authored_zero_diffusion() {
    let layers=[layer(L,2,D,S),layer(L,3,D,S),layer(L,2,D,S)];
    let seal=Some(MolarPermeance::try_new(0.0).unwrap());
    let film=MembraneDiffusion::new(AREA,&layers,&[seal;2],||false).unwrap()
        .initial_state(&[0.0,0.0,1000.0,3000.0,7000.0,5000.0,5000.0],||false).unwrap();
    let mut state=FiniteMembraneDiffusionState::new(film,[capacity(AREA*S*L),capacity(3.0*AREA*S*L)],
        [P,10000.0],||false).unwrap();
    let total=state.total_amount_mol().unwrap();let mut work=FiniteMembraneDiffusionWorkspace::new(7).unwrap();
    for _ in 0..20 {
        let next=state.step_candidate(1000.0,controls(),&mut work,||false).unwrap();
        assert_eq!(next.step().membrane.face_molar_fluxes_mol_m2_s[2],0.0);
        assert_eq!(next.step().membrane.face_molar_fluxes_mol_m2_s[5],0.0);
        state=next.into_parts().0;
    }
    near(state.reservoir_pressures_pa()[0],0.5*P,1e-10);
    near(state.reservoir_pressures_pa()[1],8750.0,1e-10);
    for &p in &state.membrane().partial_pressures_pa()[2..5] {near(p,11000.0/3.0,1e-10);}
    near(state.total_amount_mol().unwrap(),total,1e-11);
    let film=MembraneDiffusion::new(AREA,&[layer(L,4,0.0,S)],&[],||false).unwrap()
        .initial_state(&[1000.0,3000.0,7000.0,9000.0],||false).unwrap();
    let sealed=FiniteMembraneDiffusionState::new(film,[capacity(1e-12),capacity(1e-8)],[P,10000.0],||false).unwrap();
    let mut work=FiniteMembraneDiffusionWorkspace::new(4).unwrap();
    let next=sealed.step_candidate(1e8,controls(),&mut work,||false).unwrap();
    assert_eq!(next.state().reservoir_pressures_pa(),sealed.reservoir_pressures_pa());
    assert_eq!(next.state().membrane().partial_pressures_pa(),sealed.membrane().partial_pressures_pa());
    assert_eq!(next.step().membrane.face_molar_fluxes_mol_m2_s,[0.0;5]);
    assert_eq!(next.step().reservoir_amount_changes_mol,[0.0;2]);
}

#[test]
fn large_capacities_approach_prescribed_reservoir_breakthrough() {
    let film=membrane(32);let pressures=[P,1000.0];
    let mut prescribed_work=MembraneDiffusionWorkspace::new(32).unwrap();
    let prescribed=film.step_candidate(pressures,1.0,controls(),&mut prescribed_work,||false).unwrap();
    let finite=FiniteMembraneDiffusionState::new(film,[capacity(1e-4);2],pressures,||false).unwrap();
    let mut work=FiniteMembraneDiffusionWorkspace::new(32).unwrap();
    let finite=finite.step_candidate(1.0,controls(),&mut work,||false).unwrap();
    for (&a,&b) in finite.state().membrane().partial_pressures_pa().iter().zip(prescribed.state().partial_pressures_pa()) {
        near(a,b,1e-6);
    }
    assert!(finite.step().reservoir_pressures_after_pa[0]<P);
    assert!(finite.step().amount_residual_mol.abs()<1e-14);
}

#[test]
fn capacity_cache_identity_and_cancellation_never_publish_partial_inventory() {
    assert_eq!(MOLAR_PRESSURE_CAPACITY_DIMS.checked_plus(fs_qty::Pressure::DIMS),Some(Dims([0,0,0,0,0,1])));
    for value in [0.0,-1.0,f64::INFINITY,f64::NAN] {assert!(MolarPressureCapacity::try_new(value).is_err());}
    assert!(FiniteMembraneDiffusionWorkspace::new(0).is_err());
    assert!(FiniteMembraneDiffusionWorkspace::new(MAX_TRANSIENT_MEMBRANE_CELLS+1).is_err());
    let film=membrane(256);
    let a=FiniteMembraneDiffusionState::new(film.clone(),[capacity(1e-11);2],[P,0.0],||false).unwrap();
    let b=FiniteMembraneDiffusionState::new(film,[capacity(3e-11),capacity(1e-11)],[P,0.0],||false).unwrap();
    let saved=a.clone();let mut reused=FiniteMembraneDiffusionWorkspace::new(256).unwrap();let mut polls=0;
    assert!(matches!(a.step_candidate(1.0,controls(),&mut reused,||{polls+=1;polls==8}),Err(Error::Cancelled)));
    assert_eq!(a,saved);
    a.step_candidate(1.0,controls(),&mut reused,||false).unwrap();
    let from_reuse=b.step_candidate(1.0,controls(),&mut reused,||false).unwrap();
    let mut fresh=FiniteMembraneDiffusionWorkspace::new(256).unwrap();
    assert_eq!(from_reuse,b.step_candidate(1.0,controls(),&mut fresh,||false).unwrap());
    let changed_dt=b.step_candidate(2.0,controls(),&mut reused,||false).unwrap();
    assert_eq!(changed_dt,b.step_candidate(2.0,controls(),&mut fresh,||false).unwrap());
    assert!(a.step_candidate(0.0,controls(),&mut reused,||false).is_err());
    assert_eq!(a,saved);
}
