use super::*;
use crate::rk45_adaptive_checked;
use std::cell::Cell;

struct Decay { rate: Cell<f64>, calls: Cell<usize>, derivatives: Cell<usize> }
impl Decay { fn new() -> Self { Self { rate: Cell::new(-0.7), calls: Cell::new(0), derivatives: Cell::new(0) } } }
impl OdeVjp for Decay {
    fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {1}
    fn rhs(&self,_:f64,u:&[f64],o:&mut[f64]) { self.calls.set(self.calls.get()+1);o[0]=self.rate.get()*u[0]; }
    fn rhs_vjp(&self,_:f64,u:&[f64],b:&[f64],x:&mut[f64],p:&mut[f64])->Result<(),String> {
        self.derivatives.set(self.derivatives.get()+1);x[0]=self.rate.get()*b[0];p[0]=u[0]*b[0];Ok(())
    }
}
fn config(end:f64)->RecordingConfig { RecordingConfig {end,rtol:1e-10,atol:1e-12,
    controller:PiController::default(),max_workspace_components:4096} }
fn budget()->ReplayBudget {ReplayBudget{checkpoints:64,replayed_steps:100_000}}
fn recording(model:&Decay)->RecordedRk45<'_,Decay> {
    RecordedRk45::new(model,AdaptiveState::new(0.0,&[1.2],1.0),config(2.0)).unwrap()
}
fn same(a:&AdaptiveState,b:&AdaptiveState) {
    assert_eq!((a.t.to_bits(),a.h.to_bits(),a.err_prev.to_bits(),a.accepted,a.rejected),
        (b.t.to_bits(),b.h.to_bits(),b.err_prev.to_bits(),b.accepted,b.rejected));
    assert_eq!(a.u.iter().map(|x|x.to_bits()).collect::<Vec<_>>(),b.u.iter().map(|x|x.to_bits()).collect::<Vec<_>>());
}

#[test]
fn recorded_forward_and_budgeted_clone_resume_match_production() {
    let model=Decay::new();let mut a=recording(&model);let mut b=recording(&model);
    assert_eq!(a.advance(10000,10000,&mut||false).unwrap().status,RecordingStatus::ReachedEnd);
    assert_eq!(b.advance(10000,1,&mut||false).unwrap().status,RecordingStatus::RecordLimit);
    assert_eq!(b.accepted_steps(),1);let mut resumed=b.clone();
    for _ in 0..10000 {
        if resumed.advance(1,10000,&mut||false).unwrap().status==RecordingStatus::ReachedEnd {break;}
    }
    same(a.state(),resumed.state());assert_eq!(a.records,resumed.records);assert!(a.state().rejected>0);
    let mut ordinary=AdaptiveState::new(0.0,&[1.2],1.0);
    rk45_adaptive_checked(&mut ordinary,&|t,u,o|model.rhs(t,u,o),2.0,1e-10,1e-12,&PiController::default(),10000).unwrap();
    same(a.state(),&ordinary);
}

#[test]
fn checkpointed_gradient_matches_full_storage_bitwise_and_exact_decay() {
    let model=Decay::new();let mut a=recording(&model);a.advance(10000,10000,&mut||false).unwrap();
    let got=a.pullback(&[0.8],&[0.3],budget(),&mut||false).unwrap();
    let mut states=vec![a.initial.clone()];
    for r in &a.records {states.push(replay(&model,states.last().unwrap(),r.start,r.end,r.h,&mut||false).unwrap().next);}
    let (mut x,mut p)=(vec![0.8],vec![0.3]);
    for (i,r) in a.records.iter().enumerate().rev() {
        let w=replay(&model,&states[i],r.start,r.end,r.h,&mut||false).unwrap();
        let v=reverse(&model,&states[i],r.start,r.end,r.h,&x,w,&mut||false).unwrap();x=v.initial;p[0]+=v.parameters[0];
    }
    assert_eq!(got.initial,x);assert_eq!(got.parameters,p);
    assert!((got.initial[0]-0.8*(-1.4f64).exp()).abs()<1e-9);
    assert!((got.parameters[0]-(0.3+0.8*2.0*1.2*(-1.4f64).exp())).abs()<1e-9);
    assert!(got.peak_checkpoints<=fs_ad::revolve::min_budget(a.accepted_steps()));
    assert!(got.replayed_steps>=a.accepted_steps());
}

#[test]
fn recorded_discrete_adjoint_matches_independent_forward_duals() {
    use fs_ad::{Dual64, Real, gradient};
    type D=Dual64<5>;
    let model=super::super::tests::Nonlinear([0.3,-0.4,0.2]);
    let mut a=RecordedRk45::new(&model,AdaptiveState::new(0.0,&[0.8,-0.2],0.1),config(0.7)).unwrap();
    a.advance(10000,10000,&mut||false).unwrap();
    let got=a.pullback(&[0.7,-1.3],&[0.0;3],budget(),&mut||false).unwrap();
    let (_,grad)=gradient([0.8,-0.2,0.3,-0.4,0.2],|p:[D;5]| {
        let mut u=[p[0],p[1]];
        for r in &a.records {
            let mut k=[[D::zero();2];7];
            for stage in 0..7 {
                let mut y=u;
                for j in 0..stage {
                    let coefficient=super::super::super::A[stage-1][j];
                    if coefficient!=0.0 {for q in 0..2 {y[q]=D::from_f64(r.h*coefficient).mul_add(k[j][q],y[q]);}}
                }
                let t=super::super::super::stage_time(r.start,r.end,r.h,stage);
                k[stage]=[p[2]*y[0]*y[1]+p[3]*D::from_f64(t),-y[0]+p[4]*y[1]*y[1]];
            }
            for q in 0..2 {
                let mut du=D::zero();
                for j in 0..7 {du=D::from_f64(super::super::super::B5[j]).mul_add(k[j][q],du);}
                u[q]=D::from_f64(r.h).mul_add(du,u[q]);
            }
        }
        D::from_f64(0.7)*u[0]-D::from_f64(1.3)*u[1]
    });
    for (v,d) in got.initial.iter().chain(&got.parameters).zip(grad) {assert!((v-d).abs()<2e-12,"{v} != {d}");}
}

#[test]
fn changing_model_cannot_silently_change_replayed_primal() {
    let model=Decay::new();let mut a=recording(&model);a.advance(10000,10000,&mut||false).unwrap();
    model.rate.set(-0.8);
    assert!(matches!(a.pullback(&[1.0],&[0.0],budget(),&mut||false),Err(TrajectoryError::ReplayMismatch{..})));
    assert_eq!(model.derivatives.get(),0);
}

#[test]
fn reverse_limits_fail_without_returning_partial_gradients() {
    let model=Decay::new();let mut a=recording(&model);a.advance(10000,10000,&mut||false).unwrap();
    let before=a.state().clone();let calls=model.calls.get();
    let small=ReplayBudget{checkpoints:0,replayed_steps:10000};
    assert!(matches!(a.pullback(&[1.0],&[0.0],small,&mut||false),Err(TrajectoryError::CheckpointLimit{..})));
    assert_eq!(model.calls.get(),calls);
    let expected=a.pullback(&[1.0],&[0.0],budget(),&mut||false).unwrap();
    let short=ReplayBudget{checkpoints:64,replayed_steps:expected.replayed_steps-1};
    assert_eq!(a.pullback(&[1.0],&[0.0],short,&mut||false),Err(TrajectoryError::ReplayLimit));
    let exact=ReplayBudget{replayed_steps:expected.replayed_steps,..short};
    assert_eq!(a.pullback(&[1.0],&[0.0],exact,&mut||false).unwrap(),expected);same(a.state(),&before);
}

#[test]
fn cancellation_in_forward_and_reverse_preserves_recording() {
    let model=Decay::new();let mut a=recording(&model);let before=a.state().clone();
    let mut polls=0;
    assert_eq!(a.advance(10000,10000,&mut||{polls+=1;polls==5}).unwrap().status,RecordingStatus::Cancelled);
    same(a.state(),&before);assert_eq!(a.accepted_steps(),0);
    a.advance(10000,10000,&mut||false).unwrap();let completed=a.state().clone();
    let calls=Cell::new(0);
    let expected=a.pullback(&[1.0],&[0.0],budget(),&mut||{calls.set(calls.get()+1);false}).unwrap();
    for stop in [1,7,19,calls.get()/2,calls.get()] {
        let mut k=0;
        assert_eq!(a.pullback(&[1.0],&[0.0],budget(),&mut||{k+=1;k==stop}),Err(TrajectoryError::Step(AdjointError::Cancelled)));
        same(a.state(),&completed);
    }
    assert_eq!(a.pullback(&[1.0],&[0.0],budget(),&mut||false).unwrap(),expected);
}

#[test]
fn incomplete_and_empty_trajectories_have_distinct_semantics() {
    let model=Decay::new();let a=recording(&model);
    assert_eq!(a.pullback(&[1.0],&[0.0],budget(),&mut||false),Err(TrajectoryError::Incomplete));
    let a=RecordedRk45::new(&model,AdaptiveState::new(2.0,&[1.2],0.1),config(2.0)).unwrap();
    let g=a.pullback(&[2.0],&[3.0],ReplayBudget{checkpoints:0,replayed_steps:0},&mut||false).unwrap();
    assert_eq!(g.initial,vec![2.0]);assert_eq!(g.parameters,vec![3.0]);assert_eq!(g.replayed_steps,0);assert_eq!(model.calls.get(),0);
}

#[test]
fn replay_retains_h_instead_of_subtracting_rounded_endpoints() {
    let model=Decay::new();let t=1e16;
    let mut cfg=config(t+4.0);cfg.rtol=10.0;cfg.atol=10.0;
    let mut a=RecordedRk45::new(&model,AdaptiveState::new(t,&[1.2],3.0),cfg).unwrap();
    assert_eq!(a.advance(1,1,&mut||false).unwrap().status,RecordingStatus::ReachedEnd);
    assert_eq!(a.records[0].h,3.0);assert_ne!(a.records[0].h,a.records[0].end-a.records[0].start);
    let g=a.pullback(&[1.0],&[0.0],budget(),&mut||false).unwrap();
    let expected=super::super::step_vjp(&model,t,&[1.2],3.0,&[1.0],4096,&mut||false).unwrap();
    assert_eq!(g.initial,expected.initial);assert_eq!(g.parameters,expected.parameters);
}
