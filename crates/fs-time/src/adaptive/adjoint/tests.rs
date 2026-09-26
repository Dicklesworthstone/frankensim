use super::*;
use crate::{AdaptiveState, AdaptiveStatus, PiController, rk45_adaptive_checked};
use std::cell::Cell;

struct Decay(f64);
impl OdeVjp for Decay {
    fn dimension(&self) -> usize { 1 }
    fn parameter_count(&self) -> usize { 1 }
    fn rhs(&self, _: f64, u: &[f64], out: &mut [f64]) { out[0] = self.0 * u[0]; }
    fn rhs_vjp(&self, _: f64, u: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        x[0] = self.0 * seed[0]; p[0] = u[0] * seed[0]; Ok(())
    }
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol, "{a} != {b}"); }

#[test]
fn step_matches_production_and_independent_stability_polynomial() {
    let model = Decay(-0.7);
    let result = step_vjp(&model, 0.25, &[1.2], 0.4, &[2.0], 19, &mut || false).unwrap();
    let z = -0.28f64;
    let r = 1.0 + z + z.powi(2)/2.0 + z.powi(3)/6.0 + z.powi(4)/24.0 + z.powi(5)/120.0 + z.powi(6)/600.0;
    let dr = 1.0 + z + z.powi(2)/2.0 + z.powi(3)/6.0 + z.powi(4)/24.0 + z.powi(5)/100.0;
    close(result.value[0], 1.2*r, 1e-14);
    close(result.initial[0], 2.0*r, 1e-14);
    close(result.parameters[0], 2.0*1.2*0.4*dr, 1e-14);
    let mut primal = AdaptiveState::new(0.25, &[1.2], 0.4);
    let report = rk45_adaptive_checked(&mut primal, &|t,u,o| model.rhs(t,u,o), 1.0,
        1.0, 1.0, &PiController::default(), 1).unwrap();
    assert_eq!(report.status, AdaptiveStatus::StepLimit);
    assert_eq!(primal.accepted, 1);
    assert_eq!(primal.u[0].to_bits(), result.value[0].to_bits());
}

#[derive(Clone)]
pub(super) struct Nonlinear(pub(super) [f64; 3]);
impl OdeVjp for Nonlinear {
    fn dimension(&self) -> usize { 2 }
    fn parameter_count(&self) -> usize { 3 }
    fn rhs(&self, t: f64, x: &[f64], f: &mut [f64]) {
        f[0] = self.0[0]*x[0]*x[1] + self.0[1]*t;
        f[1] = -x[0] + self.0[2]*x[1]*x[1];
    }
    fn rhs_vjp(&self, t: f64, x: &[f64], b: &[f64], y: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        y[0] = self.0[0]*x[1]*b[0] - b[1];
        y[1] = self.0[0]*x[0]*b[0] + 2.0*self.0[2]*x[1]*b[1];
        p[0] = x[0]*x[1]*b[0]; p[1] = t*b[0]; p[2] = x[1]*x[1]*b[1]; Ok(())
    }
}

#[test]
fn nonlinear_nonsymmetric_transpose_matches_finite_differences() {
    let seed = [0.7, -1.3]; let x = [0.8, -0.2];
    let model = Nonlinear([0.3, -0.4, 0.2]);
    let result = step_vjp(&model, 0.5, &x, 0.2, &seed, 40, &mut || false).unwrap();
    let value = |model: &Nonlinear, u: &[f64]| {
        let w = replay(model, u, 0.5, 0.7, 0.2, &mut || false).unwrap();
        w.next.iter().zip(seed).map(|(v,b)| v*b).sum::<f64>()
    };
    let eps = 1e-6;
    for i in 0..2 {
        let (mut a, mut b) = (x,x); a[i]+=eps; b[i]-=eps;
        close(result.initial[i], (value(&model,&a)-value(&model,&b))/(2.0*eps), 2e-9);
    }
    for i in 0..3 {
        let (mut a, mut b) = (model.clone(),model.clone()); a.0[i]+=eps; b.0[i]-=eps;
        close(result.parameters[i], (value(&a,&x)-value(&b,&x))/(2.0*eps), 2e-9);
    }
}

#[test]
fn invalid_inputs_and_workspace_caps_refuse_before_rhs() {
    struct NoCall;
    impl OdeVjp for NoCall {
        fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {1}
        fn rhs(&self,_:f64,_:&[f64],_:&mut[f64]) {panic!("unexpected RHS")}
        fn rhs_vjp(&self,_:f64,_:&[f64],_:&[f64],_:&mut[f64],_:&mut[f64])->Result<(),String> {panic!("unexpected VJP")}
    }
    assert_eq!(step_vjp(&NoCall,0.0,&[1.0],0.1,&[1.0],18,&mut||false),
        Err(AdjointError::WorkspaceLimit{required:19,limit:18}));
    for h in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(step_vjp(&NoCall,0.0,&[1.0],h,&[1.0],19,&mut||false).is_err());
    }
    assert!(step_vjp(&NoCall,0.0,&[1.0],0.1,&[],19,&mut||false).is_err());
}

#[test]
fn cancellation_at_every_boundary_returns_no_partial_pullback() {
    let count=Cell::new(0);
    let expected=step_vjp(&Decay(-0.7),0.0,&[1.2],0.2,&[1.0],19,
        &mut||{count.set(count.get()+1);false}).unwrap();
    for stop in 1..=count.get() {
        let mut polls=0;
        assert_eq!(step_vjp(&Decay(-0.7),0.0,&[1.2],0.2,&[1.0],19,
            &mut||{polls+=1;polls==stop}),Err(AdjointError::Cancelled));
    }
    assert_eq!(step_vjp(&Decay(-0.7),0.0,&[1.2],0.2,&[1.0],19,&mut||false).unwrap(),expected);
}

#[test]
fn unwritten_derivative_components_are_rejected() {
    struct Bad;
    impl OdeVjp for Bad {
        fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {1}
        fn rhs(&self,_:f64,_:&[f64],out:&mut[f64]) {out[0]=0.0;}
        fn rhs_vjp(&self,_:f64,_:&[f64],_:&[f64],x:&mut[f64],_:&mut[f64])->Result<(),String> {x[0]=0.0;Ok(())}
    }
    assert_eq!(step_vjp(&Bad,0.0,&[1.0],0.1,&[1.0],19,&mut||false),Err(AdjointError::NonFiniteDerivative));
}

#[test]
fn parameter_free_model_can_pull_back_initial_conditions() {
    struct Flow;
    impl OdeVjp for Flow {
        fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {0}
        fn rhs(&self,t:f64,_:&[f64],out:&mut[f64]) {out[0]=t;}
        fn rhs_vjp(&self,_:f64,_:&[f64],_:&[f64],x:&mut[f64],p:&mut[f64])->Result<(),String> {
            x[0]=0.0;assert!(p.is_empty());Ok(())
        }
    }
    let r=step_vjp(&Flow,0.0,&[1.0],0.5,&[3.0],17,&mut||false).unwrap();
    close(r.value[0],1.125,1e-14);assert_eq!(r.initial,vec![3.0]);assert!(r.parameters.is_empty());
}
