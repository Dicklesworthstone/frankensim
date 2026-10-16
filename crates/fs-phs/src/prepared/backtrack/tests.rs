use super::*;
use crate::{Storage, QuadraticStorage};

// Synthetic mass hitting a rigid Hertz obstacle. This tests numerical work,
// not the accuracy of a large single timestep or a calibrated drum contact.
struct Hertz;
impl Storage for Hertz {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let penetration = (x[0]/0.02_f64.sqrt()).max(0.0);
        0.5*x[1]*x[1] + 1e8/2.5*penetration.powf(2.5)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 1e8/0.02_f64.sqrt()*(x[0]/0.02_f64.sqrt()).max(0.0).powf(1.5);
        out[1] = x[1];
    }
}
fn tangent(x: &[f64], d: &[f64], out: &mut [f64]) -> bool {
    out[0] = 1.5e8/0.02*(x[0]/0.02_f64.sqrt()).max(0.0).sqrt()*d[0];
    out[1] = d[1]; true
}
fn contact() -> PortHamiltonian {
    PortHamiltonian::new(2, 1, vec![0.,1.,-1.,0.], vec![0.;4], vec![0.,1.], Box::new(Hertz)).unwrap()
}
fn initial() -> [f64;2] { [-1e-4*0.02_f64.sqrt(),0.02_f64.sqrt()] }
fn quadratic() -> PortHamiltonian {
    PortHamiltonian::new(2, 1, vec![0.,1.,-1.,0.], vec![0.,0.,0.,0.7], vec![0.,1.],
        Box::new(QuadraticStorage::new(vec![9.,0.,0.,1.],2).unwrap())).unwrap()
}

#[test]
fn hertz_onset_converges_within_the_original_ten_update_budget() {
    let sys = contact(); let a = initial();
    let mut raw = StepWorkspace::new(&sys).unwrap(); raw.set_iteration_limit(10).unwrap();
    let (mut x, mut y) = ([123.;2],[456.]);
    assert!(matches!(raw.step_into_analytic(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent),
        Err(PhsError::NewtonStalled {..})));
    assert_eq!(x,[123.;2]); assert_eq!(y,[456.]);
    let mut guarded = StepWorkspace::new(&sys).unwrap();
    guarded.set_iteration_limit(10).unwrap(); guarded.set_backtracking_limit(8).unwrap();
    let record = guarded.step_into_analytic(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent).unwrap();
    assert!(guarded.backtracking_work().1 > 0);
    assert!(guarded.backtracking_work().0 <= 9*record.newton_iters);
    assert!(record.balance_residual().abs() < 1e-9);
    let mut reference = StepWorkspace::new(&sys).unwrap();
    let (mut expected, mut expected_y) = ([0.;2],[0.]);
    reference.step_into_analytic(&sys,&a,&[0.],0.001,&mut expected,&mut expected_y,&tangent).unwrap();
    for (x,y) in x.iter().zip(expected) { assert!((x-y).abs() < 1e-9); }
}

#[test]
fn full_updates_keep_reference_bits_and_driven_energy_accounting() {
    for analytic in [false,true] {
        let sys = quadratic();
        let mut raw = StepWorkspace::new(&sys).unwrap();
        let mut guarded = StepWorkspace::new(&sys).unwrap(); guarded.set_backtracking_limit(8).unwrap();
        let mut a = [0.2,-0.4]; let (mut x,mut rx,mut y,mut ry) = ([0.;2],[0.;2],[0.],[0.]);
        let h = |_: &[f64], d: &[f64], out: &mut [f64]| { out[0]=9.*d[0];out[1]=d[1];true };
        for i in 0..60 {
            let u = [if i<30 {0.8} else {0.}];
            let (r,g) = if analytic {
                (raw.step_into_analytic(&sys,&a,&u,0.003,&mut rx,&mut ry,&h).unwrap(),
                 guarded.step_into_analytic(&sys,&a,&u,0.003,&mut x,&mut y,&h).unwrap())
            } else {
                (raw.step_into(&sys,&a,&u,0.003,&mut rx,&mut ry).unwrap(),
                 guarded.step_into(&sys,&a,&u,0.003,&mut x,&mut y).unwrap())
            };
            assert_eq!(x,rx); assert_eq!(y,ry); assert_eq!(r,g);
            assert_eq!(guarded.backtracking_work().1,0);
            assert!(g.dissipated>=0. && g.balance_residual().abs()<1e-10); a=x;
        }
    }
}

#[test]
fn cancellation_at_every_solver_boundary_keeps_outputs_and_exact_retry() {
    let sys = contact(); let a = initial(); let (mut x,mut y) = ([0.;2],[0.]);
    let mut clean = StepWorkspace::new(&sys).unwrap(); clean.set_backtracking_limit(8).unwrap();
    let mut count=0;
    clean.step_into_analytic_controlled(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent,
        || {count+=1;false}).unwrap();
    let expected=x; let expected_y=y;
    for stop in 1..=count {
        let mut s = StepWorkspace::new(&sys).unwrap(); s.set_backtracking_limit(8).unwrap();
        x=[123.;2];y=[456.];let mut seen=0;
        let error=s.step_into_analytic_controlled(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent,
            || {seen+=1;seen==stop}).unwrap_err();
        assert_eq!(error,PreparedStepError::Cancelled); assert_eq!(x,[123.;2]); assert_eq!(y,[456.]);
        s.step_into_analytic(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent).unwrap();
        assert_eq!(x,expected); assert_eq!(y,expected_y);
    }
}

#[test]
fn exhausted_search_uses_the_original_refusal_and_configuration_is_transactional() {
    let sys=contact(); let a=initial(); let mut s=StepWorkspace::new(&sys).unwrap();
    s.set_backtracking_limit(1).unwrap(); assert!(s.set_backtracking_limit(17).is_err());
    assert_eq!(s.backtracking_halvings,1);
    let (mut x,mut y)=([123.;2],[456.]);
    assert!(matches!(s.step_into_analytic(&sys,&a,&[0.],0.01,&mut x,&mut y,&tangent),
        Err(PhsError::NewtonStalled {..})));
    assert_eq!(x,[123.;2]); assert_eq!(y,[456.]);
    s.set_backtracking_limit(0).unwrap();
    s.step_into_analytic(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent).unwrap();
    assert_eq!(s.backtracking_work(),(0,0));
    s.set_backtracking_limit(8).unwrap(); s.set_iteration_limit(0).unwrap();
    x=[123.;2];y=[456.];
    assert!(s.step_into_analytic(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent).is_err());
    assert_eq!(x,[123.;2]); assert_eq!(y,[456.]); assert_eq!(s.backtracking_work(),(0,0));
}

#[test]
fn condensed_and_full_tangents_solve_the_same_globalized_equation() {
    let sys=contact();let a=initial();
    let mut full=StepWorkspace::new(&sys).unwrap();full.set_backtracking_limit(8).unwrap();
    let mut condensed=StepWorkspace::new(&sys).unwrap();condensed.set_backtracking_limit(8).unwrap();
    condensed.set_condensed_pairs(&[[0,1]]).unwrap();
    let (mut x,mut c,mut y,mut cy)=([0.;2],[0.;2],[0.],[0.]);
    let r=full.step_into_analytic(&sys,&a,&[0.],0.001,&mut x,&mut y,&tangent).unwrap();
    let cr=condensed.step_into_analytic(&sys,&a,&[0.],0.001,&mut c,&mut cy,&tangent).unwrap();
    assert!(r.balance_residual().abs()<1e-9 && cr.balance_residual().abs()<1e-9);
    for (x,c) in x.iter().zip(c) {assert!((x-c).abs()<1e-8);}
    assert!(condensed.linear_solve_counts().0>0);
}
