//! Independent dense and perturbation checks of the actual tet operator.
use std::ops::ControlFlow;
use fs_topopt::{DensityElasticity, SolveBudget, SolveControl, SolveProgress, EvaluationStop};
use fs_topopt::elasticity::second_order::{ComplianceError, ComplianceGoal, ComplianceSolveOptions};
use fs_topopt::pipeline::LoadCase;

fn fixture(n: usize) -> (DensityElasticity, Vec<f64>, Vec<f64>) {
    let (mesh, positions) = fs_feec::kuhn_cube(n);
    let el = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let mut a = vec![0.0; el.n()];
    let mut b = a.clone();
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 { a[3*i+2] = -0.2; b[3*i+1] = 0.13 * (1.0 + p[2]); }
    }
    (el, a, b)
}
fn close(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    let error = a.iter().zip(b).map(|(a,b)| (a-b).abs()).fold(0.0, f64::max);
    let scale = b.iter().map(|b| b.abs()).fold(1.0, f64::max);
    assert!(error < tol*scale, "error {error:e}, scale {scale:e}, tolerance {tol:e}");
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn matvec(a: &[f64], x: &[f64]) -> Vec<f64> {
    a.chunks(x.len()).map(|row| dot(row,x)).collect()
}
// Test-only Gaussian elimination with pivoting, independent of production CG.
fn dense_solve(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n).max_by(|&i,&j| a[i*n+k].abs().total_cmp(&a[j*n+k].abs())).unwrap();
        assert!(a[pivot*n+k].abs() > 1e-15);
        for j in 0..n { a.swap(k*n+j, pivot*n+j); }
        b.swap(k,pivot);
        for i in k+1..n {
            let ratio = a[i*n+k]/a[k*n+k];
            for j in k..n { a[i*n+j] -= ratio*a[k*n+j]; }
            b[i] -= ratio*b[k];
        }
    }
    for i in (0..n).rev() {
        for j in i+1..n { b[i] -= a[i*n+j]*b[j]; }
        b[i] /= a[i*n+i];
    }
    b
}

#[test]
fn actual_tet_primal_gradient_and_hessian_match_independent_dense_oracle() {
    let (el, a, b) = fixture(1);
    let (mut reference, _, _) = fixture(1);
    let m: Vec<f64> = (0..el.cells()).map(|e| 0.7 + 0.07*e as f64).collect();
    let v: Vec<f64> = (0..el.cells()).map(|e| 0.1*(e as f64-2.0)).collect();
    reference.moduli = m.clone();
    let (matrix, _, free) = reference.assemble_dense(&vec![1.0; el.cells()]);
    let mut blocks = Vec::new();
    for e in 0..el.cells() {
        reference.moduli.fill(0.0); reference.moduli[e] = 1.0;
        blocks.push(reference.assemble_dense(&vec![1.0; el.cells()]).0);
    }
    let mut dk = vec![0.0; matrix.len()];
    for (block, &v) in blocks.iter().zip(&v) {
        for (d, &k) in dk.iter_mut().zip(block) { *d += v*k; }
    }
    let mut expected_value = 0.0;
    let mut expected_gradient = vec![0.0; el.cells()];
    let mut expected_hv = vec![0.0; el.cells()];
    for (force, weight) in [(&a, 0.7), (&b, 1.3)] {
        let rhs: Vec<f64> = free.iter().map(|&d| force[d]).collect();
        let u = dense_solve(matrix.clone(), rhs.clone());
        expected_value += weight*dot(&rhs,&u);
        let du = dense_solve(matrix.clone(), matvec(&dk,&u).iter().map(|x| -*x).collect());
        for (e, block) in blocks.iter().enumerate() {
            let ku = matvec(block,&u);
            expected_gradient[e] -= weight*dot(&u,&ku);
            expected_hv[e] -= 2.0*weight*dot(&du,&ku);
        }
    }
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let goal = ComplianceGoal::new(&el,&m,
        &[LoadCase { force: &a, weight: 0.7 },LoadCase { force: &b, weight: 1.3 }],
        ComplianceSolveOptions::default(),&mut control).unwrap();
    close(&[goal.value()], &[expected_value], 1e-9);
    close(goal.gradient(), &expected_gradient, 1e-9);
    close(&goal.hessian_vector_product(&v,&mut control).unwrap(), &expected_hv, 1e-8);
    assert_eq!(control.work().linear_solves, 4);
    assert!(goal.primal_relative_residuals().iter().all(|r| *r < 1e-10));
}

#[test]
fn heterogeneous_gradient_differences_symmetry_and_uniform_scaling() {
    let (el,a,b) = fixture(2);
    let m: Vec<f64> = (0..el.cells()).map(|e| 0.6+0.08*(e%9) as f64).collect();
    let v: Vec<f64> = (0..el.cells()).map(|e| 0.03*((e%7) as f64-3.0)).collect();
    let w: Vec<f64> = (0..el.cells()).map(|e| 0.05*((e%5) as f64-2.0)).collect();
    let loads = [LoadCase { force: &a, weight: 0.7 },LoadCase { force: &b, weight: 1.3 }];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let goal = ComplianceGoal::new(&el,&m,&loads,ComplianceSolveOptions::default(),&mut control).unwrap();
    let hv = goal.hessian_vector_product(&v,&mut control).unwrap();
    let hw = goal.hessian_vector_product(&w,&mut control).unwrap();
    close(&[dot(&w,&hv)], &[dot(&v,&hw)], 2e-8);
    assert!(dot(&v,&hv) >= -1e-9, "compliance is convex in linear moduli");
    let eps = 1e-4;
    let plus: Vec<f64> = m.iter().zip(&v).map(|(m,v)| m+eps*v).collect();
    let minus: Vec<f64> = m.iter().zip(&v).map(|(m,v)| m-eps*v).collect();
    let gp = ComplianceGoal::new(&el,&plus,&loads,ComplianceSolveOptions::default(),&mut control).unwrap();
    let gm = ComplianceGoal::new(&el,&minus,&loads,ComplianceSolveOptions::default(),&mut control).unwrap();
    let fd: Vec<f64> = gp.gradient().iter().zip(gm.gradient()).map(|(a,b)| (a-b)/(2.0*eps)).collect();
    close(&hv,&fd,2e-6);
    close(&[(gp.value()-gm.value())/(2.0*eps)], &[dot(goal.gradient(),&v)],2e-6);
    let uniform = goal.hessian_vector_product(&m,&mut control).unwrap();
    close(&uniform,&goal.gradient().iter().map(|g| -2.0*g).collect::<Vec<_>>(),1e-8);
    close(&[dot(goal.gradient(),&m)], &[-goal.value()],1e-9);
    assert_eq!(el.moduli,vec![1.0;el.cells()],"evaluations must not mutate operator scales");
}

#[test]
fn independent_opposite_loads_do_not_cancel_and_fixed_forces_do_no_work() {
    let (el,mut a,_) = fixture(1);
    for (d,&free) in el.free().iter().enumerate() { if !free { a[d] = 37.0; } }
    let b: Vec<f64> = a.iter().map(|a| -*a).collect();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut callback);
    let goal = ComplianceGoal::new(&el,&vec![1.0;el.cells()],
        &[LoadCase { force: &a,weight: 0.5 },LoadCase { force: &b,weight: 0.5 }],
        ComplianceSolveOptions::default(),&mut control).unwrap();
    assert!(goal.value() > 0.0);
    close(&[goal.value()],&[goal.case_compliances()[0]],1e-12);
    close(&[goal.value()],&[goal.case_compliances()[1]],1e-12);
    for u in goal.displacements() {
        for (d,&free) in el.free().iter().enumerate() { if !free { assert_eq!(u[d],0.0); } }
    }
    let before = control.work();
    assert_eq!(goal.hessian_vector_product(&vec![0.0;el.cells()],&mut control).unwrap(),vec![0.0;el.cells()]);
    assert_eq!(control.work(),before);
}

#[test]
fn budget_and_cancelled_tangents_leave_the_accepted_linearization_usable() {
    let (el,a,_) = fixture(2);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut callback);
    let goal = ComplianceGoal::new(&el,&vec![1.0;el.cells()],&[LoadCase { force: &a,weight: 1.0 }],
        ComplianceSolveOptions::default(),&mut control).unwrap();
    let v: Vec<f64> = (0..el.cells()).map(|e| 0.1+(e%3) as f64*0.1).collect();
    let reference = goal.hessian_vector_product(&v,&mut control).unwrap();
    let mut budget_callback = |_| ControlFlow::Continue(());
    let mut bounded = SolveControl::new(SolveBudget { per_solve_iterations: 1000,total_iterations: 1 },&mut budget_callback);
    assert!(matches!(goal.hessian_vector_product(&v,&mut bounded),
        Err(ComplianceError::Stopped(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(bounded.work().linear_iterations,1);
    let mut cancel = |p: SolveProgress| if p.work.linear_iterations > 0 {
        ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    let mut cancelled = SolveControl::new(SolveBudget::default(),&mut cancel);
    assert!(matches!(goal.hessian_vector_product(&v,&mut cancelled),
        Err(ComplianceError::Stopped(EvaluationStop::Cancelled))));
    assert!(cancelled.work().linear_iterations > 0);
    assert_eq!(goal.hessian_vector_product(&v,&mut control).unwrap(),reference);
}

#[test]
fn invalid_data_refuses_before_solves() {
    let (el,a,_) = fixture(1);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut callback);
    for bad in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        let mut m = vec![1.0;el.cells()]; m[0] = bad;
        assert!(ComplianceGoal::new(&el,&m,&[LoadCase { force: &a,weight: 1.0 }],
            ComplianceSolveOptions::default(),&mut control).is_err());
    }
    assert!(ComplianceGoal::new(&el,&vec![1.0;el.cells()],&[LoadCase { force: &a,weight: -1.0 }],
        ComplianceSolveOptions::default(),&mut control).is_err());
    assert_eq!(control.work().linear_solves,0);
}
