//! G0/G1/G4/G5: recursive correction, literal Galerkin products and real CG.
use std::cell::Cell;
use std::ops::ControlFlow;
use fs_solver::{CgState, LinearOp};
use fs_solver::op::multilevel::{
    MultilevelBudget, MultilevelControl, MultilevelError, SparseMultilevel,
    SymmetricGalerkin, sparse_galerkin,
};
use fs_sparse::{Coo, Csr, precond::{IdentityPrecond, Precond}};
struct Counted { a: Csr, applications: Cell<usize> }
impl LinearOp for Counted {
    fn n(&self) -> usize { self.a.nrows() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.applications.set(self.applications.get()+1); self.a.spmv(x,y);
    }
}
fn poisson(n: usize, scale: f64) -> Csr {
    let mut a = Coo::new(n,n);
    for i in 0..n { a.push(i,i,2.0*scale); if i>0 {a.push(i,i-1,-scale);} if i+1<n {a.push(i,i+1,-scale);} }
    a.assemble()
}
fn interpolation(n: usize) -> Csr {
    assert!(n>1 && n%2==1); let nc=(n-1)/2; let mut p=Coo::new(n,nc);
    for i in 0..n {
        if i%2==1 {p.push(i,i/2,1.0);} else {
            if i>0 {p.push(i,i/2-1,0.5);} if i/2<nc {p.push(i,i/2,0.5);}
        }
    }
    p.assemble()
}
fn hierarchy(mut n: usize) -> Vec<Csr> {
    let mut p=Vec::new(); while n>15 {p.push(interpolation(n)); n=(n-1)/2;} p
}
fn dot(x:&[f64],y:&[f64])->f64 {x.iter().zip(y).map(|(x,y)|x*y).sum()}
fn act(p:&impl Precond,r:&[f64])->Vec<f64> {let mut z=vec![0.0;r.len()];p.apply(r,&mut z);z}

#[test]
fn g1_recursive_coarse_spaces_exceed_the_old_512_coordinate_limit() {
    let op=Counted {a:poisson(2047,1.0),applications:Cell::new(0)};
    let mut poll=|_|ControlFlow::Continue(()); let mut c=MultilevelControl::new(MultilevelBudget::default(),&mut poll).unwrap();
    // Analytically P^T tridiag(-1,2,-1) P = 0.5*tridiag(-1,2,-1).
    let mg=SparseMultilevel::new(&op,&vec![0.5;op.n()],hierarchy(op.n()),poisson(1023,0.5),&mut c).unwrap();
    assert_eq!(mg.level_sizes(),[2047,1023,511,255,127,63,31,15]);
    assert_eq!(op.applications.get(),0,"setup must not probe the fine matrix");
    assert!(mg.work().matrix_entries<12*op.n());
    let rhs:Vec<f64>=(0..op.n()).map(|i|1.0+(i%13) as f64/13.0).collect();
    let mut cg=CgState::new(&op,&mg,&rhs);
    assert!(cg.run(&op,&mg,1e-10,150).converged);
    let mut residual=vec![0.0;op.n()];op.apply(&cg.x,&mut residual);
    for (r,b) in residual.iter_mut().zip(&rhs){*r-=b;}
    assert!((dot(&residual,&residual)/dot(&rhs,&rhs)).sqrt()<1e-8);
}

#[test]
fn g0_sparse_galerkin_matches_independent_column_probing() {
    let a=poisson(31,3.0);let p=interpolation(31);let mut callback=|_|ControlFlow::Continue(());
    let mut c=MultilevelControl::new(MultilevelBudget::default(),&mut callback).unwrap();
    let result=sparse_galerkin(&a,&p,&mut c).unwrap();
    for j in 0..p.ncols() {
        let mut basis=vec![0.0;p.ncols()];basis[j]=1.0;
        let mut x=vec![0.0;p.nrows()];p.spmv(&basis,&mut x);
        let mut ax=vec![0.0;p.nrows()];a.spmv(&x,&mut ax);
        for i in 0..p.ncols() {
            basis.fill(0.0);basis[i]=1.0;let mut y=vec![0.0;p.nrows()];p.spmv(&basis,&mut y);
            assert_eq!(result.get(i,j),dot(&y,&ax));
        }
    }
    assert!(c.work().galerkin_products>0);
}

#[test]
fn g0_symmetric_fixed_linear_positive_action_and_repeatable_cg() {
    let op=Counted{a:poisson(127,1.0),applications:Cell::new(0)};
    let mut callback=|_|ControlFlow::Continue(());let mut c=MultilevelControl::new(MultilevelBudget::default(),&mut callback).unwrap();
    let mg=SparseMultilevel::new(&op,&vec![0.5;op.n()],hierarchy(op.n()),poisson(63,0.5),&mut c).unwrap();
    let x:Vec<f64>=(0..op.n()).map(|i|((i*17)%23) as f64/23.0-0.5).collect();
    let y:Vec<f64>=(0..op.n()).map(|i|((i*11)%19) as f64/19.0-0.5).collect();
    let bx=act(&mg,&x);let by=act(&mg,&y);let lhs=dot(&x,&by);let rhs=dot(&y,&bx);
    assert!((lhs-rhs).abs()<1e-11*lhs.abs().max(rhs.abs()).max(1.0));assert!(dot(&x,&bx)>0.0);
    let sum:Vec<f64>=x.iter().zip(&y).map(|(x,y)|2.0*x-0.5*y).collect();let bs=act(&mg,&sum);
    for ((s,x),y) in bs.iter().zip(&bx).zip(&by){assert!((s-(2.0*x-0.5*y)).abs()<1e-10);}
    assert_eq!(bx,act(&mg,&x));
    let mut a=CgState::new(&op,&mg,&x);let mut b=CgState::new(&op,&mg,&x);
    assert!(a.run(&op,&mg,1e-11,200).converged);assert!(b.run(&op,&mg,1e-11,200).converged);
    assert_eq!(a.x,b.x);assert_eq!(a.iters,b.iters);
    let mut plain=CgState::new(&op,&IdentityPrecond,&x);assert!(plain.run(&op,&IdentityPrecond,1e-11,300).converged);
    assert!(a.iters<plain.iters);
}

#[test]
fn g4_product_entry_bottom_budgets_and_cancellation_refuse_whole_hierarchy() {
    let op=Counted{a:poisson(127,1.0),applications:Cell::new(0)};
    let mut callback=|_|ControlFlow::Continue(());
    let mut c=MultilevelControl::new(MultilevelBudget{max_galerkin_products:3,..Default::default()},&mut callback).unwrap();
    assert!(matches!(SparseMultilevel::new(&op,&vec![0.5;127],hierarchy(127),poisson(63,0.5),&mut c),Err(MultilevelError::Budget(_))));
    assert_eq!(c.work().galerkin_products,3);assert_eq!(op.applications.get(),0);
    let mut callback=|w:fs_solver::op::multilevel::MultilevelWork| if w.galerkin_products>0 {ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=MultilevelControl::new(MultilevelBudget::default(),&mut callback).unwrap();
    assert!(matches!(SparseMultilevel::new(&op,&vec![0.5;127],hierarchy(127),poisson(63,0.5),&mut c),Err(MultilevelError::Cancelled)));
    assert!(c.work().galerkin_products>0);
    for budget in [MultilevelBudget{max_matrix_entries:1,..Default::default()},MultilevelBudget{max_coarsest_dofs:1,..Default::default()}] {
        let mut callback=|_|ControlFlow::Continue(());let mut c=MultilevelControl::new(budget,&mut callback).unwrap();
        assert!(matches!(SparseMultilevel::new(&op,&vec![0.5;127],hierarchy(127),poisson(63,0.5),&mut c),Err(MultilevelError::Budget(_))));
    }
}

#[test]
fn g0_no_nonsymmetric_law_nonfinite_contribution_or_fake_bottom_pivot() {
    let mut poll=|_|ControlFlow::Continue(());let mut c=MultilevelControl::new(MultilevelBudget::default(),&mut poll).unwrap();
    let mut acc=SymmetricGalerkin::new(3,&mut c).unwrap();
    assert!(matches!(acc.add(0,0,f64::NAN,&mut c),Err(MultilevelError::Invalid(_))));
    let mut a=Coo::new(3,3);for i in 0..3 {a.push(i,i,1.0);}a.push(0,1,0.2);
    let op=Counted{a:poisson(7,1.0),applications:Cell::new(0)};
    assert!(matches!(SparseMultilevel::new(&op,&[0.5;7],vec![interpolation(7)],a.assemble(),&mut c),Err(MultilevelError::Nonsymmetric)));
    let mut a=Coo::new(3,3);for i in 0..3 {for j in 0..3 {a.push(i,j,1.0);}}
    assert!(matches!(SparseMultilevel::new(&op,&[0.5;7],vec![interpolation(7)],a.assemble(),&mut c),Err(MultilevelError::NotPositiveDefinite)));
}
