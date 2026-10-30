use super::*;
use crate::{QuadraticStorage,RelaxationBranch};

#[derive(Debug)]
struct QuadraticStrain;
impl RelaxationStrainMap for QuadraticStrain {
    fn state_dimension(&self)->usize {2}
    fn strain(&self,x:&[f64])->f64 {x[0]+0.5*x[0]*x[0]}
    fn gradient(&self,x:&[f64],g:&mut[f64]) {g.copy_from_slice(&[1.0+x[0],0.0]);}
    fn hessian_vector(&self,_:&[f64],d:&[f64],g:&mut[f64]) {g.copy_from_slice(&[d[0],0.0]);}
    fn is_quadratic_strain(&self)->bool {true}
}
fn fixed()->PortHamiltonian {
    PortHamiltonian::new(2,0,vec![0.0;4],vec![0.0;4],vec![],
        Box::new(QuadraticStorage::new(vec![0.0;4],2).unwrap())).unwrap()
}
#[test]
fn g1_nonlinear_strain_relaxes_in_original_phs_with_exact_work() {
    let sys=fixed().with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain:Rc::new(QuadraticStrain),stiffness:4.0,relaxation_time_s:0.2,
    }]).unwrap();
    let mut x=vec![0.3,0.0,0.0];let before=sys.hamiltonian(&x);let mut loss=0.0;
    let target=2.0*(0.3+0.5*0.3*0.3);let dt=0.01;
    for _ in 0..20 {
        let expected=(x[2]*(1.0-dt/0.4)+dt/0.2*target)/(1.0+dt/0.4);
        let step=crate::step(&sys,&x,&[],dt).unwrap();
        assert!((step.x[2]-expected).abs()<1e-12);
        assert_eq!(step.x[0],0.3);assert_eq!(step.x[1],0.0);
        loss+=step.dissipated;x=step.x;
        assert!((sys.hamiltonian(&x)+loss-before).abs()<1e-12);
    }
}
#[test]
fn g1_exact_geometric_remainder_retains_prestressed_cross_term() {
    let coupling=RelaxationCoupling::Mapped {strain:Rc::new(QuadraticStrain),scale:3.0};
    let (a,b)=(vec![0.12,0.0],vec![0.32,0.0]);let m=vec![0.22,0.0];
    let (za,zb,zm)=(0.07,0.11,0.09);
    let actual=coupling.mapped_remainder(&a,&b,&m,za,zb,zm).unwrap();
    let ea=3.0*(a[0]+0.5*a[0]*a[0])-za;
    let eb=3.0*(b[0]+0.5*b[0]*b[0])-zb;
    let em=3.0*(m[0]+0.5*m[0]*m[0])-zm;
    let expected=0.5*(eb*eb-ea*ea)-em*(3.0*(1.0+m[0])*(b[0]-a[0])-(zb-za));
    assert!((actual-expected).abs()<1e-15);
    // The stable remainder still resolves a tiny geometric term while two
    // large energy values would erase it in endpoint subtraction.
    let huge=1e14;
    let stable=coupling.mapped_remainder(&a,&b,&m,huge,huge,huge).unwrap();
    assert!(stable.is_finite() && stable.abs()>1e-4);
}
#[test]
fn g0_linear_storage_keeps_original_quadratic_image() {
    let sys=fixed().with_relaxation_branches(vec![RelaxationBranch {
        projection:vec![1.0,0.0],stiffness:4.0,relaxation_time_s:0.2,
    }]).unwrap();
    let x=[0.3,0.7,0.1];let e=2.0*x[0]-x[2];
    assert_eq!(sys.hamiltonian(&x).to_bits(),(0.5*e*e).to_bits());
    assert_eq!(sys.effort(&x),vec![2.0*e,0.0,-e]);
    let q=sys.storage.quadratic_part().unwrap();assert!(q.exact);
    assert_eq!(q.matrix,vec![4.0,0.0,-2.0,0.0,0.0,0.0,-2.0,0.0,1.0]);
    assert!(fixed().with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain:Rc::new(QuadraticStrain),stiffness:0.0,relaxation_time_s:0.2,
    }]).is_err());
}
