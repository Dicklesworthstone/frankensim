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

#[test]
fn g1_native_tangent_keeps_prestress_curvature_and_the_memory_row() {
    let sys = fixed().with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain: Rc::new(QuadraticStrain), stiffness: 4.0, relaxation_time_s: 0.2,
    }]).unwrap();
    let state = [0.3, 0.0, 0.17];
    let direction = [0.4, -0.7, -0.2];
    let mut actual = [f64::NAN; 3];
    assert!(sys.storage().continuous_tangent().unwrap().apply(&state, &direction, &mut actual));
    // H_arm = (2*(q+q^2/2)-z)^2/2. Differentiating its actual force
    // includes stored stress times strain curvature, besides the outer product.
    let strain_gradient = 2.0 * (1.0 + state[0]);
    let effort = 2.0 * (state[0] + 0.5 * state[0] * state[0]) - state[2];
    let change = strain_gradient * direction[0] - direction[2];
    let expected = [strain_gradient * change + 2.0 * effort * direction[0], 0.0, -change];
    for (a, b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < 2.0e-14);
    }
    assert!((actual[0] - strain_gradient * change).abs() > 0.4);
    // A pure history perturbation produces the reciprocal mechanical reaction.
    assert!(sys.storage().continuous_tangent().unwrap().apply(
        &state, &[0.0, 0.0, 1.0], &mut actual,
    ));
    assert_eq!(actual, [-strain_gradient, 0.0, 1.0]);
}

#[derive(Debug)]
struct CoupledStrain;
impl RelaxationStrainMap for CoupledStrain {
    fn state_dimension(&self) -> usize { 3 }
    fn strain(&self, x: &[f64]) -> f64 {
        x[0] + 0.5*x[0]*x[0] + 0.3*x[0]*x[1] + 0.25*x[1]*x[1]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[1.0+x[0]+0.3*x[1], 0.3*x[0]+0.5*x[1], 0.0]);
    }
    fn hessian_vector(&self, _: &[f64], v: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[v[0]+0.3*v[1], 0.3*v[0]+0.5*v[1], 0.0]);
    }
}

#[test]
fn g3_native_tangent_composes_nested_linear_and_nonlinear_material_memory() {
    let sys = fixed().with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.5], stiffness: 9.0, relaxation_time_s: 0.7,
    }]).unwrap().with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain: Rc::new(CoupledStrain), stiffness: 4.0, relaxation_time_s: 0.2,
    }]).unwrap();
    let state = [0.3, -0.2, 0.12, -0.4];
    let direction = [0.4, -0.7, -0.2, 0.3];
    let other = [-0.2, 0.1, 0.6, -0.4];
    let tangent = sys.storage().continuous_tangent().unwrap();
    let mut actual = [f64::NAN; 4];
    assert!(tangent.apply(&state, &direction, &mut actual));
    // The original effort is cubic: this independent four-point derivative
    // is exact in exact arithmetic, including both internal-history rows.
    let epsilon = 0.01;
    let gradients: Vec<_> = [-2.0, -1.0, 1.0, 2.0].iter().map(|step| {
        let x: Vec<_> = state.iter().zip(direction).map(|(x, v)| x + step*epsilon*v).collect();
        sys.effort(&x)
    }).collect();
    for k in 0..4 {
        let reference = (gradients[0][k] - 8.0*gradients[1][k]
            + 8.0*gradients[2][k] - gradients[3][k]) / (12.0*epsilon);
        assert!((actual[k] - reference).abs() < 2.0e-11 * reference.abs().max(1.0));
    }
    let mut transposed = [0.0; 4];
    assert!(tangent.apply(&state, &other, &mut transposed));
    let left: f64 = actual.iter().zip(other).map(|(a, v)| a*v).sum();
    let right: f64 = transposed.iter().zip(direction).map(|(a, v)| a*v).sum();
    assert!((left - right).abs() < 2.0e-13 * left.abs().max(1.0));
}

#[derive(Debug)]
struct IncompleteStrain { incomplete_gradient: bool }
impl RelaxationStrainMap for IncompleteStrain {
    fn state_dimension(&self) -> usize { 2 }
    fn strain(&self, x: &[f64]) -> f64 { x[0] }
    fn gradient(&self, _: &[f64], out: &mut [f64]) {
        if self.incomplete_gradient { out[0] = 1.0; }
        else { out.copy_from_slice(&[1.0, 0.0]); }
    }
    fn hessian_vector(&self, _: &[f64], _: &[f64], out: &mut [f64]) {
        if self.incomplete_gradient { out.fill(0.0); }
        else { out[0] = 0.0; }
    }
}

#[test]
fn g0_native_relaxation_tangent_refuses_missing_derivatives_and_invalid_inputs() {
    for incomplete_gradient in [false, true] {
        let sys = fixed().with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
            strain: Rc::new(IncompleteStrain { incomplete_gradient }),
            stiffness: 4.0, relaxation_time_s: 0.2,
        }]).unwrap();
        assert!(!sys.storage().continuous_tangent().unwrap().apply(
            &[0.3, 0.1, -0.2], &[0.4, 0.5, 0.6], &mut [0.0; 3],
        ));
    }
    let sys = fixed().with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0], stiffness: 4.0, relaxation_time_s: 0.2,
    }]).unwrap();
    let tangent = sys.storage().continuous_tangent().unwrap();
    let mut out = [17.0; 3];
    for (state, direction) in [
        (vec![0.0; 2], vec![0.0; 3]),
        (vec![0.0; 3], vec![0.0; 2]),
        (vec![f64::NAN, 0.0, 0.0], vec![0.0; 3]),
        (vec![0.0; 3], vec![0.0, f64::INFINITY, 0.0]),
    ] {
        assert!(!tangent.apply(&state, &direction, &mut out));
        assert_eq!(out, [17.0; 3]);
    }
    assert!(!tangent.apply(&[0.0; 3], &[0.0; 3], &mut [0.0; 2]));
}
