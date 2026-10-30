//! Manufactured source laws, resolved operating points, and original resonances.
use fs_couple::stationary_equilibrium::{solve_equilibrium,EquilibriumSolveOptions,EquilibriumSolveError};
use fs_phs::harmonic::equilibrium::{LinearizationBudget,NonlinearDissipation};
use fs_phs::harmonic::equilibrium::stationary::{StationaryCoordinates,StationaryProblem};
use fs_phs::harmonic::HarmonicBudget;
use fs_phs::{PortHamiltonian,DescriptorPortHamiltonian,Storage,QuadraticStorage,RelaxationBranch};
use fs_solver::{NewtonKrylovConfig,Globalization,LineSearchConfig};
use fs_math::c64::C64;
use std::cell::Cell;

#[derive(Clone,Copy)]
struct Spring;
impl Storage for Spring {
    fn hamiltonian(&self,x:&[f64])->f64 {0.5*x[0]*x[0]+25.0*x[0].powi(4)+0.25*x[1]*x[1]}
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        out.fill(0.0);out[0]=x[0]+100.0*x[0].powi(3);out[1]=0.5*x[1];
    }
}
fn hessian(x:&[f64],d:&[f64],out:&mut[f64])->bool {
    out.fill(0.0);out[0]=(1.0+300.0*x[0]*x[0])*d[0];out[1]=0.5*d[1];true
}
fn source()->PortHamiltonian {
    PortHamiltonian::new(2,1,vec![0.,1.,-1.,0.],vec![0.,0.,0.,0.3],
        vec![0.,1.],Box::new(Spring)).unwrap()
}
fn budget()->LinearizationBudget {
    LinearizationBudget {harmonic:HarmonicBudget {max_dimension:16,max_dense_bytes:1_000_000,
        max_dense_work:1_000_000,relative_residual_tolerance:1e-10},max_hessian_actions:16}
}
fn options()->EquilibriumSolveOptions {
    EquilibriumSolveOptions {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,
        relative_tolerance:0.0,linear_restart:8,max_linear_cycles:8,
        forcing_minimum:1e-12,forcing_maximum:1e-6,
        globalization:Globalization::LineSearch(LineSearchConfig{minimum_step:1.0/65536.0,
            ..Default::default()}),..Default::default()},relative_balance_tolerance:1e-10,
        max_iterations:40,max_residual_calls:256,max_jacobian_calls:2048,max_workspace_bytes:1_000_000}
}
fn coordinates<'a>(input:&'a[f64],states:&'a[f64],rows:&'a[f64])->StationaryCoordinates<'a> {
    StationaryCoordinates{input,state_scales:states,balance_scales:rows}
}
fn close(a:C64,b:C64,t:f64) {assert!((a-b).abs()<t*(1.0+b.abs()),"{a:?} vs {b:?}");}
fn admit_all(_: &[f64])->bool {true}
fn refuse_all(_: &[f64])->bool {false}

#[test]
fn g1_g3_source_preload_is_solved_before_its_original_resonance() {
    let system=source();let rejected=Cell::new(0);
    let domain=|x:&[f64]| {let valid=x[0].abs()<1.0;if !valid {rejected.set(rejected.get()+1);}valid};
    let problem=StationaryProblem::from_ode(&system,coordinates(&[3.],&[1.,1.],&[1.,3.]),
        budget(),&hessian,None,&domain).unwrap();
    let result=solve_equilibrium(&problem,&[0.,0.],options(),||false).unwrap();
    let image=&result.linearization;
    assert!((image.equilibrium_state()[0]-0.3).abs()<1e-11);
    assert!(image.equilibrium_state()[1].abs()<1e-11);
    assert!(rejected.get()>0,"full Newton trial must respect the source interval");
    assert!((image.tangent_storage_matrix()[0]-28.0).abs()<1e-8);
    for w in [1.,3.,4.,8.] {
        let r=image.model().factor(w,||false).unwrap().response(&[C64::ONE],||false).unwrap();
        close(r.state[0],C64::new(28.-2.*w*w,-0.3*w).recip(),1e-9);
    }
    assert!(image.maximum_scaled_residual()<=options().relative_balance_tolerance);
}

#[test]
fn g1_scaled_original_descriptor_resolves_nonzero_constraint_reaction() {
    let system=DescriptorPortHamiltonian::new(3,2,1,
        vec![0.,1.,0.,-1.,0.,2.,0.,-2.,0.],
        vec![0.,0.,0.,0.,0.3,0.,0.,0.,4.],vec![0.,1.,2.],Box::new(Spring)).unwrap();
    for states in [[1.,1.,1.],[0.01,100.,0.1]] {
        let problem=StationaryProblem::from_descriptor(&system,
            coordinates(&[1.5],&states,&[1.,3.,3.]),budget(),&hessian,None,&admit_all).unwrap();
        let solved=solve_equilibrium(&problem,&[0.2,0.1,-0.1],options(),||false).unwrap();
        let x=solved.linearization.equilibrium_state();
        assert!((x[0]-0.3).abs()<1e-9 && x[1].abs()<1e-9 && (x[2]-0.75).abs()<1e-9);
        assert_eq!(solved.linearization.tangent_storage_matrix()[8],0.0);
        let w=3.;let r=solved.linearization.model().factor(w,||false).unwrap()
            .response(&[C64::ONE],||false).unwrap();
        // G also acts on the algebraic row: the effective mechanical input is 2u.
        close(r.state[0],C64::from_re(2.)/C64::new(28.-2.*w*w,-1.3*w),1e-9);
    }
}

#[test]
fn g1_original_nonlinear_resistor_sets_bias_and_incremental_decay() {
    let system=PortHamiltonian::new(1,1,vec![0.],vec![0.2],vec![1.],
        Box::new(QuadraticStorage::new(vec![1.],1).unwrap())).unwrap();
    let action=|_:&[f64],d:&[f64],o:&mut[f64]| {o[0]=d[0];true};
    let flow=|_:&[f64],e:&[f64],o:&mut[f64]| {o[0]=e[0]+e[0].powi(3);true};
    let tangent=|_:&[f64],e:&[f64],_:&[f64],de:&[f64],o:&mut[f64]| {
        o[0]=(1.+3.*e[0]*e[0])*de[0];true};
    let loss=NonlinearDissipation{flow:&flow,tangent:&tangent,max_tangent_actions:1};
    let problem=StationaryProblem::from_ode(&system,coordinates(&[10.4],&[1.],&[10.4]),
        budget(),&action,Some(loss),&admit_all).unwrap();
    let solved=solve_equilibrium(&problem,&[0.5],options(),||false).unwrap();
    assert!((solved.linearization.equilibrium_state()[0]-2.).abs()<1e-10);
    assert!((solved.linearization.equilibrium_nonlinear_dissipated_power_w()-20.).abs()<1e-8);
    let r=solved.linearization.model().factor(3.,||false).unwrap().response(&[C64::ONE],||false).unwrap();
    close(r.state[0],C64::new(13.2,-3.).recip(),1e-10);
}

#[test]
fn g1_relaxation_history_is_solved_not_reset_or_eliminated() {
    let system=source().with_relaxation_branches(vec![RelaxationBranch{
        projection:vec![1.,0.],stiffness:9.,relaxation_time_s:0.2}]).unwrap();
    let action=|x:&[f64],d:&[f64],o:&mut[f64]| {
        hessian(x,d,o);let e=3.*d[0]-d[2];o[0]+=3.*e;o[2]=-e;true};
    let problem=StationaryProblem::from_ode(&system,coordinates(&[3.],&[1.,1.,1.],&[1.,3.,3.]),
        budget(),&action,None,&admit_all).unwrap();
    let solved=solve_equilibrium(&problem,&[0.1,0.,-0.2],options(),||false).unwrap();
    let x=solved.linearization.equilibrium_state();
    assert!((x[0]-0.3).abs()<1e-9 && (x[2]-0.9).abs()<1e-9);
    let w=4.;let s=C64::new(0.,-w);
    let modulus=(s.scale(0.2)/(C64::ONE+s.scale(0.2))).scale(9.);
    let r=solved.linearization.model().factor(w,||false).unwrap().response(&[C64::ONE],||false).unwrap();
    close(r.state[0],(C64::new(28.-2.*w*w,-0.3*w)+modulus).recip(),1e-9);
    close(r.state[2],r.state[0].scale(3.)/(C64::ONE+s.scale(0.2)),1e-9);
}

#[test]
fn g1_matrix_free_action_and_exact_point_cache_preserve_original_derivatives() {
    let system=source();let calls=Cell::new(0);
    let action=|x:&[f64],d:&[f64],o:&mut[f64]| {calls.set(calls.get()+1);hessian(x,d,o)};
    let problem=StationaryProblem::from_ode(&system,coordinates(&[3.],&[0.1,2.],&[1.,3.]),
        budget(),&action,None,&admit_all).unwrap();
    let mut scratch=problem.workspace();let z=[2.,0.3];let d=[0.7,-0.2];let mut j=[0.;2];
    problem.jacobian_action_into(&z,&d,&mut scratch,&mut j,||false).unwrap();
    assert_eq!(calls.get(),1,"one analytic action, not a dense Hessian");
    let h=1e-5;let mut plus=[0.;2];let mut minus=[0.;2];
    problem.residual_into(&[z[0]+h*d[0],z[1]+h*d[1]],&mut scratch,&mut plus,||false).unwrap();
    problem.residual_into(&[z[0]-h*d[0],z[1]-h*d[1]],&mut scratch,&mut minus,||false).unwrap();
    for i in 0..2 {assert!(((plus[i]-minus[i])/(2.*h)-j[i]).abs()<1e-8);}
    let other=StationaryProblem::from_ode(&system,coordinates(&[4.],&[0.1,2.],&[1.,3.]),
        budget(),&action,None,&admit_all).unwrap();
    let mut untouched=[99.;2];
    assert!(other.residual_into(&z,&mut scratch,&mut untouched,||false).is_err());
    assert_eq!(untouched,[99.;2]);
}

#[test]
fn g0_g4_limits_domain_refusals_and_late_cancellation_are_transactional() {
    let system=source();let problem=StationaryProblem::from_ode(&system,
        coordinates(&[3.],&[1.,1.],&[1.,3.]),budget(),&hessian,None,&admit_all).unwrap();
    let original=[0.,0.];let count=Cell::new(0);
    let baseline=solve_equilibrium(&problem,&original,options(),||{count.set(count.get()+1);false}).unwrap();
    let total=count.get();
    for cut in [1,total/2,total] {
        let count=Cell::new(0);
        assert!(solve_equilibrium(&problem,&original,options(),||{
            count.set(count.get()+1);count.get()>=cut}).is_err());
    }
    let replay=solve_equilibrium(&problem,&original,options(),||false).unwrap();
    assert_eq!(baseline.linearization.equilibrium_state(),replay.linearization.equilibrium_state());
    assert_eq!(original,[0.,0.]);
    let mut short=options();short.max_residual_calls=1;
    assert!(matches!(solve_equilibrium(&problem,&original,short,||false),Err(EquilibriumSolveError::Budget{..})));
    short=options();short.max_jacobian_calls=1;
    assert!(matches!(solve_equilibrium(&problem,&original,short,||false),Err(EquilibriumSolveError::Budget{..})));
    short=options();short.max_workspace_bytes=0;
    assert!(solve_equilibrium(&problem,&original,short,||false).is_err());
    let refused=StationaryProblem::from_ode(&system,coordinates(&[3.],&[1.,1.],&[1.,3.]),
        budget(),&hessian,None,&refuse_all).unwrap();
    assert!(solve_equilibrium(&refused,&original,options(),||false).is_err());
}

#[test]
fn g0_nonzero_unbalanced_free_flow_cannot_be_accepted_as_equilibrium() {
    let source=PortHamiltonian::new(1,1,vec![0.],vec![0.],vec![1.],
        Box::new(QuadraticStorage::new(vec![1.],1).unwrap())).unwrap();
    let action=|_:&[f64],d:&[f64],o:&mut[f64]| {o.copy_from_slice(d);true};
    let problem=StationaryProblem::from_ode(&source,coordinates(&[1.],&[1.],&[1.]),
        budget(),&action,None,&admit_all).unwrap();
    assert!(solve_equilibrium(&problem,&[0.],options(),||false).is_err());
}
