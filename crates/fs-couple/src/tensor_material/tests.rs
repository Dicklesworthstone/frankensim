use super::*;
use fs_material::{Tangent6,Voigt};
use fs_material::visco_tensor::{TensorProny,TensorPronyBranch,TensorProjectionBudget};
use fs_math::det;

fn diagonal(v:f64)->Tangent6 {core::array::from_fn(|i|core::array::from_fn(|j|if i==j {v} else {0.0}))}
fn budget()->TensorMaterialBudget {TensorMaterialBudget {max_state_coordinates:96,max_material_points:8}}
fn project(m:&TensorProny,c:&[Voigt],volume:f64)->ProjectedTensorProny {
    m.project_strain(c,volume,TensorProjectionBudget {max_coordinates:4,max_scalar_terms:20_000},||false).unwrap()
}
fn scalar(branch:bool)->TensorMaterialBody {
    let mut b=[[0.0;6];6];b[0][0]=2.0;
    let branches=if branch {vec![TensorPronyBranch::from_gram_rows(b,0.05).unwrap()]} else {vec![]};
    let m=TensorProny::new_mandel(diagonal(12.0),branches,0.2).unwrap();
    TensorMaterialBody::new(vec![2.0],vec![project(&m,&[[1.0,0.0,0.0,0.0,0.0,0.0]],1.0)],
        budget(),&CancelGate::new_clock_free()).unwrap()
}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol,"{a} vs {b}, tolerance {tol}");}

#[test]
fn elastic_limit_uses_exact_oscillator_motion_and_force_velocity_work() {
    let body=scalar(false);let gate=CancelGate::new_clock_free();
    assert_eq!(body.state_dimension(),2);
    let mut h=body.history(&[0.001],&[0.002],InitialMemory::Relaxed,&gate).unwrap();
    let w=det::sqrt(6.0);let force=0.003;let rest=force/12.0;
    for _ in 0..8 {
        let before=h.displacement_m()[0];let r=h.advance(&[force],0.03,&gate).unwrap();
        let t=h.time_s();
        let expected=rest+(0.001-rest)*det::cos(w*t)+0.002/w*det::sin(w*t);
        close(h.displacement_m()[0],expected,2e-13);
        close(r.supplied,force*(h.displacement_m()[0]-before),1e-14);
        close(h.dissipated_heat_j(),0.0,1e-15);close(h.balance_residual_j(),0.0,1e-13);
    }
}

// Independent three-state scalar Maxwell ODE: displacement, velocity, branch
// stress. Used only as a fine-step reference, not a new production integrator.
fn reference(x:[f64;3],f:f64,dt:f64)->[f64;3] {
    let rhs=|x:[f64;3]|[x[1],(f-12.0*x[0]-x[2])/2.0,4.0*x[1]-x[2]/0.05];
    let plus=|a:[f64;3],b:[f64;3],s:f64|core::array::from_fn(|i|a[i]+s*b[i]);
    let a=rhs(x);let b=rhs(plus(x,a,0.5*dt));let c=rhs(plus(x,b,0.5*dt));let d=rhs(plus(x,c,dt));
    core::array::from_fn(|i|x[i]+dt*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0)
}

#[test]
fn full_material_memory_changes_motion_and_matches_independent_maxwell_dynamics() {
    let body=scalar(true);let gate=CancelGate::new_clock_free();assert_eq!(body.state_dimension(),8);
    let mut h=body.history(&[0.001],&[0.002],InitialMemory::Instantaneous,&gate).unwrap();
    let mut x=[0.001,0.002,0.004];
    let elastic=scalar(false);let mut undamped=elastic.history(&[0.001],&[0.002],InitialMemory::Relaxed,&gate).unwrap();
    for tick in 0..12 {
        let force=if tick<3 {0.03} else {0.0};
        h.advance(&[force],0.01,&gate).unwrap();undamped.advance(&[force],0.01,&gate).unwrap();
        for _ in 0..200 {x=reference(x,force,0.01/200.0);}
        close(h.displacement_m()[0],x[0],2e-12);close(h.velocity_m_s()[0],x[1],2e-12);
        close(2.0*h.state()[2],x[2],2e-12); // A=2, resisting force=A*z.
        close(h.balance_residual_j(),0.0,2e-12);
    }
    assert!(h.dissipated_heat_j()>1e-8);
    assert!((h.displacement_m()[0]-undamped.displacement_m()[0]).abs()>1e-6);
}

fn coupled_points(cross:bool,parts:usize)->Vec<ProjectedTensorProny> {
    let mut b=[[0.0;6];6];b[0][0]=2.0;
    if cross {b[0][1]=-1.0;} else {b[1][1]=1.0;}
    let m=TensorProny::new_mandel(diagonal(12.0),
        vec![TensorPronyBranch::from_gram_rows(b,0.1).unwrap()],0.2).unwrap();
    let columns=[[1.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0]];
    (0..parts).map(|_|project(&m,&columns,1.0/parts as f64)).collect()
}

#[test]
fn anisotropic_cross_memory_survives_assembly_and_quadrature_subdivision() {
    let gate=CancelGate::new_clock_free();
    let a=TensorMaterialBody::new(vec![2.0,3.0],coupled_points(true,1),budget(),&gate).unwrap();
    let b=TensorMaterialBody::new(vec![2.0,3.0],coupled_points(true,2),budget(),&gate).unwrap();
    let diagonal=TensorMaterialBody::new(vec![2.0,3.0],coupled_points(false,1),budget(),&gate).unwrap();
    let mut ha=a.history(&[0.001,0.0],&[0.0,0.0],InitialMemory::Instantaneous,&gate).unwrap();
    let mut hb=b.history(&[0.001,0.0],&[0.0,0.0],InitialMemory::Instantaneous,&gate).unwrap();
    let mut hd=diagonal.history(&[0.001,0.0],&[0.0,0.0],InitialMemory::Instantaneous,&gate).unwrap();
    for _ in 0..5 {
        ha.advance(&[0.0,0.0],0.01,&gate).unwrap();hb.advance(&[0.0,0.0],0.01,&gate).unwrap();
        hd.advance(&[0.0,0.0],0.01,&gate).unwrap();
        for i in 0..2 {close(ha.displacement_m()[i],hb.displacement_m()[i],2e-12);}
        close(ha.dissipated_heat_j(),hb.dissipated_heat_j(),2e-12);
    }
    assert!(ha.displacement_m()[1].abs()>1e-8,"a scalar/diagonal loss replacement would erase this motion");
    close(hd.displacement_m()[1],0.0,1e-14);
}

#[test]
fn initial_history_is_explicit_and_failures_preserve_every_accepted_value() {
    let body=scalar(true);let gate=CancelGate::new_clock_free();
    let relaxed=body.history(&[0.001],&[0.0],InitialMemory::Relaxed,&gate).unwrap();
    let mut h=body.history(&[0.001],&[0.0],InitialMemory::Instantaneous,&gate).unwrap();
    close(h.stored_energy_j()-relaxed.stored_energy_j(),0.5*4.0*0.001_f64.powi(2),1e-15);
    let before=h.clone();
    for (force,dt) in [(f64::NAN,0.01),(0.0,-1.0),(0.0,f64::NAN),(1e6,0.1)] {
        assert!(h.advance(&[force],dt,&gate).is_err());
        assert_eq!(h.state(),before.state());assert_eq!(h.time_s(),before.time_s());
        assert_eq!(h.steps(),before.steps());assert_eq!(h.supplied_work_j(),before.supplied_work_j());
        assert_eq!(h.dissipated_heat_j(),before.dissipated_heat_j());
    }
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(h.advance(&[0.0],0.01,&cancelled),Err(TensorMaterialError::Cancelled)));
    let mut retry=before.clone();h.advance(&[0.0],0.01,&gate).unwrap();retry.advance(&[0.0],0.01,&gate).unwrap();
    assert_eq!(h.state(),retry.state());assert_eq!(h.supplied_work_j(),retry.supplied_work_j());
    assert!(body.history(&[0.3],&[0.0],InitialMemory::Relaxed,&gate).is_err());
    assert!(TensorMaterialBody::new(vec![0.0],body.material_points().to_vec(),budget(),&gate).is_err());
    assert!(TensorMaterialBody::new(vec![2.0],body.material_points().to_vec(),
        TensorMaterialBudget {max_state_coordinates:7,..budget()},&gate).is_err());
}
