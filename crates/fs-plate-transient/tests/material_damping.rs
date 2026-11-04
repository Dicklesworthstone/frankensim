//! Spatial constitutive loss: dense stepping, physical work and transposed actions.
#[path = "../examples/moving_plate.rs"]
#[allow(dead_code)]
mod moving;

use fs_plate::{PlateMesh, PlateModel, viscosity::{
    PlateBendingViscosity, PlateViscosityBudget, ViscousPlateSection,
}};
use fs_plate_transient::{NoPlateLoad, PlateDynamics, PlateDynamicsError};
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::{
    GeneralizedAlpha, SecondOrderOperatorWeights, SecondOrderProblem, SecondOrderState, galpha_step,
    second_order_adjoint::SecondOrderVjp,
};

fn loss<'a>(mesh: &PlateMesh, plate: &'a PlateModel) -> PlateBendingViscosity<'a> {
    let tensor = [0.04,0.01,0.005,0.01,0.08,-0.0075,0.005,-0.0075,0.02];
    PlateBendingViscosity::assemble(mesh, plate,
        vec![ViscousPlateSection::from_resultant(tensor).unwrap();mesh.tris.len()],
        PlateViscosityBudget {
            max_elements: mesh.tris.len(), max_dofs: plate.free,
            max_full_dofs: plate.dof_map.len(), max_assembly_entries:81*mesh.tris.len(),
        }, &mut || false).unwrap()
}

fn dynamics<'a>(plate: &'a PlateModel, law: &'a PlateBendingViscosity<'a>, p:[f64;5]) -> PlateDynamics<'a, NoPlateLoad> {
    PlateDynamics::new(plate, moving::dynamics_parameters(p), &NoPlateLoad,
        moving::dynamics_budget(plate), &mut || false).unwrap()
        .with_bending_viscosity(law,plate.m.nnz()+plate.k.nnz()+law.matrix().nnz(),&mut || false).unwrap()
}

fn dot(a:&[f64], b:&[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }

fn close(a:&[f64], b:&[f64], tolerance:f64) {
    let scale = b.iter().fold(1e-12_f64, |scale,v| scale.max(v.abs()));
    let error = a.iter().zip(b).fold(0.0_f64,|error,(a,b)|error.max((a-b).abs()));
    assert!(error <= tolerance*scale,"error {error:e}, scale {scale:e}");
}

fn dense_operators(plate:&PlateModel, law:&PlateBendingViscosity<'_>, p:[f64;5]) -> (Vec<f64>,Vec<f64>,Vec<f64>) {
    let m:Vec<_> = plate.m.to_dense().iter().map(|v| p[1]*v).collect();
    let k:Vec<_> = plate.k.to_dense().iter().map(|v| p[0]*v).collect();
    let c = m.iter().zip(&k).zip(law.matrix().to_dense()).map(|((m,k),c)|p[2]*m+p[3]*k+c).collect();
    (m,c,k)
}

fn dense_action(matrix:&[f64], input:&[f64]) -> Vec<f64> {
    matrix.chunks_exact(input.len()).map(|row|dot(row,input)).collect()
}

#[test]
fn full_material_coupling_reaches_primal_transpose_and_both_preconditioners() {
    let (mesh,plate) = moving::build_plate().unwrap();
    let law = loss(&mesh,&plate);
    let p = [1.1,0.9,3.0,2e-4,0.0];
    let dynamics = dynamics(&plate,&law,p);
    let (m,c,k) = dense_operators(&plate,&law,p);
    let n = plate.free;
    let v:Vec<_>=(0..n).map(|i|(0.71*i as f64).sin()).collect();
    let seed:Vec<_>=(0..n).map(|i|(0.37*i as f64).cos()).collect();
    let mut force = vec![0.0;n];
    dynamics.damping_apply(&v,&mut force);
    close(&force,&dense_action(&c,&v),2e-15);
    let mut transpose=vec![0.0;n];
    dynamics.damping_transpose_apply(&seed,&mut transpose).unwrap();
    let expected:Vec<_>=(0..n).map(|j|(0..n).map(|i|c[n*i+j]*seed[i]).sum()).collect();
    close(&transpose,&expected,2e-15);
    assert!((dot(&seed,&force)-dot(&v,&transpose)).abs()<1e-12*dot(&seed,&force).abs());
    let rho=0.6_f64;
    let alpha_m=(2.0*rho-1.0)/(1.0+rho);
    let alpha_f=rho/(1.0+rho);
    let gamma=0.5-alpha_m+alpha_f;
    let beta=0.25*(1.0-alpha_m+alpha_f).powi(2);
    let weights = SecondOrderOperatorWeights {
        mass:(1.0-alpha_m)/(beta*moving::STEP.powi(2)),
        damping:(1.0-alpha_f)*gamma/(beta*moving::STEP), tangent:1.0-alpha_f,
    };
    let expected:Vec<_>=(0..n).map(|i|seed[i]/(weights.mass*m[n*i+i]+weights.damping*c[n*i+i]+weights.tangent*k[n*i+i])).collect();
    dynamics.preconditioner_apply(&v,weights,0,0,&seed,&mut force);
    close(&force,&expected,2e-15);
    dynamics.effective_preconditioner(moving::STEP,rho,&mut || false).unwrap().apply(0,&seed,&mut force);
    close(&force,&expected,2e-15);
}

#[test]
fn material_loss_is_fixed_under_the_four_mechanical_parameter_derivatives() {
    let (mesh,plate) = moving::build_plate().unwrap();
    let law = loss(&mesh,&plate);
    let p = [1.1,0.9,3.0,2e-4,0.0];
    let n=plate.free;
    let q:Vec<_>=(0..n).map(|i|1e-3*(i as f64*0.37).sin()).collect();
    let v:Vec<_>=(0..n).map(|i|(i as f64*0.71).sin()).collect();
    let a:Vec<_>=(0..n).map(|i|(i as f64*0.23).cos()).collect();
    let seed:Vec<_>=(0..n).map(|i|(i as f64*0.41).cos()).collect();
    let mut derivative=[0.0;4];
    dynamics(&plate,&law,p).residual_parameter_vjp(&q,&v,&a,&seed,&mut derivative).unwrap();
    for (j,h) in [2e-5,2e-5,2e-3,2e-8].into_iter().enumerate() {
        let evaluate=|p| {
            let (m,c,k)=dense_operators(&plate,&law,p);
            dot(&seed,&dense_action(&m,&a))+dot(&seed,&dense_action(&c,&v))+dot(&seed,&dense_action(&k,&q))
        };
        let (mut plus,mut minus)=(p,p); plus[j]+=h;minus[j]-=h;
        let fd=(evaluate(plus)-evaluate(minus))/(2.0*h);
        assert!((derivative[j]-fd).abs()<3e-7*fd.abs(),"parameter {j}: {} != {fd}",derivative[j]);
    }
}

#[test]
fn midpoint_transient_matches_dense_steps_and_exact_material_work_balance() {
    let (mesh,plate)=moving::build_plate().unwrap();
    let law=loss(&mesh,&plate);
    // The sole dissipation is the anisotropic constitutive matrix.
    let p=[1.1,0.9,0.0,0.0,0.0];
    let model=dynamics(&plate,&law,p);
    let (m,c,k)=dense_operators(&plate,&law,p);
    let n=plate.free;
    let q:Vec<_>=(0..n).map(|i|1e-6*(0.37*i as f64).sin()).collect();
    let zero=vec![0.0;n];
    let force=dense_action(&k,&q);
    let acceleration:Vec<_>=(0..n).map(|i|-force[i]/m[n*i+i]).collect();
    let mut actual=SecondOrderState::new(0.0,&q,&zero,&acceleration);
    let mut expected=actual.clone();
    let method=moving::method(n,1.0);
    let dense=GeneralizedAlpha::new(&m,&c,&k,n,moving::STEP,1.0);
    let energy=|state:&SecondOrderState|0.5*(dot(&state.q,&dense_action(&k,&state.q))+dot(&state.v,&dense_action(&m,&state.v)));
    let initial=energy(&actual);
    let mut total_loss=0.0;
    for _ in 0..48 {
        let old=actual.clone();
        let previous=energy(&old);
        method.step(&mut actual,&model,&zero).unwrap();
        galpha_step(&dense,&mut expected.q,&mut expected.v,&mut expected.a,&zero);
        close(&actual.q,&expected.q,2e-8);
        close(&actual.v,&expected.v,2e-8);
        close(&actual.a,&expected.a,2e-8);
        let midpoint:Vec<_>=old.v.iter().zip(&actual.v).map(|(a,b)|0.5*(a+b)).collect();
        let dissipated=moving::STEP*dot(&midpoint,&dense_action(&c,&midpoint));
        assert!(dissipated>=0.0);
        assert!((previous-energy(&actual)-dissipated).abs()<3e-9*initial);
        total_loss+=dissipated;
    }
    assert!(total_loss>0.01*initial);
    assert!((initial-energy(&actual)-total_loss).abs()<1e-8*initial);
}

#[test]
fn original_model_identity_cancellation_and_additional_operator_budget_are_enforced() {
    let (mesh,plate)=moving::build_plate().unwrap();
    let law=loss(&mesh,&plate);
    let p=[1.0,1.0,0.0,0.0,0.0];
    let plain=PlateDynamics::new(&plate,moving::dynamics_parameters(p),&NoPlateLoad,
        moving::dynamics_budget(&plate),&mut || false).unwrap();
    assert!(plain.bending_viscosity().is_none());
    let cap=plate.k.nnz()+plate.m.nnz()+law.matrix().nnz();
    assert!(matches!(plain.with_bending_viscosity(&law,cap-1,&mut || false),Err(PlateDynamicsError::Budget(_))));
    assert!(matches!(plain.with_bending_viscosity(&law,cap,&mut || true),Err(PlateDynamicsError::Cancelled)));
    let different=plate.clone();
    let different_dynamics=PlateDynamics::new(&different,moving::dynamics_parameters(p),&NoPlateLoad,
        moving::dynamics_budget(&different),&mut || false).unwrap();
    assert!(matches!(different_dynamics.with_bending_viscosity(&law,cap,&mut || false),Err(PlateDynamicsError::InvalidInput(_))));
}
