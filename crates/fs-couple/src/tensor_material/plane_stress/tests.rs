use super::*;
use fs_material::{IsotropicElastic, Tangent6};
use fs_material::visco_tensor::TensorPronyBranch;
use fs_math::c64::C64;

fn law()->TensorProny {
    let eq=IsotropicElastic::new(100.0,0.27,0.1).unwrap().stiffness();
    let mut a:Tangent6=[[0.0;6];6];a[0]=[2.0,-1.0,3.0,0.7,1.0,0.2];
    let mut b:Tangent6=[[0.0;6];6];b[1]=[-1.0,2.0,0.5,1.0,-0.4,2.0];
    TensorProny::from_tensor_shear(eq,vec![TensorPronyBranch::from_gram_rows(a,0.1).unwrap(),
        TensorPronyBranch::from_gram_rows(b,0.7).unwrap()],0.1).unwrap()
}
fn budget()->TensorProjectionBudget {TensorProjectionBudget{max_coordinates:16,max_scalar_terms:1_000_000}}
fn body_budget()->TensorMaterialBudget {TensorMaterialBudget{max_state_coordinates:128,max_material_points:4}}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol*(1.0+a.abs().max(b.abs())),"{a} != {b}");}

#[test]
fn retained_plane_stress_memory_matches_the_full_complex_schur_complement() {
    let gate=CancelGate::new_clock_free();let law=law();
    let point=PlaneStressPoint::new(&law,&[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],0.2,budget(),&gate).unwrap();
    let n=3;let b=point.memory_coordinates();let m=n+b;
    assert_eq!(b,12,"both rank-deficient branches retain six coordinates");
    for w in [0.0,0.1,1.0,10.0,1000.0] {
        // a'=-D(Q_au u+Q_aa a); harmonic elimination is independent of
        // the source owner's frequency-dependent 3D plane-stress reduction.
        let matrix:Vec<_>=(0..b*b).map(|k|C64::new(
            point.rates[k/b]*point.energy[(n+k/b)*m+n+k%b],
            if k/b==k%b {w} else {0.0})).collect();
        let lu=fs_la::eigen_complex::lu_complex(&matrix,b).unwrap();
        let expected=law.plane_stress(w).unwrap().engineering();
        for col in 0..n {
            let mut a:Vec<_>=(0..b).map(|i|C64::from_re(-point.rates[i]*point.energy[(n+i)*m+col])).collect();
            lu.solve(&mut a);
            for row in 0..n {
                let actual=(0..b).fold(C64::from_re(point.energy[row*m+col]),|s,i|
                    s+a[i].scale(point.energy[row*m+n+i]));
                close(actual.re,0.2*expected[row][col].re,2e-10);
                close(actual.im,0.2*expected[row][col].im,2e-10);
            }
        }
    }
}

#[test]
fn force_driven_plane_stress_reconstructs_traction_and_closes_history_energy() {
    let gate=CancelGate::new_clock_free();let law=law();
    let columns=[[1.0,0.0,0.2],[0.0,1.0,-0.3]];
    let point=PlaneStressPoint::new(&law,&columns,0.1,budget(),&gate).unwrap();
    let body=TensorMaterialBody::from_plane_stress(vec![1.0,1.5],vec![point],body_budget(),&gate).unwrap();
    let point=&body.plane_stress_points()[0];
    assert!(body.material_points().is_empty());
    assert_eq!(body.state_dimension(),16);
    for initial in [InitialMemory::Relaxed,InitialMemory::Instantaneous] {
        let mut history=body.history(&[0.001,-0.002],&[0.0,0.0],initial,&gate).unwrap();
        let mut changed=false;
        for tick in 0..24 {
            let force=if tick<5 {[0.04,0.0]} else {[0.0,0.0]};
            let rec=history.advance(&force,0.01,&gate).unwrap();
            let memory=&history.state()[4..];
            let strain=history.plane_stress_strain(0).unwrap();
            assert_eq!(strain,point.strain(history.displacement_m(),memory).unwrap());
            let instant=point.strain(history.displacement_m(),&[0.0;12]).unwrap();
            changed|=(strain[2]-instant[2]).abs()>1e-9;
            let traction=history.plane_stress_normal_stress_pa(0).unwrap();
            assert!(traction.iter().all(|s|s.abs()<2e-11),"{traction:?}");
            assert!(rec.dissipated>=-1e-13);
            assert!(history.balance_residual_j().abs()<2e-11);
        }
        assert!(changed,"normal strain must track material memory, not a fixed Poisson map");
        assert!(history.dissipated_heat_j()>0.0);
    }
}

#[test]
fn elastic_plane_stress_and_quadrature_subdivision_keep_physical_motion() {
    let gate=CancelGate::new_clock_free();
    let material=TensorProny::from_tensor_shear(IsotropicElastic::new(120.0,0.25,0.1).unwrap().stiffness(),vec![],0.1).unwrap();
    let columns=[[2.0,0.0,0.0]];
    let full=PlaneStressPoint::new(&material,&columns,0.3,budget(),&gate).unwrap();
    let half=PlaneStressPoint::new(&material,&columns,0.15,budget(),&gate).unwrap();
    let stiffness=0.3*4.0*120.0/(1.0-0.25*0.25);
    close(full.energy_matrix()[0],stiffness,1e-12);
    let one=TensorMaterialBody::from_plane_stress(vec![2.0],vec![full],body_budget(),&gate).unwrap();
    let two=TensorMaterialBody::from_plane_stress(vec![2.0],vec![half.clone(),half],body_budget(),&gate).unwrap();
    let mut a=one.history(&[0.002],&[0.0],InitialMemory::Relaxed,&gate).unwrap();
    let mut b=two.history(&[0.002],&[0.0],InitialMemory::Relaxed,&gate).unwrap();
    for tick in 1..=12 {
        a.advance(&[0.0],0.01,&gate).unwrap();b.advance(&[0.0],0.01,&gate).unwrap();
        close(a.displacement_m()[0],0.002*det::cos(det::sqrt(stiffness/2.0)*tick as f64*0.01),1e-12);
        close(a.displacement_m()[0],b.displacement_m()[0],1e-12);
        close(a.dissipated_heat_j(),0.0,1e-12);
    }
}

#[test]
fn plane_stress_refuses_unsupported_rate_constraints_and_keeps_failed_history() {
    let gate=CancelGate::new_clock_free();let law=law();let columns=[[1.0,0.0,0.0]];
    assert!(PlaneStressPoint::new(&law,&[],1.0,budget(),&gate).is_err());
    assert!(PlaneStressPoint::new(&law,&columns,0.0,budget(),&gate).is_err());
    assert!(PlaneStressPoint::new(&law,&columns,1.0,TensorProjectionBudget{max_scalar_terms:1,..budget()},&gate).is_err());
    let kv=TensorProny::new(*law.c_eq(),Some(core::array::from_fn(|i|core::array::from_fn(|j|
        if i==j {1.0} else {0.0}))),vec![],(0.0,f64::MAX)).unwrap();
    assert!(PlaneStressPoint::new(&kv,&columns,1.0,budget(),&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(PlaneStressPoint::new(&law,&columns,1.0,budget(),&cancelled).is_err());
    let point=PlaneStressPoint::new(&law,&columns,1.0,budget(),&gate).unwrap();
    assert!(TensorMaterialBody::from_plane_stress(vec![1.0],vec![point.clone()],
        TensorMaterialBudget{max_state_coordinates:2,..body_budget()},&gate).is_err());
    let body=TensorMaterialBody::from_plane_stress(vec![1.0],vec![point],body_budget(),&gate).unwrap();
    let mut history=body.history(&[0.001],&[0.0],InitialMemory::Relaxed,&gate).unwrap();
    assert!(history.plane_stress_strain(1).is_err());
    let before=history.state().to_vec();
    for (force,dt,g) in [(vec![f64::NAN],0.01,&gate),(vec![0.0],0.01,&cancelled),(vec![0.0],-0.01,&gate)] {
        assert!(history.advance(&force,dt,g).is_err());
        assert_eq!(history.state(),before);
        assert_eq!(history.steps(),0);assert_eq!(history.supplied_work_j(),0.0);
    }
    assert!(body.history(&[1.0],&[0.0],InitialMemory::Relaxed,&gate).is_err());
    history.advance(&[0.0],0.01,&gate).unwrap();assert_eq!(history.steps(),1);
}


#[test]
fn consistent_plane_stress_matches_independent_physical_velocity_and_memory_dynamics() {
    let gate=CancelGate::new_clock_free();
    let point=PlaneStressPoint::new(&law(),&[[1.0,0.0,0.2],[0.0,1.0,-0.3]],
        0.1,budget(),&gate).unwrap();
    // Full M and anisotropic K do not share a coordinate basis.
    let mass=vec![1.4,0.45,0.45,0.9];
    let body=TensorMaterialBody::from_plane_stress_mass_matrix(2,mass.clone(),
        vec![point],body_budget(),&gate).unwrap();
    assert_eq!(body.mass_matrix_kg(),mass);
    let point=&body.plane_stress_points()[0];
    let determinant=mass[0]*mass[3]-mass[1]*mass[2];
    let inverse=[mass[3]/determinant,-mass[1]/determinant,
        -mass[2]/determinant,mass[0]/determinant];
    for memory in [InitialMemory::Relaxed,InitialMemory::Instantaneous] {
        let velocity=[0.003,-0.002];
        let mut history=body.history(&[0.001,-0.002],&velocity,memory,&gate).unwrap();
        for i in 0..2 {close(history.state()[2+i],mass[2*i]*velocity[0]+mass[2*i+1]*velocity[1],1e-15);}
        let rest=body.history(&[0.001,-0.002],&[0.0;2],memory,&gate).unwrap();
        close(history.stored_energy_j()-rest.stored_energy_j(),
            0.5*(mass[0]*velocity[0]*velocity[0]+2.0*mass[1]*velocity[0]*velocity[1]
                +mass[3]*velocity[1]*velocity[1]),1e-15);
        // Independent RK4 in physical (u,v,a), not the pH generator, augmented
        // with direct force power and each original branch's dissipated power.
        let state_len=history.state().len();
        let mut reference=history.state().to_vec();reference[2..4].copy_from_slice(&velocity);
        reference.extend([0.0,0.0]);
        for tick in 0..12 {
            let force=if tick<4 {[0.03,-0.02]} else {[0.0;2]};
            let rhs=|x:&[f64]| {
                let mut derivative=vec![0.0;x.len()];
                derivative[..2].copy_from_slice(&x[2..4]);
                let local:Vec<_>=x[..2].iter().chain(&x[4..state_len]).copied().collect();
                let d=local.len();
                let gradient:Vec<_>=(0..d).map(|i|(0..d).map(|j|point.energy[i*d+j]*local[j]).sum::<f64>()).collect();
                for i in 0..2 {derivative[2+i]=inverse[2*i]*(force[0]-gradient[0])
                    +inverse[2*i+1]*(force[1]-gradient[1]);}
                for i in 0..point.memory_coordinates() {
                    derivative[4+i]=-point.rates[i]*gradient[2+i];
                    derivative[state_len+1]+=point.rates[i]*gradient[2+i]*gradient[2+i];
                }
                derivative[state_len]=force[0]*x[2]+force[1]*x[3];
                derivative
            };
            let dt=0.002;let h=dt/100.0;
            for _ in 0..100 {
                let k1=rhs(&reference);
                let stage=|k:&[f64],scale:f64|reference.iter().zip(k).map(|(x,k)|x+scale*k).collect::<Vec<_>>();
                let k2=rhs(&stage(&k1,0.5*h));let k3=rhs(&stage(&k2,0.5*h));let k4=rhs(&stage(&k3,h));
                for i in 0..reference.len() {reference[i]+=h*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0;}
            }
            history.advance(&force,dt,&gate).unwrap();
            for (got,want) in history.displacement_m().iter().chain(&history.velocity_m_s())
                .chain(&history.state()[4..]).zip(&reference[..state_len]) {close(*got,*want,2e-11);}
            close(history.supplied_work_j(),reference[state_len],2e-11);
            close(history.dissipated_heat_j(),reference[state_len+1],2e-11);
            close(history.balance_residual_j(),0.0,2e-11);
            for traction in history.plane_stress_normal_stress_pa(0).unwrap() {close(traction,0.0,2e-11);}
        }
        assert!(history.dissipated_heat_j()>0.0);
    }
}

#[test]
fn full_plane_stress_inertia_preserves_diagonal_motion_and_refuses_invalid_mass() {
    let gate=CancelGate::new_clock_free();
    let point=PlaneStressPoint::new(&law(),&[[1.0,0.0,0.2],[0.0,1.0,-0.3]],
        0.1,budget(),&gate).unwrap();
    for mass in [vec![1.0],vec![1.0,1.0,1.0,1.0],vec![1.0,2.0,2.0,1.0],
        vec![1.0,0.1,0.2,1.0],vec![1.0,0.0,0.0,f64::NAN]] {
        assert!(TensorMaterialBody::from_plane_stress_mass_matrix(2,mass,vec![point.clone()],body_budget(),&gate).is_err());
    }
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(TensorMaterialBody::from_plane_stress_mass_matrix(2,vec![1.0,0.0,0.0,2.0],
        vec![point.clone()],body_budget(),&cancelled),Err(TensorMaterialError::Cancelled)));
    assert!(TensorMaterialBody::from_plane_stress_mass_matrix(2,vec![1.0,0.0,0.0,2.0],
        vec![point.clone()],TensorMaterialBudget{max_state_coordinates:4,..body_budget()},&gate).is_err());
    let a=TensorMaterialBody::from_plane_stress(vec![1.0,2.0],vec![point.clone()],body_budget(),&gate).unwrap();
    let b=TensorMaterialBody::from_plane_stress_mass_matrix(2,vec![1.0,0.0,0.0,2.0],vec![point],body_budget(),&gate).unwrap();
    let mut x=a.history(&[0.001,-0.002],&[0.002,0.003],InitialMemory::Relaxed,&gate).unwrap();
    let mut y=b.history(&[0.001,-0.002],&[0.002,0.003],InitialMemory::Relaxed,&gate).unwrap();
    x.advance(&[0.01,0.02],0.03,&gate).unwrap();y.advance(&[0.01,0.02],0.03,&gate).unwrap();
    for (x,y) in x.state().iter().zip(y.state()) {close(*x,*y,2e-13);}
}
