use super::*;
use super::super::super::{ImpactBody,BodyPotential,ImpactConfig,VolumeSpring,
    damping::ViscousDamper, felt::{FeltPad,KelvinBranch},
    membrane::MembranePotential, radiation::{Model,Pole}};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn config()->ImpactConfig {ImpactConfig {dt_s:2e-6,max_steps:2000,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4}}
fn body(potential:BodyPotential)->ImpactBody {
    let n=potential.count();ImpactBody {potential,initial:vec![ModalAcousticState::default();n],damping_per_s:vec![1.0;n]}
}
fn strike(loss:f64,damped:bool)->ImpactSystem {
    let (stick,w)=ImpactBody::free_mass(0.02,-0.0002,0.8).unwrap();
    let contact=Obstacle::new(vec![w,-2.0],1,2,vec![0.0],vec![1.0],2e7,1.5,
        "momentum/Hertz regression, not measured hardware".into()).unwrap().with_internal_loss(loss).unwrap();
    let pad=FeltPad {weights:vec![0.0,2.0],area_m2:0.001,thickness_m:0.006,precompression_m:0.0006,
        law:fs_material::fiber::WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),
        prior_maximum_strain:0.1,creep:vec![KelvinBranch {stiffness_n_m:1500.0,viscosity_n_s_m:8.0}]};
    let dampers=if damped {vec![ViscousDamper {weights:vec![w,-2.0],damping_n_s_m:0.1}]}else{vec![]};
    ImpactSystem::new_with_dampers(vec![stick,body(BodyPotential::Linear(vec![core::f64::consts::TAU*500.0]))],
        vec![contact],vec![pad],vec![],dampers,config()).unwrap()
}

#[test]
fn hertz_rebound_and_felt_memory_keep_the_full_dense_trajectory_and_energy_ledger() {
    let mut fast=strike(0.0,false).prepare_analytic().unwrap();
    let mut dense=strike(0.0,false).prepare_analytic().unwrap();dense.set_momentum_condensation(false).unwrap();
    assert_eq!(fast.state().len(),5);assert_eq!(fast.condensed_newton_dimension(),Some(4));
    assert_eq!(dense.condensed_newton_dimension(),None);
    let gate=CancelGate::new_clock_free();let initial=fast.stored_energy_j();let mut net=0.0;
    let mut rebound=false;let mut solves=0;
    for tick in 0..1200 {
        let force=[if tick<50 {0.1}else{0.0},0.0];
        let a=fast.step(&force,&gate).unwrap();let b=dense.step(&force,&gate).unwrap();
        assert_eq!(a.sample,b.sample);assert_eq!(a.time_s,b.time_s);
        for (a,b) in fast.state().iter().zip(dense.state()) {assert!((a-b).abs()<2e-8);}
        assert!(a.balance_residual_j.abs()<1e-8);net+=a.supplied_work_j-a.dissipated_energy_j;
        rebound|=fast.state()[1]<0.0;solves+=fast.newton_linear_solve_counts().0;
        assert_eq!(fast.newton_linear_solve_counts().1,0);
    }
    assert!(rebound && solves>0 && fast.state()[4]!=0.0);
    assert!((fast.stored_energy_j()-initial-net).abs()<1e-7);
    assert!((fast.felt_observation(0).unwrap().0-dense.felt_observation(0).unwrap().0).abs()<1e-7);
}

fn membrane()->BodyPotential {
    use fs_plate::shell::head::{TensionedDisk,TensionedDiskSpec,nonlinear::MembraneReductionBudget};
    use fs_plate::shell::profile::ProfileBudget;
    let d=TensionedDisk::new(TensionedDiskSpec {radius_m:0.17,thickness_m:0.000254,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:3000.0,radial_intervals:2,azimuths:8},
        ProfileBudget {max_nodes:100,max_triangles:100,max_feature_evaluations:0}).unwrap();
    let n=d.model.free;let mut k=vec![0.0;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=d.model.k.get(i,j);m[i*n+j]=d.model.m.get(i,j);}}
    let mut modes=fs_modal::eigh_gen_dense(&k,&m,n).unwrap();modes.truncate(3);
    let rim:Vec<_>=(0..d.mesh.nodes.len()).filter(|&i|d.model.dof_map[3*i].is_none()).collect();
    BodyPotential::Membrane(MembranePotential::from_pencil(&d.mesh,&d.section,&d.model,&modes,&rim,
        MembraneReductionBudget {max_modes:4,max_nodes:100,max_facet_pairs:1000,max_solve_entries:100000,
            relative_tolerance:1e-6},0.25).unwrap())
}
fn shell()->BodyPotential {
    use fs_plate::shell::{ShellMesh,ShellSupport,assemble_shell_sections};
    use fs_plate::shell::reduction::{ShellReduction,ReductionBudget};
    let mesh=ShellMesh::new(vec![[0.0,0.0,0.0],[0.1,0.0,0.01],[0.0,0.1,0.02],[0.08,0.09,0.035]],
        vec![[0,1,3],[0,3,2]]).unwrap();
    let sections=vec![fs_plate::PlateSection::isotropic(112.6e9,0.342,0.001,8607.0).unwrap(),
        fs_plate::PlateSection::isotropic(112.6e9,0.342,0.0007,8607.0).unwrap()];
    let model=assemble_shell_sections(&mesh,&sections,&[0,1,2],ShellSupport::Clamped).unwrap();
    let n=model.free;let mut k=vec![0.0;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=model.k.get(i,j);m[i*n+j]=model.m.get(i,j);}}
    let modes=fs_modal::eigh_gen_dense(&k,&m,n).unwrap();
    BodyPotential::Shell(ShellReduction::new(&mesh,&sections,&model,&modes,
        ReductionBudget {max_modes:12,max_facet_modes:100,relative_tolerance:1e-6}).unwrap())
}
fn geometric(curved:bool)->ImpactSystem {
    let mut b=body(if curved {shell()}else{membrane()});
    for (i,state) in b.initial.iter_mut().enumerate() {
        state.displacement_m_sqrt_kg=if curved {1e-7}else{1e-5}*((i+1) as f64).sin();
        state.velocity_m_sqrt_kg_per_s=1e-4*(i+1) as f64;
    }
    let mut cfg=config();if curved {cfg.dt_s=1e-9;}
    // This is a reciprocal gas-volume row, not a frequency tweak or extra damping.
    let volume=VolumeSpring {bulk_modulus_pa:1.4e5,volume_m3:0.01,
        areas:(0..b.initial.len()).map(|i|0.01*if i%2==0 {1.0}else{-1.0}).collect()};
    ImpactSystem::new(vec![b],vec![],vec![],vec![volume],cfg).unwrap()
}

#[test]
fn actual_curved_shell_and_stretching_head_keep_nonlinear_coupling_and_all_coordinates() {
    let gate=CancelGate::new_clock_free();
    for curved in [false,true] {
        let mut fast=geometric(curved).prepare_analytic().unwrap();
        let mut dense=geometric(curved).prepare_analytic().unwrap();dense.set_momentum_condensation(false).unwrap();
        let modes=fast.inner.modes;let dimension=fast.state().len();let mut solves=0;
        assert!(modes>1);assert_eq!(fast.condensed_newton_dimension(),Some(dimension-modes+1));
        let force=vec![0.0;modes];
        for _ in 0..30 {
            let f=fast.step(&force,&gate).unwrap();dense.step(&force,&gate).unwrap();
            assert_eq!(fast.state().len(),dimension);assert!(f.balance_residual_j.abs()<1e-8);
            for (a,b) in fast.state().iter().zip(dense.state()) {assert!((a-b).abs()<1e-8*a.abs().max(b.abs()).max(1.0));}
            solves+=fast.newton_linear_solve_counts().0;assert_eq!(fast.newton_linear_solve_counts().1,0);
        }
        assert!(solves>0);
    }
}

#[test]
fn cross_velocity_and_nonlinear_contact_loss_keep_their_original_solve_by_default() {
    let gate=CancelGate::new_clock_free();
    for (loss,damped) in [(0.0,true),(0.1,false)] {
        let mut dense=strike(loss,damped).prepare_analytic().unwrap();
        assert_eq!(dense.condensed_newton_dimension(),None);
        let mut requested=strike(loss,damped).prepare_analytic().unwrap();
        requested.set_momentum_condensation(true).unwrap();let mut fallbacks=0;
        for _ in 0..250 {
            let f=dense.step(&[0.0;2],&gate).unwrap();let r=requested.step(&[0.0;2],&gate).unwrap();
            assert_eq!(f.sample,r.sample);
            for (a,b) in dense.state().iter().zip(requested.state()) {assert!((a-b).abs()<2e-8);}
            fallbacks+=requested.newton_linear_solve_counts().1;
        }
        assert!(fallbacks>0,"a spatial or active impact resistance must not be diagonalized");
    }
}

#[test]
fn refused_steps_switching_and_repreparation_do_not_reset_felt_or_motion() {
    let gate=CancelGate::new_clock_free();let mut retry=strike(0.0,false).prepare_analytic().unwrap();
    let mut clean=strike(0.0,false).prepare_analytic().unwrap();
    for _ in 0..140 {retry.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();}
    let state=retry.state().to_vec();let history=retry.felt_history(0);let samples=retry.samples();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(retry.step(&[0.0;2],&cancelled),Err(ImpactError::Cancelled)));
    assert!(retry.step(&[f64::NAN,0.0],&gate).is_err());
    retry.set_iteration_limit(0).unwrap();assert!(retry.step(&[0.0;2],&gate).is_err());
    assert_eq!(retry.state(),state);assert_eq!(retry.felt_history(0),history);assert_eq!(retry.samples(),samples);
    retry.set_iteration_limit(50).unwrap();retry.set_momentum_condensation(false).unwrap();
    assert_eq!(retry.state(),state);retry.set_momentum_condensation(true).unwrap();
    retry.set_analytic_newton(false);assert_eq!(retry.condensed_newton_dimension(),None);
    retry.set_analytic_newton(true);retry.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();
    assert_eq!(retry.state(),clean.state());assert_eq!(retry.felt_history(0),clean.felt_history(0));
    let mut retry=retry.into_reference().prepare_analytic().unwrap();
    retry.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();
    assert_eq!(retry.state(),clean.state());assert_eq!(retry.samples(),clean.samples());
}

#[test]
fn radiation_pair_preparation_is_unchanged_and_cannot_be_accidentally_replaced() {
    let model=Model {ports:2,poles:(0..4).map(|i|Pole {omega:300.0+30.0*i as f64,zeta:0.1,
        coupling:vec![0.4,-0.2]}).collect()};
    let mut system=strike(0.0,false).with_radiation_load(&model,&[0,1],4).unwrap().prepare_analytic().unwrap();
    let state=system.state().to_vec();let history=system.felt_history(0);
    let dimension=system.condensed_newton_dimension();assert_eq!(dimension,Some(6));
    for enabled in [true,false] {
        assert!(system.set_momentum_condensation(enabled).is_err());
        assert_eq!(system.condensed_newton_dimension(),dimension);assert_eq!(system.state(),state);
        assert_eq!(system.felt_history(0),history);
    }
    system.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
    assert_eq!(system.newton_linear_solve_counts().1,0);
}
