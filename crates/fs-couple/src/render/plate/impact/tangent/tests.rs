//! Real owner assembly, not a fitted Duffing surrogate for shell or film.
use super::*;
use super::super::{ImpactBody,ImpactConfig,ImpactError,VolumeSpring,membrane::MembranePotential};
use super::super::felt::{FeltPad,KelvinBranch};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_dcontact::{ContactStorage,Obstacle};
use fs_exec::CancelGate;
use fs_material::fiber::WoolFelt;

fn config()->ImpactConfig {ImpactConfig {dt_s:2e-6,max_steps:2000,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4}}
fn body(potential:BodyPotential)->ImpactBody {
    let n=potential.count();ImpactBody {potential,initial:vec![ModalAcousticState::default();n],damping_per_s:vec![1.0;n]}
}
fn pad(weights:Vec<f64>,compression:f64)->FeltPad {
    FeltPad {area_m2:0.001,thickness_m:0.006,precompression_m:compression,weights,
        law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),prior_maximum_strain:0.1,
        creep:vec![KelvinBranch {stiffness_n_m:1500.0,viscosity_n_s_m:8.0}]}
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
fn check_tangent(s:&ImpactSystem,x:&[f64]) {
    check_tangent_at(s,x,1e-9,false);
}
fn check_tangent_at(s:&ImpactSystem,x:&[f64],h:f64,cubic_gradient:bool) {
    let n=x.len();
    let mut matrix=vec![0.0;n*n];
    for col in 0..n {
        let mut d=vec![0.0;n];d[col]=1.0;let mut hd=vec![0.0;n];
        assert!(s.hessian_vector(x,&d,&mut hd));
        let mut plus=x.to_vec();plus[col]+=h;let mut minus=x.to_vec();minus[col]-=h;
        let mut gp=vec![0.0;n];let mut gm=gp.clone();
        s.contact.gradient(&plus,&mut gp);s.contact.gradient(&minus,&mut gm);
        let mut gp2=vec![0.0;n];let mut gm2=vec![0.0;n];
        if cubic_gradient {
            plus[col]=x[col]+2.0*h;minus[col]=x[col]-2.0*h;
            s.contact.gradient(&plus,&mut gp2);s.contact.gradient(&minus,&mut gm2);
        }
        for row in 0..n {
            matrix[row*n+col]=hd[row];
            let fd=if cubic_gradient { (8.0*(gp[row]-gm[row])-(gp2[row]-gm2[row]))/(12.0*h) }
                else { (gp[row]-gm[row])/(2.0*h) };
            assert!((fd-hd[row]).abs()<1e-5*fd.abs().max(hd[row].abs()).max(1.0),
                "({row},{col}): {fd:e} != {:e}",hd[row]);
        }
    }
    for i in 0..n {for j in 0..n {
        assert!((matrix[i*n+j]-matrix[j*n+i]).abs()<1e-11*matrix[i*n+j].abs().max(matrix[j*n+i].abs()).max(1.0));
    }}
    assert!(!s.hessian_vector(&[],&[],&mut []));
    let mut bad=x.to_vec();bad[0]=f64::NAN;
    assert!(!s.hessian_vector(&bad,&vec![0.0;n],&mut vec![0.0;n]));
}
#[test]
fn exact_curved_shell_tangent_includes_geometric_and_material_stiffness() {
    let b=body(shell());let n=b.potential.count();
    let mut cfg=config();cfg.dt_s=1e-9; // Retain every test-pencil mode, no Nyquist bypass.
    let s=ImpactSystem::new(vec![b],vec![],vec![],vec![],cfg).unwrap();
    let x:Vec<_>=(0..2*n).map(|i|if i%2==0 {1e-5*((i+1) as f64).sin()}else{0.03}).collect();
    // This all-mode shell includes very stiff rotational coordinates. At
    // h=1e-9 cancellation in large gradient entries obscures small off-diagonal
    // terms. Shell H is quartic, so this four-point derivative is EXACT on its
    // cubic gradient in real arithmetic. Increase the probe spacing, not the
    // agreement tolerance; no contact/felt branch is crossed in this fixture.
    check_tangent_at(&s,&x,1e-4,true);
    // Tangent really changes with the nonlinear strain, not just frequency.
    let d=vec![1.0;2*n];let mut a=d.clone();let mut zero=d.clone();
    assert!(s.hessian_vector(&x,&d,&mut a));assert!(s.hessian_vector(&vec![0.0;2*n],&d,&mut zero));
    assert!(a.iter().zip(zero).any(|(a,b)|(a-b).abs()>1e-3));
}
#[test]
fn relaxed_membrane_tangent_retains_all_mixed_modes_and_installed_tension() {
    let b=body(membrane());let n=b.potential.count();
    let s=ImpactSystem::new(vec![b],vec![],vec![],vec![],config()).unwrap();
    let x:Vec<_>=(0..2*n).map(|i|if i%2==0 {2e-4*((i+1) as f64).sin()}else{-0.02}).collect();
    check_tangent(&s,&x);
}
#[test]
fn contact_felt_creep_and_reciprocal_volume_share_the_exact_frozen_history_tangent() {
    let contact=Obstacle::new(vec![2.0,-3.0,-1.0,2.0],2,2,vec![0.0,0.001],vec![0.4,0.6],
        2e7,1.5,"tangent regression, not calibrated contact".into()).unwrap();
    let s=ImpactSystem::new(vec![body(BodyPotential::Linear(vec![800.0,900.0]))],vec![contact],
        vec![pad(vec![1.0,-0.7],0.0003),pad(vec![0.5,-1.0],0.0018)],
        vec![VolumeSpring {bulk_modulus_pa:1.4e5,volume_m3:0.01,areas:vec![0.1,-0.08]}],config()).unwrap();
    let history=s.histories.borrow().clone();
    // Positive contact; first pad unloading and second loading, with both Kelvin tails.
    let x=[1e-4,0.2,-1e-4,-0.1,1e-5,-2e-5];check_tangent(&s,&x);
    // Detach both contacts and the first felt pad; no negative-penetration stiffness.
    check_tangent(&s,&[-0.0005,0.2,0.0,-0.1,0.0,0.0]);
    assert_eq!(*s.histories.borrow(),history);assert_eq!(s.samples(),0);
}
#[test]
fn distributed_contact_action_does_not_assume_unit_mass_or_discard_inner_tail() {
    let ob=Obstacle::new(vec![2.0,-3.0],1,2,vec![0.0],vec![0.5],100.0,1.0,"linear penalty regression".into()).unwrap();
    let diagonal=[7.0,0.5,11.0,0.25,3.0];let mut q=vec![0.0;25];
    for i in 0..5 {q[5*i+i]=diagonal[i];}
    let s=ContactStorage::new(Box::new(fs_phs::QuadraticStorage::new(q,5).unwrap()),2,vec![ob]).unwrap();
    let action=|_:&[f64],d:&[f64],o:&mut[f64]| {for i in 0..5{o[i]=diagonal[i]*d[i];}true};
    let d=[1.0,2.0,-0.5,4.0,6.0];let mut out=[0.0;5];
    for active in [false,true] {
        let x=[if active {0.1}else{-0.1},0.0,0.0,0.0,0.0];
        assert!(s.hessian_vector_with(&x,&d,&mut out,action));
        let reaction=if active {50.0*(2.0*d[0]-3.0*d[2])}else{0.0};
        assert_eq!(out,[7.0*d[0]+2.0*reaction,1.0,11.0*d[2]-3.0*reaction,1.0,18.0]);
    }
    assert!(!s.hessian_vector_with(&[0.0;5],&d,&mut out,|_,_,_|false));
}
fn strike()->ImpactSystem {
    let (stick,w)=ImpactBody::free_mass(0.02,-0.0002,0.8).unwrap();
    let c=Obstacle::new(vec![w,-2.0],1,2,vec![0.0],vec![1.0],2e7,1.5,"Hertz regression".into()).unwrap();
    ImpactSystem::new(vec![stick,body(BodyPotential::Linear(vec![core::f64::consts::TAU*500.0]))],vec![c],
        vec![pad(vec![0.0,2.0],0.0006)],vec![],config()).unwrap()
}
#[test]
fn analytic_hertz_rebound_keeps_felt_history_and_complete_work_balance() {
    let mut reference=strike().prepare().unwrap();let mut analytic=strike().prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();let initial=analytic.stored_energy_j();let mut net=0.0;let mut rebound=false;
    for tick in 0..1200 {
        let force=[if tick<50 {0.1}else{0.0},0.0];
        reference.step(&force,&gate).unwrap();let f=analytic.step(&force,&gate).unwrap();
        net+=f.supplied_work_j-f.dissipated_energy_j;rebound|=analytic.state()[1]<0.0;
        assert!(f.balance_residual_j.abs()<1e-8);
        for (a,b) in analytic.state().iter().zip(reference.state()) {
            assert!((a-b).abs()<2e-6*a.abs().max(b.abs()).max(1e-4),"{tick}: {a:e} != {b:e}");
        }
    }
    assert!(rebound && analytic.state()[4]!=0.0);
    assert!((analytic.stored_energy_j()-initial-net).abs()<1e-7);
    assert!((analytic.felt_observation(0).unwrap().0-reference.felt_observation(0).unwrap().0).abs()<1e-7);
}
#[test]
fn analytic_refusals_and_switching_preserve_motion_history_and_exact_retry() {
    let mut trial=strike().prepare_analytic().unwrap();let mut clean=strike().prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..140 {trial.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();}
    let x=trial.state().to_vec();let h=trial.felt_history(0);
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(trial.step(&[0.0;2],&cancelled),Err(ImpactError::Cancelled)));
    trial.set_iteration_limit(0).unwrap();assert!(trial.step(&[0.0;2],&gate).is_err());
    assert_eq!(trial.state(),x);assert_eq!(trial.felt_history(0),h);assert_eq!(trial.samples(),140);
    trial.set_iteration_limit(50).unwrap();trial.set_analytic_newton(false);trial.set_analytic_newton(true);
    trial.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();
    assert!(trial.state().iter().zip(clean.state()).all(|(a,b)|a.to_bits()==b.to_bits()));
    assert_eq!(trial.felt_history(0),clean.felt_history(0));
    let x=trial.state().to_vec();let h=trial.felt_history(0);
    let trial=trial.into_reference().prepare_analytic().unwrap();
    assert_eq!(trial.state(),x);assert_eq!(trial.felt_history(0),h);assert_eq!(trial.samples(),141);
}
