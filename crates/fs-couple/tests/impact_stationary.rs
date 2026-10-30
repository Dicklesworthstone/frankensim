//! Held physical loads resolved through the original impact/contact/film owners.
use fs_couple::render::plate::impact::{BodyPotential,ImpactBody,ImpactConfig,ImpactSystem};
use fs_couple::render::plate::impact::relaxation::InitialMemory;
use fs_couple::render::plate::impact::squeeze::{FilmCell,FilmChannel,FilmLimits,GapPort,ResistiveFilm};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::stationary_equilibrium::EquilibriumSolveOptions;
use fs_phs::harmonic::equilibrium::{EquilibriumPoint,LinearizationBudget};
use fs_phs::harmonic::equilibrium::stationary::StationaryCoordinates;
use fs_phs::harmonic::{HarmonicBudget,PoleOptions};
use fs_phs::RelaxationBranch;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_solver::NewtonKrylovConfig;
use std::cell::Cell;

fn config()->ImpactConfig {ImpactConfig {dt_s:1e-5,max_steps:1000,maximum_energy_j:100.,
    energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:1000.}}
fn budget()->LinearizationBudget {LinearizationBudget {harmonic:HarmonicBudget {
    max_dimension:32,max_dense_bytes:2_000_000,max_dense_work:10_000_000,
    relative_residual_tolerance:1e-10},max_hessian_actions:32}}
fn options()->EquilibriumSolveOptions {EquilibriumSolveOptions {
    newton:NewtonKrylovConfig {absolute_tolerance:1e-12,relative_tolerance:0.,
        forcing_minimum:1e-12,forcing_maximum:1e-6,linear_restart:8,max_linear_cycles:8,
        ..Default::default()},relative_balance_tolerance:1e-11,max_iterations:40,
    max_residual_calls:256,max_jacobian_calls:2048,max_workspace_bytes:2_000_000}}
fn body(w:Vec<f64>,q:&[f64],d:Vec<f64>)->ImpactBody {
    ImpactBody {potential:BodyPotential::Linear(w),damping_per_s:d,
        initial:q.iter().map(|&q|ModalAcousticState {
            displacement_m_sqrt_kg:q,velocity_m_sqrt_kg_per_s:0.,
        }).collect()}
}
fn obstacle(row:Vec<f64>,alpha:f64,chi:f64)->Obstacle {
    let n=row.len();Obstacle::new(row,1,n,vec![0.],vec![1.],1000.,alpha,
        "manufactured stationary fixture; not material qualification".into()).unwrap()
        .with_internal_loss(chi).unwrap()
}
fn one()->ImpactSystem {
    ImpactSystem::new(vec![body(vec![10.],&[0.003],vec![0.2])],
        vec![obstacle(vec![1.],1.5,0.3)],vec![],vec![],config()).unwrap()
}
fn near(a:C64,b:C64) {assert!((a-b).abs()<2e-9*(1.+b.abs()),"{a:?} != {b:?}");}

#[test]
fn g1_force_loading_determines_contact_frequency_and_decay_without_installing_a_state() {
    let system=one();let before=system.state().to_vec();
    for target in [0.01_f64,0.02,0.04] {
        let force=[100.*target+1000.*target.powf(1.5)];
        let coordinates=StationaryCoordinates {input:&force,state_scales:&[0.02,1.],balance_scales:&[1.,10.]};
        let result=system.solve_equilibrium(&[0.,0.],coordinates,budget(),32,options(),||false).unwrap();
        let image=result.linearization;let x=image.equilibrium_state();
        assert!((x[0]-target).abs()<1e-10 && x[1].abs()<1e-10);
        let k=100.+1500.*target.sqrt();let d=0.2+300.*target.powf(1.5);
        let exact=C64::new(-0.5*d,-(k-0.25*d*d).sqrt());
        let pole=image.model().resonance_pole(exact,&[C64::ONE,exact],&PoleOptions {
            spectral_scale:30.,step_tolerance:1e-12,residual_tolerance:1e-12,
            max_iterations:24,max_dimension:32,max_dense_bytes:2_000_000,max_work:10_000_000,
        },||false).unwrap();
        near(pole.laplace_pole_per_s,exact);
        assert!((image.model().nonlinear_resisting_tangent().unwrap()[3]-(d-0.2)).abs()<1e-8);
        let direct=system.linearize_equilibrium(EquilibriumPoint {state:x,input:&force,
            balance_scales:&[1.,10.],relative_tolerance:1e-11},budget(),32,||false).unwrap();
        assert_eq!(direct.tangent_storage_matrix(),image.tangent_storage_matrix());
    }
    assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
}

#[test]
fn g1_material_history_is_solved_in_original_coordinates_and_remains_in_response() {
    let system=one().with_relaxation_branches(vec![RelaxationBranch {
        projection:vec![1.,0.],stiffness:30.,relaxation_time_s:0.2,
    }],InitialMemory::Unrelaxed,1).unwrap();
    let before=system.state().to_vec();let energy=system.stored_energy_j();
    let target=0.02_f64;let force=[100.*target+1000.*target.powf(1.5)];
    let result=system.solve_equilibrium(&before,StationaryCoordinates {input:&force,
        state_scales:&[0.02,1.,0.1],balance_scales:&[1.,10.,1.]},budget(),32,options(),||false).unwrap();
    let image=result.linearization;let x=image.equilibrium_state();
    assert!((x[0]-target).abs()<1e-10);
    assert!((x[2]-30.0_f64.sqrt()*target).abs()<1e-10);
    for w in [5.,10.,20.,80.] {
        let s=C64::new(0.,-w);let divisor=C64::ONE+s.scale(0.2);
        let k=100.+1500.*target.sqrt();let d=0.2+300.*target.powf(1.5);
        let q=(s*s+s.scale(d)+C64::from_re(k)+s.scale(6.)/divisor).recip();
        let r=image.model().factor(w,||false).unwrap().response(&[C64::ONE],||false).unwrap();
        near(r.state[0],q);near(r.state[2],q.scale(30.0_f64.sqrt())/divisor);
    }
    assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
    assert_eq!(system.stored_energy_j().to_bits(),energy.to_bits());
}

#[test]
fn g3_signed_contact_and_gap_dependent_film_use_the_resolved_common_configuration() {
    let mut system=ImpactSystem::new(vec![body(vec![10.,15.],&[0.,0.],vec![0.2,0.4])],
        vec![obstacle(vec![1.,-1.],1.5,100.)],vec![],vec![],config()).unwrap();
    let gap=GapPort {reference_m:0.001,closure:vec![1.,-1.]};
    let film=ResistiveFilm::new(vec![FilmCell {area_m2:0.002,gap:gap.clone()}],
        vec![FilmChannel {from:0,to:None,width_m:0.05,length_m:0.01,gap}],2,1.8e-5,
        FilmLimits {minimum_cell_gap_m:1e-6,maximum_gap_m:0.01,maximum_pressure_pa:1000.}).unwrap();
    system=system.with_squeeze_film(film).unwrap();
    let target=[0.0001,-0.0001];let penetration=0.0002_f64;
    let contact=1000.*penetration.powf(1.5);
    let force=[100.*target[0]+contact,225.*target[1]-contact];
    let result=system.solve_equilibrium(system.state(),StationaryCoordinates {input:&force,
        state_scales:&[0.001,0.01,0.001,0.01],balance_scales:&[0.01,1.,0.01,1.]},
        budget(),32,options(),||false).unwrap();
    let image=result.linearization;let x=image.equilibrium_state();
    for i in 0..2 {assert!((x[2*i]-target[i]).abs()<1e-11);}
    let df=0.002_f64.powi(2)*12.*1.8e-5*0.01/(0.05*0.0008_f64.powi(3));
    let expected=100.*contact+df;
    let t=image.model().nonlinear_resisting_tangent().unwrap();
    for (i,j,sign) in [(1,1,1.),(1,3,-1.),(3,1,-1.),(3,3,1.)] {
        assert!((t[i*4+j]-sign*expected).abs()<1e-8);
    }
    let impossible=[0.002,0.,0.,0.];
    assert!(system.solve_equilibrium(&impossible,StationaryCoordinates {input:&force,
        state_scales:&[1.;4],balance_scales:&[1.;4]},budget(),32,options(),||false).is_err());
    assert_eq!(system.state(),[0.;4]);assert_eq!(system.samples(),0);
}

#[test]
fn g0_g4_failed_and_successful_bias_solves_do_not_change_the_accepted_time_trajectory() {
    let mut system=one();let mut reference=one();let before=system.state().to_vec();
    let force=[100.*0.02+1000.*0.02_f64.powf(1.5)];
    let coordinates=StationaryCoordinates {input:&force,state_scales:&[0.02,1.],balance_scales:&[1.,10.]};
    let count=Cell::new(0);
    let baseline=system.solve_equilibrium(&before,coordinates,budget(),32,options(),||{
        count.set(count.get()+1);false}).unwrap();
    let total=count.get();
    for cut in [1,total/2,total] {
        let count=Cell::new(0);
        assert!(system.solve_equilibrium(&before,coordinates,budget(),32,options(),||{
            count.set(count.get()+1);count.get()>=cut}).is_err());
    }
    let mut short=options();short.max_residual_calls=1;
    assert!(system.solve_equilibrium(&before,coordinates,budget(),32,short,||false).is_err());
    assert!(system.solve_equilibrium(&before,coordinates,budget(),1,options(),||false).is_err());
    let excess=[1001.];
    assert!(system.solve_equilibrium(&before,StationaryCoordinates{input:&excess,..coordinates},
        budget(),32,options(),||false).is_err());
    let replay=system.solve_equilibrium(&before,coordinates,budget(),32,options(),||false).unwrap();
    assert_eq!(baseline.linearization.equilibrium_state(),replay.linearization.equilibrium_state());
    assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
    // Even an exact zero-balance point is not an admissible linear resonance
    // when an active-branch tangent would conceal a derivative corner.
    let kink=ImpactSystem::new(vec![body(vec![10.],&[0.],vec![0.2])],
        vec![obstacle(vec![1.],1.,0.)],vec![],vec![],config()).unwrap();
    assert!(kink.solve_equilibrium(kink.state(),StationaryCoordinates {input:&[0.],
        state_scales:&[1.;2],balance_scales:&[1.;2]},budget(),32,options(),||false).is_err());
    let gate=CancelGate::new_clock_free();
    for _ in 0..8 {system.step(&force,&gate).unwrap();reference.step(&force,&gate).unwrap();}
    assert_eq!(system.state(),reference.state());assert_eq!(system.samples(),reference.samples());
}
