use super::*;
use super::super::{BodyPotential, ImpactConfig, ImpactSystem, ImpactSubstepConfig};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;
use fs_phs::Storage;

fn jaw(side: PadSide) -> CompliantJaw {
    CompliantJaw { side, mass_kg: 0.025, drag_n_s_m: 0.0,
        initial_gap_m: 0.0, initial_velocity_m_s: 0.0, thickness_m: 0.004,
        law: WoolFelt::new(30000.0, 0.2, 2.2, 3.0, 0.15, 0.7).unwrap(),
        prior_maximum_strain: 0.0, creep: vec![] }
}
fn config(steps: u64, dt_s: f64) -> ImpactConfig {
    ImpactConfig { dt_s, max_steps: steps, maximum_energy_j: 20.0,
        energy_absolute_tolerance_j: 1e-9, energy_relative_tolerance: 1e-6,
        maximum_generalized_force: 1e6 }
}
fn site(weights: Vec<f64>, area_m2: f64) -> PadSite { PadSite { weights, area_m2 } }

#[test]
fn footprint_keeps_local_contact_rows_shared_inertias_and_physical_creep_totals() {
    let mut top = jaw(PadSide::Positive);
    top.initial_gap_m = 0.0002; top.drag_n_s_m = 0.3;
    top.creep.push(KelvinBranch { stiffness_n_m: 4000.0, viscosity_n_s_m: 12.0 });
    let sites = [site(vec![0.0, 2.0, -3.0], 0.0001), site(vec![0.0, -4.0, 5.0], 0.0003)];
    let a = MovingPads::new(3, &sites, &[top, jaw(PadSide::Negative)]).unwrap();
    assert_eq!(a.bodies.len(), 2); assert_eq!(a.pads.len(), 4);
    assert_eq!(a.ports[0].coordinate, 3); assert_eq!(a.ports[1].coordinate, 4);
    assert_eq!(&a.pads[0].weights[..3], &sites[0].weights);
    assert_eq!(&a.pads[2].weights[..3], &[0.0, -2.0, 3.0]);
    assert_eq!(a.pads[0].weights[4], 0.0); assert_eq!(a.pads[2].weights[3], 0.0);
    assert_eq!(a.pads[0].weights[3], 1.0/0.025_f64.sqrt());
    assert_eq!(a.pads[0].precompression_m, -0.0002);
    assert!((a.bodies[0].damping_per_s[0]-12.0).abs() < 1e-14);
    assert_eq!(a.pads[0].creep[0].stiffness_n_m, 1000.0);
    assert_eq!(a.pads[1].creep[0].stiffness_n_m, 3000.0);
    assert_eq!(a.pads[0].creep[0].viscosity_n_s_m+a.pads[1].creep[0].viscosity_n_s_m, 12.0);
}

#[test]
fn opposed_contact_is_reciprocal_compression_only_and_not_an_average_shape() {
    let a = MovingPads::new(2, &[site(vec![1.0,1.0],0.0002),
        site(vec![1.0,-1.0],0.0002)], &[jaw(PadSide::Positive),jaw(PadSide::Negative)]).unwrap();
    let surface = ImpactBody { potential: BodyPotential::Linear(vec![0.0;2]),
        initial: vec![ModalAcousticState::default();2], damping_per_s: vec![0.0;2] };
    let mut bodies=vec![surface]; bodies.extend(a.bodies);
    let s=ImpactSystem::new(bodies,vec![],a.pads,vec![],config(1,1e-6)).unwrap();
    let mut x=s.state().to_vec(); x[2]=0.0001;
    let mut g=vec![0.0;x.len()]; s.mechanical.gradient(&x,&mut g);
    // The mean surface row [1,0] would miss this mode completely. Each side
    // instead contacts one physical site and both jaws receive a reaction.
    assert!(s.mechanical.hamiltonian(&x)>0.0); assert!(g[2]>0.0);
    assert!(g[4]>0.0 && g[6]>0.0); assert!(g[0].abs()<1e-12);
    // Uniform axial translation of surface AND jaws leaves compression fixed.
    let force_balance=g[0]-g[4]/a.ports[0].inverse_sqrt_mass
        +g[6]/a.ports[1].inverse_sqrt_mass;
    assert!(force_balance.abs()<1e-12);
    // Pull both jaws back far enough: no residual attractive contact force.
    x[4]=-0.001/a.ports[0].inverse_sqrt_mass;
    x[6]=-0.001/a.ports[1].inverse_sqrt_mass;
    s.mechanical.gradient(&x,&mut g);
    assert!(g.iter().all(|f|*f==0.0)); assert_eq!(s.mechanical.hamiltonian(&x),0.0);
}

#[test]
fn uniform_area_partition_preserves_static_force_energy_and_jaw_mass() {
    let sites_one=[site(vec![2.0],0.0004)];
    let sites_four=vec![site(vec![2.0],0.0001);4];
    let make=|sites:&[PadSite]| {
        let a=MovingPads::new(1,sites,&[jaw(PadSide::Negative)]).unwrap();
        let (surface,_)=ImpactBody::free_mass(0.25,0.0,0.0).unwrap();
        let mut bodies=vec![surface];bodies.extend(a.bodies);
        ImpactSystem::new(bodies,vec![],a.pads,vec![],config(1,1e-6)).unwrap()
    };
    let one=make(&sites_one);let four=make(&sites_four);
    let x=[-0.0001,0.04,0.00002,0.01];let mut a=[0.0;4];let mut b=[0.0;4];
    one.mechanical.gradient(&x,&mut a);four.mechanical.gradient(&x,&mut b);
    for (a,b) in a.into_iter().zip(b) {assert!((a-b).abs()<1e-12);}
    assert!((one.mechanical.hamiltonian(&x)-four.mechanical.hamiltonian(&x)).abs()<1e-14);
}

#[test]
fn squeeze_then_retract_uses_real_contact_work_loss_and_atomic_retry() {
    let make=|| {
        let mut j=jaw(PadSide::Negative);j.initial_gap_m=0.00002;j.drag_n_s_m=0.1;
        j.creep.push(KelvinBranch {stiffness_n_m:4000.0,viscosity_n_s_m:12.0});
        let a=MovingPads::new(1,&[site(vec![1.0],0.0004)],&[j]).unwrap();
        let port=a.ports[0];let (surface,_)=ImpactBody::free_mass(1.0,0.0,0.0).unwrap();
        let mut bodies=vec![surface];bodies.extend(a.bodies);
        let s=ImpactSystem::new(bodies,vec![],a.pads,vec![],config(2000,2e-6)).unwrap()
            .prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig {
                max_depth:8, max_attempts:511 }).unwrap();
        (s,port)
    };
    let (mut s,port)=make();let (mut clean,_)=make();let gate=CancelGate::new_clock_free();
    let mut work=0.0;let mut loss=0.0;let mut maximum_history=0.0_f64;
    for tick in 0..2000 {
        // Positive force closes the real gap. Negative force later retracts;
        // zero force alone is never interpreted as a teleport/release command.
        let force=if tick<800 {1.0}else{-2.0};let external=[0.0,force*port.inverse_sqrt_mass];
        if tick==900 {
            let before=s.state().to_vec();
            assert!(s.step(&[0.0,1e7],&gate).is_err());assert_eq!(s.state(),before);
        }
        let f=s.step(&external,&gate).unwrap();clean.step(&external,&gate).unwrap();
        assert_eq!(s.state(),clean.state());work+=f.supplied_work_j;loss+=f.dissipated_energy_j;
        // Inspect through the same preparation owner, not reconstructed histories.
        maximum_history=maximum_history.max(s.felt_history(0).unwrap().eps_max);
        assert!((f.stored_energy_j+loss-work).abs()<1e-7);
    }
    assert!(maximum_history>0.0);assert!(loss>0.0);
    assert!(s.state()[2]*port.inverse_sqrt_mass<0.0); // actual retraction past the initial jaw plane
    assert!(s.state()[0]>0.0); // surface received the inward contact impulse
}

#[test]
fn malformed_geometry_material_and_hidden_budget_expansion_refuse() {
    let good=[site(vec![1.0],0.0004)];let j=jaw(PadSide::Positive);
    assert!(MovingPads::new(1,&[],&[j.clone()]).is_err());
    assert!(MovingPads::new(1,&good,&[]).is_err());
    assert!(MovingPads::new(1,&good,&[j.clone(),j.clone()]).is_err());
    assert!(MovingPads::new(1,&vec![good[0].clone();5],&[j.clone()]).is_err());
    // Use the mechanical owner's current ceiling: a concurrent snare change
    // enlarged it. Adding a jaw must still consume one of those coordinates.
    assert!(MovingPads::new(MAX_IMPACT_MODES,
        &[site(vec![1.0;MAX_IMPACT_MODES],0.0004)],&[j.clone()]).is_err());
    let at_limit=MovingPads::new(MAX_IMPACT_MODES-1,
        &[site(vec![1.0;MAX_IMPACT_MODES-1],0.0004)],&[j.clone()]).unwrap();
    assert_eq!(at_limit.ports[0].coordinate,MAX_IMPACT_MODES-1);
    assert_eq!(at_limit.pads[0].weights.len(),MAX_IMPACT_MODES);
    for row in [vec![0.0],vec![f64::NAN],vec![1.0,2.0]] {
        assert!(MovingPads::new(1,&[site(row,0.0004)],&[j.clone()]).is_err());
    }
    for mass in [0.0,-1.0,f64::INFINITY] {
        assert!(MovingPads::new(1,&good,&[CompliantJaw {mass_kg:mass,..j.clone()}]).is_err());
    }
    for gap in [-0.001,f64::NAN] {
        assert!(MovingPads::new(1,&good,&[CompliantJaw {initial_gap_m:gap,..j.clone()}]).is_err());
    }
    assert!(MovingPads::new(1,&good,&[CompliantJaw {thickness_m:0.0,..j.clone()}]).is_err());
    assert!(MovingPads::new(1,&good,&[CompliantJaw {drag_n_s_m:-1.0,..j}]).is_err());
}
