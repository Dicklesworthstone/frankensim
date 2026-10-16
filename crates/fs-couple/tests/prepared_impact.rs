//! Exercise the real owners. Numerical fixtures are not calibrated instruments.
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, ImpactError, ImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::audio::ImpactSource;
use fs_couple::render::plate::impact::felt::{FeltPad, KelvinBranch};
use fs_couple::render::plate::impact::membrane::MembranePotential;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::fiber::WoolFelt;

fn config() -> ImpactConfig {
    ImpactConfig { dt_s: 2e-6, max_steps: 2000, maximum_energy_j: 20.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-7,
        maximum_generalized_force: 1e4 }
}
fn body(omega: f64, q: f64, v: f64) -> ImpactBody {
    ImpactBody { potential: BodyPotential::Linear(vec![omega]),
        initial: vec![ModalAcousticState { displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v }], damping_per_s: vec![1.0] }
}
fn pad(weights: Vec<f64>, compression: f64) -> FeltPad {
    FeltPad { area_m2: 0.001, thickness_m: 0.006, precompression_m: compression,
        weights, law: WoolFelt::new(30000.0, 0.2, 2.2, 3.0, 0.15, 0.7).unwrap(),
        prior_maximum_strain: 0.1,
        creep: vec![KelvinBranch { stiffness_n_m: 1500.0, viscosity_n_s_m: 8.0 }] }
}
fn strike(felt: bool, connected: bool) -> ImpactSystem {
    let (stick, weight) = ImpactBody::free_mass(0.02, -0.0002, 0.8).unwrap();
    let contact = Obstacle::new(vec![weight, -2.0], 1, 2, vec![0.0], vec![1.0],
        if connected { 2e7 } else { 0.0 }, 1.5, "synthetic Hertz regression".into()).unwrap();
    ImpactSystem::new(vec![stick, body(2.0*std::f64::consts::PI*500.0, 0.0, 0.0)],
        vec![contact], if felt { vec![pad(vec![0.0, 2.0], 0.0006)] } else { vec![] },
        vec![], config()).unwrap()
}
fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        assert!((a-b).abs() < 2e-7*a.abs().max(b.abs()).max(1e-4), "{a:e} != {b:e}");
    }
}

#[test]
fn prepared_hertz_strike_rebounds_and_retains_reference_felt_and_work() {
    let gate = CancelGate::new_clock_free();
    for felt in [false, true] {
        let mut reference = strike(felt, true);
        let mut prepared = strike(felt, true).prepare().unwrap();
        let initial = prepared.stored_energy_j();
        let mut losses = 0.0; let mut receiver = 0.0_f64; let mut rebound = false;
        for _ in 0..1200 {
            let a = reference.step(&[0.0; 2], &gate).unwrap();
            let b = prepared.step(&[0.0; 2], &gate).unwrap();
            assert_eq!(a.sample, b.sample); assert_eq!(a.time_s.to_bits(), b.time_s.to_bits());
            close(reference.state(), prepared.state());
            assert!((a.stored_energy_j-b.stored_energy_j).abs() < 1e-8);
            assert!(b.balance_residual_j.abs() < 1e-8);
            losses += b.dissipated_energy_j;
            receiver = receiver.max(prepared.state()[2].abs()); rebound |= prepared.state()[1] < 0.0;
        }
        assert!(receiver > 1e-6 && rebound);
        assert!((prepared.stored_energy_j()-initial+losses).abs() < 1e-7);
        if felt {
            assert!(prepared.state()[4].abs() > 0.0);
            assert!((prepared.felt_observation(0).unwrap().0-reference.felt_observation(0).unwrap().0).abs() < 1e-7);
        }
    }
}

#[test]
fn conversion_preserves_an_ongoing_strike_history_and_accepted_clock() {
    let mut reference = strike(true, true); let gate = CancelGate::new_clock_free();
    for _ in 0..180 { reference.step(&[0.0; 2], &gate).unwrap(); }
    let state = reference.state().to_vec(); let history = reference.felt_history(0);
    let mut prepared = reference.prepare().unwrap();
    assert_eq!(prepared.state(), state); assert_eq!(prepared.felt_history(0), history);
    assert_eq!(prepared.samples(), 180);
    prepared.step(&[0.0; 2], &gate).unwrap();
    let state = prepared.state().to_vec(); let history = prepared.felt_history(0);
    let returned = prepared.into_reference();
    assert_eq!(returned.state(), state); assert_eq!(returned.felt_history(0), history);
    assert_eq!(returned.samples(), 181);
}

#[test]
fn cancellation_bad_force_and_newton_refusal_do_not_poison_exact_retry() {
    let mut prepared = strike(true, true).prepare().unwrap();
    let mut clean = strike(true, true).prepare().unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..140 { prepared.step(&[0.0; 2], &gate).unwrap(); clean.step(&[0.0; 2], &gate).unwrap(); }
    let state = prepared.state().to_vec(); let history = prepared.felt_history(0);
    let cancel = CancelGate::new_clock_free(); cancel.request();
    assert!(matches!(prepared.step(&[0.0; 2], &cancel), Err(ImpactError::Cancelled)));
    for force in [&[0.0][..], &[f64::NAN,0.0][..], &[1e5,0.0][..]] {
        assert!(prepared.step(force, &gate).is_err());
    }
    prepared.set_iteration_limit(0).unwrap();
    assert!(matches!(prepared.step(&[0.0; 2], &gate), Err(ImpactError::PreparedSolve(_))));
    assert_eq!(prepared.state(), state); assert_eq!(prepared.felt_history(0), history);
    assert_eq!(prepared.samples(), 140);
    prepared.set_iteration_limit(50).unwrap();
    for count in [3, 9, 1, 27] { for _ in 0..count {
        prepared.step(&[0.0; 2], &gate).unwrap(); clean.step(&[0.0; 2], &gate).unwrap();
    }}
    assert!(prepared.state().iter().zip(clean.state()).all(|(a,b)| a.to_bits()==b.to_bits()));
    assert_eq!(prepared.felt_history(0), clean.felt_history(0));
}

#[test]
fn energy_and_felt_domain_rejections_do_not_publish_solved_candidates() {
    let mut cfg = config(); cfg.maximum_energy_j = 1e-12;
    let mut tiny = ImpactSystem::new(vec![body(200.0,0.0,0.0)],vec![],vec![],vec![],cfg)
        .unwrap().prepare().unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..2 {
        assert!(matches!(tiny.step(&[100.0],&gate),Err(ImpactError::Invalid("impact candidate exceeds finite energy limits"))));
        assert_eq!(tiny.state(), &[0.0,0.0]); assert_eq!(tiny.samples(),0);
    }
    let mut p = pad(vec![1.0],0.69999*0.006); p.creep.clear();
    let (mass,_) = ImpactBody::free_mass(1.0,0.0,1.0).unwrap();
    let mut felt = ImpactSystem::new(vec![mass],vec![],vec![p],vec![],config()).unwrap().prepare().unwrap();
    let state = felt.state().to_vec(); let history = felt.felt_history(0);
    for _ in 0..2 {
        assert!(matches!(felt.step(&[0.0],&gate),Err(ImpactError::Invalid("felt trial exceeds densification validity"))));
        assert_eq!(felt.state(),state); assert_eq!(felt.felt_history(0),history); assert_eq!(felt.samples(),0);
    }
}

#[test]
fn prepared_pressure_source_has_exact_budget_and_disconnected_silence() {
    let mut source = strike(false,false).prepare().unwrap();
    let gate = CancelGate::new_clock_free();
    assert_eq!(ImpactSource::mode_count(&source),2);
    assert_eq!(ImpactSource::sample_period_s(&source).to_bits(),config().dt_s.to_bits());
    for _ in 0..config().max_steps {
        ImpactSource::advance(&mut source,&[0.0;2],&gate).unwrap();
        assert_eq!(source.state()[2],0.0); assert_eq!(source.state()[3],0.0);
    }
    assert_eq!(ImpactSource::remaining_steps(&source),0);
    let state = source.state().to_vec();
    assert!(matches!(source.step(&[0.0;2],&gate),Err(ImpactError::Budget)));
    assert_eq!(source.state(),state);
}

#[test]
fn sealed_volume_drives_the_other_head_in_the_same_nonlinear_step() {
    let make = || ImpactSystem::new(vec![body(800.0,1e-5,0.0),body(900.0,0.0,0.0)],vec![],vec![],
        vec![VolumeSpring { bulk_modulus_pa:1.4e5,volume_m3:0.01,areas:vec![0.1,-0.1] }],config()).unwrap();
    let mut reference = make(); let mut prepared = make().prepare().unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..100 {
        reference.step(&[0.0;2],&gate).unwrap(); prepared.step(&[0.0;2],&gate).unwrap();
        close(reference.state(),prepared.state());
    }
    assert!(prepared.state()[2].abs()>1e-9);
}

// Actual geometric film assembly and eigensolve, not a frequency preset.
fn membrane(limit: f64) -> MembranePotential {
    use fs_plate::shell::head::{TensionedDisk,TensionedDiskSpec};
    use fs_plate::shell::head::nonlinear::MembraneReductionBudget;
    use fs_plate::shell::profile::ProfileBudget;
    let d = TensionedDisk::new(TensionedDiskSpec {radius_m:0.17,thickness_m:0.000254,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:3000.0,
        radial_intervals:2,azimuths:8},ProfileBudget {max_nodes:100,max_triangles:100,max_feature_evaluations:0}).unwrap();
    let n=d.model.free; let mut k=vec![0.0;n*n]; let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=d.model.k.get(i,j);m[i*n+j]=d.model.m.get(i,j);}}
    let mut modes=fs_modal::eigh_gen_dense(&k,&m,n).unwrap(); modes.truncate(1);
    let rim:Vec<_>=(0..d.mesh.nodes.len()).filter(|&i|d.model.dof_map[3*i].is_none()).collect();
    MembranePotential::from_pencil(&d.mesh,&d.section,&d.model,&modes,&rim,
        MembraneReductionBudget {max_modes:4,max_nodes:100,max_facet_pairs:1000,
            max_solve_entries:100000,relative_tolerance:1e-6},limit).unwrap()
}
#[test]
fn geometry_owned_stretching_and_its_slope_gate_survive_preparation() {
    let potential=membrane(0.2);
    let q=0.08/potential.reduction().maximum_slope(&[1.0]);
    let make=|law:MembranePotential,q,v| {
        let mut b=body(law.reduction().omegas()[0],q,v); b.potential=BodyPotential::Membrane(law);
        ImpactSystem::new(vec![b],vec![],vec![],vec![],config()).unwrap()
    };
    let mut reference=make(potential.clone(),q,0.0);
    let mut prepared=make(potential,q,0.0).prepare().unwrap(); let gate=CancelGate::new_clock_free();
    assert!(prepared.membrane_observation(0).unwrap().stretching_energy_j>0.0);
    for _ in 0..128 {
        reference.step(&[0.0],&gate).unwrap(); prepared.step(&[0.0],&gate).unwrap();
        close(reference.state(),prepared.state());
    }
    let mut strict=make(membrane(1e-7),0.0,1.0).prepare().unwrap();
    let state=strict.state().to_vec();
    assert!(matches!(strict.step(&[0.0],&gate),Err(ImpactError::Invalid("membrane state exceeds its declared slope or finite stretching-energy validity"))));
    assert_eq!(strict.state(),state); assert_eq!(strict.samples(),0);
}
