//! Physical distributed-contact checks; numerical fixtures are not measured snares.
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody};
use fs_couple::render::plate::impact::linear::{LinearImpactConfig, LinearImpactSystem};
use fs_couple::render::plate::impact::linear::wire::{HelicalWire, WireSpan, LineContact, film_shapes};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig,
    contact::{ModalContactConfig, multiple::MultiContactConfig}};
use fs_dcontact::ContactStorage;
use fs_exec::CancelGate;
use fs_phs::Storage;
use fs_plate::{ModePair, shell::{head::{TensionedDisk, TensionedDiskSpec}, profile::ProfileBudget}};

fn wire() -> WireSpan {
    WireSpan { endpoints_m: [[-0.15, 0.0], [0.15, 0.0]], linear_density_kg_m: 0.002,
        tension_n: 0.7, bending_stiffness_n_m2: 1e-6, damping_per_s: vec![2.0; 3] }
}
fn line(n: usize, k: f64) -> LineContact {
    LineContact::uniform(0.3, n, 2e-5, k, 1.5, 0.05,
        "synthetic distributed compliant wire/head contact".into()).unwrap()
}
fn config(steps: u64) -> LinearImpactConfig {
    LinearImpactConfig { sample_rate_hz: 192000, max_steps: steps, maximum_generalized_force: 1e5,
        component: ModalAcousticTimeBudget::audible_reference(),
        coupling: ModalCouplingConfig { max_modes: 128, max_connections: 8, max_setup_terms: 1000000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e5, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8 },
        contact: ModalContactConfig { max_iterations: 100, maximum_force_n: 1e4,
            maximum_penetration_m: 0.01, force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-9 },
        multiple: MultiContactConfig { max_contacts: 32, max_sweeps: 100, max_setup_terms: 1000000 } }
}
struct Zero;
impl Storage for Zero {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], g: &mut [f64]) { g.fill(0.0); }
}

#[test]
fn helix_mass_counts_metal_arc_length_without_inventing_coil_stiffness() {
    let h = HelicalWire { wire_radius_m: 0.00015, coil_radius_m: 0.00055,
        pitch_m: 0.00085, density_kg_m3: 7800.0 };
    let mut length = 0.0;
    let n = 8192;
    let point = |i: usize| { let t = i as f64/n as f64;
        [h.coil_radius_m*(std::f64::consts::TAU*t).cos(),
         h.coil_radius_m*(std::f64::consts::TAU*t).sin(), h.pitch_m*t] };
    for i in 0..n { let a=point(i); let b=point(i+1);
        length += (0..3).map(|j| (b[j]-a[j]).powi(2)).sum::<f64>().sqrt(); }
    let measured_polyline = h.density_kg_m3*std::f64::consts::PI*h.wire_radius_m.powi(2)*length/h.pitch_m;
    let analytic = h.linear_density_kg_m().unwrap();
    assert!((analytic-measured_polyline).abs()/analytic < 3e-8);
    let straight = HelicalWire { coil_radius_m: 0.0, ..h }.linear_density_kg_m().unwrap();
    assert!(analytic > 4.0*straight);
    assert!(HelicalWire { pitch_m: 0.0, ..h }.linear_density_kg_m().is_err());
}

#[test]
fn wire_frequencies_follow_tension_mass_length_and_flexure_in_the_existing_owner() {
    let w=wire();
    let initial=vec![ModalAcousticState::default();3];
    let body=w.body(initial.clone()).unwrap();
    let BodyPotential::Linear(omega)=body.potential else {panic!("linear wire expected")};
    for (i,actual) in omega.iter().enumerate() {
        let k=(i+1) as f64*std::f64::consts::PI/w.length_m();
        let expected=((w.tension_n*k*k+w.bending_stiffness_n_m2*k.powi(4))/w.linear_density_kg_m).sqrt();
        assert!((actual-expected).abs() < 1e-12*expected);
    }
    let stiff=WireSpan { tension_n:2.8, ..w.clone() }.body(initial).unwrap();
    let BodyPotential::Linear(high)=stiff.potential else {panic!()};
    assert!(high[0] > 1.99*omega[0]);
    assert!(w.body(vec![]).is_err());
}

#[test]
fn line_contact_preserves_integrated_strength_and_reciprocal_work_on_refinement() {
    let w=wire();
    for cells in [6,12,24,48] {
        let l=line(cells,5e8);
        let obstacle=w.contact(&l,&vec![vec![10.0];cells],0..1,1..4,4).unwrap();
        let storage=ContactStorage::new(Box::new(Zero),4,vec![obstacle]).unwrap();
        let x=[1e-5,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
        let penetration=10.0*x[0]-l.gaps_m[0];
        let expected=w.length_m()*l.stiffness_per_length*penetration.powf(l.alpha+1.0)/(l.alpha+1.0);
        assert!((storage.hamiltonian(&x)-expected).abs() < 1e-13*expected);
        let direction=[0.3,0.0,-0.2,0.0,0.1,0.0,0.4,0.0];
        let mut gradient=[0.0;8];storage.gradient(&x,&mut gradient);
        let eps=1e-9;let mut plus=x;let mut minus=x;
        for i in 0..8 {plus[i]+=eps*direction[i];minus[i]-=eps*direction[i];}
        let fd=(storage.hamiltonian(&plus)-storage.hamiltonian(&minus))/(2.0*eps);
        let work=gradient.iter().zip(direction).map(|(g,v)|g*v).sum::<f64>();
        assert!((fd-work).abs() < 1e-5*work.abs().max(1.0));
    }
}

#[test]
fn contact_stations_use_film_interpolation_not_nearest_node_snapping() {
    let film=TensionedDisk::new(TensionedDiskSpec {radius_m:0.1,thickness_m:0.0002,young_pa:4e9,
        poisson:0.38,density_kg_m3:1390.0,tension_n_m:1000.0,radial_intervals:3,azimuths:12},
        ProfileBudget{max_nodes:100,max_triangles:200,max_feature_evaluations:0}).unwrap();
    // An affine kinematic fixture on interior nodes, not a claimed eigenpair.
    let mut phi=vec![0.0;film.model.free];
    for (i,&(x,y)) in film.mesh.nodes.iter().enumerate() {
        if let Some(k)=film.model.dof_map[3*i] {phi[k]=0.3+0.2*x-0.1*y;}
    }
    let modes=[ModePair{lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)}];
    let points=[[0.002,0.003],[0.004,0.003],[-0.005,0.007],[0.0,0.0]];
    let rows=film_shapes(&film,&modes,&points).unwrap();
    for (row,[x,y]) in rows.iter().zip(points) {assert!((row[0]-(0.3+0.2*x-0.1*y)).abs()<1e-14);}
    assert_ne!(rows[0][0],rows[1][0]);
    assert!(film_shapes(&film,&modes,&[[0.11,0.0]]).is_err());
}

fn coupled(k:f64,steps:u64)->LinearImpactSystem {
    let w=wire();let l=line(6,k);
    let contact=w.contact(&l,&vec![vec![10.0];6],0..1,1..4,4).unwrap();
    let receiver=ImpactBody{potential:BodyPotential::Linear(vec![1000.0]),
        initial:vec![ModalAcousticState{displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:0.02}],
        damping_per_s:vec![2.0]};
    LinearImpactSystem::new(vec![receiver,w.body(vec![ModalAcousticState::default();3]).unwrap()],
        vec![contact],vec![],config(steps),&CancelGate::new_clock_free()).unwrap()
}
#[test]
fn unforced_wires_rattle_and_react_back_on_the_head_without_a_noise_source() {
    let gate=CancelGate::new_clock_free();let mut active=coupled(5e8,2048);let mut off=coupled(0.0,2048);
    let mut wire_peak=0.0_f64;let mut changed=0.0_f64;
    for _ in 0..2048 {
        let f=active.step(&[0.0;4],&gate).unwrap();off.step(&[0.0;4],&gate).unwrap();
        assert!(f.balance_residual_j.abs()<1e-8);
        wire_peak=wire_peak.max(active.state()[3].abs());
        changed=changed.max((active.state()[1]-off.state()[1]).abs());
        assert!(off.state()[2..].iter().all(|x|*x==0.0));
    }
    assert!(wire_peak>1e-5 && changed>1e-5,"wire motion and back reaction must both be physical");
}
#[test]
fn cancellation_and_refusal_keep_wire_vibration_and_the_joint_contact_clock() {
    let gate=CancelGate::new_clock_free();let mut a=coupled(5e8,64);let mut b=coupled(5e8,128);
    for _ in 0..32 {a.step(&[0.0;4],&gate).unwrap();b.step(&[0.0;4],&gate).unwrap();}
    let x=a.state().to_vec();let f=*a.frame();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(a.step(&[0.0;4],&cancelled).is_err());
    assert!(a.step(&[f64::NAN;4],&gate).is_err());
    assert_eq!(a.state(),x);assert_eq!(*a.frame(),f);
    a.extend_step_budget(128).unwrap();
    for _ in 32..128 {a.step(&[0.0;4],&gate).unwrap();b.step(&[0.0;4],&gate).unwrap();}
    assert_eq!(a.state(),b.state());assert_eq!(a.frame(),b.frame());
}
