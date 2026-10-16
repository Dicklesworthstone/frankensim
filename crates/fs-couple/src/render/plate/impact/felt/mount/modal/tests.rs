use super::*;
use crate::render::plate::impact::{BodyPotential, ImpactConfig, ImpactSystem};
use crate::render::plate::impact::felt::{KelvinBranch, mount::{Carrier, MountSite, MountedFelt}};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_material::fiber::WoolFelt;
use fs_exec::CancelGate;

fn pad(row: Vec<f64>) -> FeltPad {
    FeltPad { area_m2: 0.0002, thickness_m: 0.006, precompression_m: 0.0006,
        weights: row, law: WoolFelt::new(30000.0, 0.2, 2.2, 3.0, 0.15, 0.7).unwrap(),
        prior_maximum_strain: 0.1,
        creep: vec![KelvinBranch { stiffness_n_m: 1500.0, viscosity_n_s_m: 8.0 }] }
}
fn body(omega: Vec<f64>) -> ImpactBody {
    let n = omega.len();
    ImpactBody { potential: BodyPotential::Linear(omega),
        initial: vec![ModalAcousticState::default(); n], damping_per_s: vec![0.0; n] }
}
fn config() -> ImpactConfig {
    ImpactConfig { dt_s: 2e-6, max_steps: 1000, maximum_energy_j: 2.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-7,
        maximum_generalized_force: 1e4 }
}
fn supported() -> ImpactSystem {
    let mut support = body(vec![80.0, 110.0, 90.0]);
    support.damping_per_s = vec![1.0, 2.0, 3.0];
    let sites = vec![
        ModalMountSite { pad: pad(vec![1.0]), support_weights: vec![-2.0, -0.4, 0.3] },
        ModalMountSite { pad: pad(vec![-1.0]), support_weights: vec![2.0, -0.1, 0.5] },
    ];
    let mount = ModalMountedFelt::new(1, support, &sites).unwrap();
    let mut receiver = body(vec![500.0]);
    receiver.initial[0].velocity_m_sqrt_kg_per_s = 0.01;
    ImpactSystem::new(vec![receiver, mount.body], vec![], mount.pads, vec![], config()).unwrap()
}
#[test]
fn one_mode_support_reproduces_the_original_translating_mount() {
    let carrier = Carrier { mass_kg: 0.25, stiffness_n_m: 1600.0, resistance_n_s_m: 0.5,
        initial_position_m: 0.0001, initial_velocity_m_s: 0.02 };
    let site = MountSite { pad: pad(vec![1.0]), compression_per_carrier_m: -1.0 };
    let old = MountedFelt::new(1, carrier, &[site]).unwrap();
    let new = ModalMountedFelt::new(1, old.body.clone(), &[ModalMountSite {
        pad: pad(vec![1.0]), support_weights: vec![-2.0] }]).unwrap();
    assert_eq!(old.pads[0].weights, new.pads[0].weights);
    let build = |support: ImpactBody, pads: Vec<FeltPad>| {
        ImpactSystem::new(vec![body(vec![500.0]), support], vec![], pads, vec![], config())
            .unwrap().prepare_analytic().unwrap()
    };
    let mut a = build(old.body, old.pads); let mut b = build(new.body, new.pads);
    let gate = CancelGate::new_clock_free();
    for _ in 0..100 {
        a.step(&[0.0; 2], &gate).unwrap(); b.step(&[0.0; 2], &gate).unwrap();
        assert_eq!(a.state(), b.state());
    }
}
#[test]
fn signed_nonunit_point_rows_preserve_force_and_moment_virtual_work() {
    let row = vec![-2.0, -4.0, 7.0];
    let mount = ModalMountedFelt::new(2, body(vec![0.0, 10.0, 20.0]), &[ModalMountSite {
        pad: pad(vec![0.25, -0.5]), support_weights: row.clone() }]).unwrap();
    assert_eq!(&mount.pads[0].weights[2..], row.as_slice());
    let rates = [0.1, -0.2, 0.3, -0.4, 0.5];
    let force = 3.7;
    let speed: f64 = mount.pads[0].weights.iter().zip(rates).map(|(b,v)| b*v).sum();
    let generalized_power: f64 = mount.pads[0].weights.iter().zip(rates).map(|(b,v)| -force*b*v).sum();
    assert!((generalized_power + force*speed).abs() < 1e-13);
}
#[test]
fn coupled_motion_excites_both_support_rotations_and_keeps_the_energy_ledger() {
    let mut system = supported().prepare_analytic().unwrap();
    let initial = system.stored_energy_j(); let mut net = 0.0;
    let gate = CancelGate::new_clock_free();
    for _ in 0..400 {
        let frame = system.step(&[0.0; 4], &gate).unwrap();
        assert!(frame.dissipated_energy_j >= 0.0);
        assert!(frame.balance_residual_j.abs() < 1e-9);
        net += frame.supplied_work_j - frame.dissipated_energy_j;
    }
    assert!(system.state()[5].abs() > 1e-9 && system.state()[7].abs() > 1e-9);
    assert!(system.state()[8] != 0.0 && system.state()[9] != 0.0);
    assert!((system.stored_energy_j() - initial - net).abs() < 1e-8);
}
#[test]
fn cancellation_and_bad_force_do_not_publish_support_or_felt_memory() {
    let mut retry = supported().prepare_analytic().unwrap();
    let mut clean = supported().prepare_analytic().unwrap();
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    let before = retry.state().to_vec();
    assert!(retry.step(&[0.0; 4], &cancelled).is_err());
    assert_eq!(retry.state(), before.as_slice());
    let gate = CancelGate::new_clock_free();
    assert!(retry.step(&[f64::NAN, 0.0, 0.0, 0.0], &gate).is_err());
    assert_eq!(retry.state(), before.as_slice());
    let a = retry.step(&[0.0; 4], &gate).unwrap();
    let b = clean.step(&[0.0; 4], &gate).unwrap();
    assert_eq!(a.sample, b.sample); assert_eq!(retry.state(), clean.state());
}
#[test]
fn disconnected_rows_bad_shapes_and_negative_resistance_are_refused() {
    let mut site = ModalMountSite { pad: pad(vec![1.0]), support_weights: vec![0.0; 3] };
    assert!(ModalMountedFelt::new(1, body(vec![0.0; 3]), &[site.clone()]).is_err());
    site.support_weights[1] = 1.0;
    assert!(ModalMountedFelt::new(1, body(vec![0.0; 2]), &[site.clone()]).is_err());
    let mut support = body(vec![0.0; 3]); support.damping_per_s[0] = -1.0;
    assert!(ModalMountedFelt::new(1, support, &[site.clone()]).is_err());
    site.support_weights[0] = f64::NAN;
    assert!(ModalMountedFelt::new(1, body(vec![0.0; 3]), &[site]).is_err());
}
