//! A compliant stop can support a genuinely free assembly without a tether.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_couple::render::schedule::force::coupled::contact::{
    ContactModalSystem, ModalContact, ModalContactConfig,
};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn budget() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 16, max_connections: 8, max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 1000.0, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-10,
        energy_relative_tolerance: 1e-8,
    }
}
fn contact_budget() -> ModalContactConfig {
    ModalContactConfig {
        max_iterations: 96, maximum_force_n: 1e4, maximum_penetration_m: 0.05,
        force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-11,
    }
}
fn mass(rate: u32, kg: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_free_mass(rate, kg, 0.0, 0.0,
        ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn resonator(rate: u32) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: 800.0, damping_ratio: 0.02,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn attachment(component: usize, shape: f64) -> ModalAttachment {
    ModalAttachment { component, shapes: vec![shape] }
}
fn point(left: ModalAttachment, right: ModalAttachment, gap: f64, k: f64, weight: f64, alpha: f64)
    -> ModalContact
{
    ModalContact { left, right, law: Obstacle::new(
        vec![-1.0], 1, 1, vec![gap], vec![weight], k, alpha,
        "authored free-assembly stop".into()).unwrap().with_internal_loss(0.2).unwrap() }
}
fn free(rate: u32, kg: f64, c: ModalCouplingConfig) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![mass(rate, kg)], vec![], c, &CancelGate::new()).unwrap()
}
fn states(s: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m| m.states()).copied().collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:e} != {b:e} (tolerance {tolerance:e})");
}

#[test]
fn g1_hertz_stop_supports_free_mass_with_correct_storage_held_silence_and_release() {
    for (rate, kg, force) in [
        (24000, 0.04_f64, 0.7_f64), (48000, 0.07, 1.3), (96000, 0.13, 0.9),
    ] {
        let root = kg.sqrt();
        let mut s = free(rate, kg, budget());
        let stop = point(attachment(0, 1.0/root), attachment(0, 0.0), 0.0004, 3e6, 0.8, 1.5);
        let load = [force/root];
        // The ordinary compliance path still cannot invent a bilateral support.
        assert!(s.initialize_contact_equilibrium(&load, &[(stop.clone(), contact_budget())],
            MultiContactConfig { max_contacts: 1, max_sweeps: 64, max_setup_terms: 100_000 },
            &CancelGate::new()).is_err());
        let energy = s.initialize_contact_equilibrium_with_gauge(
            &load, &stop, contact_budget(), &[1.0], &CancelGate::new()).unwrap();
        let penetration = (force/(0.8*3e6)).powf(2.0/3.0);
        let x = s.components()[0].states()[0].displacement_m_sqrt_kg/root;
        close(x, 0.0004+penetration, 2e-13);
        close(0.8*3e6*(x-0.0004).powf(1.5), force, 2e-10);
        close(energy, force*penetration/2.5, 2e-13);
        assert_eq!(s.components()[0].modes()[0].angular_frequency_rad_s, 0.0);
        assert_eq!(s.samples_rendered(), 0); assert!(s.last_frame().is_none());
        let mut contacted = ContactModalSystem::new(s, stop, contact_budget(), &CancelGate::new()).unwrap();
        for _ in 0..64 { contacted.step(&load).unwrap(); }
        assert!(contacted.components()[0].states()[0].velocity_m_sqrt_kg_per_s.abs() < 1e-8);
        for _ in 0..64 { contacted.step(&[0.0]).unwrap(); }
        assert!(contacted.components()[0].states()[0].velocity_m_sqrt_kg_per_s < 0.0);
        assert!(contacted.components()[0].states()[0].displacement_m_sqrt_kg/root < x);
        assert!(contacted.total_energy_j().unwrap() <= energy+1e-9);
    }
}

#[test]
fn g1_stop_carries_both_unequal_masses_and_the_bilateral_rest_length() {
    for masses in [[0.04_f64, 0.09_f64], [0.07, 0.13]] {
        let roots = masses.map(f64::sqrt);
        let make = || CoupledModalSystem::new(
            vec![mass(48000, masses[0]), mass(48000, masses[1])],
            vec![ModalConnection {
                left: attachment(0, 1.0/roots[0]), right: attachment(1, 1.0/roots[1]),
                stiffness_n_m: 1200.0, damping_n_s_m: 4.0, rest_extension_m: 0.0012,
            }], budget(), &CancelGate::new()).unwrap();
        let stop = point(attachment(1, 1.0/roots[1]), attachment(1, 0.0), 0.0003, 7e6, 0.6, 2.0);
        let loads = [0.8/roots[0], 1.1/roots[1]];
        let norm = (masses[0]+masses[1]).sqrt();
        let gauge = [roots[0]/norm, roots[1]/norm];
        let mut reference = Vec::new();
        for sign in [1.0, -1.0] {
            let mut s = make();
            let energy = s.initialize_contact_equilibrium_with_gauge(
                &loads, &stop, contact_budget(), &gauge.map(|v| sign*v), &CancelGate::new()).unwrap();
            let q = states(&s);
            let x = [q[0].displacement_m_sqrt_kg/roots[0], q[1].displacement_m_sqrt_kg/roots[1]];
            let penetration = (1.9_f64/(0.6*7e6)).sqrt();
            close(x[1], 0.0003+penetration, 2e-13);
            close(x[0], x[1]+0.0012+0.8/1200.0, 2e-13);
            let spring = 1200.0*(x[0]-x[1]-0.0012);
            let reaction = 0.6*7e6*(x[1]-0.0003).powi(2);
            // Original full generalized forces, with no gauge reaction.
            close(loads[0]-spring/roots[0], 0.0, 2e-9);
            close(loads[1]+spring/roots[1]-reaction/roots[1], 0.0, 2e-9);
            close(energy, 0.8_f64.powi(2)/2400.0+1.9*penetration/3.0, 2e-13);
            if reference.is_empty() { reference = q; } else {
                for (a, b) in reference.iter().zip(q) {
                    close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg, 1e-16);
                }
            }
            let mut contacted = ContactModalSystem::new(s, stop.clone(), contact_budget(), &CancelGate::new()).unwrap();
            for _ in 0..32 { contacted.step(&loads).unwrap(); }
            for model in contacted.components() {
                assert!(model.states()[0].velocity_m_sqrt_kg_per_s.abs() < 1e-8);
            }
        }
    }
}

#[test]
fn g1_reciprocal_stop_load_reaches_the_original_elastic_receiver() {
    let root = 0.07_f64.sqrt();
    let mut s = CoupledModalSystem::new(vec![mass(48000, 0.07), resonator(48000)],
        vec![], budget(), &CancelGate::new()).unwrap();
    let stop = point(attachment(0, 1.0/root), attachment(1, 2.0), 0.0002, 2e8, 0.75, 2.0);
    let loads = [2.5/root, 0.3];
    let energy = s.initialize_contact_equilibrium_with_gauge(
        &loads, &stop, contact_budget(), &[1.0, 0.0], &CancelGate::new()).unwrap();
    let q = states(&s);
    let receiver = (0.3+2.0*2.5)/800.0_f64.powi(2);
    let penetration = (2.5_f64/(0.75*2e8)).sqrt();
    close(q[1].displacement_m_sqrt_kg, receiver, 1e-16);
    close(q[0].displacement_m_sqrt_kg/root, 0.0002+penetration+2.0*receiver, 2e-13);
    let closure = q[0].displacement_m_sqrt_kg/root-2.0*q[1].displacement_m_sqrt_kg;
    let reaction = 0.75*2e8*(closure-0.0002).powi(2);
    close(loads[0]-reaction/root, 0.0, 2e-9);
    close(800.0_f64.powi(2)*q[1].displacement_m_sqrt_kg-loads[1]-2.0*reaction, 0.0, 2e-9);
    close(energy, 0.5*800.0_f64.powi(2)*receiver*receiver+2.5*penetration/3.0, 2e-13);
}

#[test]
fn g1_original_load_cancellation_roundoff_is_admitted_without_projecting_force() {
    let kg = 0.024_f64; let root = kg.sqrt(); let force = 0.14800000000000002;
    let load = force/root; let shape = 1.0/root;
    // The reaction multiplication leaves a nonzero floating-point remainder.
    // Normalizing this remainder by itself would falsely refuse this mass.
    assert_ne!(load-shape*(load/shape), 0.0);
    let mut s = free(48000, kg, budget());
    let stop = point(attachment(0, shape), attachment(0, 0.0), 0.0001, 1e6, 1.0, 2.0);
    s.initialize_contact_equilibrium_with_gauge(
        &[load], &stop, contact_budget(), &[1.0], &CancelGate::new()).unwrap();
    let x = states(&s)[0].displacement_m_sqrt_kg/root;
    close(load-shape*1e6*(x-0.0001).powi(2), 0.0, 2e-9);
}

#[test]
fn g1_zero_resultant_keeps_an_open_pose_and_relaxes_interference_to_touching() {
    for gap in [0.001, 0.0, -0.001] {
        let mut s = free(48000, 0.04, budget());
        let stop = point(attachment(0, 5.0), attachment(0, 0.0), gap, 1e6, 1.0, 1.5);
        let energy = s.initialize_contact_equilibrium_with_gauge(
            &[0.0], &stop, contact_budget(), &[1.0], &CancelGate::new()).unwrap();
        close(states(&s)[0].displacement_m_sqrt_kg/0.2, gap.min(0.0), 1e-16);
        close(energy, 0.0, 1e-20);
        assert_eq!(states(&s)[0].velocity_m_sqrt_kg_per_s, 0.0);
    }
}

#[test]
fn g0_tension_false_nullspaces_and_contacts_blind_to_rigid_motion_refuse() {
    let stop = point(attachment(0, 5.0), attachment(0, 0.0), 0.001, 1e6, 1.0, 2.0);
    let mut s = free(48000, 0.04, budget()); let initial = states(&s);
    for load in [-1.0, -f64::MIN_POSITIVE, -f64::from_bits(1)] {
        assert!(s.initialize_contact_equilibrium_with_gauge(
            &[load], &stop, contact_budget(), &[1.0], &CancelGate::new()).is_err());
        assert_eq!(states(&s), initial);
    }
    for gauge in [vec![], vec![0.0], vec![2.0], vec![f64::NAN], vec![1.0, 0.0]] {
        assert!(s.initialize_contact_equilibrium_with_gauge(
            &[5.0], &stop, contact_budget(), &gauge, &CancelGate::new()).is_err());
        assert_eq!(states(&s), initial);
    }
    // Even a weak real ground spring breaks translation symmetry. It cannot
    // be erased by declaring a convenient contact-support direction.
    let mut anchored = CoupledModalSystem::new(vec![mass(48000, 0.04)],
        vec![ModalConnection { left: attachment(0, 5.0), right: attachment(0, 0.0),
            stiffness_n_m: 1e-12, damping_n_s_m: 0.0, rest_extension_m: 0.0 }],
        budget(), &CancelGate::new()).unwrap();
    assert!(anchored.initialize_contact_equilibrium_with_gauge(
        &[5.0], &stop, contact_budget(), &[1.0], &CancelGate::new()).is_err());
    let mut disconnected = CoupledModalSystem::new(vec![mass(48000, 0.04), mass(48000, 0.09)],
        vec![], budget(), &CancelGate::new()).unwrap();
    assert!(disconnected.initialize_contact_equilibrium_with_gauge(
        &[5.0, 0.0], &stop, contact_budget(), &[1.0, 0.0], &CancelGate::new()).is_err());
    let mut pair = CoupledModalSystem::new(vec![mass(48000, 0.04), mass(48000, 0.09)],
        vec![ModalConnection { left: attachment(0, 5.0), right: attachment(1, 1.0/0.3),
            stiffness_n_m: 1000.0, damping_n_s_m: 0.0, rest_extension_m: 0.0 }],
        budget(), &CancelGate::new()).unwrap();
    let internal = point(attachment(0, 5.0), attachment(1, 1.0/0.3), 0.0, 1e6, 1.0, 2.0);
    let norm = 0.13_f64.sqrt();
    assert!(pair.initialize_contact_equilibrium_with_gauge(
        &[5.0, 0.0], &internal, contact_budget(), &[0.2/norm, 0.3/norm], &CancelGate::new()).is_err());
    assert!(states(&pair).iter().all(|s| s.displacement_m_sqrt_kg == 0.0));
}

#[test]
fn g0_contact_caps_work_cancellation_and_energy_refusals_are_retryable() {
    let stop = point(attachment(0, 5.0), attachment(0, 0.0), 0.0003, 1e6, 1.0, 2.0);
    let c = contact_budget();
    for config in [
        ModalContactConfig { maximum_force_n: 0.5, ..c },
        ModalContactConfig { maximum_penetration_m: 0.0005, ..c },
        ModalContactConfig { max_iterations: 1, ..c },
    ] {
        let mut s = free(48000, 0.04, budget()); let initial = states(&s);
        assert!(s.initialize_contact_equilibrium_with_gauge(
            &[5.0], &stop, config, &[1.0], &CancelGate::new()).is_err());
        assert_eq!(states(&s), initial); assert_eq!(s.samples_rendered(), 0); assert!(s.last_frame().is_none());
        s.initialize_contact_equilibrium_with_gauge(
            &[0.0], &stop, config, &[1.0], &CancelGate::new()).unwrap();
    }
    let mut s = free(48000, 0.04, ModalCouplingConfig { maximum_total_energy_j: 1e-5, ..budget() });
    let initial = states(&s);
    assert!(s.initialize_contact_equilibrium_with_gauge(
        &[5.0], &stop, c, &[1.0], &CancelGate::new()).is_err());
    assert_eq!(states(&s), initial);
    let gate = CancelGate::new(); gate.request();
    assert!(matches!(s.initialize_contact_equilibrium_with_gauge(
        &[0.05], &stop, c, &[1.0], &gate), Err(ModalCouplingError::Cancelled)));
    assert_eq!(states(&s), initial);
    s.initialize_contact_equilibrium_with_gauge(
        &[0.05], &stop, c, &[1.0], &CancelGate::new()).unwrap();
    let mut short = free(48000, 0.04, ModalCouplingConfig { max_setup_terms: 6, ..budget() });
    assert!(short.initialize_contact_equilibrium_with_gauge(
        &[5.0], &stop, c, &[1.0], &CancelGate::new()).is_err());
    // Original Schur/nullspace work for one unconnected free mode is seven terms.
    free(48000, 0.04, ModalCouplingConfig { max_setup_terms: 7, ..budget() })
        .initialize_contact_equilibrium_with_gauge(&[5.0], &stop, c, &[1.0], &CancelGate::new()).unwrap();
}
