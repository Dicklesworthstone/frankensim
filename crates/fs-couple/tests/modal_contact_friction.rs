//! Shared-body friction: independent force/work checks, not producer-only totals.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use fs_couple::render::schedule::force::coupled::contact::multiple::friction::ModalFriction;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:e} != {b:e}; tolerance {tolerance:e}");
}
fn normal_config() -> ModalContactConfig {
    ModalContactConfig { max_iterations: 80, maximum_force_n: 1e5, maximum_penetration_m: 0.01,
        force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-10 }
}
fn config() -> MultiContactConfig {
    MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 4096 }
}
fn image(velocity: f64) -> ModalAcousticTimeModel {
    let modes = [800.0, 600.0].into_iter().map(|omega| ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }).collect();
    let mut model = ModalAcousticTimeModel::try_new(48000, modes, ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState::default(), ModalAcousticState {
        displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: velocity,
    }]).unwrap();
    model
}
fn attachment(component: usize, shapes: [f64; 2]) -> ModalAttachment {
    ModalAttachment { component, shapes: shapes.to_vec() }
}
fn contact(left: usize, right: usize, gap: f64) -> ModalContact {
    ModalContact { left: attachment(left, [1.0, 0.0]), right: attachment(right, [1.0, 0.0]),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], 1e5, 1.0,
            "synthetic:linear-modal-contact".into()).unwrap() }
}
fn friction(mu: f64, mixed: f64) -> ModalFriction {
    ModalFriction { left_shapes: vec![mixed, 1.0], right_shapes: vec![mixed, 1.0],
        coefficient: mu, regularization_speed_m_s: 0.01, maximum_force_n: 1e4,
        source: "synthetic:constant-mu-not-a-material-card".into() }
}
fn build(velocities: &[f64], contacts: Vec<ModalContact>, bilateral: bool, cfg: MultiContactConfig)
    -> MultiContactModalSystem
{
    let network_config = ModalCouplingConfig {
        max_modes: 8, max_connections: 4, max_setup_terms: 4096, nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6, solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8,
    };
    let links = if bilateral { vec![ModalConnection {
        left: attachment(0, [0.5, 1.0]), right: attachment(velocities.len()-1, [0.5, 1.0]),
        stiffness_n_m: 8e4, damping_n_s_m: 0.4, rest_extension_m: 0.0,
    }] } else { vec![] };
    let network = CoupledModalSystem::new(velocities.iter().map(|&v| image(v)).collect(),
        links, network_config, &CancelGate::new()).unwrap();
    MultiContactModalSystem::new(network, contacts.into_iter().map(|c| (c, normal_config())).collect(),
        cfg, &CancelGate::new()).unwrap()
}
fn pair(v: f64, gap: f64, cfg: MultiContactConfig) -> MultiContactModalSystem {
    build(&[v, 0.0], vec![contact(0, 1, gap)], false, cfg)
}
fn shared() -> MultiContactModalSystem {
    build(&[1.0, 0.2, -0.3], vec![contact(0, 1, -0.0006), contact(1, 2, -0.0006)], true, config())
}
fn states(system: &MultiContactModalSystem) -> Vec<ModalAcousticState> {
    system.components().iter().flat_map(|m| m.states().iter().copied()).collect()
}
fn displacement(state: &[ModalAcousticState], i: usize) -> f64 { state[i].displacement_m_sqrt_kg }

#[test]
fn one_step_matches_independent_closed_form_in_ramp_slide_rest_and_reversal() {
    for v in [-1.0, -0.005, 0.0, 0.005, 1.0] {
        let mut system = pair(v, -0.001, config()).with_friction(vec![Some(friction(0.3, 0.0))], &CancelGate::new()).unwrap();
        let dt = system.sample_period_s();
        let dn = (1.0-(800.0*dt).cos())/(800.0*800.0);
        let dtan = (1.0-(600.0*dt).cos())/(600.0*600.0);
        let normal = 100.0/(1.0+1e5*dn);
        let free_y = v*(600.0*dt).sin()/600.0;
        let ramp = 0.3*normal*free_y/(dt*0.01+0.3*normal*2.0*dtan);
        let expected = ramp.clamp(-0.3*normal, 0.3*normal);
        let delta_y = free_y-2.0*dtan*expected;
        let frame = system.step(&[0.0; 4]).unwrap().clone();
        let tangent = frame.friction[0].unwrap();
        close(frame.contacts[0].normal_force_n, normal, 2e-8);
        close(tangent.reaction_n, expected, 2e-8);
        close(tangent.slip_velocity_m_s, delta_y/dt, 1e-10);
        close(tangent.dissipation_j, expected*delta_y, 1e-11);
        let actual = states(&system);
        close(displacement(&actual, 0), -dn*normal, 1e-12);
        close(displacement(&actual, 2), dn*normal, 1e-12);
        close(displacement(&actual, 1), free_y-dtan*expected, 1e-12);
        close(displacement(&actual, 3), dtan*expected, 1e-12);
        close(actual[3].velocity_m_sqrt_kg_per_s, expected*(600.0*dt).sin()/600.0, 1e-11);
        assert!(tangent.dissipation_j >= 0.0);
    }
}

#[test]
fn shared_contacts_recheck_actual_motion_loads_and_independent_work() {
    let mut system = shared().with_friction(vec![Some(friction(0.3, 0.4)), Some(friction(0.5, -0.2))], &CancelGate::new()).unwrap();
    let mut normal_only = shared();
    let initial = system.total_energy_j().unwrap();
    let (mut work, mut loss, mut friction_loss, mut mixed_effect) = (0.0, 0.0, 0.0, 0.0_f64);
    let mut both_active = false;
    for sample in 0..160 {
        let old = states(&system);
        let external = [0.0, if sample < 40 { 0.5 } else { 0.0 }, 0.0, 0.0, 0.0, -0.25];
        let frame = system.step(&external).unwrap().clone();
        let reference = normal_only.step(&external).unwrap();
        let new = states(&system);
        let mut independent_work = 0.0;
        for i in 0..6 { independent_work += external[i]*(displacement(&new, i)-displacement(&old, i)); }
        close(frame.external_work_j, independent_work, 1e-12);
        let mut tangent_sum = 0.0;
        for (i, mu, mixed) in [(0, 0.3, 0.4), (1, 0.5, -0.2)] {
            let x0 = displacement(&old, 2*i)-displacement(&old, 2*i+2);
            let x1 = displacement(&new, 2*i)-displacement(&new, 2*i+2);
            let p0 = (x0+0.0006).max(0.0);
            let p1 = (x1+0.0006).max(0.0);
            let elastic = if x0 == x1 { 1e5*p0 } else { 0.5e5*(p1*p1-p0*p0)/(x1-x0) };
            close(frame.contacts[i].normal_force_n, elastic, 1e-6);
            let tangent = frame.friction[i].unwrap();
            let delta_y = mixed*(x1-x0) + displacement(&new, 2*i+1)-displacement(&old, 2*i+1)
                - (displacement(&new, 2*i+3)-displacement(&old, 2*i+3));
            let slip = delta_y/system.sample_period_s();
            let expected = mu*frame.contacts[i].normal_force_n*(slip/0.01).clamp(-1.0, 1.0);
            close(tangent.reaction_n, expected, tangent.force_tolerance_n+1e-9);
            close(tangent.dissipation_j, tangent.reaction_n*delta_y, 1e-12);
            assert!(tangent.dissipation_j >= 0.0);
            tangent_sum += tangent.dissipation_j;
            mixed_effect = mixed_effect.max((frame.contacts[i].normal_force_n-reference.contacts[i].normal_force_n).abs());
        }
        close(frame.friction_dissipation_j, tangent_sum, 1e-12);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        both_active |= frame.contacts.iter().all(|c| c.normal_force_n > 1.0)
            && frame.friction.iter().all(|c| c.unwrap().reaction_n.abs() > 0.1);
        work += frame.external_work_j;
        loss += frame.network_dissipation_j+frame.contact_dissipation_j+frame.friction_dissipation_j;
        friction_loss += frame.friction_dissipation_j;
    }
    assert!(both_active && friction_loss > 1e-4 && mixed_effect > 1e-4);
    close(system.total_energy_j().unwrap()-initial+loss-work, 0.0, 1e-7);
}

#[test]
fn zero_coefficient_is_bitwise_normal_only_and_separation_cannot_drag() {
    for separated in [false, true] {
        let gap = if separated { 1.0 } else { -0.001 };
        let mu = if separated { 0.5 } else { 0.0 };
        let mut reference = pair(1.0, gap, config());
        let mut system = pair(1.0, gap, config()).with_friction(vec![Some(friction(mu, 0.3))], &CancelGate::new()).unwrap();
        for _ in 0..100 {
            let a = reference.step(&[0.0; 4]).unwrap().clone();
            let b = system.step(&[0.0; 4]).unwrap().clone();
            assert_eq!(states(&reference), states(&system));
            assert_eq!(a.contacts, b.contacts);
            assert_eq!(a.observer_pressure_pa.to_bits(), b.observer_pressure_pa.to_bits());
            assert_eq!(a.sweeps, b.sweeps);
            assert_eq!(b.friction[0].unwrap().reaction_n, 0.0);
            assert_eq!(b.friction_dissipation_j, 0.0);
        }
    }
}

#[test]
fn reversing_the_authored_tangent_changes_no_physical_state_or_work() {
    let positive = vec![Some(friction(0.3, 0.4)), Some(friction(0.5, -0.2))];
    let mut negative = positive.clone();
    for spec in negative.iter_mut().flatten() {
        for x in spec.left_shapes.iter_mut().chain(&mut spec.right_shapes) { *x = -*x; }
    }
    let mut a = shared().with_friction(positive, &CancelGate::new()).unwrap();
    let mut b = shared().with_friction(negative, &CancelGate::new()).unwrap();
    for _ in 0..80 {
        let fa = a.step(&[0.0; 6]).unwrap().clone();
        let fb = b.step(&[0.0; 6]).unwrap().clone();
        assert_eq!(states(&a), states(&b));
        assert_eq!(fa.contacts, fb.contacts);
        for i in 0..2 {
            let (ta, tb) = (fa.friction[i].unwrap(), fb.friction[i].unwrap());
            close(ta.reaction_n, -tb.reaction_n, 1e-12);
            close(ta.slip_velocity_m_s, -tb.slip_velocity_m_s, 1e-12);
            close(ta.dissipation_j, tb.dissipation_j, 1e-12);
        }
    }
}

#[test]
fn contact_permutation_preserves_the_converged_physical_solution() {
    let mut a = shared().with_friction(vec![Some(friction(0.3, 0.4)), Some(friction(0.5, -0.2))], &CancelGate::new()).unwrap();
    let mut b = build(&[1.0, 0.2, -0.3], vec![contact(1, 2, -0.0006), contact(0, 1, -0.0006)], true, config())
        .with_friction(vec![Some(friction(0.5, -0.2)), Some(friction(0.3, 0.4))], &CancelGate::new()).unwrap();
    for _ in 0..80 {
        a.step(&[0.0; 6]).unwrap(); b.step(&[0.0; 6]).unwrap();
        for (x, y) in states(&a).iter().zip(states(&b)) {
            close(x.displacement_m_sqrt_kg, y.displacement_m_sqrt_kg, 1e-10);
            close(x.velocity_m_sqrt_kg_per_s, y.velocity_m_sqrt_kg_per_s, 1e-8);
        }
    }
}

#[test]
fn cancellation_invalid_input_and_friction_ceiling_leave_no_partial_state() {
    let make = || pair(1.0, -0.001, config()).with_friction(vec![Some(friction(0.3, 0.2))], &CancelGate::new()).unwrap();
    let (mut system, mut reference) = (make(), make());
    for _ in 0..8 { system.step(&[0.0; 4]).unwrap(); reference.step(&[0.0; 4]).unwrap(); }
    let old = states(&system);
    let report = system.last_frame().cloned();
    let gate = CancelGate::new(); gate.request();
    assert!(matches!(system.step_under_gate(&[0.0; 4], &gate), Err(ModalCouplingError::Cancelled)));
    for forces in [vec![], vec![f64::NAN; 4], vec![1e100; 4]] {
        assert!(system.step(&forces).is_err());
        assert_eq!(states(&system), old); assert_eq!(system.last_frame(), report.as_ref());
    }
    system.step_under_gate(&[0.0; 4], &CancelGate::new()).unwrap(); reference.step(&[0.0; 4]).unwrap();
    assert_eq!(states(&system), states(&reference)); assert_eq!(system.last_frame(), reference.last_frame());

    let mut spec = friction(0.3, 0.0); spec.maximum_force_n = 1e-12;
    let mut capped = pair(1.0, -0.001, config()).with_friction(vec![Some(spec)], &CancelGate::new()).unwrap();
    let old = states(&capped);
    for _ in 0..2 {
        assert!(matches!(capped.step(&[0.0; 4]), Err(ModalCouplingError::ContactSolve { .. })));
        assert_eq!(states(&capped), old); assert_eq!(capped.samples_rendered(), 0); assert!(capped.last_frame().is_none());
    }
}

#[test]
fn explicit_absence_and_redundant_maps_do_not_require_an_invertible_contact_matrix() {
    let mut system = build(&[1.0, 0.0], vec![contact(0, 1, -0.001), contact(0, 1, -0.001)], false, config())
        .with_friction(vec![Some(friction(0.3, 0.2)), Some(friction(0.5, 0.2))], &CancelGate::new()).unwrap();
    let frame = system.step(&[0.0; 4]).unwrap();
    assert!(frame.contacts.iter().all(|c| c.normal_force_n > 0.0));
    assert!(frame.friction.iter().all(|c| c.unwrap().reaction_n > 0.0));
    let mut absent = shared().with_friction(vec![None, Some(friction(0.3, 0.2))], &CancelGate::new()).unwrap();
    assert!(absent.friction_law(0).is_none()); assert!(absent.friction_law(1).is_some());
    assert!(absent.friction_law(2).is_none());
    let frame = absent.step(&[0.0; 6]).unwrap();
    assert!(frame.friction[0].is_none() && frame.friction[1].is_some());
}

#[test]
fn friction_admission_checks_shapes_source_model_budget_and_sample_zero() {
    let valid = friction(0.3, 0.2);
    for spec in [
        ModalFriction { coefficient: f64::NAN, ..valid.clone() },
        ModalFriction { coefficient: -0.1, ..valid.clone() },
        ModalFriction { regularization_speed_m_s: 0.0, ..valid.clone() },
        ModalFriction { maximum_force_n: 0.0, ..valid.clone() },
        ModalFriction { maximum_force_n: f64::INFINITY, ..valid.clone() },
        ModalFriction { source: " ".into(), ..valid.clone() },
        ModalFriction { left_shapes: vec![], ..valid.clone() },
        ModalFriction { right_shapes: vec![0.0, f64::NAN], ..valid.clone() },
        ModalFriction { left_shapes: vec![0.0; 2], right_shapes: vec![0.0; 2], ..valid.clone() },
    ] {
        assert!(pair(1.0, -0.001, config()).with_friction(vec![Some(spec)], &CancelGate::new()).is_err());
    }
    assert!(pair(1.0, -0.001, config()).with_friction(vec![], &CancelGate::new()).is_err());
    // n=4, p=1, k=0: combined original+extension setup bound is exactly 43.
    for cap in [42, 43] {
        let cfg = MultiContactConfig { max_setup_terms: cap, ..config() };
        assert_eq!(pair(1.0, -0.001, cfg).with_friction(vec![Some(valid.clone())], &CancelGate::new()).is_ok(), cap == 43);
    }
    let mut advanced = pair(1.0, -0.001, config()); advanced.step(&[0.0; 4]).unwrap();
    assert!(advanced.with_friction(vec![Some(valid.clone())], &CancelGate::new()).is_err());
    let attached = pair(1.0, -0.001, config()).with_friction(vec![Some(valid.clone())], &CancelGate::new()).unwrap();
    assert_eq!(attached.friction_law(0), Some(&valid));
    assert!(attached.with_friction(vec![Some(valid.clone())], &CancelGate::new()).is_err());
    let gate = CancelGate::new(); gate.request();
    assert!(matches!(pair(1.0, -0.001, config()).with_friction(vec![Some(valid)], &gate), Err(ModalCouplingError::Cancelled)));
}
