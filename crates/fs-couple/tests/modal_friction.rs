//! G1/G3: actual two-way motion and work, not an authored oscillator signal.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment,
    ModalConnection, ModalCouplingConfig, ModalCouplingError};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,
    FrictionModalSystem, ModalFrictionConfig, ModalFrictionError, ModalFrictionPort,
    ModalFrictionRegime};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::FrictionLaw;

fn budgets() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 4, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8 }
}
fn config() -> ModalFrictionConfig {
    ModalFrictionConfig { max_iterations: 80, maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-10, sticking_tolerance_m: 1e-12 }
}
fn image(omega: f64, q: f64, v: f64) -> ModalAcousticTimeModel {
    let mut m = ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v }]).unwrap(); m
}
fn attachment(component: usize) -> ModalAttachment {
    let mass: f64 = if component == 0 { 0.002 } else { 0.01 };
    ModalAttachment { component, shapes: vec![mass.sqrt().recip()] }
}
fn network(linked: bool, velocity: f64) -> CoupledModalSystem {
    let connections = if linked { vec![ModalConnection {
        left: attachment(0), right: attachment(1), stiffness_n_m: 1e4,
        damping_n_s_m: 0.2, rest_extension_m: 0.0,
    }] } else { vec![] };
    CoupledModalSystem::new(vec![image(800.0, 0.0, velocity), image(1100.0, 0.0, 0.0)],
        connections, budgets(), &CancelGate::new()).unwrap()
}
fn law() -> FrictionLaw { FrictionLaw::Stribeck { static_mu: 0.8, kinetic_mu: 0.4,
    characteristic_speed: 0.04, viscous_per_speed: 0.0 } }
fn system(linked: bool) -> FrictionModalSystem {
    FrictionModalSystem::new(network(linked, 0.0), ModalFrictionPort {
        left: attachment(0), right: None, law: law(),
    }, config(), &CancelGate::new()).unwrap()
}
fn states(s: &FrictionModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m| m.states().iter().copied()).collect()
}
fn drive(speed_m_s: f64, normal_force_n: f64) -> FrictionDrive {
    FrictionDrive { speed_m_s, normal_force_n }
}

#[test]
fn sticking_solves_displacement_constraint_and_connected_body_reacts_back() {
    let mut linked = system(true);
    let mut isolated = system(false);
    let d = drive(0.001, 1.0);
    let actual = *linked.step(&[0.0; 2], d).unwrap();
    let other = *isolated.step(&[0.0; 2], d).unwrap();
    assert_eq!(actual.regime, ModalFrictionRegime::Sticking);
    assert_eq!(other.regime, ModalFrictionRegime::Sticking);
    assert!(actual.traction_n > other.traction_n + 1e-7, "receiver must feed back into the traction solve");
    assert!(linked.components()[1].states()[0].displacement_m_sqrt_kg > 0.0);
    let moved = attachment(0).shapes[0] * linked.components()[0].states()[0].displacement_m_sqrt_kg;
    assert!((moved-d.speed_m_s/48000.0).abs() < 1e-13);
    assert!(actual.traction_n.abs() <= 0.8);
    assert_eq!(actual.friction_dissipation_j, 0.0);
    assert!(actual.energy_residual_j.abs() <= actual.energy_tolerance_j);
}

#[test]
fn two_body_coulomb_brake_matches_independent_forced_oscillators() {
    let mut s = FrictionModalSystem::new(network(false, 0.01), ModalFrictionPort {
        left: attachment(0), right: Some(attachment(1)),
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 },
    }, config(), &CancelGate::new()).unwrap();
    let before = states(&s);
    let energy = s.total_energy_j().unwrap();
    let f = *s.step(&[0.0; 2], drive(0.0, 1.0)).unwrap();
    assert_eq!(f.regime, ModalFrictionRegime::Sliding);
    assert!((f.traction_n + 0.4).abs() < 1e-9);
    assert_eq!(f.drive_work_j, 0.0);
    let dt = 1.0/48000.0;
    for (i, omega) in [800.0_f64, 1100.0].into_iter().enumerate() {
        let force = (if i == 0 { -0.4 } else { 0.4 }) * attachment(i).shapes[0];
        let q = before[i].velocity_m_sqrt_kg_per_s * (omega*dt).sin()/omega
            +force*(1.0-(omega*dt).cos())/(omega*omega);
        let v = before[i].velocity_m_sqrt_kg_per_s * (omega*dt).cos()
            +force*(omega*dt).sin()/omega;
        let state = s.components()[i].states()[0];
        assert!((q-state.displacement_m_sqrt_kg).abs() < 1e-13);
        assert!((v-state.velocity_m_sqrt_kg_per_s).abs() < 1e-11);
    }
    assert!(f.friction_dissipation_j > 0.0);
    assert!((s.total_energy_j().unwrap()-energy+f.friction_dissipation_j).abs() < 1e-12);
}

#[test]
fn stribeck_reversals_follow_owner_law_and_close_actual_drive_work() {
    let mut s = system(true);
    let mut work = 0.0;
    let mut loss = 0.0;
    let mut positive = false;
    let mut negative = false;
    for step in 0..400 {
        let speed = if step < 200 { 0.3 } else { -0.3 };
        let f = *s.step(&[0.02, -0.03], drive(speed, 1.0)).unwrap();
        work += f.external_work_j + f.drive_work_j;
        loss += f.network_dissipation_j + f.friction_dissipation_j;
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        assert!(f.friction_dissipation_j >= 0.0);
        if f.regime == ModalFrictionRegime::Sliding {
            let v = f.slip_distance_m / s.sample_period_s();
            let expected = law().kinetic_coefficient(v.abs()).unwrap()*v.signum();
            assert!((f.traction_n-expected).abs() < 1e-9);
            assert!((f.friction_dissipation_j-f.traction_n*f.slip_distance_m).abs() < 1e-15);
            positive |= f.traction_n > 0.0;
            negative |= f.traction_n < 0.0;
        }
    }
    assert!(positive && negative);
    assert!(s.components()[1].states()[0].displacement_m_sqrt_kg.abs() > 1e-10);
    assert!((s.total_energy_j().unwrap()+loss-work).abs() < 1e-9);
}

#[test]
fn zero_normal_load_is_bitwise_the_original_coupled_ringdown() {
    let mut reference = network(true, 0.01);
    for _ in 0..19 { reference.step(&[0.02, -0.03]).unwrap(); }
    let mut base = network(true, 0.01);
    for _ in 0..19 { base.step(&[0.02, -0.03]).unwrap(); }
    let mut s = FrictionModalSystem::new(base, ModalFrictionPort {
        left: attachment(0), right: None, law: law(),
    }, config(), &CancelGate::new()).unwrap();
    assert_eq!(s.samples_rendered(), 19);
    for _ in 0..128 {
        let f = *s.step(&[0.0; 2], drive(-3.0, 0.0)).unwrap();
        let expected = reference.step(&[0.0; 2]).unwrap().observer_pressure_pa;
        assert_eq!(f.observer_pressure_pa.to_bits(), expected.to_bits());
        assert_eq!(f.regime, ModalFrictionRegime::Released);
        assert_eq!(f.traction_n, 0.0);
        assert_eq!(f.drive_work_j, 0.0);
        for (a, b) in s.components().iter().zip(reference.components()) { assert_eq!(a.states(), b.states()); }
    }
}

#[test]
fn cancellation_and_invalid_inputs_preserve_state_frame_and_exact_retry() {
    let mut s = system(true);
    let mut direct = system(true);
    let d = drive(0.3, 1.0);
    for _ in 0..19 { s.step(&[0.0; 2], d).unwrap(); direct.step(&[0.0; 2], d).unwrap(); }
    let before = states(&s);
    let frame = s.last_frame().copied();
    let gate = CancelGate::new(); gate.request();
    assert!(matches!(s.step_under_gate(&[0.0; 2], d, &gate),
        Err(ModalFrictionError::Network(ModalCouplingError::Cancelled))));
    for bad in [drive(f64::NAN, 1.0), drive(0.3, -1.0), drive(0.3, f64::INFINITY)] {
        assert!(s.step(&[0.0; 2], bad).is_err());
    }
    for external in [vec![], vec![f64::NAN, 0.0], vec![1e100, 0.0]] {
        assert!(s.step(&external, d).is_err());
    }
    assert_eq!(states(&s), before);
    assert_eq!(s.last_frame().copied(), frame);
    s.step_under_gate(&[0.0; 2], d, &CancelGate::new()).unwrap();
    direct.step(&[0.0; 2], d).unwrap();
    assert_eq!(states(&s), states(&direct));
    assert_eq!(s.last_frame(), direct.last_frame());
}

#[test]
fn root_force_and_actual_slip_budgets_refuse_transactionally() {
    let mut one = config(); one.max_iterations = 1;
    let mut force = config(); force.maximum_force_n = 0.01;
    let mut slip = config(); slip.maximum_slip_speed_m_s = 0.001;
    for cfg in [one, force, slip] {
        let mut s = FrictionModalSystem::new(network(true, 0.0), ModalFrictionPort {
            left: attachment(0), right: None, law: law(),
        }, cfg, &CancelGate::new()).unwrap();
        let before = states(&s);
        assert!(s.step(&[0.0; 2], drive(0.3, 1.0)).is_err());
        assert_eq!(states(&s), before);
        assert_eq!(s.samples_rendered(), 0);
        assert!(s.last_frame().is_none());
        s.step(&[0.0; 2], drive(0.0, 0.0)).unwrap();
    }
}

#[test]
fn zero_coefficients_and_velocity_strengthening_remain_physical() {
    for law in [FrictionLaw::Coulomb { static_mu: 0.0, kinetic_mu: 0.0 },
        FrictionLaw::VelocityDependent { static_mu: 0.8, mu_zero: 0.4, slope_per_speed: 2.0 },
        FrictionLaw::Stribeck { static_mu: 0.8, kinetic_mu: 0.4,
            characteristic_speed: 0.04, viscous_per_speed: 2.0 }] {
        let mut s = FrictionModalSystem::new(network(true, 0.0), ModalFrictionPort {
            left: attachment(0), right: None, law,
        }, config(), &CancelGate::new()).unwrap();
        let f = *s.step(&[0.0; 2], drive(0.3, 1.0)).unwrap();
        assert_eq!(f.regime, ModalFrictionRegime::Sliding);
        assert!(f.traction_n >= 0.0 && f.friction_dissipation_j >= 0.0);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
    }
}

#[test]
fn malformed_law_or_attachment_is_refused_before_any_sample() {
    let mut bad = law();
    if let FrictionLaw::Stribeck { characteristic_speed, .. } = &mut bad { *characteristic_speed = 0.0; }
    for port in [ModalFrictionPort { left: attachment(0), right: None, law: bad },
        ModalFrictionPort { left: ModalAttachment { component: 3, shapes: vec![1.0] }, right: None, law: law() },
        ModalFrictionPort { left: attachment(0), right: Some(attachment(0)), law: law() }] {
        assert!(FrictionModalSystem::new(network(false, 0.0), port, config(), &CancelGate::new()).is_err());
    }
}

use fs_couple::render::{ControlDelta, GatedRenderOutcome, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::coupled::render::friction::{FrictionGestureConfig, FrictionModalVoice};
use fs_couple::render::schedule::force::ensemble::{EnsembleConfig, EnsembleRenderer};
use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};

fn bow_source() -> GestureSchedule {
    let value = |v, n| GestureValue::Bow { velocity_m_per_s: v, normal_force_n: n, station: 0.11 };
    GestureSchedule::try_new(700, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 }, initial: value(0.3, 1.0),
        events: vec![
            GestureEvent { time_s: 0.0, transition_s: 5.0/700.0, value: value(-0.3, 0.5) },
            GestureEvent { time_s: 3.0/700.0, transition_s: 3.0/700.0, value: value(0.25, 1.0) },
            GestureEvent { time_s: 7.0/700.0, transition_s: 0.0, value: value(0.0, 0.0) },
            GestureEvent { time_s: 9.0/700.0, transition_s: 0.0, value: value(-0.3, 1.0) },
            GestureEvent { time_s: 12.0/700.0, transition_s: 0.0, value: value(0.0, 0.0) },
        ],
    }]).unwrap()
}
fn gesture_config() -> FrictionGestureConfig {
    FrictionGestureConfig { sample_rate_hz: 48000, samples: 1024, max_block: 1024,
        max_work: 10000, max_events: 128, station_fraction: 0.11 }
}
fn performance() -> ScheduledRenderer {
    ScheduledRenderer::from_friction_gesture(system(true), vec![0.02, -0.03],
        &bow_source(), "bow", gesture_config()).unwrap()
}

#[test]
fn typed_bow_performance_matches_direct_two_way_physics_at_every_sample() {
    let source = bow_source();
    let decoded = GestureSchedule::from_canonical_bytes(&source.to_canonical_bytes()).unwrap();
    let mut direct = system(true);
    for _ in 0..19 { direct.step(&[0.02, -0.03], drive(0.3, 1.0)).unwrap(); }
    let initial_energy = direct.total_energy_j().unwrap().to_bits();
    let mut expected = Vec::new();
    for sample in 0..1024_u64 {
        let GestureValue::Bow { velocity_m_per_s, normal_force_n, .. } =
            source.sample_value("bow", sample*700/48000).unwrap() else { panic!("bow") };
        expected.push(direct.step(&[0.02, -0.03], drive(velocity_m_per_s, normal_force_n))
            .unwrap().observer_pressure_pa.to_bits());
    }
    for src in [&source, &decoded] {
        for partition in [1, 37, 256, 1024] {
            let mut physical = system(true);
            for _ in 0..19 { physical.step(&[0.02, -0.03], drive(0.3, 1.0)).unwrap(); }
            let mut render = ScheduledRenderer::from_friction_gesture(physical, vec![0.02, -0.03],
                src, "bow", gesture_config()).unwrap();
            let retained = render.context().friction_voice(0).unwrap().system();
            assert_eq!(retained.samples_rendered(), 19);
            assert_eq!(retained.total_energy_j().unwrap().to_bits(), initial_energy);
            let mut out = [0.0; 1024];
            for block in out.chunks_mut(partition) { render.block(block).unwrap(); }
            assert_eq!(out.map(f64::to_bits).as_slice(), expected);
            let retained = render.context().friction_voice(0).unwrap().system();
            assert_eq!(states(retained), states(&direct));
            assert_eq!(retained.samples_rendered(), 1043);
            assert_eq!(retained.last_frame(), direct.last_frame());
        }
    }
}

#[test]
fn shared_cancellation_preserves_unconsumed_friction_controls_and_resumes_exactly() {
    let mut render = performance();
    let mut actual = [0.0; 1024];
    render.block(&mut actual[..69]).unwrap();
    assert_eq!(render.pending_controls()[0].sample, 69);
    let pending = render.pending_controls().to_vec();
    let state = states(render.context().friction_voice(0).unwrap().system());
    let gate = CancelGate::new(); gate.request();
    let mut sentinel = [12345.0; 37];
    assert_eq!(render.render_under_gate(&gate, &mut sentinel, 37, 1).unwrap(),
        GatedRenderOutcome::Cancelled { blocks: 0 });
    assert_eq!(sentinel, [12345.0; 37]);
    assert_eq!(render.pending_controls(), pending);
    assert_eq!(states(render.context().friction_voice(0).unwrap().system()), state);
    for block in actual[69..].chunks_mut(37) { render.block(block).unwrap(); }
    let mut expected = [0.0; 1024]; performance().block(&mut expected).unwrap();
    assert_eq!(actual, expected);
}

#[test]
fn ensemble_join_reindexes_friction_controls_without_restarting_either_network() {
    let part = |prefix: usize| {
        let mut renderer = performance(); renderer.block(&mut vec![0.0; prefix]).unwrap(); renderer
    };
    let mut a = part(37); let mut b = part(69);
    let mut left = [0.0; 128]; let mut right = [0.0; 128];
    a.block(&mut left).unwrap(); b.block(&mut right).unwrap();
    let expected: Vec<u64> = left.iter().zip(right).map(|(a,b)| (a+b).to_bits()).collect();
    let mut ensemble = EnsembleRenderer::from_parts(vec![part(37), part(69)], EnsembleConfig {
        sample_rate_hz: 48000, max_block: 1024, samples: 128, max_voices: 2, max_events: 256,
    }).unwrap();
    let mut out = [0.0; 128];
    for block in out.chunks_mut(23) { ensemble.block(block).unwrap(); }
    assert_eq!(out.map(f64::to_bits).as_slice(), expected);
    assert_eq!(ensemble.renderer().context().friction_voice(0).unwrap().system().samples_rendered(), 165);
    assert_eq!(ensemble.renderer().context().friction_voice(1).unwrap().system().samples_rendered(), 197);
}

#[test]
fn gesture_binding_refuses_wrong_clocks_moving_ports_and_exhausted_budgets() {
    let source = bow_source();
    let base = gesture_config();
    for cfg in [FrictionGestureConfig { sample_rate_hz: 44100, ..base },
        FrictionGestureConfig { max_events: 0, ..base }, FrictionGestureConfig { max_work: 1, ..base },
        FrictionGestureConfig { station_fraction: 0.2, ..base }, FrictionGestureConfig { max_block: 0, ..base }] {
        assert!(ScheduledRenderer::from_friction_gesture(system(true), vec![0.0; 2], &source, "bow", cfg).is_err());
    }
    assert!(ScheduledRenderer::from_friction_gesture(system(true), vec![0.0; 2], &source, "missing", base).is_err());
    let mut track = source.tracks()[0].clone();
    track.events.push(GestureEvent { time_s: 1.0, transition_s: 0.0,
        value: GestureValue::Bow { velocity_m_per_s: 0.3, normal_force_n: 1.0, station: 0.2 } });
    let moving = GestureSchedule::try_new(700, vec![track]).unwrap();
    assert!(ScheduledRenderer::from_friction_gesture(system(true), vec![0.0; 2], &moving, "bow", base).is_err());
}

#[test]
fn mixed_control_batch_admission_cannot_partially_change_friction_inputs() {
    let make = || RenderContext::new(vec![RenderVoice::FrictionModal(Box::new(
        FrictionModalVoice::new(system(true), vec![0.02, -0.03], drive(0.3, 1.0)).unwrap(),
    ))], 128);
    let mut context = make();
    assert!(context.apply_controls(&[
        ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 100.0 },
        ControlDelta::SetFrictionDrive { voice: 0, speed_m_s: 0.3, normal_force_n: -1.0 },
    ]).is_err());
    assert!(context.control_log().is_empty());
    let mut actual = [0.0; 128]; let mut expected = [0.0; 128];
    context.block(&mut actual).unwrap(); make().block(&mut expected).unwrap();
    assert_eq!(actual, expected);
}
