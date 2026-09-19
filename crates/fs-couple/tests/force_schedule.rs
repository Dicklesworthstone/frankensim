//! Real modal-runtime regressions for physical actuator performances.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::{
    ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice,
};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::{ControlDelta, RenderError};
use fs_math::c64::C64;

fn model(rate: u32, count: usize) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(rate, (0..count).map(|k| ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * (220.0 + 110.0 * k as f64),
        damping_ratio: 0.03,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn config(block: usize) -> ForceRenderConfig {
    ForceRenderConfig { sample_rate_hz: 48_000, max_block: block,
        max_events: 64, max_controls: 128, max_projection_terms: 4096 }
}
fn voice() -> ModalForceVoice {
    ModalForceVoice::new(model(48_000, 2), vec![vec![1.0, -0.5], vec![0.25, 2.0]],
        vec![0.5, -0.25], ForceInitialization::RetainState).unwrap()
}
fn event(sample: u64, voice: usize, port: usize, force_n: f64) -> ModalForceEvent {
    ModalForceEvent { sample, voice, port, force_n }
}
fn samples(mut render: ScheduledRenderer, block: usize, count: usize) -> Vec<f64> {
    let mut out = vec![0.0; count];
    for chunk in out.chunks_mut(block) { render.block(chunk).unwrap(); }
    out
}

#[test]
fn overlapping_ports_release_independently_and_match_direct_physics_at_every_partition() {
    // Unsorted on purpose; repeated sample-32 assignments preserve input order.
    let events = vec![event(80, 0, 1, 0.0), event(13, 0, 0, 1.0),
        event(13, 0, 1, -0.5), event(17, 0, 0, 0.0), event(31, 1, 0, 0.3),
        event(32, 1, 0, 0.8), event(32, 1, 0, 0.2)];
    let mut a = model(48_000, 2);
    let mut b = model(48_000, 1);
    let mut expected = Vec::new();
    for n in 0..256 {
        let (p0, p1) = if n < 13 { (0.5, -0.25) } else if n < 17 {
            (1.0, -0.5)
        } else if n < 80 { (0.0, -0.5) } else { (0.0, 0.0) };
        let other = if n < 31 { 0.0 } else if n == 31 { 0.3 } else { 0.2 };
        expected.push(a.step(&[p0 + 0.25*p1, -0.5*p0 + 2.0*p1]).unwrap().observer_pressure_pa
            + b.step(&[1.5*other]).unwrap().observer_pressure_pa);
    }
    assert!(expected[80..].iter().any(|x| x.abs() > 1e-7));
    for block in [1, 7, 37, 64, 256] {
        let second = ModalForceVoice::new(model(48_000, 1), vec![vec![1.5]], vec![0.0],
            ForceInitialization::RetainState).unwrap();
        let render = ScheduledRenderer::from_modal_forces(vec![voice(), second], events.clone(),
            config(256)).unwrap();
        let actual = samples(render, block, 256);
        assert!(actual.iter().zip(&expected).all(|(a,b)| a.to_bits() == b.to_bits()),
            "physical force superposition must not depend on callback size {block}");
    }
}

#[test]
fn explicit_static_preload_rings_after_release_without_a_state_reset() {
    let mut direct = model(48_000, 1);
    direct.initialize_static_equilibrium(&[2.0]).unwrap();
    let expected: Vec<_> = (0..128).map(|n| {
        direct.step(&[if n < 19 { 2.0 } else { 0.0 }]).unwrap().observer_pressure_pa
    }).collect();
    assert!(expected[..19].iter().all(|v| *v == 0.0));
    assert!(expected[19..].iter().any(|v| v.abs() > 1e-6));
    let voice = ModalForceVoice::new(model(48_000, 1), vec![vec![2.0]], vec![1.0],
        ForceInitialization::StaticPreload).unwrap();
    let actual = samples(ScheduledRenderer::from_modal_forces(vec![voice],
        vec![event(19, 0, 0, 0.0)], config(128)).unwrap(), 37, 128);
    assert_eq!(actual, expected);
}

#[test]
fn same_sample_ports_are_projected_together_and_last_assignment_wins() {
    let voice = ModalForceVoice::new(model(48_000, 1), vec![vec![1.0], vec![1.0]],
        vec![0.0, 0.0], ForceInitialization::RetainState).unwrap();
    let render = ScheduledRenderer::from_modal_forces(vec![voice], vec![
        event(3, 0, 0, 10.0), event(3, 0, 0, 1.0), event(3, 0, 1, -1.0),
    ], config(16)).unwrap();
    assert!(render.pending_controls().is_empty(), "opposing held forces cancel before projection");
    assert!(samples(render, 7, 16).iter().all(|x| *x == 0.0));
}

#[test]
fn projected_forces_and_physical_port_displacements_have_the_same_work() {
    let columns = [vec![1.0, -0.5], vec![0.25, 2.0]];
    let forces = [0.75, -0.25];
    let voice = ModalForceVoice::new(model(48_000, 2), columns.to_vec(), vec![0.0; 2],
        ForceInitialization::RetainState).unwrap();
    let render = ScheduledRenderer::from_modal_forces(vec![voice],
        vec![event(0, 0, 0, forces[0]), event(0, 0, 1, forces[1])], config(32)).unwrap();
    let mut generalized = [0.0; 2];
    for event in render.pending_controls() {
        if let ControlDelta::SetModalForce { mode, force_n_per_sqrt_kg, .. } = event.delta {
            generalized[mode] = force_n_per_sqrt_kg;
        } else { panic!("physical forces must lower only to modal inputs"); }
    }
    let mut direct = model(48_000, 2);
    let frame = direct.step(&generalized).unwrap();
    let physical_work: f64 = columns.iter().zip(forces).map(|(column, force)| {
        force * column.iter().zip(direct.states()).map(|(b,q)| b*q.displacement_m_sqrt_kg).sum::<f64>()
    }).sum();
    assert!(frame.input_work_j > 0.0);
    assert!((physical_work-frame.input_work_j).abs() < 1e-14*frame.input_work_j);
}

#[test]
fn malformed_port_maps_are_refused() {
    for (columns, initial) in [
        (vec![], vec![]), (vec![vec![1.0]], vec![0.0]),
        (vec![vec![1.0, f64::NAN]], vec![0.0]),
        (vec![vec![1.0, 2.0]], vec![]),
        (vec![vec![1.0, 2.0]], vec![f64::INFINITY]),
    ] {
        assert!(ModalForceVoice::new(model(48_000, 2), columns, initial,
            ForceInitialization::RetainState).is_err());
    }
}

#[test]
fn all_events_rates_and_projection_overflow_are_checked_before_rendering() {
    for bad in [event(1, 2, 0, 1.0), event(1, 0, 2, 1.0),
        event(u64::MAX, 0, 0, 1.0), event(1, 0, 0, f64::NAN)] {
        assert!(ScheduledRenderer::from_modal_forces(vec![voice()], vec![bad], config(32)).is_err());
    }
    let wrong_rate = ModalForceVoice::new(model(44_100, 1), vec![vec![1.0]], vec![0.0],
        ForceInitialization::RetainState).unwrap();
    assert!(ScheduledRenderer::from_modal_forces(vec![voice(), wrong_rate], vec![], config(32)).is_err());
    let overflow = ModalForceVoice::new(model(48_000, 1), vec![vec![f64::MAX]], vec![0.0],
        ForceInitialization::RetainState).unwrap();
    assert!(ScheduledRenderer::from_modal_forces(vec![overflow], vec![event(3,0,0,2.0)], config(32)).is_err());
}

#[test]
fn input_expansion_and_projection_work_budgets_have_exact_boundaries() {
    let events = vec![event(3, 0, 0, 1.0), event(3, 0, 1, 1.0)];
    let exact = ForceRenderConfig { max_events: 2, max_controls: 2, max_projection_terms: 8, ..config(32) };
    assert!(ScheduledRenderer::from_modal_forces(vec![voice()], events.clone(), exact).is_ok());
    for low in [ForceRenderConfig { max_events: 1, ..exact },
        ForceRenderConfig { max_controls: 1, ..exact },
        ForceRenderConfig { max_projection_terms: 7, ..exact }] {
        assert!(matches!(ScheduledRenderer::from_modal_forces(vec![voice()], events.clone(), low),
            Err(RenderError::Sizing { .. })));
    }
}

#[test]
fn zero_force_events_are_not_necessary_to_render_an_existing_vibration() {
    let mut image = model(48_000, 1);
    image.step(&[1.0]).unwrap();
    let mut direct = image.clone();
    let voice = ModalForceVoice::new(image, vec![vec![1.0]], vec![0.0],
        ForceInitialization::RetainState).unwrap();
    let render = ScheduledRenderer::from_modal_forces(vec![voice], vec![], ForceRenderConfig {
        max_events: 0, max_controls: 0, max_projection_terms: 1, ..config(32)
    }).unwrap();
    let expected: Vec<_> = (0..32).map(|_| direct.step(&[0.0]).unwrap().observer_pressure_pa).collect();
    assert_eq!(samples(render, 7, 32), expected);
}
