//! Production modal physics is the oracle: scheduling changes only when held
//! inputs change, never the integrator, vibration state or sample ordering.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_couple::render::{
    ControlDelta, GatedRenderOutcome, ModalStringVoice, RenderContext, RenderError, RenderVoice,
};
use fs_exec::CancelGate;

const RATE: u32 = 48_000;
const SAMPLES: usize = 512;

fn model() -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        RATE,
        [220.0, 330.0]
            .into_iter()
            .map(|hz| ModalAcousticMode {
                angular_frequency_rad_s: 2.0 * core::f64::consts::PI * hz,
                damping_ratio: 0.01,
                pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}

fn context(max_block: usize) -> RenderContext {
    RenderContext::new(
        (0..2)
            .map(|_| {
                RenderVoice::ModalString(ModalStringVoice::new(model(), vec![0.0; 2]).unwrap())
            })
            .collect(),
        max_block,
    )
}

fn event(sample: u64, voice: usize, mode: usize, force: f64) -> ScheduledControl {
    ScheduledControl {
        sample,
        delta: ControlDelta::SetModalForce {
            voice,
            mode,
            force_n_per_sqrt_kg: force,
        },
    }
}

fn score() -> Vec<ScheduledControl> {
    // Deliberately unsorted, with simultaneous updates and repeated assignments.
    vec![
        event(37, 0, 0, -0.2),
        event(0, 0, 0, 0.5),
        event(127, 0, 0, 0.0),
        event(1, 1, 1, -0.3),
        event(37, 0, 0, 0.7),
        event(37, 1, 0, 0.25),
        event(127, 1, 0, 0.0),
        event(127, 1, 1, 0.0),
        event(255, 0, 1, 0.4),
        event(300, 0, 1, 0.0),
    ]
}

fn direct(events: &[ScheduledControl]) -> Vec<f64> {
    let mut models = [model(), model()];
    let mut forces = [[0.0; 2]; 2];
    (0..SAMPLES)
        .map(|sample| {
            // Independent selection: scan input order, not the scheduler's sort.
            for entry in events.iter().filter(|entry| entry.sample == sample as u64) {
                if let ControlDelta::SetModalForce {
                    voice,
                    mode,
                    force_n_per_sqrt_kg,
                } = entry.delta
                {
                    forces[voice][mode] = force_n_per_sqrt_kg;
                }
            }
            models
                .iter_mut()
                .zip(&forces)
                .fold(0.0, |sum, (model, force)| {
                    sum + model.step(force).unwrap().observer_pressure_pa
                })
        })
        .collect()
}

fn assert_bits(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (sample, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(actual.to_bits(), expected.to_bits(), "sample {sample}");
    }
}

#[test]
fn off_grid_controls_match_direct_physics_for_every_callback_partition() {
    let events = score();
    let expected = direct(&events);
    assert!(expected[128..255].iter().any(|p| p.abs() > 1e-6), "release retains ringdown");
    let mut sorted = events.clone();
    sorted.sort_by_key(|entry| entry.sample);
    for size in [1, 7, 31, 64, 257, SAMPLES] {
        let mut renderer = ScheduledRenderer::new(context(SAMPLES), events.clone(), events.len()).unwrap();
        let mut actual = vec![0.0; SAMPLES];
        for block in actual.chunks_mut(size) {
            renderer.block(block).unwrap();
        }
        assert_bits(&actual, &expected);
        assert_eq!(renderer.samples_rendered(), SAMPLES as u64);
        assert_eq!(renderer.applied_controls(), sorted);
        assert!(renderer.pending_controls().is_empty());
        assert_eq!(renderer.context().control_log().len(), events.len());
    }
}

#[test]
fn callback_end_event_is_deferred_until_the_next_sample() {
    let events = vec![event(8, 0, 0, 0.5)];
    let mut renderer = ScheduledRenderer::new(context(16), events.clone(), 1).unwrap();
    let mut first = [99.0; 8];
    renderer.block(&mut first).unwrap();
    assert!(first.iter().all(|&sample| sample == 0.0));
    assert!(renderer.applied_controls().is_empty());
    assert_eq!(renderer.pending_controls(), events);
    let mut next = [0.0; 1];
    renderer.block(&mut next).unwrap();
    assert_bits(&next, &direct(&events)[8..9]);
    assert_eq!(renderer.applied_controls(), events);
}

#[test]
fn invalid_future_inputs_refuse_during_admission() {
    for invalid in [
        event(400, 2, 0, 1.0),
        event(400, 0, 2, 1.0),
        event(400, 0, 0, f64::NAN),
        event(400, 0, 0, f64::INFINITY),
        ScheduledControl {
            sample: 400,
            delta: ControlDelta::SetBlowingPressure { voice: 0, pressure_pa: 100.0 },
        },
        event(u64::MAX, 0, 0, 1.0),
    ] {
        assert!(ScheduledRenderer::new(context(16), vec![event(0, 0, 0, 0.5), invalid], 2).is_err());
    }
}

#[test]
fn refused_callback_does_not_consume_sample_zero_controls_or_touch_output() {
    let events = vec![event(0, 0, 0, 0.5)];
    let mut renderer = ScheduledRenderer::new(context(8), events.clone(), 1).unwrap();
    assert!(matches!(renderer.block(&mut []), Err(RenderError::EmptyBlock)));
    let mut oversized = [91.0; 9];
    assert!(matches!(renderer.block(&mut oversized), Err(RenderError::Sizing { .. })));
    assert_eq!(oversized, [91.0; 9]);
    assert_eq!(renderer.samples_rendered(), 0);
    assert!(renderer.applied_controls().is_empty());
    assert!(renderer.context().control_log().is_empty());
    let mut accepted = [0.0; 8];
    renderer.block(&mut accepted).unwrap();
    assert_bits(&accepted, &direct(&events)[..8]);
}

#[test]
fn cancellation_at_an_event_boundary_preserves_pending_controls_and_replay() {
    let events = vec![event(0, 0, 0, 0.5), event(64, 0, 0, 0.0)];
    let mut renderer = ScheduledRenderer::new(context(64), events.clone(), 2).unwrap();
    let mut actual = vec![0.0; SAMPLES];
    renderer.block(&mut actual[..64]).unwrap();
    let stop = CancelGate::new();
    stop.request();
    let mut untouched = [91.0; 64];
    assert_eq!(renderer.render_under_gate(&stop, &mut untouched, 64, 1).unwrap(),
        GatedRenderOutcome::Cancelled { blocks: 0 });
    assert_eq!(untouched, [91.0; 64]);
    assert_eq!(renderer.samples_rendered(), 64);
    assert_eq!(renderer.pending_controls(), &events[1..]);
    assert_eq!(renderer.render_under_gate(&CancelGate::new(), &mut actual[64..], 64, 7).unwrap(),
        GatedRenderOutcome::Completed { blocks: 7 });
    assert_bits(&actual, &direct(&events));
    assert_eq!(renderer.applied_controls(), events);
}

#[test]
fn gated_shape_refusal_precedes_even_an_initial_event() {
    let mut renderer = ScheduledRenderer::new(context(8), vec![event(0, 0, 0, 0.5)], 1).unwrap();
    let mut out = [91.0; 8];
    for (len, count) in [(0, 1), (9, 1), (8, 2), (8, usize::MAX)] {
        assert!(renderer.render_under_gate(&CancelGate::new(), &mut out, len, count).is_err());
        assert_eq!(out, [91.0; 8]);
        assert_eq!(renderer.samples_rendered(), 0);
        assert!(renderer.applied_controls().is_empty());
    }
}

#[test]
fn advanced_context_uses_absolute_sample_times_and_refuses_the_past() {
    let mut advanced = context(32);
    advanced.block(&mut [0.0; 20]).unwrap();
    assert_eq!(advanced.samples_rendered(), 20);
    let mut renderer = ScheduledRenderer::new(advanced, vec![event(20, 0, 0, 0.5)], 1).unwrap();
    let mut out = [0.0; 32];
    renderer.block(&mut out).unwrap();
    assert_bits(&out, &direct(&[event(20, 0, 0, 0.5)])[20..52]);
    let mut recovered = renderer.into_context();
    assert_eq!(recovered.samples_rendered(), 52);
    recovered.block(&mut [0.0; 1]).unwrap();
    assert!(matches!(ScheduledRenderer::new(recovered, vec![event(52, 0, 0, 0.5)], 1),
        Err(RenderError::Control { .. })));
}

#[test]
fn empty_schedule_and_exact_event_budget_are_admitted() {
    let mut renderer = ScheduledRenderer::new(context(8), Vec::new(), 0).unwrap();
    let mut out = [91.0; 8];
    renderer.block(&mut out).unwrap();
    assert_eq!(out, [0.0; 8]);
    assert!(ScheduledRenderer::new(context(8), vec![event(0, 0, 0, 0.5)], 1).is_ok());
    assert!(matches!(ScheduledRenderer::new(context(8), vec![event(0, 0, 0, 0.5)], 0),
        Err(RenderError::Sizing { .. })));
    assert!(matches!(ScheduledRenderer::new(context(0), Vec::new(), 0),
        Err(RenderError::Sizing { .. })));
}

#[test]
fn admission_does_not_apply_controls_to_a_context() {
    let mut context = context(8);
    context.validate_controls(&[event(0, 0, 0, 1.0).delta]).unwrap();
    assert!(context.control_log().is_empty());
    let mut out = [91.0; 8];
    context.block(&mut out).unwrap();
    assert_eq!(out, [0.0; 8]);
    assert_eq!(context.samples_rendered(), 8);
}
