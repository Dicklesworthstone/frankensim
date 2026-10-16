//! Join independently scheduled physical voices without restarting their state.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::{ControlDelta, ModalStringVoice, RenderContext, RenderVoice};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_couple::render::schedule::force::ensemble::{
    EnsembleConfig, EnsemblePartOrigin, EnsembleRenderOutcome, EnsembleRenderer,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> EnsembleConfig {
    EnsembleConfig {
        sample_rate_hz: 48_000, max_block: 64, samples: 64,
        max_voices: 8, max_events: 64,
    }
}

fn part(rate: u32, initial: f64, prefix: usize, assignments: &[(u64, f64)]) -> ScheduledRenderer {
    let model = ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * 330.0,
        damping_ratio: 0.02, pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    let voice = ModalStringVoice::new(model, vec![initial]).unwrap();
    let context = RenderContext::new(vec![RenderVoice::ModalString(voice)], 128);
    let events = assignments.iter().map(|&(sample, force)| ScheduledControl {
        sample, delta: ControlDelta::SetModalForce {
            voice: 0, mode: 0, force_n_per_sqrt_kg: force,
        },
    }).collect();
    let mut renderer = ScheduledRenderer::new(context, events, 64).unwrap();
    if prefix != 0 { renderer.block(&mut vec![0.0; prefix]).unwrap(); }
    renderer
}

fn parts() -> Vec<ScheduledRenderer> {
    vec![
        part(48_000, 0.01, 7, &[(20, 0.02), (20, 0.03), (45, 0.0)]),
        part(48_000, -0.02, 11, &[(11, -0.03), (33, 0.02), (70, 0.0)]),
    ]
}

fn expected() -> Vec<f64> {
    let mut parts = parts();
    let mut a = vec![0.0; 64];
    let mut b = vec![0.0; 64];
    parts[0].block(&mut a).unwrap();
    parts[1].block(&mut b).unwrap();
    a.into_iter().zip(b).map(|(x, y)| x + y).collect()
}

#[test]
fn resumed_parts_keep_state_and_remap_controls_with_stable_same_sample_order() {
    let expected = expected().into_iter().map(f64::to_bits).collect::<Vec<_>>();
    for block in [1, 13, 37, 64] {
        let mut ensemble = EnsembleRenderer::from_parts(parts(), config()).unwrap();
        assert_eq!(ensemble.origins(), &[
            EnsemblePartOrigin { source_sample: 7, first_voice: 0, voices: 1 },
            EnsemblePartOrigin { source_sample: 11, first_voice: 1, voices: 1 },
        ]);
        assert_eq!(ensemble.renderer().pending_controls().iter().map(|e| e.sample).collect::<Vec<_>>(),
            vec![0, 13, 13, 22, 38, 59]);
        assert_eq!(ensemble.renderer().pending_controls()[0].delta, ControlDelta::SetModalForce {
            voice: 1, mode: 0, force_n_per_sqrt_kg: -0.03,
        });
        let mut out = vec![0.0; 64];
        assert_eq!(ensemble.render_under_gate(&CancelGate::new(), &mut out, block).unwrap(),
            EnsembleRenderOutcome::Completed { samples: 64 });
        assert_eq!(out.into_iter().map(f64::to_bits).collect::<Vec<_>>(), expected);
        assert_eq!(ensemble.samples_rendered(), 64);
        assert_eq!(ensemble.remaining_samples(), 0);
        assert_eq!(ensemble.renderer().applied_controls().len(), 6);
    }
}

#[test]
fn cancellation_on_a_pending_event_keeps_the_suffix_and_resumes_short_final_block() {
    let mut ensemble = EnsembleRenderer::from_parts(parts(), config()).unwrap();
    let mut out = vec![0.0; 64];
    ensemble.block(&mut out[..13]).unwrap();
    assert_eq!(ensemble.renderer().pending_controls()[0].sample, 13);
    let pending = ensemble.renderer().pending_controls().to_vec();
    let applied = ensemble.renderer().applied_controls().to_vec();
    let gate = CancelGate::new();
    gate.request();
    out[13..].fill(12345.0);
    assert_eq!(ensemble.render_under_gate(&gate, &mut out[13..], 37).unwrap(),
        EnsembleRenderOutcome::Cancelled { samples: 0 });
    assert_eq!(ensemble.samples_rendered(), 13);
    assert!(out[13..].iter().all(|&p| p == 12345.0));
    assert_eq!(ensemble.renderer().pending_controls(), pending);
    assert_eq!(ensemble.renderer().applied_controls(), applied);
    assert_eq!(ensemble.render_under_gate(&CancelGate::new(), &mut out[13..], 37).unwrap(),
        EnsembleRenderOutcome::Completed { samples: 51 });
    assert_eq!(out, expected());
}

#[test]
fn whole_window_capacity_and_empty_block_refusals_do_not_advance_any_part() {
    let mut ensemble = EnsembleRenderer::from_parts(parts(), config()).unwrap();
    let gate = CancelGate::new();
    let mut too_long = [12345.0; 65];
    assert!(ensemble.render_under_gate(&gate, &mut too_long, 32).is_err());
    assert_eq!(too_long, [12345.0; 65]);
    assert_eq!(ensemble.samples_rendered(), 0);
    assert!(ensemble.renderer().applied_controls().is_empty());
    assert!(ensemble.block(&mut []).is_err());
    assert!(ensemble.render_under_gate(&gate, &mut [12345.0; 1], 0).is_err());
    assert!(ensemble.render_under_gate(&gate, &mut [12345.0; 1], 65).is_err());
    assert_eq!(ensemble.render_under_gate(&gate, &mut [], 32).unwrap(),
        EnsembleRenderOutcome::Completed { samples: 0 });
    let mut out = [0.0; 64];
    ensemble.block(&mut out).unwrap();
    assert_eq!(out.as_slice(), expected());
    let mut untouched = [12345.0];
    assert!(ensemble.block(&mut untouched).is_err());
    assert_eq!(untouched, [12345.0]);
}

#[test]
fn admission_enforces_actual_clocks_source_capacities_and_explicit_budgets() {
    for config in [
        EnsembleConfig { sample_rate_hz: 0, ..config() },
        EnsembleConfig { sample_rate_hz: 44_100, ..config() },
        EnsembleConfig { max_voices: 1, ..config() },
        EnsembleConfig { max_events: 5, ..config() },
        EnsembleConfig { max_block: 129, ..config() },
        EnsembleConfig { samples: u64::MAX, ..config() },
    ] { assert!(EnsembleRenderer::from_parts(parts(), config).is_err()); }
    assert!(EnsembleRenderer::from_parts(Vec::new(), config()).is_err());
    let silent = ScheduledRenderer::new(RenderContext::new(vec![], 128), vec![], 0).unwrap();
    assert!(EnsembleRenderer::from_parts(vec![silent], config()).is_err());
    let mixed_rate = vec![part(48_000, 0.01, 0, &[]), part(44_100, 0.01, 0, &[])];
    assert!(EnsembleRenderer::from_parts(mixed_rate, config()).is_err());
}

#[test]
fn finite_windows_keep_later_controls_and_rejoining_never_replays_past_events() {
    let mut ensemble = EnsembleRenderer::from_parts(parts(), EnsembleConfig { samples: 13, ..config() }).unwrap();
    let mut out = vec![0.0; 64];
    ensemble.block(&mut out[..13]).unwrap();
    assert_eq!(ensemble.renderer().pending_controls()[0].sample, 13);
    let mut remaining = EnsembleRenderer::from_parts(vec![ensemble.into_renderer()],
        EnsembleConfig { samples: 51, ..config() }).unwrap();
    assert_eq!(remaining.origins()[0].source_sample, 13);
    assert_eq!(remaining.renderer().pending_controls()[0].sample, 0);
    remaining.render_under_gate(&CancelGate::new(), &mut out[13..], 37).unwrap();
    assert_eq!(out, expected());
}

fn bow_part(horizon: u64, prefix: usize) -> ScheduledRenderer {
    use fs_couple::bowed_string::{BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination};
    use fs_couple::bowed_string::runtime::BowedStringState;
    use fs_couple::bowed_string::runtime::schedule::ScheduledBowedRenderer;
    use fs_couple::stribeck_friction::StribeckFriction;
    use fs_couple::thin_plate::CompactBody;
    use fs_material::gas::{GasSpec, GasState};
    use fs_scenario::{RadiatingPlate, gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
    let config = BowedRunConfig {
        card: BowedStringCard {
            length_m: 0.65, tension_n: 60.0, linear_density_kg_m: 6e-4,
            bending_stiffness_n_m2: 0.0, viscous_bending_n_m2_s: 0.0,
            mode_count: 16, zetas: vec![0.01; 16], sample_rate_hz: 48_000,
        },
        island: FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.4, 0.04).unwrap()),
        gesture: BowGesture::admit(0.45, 3.9, 0.11).unwrap(),
        steps: 1024, subsamples: 16, listener_m: 1.0,
        termination: Termination::PlateOnePort {
            body: Box::new(CompactBody::from_radiator(RadiatingPlate {
                area_m2: 3e-3, mass_kg: 0.15, frequency_hz: 280.0, damping_ratio: 0.02,
            }).unwrap()),
            ambient: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
        },
    };
    let source = GestureSchedule::try_new(700, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow { velocity_m_per_s: 0.45, normal_force_n: 3.9, station: 0.11 },
        events: vec![],
    }]).unwrap();
    let mut state = BowedStringState::new(&config, 128).unwrap();
    for _ in 0..19 { state.step().unwrap(); }
    let bow = ScheduledBowedRenderer::new(state, &source, "bow", horizon, 10000, 64).unwrap();
    let mut part = ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::BowedString(Box::new(bow))], 128),
        vec![], 0).unwrap();
    if prefix != 0 { part.block(&mut vec![0.0; prefix]).unwrap(); }
    part
}

#[test]
fn joining_bow_and_modal_parts_preserves_nested_absolute_clocks_and_pressure() {
    let make = || vec![part(48_000, 0.01, 7, &[(20, 0.0)]), bow_part(128, 29)];
    let mut direct = make();
    let mut expected = vec![0.0; 64];
    let mut bow = vec![0.0; 64];
    direct[0].block(&mut expected).unwrap();
    direct[1].block(&mut bow).unwrap();
    assert!(bow.iter().any(|p| p.abs() > 1e-12));
    for (sum, b) in expected.iter_mut().zip(bow) { *sum += b; }
    let mut ensemble = EnsembleRenderer::from_parts(make(), config()).unwrap();
    assert_eq!(ensemble.renderer().context().bowed_performance(1).unwrap().state().samples_rendered(), 48);
    let mut actual = vec![0.0; 64];
    ensemble.render_under_gate(&CancelGate::new(), &mut actual, 37).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(ensemble.renderer().context().bowed_performance(1).unwrap().state().samples_rendered(), 112);
}

#[test]
fn the_entire_ensemble_window_must_fit_every_bowed_part_before_playback() {
    assert!(EnsembleRenderer::from_parts(vec![bow_part(64, 1)], config()).is_err());
    assert!(EnsembleRenderer::from_parts(vec![bow_part(128, 0)],
        EnsembleConfig { max_block: 129, ..config() }).is_err());
    let mut ensemble = EnsembleRenderer::from_parts(vec![bow_part(64, 0)], config()).unwrap();
    let mut out = [12345.0; 65];
    assert!(ensemble.render_under_gate(&CancelGate::new(), &mut out, 32).is_err());
    assert_eq!(out, [12345.0; 65]);
    assert_eq!(ensemble.renderer().context().bowed_performance(0).unwrap().state().samples_rendered(), 19);
    ensemble.render_under_gate(&CancelGate::new(), &mut out[..64], 37).unwrap();
    assert_eq!(ensemble.remaining_samples(), 0);
}
