//! G3 mixed pressure playback: each voice keeps its original physical stepper.
use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination,
};
use fs_couple::bowed_string::runtime::BowedStringState;
use fs_couple::bowed_string::runtime::schedule::ScheduledBowedRenderer;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::{ControlDelta, GatedRenderOutcome, ModalStringVoice, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_couple::stribeck_friction::StribeckFriction;
use fs_couple::thin_plate::CompactBody;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_scenario::RadiatingPlate;
use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};

fn config(radiation: bool, rate: u32) -> BowedRunConfig {
    BowedRunConfig {
        card: BowedStringCard {
            length_m: 0.65, tension_n: 60.0, linear_density_kg_m: 6e-4,
            bending_stiffness_n_m2: 0.0, viscous_bending_n_m2_s: 0.0,
            mode_count: 16,
            zetas: (0..16).map(|k| 1e-3 * (1.0 + 0.55 * f64::from(k))).collect(),
            sample_rate_hz: rate,
        },
        island: FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.4, 0.04).unwrap()),
        gesture: BowGesture::admit(0.45, 3.9, 0.11).unwrap(),
        steps: 1024, subsamples: 16, listener_m: 1.0,
        termination: if radiation {
            Termination::PlateOnePort {
                body: Box::new(CompactBody::from_radiator(RadiatingPlate {
                    area_m2: 3e-3, mass_kg: 0.15, frequency_hz: 280.0, damping_ratio: 0.02,
                }).unwrap()),
                ambient: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
            }
        } else { Termination::Rigid },
    }
}

fn source() -> GestureSchedule {
    let bow = |v, f, station| GestureValue::Bow { velocity_m_per_s: v, normal_force_n: f, station };
    GestureSchedule::try_new(700, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: bow(0.45, 3.9, 0.11),
        events: vec![
            GestureEvent { time_s: 0.0, transition_s: 5.0 / 700.0, value: bow(-0.3, 2.0, 0.2) },
            GestureEvent { time_s: 3.0 / 700.0, transition_s: 3.0 / 700.0, value: bow(0.25, 1.0, 0.15) },
            GestureEvent { time_s: 7.0 / 700.0, transition_s: 0.0, value: bow(0.0, 0.0, 0.15) },
            GestureEvent { time_s: 9.0 / 700.0, transition_s: 0.0, value: bow(-0.3, 2.0, 0.12) },
            GestureEvent { time_s: 12.0 / 700.0, transition_s: 0.0, value: bow(0.0, 0.0, 0.12) },
        ],
    }]).unwrap()
}

fn bow_voice(radiation: bool, rate: u32, horizon: u64, capacity: usize) -> RenderVoice {
    RenderVoice::BowedString(Box::new(ScheduledBowedRenderer::new(
        BowedStringState::new(&config(radiation, rate), capacity).unwrap(),
        &source(), "bow", horizon, 100_000, 1000,
    ).unwrap()))
}

fn modal_voice() -> RenderVoice {
    let model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * 330.0,
        damping_ratio: 0.02,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    RenderVoice::ModalString(ModalStringVoice::new(model, vec![0.01]).unwrap())
}

fn force(voice: usize, value: f64) -> ControlDelta {
    ControlDelta::SetModalForce { voice, mode: 0, force_n_per_sqrt_kg: value }
}

fn mixed(horizon: u64, capacity: usize) -> ScheduledRenderer {
    ScheduledRenderer::new(RenderContext::new(vec![
        modal_voice(), bow_voice(true, 48_000, horizon, capacity),
    ], 1024), vec![
        ScheduledControl { sample: 37, delta: force(0, 0.02) },
        ScheduledControl { sample: 511, delta: force(0, 0.0) },
    ], 2).unwrap()
}

#[test]
fn mixed_callbacks_match_independent_samplewise_pressure_and_keep_bow_ringdown() {
    let source = source();
    let mut bow = BowedStringState::new(&config(true, 48_000), 1).unwrap();
    let mut modal = RenderContext::new(vec![modal_voice()], 1);
    let mut previous = None;
    let mut expected = Vec::new();
    let mut bow_tail = Vec::new();
    for sample in 0..1024_u64 {
        let GestureValue::Bow { velocity_m_per_s, normal_force_n, station } =
            source.sample_value("bow", sample * 700 / 48_000).unwrap()
            else { panic!("typed bow value") };
        let key = [velocity_m_per_s.to_bits(), normal_force_n.to_bits(), station.to_bits()];
        if previous != Some(key) {
            bow.set_bow(velocity_m_per_s, normal_force_n, station).unwrap();
            previous = Some(key);
        }
        if sample == 37 { modal.apply_controls(&[force(0, 0.02)]).unwrap(); }
        if sample == 511 { modal.apply_controls(&[force(0, 0.0)]).unwrap(); }
        let mut pressure = [0.0];
        modal.block(&mut pressure).unwrap();
        let bow_pressure = bow.step().unwrap().radiated_pressure_pa.unwrap();
        expected.push((pressure[0] + bow_pressure).to_bits());
        if sample >= 823 { bow_tail.push(bow_pressure); }
    }
    assert!(bow_tail.iter().any(|p| p.abs() > 1e-12), "bow release retains physical ringdown");
    for partition in [1, 37, 256, 1024] {
        let mut render = mixed(1024, 1024);
        render.validate_sample_rate(48_000).unwrap();
        let mut out = vec![0.0; 1024];
        for block in out.chunks_mut(partition) { render.block(block).unwrap(); }
        assert_eq!(out.into_iter().map(f64::to_bits).collect::<Vec<_>>(), expected);
        let hosted = render.context().bowed_performance(1).unwrap();
        assert_eq!(hosted.state().samples_rendered(), 1024);
        assert_eq!(hosted.remaining_samples(), 0);
        assert_eq!(hosted.state().total_modal_energy_j().to_bits(), bow.total_modal_energy_j().to_bits());
    }
}

#[test]
fn cancellation_at_a_bow_event_preserves_both_voice_clocks_and_pending_controls() {
    let mut render = mixed(1024, 1024);
    let mut actual = vec![0.0; 1024];
    render.block(&mut actual[..69]).unwrap();
    let pending = render.context().bowed_performance(1).unwrap().pending_controls().to_vec();
    assert_eq!(pending[0].sample, 69);
    let energy = render.context().bowed_performance(1).unwrap().state().total_modal_energy_j().to_bits();
    let gate = CancelGate::new();
    gate.request();
    let mut sentinel = [12345.0; 37];
    assert_eq!(render.render_under_gate(&gate, &mut sentinel, 37, 1).unwrap(),
        GatedRenderOutcome::Cancelled { blocks: 0 });
    assert_eq!(sentinel, [12345.0; 37]);
    assert_eq!(render.samples_rendered(), 69);
    let bow = render.context().bowed_performance(1).unwrap();
    assert_eq!(bow.pending_controls(), pending);
    assert_eq!(bow.state().total_modal_energy_j().to_bits(), energy);
    for block in actual[69..].chunks_mut(37) { render.block(block).unwrap(); }
    let mut expected = vec![0.0; 1024];
    mixed(1024, 1024).block(&mut expected).unwrap();
    assert_eq!(actual, expected);
}

#[test]
fn bowed_horizon_and_capacity_refuse_before_any_voice_or_control_advances() {
    for (horizon, capacity) in [(64, 1024), (1024, 64)] {
        let mut render = mixed(horizon, capacity);
        let mut out = [12345.0; 65];
        assert!(matches!(render.block(&mut out), Err(RenderError::Sizing { .. })));
        assert_eq!(out, [12345.0; 65]);
        assert_eq!(render.samples_rendered(), 0);
        assert!(render.applied_controls().is_empty());
        assert_eq!(render.context().bowed_performance(1).unwrap().state().samples_rendered(), 0);
        let mut a = [0.0; 64];
        let mut b = [0.0; 64];
        render.block(&mut a).unwrap();
        mixed(horizon, capacity).block(&mut b).unwrap();
        assert_eq!(a, b, "shape refusal must not poison or alter either voice");
    }
}

#[test]
fn rigid_observer_wrong_control_and_mismatched_clock_are_not_silently_mixed() {
    let context = RenderContext::new(vec![modal_voice(), bow_voice(false, 48_000, 64, 64)], 64);
    assert!(ScheduledRenderer::new(context, vec![], 0).is_err());
    let mut context = RenderContext::new(vec![modal_voice(), bow_voice(true, 48_000, 64, 64)], 64);
    assert!(context.apply_controls(&[force(0, 0.5), force(1, 0.5)]).is_err());
    assert!(context.control_log().is_empty());
    assert!(context.apply_controls(&[ControlDelta::SetBlowingPressure { voice: 1, pressure_pa: 100.0 }]).is_err());
    let wrong_clock = ScheduledRenderer::new(RenderContext::new(vec![
        modal_voice(), bow_voice(true, 44_100, 64, 64),
    ], 64), vec![], 0).unwrap();
    assert!(wrong_clock.validate_sample_rate(48_000).is_err());
}
