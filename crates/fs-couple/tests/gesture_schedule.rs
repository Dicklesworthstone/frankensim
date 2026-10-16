//! Gesture-schedule e2e (music bead `frankensim-music-v8-root-3ez8g.2.3`):
//! the SAME performance driven (a) from a typed `GestureSchedule`
//! sampled at block boundaries and (b) from inline control code must be
//! BITWISE identical — schedules are replayable data, not new physics.

use fs_couple::render::ReedBoreVoice;
use fs_couple::render::schedule::{PressureGestureBinding, compile_pressure_gestures};
use fs_couple::render::schedule::{ScheduledRenderer, pressure_gesture_controls};
use fs_couple::render::{ControlDelta, RenderContext, RenderVoice};
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;
use fs_scenario::gesture::{
    GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
};

#[test]
fn signed_gesture_ramp_is_finite_through_the_public_sampler() {
    let schedule = GestureSchedule::try_new(
        4,
        vec![GestureTrack {
            id: "velocity".into(),
            target: GestureTarget::JetSpeed,
            initial: GestureValue::VelocityMPerS(-f64::MAX),
            events: vec![GestureEvent {
                time_s: 0.0,
                transition_s: 1.0,
                value: GestureValue::VelocityMPerS(f64::MAX),
            }],
        }],
    )
    .unwrap();
    let values: Vec<_> = (0..=4)
        .map(|tick| schedule.sample("velocity", tick).unwrap())
        .collect();
    assert!(values.iter().all(|value| value.is_finite()));
    assert_eq!(values[0], -f64::MAX);
    assert_eq!(values[2], 0.0);
    assert_eq!(values[4], f64::MAX);
    assert!(values.windows(2).all(|pair| pair[0] < pair[1]));
}

fn voice() -> ReedBoreVoice {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: 0.45,
        }],
    };
    let reed = BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.012,
        closing_pressure_pa: 2200.0,
        blowing_pressure_pa: 1500.0,
        attack_s: 0.01,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.35,
    };
    ReedBoreVoice::new(
        &duct,
        &gas,
        reed,
        Termination::UnflangedOpen,
        fs_couple::thin_plate::PlateBank::default(),
        1.0,
        48_000,
        4096,
        None,
    )
    .expect("voice")
}

/// The performance: pressure steps at block boundaries (the render API's
/// D17 contract — deltas between blocks; the schedule's control clock is
/// the block clock here).
fn schedule(block_rate_hz: u32) -> GestureSchedule {
    GestureSchedule::try_new(
        block_rate_hz,
        vec![GestureTrack {
            id: "blow".to_string(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(1500.0),
            events: vec![
                GestureEvent {
                    time_s: 10.0 / f64::from(block_rate_hz),
                    transition_s: 0.0,
                    value: GestureValue::PressurePa(2400.0),
                },
                GestureEvent {
                    time_s: 20.0 / f64::from(block_rate_hz),
                    transition_s: 0.0,
                    value: GestureValue::PressurePa(600.0),
                },
            ],
        }],
    )
    .expect("schedule admits")
}

#[test]
fn schedule_driven_render_is_bitwise_identical_to_inline() {
    let block_len = 480usize;
    let blocks = 30usize;
    let block_rate = 48_000 / block_len as u32; // 100 Hz control clock
    let s = schedule(block_rate);
    // Arm A: schedule-driven.
    let mut ctx_a = RenderContext::new(vec![RenderVoice::ReedBore(voice())], block_len * blocks);
    let mut out_a = vec![0.0f64; block_len * blocks];
    let mut last = f64::NAN;
    for b in 0..blocks {
        let p = s.sample("blow", b as u64).expect("sample");
        #[allow(clippy::float_cmp)] // exact change detection on a deterministic schedule
        if p != last {
            ctx_a
                .apply_controls(&[ControlDelta::SetBlowingPressure {
                    voice: 0,
                    pressure_pa: p,
                }])
                .expect("delta");
            last = p;
        }
        ctx_a
            .block(&mut out_a[b * block_len..(b + 1) * block_len])
            .expect("block");
    }
    // Arm B: the equivalent inline control code.
    let mut ctx_b = RenderContext::new(vec![RenderVoice::ReedBore(voice())], block_len * blocks);
    let mut out_b = vec![0.0f64; block_len * blocks];
    for b in 0..blocks {
        let p = if b < 10 {
            1500.0
        } else if b < 20 {
            2400.0
        } else {
            600.0
        };
        let apply = b == 0 || b == 10 || b == 20;
        if apply {
            ctx_b
                .apply_controls(&[ControlDelta::SetBlowingPressure {
                    voice: 0,
                    pressure_pa: p,
                }])
                .expect("delta");
        }
        ctx_b
            .block(&mut out_b[b * block_len..(b + 1) * block_len])
            .expect("block");
    }
    let bitwise = out_a
        .iter()
        .zip(&out_b)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    let rms = (out_a.iter().map(|x| x * x).sum::<f64>() / out_a.len() as f64).sqrt();
    assert!(
        rms > 1.0,
        "non-vacuity: the fixture must actually sound (rms {rms})"
    );
    assert!(
        bitwise,
        "schedule-driven render must be bitwise identical to inline"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"gesture-schedule-bitwise\",\"verdict\":\"pass\",\
         \"blocks\":{blocks},\"rms_pa\":{rms:.1},\"schedule_hash\":\"{}\"}}",
        s.content_hash().to_hex()
    );
}

/// G3: a physical pressure ramp on a nonintegral control/audio clock ratio
/// drives the same retained reed state regardless of the host callback size.
#[test]
fn pressure_gesture_binding_matches_samplewise_controls_across_partitions() {
    const SAMPLES: usize = 2048;
    const RATE: u32 = 700;
    let schedule = GestureSchedule::try_new(
        RATE,
        vec![GestureTrack {
            id: "pressure".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(1500.0),
            events: vec![
                GestureEvent {
                    time_s: 3.0 / f64::from(RATE),
                    transition_s: 5.0 / f64::from(RATE),
                    value: GestureValue::PressurePa(2400.0),
                },
                GestureEvent {
                    time_s: 5.0 / f64::from(RATE),
                    transition_s: 4.0 / f64::from(RATE),
                    value: GestureValue::PressurePa(600.0),
                },
            ],
        }],
    )
    .unwrap();
    let events =
        pressure_gesture_controls(&schedule, "pressure", 0, 48_000, SAMPLES as u64, 30).unwrap();
    assert_eq!(events.first().unwrap().sample, 0);
    let compiled = compile_pressure_gestures(
        &schedule,
        &[PressureGestureBinding {
            track: "pressure".into(),
            voice: 0,
        }],
        48_000,
        SAMPLES as u64,
        120,
    )
    .unwrap();
    assert_eq!(compiled, events);
    assert!(events.iter().any(|event| event.sample == 275)); // ceil(4*48000/700)
    assert!(
        pressure_gesture_controls(&schedule, "pressure", 0, 48_000, SAMPLES as u64, 29).is_err()
    );
    let mut direct = RenderContext::new(vec![RenderVoice::ReedBore(voice())], 1);
    let mut expected = vec![0.0; SAMPLES];
    for (sample, output) in expected.iter_mut().enumerate() {
        // Independent inverse clock mapping, rather than reusing the compiler's
        // event timestamps: the latest tick whose physical time has arrived.
        let tick = sample as u64 * u64::from(RATE) / 48_000;
        direct
            .apply_controls(&[ControlDelta::SetBlowingPressure {
                voice: 0,
                pressure_pa: schedule.sample("pressure", tick).unwrap(),
            }])
            .unwrap();
        direct.block(core::slice::from_mut(output)).unwrap();
    }
    assert!(expected.iter().any(|value| value.abs() > 1.0));
    for block_len in [1, 37, 480, SAMPLES] {
        let context = RenderContext::new(vec![RenderVoice::ReedBore(voice())], block_len);
        let mut renderer = ScheduledRenderer::from_pressure_gestures(
            context,
            &schedule,
            &[PressureGestureBinding {
                track: "pressure".into(),
                voice: 0,
            }],
            48_000,
            SAMPLES as u64,
            120,
            30,
        )
        .unwrap();
        let mut actual = vec![0.0; SAMPLES];
        for block in actual.chunks_mut(block_len) {
            renderer.block(block).unwrap();
        }
        assert!(
            actual
                .iter()
                .zip(&expected)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "partition {block_len}"
        );
        assert_eq!(renderer.applied_controls(), events);
    }
}

#[test]
fn pressure_gesture_binding_refuses_invalid_inputs_and_bounds_work() {
    let mut schedule = schedule(100);
    assert!(pressure_gesture_controls(&schedule, "missing", 0, 48_000, 1, 1).is_err());
    assert!(pressure_gesture_controls(&schedule, "blow", 0, 0, 1, 1).is_err());
    assert!(pressure_gesture_controls(&schedule, "blow", 0, 99, 1, 1).is_err());
    assert!(pressure_gesture_controls(&schedule, "blow", 0, 48_000, u64::MAX, 1).is_err());
    assert!(
        pressure_gesture_controls(&schedule, "blow", 0, 48_000, 0, 0)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        pressure_gesture_controls(&schedule, "blow", 0, 48_000, 1, 1)
            .unwrap()
            .len(),
        1
    );
    // The clock field is public: callers can invalidate it after admission.
    schedule.control_rate_hz = 0;
    assert!(pressure_gesture_controls(&schedule, "blow", 0, 48_000, 1, 1).is_err());
    let other = GestureSchedule::try_new(
        100,
        vec![GestureTrack {
            id: "blow".into(),
            target: GestureTarget::RestAperture,
            initial: GestureValue::LengthM(0.001),
            events: vec![],
        }],
    )
    .unwrap();
    assert!(pressure_gesture_controls(&other, "blow", 0, 48_000, 1, 1).is_err());
}

#[test]
fn pressure_performance_admits_the_actual_voice_clock_and_complete_bindings() {
    use fs_couple::render::schedule::GestureCompileError;
    let schedule = schedule(100);
    let bindings = [PressureGestureBinding {
        track: "blow".into(),
        voice: 0,
    }];
    let context = || RenderContext::new(vec![RenderVoice::ReedBore(voice())], 37);
    assert!(matches!(
        ScheduledRenderer::from_pressure_gestures(
            context(),
            &schedule,
            &bindings,
            44_100,
            480,
            100,
            10,
        ),
        Err(GestureCompileError::Render(_))
    ));
    // An empty horizon must not bypass bad destination indices.
    assert!(matches!(
        ScheduledRenderer::from_pressure_gestures(
            context(),
            &schedule,
            &[PressureGestureBinding {
                track: "blow".into(),
                voice: 1
            }],
            48_000,
            0,
            0,
            0,
        ),
        Err(GestureCompileError::Render(_))
    ));
    assert!(matches!(
        ScheduledRenderer::from_pressure_gestures(
            context(),
            &schedule,
            &bindings,
            48_000,
            480,
            100,
            0,
        ),
        Err(GestureCompileError::Render(_))
    ));
    let mut advanced = context();
    advanced.block(&mut [0.0]).unwrap();
    assert!(matches!(
        ScheduledRenderer::from_pressure_gestures(advanced, &schedule, &bindings, 48_000, 0, 0, 0,),
        Err(GestureCompileError::Invalid { .. })
    ));
    let admitted = ScheduledRenderer::from_pressure_gestures(
        context(),
        &schedule,
        &bindings,
        48_000,
        480,
        100,
        10,
    )
    .unwrap();
    assert_eq!(admitted.samples_rendered(), 0);
    assert!(admitted.applied_controls().is_empty());
    assert_eq!(admitted.pending_controls().len(), 1);
}

/// G0: exact analytical values across ramp-to-ramp and ramp-to-step interruption.
#[test]
fn interrupted_gesture_ramps_start_from_the_current_value() {
    let schedule = GestureSchedule::try_new(
        8,
        vec![GestureTrack {
            id: "pressure".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(10.0),
            events: vec![
                GestureEvent {
                    time_s: 0.0,
                    transition_s: 1.0,
                    value: GestureValue::PressurePa(50.0),
                },
                GestureEvent {
                    time_s: 0.25,
                    transition_s: 0.5,
                    value: GestureValue::PressurePa(4.0),
                },
                GestureEvent {
                    time_s: 0.5,
                    transition_s: 0.0,
                    value: GestureValue::PressurePa(30.0),
                },
                GestureEvent {
                    time_s: 0.75,
                    transition_s: 0.25,
                    value: GestureValue::PressurePa(14.0),
                },
            ],
        }],
    )
    .unwrap();
    let decoded = GestureSchedule::from_canonical_bytes(&schedule.to_canonical_bytes()).unwrap();
    for source in [&schedule, &decoded] {
        for (tick, expected) in [10.0, 15.0, 20.0, 16.0, 30.0, 30.0, 30.0, 22.0, 14.0]
            .into_iter()
            .enumerate()
        {
            assert_eq!(
                source.sample("pressure", tick as u64).unwrap(),
                expected,
                "tick {tick}"
            );
        }
    }
}
