//! Gesture-schedule e2e (music bead `frankensim-music-v8-root-3ez8g.2.3`):
//! the SAME performance driven (a) from a typed `GestureSchedule`
//! sampled at block boundaries and (b) from inline control code must be
//! BITWISE identical — schedules are replayable data, not new physics.

use fs_couple::render::ReedBoreVoice;
use fs_couple::render::{ControlDelta, RenderContext, RenderVoice};
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;
use fs_scenario::gesture::{
    GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
};

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
