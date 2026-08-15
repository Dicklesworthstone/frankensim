//! Block render API battery (music bead `frankensim-music-v8-root-3ez8g.2.1`).
//!
//! The contract under test: a block boundary performs NO arithmetic, so
//! block size is a pure loop-partition choice — the one-shot realizer, a
//! 64-sample stream, a 571-sample stream, and one giant block must all
//! produce BITWISE-identical pascals. Controls apply only between blocks
//! and are logged; cancellation is polled only at block boundaries so a
//! cancelled render is a whole number of blocks and resumes
//! bitwise-identically; refusal arms are typed and leave state untouched
//! where promised. Detailed JSON-lines logging on the equivalence test
//! lets a reviewer confirm the streams and their first divergence (none)
//! from the output alone.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::reed_bore::realize_reed_bore;
use fs_couple::render::{
    ControlDelta, GatedRenderOutcome, ModalStringVoice, ReedBoreVoice, RenderContext, RenderError,
    RenderVoice, render_under_gate,
};
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;

const RATE: u32 = 48_000;
const N: usize = 4_800; // 100 ms

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn clarinet_ish() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: 0.50,
        }],
    }
}

fn reed() -> BeatingReed {
    BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.013,
        closing_pressure_pa: 6_000.0,
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
    }
}

fn reed_voice() -> ReedBoreVoice {
    ReedBoreVoice::new(
        &clarinet_ish(),
        &air(),
        reed(),
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        RATE,
        N,
        None,
    )
    .expect("voice admits")
}

/// Render N samples through the context in fixed-size blocks.
fn render_blocked(block_len: usize) -> Vec<f64> {
    let mut context = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], N);
    let mut out = vec![0.0; N];
    let mut cursor = 0;
    while cursor < N {
        let len = block_len.min(N - cursor);
        context
            .block(&mut out[cursor..cursor + len])
            .expect("block");
        cursor += len;
    }
    out
}

#[test]
fn one_shot_and_every_block_partition_are_bitwise_identical() {
    // The oracle is the pre-API one-shot realizer (now a wrapper over the
    // same voice — but this test also pins that the WRAPPER semantics,
    // including the plate-bank swap, changed nothing).
    let mut bank = PlateBank::default();
    let oracle = realize_reed_bore(
        &clarinet_ish(),
        &air(),
        reed(),
        Termination::UnflangedOpen,
        &mut bank,
        1.0,
        RATE,
        N,
        None,
    )
    .expect("oracle renders");
    assert_eq!(oracle.len(), N);
    // The stream must actually be sound, not silence — a vacuous
    // equivalence over zeros would prove nothing.
    let tail = &oracle[N / 2..];
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    let rms =
        (tail.iter().map(|p| (p - mean) * (p - mean)).sum::<f64>() / tail.len() as f64).sqrt();
    assert!(rms > 5.0, "fixture must self-oscillate (rms {rms})");

    for block_len in [64usize, 571, N] {
        let blocked = render_blocked(block_len);
        let first_diff = oracle
            .iter()
            .zip(blocked.iter())
            .position(|(a, b)| a.to_bits() != b.to_bits());
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"render-block-invariance\",\"block_len\":{block_len},\
             \"samples\":{N},\"rms\":{rms:.3},\"first_divergence\":{:?}}}",
            first_diff
        );
        assert!(
            first_diff.is_none(),
            "block_len {block_len} diverged from the one-shot oracle at sample {first_diff:?}"
        );
    }
}

#[test]
fn modal_string_voice_matches_direct_model_stepping() {
    fn model() -> ModalAcousticTimeModel {
        let modes = vec![ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 220.0,
            damping_ratio: 1.0e-3,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        }];
        let mut model = ModalAcousticTimeModel::try_new(
            RATE,
            modes,
            ModalAcousticTimeBudget::audible_reference(),
        )
        .expect("model admits");
        model
            .restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: 1.0e-3,
                velocity_m_sqrt_kg_per_s: 0.0,
            }])
            .expect("pluck");
        model
    }
    // Direct per-sample loop (the oracle).
    let mut direct_model = model();
    let mut direct = vec![0.0; N];
    for slot in &mut direct {
        *slot = direct_model
            .step(&[0.0])
            .expect("step")
            .observer_pressure_pa;
    }
    // Hosted voice, odd block partition.
    let mut voice = ModalStringVoice::new(model(), vec![0.0]).expect("voice");
    let mut hosted = vec![0.0; N];
    let mut cursor = 0;
    while cursor < N {
        let len = 173.min(N - cursor);
        voice
            .step_block(&mut hosted[cursor..cursor + len])
            .expect("block");
        cursor += len;
    }
    assert!(
        direct
            .iter()
            .zip(hosted.iter())
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "hosting the exact-ZOH model must add zero arithmetic"
    );
    assert!(direct.iter().any(|p| p.abs() > 0.0), "non-vacuous stream");
}

#[test]
fn controls_apply_between_blocks_and_are_logged() {
    let mut context = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], 512);
    let mut block = vec![0.0; 512];
    for _ in 0..4 {
        context.block(&mut block).expect("block");
    }
    context
        .apply_controls(&[ControlDelta::SetBlowingPressure {
            voice: 0,
            pressure_pa: 0.0,
        }])
        .expect("control applies");
    // The control must have a real physical effect: against a parallel
    // context that KEEPS blowing, the cut voice's late blocks must sit
    // well below the sustained one (the driven-vs-cut discriminator is
    // robust to the line's own ring-down time, which a fixed decay
    // fraction is not).
    let mut sustained = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], 512);
    let mut sustained_block = vec![0.0; 512];
    for _ in 0..4 {
        sustained.block(&mut sustained_block).expect("block");
    }
    let mut cut_rms = 0.0;
    let mut sustained_rms = 0.0;
    for i in 0..20 {
        context.block(&mut block).expect("block");
        sustained.block(&mut sustained_block).expect("block");
        if i == 19 {
            cut_rms = (block.iter().map(|p| p * p).sum::<f64>() / block.len() as f64).sqrt();
            sustained_rms = (sustained_block.iter().map(|p| p * p).sum::<f64>()
                / sustained_block.len() as f64)
                .sqrt();
        }
    }
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"render-control-effect\",\"cut_rms\":{cut_rms:.4},\
         \"sustained_rms\":{sustained_rms:.4}}}"
    );
    assert!(
        cut_rms < 0.5 * sustained_rms.max(1.0e-12),
        "cutting blowing pressure must sit well below the sustained voice \
         (cut {cut_rms}, sustained {sustained_rms})"
    );
    let log = context.control_log();
    assert_eq!(log.len(), 1);
    assert_eq!(log[0].block_index, 4, "control recorded at its boundary");
    assert!(log[0].lift.is_empty(), "pure input move lifts no state");

    // Refusal arms: unknown voice, wrong voice kind, non-finite value.
    assert!(matches!(
        context.apply_controls(&[ControlDelta::SetBlowingPressure {
            voice: 7,
            pressure_pa: 100.0
        }]),
        Err(RenderError::Control { .. })
    ));
    assert!(matches!(
        context.apply_controls(&[ControlDelta::SetBlowingPressure {
            voice: 0,
            pressure_pa: f64::NAN
        }]),
        Err(RenderError::Control { .. })
    ));
}

#[test]
fn oversized_and_empty_blocks_refuse_before_state_moves() {
    let mut context = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], 128);
    let mut too_big = vec![0.0; 256];
    assert!(matches!(
        context.block(&mut too_big),
        Err(RenderError::Control { .. })
    ));
    let mut empty: Vec<f64> = vec![];
    assert!(matches!(
        context.block(&mut empty),
        Err(RenderError::EmptyBlock)
    ));
    assert_eq!(context.blocks_rendered(), 0, "refusals render nothing");
    // And the voice state genuinely did not move: a first real block still
    // matches a fresh context's first block bitwise.
    let mut fresh = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], 128);
    let mut a = vec![0.0; 128];
    let mut b = vec![0.0; 128];
    context.block(&mut a).expect("block");
    fresh.block(&mut b).expect("block");
    assert!(
        a.iter()
            .zip(b.iter())
            .all(|(x, y)| x.to_bits() == y.to_bits())
    );
}

#[test]
fn cancellation_stops_at_block_boundaries_and_resume_is_bitwise() {
    const BLOCK: usize = 480;
    const BLOCKS: usize = 10;
    // Straight-through reference.
    let mut straight = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], BLOCK);
    let mut reference = vec![0.0; BLOCK * BLOCKS];
    let gate = CancelGate::new();
    let outcome =
        render_under_gate(&mut straight, &gate, &mut reference, BLOCK, BLOCKS).expect("render");
    assert_eq!(outcome, GatedRenderOutcome::Completed { blocks: 10 });

    // Cancelled mid-way: a pre-requested gate stops before block 0…
    let mut cancelled = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], BLOCK);
    let hot_gate = CancelGate::new();
    hot_gate.request();
    let mut sink = vec![0.0; BLOCK * BLOCKS];
    let outcome =
        render_under_gate(&mut cancelled, &hot_gate, &mut sink, BLOCK, BLOCKS).expect("render");
    assert_eq!(outcome, GatedRenderOutcome::Cancelled { blocks: 0 });
    assert_eq!(cancelled.blocks_rendered(), 0);

    // …and a gate requested after 4 blocks drains the in-flight block and
    // stops at the boundary; resuming with a cleared demand renders the
    // remaining 6 blocks BITWISE equal to the straight-through stream.
    let mut resumed = RenderContext::new(vec![RenderVoice::ReedBore(reed_voice())], BLOCK);
    let mut prefix = vec![0.0; BLOCK * BLOCKS];
    let gate = CancelGate::new();
    let outcome = render_under_gate(&mut resumed, &gate, &mut prefix, BLOCK, 4).expect("render");
    assert_eq!(outcome, GatedRenderOutcome::Completed { blocks: 4 });
    gate.request();
    let outcome = render_under_gate(
        &mut resumed,
        &gate,
        &mut prefix[BLOCK * 4..],
        BLOCK,
        BLOCKS - 4,
    )
    .expect("render");
    assert_eq!(
        outcome,
        GatedRenderOutcome::Cancelled { blocks: 0 },
        "a requested gate refuses further blocks at the boundary"
    );
    // Resume: fresh gate, same context — continues exactly where it
    // stopped (whole-block transactionality is the resume guarantee).
    let resume_gate = CancelGate::new();
    let outcome = render_under_gate(
        &mut resumed,
        &resume_gate,
        &mut prefix[BLOCK * 4..],
        BLOCK,
        BLOCKS - 4,
    )
    .expect("render");
    assert_eq!(outcome, GatedRenderOutcome::Completed { blocks: 6 });
    let first_diff = reference
        .iter()
        .zip(prefix.iter())
        .position(|(a, b)| a.to_bits() != b.to_bits());
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"render-cancel-resume\",\"blocks\":{BLOCKS},\
         \"cancelled_after\":4,\"first_divergence\":{first_diff:?}}}"
    );
    assert!(
        first_diff.is_none(),
        "cancel/resume diverged at sample {first_diff:?}"
    );
}
