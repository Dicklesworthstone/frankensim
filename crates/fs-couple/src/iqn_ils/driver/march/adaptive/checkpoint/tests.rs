use super::*;
use super::super::{AdaptiveSettings, AdaptiveReport};
use super::super::super::super::{BalanceControl, CouplingControls, CouplingTrial,
    InterfaceControl, StepInterval, IqnIlsConfig};
use std::cell::Cell;

const MODEL: &[u8] = b"decay-v1:rate=1;initial=1;distance=absolute/1e-4;codec=f64-le";
const CAP: usize = 4096;
fn make() -> AdaptiveEvolution<f64> {
    AdaptiveEvolution::new(1.0, vec![1.0], 0.0, vec![0.5, 1.0],
        CouplingControls { max_evaluations: 8, relaxation: 0.5,
            method: CouplingMethod::IqnIls(IqnIlsConfig::default()),
            interfaces: vec![InterfaceControl { scale: 1.0, absolute_tolerance: 1e-12,
                relative_tolerance: 0.0 }],
            balances: vec![BalanceControl { name: "decay-equation".into(), absolute_tolerance: 1e-12 }],
        }, AdaptiveSettings { initial_step_s: 0.5, minimum_step_s: 1e-5,
            maximum_step_s: 0.5, method_order: 1 }).unwrap()
}
fn go(run: &mut AdaptiveEvolution<f64>, attempts: usize) -> AdaptiveReport {
    run.advance(attempts, &mut |old: &f64, dt: StepInterval, _| {
        let next = *old / (1.0 + dt.duration_s());
        Ok::<_, String>(CouplingTrial { state: next, image: vec![next], balance_residuals: vec![0.0] })
    }, &mut |_, a, b, _| Ok((a - b).abs() / 1e-4), &mut || false).unwrap()
}
fn snapshot(run: &AdaptiveEvolution<f64>) -> Vec<u8> {
    run.checkpoint_bytes(MODEL, &run.state().to_le_bytes(), CAP).unwrap()
}
fn decode(bytes: &[u8]) -> Result<f64, String> {
    let value = f64::from_le_bytes(bytes.try_into().map_err(|_| "wrong scalar length")?);
    if !value.is_finite() { return Err("nonfinite physical state".into()); } Ok(value)
}
fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 32;
    let hash = hash_domain(DOMAIN, &bytes[..end]);
    bytes[end..].copy_from_slice(hash.as_bytes());
}

#[test]
fn rejected_attempt_checkpoint_preserves_the_shrunken_next_duration() {
    let mut original = make();
    let first = go(&mut original, 1);
    assert_eq!(first.rejected, 1);
    assert_eq!(original.time_s(), 0.0);
    assert!(original.next_step_s() < 0.5);
    let bytes = snapshot(&original);
    let mut restored = make();
    restored.restore_checkpoint(&bytes, MODEL, CAP, decode).unwrap();
    assert_eq!(restored, original);
    let continuation = go(&mut original, 4096);
    assert!(continuation.complete);
    assert_eq!(go(&mut restored, 4096), continuation);
    assert_eq!(restored.state().to_bits(), original.state().to_bits());
    assert_eq!(snapshot(&restored), snapshot(&original));
}

#[test]
fn fresh_process_style_restore_replays_every_attempt_and_mandatory_boundary() {
    let mut uninterrupted = make();
    let mut restored = make();
    let mut hit_boundary = false;
    for _ in 0..4096 {
        let expected = go(&mut uninterrupted, 1);
        assert_eq!(go(&mut restored, 1), expected);
        let bytes = snapshot(&restored);
        let mut fresh = make();
        fresh.restore_checkpoint(&bytes, MODEL, CAP, decode).unwrap();
        assert_eq!(fresh, uninterrupted);
        hit_boundary |= fresh.time_s() == 0.5;
        restored = fresh;
        if expected.complete { break; }
    }
    assert!(hit_boundary && restored.is_complete());
    assert_eq!(go(&mut restored, 1).attempts, 0);
}

#[test]
fn wrong_model_or_policy_refuses_before_decoding_and_preserves_destination() {
    let bytes = snapshot(&make());
    let called = Cell::new(false);
    let mut other = make();
    let before = other.clone();
    assert_eq!(other.restore_checkpoint(&bytes, b"different-physics", CAP, |_| {
        called.set(true); Ok(0.0)
    }), Err(CheckpointError::Model));
    assert!(!called.get()); assert_eq!(other, before);
    for which in 0..11 {
        let mut other = make();
        match which {
            0 => other.settings.initial_step_s = 0.25,
            1 => other.settings.minimum_step_s = 2e-5,
            2 => other.settings.maximum_step_s = 1.0,
            3 => other.settings.method_order = 2,
            4 => other.base.times_s[1] = 0.6,
            5 => other.base.controls.max_evaluations = 16,
            6 => other.base.controls.relaxation = 0.75,
            7 => other.base.controls.interfaces[0].absolute_tolerance = 1e-10,
            8 => other.base.controls.interfaces[0].scale = 2.0,
            9 => other.base.controls.balances[0].name = "other-equation".into(),
            _ => other.base.controls.method = CouplingMethod::RelaxedPicard,
        }
        let before = other.clone();
        assert_eq!(other.restore_checkpoint(&bytes, MODEL, CAP, |_| {
            called.set(true); Ok(0.0)
        }), Err(CheckpointError::Policy));
        assert!(!called.get()); assert_eq!(other, before);
    }
}

#[test]
fn truncation_corruption_and_trailing_data_never_reach_the_state_decoder() {
    let bytes = snapshot(&make());
    let mut bad_inputs: Vec<Vec<u8>> = (0..bytes.len()).map(|n| bytes[..n].to_vec()).collect();
    for offset in [0, 8, 40, 72, 104, 128, bytes.len()-1] {
        let mut corrupted = bytes.clone(); corrupted[offset] ^= 1; bad_inputs.push(corrupted);
    }
    let mut extra = bytes.clone(); extra.push(0); bad_inputs.push(extra);
    for bad in bad_inputs {
        let mut target = make(); let before = target.clone();
        assert!(target.restore_checkpoint(&bad, MODEL, CAP, |_| panic!("invalid bytes decoded")).is_err());
        assert_eq!(target, before);
    }
}

#[test]
fn valid_digest_cannot_hide_invalid_cursor_interface_or_payload_lengths() {
    let mut run = make();
    while run.time_s() == 0.0 { go(&mut run, 1); }
    let bytes = snapshot(&run);
    for (offset, word) in [(72, u64::MAX), (80, f64::NAN.to_bits()),
        (80, 1.0_f64.to_bits()), (88, 0_u64), (96, 0_u64), (104, u64::MAX),
        (112, f64::INFINITY.to_bits()), (120, u64::MAX)]
    {
        let mut bad = bytes.clone(); bad[offset..offset+8].copy_from_slice(&word.to_le_bytes());
        reseal(&mut bad);
        let mut target = make(); let before = target.clone();
        assert!(target.restore_checkpoint(&bad, MODEL, CAP, |_| panic!("malformed metadata decoded")).is_err());
        assert_eq!(target, before);
    }
}

#[test]
fn state_refusal_is_atomic_and_keeps_the_owner_diagnostic() {
    let bytes = snapshot(&make());
    let mut target = make(); go(&mut target, 30); let before = target.clone();
    let result = target.restore_checkpoint(&bytes, MODEL, CAP, |_| Err("wrong field mesh".into()));
    assert_eq!(result, Err(CheckpointError::State("wrong field mesh".into())));
    assert_eq!(target, before);
}

#[test]
fn limits_include_the_complete_envelope_and_model_binding() {
    let run = make(); let bytes = snapshot(&run);
    assert_eq!(bytes.len(), 168);
    assert_eq!(run.checkpoint_bytes(MODEL, &[0;8], 167), Err(CheckpointError::Limit));
    assert_eq!(run.checkpoint_bytes(&[], &[0;8], CAP), Err(CheckpointError::Model));
    assert_eq!(run.checkpoint_bytes(&vec![0; MAX_MODEL_BINDING_BYTES+1], &[0;8], CAP), Err(CheckpointError::Model));
    let mut target = make();
    assert_eq!(target.restore_checkpoint(&bytes, MODEL, bytes.len()-1, |_| panic!("over budget")),
        Err(CheckpointError::Limit));
}

#[test]
fn state_payload_preserves_signed_zero_and_subnormal_bits() {
    for bits in [(-0.0_f64).to_bits(), 1_u64] {
        let mut run = make(); run.base.state = f64::from_bits(bits);
        let bytes = snapshot(&run); let mut target = make();
        target.restore_checkpoint(&bytes, MODEL, CAP, decode).unwrap();
        assert_eq!(target.state().to_bits(), bits);
        assert_eq!(snapshot(&target), bytes);
    }
}
