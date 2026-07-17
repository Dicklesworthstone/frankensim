//! BMI2 Morton capsule: PDEP/PEXT bit interleave on x86-64.
//!
//! One instruction per axis instead of a six-step magic-bits cascade; on
//! Zen 3+ and Intel parts PDEP/PEXT are 3-cycle ops (pre-Zen 3 AMD
//! microcode is slow — the dispatch layer only routes here when the CPU
//! advertises BMI2, and correctness never depends on speed). The G0
//! equivalence battery pins this capsule bit-for-bit to the magic-bits
//! reference twin in `super`.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use core::arch::x86_64::{_pdep_u64, _pext_u64};

const MASK_X: u64 = 0x1249_2492_4924_9249;
const MASK_Y: u64 = MASK_X << 1;
const MASK_Z: u64 = MASK_X << 2;

/// Interleave via PDEP. Called through the dispatch table only after
/// `is_x86_feature_detected!("bmi2")` returned true.
pub(super) fn encode(x: u32, y: u32, z: u32) -> u64 {
    // SAFETY: the dispatch layer (super::fns) selects this function only
    // when BMI2 is detected at runtime, so the target-feature contract of
    // `encode_bmi2` is met. Inputs are plain integers; no memory access.
    unsafe { encode_bmi2(x, y, z) }
}

/// De-interleave via PEXT. Same dispatch contract as [`encode`].
pub(super) fn decode(code: u64) -> (u32, u32, u32) {
    // SAFETY: as in `encode` — BMI2 presence established by the dispatcher.
    unsafe { decode_bmi2(code) }
}

#[target_feature(enable = "bmi2")]
unsafe fn encode_bmi2(x: u32, y: u32, z: u32) -> u64 {
    // `_pdep_u64` is register-only and safe INSIDE this target_feature
    // context (target-feature 1.1); the unsafety lives at the call sites
    // above, which the dispatcher guards with runtime detection.
    _pdep_u64(u64::from(x), MASK_X)
        | _pdep_u64(u64::from(y), MASK_Y)
        | _pdep_u64(u64::from(z), MASK_Z)
}

#[target_feature(enable = "bmi2")]
unsafe fn decode_bmi2(code: u64) -> (u32, u32, u32) {
    // Register-only PEXT; same contract as `encode_bmi2`.
    (
        _pext_u64(code, MASK_X) as u32,
        _pext_u64(code, MASK_Y) as u32,
        _pext_u64(code, MASK_Z) as u32,
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn capsule_matches_the_magic_reference_when_available() {
        if !std::arch::is_x86_feature_detected!("bmi2") {
            let mut emitter =
                fs_obs::Emitter::new("fs-substrate/morton-bmi2", "capsule-equivalence/outcome");
            let event = emitter.emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::CapabilityDomainDecision {
                    capability: "morton-bmi2".to_string(),
                    domain: "current-host".to_string(),
                    decision: fs_obs::CapabilityDecision::Refused,
                    detail:
                        "CPU does not advertise BMI2; capsule equivalence execution was skipped"
                            .to_string(),
                },
                None,
            );
            fs_obs::lint_failure_record(&event)
                .expect("BMI2 capability refusal must be replayable");
            let line = event.to_jsonl();
            fs_obs::validate_line(&line)
                .expect("BMI2 capability refusal must use the fs-obs wire schema");
            eprintln!("{line}");
            return;
        }
        for x in [0u32, 1, 7, 63, 1 << 20, (1 << 21) - 1] {
            for y in [0u32, 2, 30, (1 << 21) - 1] {
                for z in [0u32, 5, 511, (1 << 21) - 1] {
                    let want = crate::morton::encode_magic(x, y, z);
                    assert_eq!(super::encode(x, y, z), want);
                    assert_eq!(super::decode(want), (x, y, z));
                }
            }
        }
    }
}
