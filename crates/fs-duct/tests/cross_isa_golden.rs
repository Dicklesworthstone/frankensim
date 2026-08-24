//! Cross-ISA determinism goldens for the `fs-duct` TMM sweep core (bead
//! `frankensim-music-v8-root-3ez8g.13.4`, G5 audit).
//!
//! The TMM evaluation path (`cis`/`exp_i`, ZK losses, terminations) is
//! built exclusively from `fs_math::det` IEEE-754 arithmetic, so the sweep
//! digests below must be bit-identical on the two reference ISA families
//! (aarch64-apple, x86-64-linux) in BOTH build modes. A mismatch is a
//! golden event under the audit protocol: bisect stage-wise, name the
//! platform-libm hazard, route it through `det::` in the same commit —
//! never loosen or re-pin silently.

use fs_duct::{Duct, HoleState, LossModel, Segment, Termination, impedance_peaks, impedance_sweep};
use fs_material::gas::{GasSpec, GasState};

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fold(acc: u64, v: f64) -> u64 {
    v.to_bits()
        .to_le_bytes()
        .iter()
        .fold(acc, |a, &b| (a ^ u64::from(b)).wrapping_mul(FNV_PRIME))
}

fn fold_u64(acc: u64, v: u64) -> u64 {
    (acc ^ v).wrapping_mul(FNV_PRIME)
}

/// Authored four-hole cylinder (not corpus-derived: this golden pins the
/// kernel arithmetic, not the corpus registration seam).
fn four_hole(states: [HoleState; 4]) -> Duct {
    let bore = 6.0e-3;
    let total_length = 0.6;
    let positions = [0.12, 0.24, 0.36, 0.48];
    let mut segments = Vec::new();
    let mut cursor = 0.0;
    for (index, state) in states.iter().enumerate() {
        segments.push(Segment::Cylinder {
            radius: bore,
            length: positions[index] - cursor,
        });
        segments.push(Segment::ToneHole {
            hole_radius: 3.0e-3,
            chimney_height: 4.0e-3,
            bore_radius: bore,
            state: *state,
        });
        cursor = positions[index];
    }
    segments.push(Segment::Cylinder {
        radius: bore,
        length: total_length - cursor,
    });
    Duct { segments }
}

fn sweep_digest(duct: &Duct, air: &GasState, peak_count_out: &mut usize) -> u64 {
    let sweep = impedance_sweep(
        duct,
        air,
        2.0 * core::f64::consts::PI * 150.0,
        2.0 * core::f64::consts::PI * 1000.0,
        12_000,
        LossModel::WideTube,
        Termination::UnflangedOpen,
    )
    .expect("sweep");
    let mut acc = FNV_OFFSET;
    for response in &sweep {
        acc = fold(acc, response.omega);
        acc = fold(acc, response.impedance.re);
        acc = fold(acc, response.impedance.im);
        acc = fold(acc, response.min_shear_number);
        acc = fold(acc, response.mouth_ka);
    }
    let peaks = impedance_peaks(&sweep);
    *peak_count_out = peaks.len();
    for peak in peaks {
        acc = fold(acc, sweep[peak].omega);
    }
    acc
}

/// Combined golden over two fingerings' full sweeps plus their peak
/// counts (counts mixed as integers — no float round trip).
/// Verified bit-identical aarch64-apple and x86_64-linux (debug) on
/// 2026-08-23, bead frankensim-music-v8-root-3ez8g.13.4.
const GOLDEN_HASH: u64 = 0xd6e9_724c_5414_cf8d;

#[test]
fn tmm_sweep_digest_is_cross_isa_golden() {
    use HoleState::{Closed as X, Open as O};
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let mut closed_peaks = 0usize;
    let mut open_peaks = 0usize;
    let all_closed = sweep_digest(&four_hole([X, X, X, X]), &air, &mut closed_peaks);
    let three_closed = sweep_digest(&four_hole([X, X, X, O]), &air, &mut open_peaks);
    // The two fingerings must stay acoustically DISTINCT inside the golden:
    // identical digests across different geometry would mean the sweep is
    // not actually exercising the tone-hole chain.
    assert!(
        all_closed != three_closed && closed_peaks > 0 && open_peaks > 0,
        "fingerings must produce distinct sweeps with visible peaks"
    );
    let acc = fold_u64(fold_u64(all_closed, three_closed), {
        let c = u64::try_from(closed_peaks).expect("peak count");
        let o = u64::try_from(open_peaks).expect("peak count");
        (c << 32) | o
    });
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"cross-isa-tmm\",\"arch\":\"{}\",\"profile\":\"{}\",\
         \"digest\":\"{acc:#018x}\",\"verdict\":\"golden-check\"}}",
        std::env::consts::ARCH,
        if cfg!(debug_assertions) { "debug" } else { "release" },
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "TMM sweep bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — cross-ISA golden \
         event: bisect stage-wise, name the hazard, route through det:: in the same commit"
    );
}
