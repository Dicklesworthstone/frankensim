//! Cross-ISA determinism goldens for the `fs-couple` exact-ZOH modal
//! acoustic render (bead `frankensim-music-v8-root-3ez8g.13.4`).
//!
//! The exact zero-order-hold modal runtime is the crate's render kernel:
//! instrument fixtures become GOLDEN AUDIO artifacts, so its full
//! trajectory must be bit-identical on both reference ISA families in
//! both build modes. The digest below folds every state (q, v) of every
//! mode at every sample of a pinned three-mode pluck — the same fixture
//! family as the committed bakeoff receipt, without receipt machinery.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fold(acc: u64, v: f64) -> u64 {
    v.to_bits()
        .to_le_bytes()
        .iter()
        .fold(acc, |a, &b| (a ^ u64::from(b)).wrapping_mul(FNV_PRIME))
}

/// Shared string card: a light steel-ish string, three retained modes
/// (mirrors the bakeoff fixture so a drift localizes to the runtime).
const LENGTH_M: f64 = 0.65;
const TENSION_N: f64 = 60.0;
const LIN_DENSITY_KG_M: f64 = 6.0e-4;
const MODES: usize = 3;
const ZETAS: [f64; MODES] = [1.0e-3, 1.5e-3, 2.0e-3];
const Q0: f64 = 1.0e-3;
const SAMPLE_RATE_HZ: u32 = 48_000;
const STEPS: usize = 2_400;

/// Verified bit-identical aarch64-apple and x86_64-linux (debug) on
/// 2026-08-23, bead frankensim-music-v8-root-3ez8g.13.4.
const GOLDEN_HASH: u64 = 0x4323_20b3_c2bf_06d9;

#[test]
fn exact_zoh_render_trajectory_is_cross_isa_golden() {
    let wave_speed = (TENSION_N / LIN_DENSITY_KG_M).sqrt();
    let omegas: [f64; MODES] =
        core::array::from_fn(|k| (k + 1) as f64 * core::f64::consts::PI * wave_speed / LENGTH_M);
    let modes = omegas
        .iter()
        .zip(ZETAS)
        .map(|(&omega, zeta)| ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        })
        .collect::<Vec<_>>();
    let mut model = ModalAcousticTimeModel::try_new(
        SAMPLE_RATE_HZ,
        modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .expect("string modes admit under the Nyquist guard");
    let plucked = (0..MODES)
        .map(|_| ModalAcousticState {
            displacement_m_sqrt_kg: Q0,
            velocity_m_sqrt_kg_per_s: 0.0,
        })
        .collect::<Vec<_>>();
    model
        .restore_states(&plucked)
        .expect("finite pluck states restore");

    let zero_force = vec![0.0; MODES];
    let mut acc = FNV_OFFSET;
    for _ in 0..STEPS {
        model.step(&zero_force).expect("free decay step");
        for state in model.states() {
            acc = fold(acc, state.displacement_m_sqrt_kg);
            acc = fold(acc, state.velocity_m_sqrt_kg_per_s);
        }
    }

    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"cross-isa-zoh-render\",\"arch\":\"{}\",\
         \"profile\":\"{}\",\"digest\":\"{acc:#018x}\",\"verdict\":\"golden-check\"}}",
        std::env::consts::ARCH,
        if cfg!(debug_assertions) {
            "debug"
        } else {
            "release"
        },
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "ZOH render trajectory bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — \
         cross-ISA golden event: bisect stage-wise, name the hazard, route through \
         det:: in the same commit"
    );
}
