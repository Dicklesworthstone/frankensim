//! Fusion #2 receipt lane (music bead `frankensim-2s4i5`).
//!
//! The aperture-junction fast mode (`ReedSolverMode::FastNewton`,
//! guarded analytic-Jacobian island Newton landed in 6a3e6afd) is a
//! DECLARED FAST MODE: not bitwise-equal to the certification-default
//! strict bisection, priced by measured budget rows, and bounded by a
//! deviation receipt. This lane emits exactly those receipts:
//!
//! 1. Junction-level root-agreement receipt — owned by the inline
//!    `fast_mode_tests::
//!    fast_newton_meets_the_strict_residual_contract_within_conditioning_band`
//!    battery (64-point stimulus sweep): every Newton root satisfies the
//!    strict path's own residual contract, and positional deviation is
//!    bounded by per-sample solver conditioning (flow tolerance over
//!    |Jacobian|) rather than a fixed number.
//!    A full RENDER is a self-oscillating feedback loop: measured this
//!    session, sample-wise |p_fast − p_strict| reaches O(kPa) within a
//!    quarter-second because tiny per-sample root differences amplify
//!    through the bore feedback — phase divergence, not solver error.
//!    This lane therefore asserts NO waveform equality.
//! 2. Fallback-hit-rate receipt — the fraction of solver-work samples
//!    handed to the strict path, printed and bounded structurally (the
//!    fast path must win the majority on nominal stimuli; a higher
//!    rate means the guard is mis-tuned and this test SHOULD fail so
//!    the number becomes a finding, not folklore).
//! 3. Before/after budget rows — median wall-clock per render for both
//!    modes on a quiet host, emitted as canonical tab-separated rows
//!    with build profile stamped (debug rows are diagnostics forever,
//!    mirroring the render_budget_lane admissibility doctrine).

use fs_couple::reed_bore::{FastSolveStats, ReedSolverMode};
use fs_couple::render::ReedBoreVoice;
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;

const RATE: u32 = 48_000;
const SAMPLES: usize = 48_000 / 4;
/// Structural bound: Newton must resolve the majority of solver-work
/// samples on the nominal fixture.
const FALLBACK_RATE_BAND: f64 = 0.5;

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn clarinet_reed() -> BeatingReed {
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

fn duct() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: 0.50,
        }],
    }
}

fn voice_in_mode(mode: ReedSolverMode) -> ReedBoreVoice {
    let mut voice = ReedBoreVoice::new(
        &duct(),
        &air20(),
        clarinet_reed(),
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        RATE,
        SAMPLES,
        None,
    )
    .expect("voice admits");
    voice.set_solver_mode(mode);
    voice
}

fn render(mode: ReedSolverMode) -> (Vec<f64>, FastSolveStats) {
    let mut voice = voice_in_mode(mode);
    let mut hist = vec![0.0; SAMPLES];
    voice.step_block(&mut hist).expect("block renders");
    let stats = voice.fast_solver_stats();
    (hist, stats)
}

fn median_of(mut times: Vec<f64>) -> f64 {
    times.sort_by(f64::total_cmp);
    times[times.len() / 2]
}

#[test]
fn fusion_receipts_deviation_fallback_and_budget() {
    let build_profile = if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    };

    // --- Render-level equivalence scope --------------------------------
    // Both modes must admit and render the identical sample count on
    // the nominal fixture. Waveform equality is NOT asserted: the loop
    // amplifies per-sample root differences (see header). Junction-level
    // root agreement lives in reed_bore's fast_mode_tests battery.
    let (_p_strict, _) = render(ReedSolverMode::Strict);
    let (p_fast, stats) = render(ReedSolverMode::FastNewton);
    assert_eq!(p_fast.len(), SAMPLES);

    // --- 2. Fallback-hit-rate receipt ---------------------------------
    let rate = stats.fallback_rate();
    println!(
        "frankensim-fusion-receipt-v1\nkind\tfallback-rate\nnewton-samples\t{}\nfallback-samples\t{}\nrate\t{:e}\nband\t{:e}",
        stats.newton_samples, stats.fallback_samples, rate, FALLBACK_RATE_BAND
    );
    assert!(
        rate <= FALLBACK_RATE_BAND,
        "fallback rate {rate} exceeds the structural {FALLBACK_RATE_BAND} \
         band: the guard is routing most samples to bisection, which is a \
         finding, not a pass"
    );

    // --- 3. Before/after budget rows ----------------------------------
    // Fresh voices per rep (state carries across blocks), quiet-host
    // caveat applies exactly as in render_budget_lane.
    let reps = 5;
    let mut strict_times = Vec::with_capacity(reps);
    let mut fast_times = Vec::with_capacity(reps);
    for _ in 0..reps {
        let mut v = voice_in_mode(ReedSolverMode::Strict);
        let mut hist = vec![0.0; SAMPLES];
        let start = std::time::Instant::now();
        v.step_block(&mut hist).expect("strict timed block");
        strict_times.push(start.elapsed().as_secs_f64());

        let mut v = voice_in_mode(ReedSolverMode::FastNewton);
        let mut hist = vec![0.0; SAMPLES];
        let start = std::time::Instant::now();
        v.step_block(&mut hist).expect("fast timed block");
        fast_times.push(start.elapsed().as_secs_f64());
    }
    let strict_median = median_of(strict_times);
    let fast_median = median_of(fast_times);
    println!(
        "frankensim-fusion-budget-row-v1\nmode\tstrict\nprofile\t{}\nsamples\t{}\nmedian-sec\t{:e}\n",
        build_profile, SAMPLES, strict_median
    );
    println!(
        "frankensim-fusion-budget-row-v1\nmode\tfast-newton\nprofile\t{}\nsamples\t{}\nmedian-sec\t{:e}\nratio-vs-strict\t{:e}\n",
        build_profile,
        SAMPLES,
        fast_median,
        strict_median / fast_median.max(f64::MIN_POSITIVE)
    );
}
