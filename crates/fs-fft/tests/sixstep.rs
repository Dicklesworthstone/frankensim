//! Six-step conformance battery (bead 27d3; runs under
//! `frontier-sixstep`). The cache-blocked six-step path is correct and
//! golden-frozen but MEASURED SLOWER than the stage walk on M4 (scalar
//! transposes vs prefetch-friendly stage streams), so it ships
//! feature-gated per the Ambition-Tag rule until the SIMD-tiled
//! transpose capsule makes it win. This battery pins: the enabled
//! dispatch predicate, cross-path value agreement, transform laws at a
//! six-step size, and the six-step golden (registered in
//! golden-couplings.json against `fs-fft:transform-bits`).

use fs_fft::{C64, Fft};
use fs_math::det;

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

#[test]
fn sixstep_dispatch_is_a_pure_function_of_n() {
    // The enabled bit contract, pinned: large even-log₂ sizes take
    // six-step; everything else stays on the stage walk.
    for n in [1usize << 16, 1 << 18, 1 << 20, 1 << 22] {
        assert!(
            Fft::takes_sixstep(n),
            "n=2^{} must dispatch six-step",
            n.ilog2()
        );
    }
    for n in [128usize, 1 << 15, 1 << 17, 1 << 21] {
        assert!(
            !Fft::takes_sixstep(n),
            "n=2^{} must stay on stages",
            n.ilog2()
        );
    }
}

#[test]
fn sixstep_agrees_with_the_stage_path() {
    // Same values, different summation order: elementwise agreement to
    // 1e-12 relative (NOT bitwise — that is the point of the per-path
    // goldens).
    let mut seed = 0x6_57E9_u64;
    for e in [16usize, 18] {
        let n = 1usize << e;
        let plan = Fft::new(n);
        let x: Vec<C64> = (0..n)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let mut six = x.clone();
        let mut scratch = vec![C64::default(); n];
        plan.forward(&mut six, &mut scratch);
        let mut staged = x.clone();
        plan.forward_via_stages(&mut staged, &mut scratch);
        let scale: f64 = staged.iter().map(|v| v.norm_sq()).sum::<f64>().sqrt();
        for (k, (a, b)) in six.iter().zip(&staged).enumerate() {
            let d = ((a.re - b.re).powi(2) + (a.im - b.im).powi(2)).sqrt();
            assert!(
                d <= 1e-12 * scale,
                "n=2^{e} bin {k}: six-step {a:?} vs stages {b:?} (d={d:.3e})"
            );
        }
        // Round-trip through the dispatched path.
        let mut back = six.clone();
        plan.inverse(&mut back, &mut scratch);
        for (k, (b, x0)) in back.iter().zip(&x).enumerate() {
            let d = ((b.re - x0.re).powi(2) + (b.im - x0.im).powi(2)).sqrt();
            assert!(d <= 1e-9, "n=2^{e} round-trip bin {k} off by {d:.3e}");
        }
    }
}

#[test]
fn sixstep_impulse_shift_and_parseval() {
    let n = 1usize << 16;
    let plan = Fft::new(n);
    let mut scratch = vec![C64::default(); n];
    // Impulse at 0 → constant spectrum.
    let mut x = vec![C64::default(); n];
    x[0] = C64::new(1.0, 0.0);
    plan.forward(&mut x, &mut scratch);
    for (k, v) in x.iter().enumerate() {
        assert!(
            (v.re - 1.0).abs() < 1e-12 && v.im.abs() < 1e-12,
            "impulse spectrum bin {k}: {v:?}"
        );
    }
    // Parseval + circular-shift theorem on a random vector.
    let mut seed = 0x9A55_u64;
    let x0: Vec<C64> = (0..n)
        .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
        .collect();
    let mut fx = x0.clone();
    plan.forward(&mut fx, &mut scratch);
    let time: f64 = x0.iter().map(|v| v.norm_sq()).sum();
    let freq: f64 = fx.iter().map(|v| v.norm_sq()).sum::<f64>() / n as f64;
    assert!(
        ((time - freq) / time).abs() < 1e-12,
        "Parseval: {time} vs {freq}"
    );
    let shift = 12_345usize;
    let mut xs: Vec<C64> = (0..n).map(|j| x0[(j + shift) % n]).collect();
    plan.forward(&mut xs, &mut scratch);
    let scale: f64 = time.sqrt();
    for k in (0..n).step_by(997) {
        // x[(j+s) mod n] ⇒ X[k]·w⁻ᵏˢ with w = e^(−2πi/n): build the
        // conjugated twiddle from the same strict kernels the plan uses.
        let ks = (k * shift) % n;
        let theta = -2.0 * std::f64::consts::PI * (ks as f64) / (n as f64);
        let w = C64::new(det::cos(theta), -det::sin(theta));
        let want = C64::new(
            fx[k].re.mul_add(w.re, -(fx[k].im * w.im)),
            fx[k].re.mul_add(w.im, fx[k].im * w.re),
        );
        let d = ((xs[k].re - want.re).powi(2) + (xs[k].im - want.im).powi(2)).sqrt();
        assert!(d < 1e-10 * scale, "shift theorem bin {k}: off {d:.3e}");
    }
}

/// JUSTIFIED FREEZE (bead 27d3, six-step slice): NEW bit territory — the
/// n=128 stage golden is untouched by construction (default dispatch is
/// pinned in-lib; the enabled dispatch is pinned above). Recorded on
/// aarch64-apple (M4 Pro), identical in debug and release; registered
/// against fs-fft:transform-bits=1 in golden-couplings.json.
const SIXSTEP_GOLDEN_HASH: u64 = 0x79aa_108f_a517_012f;

#[test]
fn sixstep_golden_hash() {
    let n = 1usize << 16;
    let plan = Fft::new(n);
    let mut scratch = vec![C64::default(); n];
    let mut seed = 0x6606_u64;
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut x: Vec<C64> = (0..n)
        .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
        .collect();
    plan.forward(&mut x, &mut scratch);
    for v in &x {
        for b in v.re.to_bits().to_le_bytes() {
            acc ^= u64::from(b);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
        for b in v.im.to_bits().to_le_bytes() {
            acc ^= u64::from(b);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    println!(
        "{{\"suite\":\"fs-fft\",\"case\":\"sixstep-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, SIXSTEP_GOLDEN_HASH,
        "six-step bits changed: {acc:#018x} vs {SIXSTEP_GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}

/// RELATIVE perf instrument (wall-clock; run explicitly in release):
/// `cargo test -p fs-fft --release --features frontier-sixstep --test sixstep -- --ignored --nocapture`
///
/// Measures the six-step path against the stage walk back-to-back and
/// INTERLEAVED (ambient load hits both sides ~equally, so the RATIO is
/// meaningful on a shared host even when absolute numbers are not).
/// Reports only — the frontier-sixstep default flip additionally
/// requires a quiet-host win per the 27d3 registry note.
#[test]
#[ignore = "wall-clock comparison lane: run explicitly in release with --ignored"]
fn sixstep_vs_stage_walk_relative_throughput() {
    let n = 1usize << 20;
    let plan = Fft::new(n);
    let mut seed = 0x27d3;
    let signal: Vec<C64> = (0..n)
        .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
        .collect();
    let mut buf = signal.clone();
    let mut scratch = vec![C64::default(); n];
    // Warm both paths once.
    buf.copy_from_slice(&signal);
    plan.forward_via_stages(&mut buf, &mut scratch);
    buf.copy_from_slice(&signal);
    plan.forward(&mut buf, &mut scratch); // dispatched = six-step here
    // Interleaved best-of-5 per side.
    let mut best_stage = f64::INFINITY;
    let mut best_six = f64::INFINITY;
    for _ in 0..5 {
        buf.copy_from_slice(&signal);
        let t0 = std::time::Instant::now();
        plan.forward_via_stages(&mut buf, &mut scratch);
        best_stage = best_stage.min(t0.elapsed().as_secs_f64());
        buf.copy_from_slice(&signal);
        let t1 = std::time::Instant::now();
        plan.forward(&mut buf, &mut scratch);
        best_six = best_six.min(t1.elapsed().as_secs_f64());
    }
    let ratio = best_stage / best_six; // > 1 means six-step is FASTER
    println!(
        "{{\"metric\":\"sixstep-vs-stage\",\"n\":{n},\"stage_s\":{best_stage:.6},\
         \"sixstep_s\":{best_six:.6},\"ratio\":{ratio:.3},\"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    // No assertion on the ratio: this is a reporting instrument. The
    // paths must still agree (the correctness batteries above gate it).
}
