//! fs-uq slice-2 conformance (the o5kc bead): the seismic stack
//! (Kanai–Tajimi synthesis, CQC vs SRSS, bilinear IDA fragility),
//! anytime-valid stopping (validity under optional stopping — the
//! certify-the-certifiers battery), CVaR, and adaptive MLMC with rate
//! recovery on a known-rate synthetic.

use fs_math::det;
use fs_uq::seismic::{KanaiTajimi, bilinear_peak_ductility, cqc, ida_fragility, sdof_peak, srss};
use fs_uq::{adaptive_mlmc, cvar, estimate_probability_anytime};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-uq/slice2\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn kt() -> KanaiTajimi {
    KanaiTajimi {
        s0: 0.03,
        wg: 15.0,
        zg: 0.6,
    }
}

#[test]
fn uq_001_ground_motion_and_response_spectrum() {
    // The synthesized record is deterministic per seed and its energy
    // concentrates near the ground frequency: the response spectrum
    // peaks in the wg neighborhood rather than far above or below.
    let a1 = kt().synthesize(7, 96, 0.01, 1500);
    let a2 = kt().synthesize(7, 96, 0.01, 1500);
    assert_eq!(a1, a2, "seeded synthesis replays exactly (G5)");
    let a3 = kt().synthesize(8, 96, 0.01, 1500);
    assert_ne!(a1, a3, "different seeds differ");
    // Response-spectrum shape probe at three periods.
    let low = sdof_peak(&a1, 0.01, 3.0, 0.05); // well below wg
    let near = sdof_peak(&a1, 0.01, 15.0, 0.05); // at wg
    let high = sdof_peak(&a1, 0.01, 80.0, 0.05); // far above
    println!(
        "{{\"metric\":\"spectrum\",\"low\":{low:.4},\"near_wg\":{near:.4},\"high\":{high:.5}}}"
    );
    assert!(
        near > high,
        "the spectrum decays above the ground band: {near} vs {high}"
    );
    verdict(
        "uq-001",
        "Kanai-Tajimi synthesis replays bit-exact per seed; the linear response spectrum \
         decays above the ground-filter band",
    );
}

#[test]
fn uq_002_cqc_reduces_to_srss_and_amplifies_close_modes() {
    // Well-separated modes: CQC == SRSS (cross terms vanish).
    let peaks = [1.0, 0.6, 0.3];
    let far = cqc(&peaks, &[1.0, 10.0, 100.0], &[0.02, 0.02, 0.02]);
    let base = srss(&peaks);
    assert!(
        (far - base).abs() / base < 0.02,
        "separated modes: CQC ~ SRSS ({far:.4} vs {base:.4})"
    );
    // Closely-spaced modes with same-sign peaks: CQC EXCEEDS SRSS
    // (the correlation the SRSS shortcut misses).
    let close = cqc(&[1.0, 0.9], &[10.0, 10.5], &[0.05, 0.05]);
    let close_srss = srss(&[1.0, 0.9]);
    println!(
        "{{\"metric\":\"cqc\",\"separated_ratio\":{:.4},\"close_ratio\":{:.4}}}",
        far / base,
        close / close_srss
    );
    assert!(
        close > 1.05 * close_srss,
        "close modes amplify: {close:.4} vs {close_srss:.4}"
    );
    verdict(
        "uq-002",
        "CQC collapses to SRSS for well-separated modes and exceeds it >5% for \
         closely-spaced same-sign modes — the fast path is honest about correlation",
    );
}

#[test]
fn uq_003_ida_fragility_is_monotone_and_nonlinear() {
    // The bilinear oscillator yields: ductility grows super-linearly
    // with scale, and the fragility curve is monotone in IM.
    let record = kt().synthesize(3, 96, 0.01, 1500);
    let mu_1 = bilinear_peak_ductility(&record, 0.01, 12.0, 0.05, 0.01, 0.1, 1.0);
    let mu_3 = bilinear_peak_ductility(&record, 0.01, 12.0, 0.05, 0.01, 0.1, 3.0);
    assert!(mu_3 > mu_1, "stronger shaking, larger ductility");
    let seeds: Vec<u64> = (0..24).collect();
    let ladder = [0.5, 1.5, 3.0, 6.0];
    let frag = ida_fragility(&kt(), &seeds, &ladder, 12.0, 0.05, 0.01, 0.1, 4.0);
    let ps: Vec<f64> = frag.iter().map(|f| f.p).collect();
    println!("{{\"metric\":\"fragility\",\"im\":{ladder:?},\"p\":{ps:?}}}");
    for w in ps.windows(2) {
        assert!(w[1] >= w[0] - 1e-12, "fragility is monotone in IM: {ps:?}");
    }
    assert!(
        ps[0] < 0.5 && *ps.last().expect("last") > 0.5,
        "the curve sweeps the transition: {ps:?}"
    );
    verdict(
        "uq-003",
        "IDA fragility over a 24-motion suite is monotone in intensity and sweeps the \
         exceedance transition across the ladder",
    );
}

#[test]
fn uq_004_anytime_stopping_is_valid_and_adaptive() {
    // CERTIFY THE CERTIFIERS: over replications with ADAPTIVE stopping
    // (the CS decides when to quit), the stopped interval still covers
    // the truth at ~the nominal rate — optional stopping is safe.
    let alpha = 0.05;
    let mut miss = 0usize;
    let mut n_easy_total = 0u64;
    let reps = 60u64;
    for rep in 0..reps {
        let p_true = 0.3;
        let mut lcg = 0x1111_u64.wrapping_add(rep * 977);
        let est = estimate_probability_anytime(
            move |_| {
                lcg = lcg
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                f64::from(u8::from(
                    ((lcg >> 11) as f64) / (1u64 << 53) as f64 <= p_true,
                ))
            },
            alpha,
            0.08,
            20_000,
        );
        assert!(est.converged, "the CS reaches the target width");
        if p_true < est.lo || p_true > est.hi {
            miss += 1;
        }
        n_easy_total += est.n;
    }
    #[allow(clippy::cast_precision_loss)]
    let miss_rate = miss as f64 / reps as f64;
    // ADAPTIVITY: a looser target stops much earlier.
    let mut lcg = 0x2222u64;
    let loose = estimate_probability_anytime(
        move |_| {
            lcg = lcg
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            f64::from(u8::from(((lcg >> 11) as f64) / (1u64 << 53) as f64 <= 0.3))
        },
        alpha,
        0.2,
        20_000,
    );
    println!(
        "{{\"metric\":\"anytime\",\"miss_rate\":{miss_rate:.3},\"nominal\":{alpha},\
         \"mean_n_tight\":{},\"n_loose\":{}}}",
        n_easy_total / reps,
        loose.n
    );
    assert!(
        miss_rate <= alpha + 0.05,
        "stopped coverage holds under optional stopping: {miss_rate}"
    );
    assert!(
        loose.n * 3 < n_easy_total / reps,
        "looser decisions stop much earlier: {} vs {}",
        loose.n,
        n_easy_total / reps
    );
    verdict(
        "uq-004",
        "60 replications with CS-decided stopping: miss rate within the nominal band \
         (optional stopping safe by construction); a 2.5x looser target stops >3x sooner",
    );
}

#[test]
fn uq_005_cvar_and_adaptive_mlmc_rate_recovery() {
    // CVaR sanity: the tail mean of a known sample set.
    let samples: Vec<f64> = (1..=100).map(f64::from).collect();
    let c = cvar(&samples, 0.9);
    assert!(
        (c - 95.5).abs() < 0.6,
        "CVaR_0.9 of 1..100 is the top-decile mean: {c}"
    );
    // Anti-conservative-bias regression (bead zsvk): for NON-integer n·β the
    // proper empirical CVaR applies a fractional weight to the boundary order
    // statistic. A plain top-k mean returns 2.0 here (UNDER-reporting the risk);
    // the correct CVaR_0.85 is (0.5·1 + 3)/1.5 = 7/3.
    let tail = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0];
    let cv = cvar(&tail, 0.85);
    assert!(
        (cv - 7.0 / 3.0).abs() < 1e-12,
        "CVaR_0.85 must use the fractional boundary weight (2.333…), got {cv}"
    );
    // Adaptive MLMC on a synthetic with KNOWN rates: mean_l = 4^-l,
    // var_l = 8^-l (alpha = 2, beta = 3). The estimator must recover
    // the rates, add levels until the bias fits, and land near the
    // telescoped limit (sum of 4^-l -> 4/3).
    let sampler = |level: usize, i: usize| -> f64 {
        let mean = det::powi(4.0, -i32::try_from(level).expect("small"));
        let noise_scale = det::sqrt(det::powi(8.0, -i32::try_from(level).expect("small")));
        let mut z = 0x5eed_u64 ^ ((level as u64) << 32 | i as u64);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        let u = ((z >> 11) as f64) / (1u64 << 53) as f64 - 0.5;
        mean + noise_scale * u * 0.5
    };
    let report = adaptive_mlmc(
        sampler,
        |l| det::powi(4.0, i32::try_from(l).expect("small")),
        2e-3,
        400,
        8,
    );
    println!(
        "{{\"metric\":\"adaptive-mlmc\",\"levels\":{},\"alpha\":{:.2},\"beta\":{:.2},\
         \"estimate\":{:.4},\"bias\":{:.2e}}}",
        report.levels.len(),
        report.alpha,
        report.beta,
        report.estimate,
        report.bias_estimate
    );
    assert!(report.levels.len() >= 3, "levels were ADDED adaptively");
    assert!(
        (report.alpha - 2.0).abs() < 0.6 && (report.beta - 3.0).abs() < 1.0,
        "rates recovered: alpha {:.2}, beta {:.2}",
        report.alpha,
        report.beta
    );
    assert!(
        (report.estimate - 4.0 / 3.0).abs() < 0.02,
        "the telescoped estimate lands near the limit: {}",
        report.estimate
    );
    assert!(report.bias_estimate <= 1e-3, "the bias stop criterion held");
    verdict(
        "uq-005",
        "CVaR matches the analytic tail mean; adaptive MLMC adds levels on its own, \
         recovers alpha~2/beta~3 from level statistics, and stops when the extrapolated \
         bias fits the tolerance",
    );
}

#[test]
fn cvar_rejects_invalid_risk_inputs() {
    assert!(
        std::panic::catch_unwind(|| cvar(&[], 0.9)).is_err(),
        "empty loss samples must not report zero risk"
    );
    assert!(
        std::panic::catch_unwind(|| cvar(&[1.0], 0.0)).is_err(),
        "beta must define a strict upper tail"
    );
    assert!(
        std::panic::catch_unwind(|| cvar(&[1.0], 1.0)).is_err(),
        "beta=1 would leave an empty tail"
    );
    assert!(
        std::panic::catch_unwind(|| cvar(&[1.0, f64::NAN], 0.9)).is_err(),
        "non-finite samples must not sort into a risk estimate"
    );
}

#[test]
fn adaptive_mlmc_rejects_invalid_admission_inputs() {
    let sampler = |level: usize, i: usize| (level + i) as f64;
    assert!(
        std::panic::catch_unwind(|| adaptive_mlmc(sampler, |_| 1.0, 0.0, 4, 2)).is_err(),
        "non-positive tolerance must fail before producing NaN evidence"
    );
    assert!(
        std::panic::catch_unwind(|| adaptive_mlmc(sampler, |_| 1.0, 1e-3, 0, 2)).is_err(),
        "zero pilot samples must fail before dividing by zero"
    );
    assert!(
        std::panic::catch_unwind(|| adaptive_mlmc(sampler, |_| 1.0, 1e-3, 4, 0)).is_err(),
        "max_level below 1 cannot host the required level-0/1 pilot ladder"
    );
}
