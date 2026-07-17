//! Conformal-hardening conformance (the 7tv.9 bead; runs under
//! `conformal-hardening`). Acceptance: per-bucket coverage verified
//! empirically across regime shifts (Mondrian catches what marginal
//! hides); drift e-tests detect seeded covariate shift within a
//! documented sample budget with NO false alarm on the null (the
//! certify-the-certifiers battery); FCR budgeting correct under
//! simultaneous monitoring; and the adversarial test that matters
//! most — a surrogate-exploiting optimizer does NOT achieve silent
//! invalid coverage.
#![cfg(feature = "conformal-hardening")]

use fs_eproc::hardening::{
    BucketBand, CoverageClaim, DriftMonitor, ExchangeabilityCard, MondrianConformal,
    admission_alpha, fcr_flag,
};

const SUITE: &str = "fs-eproc/hardening";
const CH_001_INPUT_SEED: u64 = 0xa11ce;
const CH_002_INPUT_SEED: u64 = 0xd21f7;
const CH_003_INPUT_SEED: u64 = 0xfc12;
const CH_004_INPUT_SEED: u64 = 0x0b71;
const CH_005_INPUT_SEED: u64 = 0x5eed;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("hardening verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("hardening verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn measurement(case: &str, name: &str, json: String) {
    let identity = format!("{case}/measurement");
    let mut emitter = fs_obs::Emitter::new(SUITE, &identity);
    let line = emitter
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: name.to_string(),
                json,
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("hardening measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn finite_json(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "null".to_string()
    }
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / (1u64 << 53) as f64
    }

    /// Approximate standard normal (sum of 12 uniforms − 6).
    fn gauss(&mut self) -> f64 {
        (0..12).map(|_| self.next()).sum::<f64>() - 6.0
    }
}

#[test]
fn ch_001_mondrian_catches_what_marginal_hides() {
    // Two regimes: "smooth" residuals ~ 0.1|N|, "shock" ~ 0.5|N|.
    // The pooled marginal band undercovers the shock regime; the
    // Mondrian per-bucket bands hold coverage in BOTH.
    let mut rng = Lcg(CH_001_INPUT_SEED);
    let mut cal = MondrianConformal::new(50);
    for _ in 0..400 {
        cal.add("smooth", 0.1 * rng.gauss().abs());
        cal.add("shock", 0.5 * rng.gauss().abs());
    }
    let alpha = 0.1;
    let BucketBand::Calibrated {
        half_width: w_marg, ..
    } = cal.marginal_band(alpha)
    else {
        panic!("marginal calibrated")
    };
    let BucketBand::Calibrated {
        half_width: w_shock,
        ..
    } = cal.band("shock", alpha)
    else {
        panic!("shock calibrated")
    };
    // Fresh test draws per regime: empirical coverage.
    let cover = |w: f64, scale: f64, rng: &mut Lcg| -> f64 {
        let n = 4000;
        let hits = (0..n).filter(|_| scale * rng.gauss().abs() <= w).count();
        f64::from(u32::try_from(hits).expect("count")) / f64::from(n)
    };
    let marg_on_shock = cover(w_marg, 0.5, &mut rng);
    let mondrian_on_shock = cover(w_shock, 0.5, &mut rng);
    let nominal_json = finite_json(1.0 - alpha);
    let marginal_json = finite_json(marg_on_shock);
    let mondrian_json = finite_json(mondrian_on_shock);
    measurement(
        "ch-001",
        "bucket-coverage",
        format!(
            "{{\"metric\":\"bucket-coverage\",\"nominal\":{nominal_json},\
             \"marginal_on_shock\":{marginal_json},\
             \"mondrian_on_shock\":{mondrian_json},\
             \"input_seed\":{CH_001_INPUT_SEED}}}"
        ),
    );
    assert!(
        marg_on_shock < 0.86,
        "the pooled band SILENTLY undercovers the shock regime: {marg_on_shock}"
    );
    assert!(
        mondrian_on_shock >= 0.88,
        "the Mondrian bucket holds its guarantee: {mondrian_on_shock}"
    );
    // An unseen regime REFUSES rather than pretending.
    assert!(
        matches!(
            cal.band("cryogenic", alpha),
            BucketBand::Refused { have: 0, .. }
        ),
        "unseen bucket refuses"
    );
    verdict(
        "ch-001",
        true,
        "marginal band covers only ~0.8 on the high-noise regime while the Mondrian \
         bucket holds ~0.9 nominal; unseen regimes refuse instead of extrapolating",
        CH_001_INPUT_SEED,
    );
}

#[test]
fn ch_002_drift_etest_detects_within_budget_no_false_alarm() {
    let mut rng = Lcg(CH_002_INPUT_SEED);
    let train: Vec<f64> = (0..500).map(|_| rng.gauss()).collect();
    // NULL: same distribution — no false alarm across 2000 samples
    // (anytime validity, adversarial-stopping-proof by construction).
    let mut null_mon = DriftMonitor::new(train.clone(), 0.05);
    let mut fired = false;
    for _ in 0..2000 {
        fired |= null_mon.observe(rng.gauss()).drifted;
    }
    assert!(!fired, "no false alarm on the null over 2000 samples");
    // SEEDED covariate shift (+1 sigma): detected within the
    // documented budget of 200 samples; escalation shrinks validity.
    let mut mon = DriftMonitor::new(train, 0.05);
    let mut at = 0u64;
    let mut scale = 1.0f64;
    for _ in 0..400 {
        let v = mon.observe(rng.gauss() + 1.0);
        if v.drifted {
            at = v.samples_at_detection;
            scale = v.validity_scale;
            break;
        }
    }
    let scale_json = finite_json(scale);
    measurement(
        "ch-002",
        "drift-detection",
        format!(
            "{{\"metric\":\"drift-detection\",\"shift\":\"+1sigma\",\
             \"detected_at\":{at},\"validity_scale\":{scale_json},\
             \"input_seed\":{CH_002_INPUT_SEED}}}"
        ),
    );
    assert!(
        at > 0 && at <= 200,
        "detected within the 200-sample budget: {at}"
    );
    assert!(scale < 1.0, "the validity domain SHRINKS on drift: {scale}");
    verdict(
        "ch-002",
        true,
        "1-sigma covariate shift detected within 200 samples; zero false alarms across \
         2000 null samples; escalation shrinks the surrogate validity domain",
        CH_002_INPUT_SEED,
    );
}

#[test]
fn ch_003_fcr_budget_flags_broken_claims() {
    // 60 simultaneous coverage claims at alpha = 0.1; 10 are BROKEN
    // (true miss rate 0.4). The e-BH budget pass flags mostly the
    // broken ones (FDR controlled by the e-BH guarantee).
    let mut rng = Lcg(CH_003_INPUT_SEED);
    let mut claims: Vec<CoverageClaim> = (0..60)
        .map(|k| CoverageClaim::new(&format!("claim-{k:02}"), 0.1))
        .collect();
    for (k, claim) in claims.iter_mut().enumerate() {
        let miss = if k < 10 { 0.4 } else { 0.08 };
        for _ in 0..300 {
            claim.observe(rng.next() >= miss);
        }
    }
    let flagged = fcr_flag(&claims, 0.05);
    let true_pos = flagged.iter().filter(|&&i| i < 10).count();
    let false_pos = flagged.len() - true_pos;
    measurement(
        "ch-003",
        "fcr-budget",
        format!(
            "{{\"metric\":\"fcr-budget\",\"flagged\":{},\"true_pos\":{true_pos},\
             \"false_pos\":{false_pos},\"input_seed\":{CH_003_INPUT_SEED}}}",
            flagged.len()
        ),
    );
    assert!(true_pos >= 8, "most broken claims flagged: {true_pos}/10");
    assert!(
        false_pos <= 2,
        "false flags within the budget: {false_pos} of {}",
        flagged.len()
    );
    // Admission math: the per-claim reservation under a total budget.
    let per = admission_alpha(0.1, 1000);
    assert!((per - 1e-4).abs() < 1e-12, "Bonferroni reservation: {per}");
    verdict(
        "ch-003",
        true,
        "e-BH over 60 simultaneous miscoverage e-processes flags 8+ of 10 broken claims \
         with <=2 false flags; admission reserves budget/k per claim",
        CH_003_INPUT_SEED,
    );
}

#[test]
fn ch_004_adversarial_optimizer_shift() {
    // THE TEST THAT MATTERS MOST: a surrogate calibrated on region A
    // (x in [0,1], residual scale 0.1) is exploited by an optimizer
    // walking toward region B (x in [2,3], residual scale 0.6, where
    // the surrogate errs). WITHOUT hardening the A-calibrated band
    // silently undercovers on B. WITH hardening: the drift monitor
    // fires DURING the walk and the exchangeability card names the
    // policy — no silent invalid coverage.
    let mut rng = Lcg(CH_004_INPUT_SEED);
    let mut cal = MondrianConformal::new(50);
    let train_x: Vec<f64> = (0..300).map(|_| rng.next()).collect();
    for _ in 0..300 {
        cal.add("region-A", 0.1 * rng.gauss().abs());
    }
    let alpha = 0.1;
    let BucketBand::Calibrated {
        half_width: w_a, ..
    } = cal.band("region-A", alpha)
    else {
        panic!("calibrated")
    };
    // The optimizer's walk: candidates drift from x~0.5 to x~2.5.
    let mut mon = DriftMonitor::new(train_x, 0.05);
    let mut fired_at = 0u64;
    let mut naive_hits = 0u32;
    let mut naive_total = 0u32;
    for step in 0..300 {
        let t = f64::from(step) / 300.0;
        let x = 0.5 + 2.0 * t + 0.1 * rng.gauss();
        let v = mon.observe(x);
        if v.drifted && fired_at == 0 {
            fired_at = v.samples_at_detection;
        }
        // In the drifted region the true residual scale is 0.6.
        if t > 0.5 {
            naive_total += 1;
            if 0.6 * rng.gauss().abs() <= w_a {
                naive_hits += 1;
            }
        }
    }
    let naive_cov = f64::from(naive_hits) / f64::from(naive_total);
    let naive_cov_json = finite_json(naive_cov);
    let nominal_json = finite_json(1.0 - alpha);
    measurement(
        "ch-004",
        "adversarial-shift",
        format!(
            "{{\"metric\":\"adversarial-shift\",\"naive_coverage_on_B\":{naive_cov_json},\
             \"nominal\":{nominal_json},\"drift_fired_at\":{fired_at},\
             \"input_seed\":{CH_004_INPUT_SEED}}}"
        ),
    );
    assert!(
        naive_cov < 0.6,
        "the naive band's coverage silently collapses on region B: {naive_cov}"
    );
    assert!(
        fired_at > 0 && fired_at <= 250,
        "the drift monitor fires during the walk: {fired_at}"
    );
    // And the B bucket refuses until calibrated there.
    assert!(
        matches!(cal.band("region-B", alpha), BucketBand::Refused { .. }),
        "region-B refuses without calibration mass"
    );
    let card = ExchangeabilityCard {
        bucketing: "regime-class".to_string(),
        drift_alpha: 0.05,
        fcr_budget: 0.05,
        refresh_policy: "recalibrate-on-drift".to_string(),
    };
    assert!(card.to_json().contains("recalibrate-on-drift"));
    verdict(
        "ch-004",
        true,
        "the exploiting optimizer collapses naive coverage to ~0.5 on region B, but the \
         drift e-test fires mid-walk and the B bucket refuses — no silent invalid \
         coverage; the exchangeability card declares the policy",
        CH_004_INPUT_SEED,
    );
}

#[test]
fn ch_005_g5_determinism() {
    // All monitors replay bit-equal from the same stream.
    let run = || -> (f64, u64) {
        let mut rng = Lcg(CH_005_INPUT_SEED);
        let train: Vec<f64> = (0..200).map(|_| rng.gauss()).collect();
        let mut mon = DriftMonitor::new(train, 0.05);
        let mut last = 0.0;
        for _ in 0..300 {
            let v = mon.observe(rng.gauss() + 0.8);
            last = v.validity_scale;
        }
        (last, mon.observe(0.0).samples_at_detection)
    };
    let (a1, d1) = run();
    let (a2, d2) = run();
    assert!(a1.to_bits() == a2.to_bits() && d1 == d2, "bit-equal replay");
    verdict(
        "ch-005",
        true,
        "drift monitors replay bit-equal from the same stream (G5)",
        CH_005_INPUT_SEED,
    );
}

#[test]
fn ch_006_undercalibrated_conformal_band_is_infinite_not_undercovering() {
    // n=10 residuals, alpha=0.05: k = ceil(11*0.95) = 11 > 10, so NO finite
    // residual achieves 95% coverage — the honest split-conformal band is
    // INFINITE. Capping at the max residual would cover only 10/11 = 0.909
    // < 0.95, a false coverage certificate (bead q2tf).
    let mut cal = MondrianConformal::new(1);
    for i in 0..10 {
        cal.add("b", f64::from(i));
    }
    let band = cal.band("b", 0.05);
    let BucketBand::Calibrated { half_width, n } = band else {
        panic!("expected an infinite Calibrated band, got {band:?}");
    };
    assert_eq!(n, 10);
    assert!(
        half_width.is_infinite(),
        "under-calibrated band must be infinite, got {half_width}"
    );
    // Enough data (n=25, k = ceil(26*0.95) = 25 <= 25): a finite band.
    let mut big = MondrianConformal::new(1);
    for i in 0..25 {
        big.add("b", f64::from(i));
    }
    assert!(matches!(
        big.band("b", 0.05),
        BucketBand::Calibrated { half_width, .. } if half_width.is_finite()
    ));
    // Default construction is non-degenerate and refuses an empty bucket.
    assert!(matches!(
        MondrianConformal::default().band("empty", 0.1),
        BucketBand::Refused { have: 0, need: 1 }
    ));
}

#[test]
fn ch_007_malformed_validity_inputs_fail_closed() {
    assert!(std::panic::catch_unwind(|| MondrianConformal::new(0)).is_err());
    let mut cal = MondrianConformal::new(1);
    assert!(
        std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
            cal.add("bad", f64::NAN);
        }))
        .is_err()
    );
    assert!(cal.bucket_names().is_empty());
    assert!(std::panic::catch_unwind(|| cal.band("missing", 1.0)).is_err());
    assert!(std::panic::catch_unwind(|| DriftMonitor::new(Vec::new(), 0.05)).is_err());
    assert!(std::panic::catch_unwind(|| DriftMonitor::new(vec![f64::INFINITY], 0.05)).is_err());
    let mut monitor = DriftMonitor::new(vec![0.0, 1.0], 0.05);
    assert!(
        std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
            let _ = monitor.observe(f64::NAN);
        }))
        .is_err()
    );
    assert!(std::panic::catch_unwind(|| admission_alpha(0.1, 0)).is_err());
    assert!(std::panic::catch_unwind(|| admission_alpha(f64::NAN, 1)).is_err());
}

#[test]
fn ch_008_exchangeability_card_serializes_canonical_json_strings() {
    let card = ExchangeabilityCard {
        bucketing: "region-\"A\"\\cold".to_string(),
        drift_alpha: 0.0,
        fcr_budget: 0.05,
        refresh_policy: "line1\nline2".to_string(),
    };
    assert_eq!(
        card.to_json(),
        "{\"bucketing\":\"region-\\\"A\\\"\\\\cold\",\"drift_alpha\":0,\"fcr_budget\":0.05,\"refresh_policy\":\"line1\\nline2\"}"
    );
}
