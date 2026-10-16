//! Perf-regression-CI conformance (the fz2.4 bead): gate arithmetic,
//! change-point calibration on synthetic series (zero false alarms at
//! the declared confidence), seeded-regression attribution (the red
//! arrives WITH its flame-graph-level diagnosis), and the dashboard
//! one-liner answering the canonical question in one call.

use std::collections::BTreeMap;

use fs_roofline::regress::{
    Cusum, GateSpec, GateVerdict, MAX_REGRESSION_DASHBOARD_NIGHTS,
    MAX_REGRESSION_DASHBOARD_PHASE_OBSERVATIONS, MAX_REGRESSION_HISTORY_NIGHTS,
    MAX_REGRESSION_KERNELS, MAX_REGRESSION_PHASES_PER_NIGHT, MAX_REGRESSION_SERIES_SAMPLES, Night,
    gate, slower_this_month, standardize,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-roofline/regress\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn unit(seed: u64, k: u64) -> f64 {
    let mut z = seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(k + 1);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

fn gauss(seed: u64, k: u64) -> f64 {
    (0..12).map(|j| unit(seed, k * 12 + j)).sum::<f64>() - 6.0
}

/// A stable kernel: attainment ~ N(0.72, 0.01) with steady phases.
fn stable_night(night: u64, seed: u64) -> fs_roofline::regress::Night {
    let mut phases = BTreeMap::new();
    phases.insert(
        "assemble".to_string(),
        0.30 + 0.005 * gauss(seed, night * 3),
    );
    phases.insert(
        "solve".to_string(),
        0.55 + 0.005 * gauss(seed, night * 3 + 1),
    );
    phases.insert(
        "reduce".to_string(),
        0.15 + 0.003 * gauss(seed, night * 3 + 2),
    );
    fs_roofline::regress::Night {
        night,
        attainment: 0.72 + 0.01 * gauss(seed, night * 7 + 5),
        phases,
    }
}

#[test]
fn rg_001_noise_robustness_zero_false_alarms() {
    // 60 nights of stable code x 20 independent kernels: ZERO alarms
    // from the gate AND the CUSUM at the declared settings — thermal
    // jitter must not cry wolf.
    let mut gate_alarms = 0usize;
    let mut cusum_alarms = 0usize;
    for kernel in 0..20u64 {
        let history: Vec<_> = (0..60).map(|n| stable_night(n, 0xace + kernel)).collect();
        for t in 10..60 {
            if let GateVerdict::Red { .. } = gate(&history[..=t], GateSpec::default()) {
                gate_alarms += 1;
            }
        }
        let xs: Vec<f64> = history.iter().map(|n| n.attainment).collect();
        let z = standardize(&xs, 8).expect("bounded history");
        if Cusum::default().first_alarm(&z).is_some() {
            cusum_alarms += 1;
        }
    }
    println!(
        "{{\"metric\":\"false-alarms\",\"nights\":60,\"kernels\":20,\
         \"gate_alarms\":{gate_alarms},\"cusum_alarms\":{cusum_alarms}}}"
    );
    assert_eq!(gate_alarms, 0, "the 4-sigma band never cries wolf");
    assert_eq!(cusum_alarms, 0, "the CUSUM never cries wolf on stable code");
    verdict(
        "rg-001",
        "20 kernels x 60 stable nights: zero gate alarms and zero CUSUM alarms at the \
         declared settings — dispersion-aware bands hold",
    );
}

#[test]
fn rg_002_seeded_regression_red_with_attribution() {
    // Night 30 de-tunes the SOLVE phase (2x slower): the gate must go
    // RED with `solve` as the top attribution — the regression arrives
    // with its own diagnosis.
    let mut history: Vec<_> = (0..30).map(|n| stable_night(n, 0xbead)).collect();
    let mut bad = stable_night(30, 0xbead);
    bad.attainment = 0.48; // the de-tuned kernel's roofline drop
    *bad.phases.get_mut("solve").expect("solve") *= 2.0;
    history.push(bad);
    let verdict_ = gate(&history, GateSpec::default());
    let GateVerdict::Red { z, attribution } = verdict_ else {
        panic!("the seeded regression must gate red: {verdict_:?}")
    };
    println!(
        "{{\"metric\":\"seeded-regression\",\"z\":{z:.1},\"top\":\"{}\",\
         \"shares\":[{:.3},{:.3}]}}",
        attribution[0].0, attribution[0].1, attribution[0].2
    );
    assert!(z < -4.0, "far outside the band: z = {z:.1}");
    assert_eq!(
        attribution[0].0, "solve",
        "the flame-graph diff names the phase"
    );
    assert!(
        attribution[0].2 > attribution[0].1 + 0.1,
        "the share growth is visible: {:.3} -> {:.3}",
        attribution[0].1,
        attribution[0].2
    );
    verdict(
        "rg-002",
        "the de-tuned solve phase gates red at z < -4 with `solve` ranked first in the \
         flame-graph-equivalent attribution",
    );
}

#[test]
fn rg_003_cusum_catches_the_slow_drift() {
    // A 0.3-sigma-per-night drift never trips the single-night gate
    // but MUST trip the CUSUM within the month — the complementary
    // detector pair.
    let mut xs: Vec<f64> = (0..20).map(|n| 0.72 + 0.01 * gauss(0xd1f7, n)).collect();
    for n in 20..50u64 {
        let drift = 0.003 * (n - 19) as f64;
        xs.push(0.72 - drift + 0.01 * gauss(0xd1f7, n));
    }
    let z = standardize(&xs, 8).expect("bounded history");
    let single_night_reds = z.iter().skip(10).filter(|&&v| v < -4.0).count();
    let alarm = Cusum::default().first_alarm(&z);
    println!(
        "{{\"metric\":\"drift\",\"single_night_reds\":{single_night_reds},\
         \"cusum_alarm_at\":{alarm:?}}}"
    );
    let at = alarm.expect("the CUSUM catches the drift");
    assert!(at < 45, "caught within the month: night {at}");
    verdict(
        "rg-003",
        "a 0.3-sigma/night drift trips the CUSUM mid-month — the change-point detector \
         covers what the single-night band cannot",
    );
}

#[test]
fn rg_004_dashboard_one_liner() {
    // Three kernels: one regressed (12% drop, reduce-phase bloat), two
    // stable. The canonical question answers in ONE call, ranked, with
    // the why attached.
    let mut kernels = BTreeMap::new();
    kernels.insert(
        "gemm".to_string(),
        (0..30).map(|n| stable_night(n, 1)).collect::<Vec<_>>(),
    );
    kernels.insert(
        "spmv".to_string(),
        (0..30).map(|n| stable_night(n, 2)).collect::<Vec<_>>(),
    );
    let mut regressed: Vec<_> = (0..30).map(|n| stable_night(n, 3)).collect();
    for night in regressed.iter_mut().skip(20) {
        night.attainment *= 0.86;
        *night.phases.get_mut("reduce").expect("reduce") *= 3.0;
    }
    kernels.insert("fft".to_string(), regressed);
    let report = slower_this_month(&kernels, 5.0).expect("finite non-negative threshold");
    println!("{{\"metric\":\"dashboard\",\"report\":{report:?}}}");
    assert_eq!(report.len(), 1, "only the regressed kernel is named");
    assert_eq!(report[0].0, "fft");
    assert!(
        report[0].1 > 10.0,
        "the drop percentage is right: {:.1}",
        report[0].1
    );
    assert_eq!(report[0].2, "reduce", "and the why names the bloated phase");
    verdict(
        "rg-004",
        "'what got slower this month, and why' answers in one call: fft, ~13%, reduce — \
         stable kernels stay unnamed",
    );
}

/// rg-005 (bead fz2.4.1): malformed evidence FAILS CLOSED — NaN or
/// infinite attainment, negative attainment, non-finite or negative
/// phase durations, and unusable specs all yield Invalid, never Green,
/// each with a diagnosis.
#[test]
fn rg_005_malformed_evidence_fails_closed() {
    let good = |night: u64| Night {
        night,
        attainment: 0.8,
        phases: BTreeMap::from([("solve".to_string(), 1.0), ("io".to_string(), 0.5)]),
    };
    let mut history: Vec<Night> = (0..12).map(good).collect();
    let spec = GateSpec::default();
    // Baseline sanity: the clean history gates Green.
    assert!(matches!(gate(&history, spec), GateVerdict::Green { .. }));
    for (label, poison) in [
        ("nan-attainment", f64::NAN),
        ("inf-attainment", f64::INFINITY),
        ("neg-inf-attainment", f64::NEG_INFINITY),
        ("negative-attainment", -0.25),
    ] {
        let mut h = history.clone();
        h[11].attainment = poison;
        let v = gate(&h, spec);
        assert!(
            matches!(v, GateVerdict::Invalid { .. }),
            "{label}: expected Invalid, got {v:?}"
        );
    }
    // Poison BURIED in the baseline (not the newest night) also refuses.
    history[3].phases.insert("solve".to_string(), f64::NAN);
    assert!(matches!(gate(&history, spec), GateVerdict::Invalid { .. }));
    history[3].phases.insert("solve".to_string(), -2.0);
    assert!(matches!(gate(&history, spec), GateVerdict::Invalid { .. }));
    // Unusable specs.
    let clean: Vec<Night> = (0..12).map(good).collect();
    for bad_spec in [
        GateSpec {
            k_sigma: f64::NAN,
            min_baseline: 8,
        },
        GateSpec {
            k_sigma: 0.0,
            min_baseline: 8,
        },
        GateSpec {
            k_sigma: -1.0,
            min_baseline: 8,
        },
        GateSpec {
            k_sigma: f64::INFINITY,
            min_baseline: 8,
        },
        GateSpec {
            k_sigma: 4.0,
            min_baseline: 1,
        },
    ] {
        assert!(
            matches!(gate(&clean, bad_spec), GateVerdict::Invalid { .. }),
            "spec {bad_spec:?} must be refused"
        );
    }
    println!(
        "{{\"suite\":\"fs-roofline/regress\",\"case\":\"rg-005\",\"verdict\":\"pass\",\
         \"detail\":\"NaN/inf/negative fields and unusable specs all Invalid, never Green\"}}"
    );
}

/// rg-006 (bead fz2.4.1): METAMORPHIC — phase durations are shares, so
/// rescaling every phase by a constant (a time-unit change, seconds to
/// milliseconds) preserves the verdict AND the attribution ranking.
#[test]
fn rg_006_time_unit_invariance() {
    let mk = |night: u64, att: f64, solve: f64, reduce: f64, scale: f64| Night {
        night,
        attainment: att,
        phases: BTreeMap::from([
            ("solve".to_string(), solve * scale),
            ("reduce".to_string(), reduce * scale),
        ]),
    };
    for scale in [1.0f64, 1000.0] {
        let mut history: Vec<Night> = (0..14)
            .map(|t| mk(t, 0.80 + 0.001 * (t % 3) as f64, 2.0, 1.0, scale))
            .collect();
        // The regressed night: attainment collapses, reduce blows up.
        history.push(mk(14, 0.40, 2.0, 4.0, scale));
        let v = gate(&history, GateSpec::default());
        let GateVerdict::Red { attribution, .. } = &v else {
            panic!("expected Red at scale {scale}, got {v:?}");
        };
        assert_eq!(
            attribution[0].0, "reduce",
            "top offender is scale-invariant (scale {scale})"
        );
    }
    println!(
        "{{\"suite\":\"fs-roofline/regress\",\"case\":\"rg-006\",\"verdict\":\"pass\",\
         \"detail\":\"verdict and attribution ranking invariant under time-unit rescaling\"}}"
    );
}

/// rg-007 (bead fz2.4.1): poisoned trend/CUSUM state fails closed —
/// non-finite residuals alarm at their index instead of resetting the
/// shortfall; standardize maps poisoned history to −∞ from the first
/// bad index; an invalid detector spec cannot certify quiet; and
/// slower_this_month flags the poisoned kernel loudest instead of
/// skipping it.
#[test]
fn rg_007_poison_never_enters_state() {
    // NaN in the residual stream: alarm AT the poison, not suppression.
    let mut z = vec![0.1f64; 30];
    z[17] = f64::NAN;
    assert_eq!(Cusum::default().first_alarm(&z), Some(17));
    // Clean quiet stream stays quiet.
    assert_eq!(Cusum::default().first_alarm(&vec![0.1f64; 30]), None);
    // Invalid detector: cannot certify quiet.
    let bad = Cusum {
        k: f64::NAN,
        h: 8.0,
    };
    assert_eq!(bad.first_alarm(&[0.0, 0.0]), Some(0));
    // standardize: poison propagates as -inf from the first bad index.
    let zs =
        standardize(&[1.0, 1.0, 1.0, f64::INFINITY, 1.0], 2).expect("bounded poisoned history");
    assert!(zs[..3].iter().all(|v| v.is_finite()));
    assert!(zs[3] == f64::NEG_INFINITY && zs[4] == f64::NEG_INFINITY);
    // ...and the -inf stream alarms downstream.
    assert_eq!(Cusum::default().first_alarm(&zs), Some(3));
    // slower_this_month: poisoned kernel flagged first with INVALID why.
    let good_hist: Vec<Night> = (0..14)
        .map(|t| Night {
            night: t,
            attainment: 0.8,
            phases: BTreeMap::from([("solve".to_string(), 1.0)]),
        })
        .collect();
    let mut poisoned = good_hist.clone();
    poisoned[9].attainment = f64::NAN;
    let kernels = BTreeMap::from([
        ("clean".to_string(), good_hist),
        ("poisoned".to_string(), poisoned),
    ]);
    let report = slower_this_month(&kernels, 5.0).expect("finite non-negative threshold");
    assert_eq!(
        report.len(),
        1,
        "clean kernel has no drop; poisoned is flagged"
    );
    assert_eq!(report[0].0, "poisoned");
    assert!(report[0].1.is_infinite() && report[0].2.starts_with("INVALID"));
    println!(
        "{{\"suite\":\"fs-roofline/regress\",\"case\":\"rg-007\",\"verdict\":\"pass\",\
         \"detail\":\"poison alarms instead of suppressing; invalid kernels flagged loudest\"}}"
    );
}

/// rg-008: logical time and dashboard configuration are evidence, not hints.
/// Duplicate/reversed nights and invalid thresholds must fail closed, while the
/// inclusive zero bound and the dashboard's strict `drop > floor` boundary
/// remain usable and deterministic.
#[test]
fn rg_008_history_order_and_thresholds_fail_closed() {
    let history = |attainment: f64| {
        (0..14)
            .map(|night| Night {
                night,
                attainment,
                phases: BTreeMap::from([("solve".to_string(), 1.0)]),
            })
            .collect::<Vec<_>>()
    };

    let mut duplicate = history(0.8);
    duplicate[8].night = duplicate[7].night;
    let duplicate_gate = gate(&duplicate, GateSpec::default());
    assert!(
        matches!(duplicate_gate, GateVerdict::Invalid { .. }),
        "duplicate logical time must not be treated as chronological: {duplicate_gate:?}"
    );

    let mut reversed = history(0.8);
    reversed.swap(7, 8);
    let reversed_gate = gate(&reversed, GateSpec::default());
    assert!(
        matches!(reversed_gate, GateVerdict::Invalid { .. }),
        "reversed logical time must not be treated as chronological: {reversed_gate:?}"
    );

    let kernels = BTreeMap::from([
        ("duplicate".to_string(), duplicate),
        ("reversed".to_string(), reversed),
    ]);
    let invalid_rows = slower_this_month(&kernels, 0.0).expect("zero is a valid threshold");
    assert_eq!(invalid_rows.len(), 2);
    assert!(invalid_rows.iter().all(|row| {
        row.1.is_infinite()
            && row
                .2
                .starts_with("INVALID: logical night must increase strictly")
    }));

    for threshold in [f64::NAN, f64::INFINITY, -f64::MIN_POSITIVE] {
        let error = slower_this_month(&BTreeMap::new(), threshold)
            .expect_err("non-finite or negative threshold must fail closed");
        assert!(error.reason().contains("pct_floor"));
    }

    let mut exact = history(1.0);
    for night in exact.iter_mut().skip(7) {
        night.attainment = 0.5;
    }
    let exact = BTreeMap::from([("exact-bound".to_string(), exact)]);
    assert!(
        slower_this_month(&exact, 50.0)
            .expect("exact finite threshold")
            .is_empty(),
        "the documented 'more than pct_floor' comparison is strict at equality"
    );

    println!(
        "{{\"suite\":\"fs-roofline/regress\",\"case\":\"rg-008\",\"verdict\":\"pass\",\
         \"detail\":\"logical time is strictly increasing and dashboard thresholds fail closed\"}}"
    );
}

#[test]
fn rg_009_extreme_finite_values_and_sparse_phases_remain_sound() {
    let extreme: Vec<Night> = (0..14)
        .map(|night| Night {
            night,
            attainment: f64::MAX,
            phases: BTreeMap::from([
                ("left".to_string(), f64::MAX),
                ("right".to_string(), f64::MAX),
            ]),
        })
        .collect();
    assert!(matches!(
        gate(&extreme, GateSpec::default()),
        GateVerdict::Green { z } if z == 0.0
    ));
    assert!(
        slower_this_month(&BTreeMap::from([("extreme".to_string(), extreme)]), 0.0)
            .expect("bounded extreme history")
            .is_empty(),
        "finite values must not overflow into a silent NaN trend"
    );
    assert_eq!(
        standardize(&[f64::MAX, f64::MAX, f64::MAX], 1).expect("stable extreme series is bounded"),
        vec![0.0, 0.0, 0.0],
        "a stable extreme series must not underflow its normalized floor into 0/0 poison"
    );

    let mut sparse: Vec<Night> = (0..8)
        .map(|night| Night {
            night,
            attainment: 1.0,
            phases: if night == 0 {
                BTreeMap::from([("common".to_string(), 1.0), ("rare".to_string(), 9.0)])
            } else {
                BTreeMap::from([("common".to_string(), 1.0)])
            },
        })
        .collect();
    sparse.push(Night {
        night: 8,
        attainment: 0.0,
        phases: BTreeMap::from([("common".to_string(), 4.0), ("rare".to_string(), 6.0)]),
    });
    let GateVerdict::Red { attribution, .. } = gate(&sparse, GateSpec::default()) else {
        panic!("sparse-phase regression must gate red");
    };
    assert_eq!(attribution[0].0, "rare");
    assert_eq!(
        attribution[0].1.to_bits(),
        0.0_f64.to_bits(),
        "seven absent nights contribute positive zero"
    );
    assert!((attribution[0].2 - 0.6).abs() < 1e-12);
}

#[test]
fn rg_010_public_regression_inputs_are_bounded() {
    let night = |night| Night {
        night,
        attainment: 1.0,
        phases: BTreeMap::new(),
    };
    let oversized_history: Vec<_> = (0..=MAX_REGRESSION_HISTORY_NIGHTS as u64)
        .map(night)
        .collect();
    assert!(matches!(
        gate(&oversized_history, GateSpec::default()),
        GateVerdict::Invalid { .. }
    ));

    let mut too_many_phases = night(0);
    too_many_phases.phases = (0..=MAX_REGRESSION_PHASES_PER_NIGHT)
        .map(|index| (format!("phase-{index}"), 1.0))
        .collect();
    assert!(matches!(
        gate(&[too_many_phases], GateSpec::default()),
        GateVerdict::Invalid { .. }
    ));

    let kernels: BTreeMap<_, _> = (0..=MAX_REGRESSION_KERNELS)
        .map(|index| (format!("kernel-{index}"), Vec::new()))
        .collect();
    assert!(slower_this_month(&kernels, 0.0).is_err());
    assert!(standardize(&vec![0.0; MAX_REGRESSION_SERIES_SAMPLES + 1], 8).is_err());

    for min_baseline in [MAX_REGRESSION_HISTORY_NIGHTS, usize::MAX] {
        assert!(matches!(
            gate(
                &[],
                GateSpec {
                    k_sigma: 4.0,
                    min_baseline,
                }
            ),
            GateVerdict::Invalid { .. }
        ));
    }
    assert!(matches!(
        gate(
            &[],
            GateSpec {
                k_sigma: 4.0,
                min_baseline: MAX_REGRESSION_HISTORY_NIGHTS - 1,
            }
        ),
        GateVerdict::Green { z: 0.0 }
    ));

    let mut dashboard_nights = BTreeMap::new();
    for kernel in 0..(MAX_REGRESSION_DASHBOARD_NIGHTS / MAX_REGRESSION_HISTORY_NIGHTS) {
        dashboard_nights.insert(
            format!("history-{kernel}"),
            (0..MAX_REGRESSION_HISTORY_NIGHTS as u64)
                .map(night)
                .collect(),
        );
    }
    dashboard_nights.insert("history-overflow".to_string(), vec![night(0)]);
    let error = slower_this_month(&dashboard_nights, 0.0)
        .expect_err("dashboard-wide night limit+1 must refuse");
    assert!(error.reason().contains("dashboard night count"));

    let phases: BTreeMap<_, _> = (0..MAX_REGRESSION_PHASES_PER_NIGHT)
        .map(|index| (format!("phase-{index}"), 1.0))
        .collect();
    let full_phase_nights =
        MAX_REGRESSION_DASHBOARD_PHASE_OBSERVATIONS / MAX_REGRESSION_PHASES_PER_NIGHT;
    let phase_history: Vec<_> = (0..full_phase_nights as u64)
        .map(|night| Night {
            night,
            attainment: 1.0,
            phases: phases.clone(),
        })
        .collect();
    let phase_overflow = Night {
        night: 0,
        attainment: 1.0,
        phases: BTreeMap::from([("phase-overflow".to_string(), 1.0)]),
    };
    let error = slower_this_month(
        &BTreeMap::from([
            ("phase-full".to_string(), phase_history),
            ("phase-overflow".to_string(), vec![phase_overflow]),
        ]),
        0.0,
    )
    .expect_err("dashboard-wide phase observation limit+1 must refuse");
    assert!(error.reason().contains("dashboard phase observation count"));
}

#[test]
fn rg_011_tiny_regressions_and_extreme_improvements_are_scale_sound() {
    let phase = || BTreeMap::from([("solve".to_string(), 1.0)]);
    let mut gate_history: Vec<_> = (0..8)
        .map(|night| Night {
            night,
            attainment: 1e-20,
            phases: phase(),
        })
        .collect();
    gate_history.push(Night {
        night: 8,
        attainment: 0.0,
        phases: phase(),
    });
    assert!(matches!(
        gate(&gate_history, GateSpec::default()),
        GateVerdict::Red { z, .. } if z == f64::NEG_INFINITY
    ));

    let tiny: Vec<_> = (0..14)
        .map(|night| Night {
            night,
            attainment: if night < 7 { 1e-20 } else { 0.0 },
            phases: phase(),
        })
        .collect();
    let tiny_report = slower_this_month(&BTreeMap::from([("tiny".to_string(), tiny)]), 99.0)
        .expect("tiny finite values are valid evidence");
    assert_eq!(tiny_report.len(), 1);
    assert_eq!(tiny_report[0].0, "tiny");
    assert_eq!(tiny_report[0].1.to_bits(), 100.0_f64.to_bits());

    let improvement: Vec<_> = (0..14)
        .map(|night| Night {
            night,
            attainment: if night < 7 { 0.0 } else { f64::MAX },
            phases: phase(),
        })
        .collect();
    assert!(
        slower_this_month(
            &BTreeMap::from([("extreme-improvement".to_string(), improvement)]),
            0.0,
        )
        .expect("extreme finite improvement is valid evidence")
        .is_empty(),
        "an improvement must not overflow into an INVALID regression"
    );
}

#[test]
fn rg_012_dashboard_attribution_uses_its_opening_and_trailing_windows() {
    let history: Vec<_> = (0..30)
        .map(|night| {
            let regressed = night >= 7;
            Night {
                night,
                attainment: if regressed { 0.5 } else { 1.0 },
                phases: if regressed {
                    BTreeMap::from([("base".to_string(), 1.0), ("slow".to_string(), 9.0)])
                } else {
                    BTreeMap::from([("base".to_string(), 9.0), ("slow".to_string(), 1.0)])
                },
            }
        })
        .collect();
    let report = slower_this_month(
        &BTreeMap::from([("sustained-shift".to_string(), history)]),
        5.0,
    )
    .expect("bounded sustained regression");
    assert_eq!(report.len(), 1);
    assert_eq!(report[0].0, "sustained-shift");
    assert!((report[0].1 - 50.0).abs() < f64::EPSILON);
    assert_eq!(
        report[0].2, "slow",
        "the trailing-window phase increase must not be diluted by regressed middle nights"
    );
}

#[test]
fn rg_013_standardization_is_scale_invariant_without_false_improvement_alarm() {
    let canonical = [1.0, 1.01, 0.99, 1.005, 0.995, 1.002, 0.998, 1.0, 0.5, 0.5];
    let tiny: Vec<_> = canonical.iter().map(|value| value * 1e-200).collect();
    let canonical_z = standardize(&canonical, 8).expect("bounded canonical series");
    let tiny_z = standardize(&tiny, 8).expect("bounded rescaled series");
    for (canonical, tiny) in canonical_z.iter().zip(&tiny_z) {
        let tolerance = 1e-10 * canonical.abs().max(1.0);
        assert!(
            (canonical - tiny).abs() <= tolerance,
            "standardized score changed under positive rescaling: {canonical} vs {tiny}"
        );
    }
    let canonical_alarm = Cusum::default().first_alarm(&canonical_z);
    assert_eq!(canonical_alarm, Some(8));
    assert_eq!(Cusum::default().first_alarm(&tiny_z), canonical_alarm);

    let improvement =
        standardize(&[1.0, 1.0, 1.0, f64::MAX], 2).expect("bounded extreme improvement series");
    assert_eq!(improvement[3].to_bits(), f64::MAX.to_bits());
    assert_eq!(
        Cusum::default().first_alarm(&improvement),
        None,
        "a positive zero-dispersion improvement must not become poison or shortfall"
    );

    let decline =
        standardize(&[1.0, 1.0, 1.0, 0.0], 2).expect("bounded zero-dispersion decline series");
    assert_eq!(decline[3].to_bits(), (-f64::MAX).to_bits());
    assert_eq!(Cusum::default().first_alarm(&decline), Some(3));
}
