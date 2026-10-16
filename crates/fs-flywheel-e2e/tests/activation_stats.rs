//! Activation-statistics battery (bead frankensim-sj31i.18): the
//! flywheel checkpoint can activate only from preregistered,
//! holdout-only, stratified evidence with exact denominators and
//! anytime-valid bounds. Rides the crate's `flywheel-e2e` feature gate
//! (the crate is feature-gated at the root); the statistical authority
//! itself has no loop dependencies.

use fs_flywheel_e2e::activation::{
    ActivationRefusal, ActivationThresholds, Outcome, Partition, Sample, SamplingFrame,
    StratumSpec, adjudicate,
};

fn thresholds() -> ActivationThresholds {
    ActivationThresholds {
        accept_rate_floor: 0.30,
        savings_floor: 1.5,
        confidence_delta: 0.05,
        betting_fraction: 0.5,
    }
}

fn frame(min_samples: u32) -> SamplingFrame {
    SamplingFrame {
        study: "flywheel-activation-study",
        seed: 0x5A31_0018,
        strata: vec![
            StratumSpec {
                name: "elliptic-low",
                min_samples,
            },
            StratumSpec {
                name: "elliptic-high",
                min_samples,
            },
        ],
        thresholds: thresholds(),
    }
}

fn winning_sample(key: u64, stratum: &'static str) -> Sample {
    Sample {
        key,
        stratum,
        partition: Partition::Holdout,
        outcome: Outcome::WarmStarted { cold: 30, warm: 10 },
    }
}

fn losing_sample(key: u64, stratum: &'static str) -> Sample {
    Sample {
        key,
        stratum,
        partition: Partition::Holdout,
        outcome: Outcome::ColdSolve,
    }
}

/// Strong evidence in both strata: mostly warm-started 3x savings with
/// a few cold solves.
fn strong_evidence(per_stratum: u64) -> Vec<Sample> {
    let mut samples = Vec::new();
    for (base, stratum) in [(0u64, "elliptic-low"), (1_000_000, "elliptic-high")] {
        for k in 0..per_stratum {
            let sample = if k % 8 == 7 {
                losing_sample(base + k, stratum)
            } else {
                winning_sample(base + k, stratum)
            };
            samples.push(sample);
        }
    }
    samples
}

#[test]
fn empty_and_below_minimum_evidence_refuses() {
    let frame = frame(12);
    assert!(matches!(
        adjudicate(&frame, &[]),
        Err(ActivationRefusal::InsufficientSamples {
            stratum: "elliptic-low",
            have: 0,
            need: 12,
        })
    ));
    // One below the preregistered minimum in one stratum refuses even
    // when the other stratum is rich.
    let mut samples = strong_evidence(11);
    samples.extend((11..40).map(|k| winning_sample(1_000_000 + k, "elliptic-high")));
    assert!(matches!(
        adjudicate(&frame, &samples),
        Err(ActivationRefusal::InsufficientSamples {
            stratum: "elliptic-low",
            have: 11,
            need: 12,
        })
    ));
}

#[test]
fn strong_holdout_evidence_activates_with_replayable_claim() {
    let frame = frame(12);
    let samples = strong_evidence(48);
    let report = adjudicate(&frame, &samples).expect("structurally sound evidence");
    assert!(
        report.activated,
        "strong evidence must activate: {report:?}"
    );
    for stratum in &report.strata {
        assert!(stratum.accept_rejected && stratum.savings_rejected);
        assert!(stratum.accept_e_max >= 20.0);
        assert!(stratum.savings_e_max >= 20.0);
        assert_eq!(stratum.samples, 48);
    }
    // Replay identity: the same frame and evidence give a bit-identical
    // report, bound to the frame identity.
    let replay = adjudicate(&frame, &samples).expect("replay");
    assert_eq!(report, replay);
    assert_eq!(report.frame_identity, frame.identity());
}

#[test]
fn correlated_duplicates_cannot_inflate_the_denominator() {
    let frame = frame(2);
    let mut samples = strong_evidence(4);
    samples.push(samples[0]);
    assert!(matches!(
        adjudicate(&frame, &samples),
        Err(ActivationRefusal::DuplicateSample { key }) if key == samples[0].key
    ));
}

#[test]
fn holdout_leakage_is_a_typed_refusal() {
    let frame = frame(2);
    let mut samples = strong_evidence(4);
    samples[3].partition = Partition::Development;
    let leaked = samples[3].key;
    assert!(matches!(
        adjudicate(&frame, &samples),
        Err(ActivationRefusal::HoldoutLeakage { key }) if key == leaked
    ));
}

#[test]
fn zero_warm_denominators_and_unmeasured_savings_refuse() {
    let loose_frame = frame(2);
    // The +INFINITY hole: a warm start with zero warm work.
    let mut samples = strong_evidence(4);
    samples[0].outcome = Outcome::WarmStarted { cold: 30, warm: 0 };
    let zero = samples[0].key;
    assert!(matches!(
        adjudicate(&loose_frame, &samples),
        Err(ActivationRefusal::ZeroWarmDenominator { key }) if key == zero
    ));

    // All-outright stratum: accepts are perfect but savings were never
    // measured; unmeasured savings must refuse, not activate.
    let outright: Vec<Sample> = (0..16)
        .map(|k| Sample {
            key: k,
            stratum: "elliptic-low",
            partition: Partition::Holdout,
            outcome: Outcome::AcceptedOutright,
        })
        .chain((0..16).map(|k| winning_sample(1_000_000 + k, "elliptic-high")))
        .collect();
    let strict_frame = frame(8);
    assert!(matches!(
        adjudicate(&strict_frame, &outright),
        Err(ActivationRefusal::NoWarmStartEvidence {
            stratum: "elliptic-low"
        })
    ));
}

#[test]
fn adverse_strata_cannot_hide_inside_a_favorable_pool() {
    // The favorable stratum is overwhelming; the adverse stratum meets
    // its minimum but its evidence is weak (distribution shift). The
    // pooled rate would sail past every floor — stratified adjudication
    // still refuses to activate.
    let frame = frame(12);
    let mut samples = strong_evidence(0);
    samples.extend((0..200).map(|k| winning_sample(k, "elliptic-low")));
    for k in 0..16u64 {
        let sample = if k % 4 == 0 {
            winning_sample(1_000_000 + k, "elliptic-high")
        } else {
            losing_sample(1_000_000 + k, "elliptic-high")
        };
        samples.push(sample);
    }
    let report = adjudicate(&frame, &samples).expect("structurally sound evidence");
    assert!(
        !report.activated,
        "an adverse stratum must block activation: {report:?}"
    );
    let adverse = report
        .strata
        .iter()
        .find(|s| s.stratum == "elliptic-high")
        .expect("adverse stratum reported");
    assert!(!adverse.accept_rejected || !adverse.savings_rejected);
    let favorable = report
        .strata
        .iter()
        .find(|s| s.stratum == "elliptic-low")
        .expect("favorable stratum reported");
    assert!(favorable.accept_rejected && favorable.savings_rejected);
}

#[test]
fn optional_stopping_cannot_manufacture_or_lose_a_rejection() {
    // The e-process rejection is a running maximum: once the wealth
    // crosses 1/δ at some prefix, later adverse samples cannot undo the
    // rejection — and conversely, weak evidence never crosses no matter
    // where you stop.
    let frame = frame(12);
    let mut samples = strong_evidence(48);
    // Append a long adverse tail to both strata AFTER the strong runs.
    for k in 0..64u64 {
        samples.push(losing_sample(5_000_000 + k, "elliptic-low"));
        samples.push(losing_sample(6_000_000 + k, "elliptic-high"));
    }
    // Keys sort after the strong evidence, so the canonical path sees
    // the strong prefix first; the tail must not un-reject it.
    let report = adjudicate(&frame, &samples).expect("structurally sound evidence");
    assert!(
        report.activated,
        "a crossed e-threshold is retained under optional continuation: {report:?}"
    );

    // Weak evidence at exactly the minimum count stays un-activated —
    // small-n flukes cannot pass the anytime bound.
    let weak: Vec<Sample> = (0..12)
        .map(|k| {
            if k % 2 == 0 {
                winning_sample(k, "elliptic-low")
            } else {
                losing_sample(k, "elliptic-low")
            }
        })
        .chain((0..12).map(|k| {
            if k % 2 == 0 {
                winning_sample(1_000_000 + k, "elliptic-high")
            } else {
                losing_sample(1_000_000 + k, "elliptic-high")
            }
        }))
        .collect();
    let report = adjudicate(&frame, &weak).expect("structurally sound evidence");
    assert!(
        !report.activated,
        "borderline small-n evidence must not activate: {report:?}"
    );
}

#[test]
fn frame_identity_binds_every_preregistered_field() {
    let base = frame(12);
    let base_id = base.identity();
    let mut seed = base.clone();
    seed.seed += 1;
    let mut floor = base.clone();
    floor.thresholds.savings_floor = 1.6;
    let mut delta = base.clone();
    delta.thresholds.confidence_delta = 0.01;
    let mut strata = base.clone();
    strata.strata[0].min_samples += 1;
    for (what, variant) in [
        ("seed", seed),
        ("savings floor", floor),
        ("confidence", delta),
        ("stratum minimum", strata),
    ] {
        assert_ne!(
            variant.identity(),
            base_id,
            "{what} must be identity-bearing"
        );
    }

    // Malformed frames refuse before any evidence is read.
    let mut bad = frame(12);
    bad.thresholds.confidence_delta = 0.0;
    assert!(matches!(
        adjudicate(&bad, &[]),
        Err(ActivationRefusal::MalformedThresholds)
    ));
    let mut empty = frame(12);
    empty.strata.clear();
    assert!(matches!(
        adjudicate(&empty, &[]),
        Err(ActivationRefusal::MalformedFrame)
    ));
    let unknown = [winning_sample(7, "unregistered-stratum")];
    assert!(matches!(
        adjudicate(&frame(1), &unknown),
        Err(ActivationRefusal::UnknownStratum { key: 7 })
    ));
}
