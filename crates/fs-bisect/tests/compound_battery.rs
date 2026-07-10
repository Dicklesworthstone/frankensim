//! Failure-compounding acceptance battery (bead 6nb.9): seed the workflow
//! with (a) a deliberately broken cross-crate golden modeled EXACTLY on the
//! real powi incident (bead 4xnt) and (b) a falsifier hit on a wrong
//! certificate constant; both must produce minimized replayable cases,
//! neighborhood boundary evidence, permanent regression families with
//! tracking references, and a content-addressed manifest whose hash is
//! frozen as a golden (identical in both build modes and on both ISAs —
//! integer/`to_bits` arithmetic only).

use fs_bisect::compound::{
    Canon, CompoundError, FailureCase, InvariantClass, RegressionFamily, Shrink, compound, fnv64,
    minimize,
};

// ---- Scenario (a): the powi golden break, faithfully modeled ----
//
// Two integer-power implementations with different rounding orders:
// sequential (one rounding per multiply) vs square-and-multiply. They agree
// bitwise through exponent 3 and diverge from 4 — the exact mechanism that
// made the rand_nla golden build-mode-dependent.

fn pow_sequential(x: f64, k: u32) -> f64 {
    let mut p = 1.0f64;
    for _ in 0..k {
        p *= x;
    }
    p
}

fn pow_squaring(x: f64, k: u32) -> f64 {
    let mut b = k;
    let mut a = x;
    let mut r = 1.0f64;
    loop {
        if b & 1 == 1 {
            r *= a;
        }
        b /= 2;
        if b == 0 {
            break;
        }
        a *= a;
    }
    r
}

/// A golden fixture: a sweep of exponents whose combined bits are hashed.
#[derive(Debug, Clone, PartialEq)]
struct Sweep {
    base: f64,
    exponents: Vec<u32>,
}

impl Canon for Sweep {
    fn canon(&self, out: &mut Vec<u8>) {
        self.base.canon(out);
        let exps: Vec<i64> = self.exponents.iter().map(|&e| i64::from(e)).collect();
        exps.canon(out);
    }
}

impl Shrink for Sweep {
    /// Drop-one (leftmost first), then decrement-one — so minimization
    /// walks down to the exact divergence boundary, not just any witness.
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut out = Vec::new();
        if self.exponents.len() > 1 {
            for i in 0..self.exponents.len() {
                let mut exps = self.exponents.clone();
                exps.remove(i);
                out.push(Sweep {
                    base: self.base,
                    exponents: exps,
                });
            }
        }
        for i in 0..self.exponents.len() {
            if self.exponents[i] > 0 {
                let mut exps = self.exponents.clone();
                exps[i] -= 1;
                out.push(Sweep {
                    base: self.base,
                    exponents: exps,
                });
            }
        }
        out
    }
}

/// The "golden drifted" predicate: reference-chain bits != suspect-chain bits.
fn golden_breaks(s: &Sweep) -> bool {
    let feed = |f: &dyn Fn(f64, u32) -> f64| -> u64 {
        let mut bytes = Vec::new();
        for &k in &s.exponents {
            f(s.base, k).canon(&mut bytes);
        }
        fnv64(&bytes)
    };
    feed(&pow_sequential) != feed(&pow_squaring)
}

fn powi_case() -> FailureCase<Sweep> {
    FailureCase {
        id: "powi-order-divergence".to_string(),
        seed: 0xFEED,
        input: Sweep {
            base: 0.7,
            exponents: vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
        },
        invariant: InvariantClass::BuildModeDeterminism,
        contract: "fs-la::rand_nla golden (modeled)".to_string(),
        detail: "sequential vs square-multiply power chains disagree".to_string(),
    }
}

#[test]
fn powi_model_minimizes_to_the_exact_boundary() {
    // Precondition sanity: for base 0.7 the two orders agree through 6 and
    // diverge from 7 onward (an IEEE fact of the two chains — the same
    // mechanism as the real incident, whose lowerings diverged from 4).
    for k in 0..=6u32 {
        assert_eq!(
            pow_sequential(0.7, k).to_bits(),
            pow_squaring(0.7, k).to_bits(),
            "orders must agree at k={k}"
        );
    }
    assert_ne!(
        pow_sequential(0.7, 7).to_bits(),
        pow_squaring(0.7, 7).to_bits(),
        "orders must diverge at k=7"
    );

    let report = compound(
        powi_case(),
        &golden_breaks,
        &|min: &Sweep| {
            (1..=10u32)
                .map(|k| {
                    (
                        format!("k={k}"),
                        Sweep {
                            base: min.base,
                            exponents: vec![k],
                        },
                    )
                })
                .collect()
        },
        vec![
            "frankensim-epic-gauntlet-6nb.9".to_string(),
            "frankensim-powi-build-mode-determinism-4xnt".to_string(),
        ],
        Some("forbid variable-exponent f64::powi in deterministic paths (check-powi)".to_string()),
        1000,
    )
    .expect("the seeded golden break must minimize");

    // Minimized to the EXACT boundary: one exponent, value 7.
    assert!(report.converged, "minimization must reach a fixpoint");
    assert_eq!(report.case.input.exponents, vec![7], "boundary is k=7");
    // Neighborhood: 1..=6 pass, 7..=10 fail — region evidence, sharp edge.
    let failing: Vec<&str> = report
        .neighborhood
        .probes
        .iter()
        .filter(|p| p.fails)
        .map(|p| p.label.as_str())
        .collect();
    assert_eq!(failing, ["k=7", "k=8", "k=9", "k=10"]);
    // The family: minimum first, then every failing neighbor; tracked.
    assert_eq!(report.family.members.len(), 5);
    assert_eq!(report.family.members[0].0, "minimized");
    assert!(
        !report.family.tracking.is_empty(),
        "no paper trail, no family"
    );
    // Replay: every member still fails under the suspect implementation...
    let live = report.family.replay(&golden_breaks);
    assert!(live.now_passing.is_empty(), "family must be live: {live:?}");
    // ...and the SAME family goes fully stale once the bug is "fixed"
    // (both chains sequential) — stale detection is the point of replay.
    let fixed = |_: &Sweep| false;
    let stale = report.family.replay(&fixed);
    assert!(stale.still_failing.is_empty());
    assert_eq!(stale.now_passing.len(), 5);
}

/// Recorded on aarch64-apple (M4 Pro); must be identical in debug and
/// release and on x86-64 (integer/to_bits arithmetic only).
const POWI_FAMILY_MANIFEST_HASH: u64 = 0x9b2d_3f23_3704_8523;

#[test]
fn powi_family_manifest_is_content_addressed_and_frozen() {
    let report = compound(
        powi_case(),
        &golden_breaks,
        &|min: &Sweep| {
            (1..=8u32)
                .map(|k| {
                    (
                        format!("k={k}"),
                        Sweep {
                            base: min.base,
                            exponents: vec![k],
                        },
                    )
                })
                .collect()
        },
        vec!["frankensim-epic-gauntlet-6nb.9".to_string()],
        None,
        1000,
    )
    .expect("must minimize");
    let manifest = report.family.manifest();
    // The manifest carries the hash in its trailer and is replay-complete.
    assert!(manifest.contains("\"family\":\"powi-order-divergence\""));
    assert_eq!(manifest.lines().count(), 2 + report.family.members.len());
    println!(
        "{{\"suite\":\"fs-bisect\",\"case\":\"compound-manifest\",\"verdict\":\"info\",\"detail\":\"{:#018x}\"}}",
        report.content_hash
    );
    assert_eq!(
        report.content_hash, POWI_FAMILY_MANIFEST_HASH,
        "family bits changed: {:#018x} vs {POWI_FAMILY_MANIFEST_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)",
        report.content_hash
    );
}

// ---- Scenario (b): a falsifier hit on a wrong certificate constant ----
//
// A toy certificate claims |Σ_{k≤n} 1/k² − π²/6| ≤ 1/(2n). The true tail
// is ~1/n, so the claim is wrong for every n — a systematic constant
// error, exactly what an independent falsifier exists to catch.

#[derive(Debug, Clone, PartialEq)]
struct TailClaim {
    n: u64,
}

impl Canon for TailClaim {
    fn canon(&self, out: &mut Vec<u8>) {
        self.n.canon(out);
    }
}

impl Shrink for TailClaim {
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut out = Vec::new();
        if self.n > 1 {
            out.push(TailClaim { n: self.n / 2 });
            out.push(TailClaim { n: self.n - 1 });
        }
        out
    }
}

fn falsifier_refutes(c: &TailClaim) -> bool {
    let n = c.n;
    let mut s = 0.0f64;
    for k in 1..=n {
        let kf = k as f64;
        s += 1.0 / (kf * kf);
    }
    let truth = std::f64::consts::PI * std::f64::consts::PI / 6.0;
    let claimed_bound = 1.0 / (2.0 * n as f64);
    (s - truth).abs() > claimed_bound
}

#[test]
fn falsifier_hit_compounds_into_a_family() {
    let report = compound(
        FailureCase {
            id: "basel-tail-constant".to_string(),
            seed: 0,
            input: TailClaim { n: 4096 },
            invariant: InvariantClass::CertificateForgery,
            contract: "toy tail certificate (modeled)".to_string(),
            detail: "claimed 1/(2n) tail bound; true tail ~ 1/n".to_string(),
        },
        &falsifier_refutes,
        &|_min: &TailClaim| {
            [1u64, 2, 8, 64, 1024]
                .iter()
                .map(|&n| (format!("n={n}"), TailClaim { n }))
                .collect()
        },
        vec!["frankensim-epic-gauntlet-6nb.9".to_string()],
        Some("bound constants need a proof or a falsifier-passing margin".to_string()),
        100,
    )
    .expect("the falsifier hit must minimize");
    // Systematic error ⇒ minimizes all the way down to n = 1 and the whole
    // neighborhood fails: region evidence, not a point.
    assert!(report.converged);
    assert_eq!(report.case.input, TailClaim { n: 1 });
    assert_eq!(
        report.neighborhood.failing,
        report.neighborhood.probes.len(),
        "systematic constant error: every probe must fail"
    );
    assert!(report.family.recommended_admission.is_some());
    let live = report.family.replay(&falsifier_refutes);
    assert!(live.now_passing.is_empty());
}

// ---- G0 units: determinism, refusal, canon integrity ----

#[test]
fn minimize_is_deterministic_and_refuses_non_failures() {
    let case = powi_case();
    let a = minimize("a", &case.input, &golden_breaks, 1000).expect("fails");
    let b = minimize("b", &case.input, &golden_breaks, 1000).expect("fails");
    let canon = |s: &Sweep| {
        let mut v = Vec::new();
        s.canon(&mut v);
        v
    };
    assert_eq!(
        canon(&a.minimized),
        canon(&b.minimized),
        "bitwise-identical minimum"
    );
    assert_eq!(a.steps, b.steps);
    assert_eq!(a.tried, b.tried);
    // A passing input is a typed refusal, never a fake minimum.
    let passing = Sweep {
        base: 0.7,
        exponents: vec![1, 2, 3],
    };
    assert_eq!(
        minimize("p", &passing, &golden_breaks, 10).unwrap_err(),
        CompoundError::NotFailing {
            id: "p".to_string()
        }
    );
}

#[test]
fn canon_encoding_resists_concatenation_collisions() {
    let h = |parts: &[&str]| {
        let mut v = Vec::new();
        for p in parts {
            p.canon(&mut v);
        }
        fnv64(&v)
    };
    assert_ne!(h(&["ab", "c"]), h(&["a", "bc"]));
    let mut v1 = Vec::new();
    vec![1u64, 2].canon(&mut v1);
    let mut v2 = Vec::new();
    vec![1u64].canon(&mut v2);
    2u64.canon(&mut v2);
    assert_ne!(fnv64(&v1), fnv64(&v2), "length prefixes must separate");
}

#[test]
fn content_hash_is_sensitive_to_every_field() {
    let base = RegressionFamily {
        name: "f".to_string(),
        invariant: InvariantClass::GoldenDrift,
        members: vec![("m".to_string(), 1u64)],
        tracking: vec!["t".to_string()],
        recommended_admission: None,
    };
    let h0 = base.content_hash();
    let mut renamed = base.clone();
    renamed.name = "g".to_string();
    assert_ne!(h0, renamed.content_hash());
    let mut reclassed = base.clone();
    reclassed.invariant = InvariantClass::EnclosureViolation;
    assert_ne!(h0, reclassed.content_hash());
    let mut remembered = base.clone();
    remembered.members[0].1 = 2u64;
    assert_ne!(h0, remembered.content_hash());
    let mut retracked = base.clone();
    retracked.tracking.push("u".to_string());
    assert_ne!(h0, retracked.content_hash());
    let mut readmitted = base;
    readmitted.recommended_admission = Some("rule".to_string());
    assert_ne!(h0, readmitted.content_hash());
}
