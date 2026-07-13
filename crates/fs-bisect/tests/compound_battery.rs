//! Failure-compounding acceptance battery (bead 6nb.9): seed the workflow
//! with (a) a deliberately broken cross-crate golden modeled EXACTLY on the
//! real powi incident (bead 4xnt) and (b) a falsifier hit on a wrong
//! certificate constant; both must produce minimized replayable cases,
//! neighborhood boundary evidence, permanent regression families with
//! tracking references, and a content-addressed manifest whose hash is
//! frozen as a golden (identical in both build modes on the recorded aarch64
//! host; cross-ISA reproduction remains pending).

use fs_bisect::compound::{
    Canon, CanonWriter, CompoundError, FailureCase, FamilyProvenance, InvariantClass,
    MAX_CANONICAL_MEMBER_BYTES, MAX_IDENTIFIER_BYTES, MAX_MINIMIZE_EVALUATIONS, MAX_MINIMIZE_STEPS,
    MAX_NEIGHBOR_PROBES, MAX_SHRINK_CANDIDATES_PER_STEP, RegressionFamily, Shrink, canonical_bytes,
    compound, minimize, probe_neighborhood,
};

fn modeled_golden_hash(bytes: &[u8]) -> u64 {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    for &byte in bytes {
        acc ^= u64::from(byte);
        acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
    acc
}

fn provenance(seed: u64) -> FamilyProvenance {
    FamilyProvenance::new(
        seed,
        "fs-bisect::compound".to_string(),
        "seeded regression".to_string(),
    )
    .expect("valid provenance")
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TestMember(u64);

impl Shrink for TestMember {
    fn shrink_candidates(&self) -> Vec<Self> {
        Vec::new()
    }
}

impl Canon for TestMember {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test-member";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        self.0.canon(out)
    }
}

fn test_family_with_provenance(
    name: &str,
    invariant: InvariantClass,
    member: u64,
    tracking: Vec<String>,
    admission: Option<String>,
    provenance: &FamilyProvenance,
) -> RegressionFamily<TestMember> {
    compound(
        FailureCase {
            id: name.to_string(),
            seed: provenance.seed(),
            input: TestMember(member),
            invariant,
            contract: provenance.contract().to_string(),
            detail: provenance.detail().to_string(),
        },
        &|_| true,
        &|_| Vec::new(),
        tracking,
        admission,
        1,
    )
    .expect("valid family")
    .family
}

fn test_family(
    name: &str,
    invariant: InvariantClass,
    member: u64,
    tracking: Vec<String>,
    admission: Option<String>,
) -> RegressionFamily<TestMember> {
    test_family_with_provenance(name, invariant, member, tracking, admission, &provenance(7))
}

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
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.sweep";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        self.base.canon(out)?;
        let exps: Vec<i64> = self.exponents.iter().map(|&e| i64::from(e)).collect();
        exps.canon(out)
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
            bytes.extend_from_slice(
                &canonical_bytes(&f(s.base, k)).expect("fixed-width f64 canonical payload"),
            );
        }
        modeled_golden_hash(&bytes)
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
        Some(
            "forbid variable-exponent integer powers in deterministic paths (4xnt lint)"
                .to_string(),
        ),
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
    assert_eq!(report.family.member_count(), 5);
    assert_eq!(report.family.member_label(0), Some("minimized"));
    assert!(
        !report.family.tracking().is_empty(),
        "no paper trail, no family"
    );
    // Replay: every member still fails under the suspect implementation...
    let live = report
        .family
        .replay(&golden_breaks)
        .expect("sealed family replays");
    assert!(live.now_passing.is_empty(), "family must be live: {live:?}");
    // ...and the SAME family goes fully stale once the bug is "fixed"
    // (both chains sequential) — stale detection is the point of replay.
    let fixed = |_: &Sweep| false;
    let stale = report.family.replay(&fixed).expect("sealed family replays");
    assert!(stale.still_failing.is_empty());
    assert_eq!(stale.now_passing.len(), 5);
}

/// Recorded on aarch64-apple (M4 Pro); reproduced in debug and release.
/// Cross-ISA reproduction is pending and is not claimed by this fixture.
const POWI_FAMILY_CONTENT_HASH: &str =
    "6aed2aab4250ca30b657e6e016f628eb5ba33e09c3e49b732156a5058eb9141f";
const POWI_FAMILY_MANIFEST_HASH: &str =
    "43bd1ebddb606eb5a156ca06c642cf03ddd06770223584c5ea45ae031d0cf6b2";

fn assert_unique_json_object(line: &str) {
    use serde::Deserializer as _;
    use serde::de::{Error as _, MapAccess, Visitor};
    use std::collections::BTreeSet;

    struct UniqueKeys;

    impl<'de> Visitor<'de> for UniqueKeys {
        type Value = ();

        fn expecting(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            formatter.write_str("a JSON object with unique keys")
        }

        fn visit_map<M: MapAccess<'de>>(self, mut map: M) -> Result<(), M::Error> {
            let mut keys = BTreeSet::new();
            while let Some(key) = map.next_key::<String>()? {
                if !keys.insert(key.clone()) {
                    return Err(M::Error::custom(format!("duplicate key {key:?}")));
                }
                let _: serde_json::Value = map.next_value()?;
            }
            Ok(())
        }
    }

    let mut parser = serde_json::Deserializer::from_str(line);
    parser
        .deserialize_map(UniqueKeys)
        .expect("manifest line is valid JSON with unique keys");
    parser.end().expect("manifest line has no trailing syntax");
}

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
    for line in manifest.lines() {
        assert_unique_json_object(line);
    }
    // The manifest carries its member codec and content hash explicitly.
    assert!(manifest.contains("\"family\":\"powi-order-divergence\""));
    assert!(manifest.contains("\"member_type\":\"org.frankensim.fs-bisect.test.sweep\""));
    assert!(manifest.contains("\"member_schema_version\":1"));
    assert_eq!(manifest.lines().count(), 2 + report.family.member_count());
    assert_eq!(
        manifest.lines().next_back(),
        Some(format!("{{\"content_hash\":\"{}\"}}", report.content_hash).as_str())
    );
    let manifest_hash = fs_blake3::hash_bytes(manifest.as_bytes()).to_hex();
    println!(
        "{{\"suite\":\"fs-bisect\",\"case\":\"compound-manifest\",\"verdict\":\"info\",\"content_hash\":\"{}\",\"manifest_hash\":\"{manifest_hash}\"}}",
        report.content_hash,
    );
    assert_eq!(
        report.content_hash.to_hex(),
        POWI_FAMILY_CONTENT_HASH,
        "family bits changed: {} vs {POWI_FAMILY_CONTENT_HASH} — bump only with \
         semantic justification (golden-evidence policy)",
        report.content_hash
    );
    assert_eq!(
        manifest_hash, POWI_FAMILY_MANIFEST_HASH,
        "manifest bytes changed without a deliberate golden re-freeze"
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
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.tail-claim";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        self.n.canon(out)
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
    assert!(report.family.recommended_admission().is_some());
    let live = report
        .family
        .replay(&falsifier_refutes)
        .expect("sealed family replays");
    assert!(live.now_passing.is_empty());
}

// ---- G0 units: determinism, refusal, canon integrity ----

#[test]
fn minimize_is_deterministic_and_refuses_non_failures() {
    let case = powi_case();
    let a = minimize("a", &case.input, &golden_breaks, 1000).expect("fails");
    let b = minimize("b", &case.input, &golden_breaks, 1000).expect("fails");
    let canon = |s: &Sweep| canonical_bytes(s).expect("bounded sweep canon");
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
        let owned: Vec<String> = parts.iter().map(|part| (*part).to_string()).collect();
        fs_blake3::hash_bytes(&canonical_bytes(&owned).expect("bounded strings"))
    };
    assert_ne!(h(&["ab", "c"]), h(&["a", "bc"]));
    let v1 = canonical_bytes(&vec![1u64, 2]).expect("bounded vector");
    let v2 = canonical_bytes(&(vec![1u64], 2u64)).expect("bounded tuple");
    assert_ne!(
        fs_blake3::hash_bytes(&v1),
        fs_blake3::hash_bytes(&v2),
        "length prefixes must separate"
    );
}

#[test]
fn compound_family_name_moves_identity() {
    let left = test_family(
        "family-a",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking".to_string()],
        None,
    );
    let right = test_family(
        "family-b",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking".to_string()],
        None,
    );
    assert_ne!(left.content_hash(), right.content_hash());
}

#[test]
fn compound_family_invariant_moves_identity() {
    let left = test_family(
        "family",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking".to_string()],
        None,
    );
    let right = test_family(
        "family",
        InvariantClass::EnclosureViolation,
        1,
        vec!["tracking".to_string()],
        None,
    );
    assert_ne!(left.content_hash(), right.content_hash());
}

#[test]
fn compound_family_seed_moves_identity() {
    let hash = |seed| {
        test_family_with_provenance(
            "family",
            InvariantClass::GoldenDrift,
            1,
            vec!["tracking".to_string()],
            None,
            &provenance(seed),
        )
        .content_hash()
    };
    assert_ne!(hash(7), hash(8));
}

#[test]
fn compound_family_contract_moves_identity() {
    let hash = |contract: &str| {
        test_family_with_provenance(
            "family",
            InvariantClass::GoldenDrift,
            1,
            vec!["tracking".to_string()],
            None,
            &FamilyProvenance::new(7, contract.to_string(), "seeded regression".to_string())
                .expect("valid provenance"),
        )
        .content_hash()
    };
    assert_ne!(hash("fs-bisect::compound"), hash("fs-bisect::other"));
}

#[test]
fn compound_family_detail_moves_identity() {
    let hash = |detail: &str| {
        test_family_with_provenance(
            "family",
            InvariantClass::GoldenDrift,
            1,
            vec!["tracking".to_string()],
            None,
            &FamilyProvenance::new(7, "fs-bisect::compound".to_string(), detail.to_string())
                .expect("valid provenance"),
        )
        .content_hash()
    };
    assert_ne!(hash("diagnosis-a"), hash("diagnosis-b"));
}

#[test]
fn compound_family_member_labels_move_identity() {
    let family = |label: &str| {
        compound(
            FailureCase {
                id: "family".to_string(),
                seed: 7,
                input: TestMember(1),
                invariant: InvariantClass::GoldenDrift,
                contract: "fs-bisect::compound".to_string(),
                detail: "seeded regression".to_string(),
            },
            &|_| true,
            &|_| vec![(label.to_string(), TestMember(2))],
            vec!["tracking".to_string()],
            None,
            1,
        )
        .expect("valid family")
        .content_hash
    };
    assert_ne!(family("neighbor-a"), family("neighbor-b"));
}

#[test]
fn compound_family_member_bytes_move_identity() {
    let left = test_family(
        "family",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking".to_string()],
        None,
    );
    let right = test_family(
        "family",
        InvariantClass::GoldenDrift,
        2,
        vec!["tracking".to_string()],
        None,
    );
    assert_ne!(left.content_hash(), right.content_hash());
}

#[test]
fn compound_family_tracking_moves_identity() {
    let left = test_family(
        "family",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking-a".to_string()],
        None,
    );
    let right = test_family(
        "family",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking-b".to_string()],
        None,
    );
    assert_ne!(left.content_hash(), right.content_hash());
}

#[test]
fn compound_family_admission_moves_identity() {
    let left = test_family(
        "family",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking".to_string()],
        None,
    );
    let right = test_family(
        "family",
        InvariantClass::GoldenDrift,
        1,
        vec!["tracking".to_string()],
        Some("admit only after replay".to_string()),
    );
    assert_ne!(left.content_hash(), right.content_hash());
}

#[test]
fn content_hash_is_sensitive_to_every_field() {
    let base = test_family(
        "f",
        InvariantClass::GoldenDrift,
        1,
        vec!["t".to_string()],
        None,
    );
    let h0 = base.content_hash();
    assert_ne!(
        h0,
        test_family(
            "g",
            InvariantClass::GoldenDrift,
            1,
            vec!["t".to_string()],
            None,
        )
        .content_hash()
    );
    assert_ne!(
        h0,
        test_family(
            "f",
            InvariantClass::EnclosureViolation,
            1,
            vec!["t".to_string()],
            None,
        )
        .content_hash()
    );
    assert_ne!(
        h0,
        test_family(
            "f",
            InvariantClass::GoldenDrift,
            2,
            vec!["t".to_string()],
            None,
        )
        .content_hash()
    );
    assert_ne!(
        h0,
        test_family(
            "f",
            InvariantClass::GoldenDrift,
            1,
            vec!["t".to_string(), "u".to_string()],
            None,
        )
        .content_hash()
    );
    assert_ne!(
        h0,
        test_family(
            "f",
            InvariantClass::GoldenDrift,
            1,
            vec!["t".to_string()],
            Some("rule".to_string()),
        )
        .content_hash()
    );
}

#[test]
fn content_hash_is_sensitive_to_provenance_fields() {
    let hash = |provenance| {
        test_family_with_provenance(
            "f",
            InvariantClass::GoldenDrift,
            1,
            vec!["t".to_string()],
            None,
            &provenance,
        )
        .content_hash()
    };
    let h0 = hash(provenance(7));
    assert_ne!(h0, hash(provenance(8)), "seed is semantic");
    assert_ne!(
        h0,
        hash(
            FamilyProvenance::new(
                7,
                "fs-bisect::other-contract".to_string(),
                "seeded regression".to_string(),
            )
            .unwrap()
        ),
        "contract is semantic"
    );
    assert_ne!(
        h0,
        hash(
            FamilyProvenance::new(
                7,
                "fs-bisect::compound".to_string(),
                "different diagnosis".to_string(),
            )
            .unwrap()
        ),
        "detail is semantic"
    );
}

#[derive(Clone)]
struct WideShrink;

impl Shrink for WideShrink {
    fn shrink_candidates(&self) -> Vec<Self> {
        vec![Self; MAX_SHRINK_CANDIDATES_PER_STEP + 1]
    }
}

impl Canon for WideShrink {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.wide-shrink";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        1u64.canon(out)
    }
}

#[derive(Clone)]
struct HugeCanon;

impl Canon for HugeCanon {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.huge-canon";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.repeat(0, MAX_CANONICAL_MEMBER_BYTES + 1)
    }
}

impl Shrink for HugeCanon {
    fn shrink_candidates(&self) -> Vec<Self> {
        Vec::new()
    }
}

#[derive(Clone)]
struct EmptyCanon;

impl Canon for EmptyCanon {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.empty-canon";

    fn canon(&self, _out: &mut CanonWriter) -> Result<(), CompoundError> {
        Ok(())
    }
}

impl Shrink for EmptyCanon {
    fn shrink_candidates(&self) -> Vec<Self> {
        Vec::new()
    }
}

#[test]
fn work_envelopes_refuse_at_limit_plus_one() {
    assert!(matches!(
        minimize("wide", &WideShrink, &|_| true, 1),
        Err(CompoundError::LimitExceeded {
            resource: "shrink_candidates_per_step",
            requested,
            max: MAX_SHRINK_CANDIDATES_PER_STEP,
        }) if requested == MAX_SHRINK_CANDIDATES_PER_STEP + 1
    ));
    assert!(matches!(
        minimize(
            "steps",
            &WideShrink,
            &|_| true,
            MAX_MINIMIZE_STEPS + 1,
        ),
        Err(CompoundError::LimitExceeded {
            resource: "minimize_steps",
            requested,
            max: MAX_MINIMIZE_STEPS,
        }) if requested == MAX_MINIMIZE_STEPS + 1
    ));
    let oversized_id = "x".repeat(MAX_IDENTIFIER_BYTES + 1);
    assert!(matches!(
        minimize(&oversized_id, &WideShrink, &|_| true, 0),
        Err(CompoundError::LimitExceeded {
            resource: "case_id",
            requested,
            max: MAX_IDENTIFIER_BYTES,
        }) if requested == MAX_IDENTIFIER_BYTES + 1
    ));
}

#[derive(Clone)]
struct EvaluationBudgetInput {
    generation: usize,
    passing_candidate: bool,
}

impl Shrink for EvaluationBudgetInput {
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut candidates = vec![
            Self {
                generation: self.generation + 1,
                passing_candidate: true,
            };
            MAX_SHRINK_CANDIDATES_PER_STEP - 1
        ];
        candidates.push(Self {
            generation: self.generation + 1,
            passing_candidate: false,
        });
        candidates
    }
}

impl Canon for EvaluationBudgetInput {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.evaluation-budget-input";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        u64::try_from(self.generation)
            .map_err(|_| CompoundError::LimitExceeded {
                resource: "test_generation",
                requested: usize::MAX,
                max: u64::MAX as usize,
            })?
            .canon(out)?;
        self.passing_candidate.canon(out)
    }
}

#[derive(Debug, Clone)]
struct OneStep(u8);

impl Shrink for OneStep {
    fn shrink_candidates(&self) -> Vec<Self> {
        (self.0 > 0).then(|| Self(self.0 - 1)).into_iter().collect()
    }
}

impl Canon for OneStep {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.one-step";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(self.0)
    }
}

#[test]
fn aggregate_evaluation_ceiling_counts_the_seed_call() {
    use std::cell::Cell;

    let calls = Cell::new(0usize);
    let result = minimize(
        "evaluation-budget",
        &EvaluationBudgetInput {
            generation: 0,
            passing_candidate: false,
        },
        &|candidate| {
            calls.set(calls.get() + 1);
            !candidate.passing_candidate
        },
        MAX_MINIMIZE_STEPS,
    );
    assert!(matches!(
        result,
        Err(CompoundError::LimitExceeded {
            resource: "minimize_evaluations",
            requested,
            max: MAX_MINIMIZE_EVALUATIONS,
        }) if requested == MAX_MINIMIZE_EVALUATIONS + 1
    ));
    assert_eq!(
        calls.get(),
        MAX_MINIMIZE_EVALUATIONS,
        "the hard ceiling covers every predicate evaluation"
    );
}

#[test]
fn permanent_family_refuses_budget_limited_minimization() {
    use std::cell::Cell;

    let neighbor_calls = Cell::new(0usize);
    let result = compound(
        FailureCase {
            id: "budget-limited".to_string(),
            seed: 0,
            input: OneStep(1),
            invariant: InvariantClass::GoldenDrift,
            contract: "fs-bisect::compound".to_string(),
            detail: "zero-step minimization budget".to_string(),
        },
        &|_| true,
        &|_| {
            neighbor_calls.set(neighbor_calls.get() + 1);
            Vec::new()
        },
        vec!["frankensim-j3q2".to_string()],
        None,
        0,
    );
    assert_eq!(
        result.unwrap_err(),
        CompoundError::MinimizationIncomplete {
            id: "budget-limited".to_string(),
            steps: 0,
            evaluations: 2,
        }
    );
    assert_eq!(neighbor_calls.get(), 0, "incomplete evidence cannot probe");
}

#[test]
fn exact_accepted_step_budget_still_checks_the_reached_fixpoint() {
    let report = minimize("exact-step-budget", &OneStep(1), &|_| true, 1)
        .expect("one accepted step reaches a fixpoint");
    assert!(report.converged);
    assert_eq!(report.steps, 1);
    assert_eq!(report.tried, 2);
    assert_eq!(report.minimized.0, 0);

    let already_minimal = minimize("zero-step-fixpoint", &OneStep(0), &|_| true, 0)
        .expect("zero accepted steps suffice for an existing fixpoint");
    assert!(already_minimal.converged);
    assert_eq!(already_minimal.tried, 1);
}

#[test]
fn codec_type_domain_separates_identical_payload_bytes() {
    #[derive(Clone)]
    struct CodecA(u64);
    #[derive(Clone)]
    struct CodecB(u64);
    #[derive(Clone)]
    struct CodecV2(u64);

    impl Shrink for CodecA {
        fn shrink_candidates(&self) -> Vec<Self> {
            Vec::new()
        }
    }
    impl Shrink for CodecB {
        fn shrink_candidates(&self) -> Vec<Self> {
            Vec::new()
        }
    }
    impl Shrink for CodecV2 {
        fn shrink_candidates(&self) -> Vec<Self> {
            Vec::new()
        }
    }
    impl Canon for CodecA {
        const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.codec-a";

        fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
            self.0.canon(out)
        }
    }
    impl Canon for CodecB {
        const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.codec-b";

        fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
            self.0.canon(out)
        }
    }
    impl Canon for CodecV2 {
        const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.codec-a";
        const SCHEMA_VERSION: u32 = 2;

        fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
            self.0.canon(out)
        }
    }

    fn family_hash<I: Shrink + Canon>(input: I) -> fs_blake3::ContentHash {
        compound(
            FailureCase {
                id: "codec-domain".to_string(),
                seed: 0,
                input,
                invariant: InvariantClass::GoldenDrift,
                contract: "fs-bisect::compound".to_string(),
                detail: "identical payload, distinct codec".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            vec!["frankensim-j3q2".to_string()],
            None,
            1,
        )
        .expect("codec-domain family")
        .content_hash
    }

    assert_eq!(
        canonical_bytes(&CodecA(7)).unwrap(),
        canonical_bytes(&CodecB(7)).unwrap(),
        "the regression requires identical payload bytes"
    );
    let a = family_hash(CodecA(7));
    assert_ne!(a, family_hash(CodecB(7)), "codec id is semantic");
    assert_ne!(a, family_hash(CodecV2(7)), "codec version is semantic");
}

#[test]
fn replay_refuses_codec_schema_drift_before_predicate_work() {
    use std::cell::Cell;
    use std::sync::atomic::{AtomicU8, Ordering};

    static CHILD_SCHEMA: AtomicU8 = AtomicU8::new(1);

    #[derive(Clone)]
    struct MutableSchema;

    impl Shrink for MutableSchema {
        fn shrink_candidates(&self) -> Vec<Self> {
            Vec::new()
        }
    }
    impl Canon for MutableSchema {
        const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.mutable-schema";

        fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
            7u64.canon(out)
        }

        fn canon_child_schemas(out: &mut CanonWriter) -> Result<(), CompoundError> {
            out.push(CHILD_SCHEMA.load(Ordering::SeqCst))
        }
    }

    CHILD_SCHEMA.store(1, Ordering::SeqCst);
    let family = compound(
        FailureCase {
            id: "schema-drift".to_string(),
            seed: 0,
            input: MutableSchema,
            invariant: InvariantClass::GoldenDrift,
            contract: "fs-bisect::compound".to_string(),
            detail: "stateful codec schema".to_string(),
        },
        &|_| true,
        &|_| Vec::new(),
        vec!["frankensim-j3q2".to_string()],
        None,
        1,
    )
    .expect("initial schema is admitted")
    .family;
    CHILD_SCHEMA.store(2, Ordering::SeqCst);
    let predicate_calls = Cell::new(0usize);
    assert_eq!(
        family
            .replay(&|_| {
                predicate_calls.set(predicate_calls.get() + 1);
                true
            })
            .unwrap_err(),
        CompoundError::ReplayIdentityDrift { member: None }
    );
    assert_eq!(predicate_calls.get(), 0);
    CHILD_SCHEMA.store(1, Ordering::SeqCst);
}

#[test]
fn neighborhood_count_and_labels_are_bounded_and_unambiguous() {
    use std::cell::Cell;

    let at_limit: Vec<(String, u64)> = (0..MAX_NEIGHBOR_PROBES)
        .map(|index| (format!("n-{index}"), index as u64))
        .collect();
    assert_eq!(
        probe_neighborhood(&at_limit, &|_| false)
            .expect("exact neighborhood cap")
            .probes
            .len(),
        MAX_NEIGHBOR_PROBES
    );
    let mut over_limit = at_limit;
    over_limit.push(("over".to_string(), 0));
    assert!(matches!(
        probe_neighborhood(&over_limit, &|_| false),
        Err(CompoundError::LimitExceeded {
            resource: "neighbor_probes",
            requested,
            max: MAX_NEIGHBOR_PROBES,
        }) if requested == MAX_NEIGHBOR_PROBES + 1
    ));
    let duplicate = vec![("same".to_string(), 1), ("same".to_string(), 2)];
    assert!(matches!(
        probe_neighborhood(&duplicate, &|_| true),
        Err(CompoundError::DuplicateIdentity {
            field: "neighbor_label",
            ..
        })
    ));
    let calls = Cell::new(0usize);
    assert!(matches!(
        probe_neighborhood(&[("minimized".to_string(), 1u64)], &|_| {
            calls.set(calls.get() + 1);
            true
        }),
        Err(CompoundError::DuplicateIdentity {
            field: "neighbor_label",
            ..
        })
    ));
    assert_eq!(calls.get(), 0, "reserved labels refuse before work");
}

#[test]
fn family_construction_seals_authority_fields() {
    assert!(matches!(
        compound(
            FailureCase {
                id: "double--separator".to_string(),
                seed: 0,
                input: TestMember(1),
                invariant: InvariantClass::GoldenDrift,
                contract: "fs-bisect::compound".to_string(),
                detail: "ambiguous family separator".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            vec!["frankensim-j3q2".to_string()],
            None,
            1,
        ),
        Err(CompoundError::InvalidField {
            field: "case_id",
            ..
        })
    ));
    assert!(matches!(
        compound(
            FailureCase {
                id: "Not_Kebab".to_string(),
                seed: 0,
                input: TestMember(1),
                invariant: InvariantClass::GoldenDrift,
                contract: "fs-bisect::compound".to_string(),
                detail: "invalid family".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            vec!["frankensim-j3q2".to_string()],
            None,
            1,
        ),
        Err(CompoundError::InvalidField {
            field: "case_id",
            ..
        })
    ));
    assert!(matches!(
        compound(
            FailureCase {
                id: "untracked".to_string(),
                seed: 0,
                input: TestMember(1),
                invariant: InvariantClass::GoldenDrift,
                contract: "fs-bisect::compound".to_string(),
                detail: "untracked family".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            Vec::new(),
            None,
            1,
        ),
        Err(CompoundError::InvalidField {
            field: "tracking",
            ..
        })
    ));
    assert!(matches!(
        compound(
            FailureCase {
                id: "reserved".to_string(),
                seed: 0,
                input: TestMember(1),
                invariant: InvariantClass::Other("golden-drift".to_string()),
                contract: "fs-bisect::compound".to_string(),
                detail: "reserved invariant".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            vec!["t".to_string()],
            None,
            1,
        ),
        Err(CompoundError::InvalidField {
            field: "invariant",
            ..
        })
    ));
}

#[test]
fn family_construction_bounds_canonical_payloads() {
    assert!(matches!(
        compound(
            FailureCase {
                id: "huge".to_string(),
                seed: 0,
                input: HugeCanon,
                invariant: InvariantClass::Other("custom-bound".to_string()),
                contract: "fs-bisect::compound".to_string(),
                detail: "huge canonical payload".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            vec!["frankensim-j3q2".to_string()],
            None,
            1,
        ),
        Err(CompoundError::LimitExceeded {
            resource: "canonical_member_bytes",
            requested,
            max: MAX_CANONICAL_MEMBER_BYTES,
        }) if requested == MAX_CANONICAL_MEMBER_BYTES + 1
    ));
    assert!(matches!(
        compound(
            FailureCase {
                id: "empty-canon".to_string(),
                seed: 0,
                input: EmptyCanon,
                invariant: InvariantClass::Other("custom-bound".to_string()),
                contract: "fs-bisect::compound".to_string(),
                detail: "empty canonical payload".to_string(),
            },
            &|_| true,
            &|_| Vec::new(),
            vec!["frankensim-j3q2".to_string()],
            None,
            1,
        ),
        Err(CompoundError::InvalidField {
            field: "member_canon",
            ..
        })
    ));
}

#[test]
fn invalid_family_authority_refuses_before_predicate_work() {
    use std::cell::Cell;

    let calls = Cell::new(0usize);
    let neighbor_calls = Cell::new(0usize);
    let result = compound(
        FailureCase {
            id: "preflight".to_string(),
            seed: 0,
            input: WideShrink,
            invariant: InvariantClass::GoldenDrift,
            contract: "fs-bisect::compound".to_string(),
            detail: "seeded failure".to_string(),
        },
        &|_| {
            calls.set(calls.get() + 1);
            true
        },
        &|_| {
            neighbor_calls.set(neighbor_calls.get() + 1);
            Vec::new()
        },
        Vec::new(),
        None,
        0,
    );
    assert!(matches!(
        result,
        Err(CompoundError::InvalidField {
            field: "tracking",
            ..
        })
    ));
    assert_eq!(calls.get(), 0, "authority preflight must precede work");
    assert_eq!(neighbor_calls.get(), 0, "invalid authority cannot probe");
}

#[test]
fn manifest_escapes_fields_and_hashes_the_canonical_snapshot() {
    use std::cell::Cell;
    use std::rc::Rc;

    #[derive(Clone)]
    struct MutableCanon(Rc<Cell<u64>>);

    impl Shrink for MutableCanon {
        fn shrink_candidates(&self) -> Vec<Self> {
            Vec::new()
        }
    }

    impl Canon for MutableCanon {
        const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.mutable-canon";

        fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
            self.0.get().canon(out)
        }
    }

    let live_value = Rc::new(Cell::new(7));
    let mutation_handle = Rc::clone(&live_value);
    let family = compound(
        FailureCase {
            id: "family-escape".to_string(),
            seed: 11,
            input: MutableCanon(live_value),
            invariant: InvariantClass::Other("custom-\"\\".to_string()),
            contract: "contract \"quoted\"".to_string(),
            detail: "detail line one\nline two".to_string(),
        },
        &|_| true,
        &|_| {
            vec![(
                "member-\"\\".to_string(),
                MutableCanon(Rc::new(Cell::new(9))),
            )]
        },
        vec!["bead-\"\\".to_string()],
        Some("line one\n\"line two\" \\".to_string()),
        1,
    )
    .expect("escapable visible identifiers")
    .family;
    let before = family.content_hash();
    mutation_handle.set(8);
    assert_eq!(
        family.content_hash(),
        before,
        "content identity must use the sealed construction-time canonical bytes"
    );
    let predicate_calls = Cell::new(0usize);
    assert_eq!(
        family
            .replay(&|_| {
                predicate_calls.set(predicate_calls.get() + 1);
                true
            })
            .unwrap_err(),
        CompoundError::ReplayIdentityDrift {
            member: Some("minimized".to_string())
        }
    );
    assert_eq!(
        predicate_calls.get(),
        0,
        "identity preflight must precede predicate work"
    );
    let manifest = family.manifest();
    let header = manifest.lines().next().expect("header");
    assert!(header.contains("\\\""), "quotes escaped: {header}");
    assert!(header.contains("\\\\"), "backslashes escaped: {header}");
    assert!(header.contains("\\n"), "newlines escaped: {header}");
    assert!(!header.contains('\n'), "one JSON object per line");
    assert!(manifest.ends_with('\n'));
}

#[derive(Clone)]
struct CallbackMutable(std::rc::Rc<std::cell::Cell<u64>>);

impl Shrink for CallbackMutable {
    fn shrink_candidates(&self) -> Vec<Self> {
        Vec::new()
    }
}

impl Canon for CallbackMutable {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.callback-mutable";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        self.0.get().canon(out)
    }
}

#[derive(Clone)]
struct MutatingShrink(std::rc::Rc<std::cell::Cell<u64>>);

impl Shrink for MutatingShrink {
    fn shrink_candidates(&self) -> Vec<Self> {
        self.0.set(self.0.get() + 1);
        Vec::new()
    }
}

impl Canon for MutatingShrink {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.mutating-shrink";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        self.0.get().canon(out)
    }
}

#[derive(Clone)]
struct DetachedCandidate {
    rank: u8,
    state: std::rc::Rc<std::cell::Cell<u64>>,
}

impl Shrink for DetachedCandidate {
    fn shrink_candidates(&self) -> Vec<Self> {
        (self.rank > 0)
            .then(|| Self {
                rank: self.rank - 1,
                state: std::rc::Rc::new(std::cell::Cell::new(self.state.get())),
            })
            .into_iter()
            .collect()
    }
}

impl Canon for DetachedCandidate {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.detached-candidate";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(self.rank)?;
        self.state.get().canon(out)
    }
}

#[test]
fn predicates_cannot_mutate_direct_arguments_or_replay_members() {
    use std::cell::Cell;
    use std::rc::Rc;

    let state = Rc::new(Cell::new(1));
    assert!(matches!(
        minimize(
            "mutating-minimizer",
            &CallbackMutable(Rc::clone(&state)),
            &|input| {
                input.0.set(2);
                true
            },
            0,
        ),
        Err(CompoundError::CallbackIdentityDrift {
            phase: "minimize",
            ..
        })
    ));

    state.set(1);
    let family = compound(
        FailureCase {
            id: "mutating-replay".to_string(),
            seed: 0,
            input: CallbackMutable(Rc::clone(&state)),
            invariant: InvariantClass::GoldenDrift,
            contract: "fs-bisect::compound".to_string(),
            detail: "callback-time interior mutation".to_string(),
        },
        &|_| true,
        &|_| Vec::new(),
        vec!["frankensim-j3q2".to_string()],
        None,
        0,
    )
    .expect("stable construction")
    .family;
    assert!(matches!(
        family.replay(&|input| {
            input.0.set(2);
            true
        }),
        Err(CompoundError::ReplayIdentityDrift {
            member: Some(member)
        }) if member == "minimized"
    ));
}

#[test]
fn minimizer_callbacks_cannot_mutate_retained_witnesses() {
    use std::cell::Cell;
    use std::rc::Rc;

    let shrink_state = Rc::new(Cell::new(1));
    assert!(matches!(
        minimize(
            "mutating-shrink-generator",
            &MutatingShrink(Rc::clone(&shrink_state)),
            &|_| true,
            0,
        ),
        Err(CompoundError::CallbackIdentityDrift {
            phase: "shrink_candidates",
            ..
        })
    ));

    let retained_state = Rc::new(Cell::new(7));
    assert!(matches!(
        minimize(
            "mutating-retained-witness",
            &DetachedCandidate {
                rank: 1,
                state: Rc::clone(&retained_state),
            },
            &|input| {
                if input.rank == 0 {
                    retained_state.set(8);
                }
                true
            },
            1,
        ),
        Err(CompoundError::CallbackIdentityDrift {
            phase: "minimize",
            ..
        })
    ));
}

#[test]
fn neighborhood_callbacks_cannot_mutate_other_members() {
    use std::cell::Cell;
    use std::rc::Rc;

    let state = Rc::new(Cell::new(1));
    assert!(matches!(
        compound(
            FailureCase {
                id: "mutating-neighbor-generator".to_string(),
                seed: 0,
                input: CallbackMutable(Rc::clone(&state)),
                invariant: InvariantClass::GoldenDrift,
                contract: "fs-bisect::compound".to_string(),
                detail: "neighbor callback mutation".to_string(),
            },
            &|_| true,
            &|input| {
                input.0.set(2);
                Vec::new()
            },
            vec!["frankensim-j3q2".to_string()],
            None,
            0,
        ),
        Err(CompoundError::CallbackIdentityDrift {
            phase: "neighbors_of",
            ..
        })
    ));

    let first_neighbor = Rc::new(Cell::new(10));
    let second_neighbor = Rc::new(Cell::new(20));
    assert!(matches!(
        probe_neighborhood(
            &[
                (
                    "first".to_string(),
                    CallbackMutable(Rc::clone(&first_neighbor)),
                ),
                (
                    "second".to_string(),
                    CallbackMutable(Rc::clone(&second_neighbor)),
                ),
            ],
            &|input| {
                if input.0.get() == 20 {
                    first_neighbor.set(11);
                }
                true
            },
        ),
        Err(CompoundError::CallbackIdentityDrift {
            phase: "neighborhood",
            identity,
        }) if identity == "first"
    ));
}

#[test]
fn incomplete_canon_is_an_explicit_replay_no_claim() {
    use std::cell::Cell;
    use std::rc::Rc;

    #[derive(Clone)]
    struct IncompleteCanon {
        encoded: u64,
        hidden_semantic: Rc<Cell<bool>>,
    }

    impl Shrink for IncompleteCanon {
        fn shrink_candidates(&self) -> Vec<Self> {
            Vec::new()
        }
    }

    impl Canon for IncompleteCanon {
        const TYPE_ID: &'static str = "org.frankensim.fs-bisect.test.intentionally-incomplete";

        fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
            self.encoded.canon(out)
        }
    }

    let hidden = Rc::new(Cell::new(true));
    let family = compound(
        FailureCase {
            id: "incomplete-codec-no-claim".to_string(),
            seed: 0,
            input: IncompleteCanon {
                encoded: 7,
                hidden_semantic: Rc::clone(&hidden),
            },
            invariant: InvariantClass::Other("codec-trust-boundary".to_string()),
            contract: "fs-bisect::Canon completeness".to_string(),
            detail: "fixture intentionally omits predicate state".to_string(),
        },
        &|input| input.hidden_semantic.get(),
        &|_| Vec::new(),
        vec!["frankensim-j3q2".to_string()],
        None,
        0,
    )
    .expect("stable but incomplete codec is caller-owned")
    .family;
    hidden.set(false);
    let replay = family
        .replay(&|input| input.hidden_semantic.get())
        .expect("omitted fields are outside the authentication claim");
    assert_eq!(replay.now_passing, vec!["minimized"]);
}
