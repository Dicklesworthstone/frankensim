//! Three-color schema conformance (bead qmao.1; CONTRACT.md addendum):
//! composition totality (G0), regime-exit auto-demotion, the adversarial
//! LAUNDERING gauntlet (G3, a Certifying-the-Certifiers gate), the
//! waiver-in-provenance path, the fs-evidence bridge, and determinism.
//! JSON-line verdicts; seeded cases carry seeds.

use fs_evidence::{
    Color, ColorError, ColorRank, IntervalOp, ModelEvidence, NumericalCertificate, ValidityDomain,
    check_regime, color_of, compose, verified_from,
};
use fs_ledger::{
    ColorGraph, ColorStructureRejection, ColorWriteError, MAX_COLOR_PARENTS,
    NoSourceOriginVerifier, NoWaiverVerifier, PolicyDecision, SourceOrigin, SourceOriginRejection,
    SourceOriginRequest, SourceOriginVerifier, WAIVER_SCOPE_COLOR_UPGRADE,
    WAIVER_SCOPE_SOURCE_COLOR, Waiver, WaiverDependency, WaiverGrant, WaiverRejection,
    WaiverVerifier, hash_bytes,
};
use std::collections::BTreeMap as KeyMap;

/// TEST-ONLY keyed MAC over FNV — NOT cryptography; it stands in for a
/// Franken-compliant signature capability so the authorization
/// plumbing (binding, expiry, rotation, tamper) can be exercised.
struct TestVerifier {
    keys: KeyMap<String, u64>,
}

fn test_mac(secret: u64, payload: &[u8]) -> Vec<u8> {
    let mut acc = 0xcbf2_9ce4_8422_2325u64 ^ secret;
    for &b in payload {
        acc ^= u64::from(b);
        acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
    acc.to_le_bytes().to_vec()
}

impl WaiverVerifier for TestVerifier {
    fn verify(&self, key_id: &str, payload: &[u8], signature: &[u8]) -> PolicyDecision {
        let accepted = self
            .keys
            .get(key_id)
            .is_some_and(|&secret| test_mac(secret, payload) == signature);
        let mut policy = b"fs-ledger/color-battery/waiver-policy/v1".to_vec();
        for (key, secret) in &self.keys {
            let key_len = u64::try_from(key.len()).expect("test key length fits u64");
            policy.extend_from_slice(&key_len.to_le_bytes());
            policy.extend_from_slice(key.as_bytes());
            policy.extend_from_slice(&secret.to_le_bytes());
        }
        let fingerprint = hash_bytes(&policy);
        if accepted {
            PolicyDecision::accept(fingerprint)
        } else {
            PolicyDecision::reject(fingerprint)
        }
    }
}

#[derive(Clone)]
struct TestSourceVerifier {
    accepted_request: Vec<u8>,
    policy_fingerprint: fs_ledger::ContentHash,
}

impl TestSourceVerifier {
    fn authorizing(name: &str, color: &Color, origin: &SourceOrigin) -> Self {
        let accepted_request = SourceOriginRequest::new(name, color, origin).canonical_bytes();
        let mut policy = b"fs-ledger/color-battery/source-policy/v1".to_vec();
        policy.extend_from_slice(&accepted_request);
        Self {
            accepted_request,
            policy_fingerprint: hash_bytes(&policy),
        }
    }
}

impl SourceOriginVerifier for TestSourceVerifier {
    fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision {
        if self.accepted_request == request.canonical_bytes() {
            PolicyDecision::accept(self.policy_fingerprint)
        } else {
            PolicyDecision::reject(self.policy_fingerprint)
        }
    }
}

struct PanickingSourceVerifier;

impl SourceOriginVerifier for PanickingSourceVerifier {
    fn verify(&self, _request: &SourceOriginRequest<'_>) -> PolicyDecision {
        panic!("hostile source verifier")
    }
}

struct PanickingWaiverVerifier;

impl WaiverVerifier for PanickingWaiverVerifier {
    fn verify(&self, _key_id: &str, _payload: &[u8], _signature: &[u8]) -> PolicyDecision {
        panic!("hostile waiver verifier")
    }
}

fn admitted_source(
    graph: &mut ColorGraph,
    name: &str,
    color: &Color,
    origin: SourceOrigin,
) -> Result<u64, ColorWriteError> {
    let verifier = TestSourceVerifier::authorizing(name, color, &origin);
    graph.source_with_origin(name, color, origin, &verifier)
}

fn signed_grant(
    secret: u64,
    key_id: &str,
    name: &str,
    color: &Color,
    parent_hashes: Vec<fs_ledger::ContentHash>,
    op: IntervalOp,
    expires_day: u32,
) -> WaiverGrant {
    let mut grant = WaiverGrant {
        annotation: Waiver {
            id: "WVR-2026-041".to_string(),
            signer: "chief-engineer".to_string(),
            reason: "surrogate validated offline against holdout campaign 7".to_string(),
        },
        key_id: key_id.to_string(),
        scope: WAIVER_SCOPE_COLOR_UPGRADE.to_string(),
        node_name: name.to_string(),
        claimed_color: color.canonical_bytes(),
        parent_hashes,
        expires_day,
        signature: Vec::new(),
    };
    grant.signature = test_mac(secret, &grant.signing_payload(op));
    grant
}

fn signed_source_grant(
    secret: u64,
    key_id: &str,
    name: &str,
    color: &Color,
    expires_day: u32,
) -> WaiverGrant {
    let mut grant = WaiverGrant {
        annotation: Waiver {
            id: "WVR-SOURCE-2026-001".to_string(),
            signer: "chief-engineer".to_string(),
            reason: "authenticated exceptional source admission".to_string(),
        },
        key_id: key_id.to_string(),
        scope: WAIVER_SCOPE_SOURCE_COLOR.to_string(),
        node_name: name.to_string(),
        claimed_color: color.canonical_bytes(),
        parent_hashes: Vec::new(),
        expires_day,
        signature: Vec::new(),
    };
    grant.signature = test_mac(secret, &grant.signing_payload_source());
    grant
}

fn structural_refusal(result: Result<u64, ColorWriteError>) -> ColorStructureRejection {
    match result {
        Err(ColorWriteError::InvalidColor { rejection }) => rejection,
        other => panic!("expected structural color refusal, got {other:?}"),
    }
}

fn write_source(graph: &mut ColorGraph, name: &str, color: Color) -> u64 {
    match &color {
        Color::Estimated { .. } => graph.source(name, color).expect("Estimated source"),
        Color::Verified { lo, hi } => admitted_source(
            graph,
            name,
            &color,
            SourceOrigin::Certificate {
                producer: "fs-ledger/color-battery".to_string(),
                certificate_hash: hash_bytes(name.as_bytes()),
                certificate: NumericalCertificate::enclosure(*lo, *hi),
            },
        )
        .expect("Verified source certificate"),
        Color::Validated { regime, dataset } => admitted_source(
            graph,
            name,
            &color,
            SourceOrigin::Anchoring {
                dataset_id: dataset.clone(),
                content_hash: hash_bytes(dataset.as_bytes()),
                regime: regime.clone(),
            },
        )
        .expect("Validated source anchor"),
    }
}

fn test_hex(bytes: &[u8]) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

fn decode_hex(hex: &str) -> Option<Vec<u8>> {
    let (pairs, remainder) = hex.as_bytes().as_chunks::<2>();
    if !remainder.is_empty() {
        return None;
    }
    pairs
        .iter()
        .map(|pair| {
            let digits = core::str::from_utf8(pair).ok()?;
            u8::from_str_radix(digits, 16).ok()
        })
        .collect()
}

fn json_string_field<'a>(row: &'a str, key: &str) -> Option<&'a str> {
    let marker = format!("\"{key}\":\"");
    let value = row.get(row.find(&marker)? + marker.len()..)?;
    value.get(..value.find('"')?)
}
use std::collections::BTreeMap;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ledger/colors\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn below(&mut self, n: u64) -> u64 {
        (self.next() >> 32) % n
    }
}

fn rand_color(rng: &mut Lcg) -> Color {
    match rng.below(3) {
        0 => {
            let lo = rng.unit() * 2.0 - 1.0;
            Color::Verified {
                lo,
                hi: lo + rng.unit(),
            }
        }
        1 => Color::Validated {
            regime: ValidityDomain::unconstrained().with(
                "reynolds",
                1e3 * (1.0 + rng.unit()),
                1e5 * (1.0 + rng.unit()),
            ),
            dataset: format!("ds-{}", rng.below(4)),
        },
        _ => Color::Estimated {
            estimator: format!("est-{}", rng.below(4)),
            dispersion: rng.unit(),
        },
    }
}

/// col-001 — G0 totality: every color pair composes to a defined,
/// CONSERVATIVE result (rank = min of operand ranks, never higher);
/// verified bounds compose per op; validated regimes intersect;
/// estimated dispersions add.
#[test]
fn col_001_composition_totality() {
    let mut rng = Lcg(0x1001_2026_0707_0041);
    let mut total_ok = true;
    let mut conservative_ok = true;
    for _ in 0..600 {
        let (a, b) = (rand_color(&mut rng), rand_color(&mut rng));
        let op = match rng.below(3) {
            0 => IntervalOp::Add,
            1 => IntervalOp::Mul,
            _ => IntervalOp::Hull,
        };
        let out = compose(&a, &b, op);
        conservative_ok &= out.rank() <= a.rank().min(b.rank());
        total_ok &= matches!(
            out,
            Color::Verified { .. } | Color::Validated { .. } | Color::Estimated { .. }
        );
    }
    // Verified interval arithmetic spot checks.
    let v = |lo: f64, hi: f64| Color::Verified { lo, hi };
    let add = compose(&v(1.0, 2.0), &v(10.0, 20.0), IntervalOp::Add);
    let mul = compose(&v(-1.0, 2.0), &v(3.0, 4.0), IntervalOp::Mul);
    let verified_math = matches!(
        add,
        Color::Verified { lo, hi }
            if lo.to_bits() == 11.0_f64.next_down().to_bits()
                && hi.to_bits() == 22.0_f64.next_up().to_bits()
    ) && matches!(
        mul,
        Color::Verified { lo, hi }
            if lo.to_bits() == (-4.0_f64).next_down().to_bits()
                && hi.to_bits() == 8.0_f64.next_up().to_bits()
    );
    // Regime intersection: both anchors must hold.
    let val = |lo: f64, hi: f64| Color::Validated {
        regime: ValidityDomain::unconstrained().with("re", lo, hi),
        dataset: "wind-tunnel-a".to_string(),
    };
    let both = compose(&val(1e3, 1e5), &val(5e3, 5e5), IntervalOp::Add);
    let intersected = matches!(&both, Color::Validated { regime, .. }
        if regime.bounds()["re"] == (5e3, 1e5));
    // Estimated absorbs everything, dispersion conservatively.
    let est = compose(
        &Color::Estimated {
            estimator: "koopman".to_string(),
            dispersion: 0.1,
        },
        &v(0.0, 1.0),
        IntervalOp::Add,
    );
    let absorbed = matches!(est, Color::Estimated { dispersion, .. } if dispersion >= 0.1);
    verdict(
        "col-001",
        total_ok && conservative_ok && verified_math && intersected && absorbed,
        "all 600 random pairs compose totally with rank = min (never higher), \
         verified interval add/mul outward-round to true enclosures, validated \
         regimes INTERSECT, and \
         estimated absorbs everything with conservative dispersion; \
         seed 0x1001_2026_0707_0041",
    );
}

/// col-002 — regime-exit AUTO-DEMOTION: validated survives inside its
/// regime, demotes to estimated (with the flag naming dataset, axis,
/// value) the moment any axis exits or goes unreported.
#[test]
fn col_002_regime_demotion() {
    let validated = Color::Validated {
        regime: ValidityDomain::unconstrained()
            .with("reynolds", 1e3, 1e5)
            .with("mach", 0.0, 0.3),
        dataset: "wind-tunnel-a".to_string(),
    };
    let inside: BTreeMap<String, f64> =
        [("reynolds".to_string(), 5e4), ("mach".to_string(), 0.2)].into();
    let (c1, d1) = check_regime(&validated, &inside);
    let stays = c1 == validated && d1.is_none();

    let outside: BTreeMap<String, f64> =
        [("reynolds".to_string(), 2e5), ("mach".to_string(), 0.2)].into();
    let (c2, d2) = check_regime(&validated, &outside);
    let demotes = matches!(&c2, Color::Estimated { estimator, dispersion }
        if estimator.contains("regime-exit") && estimator.contains("wind-tunnel-a")
            && dispersion.is_infinite())
        && d2
            .as_ref()
            .is_some_and(|d| d.axis == "reynolds" && (d.value - 2e5).abs() < 1.0);

    let unreported: BTreeMap<String, f64> = [("reynolds".to_string(), 5e4)].into();
    let (c3, d3) = check_regime(&validated, &unreported);
    let unreported_demotes =
        matches!(c3, Color::Estimated { .. }) && d3.is_some_and(|d| d.axis == "mach");

    // Verified and estimated pass through untouched.
    let v = Color::Verified { lo: 0.0, hi: 1.0 };
    let (cv, dv) = check_regime(&v, &outside);
    let passthrough = cv == v && dv.is_none();

    verdict(
        "col-002",
        stays && demotes && unreported_demotes && passthrough,
        "validated survives inside its regime, AUTO-DEMOTES to estimated (infinite \
         dispersion, flag naming wind-tunnel-a/reynolds/2e5) on exit, demotes on an \
         UNREPORTED axis, and verified/estimated pass through untouched",
    );
}

/// col-003 — the LAUNDERING gauntlet (G3, security-critical): every
/// adversarial attempt to upgrade a color fails the type check with
/// the capping parents named.
#[test]
#[allow(clippy::too_many_lines)] // the gauntlet's five doors are one story
fn col_003_laundering_gauntlet() {
    // (a) The constructor door: estimate/no-claim certs refuse verified.
    let est_cert = NumericalCertificate::estimate(0.0, 1.0);
    let door = verified_from(&est_cert);
    let door_refuses = matches!(&door, Err(ColorError::LaunderingRefused { actual }) if *actual == "estimate")
        && door.unwrap_err().to_string().contains("waiver");
    let no_claim = verified_from(&NumericalCertificate::no_claim());
    let no_claim_refuses = no_claim.is_err();
    // Enclosures pass the door.
    let ok_door = verified_from(&NumericalCertificate::enclosure(0.0, 0.5)).is_ok();

    // (b) The write gate: estimated parent caps every claim.
    let state = BTreeMap::new();
    let mut g = ColorGraph::new();
    let clean = write_source(&mut g, "fem-bound", Color::Verified { lo: 0.9, hi: 1.1 });
    let dirty = write_source(
        &mut g,
        "surrogate-drag",
        Color::Estimated {
            estimator: "deeponet-v2".to_string(),
            dispersion: 0.05,
        },
    );
    let attempt = g.derive(
        "total-drag",
        &[clean, dirty],
        IntervalOp::Add,
        Some(Color::Verified { lo: 0.0, hi: 2.0 }),
        &state,
        None,
    );
    let gate_refuses = matches!(&attempt,
        Err(ColorWriteError::LaunderingRefused { claimed: ColorRank::Verified,
            derived: ColorRank::Estimated, offending_parents })
            if offending_parents.contains(&dirty));

    // (c) Claiming validated over estimated parents refuses too.
    let attempt2 = g.derive(
        "calibrated-drag",
        &[dirty],
        IntervalOp::Hull,
        Some(Color::Validated {
            regime: ValidityDomain::unconstrained().with("re", 1e3, 1e5),
            dataset: "wishful".to_string(),
        }),
        &state,
        None,
    );
    let gate_refuses2 = attempt2.is_err();

    // (d) Post-demotion upgrade attempts refuse: a validated parent
    // whose regime the state has exited caps at estimated.
    let val = write_source(
        &mut g,
        "turbulence-closure",
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
            dataset: "wind-tunnel-a".to_string(),
        },
    );
    let exited: BTreeMap<String, f64> = [("reynolds".to_string(), 9e5)].into();
    let attempt3 = g.derive(
        "lift-coefficient",
        &[val],
        IntervalOp::Hull,
        Some(Color::Validated {
            regime: ValidityDomain::unconstrained().with("reynolds", 1e3, 1e6),
            dataset: "wind-tunnel-a".to_string(),
        }),
        &exited,
        None,
    );
    let post_demotion_refuses = attempt3.is_err();

    // (e) Seeded adversarial pipelines: random DAGs, random upgrade
    // attempts — every single one must refuse.
    let mut rng = Lcg(0x1001_2026_0707_0043);
    let mut attempts = 0u32;
    let mut refusals = 0u32;
    for _ in 0..60 {
        let mut gg = ColorGraph::new();
        let mut ids = Vec::new();
        for k in 0..5 {
            ids.push(write_source(
                &mut gg,
                &format!("s{k}"),
                rand_color(&mut rng),
            ));
        }
        for k in 0..6 {
            let a = ids[rng.below(ids.len() as u64) as usize];
            let b = ids[rng.below(ids.len() as u64) as usize];
            let derived_rank = {
                let (ca, _) =
                    check_regime(gg.node(a).expect("a").declared_color_unverified(), &state);
                let (cb, _) =
                    check_regime(gg.node(b).expect("b").declared_color_unverified(), &state);
                compose(&ca, &cb, IntervalOp::Hull).rank()
            };
            // Claim strictly ABOVE what the parents support.
            let claim = match derived_rank {
                ColorRank::Verified => None, // nothing to launder
                ColorRank::Validated | ColorRank::Estimated => {
                    Some(Color::Verified { lo: 0.0, hi: 1.0 })
                }
            };
            if let Some(c) = claim {
                attempts += 1;
                if gg
                    .derive(
                        &format!("d{k}"),
                        &[a, b],
                        IntervalOp::Hull,
                        Some(c),
                        &state,
                        None,
                    )
                    .is_err()
                {
                    refusals += 1;
                }
            } else if let Ok(id) = gg.derive(
                &format!("d{k}"),
                &[a, b],
                IntervalOp::Hull,
                None,
                &state,
                None,
            ) {
                ids.push(id);
            }
        }
    }
    verdict(
        "col-003",
        door_refuses
            && no_claim_refuses
            && ok_door
            && gate_refuses
            && gate_refuses2
            && post_demotion_refuses
            && attempts > 100
            && refusals == attempts,
        &format!(
            "every laundering path refuses: the constructor door (estimate/no-claim \
             -> verified), the write gate (estimated parent caps all claims, capping \
             parents NAMED), validated-over-estimated, post-demotion re-claims, and \
             {refusals}/{attempts} seeded adversarial upgrade attempts; \
             seed 0x1001_2026_0707_0043"
        ),
    );
}

/// col-004 — the AUTHENTICATED waiver path (qmao.1.1): only a grant
/// bound to this node, lineage, color, and scope, unexpired, with a
/// signature the verifier accepts, authorizes an upgrade; the grant
/// travels in the provenance hash; annotations alone authorize
/// nothing; tamper/replay/expiry/rotation/no-crypto all fail closed.
#[test]
#[allow(clippy::too_many_lines)] // tamper, replay, expiry, rotation, and no-crypto form one auth story
fn col_004_waiver_in_provenance() {
    let state = BTreeMap::new();
    let verifier = TestVerifier {
        keys: KeyMap::from([("release-key-1".to_string(), 0x5EC2E7u64)]),
    };
    let claimed = Color::Verified { lo: 0.0, hi: 1.0 };
    let fresh = || {
        let mut g = ColorGraph::new();
        let dirty = write_source(
            &mut g,
            "surrogate",
            Color::Estimated {
                estimator: "pod-deim".to_string(),
                dispersion: 0.02,
            },
        );
        (g, dirty)
    };
    // An annotation via derive() authorizes NOTHING.
    let (mut g0, d0) = fresh();
    let annotated = g0.derive(
        "release-metric",
        &[d0],
        IntervalOp::Hull,
        Some(claimed.clone()),
        &state,
        Some(Waiver {
            id: "WVR-2026-041".to_string(),
            signer: "chief-engineer".to_string(),
            reason: "trust me".to_string(),
        }),
    );
    assert!(
        matches!(annotated, Err(ColorWriteError::LaunderingRefused { .. })),
        "caller-created strings cannot authorize promotion"
    );
    // The authenticated grant DOES authorize, deterministically.
    let (mut g, dirty) = fresh();
    let lineage = vec![g.node(dirty).expect("dirty source").hash()];
    let grant = signed_grant(
        0x5EC2E7,
        "release-key-1",
        "release-metric",
        &claimed,
        lineage.clone(),
        IntervalOp::Hull,
        400,
    );
    let id = g
        .derive_waived(
            "release-metric",
            &[dirty],
            IntervalOp::Hull,
            claimed.clone(),
            &state,
            grant.clone(),
            &verifier,
            200,
        )
        .expect("authenticated grant authorizes");
    let node = g.node(id).expect("authorized node").clone();
    assert!(node.grant().is_some() && node.waiver().is_some());
    assert!(
        g.rows()
            .iter()
            .any(|r| r.contains("release-key-1") && r.contains("\"authorized\":true"))
    );
    // Re-verifiable FROM the ledger: the stored grant's payload still
    // authenticates under the same verifier.
    let stored = node.grant().expect("stored grant");
    assert!(
        verifier
            .verify(
                &stored.key_id,
                &stored.signing_payload(IntervalOp::Hull),
                &stored.signature
            )
            .accepted()
    );
    // Provenance: the same write without the grant hashes differently.
    let (mut g2, d2) = fresh();
    let plain = g2
        .derive(
            "release-metric",
            &[d2],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("plain");
    assert_ne!(
        g2.node(plain).expect("plain node").hash().to_hex(),
        node.hash().to_hex()
    );
    // REFUSALS, each with its structured reason:
    let refusal = |grant: WaiverGrant, verifier: &dyn WaiverVerifier, day: u32| {
        let (mut gx, dx) = fresh();
        match gx.derive_waived(
            "release-metric",
            &[dx],
            IntervalOp::Hull,
            claimed.clone(),
            &state,
            grant,
            verifier,
            day,
        ) {
            Err(ColorWriteError::WaiverRefused { rejection }) => rejection,
            other => panic!("expected refusal, got {other:?}"),
        }
    };
    // Tampered payload (signature no longer matches).
    let mut tampered = grant.clone();
    tampered.annotation.reason = "edited after signing".to_string();
    assert!(matches!(
        refusal(tampered, &verifier, 200),
        WaiverRejection::VerifierRefused { .. }
    ));
    // Replay to another node name.
    let mut wrong_node = grant.clone();
    wrong_node.node_name = "other-metric".to_string();
    assert_eq!(
        refusal(wrong_node, &verifier, 200),
        WaiverRejection::NodeMismatch
    );
    // Replay to a different lineage.
    let mut wrong_lineage = grant.clone();
    wrong_lineage.parent_hashes = vec![];
    wrong_lineage.signature = test_mac(0x5EC2E7, &wrong_lineage.signing_payload(IntervalOp::Hull));
    assert_eq!(
        refusal(wrong_lineage, &verifier, 200),
        WaiverRejection::LineageMismatch
    );
    // A signature for one exact interval cannot authorize another interval of
    // the same rank, even when both have the same human-readable color name.
    let different_interval = Color::Verified {
        lo: 0.0,
        hi: 1.0f64.next_up(),
    };
    let (mut exact_graph, exact_parent) = fresh();
    assert!(matches!(
        exact_graph.derive_waived(
            "release-metric",
            &[exact_parent],
            IntervalOp::Hull,
            different_interval,
            &state,
            grant.clone(),
            &verifier,
            200,
        ),
        Err(ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::ColorMismatch
        })
    ));
    // Wrong color, wrong scope, expiry.
    let mut wrong_color = grant.clone();
    wrong_color.claimed_color = Color::Validated {
        regime: ValidityDomain::unconstrained(),
        dataset: "wrong".to_string(),
    }
    .canonical_bytes();
    assert_eq!(
        refusal(wrong_color, &verifier, 200),
        WaiverRejection::ColorMismatch
    );
    let mut wrong_scope = grant.clone();
    wrong_scope.scope = "deploy".to_string();
    assert_eq!(
        refusal(wrong_scope, &verifier, 200),
        WaiverRejection::ScopeMismatch
    );
    assert_eq!(
        refusal(grant.clone(), &verifier, 401),
        WaiverRejection::Expired
    );
    // Wrong key + KEY ROTATION (old key removed from the verifier).
    let mut wrong_key = grant.clone();
    wrong_key.key_id = "release-key-2".to_string();
    assert!(matches!(
        refusal(wrong_key, &verifier, 200),
        WaiverRejection::VerifierRefused { .. }
    ));
    let rotated = TestVerifier {
        keys: KeyMap::from([("release-key-2".to_string(), 0x0FF1CEu64)]),
    };
    assert!(matches!(
        refusal(grant.clone(), &rotated, 200),
        WaiverRejection::VerifierRefused { .. }
    ));
    // Delimiter injection: adversarial names cannot collide the
    // length-prefixed payload of a DIFFERENT legitimate grant.
    let evil = signed_grant(
        0x5EC2E7,
        "release-key-1",
        "release-metric\ninjected",
        &claimed,
        lineage,
        IntervalOp::Hull,
        400,
    );
    assert_ne!(
        evil.signing_payload(IntervalOp::Hull),
        grant.signing_payload(IntervalOp::Hull)
    );
    // No-crypto no-claim: the in-tree default verifier refuses all.
    assert!(matches!(
        refusal(grant, &NoWaiverVerifier, 200),
        WaiverRejection::VerifierRefused { .. }
    ));
    verdict(
        "col-004",
        true,
        "annotations never authorize; the authenticated grant does (and re-verifies \
         from the ledger); tamper, node/lineage replay, color/scope mismatch, expiry, \
         wrong key, rotation, injection, and no-crypto all fail closed",
    );
}

/// col-005 — the fs-evidence bridge: existing receipts color honestly
/// (enclosure -> verified; carded model with bounded validity ->
/// validated with THAT regime; estimates -> estimated).
#[test]
fn col_005_evidence_bridge() {
    let verified = color_of(
        &NumericalCertificate::enclosure(0.9, 1.1),
        &ModelEvidence::none(),
    );
    // Bounds pass through by BITS (no arithmetic on this path).
    let v_ok = matches!(verified, Color::Verified { lo, hi }
        if lo.to_bits() == 0.9f64.to_bits() && hi.to_bits() == 1.1f64.to_bits());
    let validated = color_of(
        &NumericalCertificate::estimate(0.0, 1.0),
        &ModelEvidence {
            cards: vec!["k-epsilon".to_string()],
            assumptions: vec![],
            validity: ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
            discrepancy_rel: 0.03,
            in_domain: true,
        },
    );
    let val_ok = matches!(&validated, Color::Validated { regime, dataset }
        if dataset == "k-epsilon" && regime.bounds().contains_key("reynolds"));
    let estimated = color_of(
        &NumericalCertificate::estimate(0.0, 1.0),
        &ModelEvidence::none(),
    );
    let est_ok = matches!(estimated, Color::Estimated { .. });
    verdict(
        "col-005",
        v_ok && val_ok && est_ok,
        "existing fs-evidence receipts color honestly: enclosures become verified \
         with their bounds, carded models with bounded validity become validated \
         with THAT regime, uncarded estimates stay estimated",
    );
}

/// col-006 — determinism: identical write sequences give bitwise
/// identical rows and hashes.
#[test]
fn col_006_determinism() {
    let build = || -> Vec<String> {
        let state: BTreeMap<String, f64> = [("re".to_string(), 5e4)].into();
        let mut g = ColorGraph::new();
        let a = write_source(&mut g, "a", Color::Verified { lo: 0.0, hi: 1.0 });
        let b = write_source(
            &mut g,
            "b",
            Color::Validated {
                regime: ValidityDomain::unconstrained().with("re", 1e3, 1e5),
                dataset: "ds".to_string(),
            },
        );
        let c = g
            .derive("c", &[a, b], IntervalOp::Add, None, &state, None)
            .expect("c");
        let _ = c;
        g.rows().to_vec()
    };
    let (r1, r2) = (build(), build());
    verdict(
        "col-006",
        r1 == r2 && !r1.is_empty(),
        "identical write sequences produce bitwise-identical rows and provenance \
         hashes",
    );
}

/// col-006b — provenance identity uses the exact Color payload, not its rounded
/// JSON rendering. Sub-render-resolution interval changes must change the hash.
#[test]
fn col_006b_bit_exact_color_identity() {
    let first = Color::Verified { lo: 1.0, hi: 2.0 };
    let second = Color::Verified {
        lo: 1.0f64.next_up(),
        hi: 2.0,
    };
    assert_eq!(
        first.payload_json(),
        second.payload_json(),
        "the regression needs two values hidden by the display precision"
    );
    assert_ne!(first.canonical_bytes(), second.canonical_bytes());

    let mut first_graph = ColorGraph::new();
    let first_id = write_source(&mut first_graph, "same-name", first);
    let mut second_graph = ColorGraph::new();
    let second_id = write_source(&mut second_graph, "same-name", second);
    assert_ne!(
        first_graph.node(first_id).expect("first").hash(),
        second_graph.node(second_id).expect("second").hash(),
        "ledger identity must include the bit-exact color payload"
    );
}

/// col-008 — an ordinary explicit claim must be the exact algebra result;
/// equal rank is not permission to narrow intervals or invent a different
/// lower-rank payload.
#[test]
fn col_008_ordinary_claim_must_equal_the_derived_color() {
    let mut graph = ColorGraph::new();
    let source_color = Color::Verified { lo: 0.0, hi: 1.0 };
    let parent = write_source(&mut graph, "verified-source", source_color.clone());

    let exact = graph.derive(
        "exact-claim",
        &[parent],
        IntervalOp::Hull,
        Some(source_color),
        &BTreeMap::new(),
        None,
    );
    assert!(
        exact.is_ok(),
        "the bit-exact derived color remains admissible"
    );

    let narrowed = graph.derive(
        "narrowed-claim",
        &[parent],
        IntervalOp::Hull,
        Some(Color::Verified { lo: 0.25, hi: 0.75 }),
        &BTreeMap::new(),
        None,
    );
    assert!(matches!(
        narrowed,
        Err(ColorWriteError::ClaimMismatch { .. })
    ));

    let invented_weakening = graph.derive(
        "invented-weaker-claim",
        &[parent],
        IntervalOp::Hull,
        Some(Color::Estimated {
            estimator: "unsupported-relabel".to_string(),
            dispersion: 0.0,
        }),
        &BTreeMap::new(),
        None,
    );
    assert!(matches!(
        invented_weakening,
        Err(ColorWriteError::ClaimMismatch { .. })
    ));

    verdict(
        "col-008",
        true,
        "ordinary derive admits only the exact canonical color; same-rank interval \
         narrowing and unmodeled rank weakening both refuse",
    );
}

/// col-009 — operation identity is part of both node provenance and signed
/// waiver authority. Add and Mul cannot collide or reuse one grant.
#[test]
fn col_009_operation_identity_and_waivers_are_separated() {
    let build = |op| {
        let mut graph = ColorGraph::new();
        let parent = write_source(
            &mut graph,
            "same-parent",
            Color::Verified { lo: 0.0, hi: 1.0 },
        );
        let node = graph
            .derive(
                "same-derived-node",
                &[parent],
                op,
                None,
                &BTreeMap::new(),
                None,
            )
            .expect("single-parent derivation");
        let written = graph.node(node).expect("derived node");
        (written.declared_color_unverified().clone(), written.hash())
    };
    let (add_color, add_hash) = build(IntervalOp::Add);
    let (mul_color, mul_hash) = build(IntervalOp::Mul);
    assert_eq!(
        add_color, mul_color,
        "the hash test needs equal color payloads"
    );
    assert_ne!(
        add_hash, mul_hash,
        "the operation must change node identity"
    );

    let secret = 0x5EC2E7;
    let verifier = TestVerifier {
        keys: KeyMap::from([("operation-key".to_string(), secret)]),
    };
    let mut graph = ColorGraph::new();
    let parent = write_source(
        &mut graph,
        "estimated-parent",
        Color::Estimated {
            estimator: "rom".to_string(),
            dispersion: 0.1,
        },
    );
    let claimed = Color::Verified { lo: 0.0, hi: 1.0 };
    let grant = signed_grant(
        secret,
        "operation-key",
        "operation-bound-node",
        &claimed,
        vec![graph.node(parent).expect("parent").hash()],
        IntervalOp::Add,
        400,
    );
    let wrong_operation = graph.derive_waived(
        "operation-bound-node",
        &[parent],
        IntervalOp::Mul,
        claimed.clone(),
        &BTreeMap::new(),
        grant.clone(),
        &verifier,
        200,
    );
    assert!(matches!(
        wrong_operation,
        Err(ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::VerifierRefused { .. }
        })
    ));
    let authorized = graph
        .derive_waived(
            "operation-bound-node",
            &[parent],
            IntervalOp::Add,
            claimed,
            &BTreeMap::new(),
            grant,
            &verifier,
            200,
        )
        .expect("the grant authorizes only Add");
    assert_eq!(
        graph.node(authorized).expect("authorized").operation(),
        Some(IntervalOp::Add)
    );

    verdict(
        "col-009",
        true,
        "v4 node identity separates equal-payload Add and Mul nodes; v3 grant \
         signatures authorize exactly one operation",
    );
}

/// col-010 — choosing the waived API is itself an authorization claim. A
/// forged grant must refuse even when the claimed color needs no rank upgrade.
#[test]
fn col_010_no_upgrade_waiver_is_still_authenticated() {
    let verifier = TestVerifier {
        keys: KeyMap::from([("release-key".to_string(), 0xA11CE)]),
    };
    let mut graph = ColorGraph::new();
    let claimed = Color::Verified { lo: 0.0, hi: 1.0 };
    let parent = write_source(&mut graph, "verified-parent", claimed.clone());
    let forged = signed_grant(
        0xBAD,
        "release-key",
        "no-upgrade-node",
        &claimed,
        vec![graph.node(parent).expect("parent").hash()],
        IntervalOp::Hull,
        400,
    );
    let result = graph.derive_waived(
        "no-upgrade-node",
        &[parent],
        IntervalOp::Hull,
        claimed,
        &BTreeMap::new(),
        forged,
        &verifier,
        200,
    );
    assert!(matches!(
        result,
        Err(ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::VerifierRefused { .. }
        })
    ));
    assert_eq!(graph.nodes().len(), 1, "a refused grant writes no node");

    verdict(
        "col-010",
        true,
        "derive_waived authenticates every grant before writing, including exact \
         no-upgrade claims",
    );
}

/// col-011 — persisted rows carry the exact signed authorization material,
/// allowing an independent verifier to replay the decision. Tampering either
/// serialized payload or signature refuses.
#[test]
#[allow(clippy::too_many_lines)] // one persisted-grant round-trip and tamper matrix
fn col_011_serialized_grant_reverifies_and_tamper_refuses() {
    let secret = 0x51A1ED;
    let verifier = TestVerifier {
        keys: KeyMap::from([("persistence-key".to_string(), secret)]),
    };
    let mut graph = ColorGraph::new();
    let parent = write_source(
        &mut graph,
        "persisted-parent",
        Color::Estimated {
            estimator: "surrogate".to_string(),
            dispersion: 0.2,
        },
    );
    let parent_hash = graph.node(parent).expect("parent").hash();
    let claimed = Color::Verified { lo: -1.0, hi: 1.0 };
    let grant = signed_grant(
        secret,
        "persistence-key",
        "persisted-waiver",
        &claimed,
        vec![parent_hash],
        IntervalOp::Add,
        400,
    );
    let node = graph
        .derive_waived(
            "persisted-waiver",
            &[parent],
            IntervalOp::Add,
            claimed.clone(),
            &BTreeMap::new(),
            grant,
            &verifier,
            200,
        )
        .expect("valid grant writes");
    let row = graph
        .rows()
        .iter()
        .find(|row| row.contains(&format!("\"node\":{node}")))
        .expect("persisted color-write row");

    let serialized_key = json_string_field(row, "key_id").expect("serialized key id");
    let serialized_payload =
        decode_hex(json_string_field(row, "signing_payload_hex").expect("serialized payload"))
            .expect("payload hex");
    let serialized_signature =
        decode_hex(json_string_field(row, "signature_hex").expect("serialized signature"))
            .expect("signature hex");
    let stored_grant = graph
        .node(node)
        .expect("stored node")
        .grant()
        .expect("stored grant");
    assert_eq!(
        serialized_payload,
        stored_grant.signing_payload(IntervalOp::Add)
    );
    assert_ne!(
        serialized_payload,
        stored_grant.signing_payload(IntervalOp::Mul)
    );
    assert!(
        verifier
            .verify(serialized_key, &serialized_payload, &serialized_signature)
            .accepted()
    );
    assert!(row.contains("\"schema_version\":5"));
    assert!(row.contains("\"payload_version\":3"));
    assert!(row.contains("\"operation\":\"add\""));
    assert_eq!(
        json_string_field(row, "node_name"),
        Some("persisted-waiver")
    );
    let claimed_color_hex = test_hex(&claimed.canonical_bytes());
    assert_eq!(
        json_string_field(row, "claimed_color_hex"),
        Some(claimed_color_hex.as_str())
    );
    assert!(row.contains(&parent_hash.to_hex()));
    fs_ledger::Ledger::open(":memory:")
        .expect("validation ledger")
        .append_event(&fs_ledger::EventRow {
            session: None,
            t: 0,
            kind: "serialized-color-grant",
            payload: Some(row),
        })
        .expect("serialized grant row is strict JSON");

    let mut tampered_payload = serialized_payload.clone();
    tampered_payload[0] ^= 1;
    assert!(
        !verifier
            .verify(serialized_key, &tampered_payload, &serialized_signature)
            .accepted()
    );
    let mut tampered_signature = serialized_signature.clone();
    tampered_signature[0] ^= 1;
    assert!(
        !verifier
            .verify(serialized_key, &serialized_payload, &tampered_signature)
            .accepted()
    );

    verdict(
        "col-011",
        true,
        "schema-v2 rows persist the v3 payload, signature, node, canonical color, \
         operation, and parent hashes; replay verifies and tampering refuses",
    );
}

/// col-012 — positive leaves cannot assert their own color. Typed source
/// origins rederive the complete claim and bind every origin field into
/// provenance; mismatched or forged origins fail closed.
#[test]
#[allow(clippy::too_many_lines)] // one forged-source matrix across both positive color kinds
fn col_012_sources_are_minted_from_origin_evidence() {
    let mut graph = ColorGraph::new();
    let estimate = graph
        .source(
            "surrogate",
            Color::Estimated {
                estimator: "rom".to_string(),
                dispersion: 0.2,
            },
        )
        .expect("Estimated is the only direct source color");
    assert!(graph.node(estimate).is_some());
    for estimator in ["", "  ", "pending", "UNKNOWN", " rom", "rom "] {
        assert!(matches!(
            graph.source(
                "bad-estimator",
                Color::Estimated {
                    estimator: estimator.to_string(),
                    dispersion: 0.2,
                },
            ),
            Err(ColorWriteError::InvalidEstimatedSource {
                field: "estimator",
                ..
            })
        ));
    }
    for dispersion in [f64::NAN, -0.1, f64::NEG_INFINITY] {
        assert!(matches!(
            graph.source(
                "bad-dispersion",
                Color::Estimated {
                    estimator: "rom".to_string(),
                    dispersion,
                },
            ),
            Err(ColorWriteError::InvalidEstimatedSource {
                field: "dispersion",
                ..
            })
        ));
    }
    graph
        .source(
            "no-spread-claim",
            Color::Estimated {
                estimator: "rom".to_string(),
                dispersion: f64::INFINITY,
            },
        )
        .expect("positive infinity is the explicit no-spread-claim sentinel");
    assert!(graph.node(u64::MAX).is_none(), "invalid ids are checked");

    let verified = Color::Verified { lo: -2.0, hi: 3.0 };
    assert!(matches!(
        graph.source("forged-verified", verified.clone()),
        Err(ColorWriteError::SourceOriginRequired {
            rank: ColorRank::Verified
        })
    ));
    let unauthenticated_origin = SourceOrigin::Certificate {
        producer: "fs-ivl/enclosure".to_string(),
        certificate_hash: hash_bytes(b"unauthenticated certificate"),
        certificate: NumericalCertificate::enclosure(-2.0, 3.0),
    };
    assert!(matches!(
        graph.source_with_origin(
            "unauthenticated",
            &verified,
            unauthenticated_origin,
            &NoSourceOriginVerifier,
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::VerifierRefused { .. }
        })
    ));
    let verified_id = admitted_source(
        &mut graph,
        "certified",
        &verified,
        SourceOrigin::Certificate {
            producer: "fs-ivl/enclosure".to_string(),
            certificate_hash: hash_bytes(b"certified enclosure"),
            certificate: NumericalCertificate::enclosure(-2.0, 3.0),
        },
    )
    .expect("matching enclosure mints Verified");
    assert_eq!(
        graph
            .node(verified_id)
            .expect("verified")
            .declared_color_unverified(),
        &verified
    );
    let authorized_origin = SourceOrigin::Certificate {
        producer: "fs-ivl/enclosure".to_string(),
        certificate_hash: hash_bytes(b"origin-bound enclosure"),
        certificate: NumericalCertificate::enclosure(-2.0, 3.0),
    };
    let exact_verifier =
        TestSourceVerifier::authorizing("origin-bound", &verified, &authorized_origin);
    assert!(matches!(
        graph.source_with_origin(
            "origin-bound",
            &verified,
            SourceOrigin::Certificate {
                producer: "forged-producer".to_string(),
                certificate_hash: hash_bytes(b"origin-bound enclosure"),
                certificate: NumericalCertificate::enclosure(-2.0, 3.0),
            },
            &exact_verifier,
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::VerifierRefused { .. }
        })
    ));
    assert!(matches!(
        graph.source_with_origin(
            "different-node",
            &verified,
            authorized_origin,
            &exact_verifier,
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::VerifierRefused { .. }
        })
    ));
    assert!(matches!(
        admitted_source(
            &mut graph,
            "mismatched-cert",
            &verified,
            SourceOrigin::Certificate {
                producer: "fs-ivl/enclosure".to_string(),
                certificate_hash: hash_bytes(b"mismatched enclosure"),
                certificate: NumericalCertificate::enclosure(-2.0, 4.0),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::CertificateMismatch
        })
    ));
    assert!(matches!(
        admitted_source(
            &mut graph,
            "padded-producer",
            &Color::Verified { lo: 0.0, hi: 1.0 },
            SourceOrigin::Certificate {
                producer: " fs-ivl/enclosure".to_string(),
                certificate_hash: hash_bytes(b"padded producer enclosure"),
                certificate: NumericalCertificate::enclosure(0.0, 1.0),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::BlankProducer
        })
    ));
    assert!(matches!(
        admitted_source(
            &mut graph,
            "estimated-cert",
            &verified,
            SourceOrigin::Certificate {
                producer: "surrogate".to_string(),
                certificate_hash: hash_bytes(b"estimated certificate"),
                certificate: NumericalCertificate::estimate(-2.0, 3.0),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::CertificateRefused { .. }
        })
    ));
    assert!(matches!(
        admitted_source(
            &mut graph,
            "placeholder-producer",
            &Color::Verified { lo: 0.0, hi: 1.0 },
            SourceOrigin::Certificate {
                producer: "pending".to_string(),
                certificate_hash: hash_bytes(b"placeholder producer enclosure"),
                certificate: NumericalCertificate::enclosure(0.0, 1.0),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::BlankProducer
        })
    ));

    let regime = ValidityDomain::unconstrained().with("re", 1e3, 1e5);
    let validated = Color::Validated {
        regime: regime.clone(),
        dataset: "wind-tunnel-a".to_string(),
    };
    assert!(matches!(
        graph.source("forged-validated", validated.clone()),
        Err(ColorWriteError::SourceOriginRequired {
            rank: ColorRank::Validated
        })
    ));
    let anchor = SourceOrigin::Anchoring {
        dataset_id: "wind-tunnel-a".to_string(),
        content_hash: hash_bytes(b"campaign-a"),
        regime: regime.clone(),
    };
    let validated_id = admitted_source(&mut graph, "anchored", &validated, anchor.clone())
        .expect("matching anchor mints Validated");
    assert_eq!(
        graph
            .node(validated_id)
            .expect("validated")
            .declared_color_unverified(),
        &validated
    );
    let anchor_verifier = TestSourceVerifier::authorizing("anchor-bound", &validated, &anchor);
    assert!(matches!(
        graph.source_with_origin(
            "anchor-bound",
            &validated,
            SourceOrigin::Anchoring {
                dataset_id: "wind-tunnel-a".to_string(),
                content_hash: hash_bytes(b"campaign-a-tampered"),
                regime: regime.clone(),
            },
            &anchor_verifier,
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::VerifierRefused { .. }
        })
    ));
    let wrong_regime = SourceOrigin::Anchoring {
        dataset_id: "wind-tunnel-a".to_string(),
        content_hash: hash_bytes(b"campaign-a"),
        regime: ValidityDomain::unconstrained().with("re", 1e3, 2e5),
    };
    assert!(matches!(
        admitted_source(&mut graph, "wrong-regime", &validated, wrong_regime),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::RegimeMismatch
        })
    ));
    let empty_regime_color = Color::Validated {
        regime: ValidityDomain::unconstrained(),
        dataset: "wind-tunnel-a".to_string(),
    };
    assert!(matches!(
        admitted_source(
            &mut graph,
            "empty-regime",
            &empty_regime_color,
            SourceOrigin::Anchoring {
                dataset_id: "wind-tunnel-a".to_string(),
                content_hash: hash_bytes(b"campaign-a"),
                regime: ValidityDomain::unconstrained(),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::InvalidRegime { .. }
        })
    ));
    for bad_regime in [
        ValidityDomain::unconstrained().with("pending", 0.0, 1.0),
        ValidityDomain::unconstrained().with(" re ", 0.0, 1.0),
        ValidityDomain::unconstrained().with("re", f64::NAN, 1.0),
        ValidityDomain::unconstrained().with("re", 0.0, f64::INFINITY),
    ] {
        let bad_color = Color::Validated {
            regime: bad_regime.clone(),
            dataset: "wind-tunnel-a".to_string(),
        };
        assert!(matches!(
            admitted_source(
                &mut graph,
                "bad-regime",
                &bad_color,
                SourceOrigin::Anchoring {
                    dataset_id: "wind-tunnel-a".to_string(),
                    content_hash: hash_bytes(b"campaign-a"),
                    regime: bad_regime,
                },
            ),
            Err(ColorWriteError::SourceOriginRefused {
                rejection: SourceOriginRejection::InvalidRegime { .. }
            })
        ));
    }
    let placeholder_dataset_color = Color::Validated {
        regime: regime.clone(),
        dataset: "pending".to_string(),
    };
    assert!(matches!(
        admitted_source(
            &mut graph,
            "placeholder-dataset",
            &placeholder_dataset_color,
            SourceOrigin::Anchoring {
                dataset_id: "pending".to_string(),
                content_hash: hash_bytes(b"campaign-a"),
                regime: regime.clone(),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::BlankDataset
        })
    ));
    let padded_dataset_color = Color::Validated {
        regime: regime.clone(),
        dataset: "wind-tunnel-a ".to_string(),
    };
    assert!(matches!(
        admitted_source(
            &mut graph,
            "padded-dataset",
            &padded_dataset_color,
            SourceOrigin::Anchoring {
                dataset_id: "wind-tunnel-a ".to_string(),
                content_hash: hash_bytes(b"campaign-a"),
                regime: regime.clone(),
            },
        ),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::BlankDataset
        })
    ));

    let mut changed_anchor = ColorGraph::new();
    let changed = admitted_source(
        &mut changed_anchor,
        "anchored",
        &validated,
        SourceOrigin::Anchoring {
            dataset_id: "wind-tunnel-a".to_string(),
            content_hash: hash_bytes(b"campaign-a-tampered"),
            regime,
        },
    )
    .expect("a different real artifact is admissible but has a different identity");
    assert_ne!(
        graph.node(validated_id).expect("original").hash(),
        changed_anchor.node(changed).expect("changed").hash(),
        "the anchoring artifact hash is provenance, not an unbound annotation"
    );
    graph.verify_replay().expect("source origins replay");
}

/// col-013 — source waivers use their own scope and v4 signing payload.
/// A derived grant cannot be replayed at a leaf; payload or signature
/// tampering refuses; the persisted row independently re-verifies.
#[test]
fn col_013_source_waiver_is_separate_and_replayable() {
    let secret = 0x50A7CE;
    let verifier = TestVerifier {
        keys: KeyMap::from([("source-key".to_string(), secret)]),
    };
    let color = Color::Verified { lo: 4.0, hi: 5.0 };
    let grant = signed_source_grant(secret, "source-key", "exceptional-source", &color, 400);
    assert_ne!(
        grant.signing_payload_source(),
        grant.signing_payload(IntervalOp::Hull),
        "source and derive authority are domain-separated"
    );
    let mut graph = ColorGraph::new();
    let id = graph
        .source_waived(
            "exceptional-source",
            color.clone(),
            grant.clone(),
            &verifier,
            200,
        )
        .expect("valid source grant");
    let row = graph
        .rows()
        .iter()
        .find(|row| row.contains(&format!("\"node\":{id}")))
        .expect("source row");
    let payload =
        decode_hex(json_string_field(row, "signing_payload_hex").expect("source payload field"))
            .expect("source payload hex");
    let signature =
        decode_hex(json_string_field(row, "signature_hex").expect("source signature field"))
            .expect("source signature hex");
    assert_eq!(payload, grant.signing_payload_source());
    assert!(row.contains("\"payload_version\":4"));
    assert!(row.contains("\"schema_version\":5"));
    assert!(
        verifier
            .verify("source-key", &payload, &signature)
            .accepted()
    );
    graph.verify_replay().expect("source grant fields replay");

    let mut tampered = grant.clone();
    tampered.annotation.reason.push_str(" edited");
    assert!(matches!(
        ColorGraph::new().source_waived(
            "exceptional-source",
            color.clone(),
            tampered,
            &verifier,
            200,
        ),
        Err(ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::VerifierRefused { .. }
        })
    ));
    let derive_grant = signed_grant(
        secret,
        "source-key",
        "exceptional-source",
        &color,
        Vec::new(),
        IntervalOp::Hull,
        400,
    );
    assert!(matches!(
        ColorGraph::new().source_waived("exceptional-source", color, derive_grant, &verifier, 200,),
        Err(ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::ScopeMismatch
        })
    ));
}

/// col-014 — every parent's regime exit survives in canonical parent-list
/// order and is part of the replayable node identity.
#[test]
fn col_014_all_regime_demotions_are_preserved() {
    let mut graph = ColorGraph::new();
    let first = write_source(
        &mut graph,
        "first-anchor",
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("mach", 0.0, 0.8),
            dataset: "tunnel-mach".to_string(),
        },
    );
    let second = write_source(
        &mut graph,
        "second-anchor",
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("re", 1e3, 1e5),
            dataset: "tunnel-re".to_string(),
        },
    );
    let state: BTreeMap<String, f64> = [("mach".to_string(), 1.2), ("re".to_string(), 5e5)].into();
    let derived = graph
        .derive(
            "both-exited",
            &[first, second],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("both parents demote");
    let demotions = graph.node(derived).expect("derived").demotions();
    assert_eq!(demotions.len(), 2);
    assert_eq!(
        (demotions[0].parent_index(), demotions[0].parent_id()),
        (0, first)
    );
    assert_eq!(
        (demotions[1].parent_index(), demotions[1].parent_id()),
        (1, second)
    );
    assert_eq!(demotions[0].reason().axis, "mach");
    assert_eq!(demotions[1].reason().axis, "re");
    let rows: Vec<&str> = graph
        .rows()
        .iter()
        .filter(|row| row.contains("\"event\":\"demotion\""))
        .map(String::as_str)
        .collect();
    assert_eq!(rows.len(), 2);
    assert!(rows[0].contains("\"parent_index\":0"));
    assert!(rows[1].contains("\"parent_index\":1"));
    let duplicate = graph
        .derive(
            "duplicate-parent-exit",
            &[first, first],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("both occurrences of one parent demote");
    let duplicate_demotions = graph.node(duplicate).expect("duplicate").demotions();
    assert_eq!(duplicate_demotions.len(), 2);
    assert_eq!(duplicate_demotions[0].parent_index(), 0);
    assert_eq!(duplicate_demotions[1].parent_index(), 1);
    assert_eq!(duplicate_demotions[0].parent_id(), first);
    assert_eq!(duplicate_demotions[1].parent_id(), first);

    let mut changed = ColorGraph::new();
    let changed_first = write_source(
        &mut changed,
        "first-anchor",
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("mach", 0.0, 0.8),
            dataset: "tunnel-mach".to_string(),
        },
    );
    let changed_second = write_source(
        &mut changed,
        "second-anchor",
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("re", 1e3, 1e5),
            dataset: "tunnel-re".to_string(),
        },
    );
    let changed_state: BTreeMap<String, f64> =
        [("mach".to_string(), 1.3), ("re".to_string(), 6e5)].into();
    let changed_derived = changed
        .derive(
            "both-exited",
            &[changed_first, changed_second],
            IntervalOp::Hull,
            None,
            &changed_state,
            None,
        )
        .expect("changed exits");
    assert_ne!(
        graph.node(derived).expect("original").hash(),
        changed.node(changed_derived).expect("changed").hash(),
        "demotion values participate in the node hash"
    );
    graph.verify_replay().expect("all demotions replay");
}

/// col-015 — authentication authorizes policy, never malformed epistemic
/// payloads. Both source and derived waiver doors reject every malformed
/// positive color; the derived door also rejects malformed estimates.
#[test]
#[allow(clippy::too_many_lines)] // all color variants form one structural-invariant matrix
fn col_015_waivers_cannot_authorize_malformed_colors() {
    let secret = 0xC010_5A7Eu64;
    let key_id = "structure-key";
    let verifier = TestVerifier {
        keys: KeyMap::from([(key_id.to_string(), secret)]),
    };
    let valid_regime = ValidityDomain::unconstrained().with("re", 1e3, 1e5);
    let malformed_positive = vec![
        (
            "verified-nan",
            Color::Verified {
                lo: f64::NAN,
                hi: 1.0,
            },
        ),
        ("verified-inverted", Color::Verified { lo: 2.0, hi: 1.0 }),
        (
            "validated-blank-dataset",
            Color::Validated {
                regime: valid_regime.clone(),
                dataset: String::new(),
            },
        ),
        (
            "validated-placeholder-dataset",
            Color::Validated {
                regime: valid_regime.clone(),
                dataset: "pending".to_string(),
            },
        ),
        (
            "validated-padded-dataset",
            Color::Validated {
                regime: valid_regime.clone(),
                dataset: "tunnel-a ".to_string(),
            },
        ),
        (
            "validated-empty-regime",
            Color::Validated {
                regime: ValidityDomain::unconstrained(),
                dataset: "tunnel-a".to_string(),
            },
        ),
        (
            "validated-padded-axis",
            Color::Validated {
                regime: ValidityDomain::unconstrained().with(" re", 1e3, 1e5),
                dataset: "tunnel-a".to_string(),
            },
        ),
        (
            "validated-nonfinite-regime",
            Color::Validated {
                regime: ValidityDomain::unconstrained().with("re", f64::NAN, 1e5),
                dataset: "tunnel-a".to_string(),
            },
        ),
        (
            "validated-infinite-regime",
            Color::Validated {
                regime: ValidityDomain::unconstrained().with("re", 1e3, f64::INFINITY),
                dataset: "tunnel-a".to_string(),
            },
        ),
    ];

    for (label, color) in &malformed_positive {
        let name = format!("source-{label}");
        let grant = signed_source_grant(secret, key_id, &name, color, 400);
        let rejection = structural_refusal(ColorGraph::new().source_waived(
            &name,
            color.clone(),
            grant,
            &verifier,
            200,
        ));
        assert!(
            matches!(
                rejection,
                ColorStructureRejection::InvalidIdentity { .. }
                    | ColorStructureRejection::InvalidVerifiedInterval { .. }
                    | ColorStructureRejection::InvalidValidatedRegime { .. }
            ),
            "unexpected source refusal for {label}: {rejection:?}"
        );
    }

    let mut graph = ColorGraph::new();
    let parent = graph
        .source(
            "valid-parent",
            Color::Estimated {
                estimator: "rom".to_string(),
                dispersion: 0.2,
            },
        )
        .expect("valid parent");
    let lineage = vec![graph.node(parent).expect("parent").hash()];
    let malformed_estimates = vec![
        (
            "estimated-blank",
            Color::Estimated {
                estimator: String::new(),
                dispersion: 0.2,
            },
        ),
        (
            "estimated-placeholder",
            Color::Estimated {
                estimator: "pending".to_string(),
                dispersion: 0.2,
            },
        ),
        (
            "estimated-padded",
            Color::Estimated {
                estimator: " rom".to_string(),
                dispersion: 0.2,
            },
        ),
        (
            "estimated-nan",
            Color::Estimated {
                estimator: "rom".to_string(),
                dispersion: f64::NAN,
            },
        ),
        (
            "estimated-negative",
            Color::Estimated {
                estimator: "rom".to_string(),
                dispersion: -0.1,
            },
        ),
    ];
    for (label, color) in malformed_positive.into_iter().chain(malformed_estimates) {
        let name = format!("derived-{label}");
        let grant = signed_grant(
            secret,
            key_id,
            &name,
            &color,
            lineage.clone(),
            IntervalOp::Hull,
            400,
        );
        structural_refusal(graph.derive_waived(
            &name,
            &[parent],
            IntervalOp::Hull,
            color,
            &BTreeMap::new(),
            grant,
            &verifier,
            200,
        ));
    }
    assert!(
        graph.node(parent + 1).is_none(),
        "no structurally invalid waived node was appended"
    );
    graph.verify_replay().expect("valid prefix still replays");
    verdict(
        "col-015",
        true,
        "valid source and derive signatures cannot authorize NaN/inverted intervals, malformed \
         validated identities/regimes, or malformed estimates",
    );
}

/// col-015b — ordinary composition cannot append a color that its own replay
/// validator rejects, and public 64-bit parent ids never truncate into a
/// platform-sized index.
#[test]
fn col_015b_overflow_and_wide_parent_ids_fail_before_append() {
    let mut graph = ColorGraph::new();
    let left = write_source(
        &mut graph,
        "overflow-left",
        Color::Verified {
            lo: f64::MAX,
            hi: f64::MAX,
        },
    );
    let right = write_source(
        &mut graph,
        "overflow-right",
        Color::Verified {
            lo: f64::MAX,
            hi: f64::MAX,
        },
    );
    let node_count = graph.nodes().len();
    let row_count = graph.rows().len();
    assert!(matches!(
        graph.derive(
            "overflowed-sum",
            &[left, right],
            IntervalOp::Add,
            None,
            &BTreeMap::new(),
            None,
        ),
        Err(ColorWriteError::InvalidColor {
            rejection: ColorStructureRejection::InvalidVerifiedInterval { .. }
        })
    ));
    assert_eq!(graph.nodes().len(), node_count);
    assert_eq!(graph.rows().len(), row_count);
    graph.verify_replay().expect("valid prefix still replays");

    assert!(matches!(
        graph.derive(
            "wide-parent",
            &[u64::MAX],
            IntervalOp::Hull,
            None,
            &BTreeMap::new(),
            None,
        ),
        Err(ColorWriteError::UnknownParent { id: u64::MAX })
    ));
    let claimed = Color::Verified { lo: 0.0, hi: 1.0 };
    let grant = signed_grant(
        7,
        "wide-parent-key",
        "wide-parent-waived",
        &claimed,
        Vec::new(),
        IntervalOp::Hull,
        400,
    );
    assert!(matches!(
        graph.derive_waived(
            "wide-parent-waived",
            &[u64::MAX],
            IntervalOp::Hull,
            claimed,
            &BTreeMap::new(),
            grant,
            &NoWaiverVerifier,
            200,
        ),
        Err(ColorWriteError::UnknownParent { id: u64::MAX })
    ));
}

/// col-015c — source provenance distinguishes both the retained certificate
/// artifact and the policy that authenticated it. Hostile verifier callbacks
/// fail closed without appending a partial node.
#[test]
fn col_015c_source_artifact_and_policy_are_hash_bound() {
    let color = Color::Verified { lo: 1.0, hi: 2.0 };
    let origin_a = SourceOrigin::Certificate {
        producer: "fs-ivl/enclosure".to_string(),
        certificate_hash: hash_bytes(b"certificate artifact A"),
        certificate: NumericalCertificate::enclosure(1.0, 2.0),
    };
    let origin_b = SourceOrigin::Certificate {
        producer: "fs-ivl/enclosure".to_string(),
        certificate_hash: hash_bytes(b"certificate artifact B"),
        certificate: NumericalCertificate::enclosure(1.0, 2.0),
    };
    let policy_a = hash_bytes(b"source trust policy A");
    let policy_b = hash_bytes(b"source trust policy B");
    let request_a = SourceOriginRequest::new("certified", &color, &origin_a).canonical_bytes();
    let request_b = SourceOriginRequest::new("certified", &color, &origin_b).canonical_bytes();

    let write = |origin: SourceOrigin, request: Vec<u8>, policy_fingerprint| {
        let verifier = TestSourceVerifier {
            accepted_request: request,
            policy_fingerprint,
        };
        let mut graph = ColorGraph::new();
        let id = graph
            .source_with_origin("certified", &color, origin, &verifier)
            .expect("authorized source");
        graph.verify_replay().expect("source replays");
        (graph, id)
    };

    let (graph_a, id_a) = write(origin_a.clone(), request_a.clone(), policy_a);
    let (artifact_changed, artifact_id) = write(origin_b, request_b, policy_a);
    let (policy_changed, policy_id) = write(origin_a.clone(), request_a, policy_b);
    let node_a = graph_a.node(id_a).expect("source A");
    assert_eq!(node_a.origin_policy_fingerprint(), Some(policy_a));
    assert_ne!(
        node_a.hash(),
        artifact_changed
            .node(artifact_id)
            .expect("artifact B")
            .hash()
    );
    assert_ne!(
        node_a.hash(),
        policy_changed.node(policy_id).expect("policy B").hash()
    );
    let row = graph_a.rows().last().expect("source row");
    assert!(row.contains(&hash_bytes(b"certificate artifact A").to_hex()));
    assert!(row.contains(&policy_a.to_hex()));

    let mut panicked = ColorGraph::new();
    assert!(matches!(
        panicked.source_with_origin("certified", &color, origin_a, &PanickingSourceVerifier),
        Err(ColorWriteError::SourceOriginRefused {
            rejection: SourceOriginRejection::VerifierPanicked
        })
    ));
    assert!(panicked.nodes().is_empty() && panicked.rows().is_empty());

    let grant = signed_source_grant(11, "panic-key", "waived", &color, 400);
    assert!(matches!(
        panicked.source_waived("waived", color, grant, &PanickingWaiverVerifier, 200,),
        Err(ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::VerifierPanicked
        })
    ));
    assert!(panicked.nodes().is_empty() && panicked.rows().is_empty());
}

/// col-016 — a waiver stays visible through ordinary and waived descendants.
/// Duplicate parents do not duplicate dependencies; independent branches form
/// a canonical union; a node's own grant enters its children's closure.
#[test]
#[allow(clippy::too_many_lines)] // one transitive authority-closure graph and its assertions
fn col_016_waiver_dependencies_propagate_transitively() {
    let secret = 0x7A17_2026u64;
    let key_id = "lineage-key";
    let verifier = TestVerifier {
        keys: KeyMap::from([(key_id.to_string(), secret)]),
    };
    let state = BTreeMap::new();
    let mut graph = ColorGraph::new();

    let color_a = Color::Verified { lo: 1.0, hi: 2.0 };
    let mut grant_a = signed_source_grant(secret, key_id, "waived-a", &color_a, 400);
    grant_a.annotation.id = "waiver-a".to_string();
    grant_a.signature = test_mac(secret, &grant_a.signing_payload_source());
    let source_a = graph
        .source_waived("waived-a", color_a, grant_a, &verifier, 200)
        .expect("first waived source");
    assert!(graph.node(source_a).expect("source A").depends_on_waiver());
    assert!(
        graph
            .node(source_a)
            .expect("source A")
            .waiver_dependencies()
            .is_empty(),
        "a direct grant is not duplicated as an inherited dependency"
    );
    let first_generation = graph
        .derive(
            "first-generation",
            &[source_a],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("ordinary child retains source waiver");
    assert!(
        graph
            .node(first_generation)
            .expect("first generation")
            .depends_on_waiver()
    );
    let duplicate_parent = graph
        .derive(
            "duplicate-parent",
            &[first_generation, first_generation],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("duplicate parent dependency is deduplicated");
    assert_eq!(
        graph
            .node(duplicate_parent)
            .expect("duplicate child")
            .waiver_dependencies()
            .len(),
        1
    );

    let color_b = Color::Verified { lo: 3.0, hi: 4.0 };
    let mut grant_b = signed_source_grant(secret, key_id, "waived-b", &color_b, 400);
    grant_b.annotation.id = "waiver-b".to_string();
    grant_b.signature = test_mac(secret, &grant_b.signing_payload_source());
    let source_b = graph
        .source_waived("waived-b", color_b, grant_b, &verifier, 200)
        .expect("second waived source");
    let merged = graph
        .derive(
            "merged-waivers",
            &[duplicate_parent, source_b],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("ordinary merge retains both source waivers");
    assert_eq!(
        graph
            .node(merged)
            .expect("merged")
            .waiver_dependencies()
            .iter()
            .map(WaiverDependency::authorizing_node)
            .collect::<Vec<_>>(),
        vec![source_a, source_b]
    );

    let merged_color = graph
        .node(merged)
        .expect("merged")
        .declared_color_unverified()
        .clone();
    let mut grant_c = signed_grant(
        secret,
        key_id,
        "waived-derived",
        &merged_color,
        vec![graph.node(merged).expect("merged").hash()],
        IntervalOp::Hull,
        400,
    );
    grant_c.annotation.id = "waiver-c".to_string();
    grant_c.signature = test_mac(secret, &grant_c.signing_payload(IntervalOp::Hull));
    let waived_derived = graph
        .derive_waived(
            "waived-derived",
            &[merged],
            IntervalOp::Hull,
            merged_color,
            &state,
            grant_c,
            &verifier,
            200,
        )
        .expect("waived child retains inherited dependencies and its own grant");
    let final_node = graph
        .derive(
            "ordinary-grandchild",
            &[waived_derived],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("ordinary grandchild inherits all waivers");
    let final_node = graph.node(final_node).expect("final node");
    assert!(
        final_node.grant().is_none() && final_node.waiver().is_none(),
        "the node did not need a new waiver"
    );
    assert!(final_node.depends_on_waiver());
    assert_eq!(
        final_node
            .waiver_dependencies()
            .iter()
            .map(WaiverDependency::authorizing_node)
            .collect::<Vec<_>>(),
        vec![source_a, source_b, waived_derived]
    );
    assert_eq!(
        final_node
            .waiver_dependencies()
            .iter()
            .map(WaiverDependency::operation)
            .collect::<Vec<_>>(),
        vec![None, None, Some(IntervalOp::Hull)]
    );
    assert_eq!(
        final_node
            .waiver_dependencies()
            .iter()
            .map(|dependency| dependency.grant().annotation.id.as_str())
            .collect::<Vec<_>>(),
        vec!["waiver-a", "waiver-b", "waiver-c"]
    );
    let policy = graph
        .node(source_a)
        .expect("source A")
        .waiver_policy_fingerprint()
        .expect("direct source waiver policy");
    assert_eq!(
        final_node
            .waiver_dependencies()
            .iter()
            .map(WaiverDependency::policy_fingerprint)
            .collect::<Vec<_>>(),
        vec![policy, policy, policy]
    );
    assert_eq!(
        final_node
            .waiver_dependencies()
            .iter()
            .map(WaiverDependency::admission_day)
            .collect::<Vec<_>>(),
        vec![200, 200, 200]
    );
    let row = graph
        .rows()
        .iter()
        .find(|row| row.contains("\"name\":\"ordinary-grandchild\""))
        .expect("final color row");
    assert!(row.contains("\"schema_version\":5"));
    assert!(row.contains("\"waiver_dependencies\":["));
    assert!(row.contains("\"admission_day\":200"));
    assert!(row.contains(&policy.to_hex()));
    assert!(row.contains("waiver-a") && row.contains("waiver-b") && row.contains("waiver-c"));
    let clean = graph
        .source(
            "clean-estimate",
            Color::Estimated {
                estimator: "independent-rom".to_string(),
                dispersion: 0.1,
            },
        )
        .expect("clean source");
    assert!(
        !graph.node(clean).expect("clean").depends_on_waiver(),
        "clean nodes must remain distinguishable from waived dependencies"
    );
    graph
        .verify_replay()
        .expect("multi-generation waiver union replays");
}

/// col-017 — authenticated metadata is structural data, not something an
/// injected verifier may bless. Waiver taint also blocks the ordinary
/// scientific-color accessor, and graph fan-in is bounded before append.
#[test]
#[allow(clippy::too_many_lines)] // one adversarial authority-shape matrix and taint proof
fn col_017_waiver_shape_taint_and_resource_limits_fail_closed() {
    let color = Color::Verified { lo: 1.0, hi: 2.0 };
    let valid = signed_source_grant(17, "release-key", "waived-source", &color, 400);

    let refusal = |grant: WaiverGrant| {
        let mut graph = ColorGraph::new();
        let result = graph.source_waived(
            "waived-source",
            color.clone(),
            grant,
            &PanickingWaiverVerifier,
            200,
        );
        assert!(graph.nodes().is_empty(), "a refused grant must not append");
        match result {
            Err(ColorWriteError::WaiverRefused { rejection }) => rejection,
            other => panic!("expected structural waiver refusal, got {other:?}"),
        }
    };

    let mut malformed = valid.clone();
    malformed.key_id.clear();
    assert!(matches!(
        refusal(malformed),
        WaiverRejection::InvalidField {
            field: "key_id",
            reason: "blank"
        }
    ));

    let mut malformed = valid.clone();
    malformed.annotation.id = "TODO".to_string();
    assert!(matches!(
        refusal(malformed),
        WaiverRejection::InvalidField {
            field: "waiver_id",
            reason: "placeholder"
        }
    ));

    let mut malformed = valid.clone();
    malformed.annotation.signer = "chief\u{202e}reenigne".to_string();
    assert!(matches!(
        refusal(malformed),
        WaiverRejection::InvalidField {
            field: "signer",
            reason: "invalid-character"
        }
    ));

    let mut malformed = valid.clone();
    malformed.annotation.reason = "approved\nwithout a stable record".to_string();
    assert!(matches!(
        refusal(malformed),
        WaiverRejection::InvalidField {
            field: "reason",
            reason: "control-character"
        }
    ));

    let mut malformed = valid.clone();
    malformed.signature.clear();
    assert!(matches!(
        refusal(malformed),
        WaiverRejection::InvalidField {
            field: "signature",
            reason: "blank"
        }
    ));

    let verifier = TestVerifier {
        keys: KeyMap::from([("release-key".to_string(), 17)]),
    };
    let mut graph = ColorGraph::new();
    let waived = graph
        .source_waived("waived-source", color, valid, &verifier, 200)
        .expect("well-formed signed grant");
    let waived = graph.node(waived).expect("waived node");
    assert!(matches!(
        waived.declared_color_unverified(),
        Color::Verified { .. }
    ));
    assert!(
        waived.scientific_color().is_none(),
        "waiver taint must not detach from the scientific accessor"
    );

    let clean = graph
        .source(
            "clean-estimate",
            Color::Estimated {
                estimator: "rom".to_string(),
                dispersion: 0.1,
            },
        )
        .expect("clean source");
    assert!(
        graph
            .node(clean)
            .expect("clean node")
            .scientific_color()
            .is_some()
    );

    let oversized_parents = vec![clean; MAX_COLOR_PARENTS + 1];
    assert!(matches!(
        graph.derive(
            "oversized-fan-in",
            &oversized_parents,
            IntervalOp::Hull,
            None,
            &BTreeMap::new(),
            None,
        ),
        Err(ColorWriteError::ResourceLimitExceeded {
            resource: "parents",
            limit: MAX_COLOR_PARENTS,
            actual
        }) if actual == MAX_COLOR_PARENTS + 1
    ));
    graph.verify_replay().expect("accepted prefix replays");
}

/// col-007 — color rows remain strict JSON when fail-closed demotion emits
/// non-finite sentinels and caller-controlled metadata contains JSON syntax or
/// control characters. Validation goes through the ledger's SQLite `json_valid`
/// path, the same parser that enforces persisted payloads.
#[test]
fn col_007_color_rows_are_strict_json_under_hostile_metadata() {
    let build = || -> Vec<String> {
        let hostile = "meta\"\\\n\r\t\u{0007}";
        let axis = "re-json".to_string();
        let dataset = "anchors-json".to_string();
        let mut graph = ColorGraph::new();
        let validated = write_source(
            &mut graph,
            &format!("validated-{hostile}"),
            Color::Validated {
                regime: ValidityDomain::unconstrained().with(&axis, 1.0, 10.0),
                dataset,
            },
        );
        let state: BTreeMap<String, f64> = [(axis, f64::NAN)].into();
        graph
            .derive(
                &format!("demoted-{hostile}"),
                &[validated],
                IntervalOp::Hull,
                None,
                &state,
                None,
            )
            .expect("non-finite state demotes instead of refusing the write");

        let estimated = write_source(
            &mut graph,
            &format!("estimated-{hostile}"),
            Color::Estimated {
                estimator: "surrogate-json".to_string(),
                dispersion: f64::INFINITY,
            },
        );
        // Hostile ANNOTATION strings must render as strict JSON; the
        // claim stays within the derivation cap (annotations no longer
        // authorize upgrades — qmao.1.1).
        graph
            .derive(
                &format!("waived-{hostile}"),
                &[estimated],
                IntervalOp::Hull,
                None,
                &BTreeMap::new(),
                Some(Waiver {
                    id: format!("id-{hostile}"),
                    signer: format!("signer-{hostile}"),
                    reason: format!("reason-{hostile}"),
                }),
            )
            .expect("annotated write");
        graph.rows().to_vec()
    };

    let rows = build();
    let deterministic = rows == build();
    let ledger = fs_ledger::Ledger::open(":memory:").expect("open validation ledger");
    let mut parser_accepts_every_row = true;
    for (index, row) in rows.iter().enumerate() {
        parser_accepts_every_row &= ledger
            .append_event(&fs_ledger::EventRow {
                session: None,
                t: i64::try_from(index).expect("small row index"),
                kind: "color-json-validation",
                payload: Some(row),
            })
            .is_ok();
    }
    let no_raw_controls = rows.iter().all(|row| !row.chars().any(char::is_control));
    let sentinels_and_escapes_present = rows.iter().any(|row| row.contains("non-finite:NaN"))
        && rows.iter().any(|row| row.contains("non-finite:inf"))
        && rows.iter().any(|row| row.contains(r#"\""#))
        && rows.iter().any(|row| row.contains(r"\\"))
        && rows.iter().any(|row| row.contains(r"\n"))
        && rows.iter().any(|row| row.contains(r"\u0007"));

    verdict(
        "col-007",
        deterministic
            && parser_accepts_every_row
            && no_raw_controls
            && sentinels_and_escapes_present,
        "SQLite json_valid accepts every deterministic color/demotion/waiver row; \
         NaN and infinity are tagged strings and hostile metadata is escaped",
    );
}

/// v3 migration regression (bead lmp4.3): the speculation extension
/// table exists, round-trips the four solve-node fields, and every
/// pre-existing table still answers queries (nothing broke).
#[test]
fn speculation_schema_migration() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("open");
    assert_eq!(ledger.schema_version().expect("version"), 3);
    let body = "{\"proposer_id\":\"neighbor-extrapolation\",\"accepted\":true,\
                \"bound\":3.2e-4,\"iterations_saved\":4}";
    ledger
        .put_extension(fs_ledger::ExtensionTable::Speculation, "solve-op-17", body)
        .expect("put");
    let back = ledger
        .get_extension(fs_ledger::ExtensionTable::Speculation, "solve-op-17")
        .expect("get")
        .expect("present");
    assert!(back.contains("iterations_saved"), "{back}");
    // Existing tables unbroken.
    for table in fs_ledger::ALL_TABLES {
        let _ = ledger.table_count(table).expect("old queries still work");
    }
}
