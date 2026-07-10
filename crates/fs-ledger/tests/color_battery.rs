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
    ColorGraph, ColorWriteError, NoWaiverVerifier, WAIVER_SCOPE_COLOR_UPGRADE, Waiver, WaiverGrant,
    WaiverRejection, WaiverVerifier,
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
    fn verify(&self, key_id: &str, payload: &[u8], signature: &[u8]) -> bool {
        self.keys
            .get(key_id)
            .is_some_and(|&secret| test_mac(secret, payload) == signature)
    }
}

fn signed_grant(
    secret: u64,
    key_id: &str,
    name: &str,
    color: &Color,
    parent_hashes: Vec<fs_ledger::ContentHash>,
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
        claimed_color: color.name().to_string(),
        parent_hashes,
        expires_day,
        signature: Vec::new(),
    };
    grant.signature = test_mac(secret, &grant.signing_payload());
    grant
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
    let clean = g.source("fem-bound", Color::Verified { lo: 0.9, hi: 1.1 });
    let dirty = g.source(
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
    let val = g.source(
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
            ids.push(gg.source(&format!("s{k}"), rand_color(&mut rng)));
        }
        for k in 0..6 {
            let a = ids[rng.below(ids.len() as u64) as usize];
            let b = ids[rng.below(ids.len() as u64) as usize];
            let derived_rank = {
                let (ca, _) = check_regime(&gg.node(a).color.clone(), &state);
                let (cb, _) = check_regime(&gg.node(b).color.clone(), &state);
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
        let dirty = g.source(
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
    let lineage = vec![g.node(dirty).hash];
    let grant = signed_grant(
        0x5EC2E7,
        "release-key-1",
        "release-metric",
        &claimed,
        lineage.clone(),
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
    let node = g.node(id).clone();
    assert!(node.grant.is_some() && node.waiver.is_some());
    assert!(
        g.rows()
            .iter()
            .any(|r| r.contains("release-key-1") && r.contains("\"authorized\":true"))
    );
    // Re-verifiable FROM the ledger: the stored grant's payload still
    // authenticates under the same verifier.
    let stored = node.grant.as_ref().unwrap();
    assert!(verifier.verify(&stored.key_id, &stored.signing_payload(), &stored.signature));
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
    assert_ne!(g2.node(plain).hash.to_hex(), node.hash.to_hex());
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
    assert_eq!(
        refusal(tampered, &verifier, 200),
        WaiverRejection::BadSignature
    );
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
    wrong_lineage.signature = test_mac(0x5EC2E7, &wrong_lineage.signing_payload());
    assert_eq!(
        refusal(wrong_lineage, &verifier, 200),
        WaiverRejection::LineageMismatch
    );
    // Wrong color, wrong scope, expiry.
    let mut wrong_color = grant.clone();
    wrong_color.claimed_color = "validated".to_string();
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
    assert_eq!(
        refusal(wrong_key, &verifier, 200),
        WaiverRejection::BadSignature
    );
    let rotated = TestVerifier {
        keys: KeyMap::from([("release-key-2".to_string(), 0x0FF1CEu64)]),
    };
    assert_eq!(
        refusal(grant.clone(), &rotated, 200),
        WaiverRejection::BadSignature
    );
    // Delimiter injection: adversarial names cannot collide the
    // length-prefixed payload of a DIFFERENT legitimate grant.
    let evil = signed_grant(
        0x5EC2E7,
        "release-key-1",
        "release-metric\ninjected",
        &claimed,
        lineage,
        400,
    );
    assert_ne!(evil.signing_payload(), grant.signing_payload());
    // No-crypto no-claim: the in-tree default verifier refuses all.
    assert_eq!(
        refusal(grant, &NoWaiverVerifier, 200),
        WaiverRejection::BadSignature
    );
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
        let a = g.source("a", Color::Verified { lo: 0.0, hi: 1.0 });
        let b = g.source(
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

/// col-007 — color rows remain strict JSON when fail-closed demotion emits
/// non-finite sentinels and caller-controlled metadata contains JSON syntax or
/// control characters. Validation goes through the ledger's SQLite `json_valid`
/// path, the same parser that enforces persisted payloads.
#[test]
fn col_007_color_rows_are_strict_json_under_hostile_metadata() {
    let build = || -> Vec<String> {
        let hostile = "meta\"\\\n\r\t\u{0007}";
        let axis = format!("Re-{hostile}");
        let dataset = format!("anchors-{hostile}");
        let mut graph = ColorGraph::new();
        let validated = graph.source(
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

        let estimated = graph.source(
            &format!("estimated-{hostile}"),
            Color::Estimated {
                estimator: format!("surrogate-{hostile}"),
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
