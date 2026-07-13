//! EPISTEMIC TYPE-SYSTEM E2E SUITE (bead qmao.10): the Layer-2 battery
//! run as ONE script — laundering, falsifier economy, Goodhart guard,
//! objective epistemics, and the evidence-package round-trip — with
//! every stage logging its enumerated fields as STRUCTURED LEDGER
//! EVENTS (kind `epi-e2e`), not stdout prose. This suite is the
//! artifact of record that the type system FAILS SAFE, not just
//! correct.

use std::collections::BTreeMap;

use fs_evidence::falsify::{ClaimContext, FalsifierHistory, FalsifierRegistry};
use fs_evidence::{Color, IntervalOp, ValidityDomain};
use fs_ledger::{
    ColorGraph, EventRow, Ledger, PolicyDecision, SourceOrigin, SourceOriginRequest,
    SourceOriginVerifier, Waiver, hash_bytes,
};
use fs_opt::{
    DeltaPerturbationStep, Endpoint, EscalationKind, EscalationStep, GoodhartGuard, GuardStatus,
    StepOutcome,
};
use fs_package::{Claim, EvidencePackage, Provenance, SemanticWitness};

struct FixtureSourceVerifier(Vec<u8>);

impl SourceOriginVerifier for FixtureSourceVerifier {
    fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision {
        let fingerprint = hash_bytes(b"fs-checker/epi-ledger-source-policy/v1");
        if self.0 == request.canonical_bytes() {
            PolicyDecision::accept(fingerprint)
        } else {
            PolicyDecision::reject(fingerprint)
        }
    }
}

struct PackageSourceVerifier;
struct PackageSignatureVerifier;

fn drag_interval_witness() -> SemanticWitness {
    let mut payload = Vec::new();
    payload.extend_from_slice(&1_u32.to_le_bytes());
    payload.push(1);
    payload.extend_from_slice(&0.31_f64.to_bits().to_le_bytes());
    payload.extend_from_slice(&0.33_f64.to_bits().to_le_bytes());
    payload.extend_from_slice(&0_u32.to_le_bytes());
    SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, payload)
}

/// Counts dispatches; used to prove the checker NEVER hands
/// attacker-selected certificate addresses to a trusted capability
/// before the package content-address authenticates (bead x2ch).
struct CountingSourceVerifier {
    calls: std::cell::Cell<usize>,
}

impl fs_checker::SourceCertificateVerifier for CountingSourceVerifier {
    fn verify(
        &self,
        _request: &fs_checker::SourceCertificateRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        self.calls.set(self.calls.get() + 1);
        fs_checker::VerificationDecision::reject(hash_bytes(b"counting-policy"))
    }
}

impl fs_checker::SourceCertificateVerifier for PackageSourceVerifier {
    fn verify(
        &self,
        request: &fs_checker::SourceCertificateRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        let witness = drag_interval_witness();
        let accepted = request.package_provenance.code_version == "frankensim@3fab970"
            && request.package_provenance.constellation_lock == "lock-digest-77"
            && request.package_root != fs_checker::ContentHash([0; 32])
            && request.claim_index == 0
            && request.claim_id == "drag"
            && request.statement == "Cd in [0.31, 0.33] at Re 2e5"
            && request.claim_subject_hash != fs_checker::ContentHash([0; 32])
            && request.lo.to_bits() == 0.31f64.to_bits()
            && request.hi.to_bits() == 0.33f64.to_bits()
            && request.producer == "test-solver/cert"
            && request.certificate_hash == witness.content_hash()
            && request.semantic_witness == Some(&witness);
        let fingerprint = hash_bytes(b"fs-checker:epi-package-source-policy:v1");
        if accepted {
            fs_checker::VerificationDecision::accept(fingerprint)
        } else {
            fs_checker::VerificationDecision::reject(fingerprint)
        }
    }
}

impl fs_checker::SignatureVerifier for PackageSignatureVerifier {
    fn verify(
        &self,
        request: &fs_checker::SignatureRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        let fingerprint = hash_bytes(b"fs-checker:epi-package-signature-policy:v1");
        if request.signature
            == format!("fs-checker:epi-subject:{}", request.subject_hash().to_hex())
            && request.purpose == fs_checker::SignaturePurpose::PackageRootAttestation
        {
            fs_checker::VerificationDecision::accept(fingerprint)
        } else {
            fs_checker::VerificationDecision::reject(fingerprint)
        }
    }
}
use fs_robust::{ColoredObjective, RobustError, fragility_curve};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"epi-e2e\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The stage logger: structured ledger events, never stdout.
struct StageLog {
    ledger: Ledger,
    t: i64,
}

impl StageLog {
    fn log(&mut self, stage: &str, payload: &str) {
        self.t += 1;
        self.ledger
            .append_event(&EventRow {
                session: None,
                t: self.t,
                kind: "epi-e2e",
                payload: Some(&format!("{{\"stage\":\"{stage}\",{payload}}}")),
            })
            .expect("event");
    }
}

/// A mock escalation step for ladder slots whose machinery lives in
/// other crates (rung k+1, cross-rep, estimator independence): the e2e
/// battery exercises the LADDER SEMANTICS with controlled outcomes.
struct MockStep {
    kind: EscalationKind,
    outcome: StepOutcome,
}

impl EscalationStep for MockStep {
    fn kind(&self) -> EscalationKind {
        self.kind
    }

    fn evaluate(&self, _endpoint: &Endpoint) -> StepOutcome {
        self.outcome.clone()
    }
}

struct MacVerifier;

fn mac(payload: &[u8]) -> Vec<u8> {
    let mut acc = 0xcbf2_9ce4_8422_2325u64 ^ 0xE2E;
    for &byte in payload {
        acc ^= u64::from(byte);
        acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
    acc.to_le_bytes().to_vec()
}

impl fs_ledger::WaiverVerifier for MacVerifier {
    fn verify(&self, key_id: &str, payload: &[u8], signature: &[u8]) -> PolicyDecision {
        let fingerprint = hash_bytes(b"fs-checker/epi-ledger-waiver-policy/v1");
        if key_id == "epi-key" && mac(payload) == signature {
            PolicyDecision::accept(fingerprint)
        } else {
            PolicyDecision::reject(fingerprint)
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn epi_e2e_battery() {
    let dir = std::env::temp_dir().join(format!("epi-e2e-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let mut log = StageLog {
        ledger: Ledger::open(dir.join("e2e.led").to_str().expect("utf8")).expect("ledger"),
        t: 0,
    };

    // ---- STAGE 1: the laundering battery --------------------------------
    let mut graph = ColorGraph::new();
    let state_in: BTreeMap<String, f64> = [("Re".to_string(), 2.0e5)].into();
    let state_out: BTreeMap<String, f64> = [("Re".to_string(), 5.0e5)].into();
    let surrogate = graph
        .source(
            "drag-surrogate",
            Color::Estimated {
                estimator: "fno-v1".to_string(),
                dispersion: 0.1,
            },
        )
        .expect("surrogate is Estimated");
    let anchor_regime = ValidityDomain::unconstrained().with("Re", 1.0e5, 3.0e5);
    let anchor_color = Color::Validated {
        regime: anchor_regime.clone(),
        dataset: "tunnel-2026".to_string(),
    };
    let anchor_origin = SourceOrigin::Anchoring {
        dataset_id: "tunnel-2026".to_string(),
        content_hash: hash_bytes(b"tunnel-2026 fixture"),
        regime: anchor_regime,
    };
    let anchor_verifier = FixtureSourceVerifier(
        SourceOriginRequest::new("tunnel-anchor", &anchor_color, &anchor_origin).canonical_bytes(),
    );
    let anchored = graph
        .source_with_origin(
            "tunnel-anchor",
            &anchor_color,
            anchor_origin,
            &anchor_verifier,
        )
        .expect("dataset anchor mints Validated");
    // Adversarial upgrade: REFUSED at write time.
    let refusal = graph
        .derive(
            "polished",
            &[surrogate, anchored],
            IntervalOp::Add,
            Some(Color::Verified { lo: 0.0, hi: 1.0 }),
            &state_in,
            None,
        )
        .expect_err("laundering must fail at write time");
    log.log(
        "laundering",
        &format!("\"event\":\"refusal\",\"rule\":\"no-outrank\",\"detail\":\"{refusal}\""),
    );
    // Validated OUT of regime: auto-demotion to estimated.
    let demoted = graph
        .derive(
            "out-of-regime",
            &[anchored],
            IntervalOp::Hull,
            None,
            &state_out,
            None,
        )
        .expect("derivation runs");
    let node = graph.node(demoted).expect("demoted node");
    assert!(
        matches!(node.declared_color_unverified(), Color::Estimated { .. }),
        "regime exit demotes: {:?}",
        node.declared_color_unverified()
    );
    assert_eq!(node.demotions().len(), 1, "the demotion event is recorded");
    log.log(
        "laundering",
        "\"event\":\"auto-demotion\",\"input\":\"validated\",\"state\":\"Re=5e5\",\
         \"output\":\"estimated\"",
    );
    // The waiver path (qmao.1.1): an AUTHENTICATED grant — bound to
    // node, lineage, color, and scope, unexpired, verifier-accepted —
    // authorizes the upgrade; a bare annotation would be refused.
    let claimed_color = Color::Validated {
        regime: ValidityDomain::unconstrained().with("Re", 1.0e5, 3.0e5),
        dataset: "engineer-judgment".to_string(),
    };
    let mut grant = fs_ledger::WaiverGrant {
        annotation: Waiver {
            id: "MEMO-42".to_string(),
            signer: "chief-engineer".to_string(),
            reason: "surrogate validated offline against tunnel run 9".to_string(),
        },
        key_id: "epi-key".to_string(),
        scope: fs_ledger::WAIVER_SCOPE_COLOR_UPGRADE.to_string(),
        node_name: "waived-upgrade".to_string(),
        claimed_color: claimed_color.canonical_bytes(),
        parent_hashes: vec![graph.node(surrogate).expect("surrogate").hash()],
        expires_day: 400,
        signature: Vec::new(),
    };
    grant.signature = mac(&grant.signing_payload(IntervalOp::Hull));
    let waived = graph
        .derive_waived(
            "waived-upgrade",
            &[surrogate],
            IntervalOp::Hull,
            claimed_color,
            &state_in,
            grant,
            &MacVerifier,
            200,
        )
        .expect("an authenticated grant authorizes the upgrade");
    let wnode = graph.node(waived).expect("waived node");
    assert_eq!(
        wnode.waiver().map(|w| w.signer.as_str()),
        Some("chief-engineer"),
        "the waiver travels with the node"
    );
    log.log(
        "laundering",
        "\"event\":\"waiver\",\"signer\":\"chief-engineer\",\"id\":\"MEMO-42\"",
    );
    verdict(
        "stage-1",
        "laundering refused; regime exit auto-demotes; waiver travels signed",
    );

    // ---- STAGE 2: the falsifier economy ---------------------------------
    let registry = FalsifierRegistry::standard();
    let blocked = registry.ship_gate(&["adjoint-gradient", "unpaired-novel-cert"]);
    assert_eq!(blocked.len(), 1, "no falsifier, no ship");
    log.log(
        "falsifier",
        "\"event\":\"ship-gate\",\"blocked\":\"unpaired-novel-cert\"",
    );
    let mut history = FalsifierHistory::new();
    // Seeded claims: high-consequence high-doubt, low-low, and COLD START.
    let claims = [
        ClaimContext {
            class: "conservation".to_string(),
            regime: "Re-2e5".to_string(),
            consequence: 10.0,
        },
        ClaimContext {
            class: "watertightness".to_string(),
            regime: "Re-2e5".to_string(),
            consequence: 1.0,
        },
        ClaimContext {
            class: "brand-new-cert".to_string(),
            regime: "Re-2e5".to_string(),
            consequence: 5.0,
        },
    ];
    let budget = fs_evidence::falsify::allocate_budget(100.0, &claims, &history);
    let get = |c: &str| {
        claims
            .iter()
            .zip(&budget)
            .find(|(cl, _)| cl.class == c)
            .map(|(_, v)| *v)
            .expect("allocated")
    };
    assert!(
        get("brand-new-cert") > get("watertightness"),
        "cold start carries max doubt: {budget:?}"
    );
    log.log(
        "falsifier",
        &format!(
            "\"event\":\"budget\",\"conservation\":{:.2},\"watertightness\":{:.2},\
             \"cold-start\":{:.2}",
            get("conservation"),
            get("watertightness"),
            get("brand-new-cert")
        ),
    );
    // A falsifier HIT: the tombstone + estimator bug report auto-create.
    let (tombstone, bug) = history.record_hit(
        &fs_evidence::falsify::FalsifierHit {
            class: "conservation".to_string(),
            regime: "Re-2e5".to_string(),
            falsifier: "global-flux-audit".to_string(),
            detail: "flux imbalance 3.2e-2 on the independent quadrature".to_string(),
        },
        42.0,
    );
    assert!(tombstone.json.contains("conservation"));
    assert!(bug.json.contains("global-flux-audit"));
    let (hits, spend, _) = history.yield_of("conservation");
    assert_eq!(hits, 1);
    log.log(
        "falsifier",
        &format!(
            "\"event\":\"hit\",\"class\":\"conservation\",\"yield_hits\":{hits},\
             \"spend_s\":{spend}"
        ),
    );
    verdict(
        "stage-2",
        "ship gate blocks; cold start dominates budget; hits mint tombstones",
    );

    // ---- STAGE 3: the Goodhart guard ------------------------------------
    // A seeded discretization EXPLOIT: the reported optimum lives in a
    // crack of the honest objective.
    let honest = |x: &[f64]| {
        let d: f64 = x.iter().map(|v| (v - 0.5).abs()).sum();
        if d < 1e-6 { -50.0 } else { d }
    };
    let full_guard = |obj: fn(&[f64]) -> f64| {
        GoodhartGuard::new()
            .with_step(Box::new(MockStep {
                kind: EscalationKind::RungKPlus1,
                outcome: StepOutcome::Passed,
            }))
            .with_step(Box::new(MockStep {
                kind: EscalationKind::CrossRepresentation,
                outcome: StepOutcome::Passed,
            }))
            .with_step(Box::new(DeltaPerturbationStep::new(0.01, 1e-9, 1.0, obj)))
            .with_step(Box::new(MockStep {
                kind: EscalationKind::EstimatorIndependence,
                outcome: StepOutcome::Passed,
            }))
    };
    let exploit = full_guard(honest).evaluate(&Endpoint::new("exploit", vec![0.5, 0.5], -50.0));
    assert_eq!(exploit.status, GuardStatus::Failed, "the exploit is vetoed");
    assert!(!exploit.findings.is_empty(), "the veto is treasure");
    let smooth = |x: &[f64]| x.iter().map(|v| (v - 0.5) * (v - 0.5)).sum::<f64>();
    let genuine = full_guard(smooth).evaluate(&Endpoint::new("genuine", vec![0.5, 0.5], 0.0));
    assert_eq!(
        genuine.status,
        GuardStatus::Cleared,
        "no false veto on the full ladder"
    );
    assert!(genuine.is_honored());
    // The unavailable-step case: one slot cannot run → PROVISIONAL.
    let partial = GoodhartGuard::new()
        .with_step(Box::new(MockStep {
            kind: EscalationKind::RungKPlus1,
            outcome: StepOutcome::NotPerformed {
                reason: "no coarser rung registered for this kernel".to_string(),
            },
        }))
        .with_step(Box::new(DeltaPerturbationStep::new(
            0.01, 1e-9, 1.0, smooth,
        )))
        .evaluate(&Endpoint::new("partial", vec![0.5, 0.5], 0.0));
    assert_eq!(partial.status, GuardStatus::Provisional);
    assert!(!partial.is_honored(), "provisional is not honored");
    assert!(partial.findings.is_empty(), "and carries no false veto");
    // Catch-rate over a small endpoint population.
    let endpoints = [
        ("e1", -50.0, true),
        ("e2", 0.0, false),
        ("e3", -50.0, true),
        ("e4", 0.0, false),
    ];
    let mut caught = 0usize;
    for (label, obj_val, is_exploit) in &endpoints {
        let g = if *is_exploit {
            full_guard(honest).evaluate(&Endpoint::new(*label, vec![0.5, 0.5], *obj_val))
        } else {
            full_guard(smooth).evaluate(&Endpoint::new(*label, vec![0.5, 0.5], *obj_val))
        };
        if *is_exploit && g.status == GuardStatus::Failed {
            caught += 1;
        }
    }
    assert_eq!(caught, 2, "catch-rate 2/2 on seeded exploits");
    log.log(
        "goodhart",
        "\"event\":\"catch-rate\",\"exploits\":2,\"caught\":2,\"false_vetoes\":0,\
         \"provisional_on_unavailable\":true",
    );
    verdict(
        "stage-3",
        "exploit vetoed, genuine cleared on the full ladder, unavailable step stays provisional",
    );

    // ---- STAGE 4: objective epistemics ----------------------------------
    let uncolored = ColoredObjective::new("naive", vec![1.0, 2.0], vec![]);
    assert!(
        matches!(
            uncolored.headline_color(),
            Err(RobustError::UncoloredObjective { .. })
        ),
        "no optimization against an un-colored objective"
    );
    // A verified solve under an ESTIMATED hazard: weakest input wins.
    let mixed = ColoredObjective::new(
        "bracket-v3",
        vec![1.0, 1.4, 0.9, 2.2],
        vec![
            Color::Verified { lo: 0.9, hi: 1.1 },
            Color::Estimated {
                estimator: "hazard-model-v2".to_string(),
                dispersion: 0.3,
            },
        ],
    );
    let headline = mixed.headline_color().expect("colored");
    assert!(
        matches!(headline, Color::Estimated { .. }),
        "verified solve under estimated hazard = ESTIMATED headline: {headline:?}"
    );
    // The seismic deliverable: a COLORED fragility curve.
    let fragility = fragility_curve(
        &[1.0, 1.2, 1.4, 1.6, 1.8],
        &[0.5, 1.3, 2.0],
        headline.clone(),
    )
    .expect("curve");
    assert_eq!(fragility.curve.len(), 3);
    assert!((fragility.curve[1].prob_failure - 0.4).abs() < 1e-12);
    assert!(matches!(fragility.color, Color::Estimated { .. }));
    log.log(
        "objective",
        "\"event\":\"weakest-input\",\"inputs\":\"verified+estimated\",\
         \"headline\":\"estimated\",\"fragility_points\":3",
    );
    verdict(
        "stage-4",
        "uncolored refused; weakest-input headline; colored fragility curve",
    );

    // ---- STAGE 5: the evidence-package round-trip -----------------------
    let unsigned = EvidencePackage::new(Provenance::new("frankensim@3fab970", "lock-digest-77"))
        .with_claim(Claim::from_portable_certificate(
            "drag",
            "Cd in [0.31, 0.33] at Re 2e5",
            0.31,
            0.33,
            "test-solver/cert",
            drag_interval_witness(),
        ))
        .with_claim(Claim::estimated(
            "fatigue",
            "life > 1e7 cycles per hazard-model-v2",
            "hazard-model-v2",
            0.2,
        ));
    let root = unsigned.try_merkle_root().expect("bounded fixture root");
    let signature_subject = fs_checker::signature_subject_hash(
        root,
        fs_checker::SignaturePurpose::PackageRootAttestation,
    );
    let package = unsigned.signed(format!(
        "fs-checker:epi-subject:{}",
        signature_subject.to_hex()
    ));
    let package_source_verifier = PackageSourceVerifier;
    let package_signature_verifier = PackageSignatureVerifier;
    let capabilities = fs_checker::VerificationCapabilities::deny_all()
        .with_source_certificates(&package_source_verifier);
    let report = fs_checker::check_with_capabilities(
        &package,
        None,
        Some(&package_signature_verifier),
        &capabilities,
    );
    assert!(
        report.passed(),
        "the exact-subject fixture capability accepts without rerunning the solver"
    );
    assert!(matches!(
        report.signature(),
        fs_checker::SignatureStatus::Authenticated(_)
    ));
    assert_eq!(
        report.semantic_status(),
        fs_checker::SemanticStatus::Verified
    );
    assert_eq!(
        report.origin_status(),
        fs_checker::OriginStatus::Authenticated
    );
    let json_report = fs_checker::check_json_with_capabilities(
        &package.to_json().expect("portable package JSON"),
        Some(root),
        Some(&package_signature_verifier),
        &capabilities,
    );
    assert_eq!(
        json_report, report,
        "portable witness survives JSON round-trip"
    );
    let pie = report.render_pie();
    assert!(
        pie.contains("verified") && pie.contains("estimated"),
        "the pie renders: {pie}"
    );
    // TAMPER: polish the estimated claim after signing.
    let tampered = EvidencePackage::new(Provenance::new("frankensim@3fab970", "lock-digest-77"))
        .with_claim(Claim::from_certificate(
            "drag",
            "Cd in [0.31, 0.33] at Re 2e5",
            0.31,
            0.33,
            "test-solver/cert",
            "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
        ))
        .with_claim(Claim::from_certificate(
            "fatigue",
            "life > 1e7 cycles per hazard-model-v2",
            0.0,
            1.0,
            "tampered/certificate",
            "abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789",
        ))
        .signed(format!(
            "fs-checker:epi-subject:{}",
            signature_subject.to_hex()
        ));
    // DELIBERATE CONTRACT UPDATE (bead x2ch): before the authenticated
    // source-certificate-refused finding here — which required the
    // checker to hand the ATTACKER-SELECTED certificate hash to the
    // trusted capability before the package authenticated. The checker
    // preflight boundary (correctly) refuses to dispatch any injected
    // capability on unauthenticated bytes, so the contract is now the
    // STRONGER one: tamper fails closed on the content root alone, and
    // the capability is provably NEVER dispatched.
    let counting_verifier = CountingSourceVerifier {
        calls: std::cell::Cell::new(0),
    };
    let counting_capabilities = fs_checker::VerificationCapabilities::deny_all()
        .with_source_certificates(&counting_verifier);
    let bad = fs_checker::check_with_capabilities(
        &tampered,
        Some(root),
        Some(&package_signature_verifier),
        &counting_capabilities,
    );
    assert!(!bad.passed(), "tampering fails closed");
    assert!(
        bad.findings()
            .iter()
            .any(|finding| finding.kind == "content-address-mismatch"),
        "claim substitution must break the expected content root: {:?}",
        bad.findings()
    );
    assert!(
        !bad.findings()
            .iter()
            .any(|finding| finding.kind == "source-certificate-refused"),
        "per-claim certificate authentication must NOT run before the \
         content address authenticates (capability trust boundary): {:?}",
        bad.findings()
    );
    assert_eq!(
        counting_verifier.calls.get(),
        0,
        "the attacker-selected certificate address was dispatched to a \
         trusted capability on unauthenticated bytes"
    );
    let named = bad.findings().iter().any(|f| f.detail.contains("root"));
    assert!(
        named,
        "the failure is localized to the root: {:?}",
        bad.findings()
    );
    log.log(
        "package",
        "\"event\":\"round-trip\",\"claims\":2,\"reverified\":true,\
         \"portable_semantics\":\"verified\",\"tamper_localized\":true",
    );
    verdict(
        "stage-5",
        "round-trip re-verifies solver-free; tamper fails closed, localized",
    );

    // ---- The suite's own artifact: the structured log ------------------
    let events = log.ledger.table_count("events").expect("count");
    assert_eq!(events, 9, "every stage logged its enumerated fields");
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "epi-e2e",
        "all five stages green in one script; 9 structured ledger events; the type \
         system fails safe",
    );
}
