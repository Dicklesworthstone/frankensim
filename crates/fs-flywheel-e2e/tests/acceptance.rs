//! EPISTEMIC-ENGINE ACCEPTANCE (bead xpck.8): the TOP-LEVEL runnable
//! integration gate for the whole addendum — a declarative query becomes a
//! colored, priced, auditable answer and crosses the solver-free checker's
//! typed capability boundary. The in-test certificate and signature policies
//! are deterministic fixtures, not vendor-independent scientific or
//! cryptographic verification.
//!
//! The path: admission (typed, teaching refusals) → flywheel discharge
//! (planner + cache) → anytime colored answer (+ the VoI-priced hint)
//! → fixture-authenticated evidence package → SOLVER-FREE policy re-check →
//! G5 whole-path replay → the laundering invariant at every hop.
#![cfg(feature = "flywheel-e2e")]

use fs_evidence::{Color, IntervalOp, compose};
use fs_ir::planner::{CostTable, MemCache, PlanError, ProblemFamily};
use fs_ir::{admission, sexpr};
use fs_package::{Claim, EvidencePackage, Provenance};

const CERTIFICATE_HASH: &str = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

// These fixtures prove exact request binding, policy plumbing, and fail-closed
// substitution behavior. They do not resolve retained certificate artifacts
// or establish an independent trust root.

struct ExactCertificate {
    provenance: Provenance,
    claim_index: usize,
    claim_id: String,
    statement: String,
    lo: f64,
    hi: f64,
    producer: String,
}

impl ExactCertificate {
    fn matches(&self, request: &fs_checker::SourceCertificateRequest<'_>) -> bool {
        request.package_provenance == &self.provenance
            && request.claim_index == self.claim_index
            && request.claim_id == self.claim_id
            && request.statement == self.statement
            && request.lo.to_bits() == self.lo.to_bits()
            && request.hi.to_bits() == self.hi.to_bits()
            && request.producer == self.producer
            && request.certificate_hash.to_hex() == CERTIFICATE_HASH
    }
}

struct ExactCertificateVerifier(Vec<ExactCertificate>);

impl ExactCertificateVerifier {
    fn policy_fingerprint(&self) -> fs_checker::ContentHash {
        fn push(bytes: &mut Vec<u8>, value: &[u8]) {
            let len = u64::try_from(value.len()).expect("fixture field length fits u64");
            bytes.extend_from_slice(&len.to_le_bytes());
            bytes.extend_from_slice(value);
        }

        let mut bytes = b"fs-flywheel-e2e:exact-certificate-policy:v1".to_vec();
        for certificate in &self.0 {
            push(&mut bytes, certificate.provenance.code_version.as_bytes());
            push(
                &mut bytes,
                certificate.provenance.constellation_lock.as_bytes(),
            );
            let claim_index =
                u64::try_from(certificate.claim_index).expect("fixture claim index fits u64");
            bytes.extend_from_slice(&claim_index.to_le_bytes());
            push(&mut bytes, certificate.claim_id.as_bytes());
            push(&mut bytes, certificate.statement.as_bytes());
            bytes.extend_from_slice(&certificate.lo.to_bits().to_le_bytes());
            bytes.extend_from_slice(&certificate.hi.to_bits().to_le_bytes());
            push(&mut bytes, certificate.producer.as_bytes());
            push(&mut bytes, CERTIFICATE_HASH.as_bytes());
        }
        fs_ledger::hash_bytes(&bytes)
    }
}

impl fs_checker::SourceCertificateVerifier for ExactCertificateVerifier {
    fn verify(
        &self,
        request: &fs_checker::SourceCertificateRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        let fingerprint = self.policy_fingerprint();
        if self
            .0
            .iter()
            .any(|certificate| certificate.matches(request))
        {
            fs_checker::VerificationDecision::accept(fingerprint)
        } else {
            fs_checker::VerificationDecision::reject(fingerprint)
        }
    }
}

fn certificate_claim(
    provenance: &Provenance,
    claim_index: usize,
    claim_id: impl Into<String>,
    statement: impl Into<String>,
    lo: f64,
    hi: f64,
    producer: impl Into<String>,
) -> (Claim, ExactCertificate) {
    let claim_id = claim_id.into();
    let statement = statement.into();
    let producer = producer.into();
    let certificate = ExactCertificate {
        provenance: provenance.clone(),
        claim_index,
        claim_id: claim_id.clone(),
        statement: statement.clone(),
        lo,
        hi,
        producer: producer.clone(),
    };
    let claim = Claim::from_certificate(claim_id, statement, lo, hi, producer, CERTIFICATE_HASH);
    (claim, certificate)
}

struct ExactRootSignatureVerifier {
    domain: &'static str,
}

impl fs_checker::SignatureVerifier for ExactRootSignatureVerifier {
    fn verify(
        &self,
        request: &fs_checker::SignatureRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        let fingerprint = fs_ledger::hash_bytes(
            format!("fs-flywheel-e2e:signature-policy:v1:{}", self.domain).as_bytes(),
        );
        if request.signature == format!("{}:{}", self.domain, request.subject_hash().to_hex())
            && request.purpose == fs_checker::SignaturePurpose::PackageRootAttestation
        {
            fs_checker::VerificationDecision::accept(fingerprint)
        } else {
            fs_checker::VerificationDecision::reject(fingerprint)
        }
    }
}

fn signed_fixture(package: EvidencePackage, domain: &'static str) -> EvidencePackage {
    let root = package.try_merkle_root().expect("bounded fixture root");
    let subject = fs_checker::signature_subject_hash(
        root,
        fs_checker::SignaturePurpose::PackageRootAttestation,
    );
    package.signed(format!("{domain}:{}", subject.to_hex()))
}

/// A deliberately corrupted content root (one byte flipped): the v4
/// 32-byte replacement for the old `root ^ 0xdead` tamper idiom.
fn flip(root: fs_package::ContentHash) -> fs_package::ContentHash {
    let mut bytes = *root.as_bytes();
    bytes[0] ^= 0xde;
    fs_package::ContentHash(bytes)
}
use fs_verify::fem1d::Poly;
use std::collections::BTreeMap;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-flywheel-e2e/acceptance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The user's wedge query as a declarative study (the CHT wedge with a
/// priced budget — "answer for under the budget or teach me why not").
const WEDGE: &str = r#"(study "cht-wedge-acceptance"
  (seed 0x5EED0008) (versions (constellation :lock "2026-07"))
  (budget (wall 1h) (mem 32GiB) (qoi-rel-error 1e-2))
  (let wedge (frep (revolve (cheb-profile "wedge.chb"))))
  (let field (flux.free-surface-lbm wedge
               (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
               (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth field :at lip :modes (1 .. 4))))
  (ascent.optimize J :over wedge :method (lbfgs :m 7)
    :until (any (grad-norm 1e-4) (budget-exhausted))
    :emit (ledger report)))"#;

/// The ill-posed variant: the same study demanding wall and memory far
/// beyond the session capability token (capability infeasibility — the
/// check that needs no cost model).
const WEDGE_ILL: &str = r#"(study "cht-wedge-illposed"
  (seed 0x5EED0009) (versions (constellation :lock "2026-07"))
  (capability :wall 100h :mem 512GiB)
  (budget (wall 1h) (mem 32GiB) (qoi-rel-error 1e-2))
  (let wedge (frep (revolve (cheb-profile "wedge.chb"))))
  (let field (flux.free-surface-lbm wedge
               (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
               (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth field :at lip :modes (1 .. 4))))
  (ascent.optimize J :over wedge :method (lbfgs :m 7)
    :until (any (grad-norm 1e-4) (budget-exhausted))
    :emit (ledger report)))"#;

fn admission_cx() -> admission::AdmissionContext<'static> {
    let token = admission::SessionCapability {
        ops: vec![
            "flux.*".to_string(),
            "ascent.*".to_string(),
            "frep".to_string(),
            "xform.*".to_string(),
        ],
        cores: 32,
        mem_bytes: 64 * 1024 * 1024 * 1024,
        wall_s: 7200.0,
    };
    admission::AdmissionContext {
        router: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(token),
        regime: None,
        regime_policy: admission::RegimePolicy::Warn,
    }
}

fn steep_family() -> Result<ProblemFamily, PlanError> {
    let mut c = vec![0.0; 6];
    c[1] = 0.2;
    c[2] = -0.2;
    c[4] = 1.0;
    c[5] = -1.0;
    ProblemFamily::new(Poly(c), "cht-wedge-acceptance")
}

const RUNGS: [usize; 4] = [12, 24, 48, 96];

#[test]
fn ac_001_admission_teaches_in_milliseconds() {
    let cx = admission_cx();
    let node = sexpr::parse(WEDGE).expect("well-formed study");
    let start = std::time::Instant::now();
    let report = admission::admit(&node, &cx);
    let ok_ms = start.elapsed().as_millis();
    assert!(report.admitted, "the well-posed wedge query is admitted");
    // The ILL-POSED variant is refused fast, with ranked teaching fixes.
    let bad = sexpr::parse(WEDGE_ILL).expect("parses");
    let start = std::time::Instant::now();
    let refusal = admission::admit(&bad, &cx);
    let bad_ms = start.elapsed().as_millis();
    assert!(!refusal.admitted, "the infeasible budget is refused");
    let has_fix = refusal.findings.iter().any(|f| !f.fixes.is_empty());
    assert!(has_fix, "the refusal carries ranked fixes (teaching)");
    println!(
        "{{\"metric\":\"admission\",\"ok_ms\":{ok_ms},\"refusal_ms\":{bad_ms},\
         \"findings\":{}}}",
        refusal.findings.len()
    );
    assert!(bad_ms < 100, "refusal in milliseconds: {bad_ms} ms");
    verdict(
        "ac-001",
        "the wedge query admits; the infeasible variant refuses in milliseconds with \
         ranked teaching fixes",
    );
}

#[test]
fn ac_002_flywheel_discharge_and_anytime_answer() -> Result<(), PlanError> {
    use fs_ir::anytime::run_anytime;
    let family = steep_family()?;
    let tol = 6e-3;
    let ladder = [30.0, 90.0, 400.0];
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0)?;
    let report = run_anytime(&family, 1.0, tol, &ladder, &RUNGS, &mut cache, &mut costs)?;
    // ANYTIME: an immediate colored interval that tightens.
    assert!(!report.trajectory.is_empty(), "an immediate answer exists");
    for step in &report.trajectory {
        assert!(
            matches!(step.color, Color::Verified { .. }),
            "every step is a CERTIFIED interval"
        );
    }
    for w in report.trajectory.windows(2) {
        assert!(w[1].bound <= w[0].bound + 1e-12, "monotone tightening");
    }
    // The flywheel reuse: the repeat query discharges from cache at the
    // smallest budget (the cheap-query loop).
    let again = run_anytime(&family, 1.0, tol, &[5.0], &RUNGS, &mut cache, &mut costs)?;
    assert!(
        again.refusal.is_none() && again.trajectory.last().expect("step").discharged,
        "the repeat query is a cache hit within a 5-cell budget"
    );
    // TEACHING REFUSAL: an impossible tolerance returns the achieved
    // interval, the priced gap, and the no-point-estimate clause.
    let mut cache2 = MemCache::default();
    let mut costs2 = CostTable::new(200.0)?;
    let refused = run_anytime(
        &family,
        1.0,
        1e-9,
        &[60.0],
        &RUNGS,
        &mut cache2,
        &mut costs2,
    )?;
    let note = refused.refusal.expect("the refusal note");
    assert!(
        note.contains("achieved a certified") && note.contains("No best-effort point estimate"),
        "the refusal teaches: {note}"
    );
    // A budget below the first solve produces the distinct no-answer
    // refusal. It must not fabricate a trajectory interval or evidence color.
    let mut cache3 = MemCache::default();
    let mut costs3 = CostTable::new(200.0)?;
    let no_answer = run_anytime(&family, 1.0, tol, &[1.0], &RUNGS, &mut cache3, &mut costs3)?;
    let no_answer_note = no_answer
        .refusal
        .expect("an unfunded valid query carries a teaching refusal");
    assert!(
        no_answer.trajectory.is_empty()
            && no_answer_note.contains("without a certified interval")
            && no_answer_note.contains("No best-effort point estimate"),
        "no budget means no fabricated answer: {no_answer_note}"
    );
    println!(
        "{{\"metric\":\"anytime\",\"steps\":{},\"final_bound\":{:.3e},\
         \"cache_hit_budget\":5}}",
        report.trajectory.len(),
        report.trajectory.last().expect("step").bound
    );
    verdict(
        "ac-002",
        "immediate certified interval, monotone tightening, cache-hit repeat at a 5-cell \
         budget, plus distinct best-certified and no-certified-answer teaching refusals",
    );
    Ok(())
}

#[test]
#[allow(clippy::too_many_lines)] // one auditable query-to-package fixture
fn ac_003_package_recheck_solver_free_and_voi_hint() -> Result<(), PlanError> {
    use fs_ir::anytime::run_anytime;
    use fs_plan::voi::{LiveDecision, Probe, ProbeKind, UncertaintyNode, rank_purchases};
    // Discharge the query, wrap the answer, and exercise the STANDALONE
    // checker's certificate-policy, composition, and content-address paths
    // with zero solver dependency.
    let family = steep_family()?;
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0)?;
    let report = run_anytime(&family, 1.0, 6e-3, &[400.0], &RUNGS, &mut cache, &mut costs)?;
    let last = report.trajectory.last().expect("answer");
    let bound = last.bound;
    // The VoI-priced hint (Proposal C riding the answer).
    let margin = move |v: &[f64]| v[0] - 5e-3;
    let decision = LiveDecision {
        margin: &margin,
        arity: 1,
    };
    let nodes = vec![UncertaintyNode {
        name: "qoi-bound".to_string(),
        lo: 0.0,
        hi: bound.max(1e-6) * 2.0,
        nominal: bound,
    }];
    let menu = vec![Probe {
        name: "climb-final-rung".to_string(),
        target: "qoi-bound".to_string(),
        cost: 12.0,
        shrink: 0.25,
        kind: ProbeKind::Computational,
    }];
    let ranked = rank_purchases(&decision, &nodes, &menu, 32).expect("valid bounded VoI request");
    let hint = fs_plan::voi::hint_for_query(&ranked).expect("valid VoI hint");
    // The package: colored claims, fixture-authenticated, Merkle-rooted.
    let Color::Verified { lo, hi } = last.color else {
        panic!("the wedge trajectory ends verified");
    };
    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    let (qoi_claim, qoi_certificate) = certificate_claim(
        &provenance,
        0,
        "wedge-qoi-interval",
        format!("certified half-width {bound:.3e} at tol 6e-3"),
        lo,
        hi,
        "fs-wedge/dwr-certifier",
    );
    let pkg = signed_fixture(
        EvidencePackage::new(provenance)
            .with_claim(qoi_claim)
            .with_claim(Claim::estimated(
                "voi-hint",
                hint.clone(),
                "voi-myopic",
                1.0,
            )),
        "acceptance-gate",
    );
    // SOLVER-FREE FIXTURE RE-CHECK: fs-checker exercises exact certificate,
    // signature-policy, composition, and root bindings without solver deps.
    let source_verifier = ExactCertificateVerifier(vec![qoi_certificate]);
    let signature_verifier = ExactRootSignatureVerifier {
        domain: "acceptance-gate",
    };
    let capabilities =
        fs_checker::VerificationCapabilities::deny_all().with_source_certificates(&source_verifier);
    let check =
        fs_checker::check_with_capabilities(&pkg, None, Some(&signature_verifier), &capabilities);
    assert!(
        check.passed(),
        "the package passes the solver-free fixture policy"
    );
    assert!(matches!(
        check.signature(),
        fs_checker::SignatureStatus::Authenticated(_)
    ));
    let root = pkg.try_merkle_root().expect("bounded fixture root");
    assert!(
        fs_checker::check_with_capabilities(
            &pkg,
            Some(root),
            Some(&signature_verifier),
            &capabilities,
        )
        .passed(),
        "the content address matches"
    );
    assert!(
        !fs_checker::check_with_capabilities(
            &pkg,
            Some(flip(root)),
            Some(&signature_verifier),
            &capabilities,
        )
        .passed(),
        "a tampered root fails the independent checker code path"
    );
    let pie = check.render_pie();
    println!(
        "{{\"metric\":\"package\",\"root\":\"{root}\",\"hint\":\"{hint}\",\
         \"pie\":\"{}\"}}",
        pie.replace('"', "'").replace('\n', " | ")
    );
    verdict(
        "ac-003",
        "the colored answer enters a fixture-authenticated Merkle-rooted package carrying \
         its VoI-priced hint; the standalone checker exercises solver-free capability \
         binding and catches a tampered root",
    );
    Ok(())
}

#[test]
fn ac_004_g5_whole_path_replay() -> Result<(), PlanError> {
    use fs_ir::anytime::run_anytime;
    let run = || -> Result<(Vec<u64>, fs_package::ContentHash), PlanError> {
        let family = steep_family()?;
        let mut cache = MemCache::default();
        let mut costs = CostTable::new(200.0)?;
        let report = run_anytime(
            &family,
            1.0,
            6e-3,
            &[30.0, 400.0],
            &RUNGS,
            &mut cache,
            &mut costs,
        )?;
        let bits: Vec<u64> = report
            .trajectory
            .iter()
            .map(|s| s.bound.to_bits())
            .collect();
        let pkg = signed_fixture(
            EvidencePackage::new(Provenance::new("acceptance-e2e", "Cargo.lock")).with_claim({
                let Color::Verified { lo, hi } =
                    report.trajectory.last().expect("step").color.clone()
                else {
                    panic!("the replay trajectory ends verified");
                };
                Claim::from_certificate(
                    "wedge-qoi-interval",
                    "replay claim",
                    lo,
                    hi,
                    "fs-wedge/dwr-certifier",
                    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
                )
            }),
            "acceptance-gate/ac-004",
        );
        Ok((bits, pkg.try_merkle_root().expect("bounded fixture root")))
    };
    let (bits_a, root_a) = run()?;
    let (bits_b, root_b) = run()?;
    assert_eq!(bits_a, bits_b, "the trajectory replays bit-exact");
    assert_eq!(root_a, root_b, "the artifact hash replays exactly");
    verdict(
        "ac-004",
        "the whole path — discharge, trajectory, package root — replays bit-exact (G5)",
    );
    Ok(())
}

#[test]
fn ac_005_laundering_invariant_across_the_path() {
    // An ESTIMATED intermediate anywhere in the composition can never
    // surface as VERIFIED — checked at the color algebra AND at the
    // package layer.
    let estimated = Color::Estimated {
        estimator: "dwr-guess".to_string(),
        dispersion: 0.1,
    };
    let verified = Color::Verified { lo: 1.0, hi: 1.1 };
    let composed = compose(&verified, &estimated, IntervalOp::Add);
    assert!(
        !matches!(composed, Color::Verified { .. }),
        "weakest-input: verified x estimated is NOT verified: {composed:?}"
    );
    // At the package layer the breakdown keeps them apart — an audit
    // sees exactly how much of the answer is estimated.
    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    let (hard_claim, hard_certificate) = certificate_claim(
        &provenance,
        0,
        "hard",
        "certified part",
        1.0,
        1.1,
        "test-solver/cert",
    );
    let pkg = signed_fixture(
        EvidencePackage::new(provenance)
            .with_claim(hard_claim)
            .with_claim(Claim::estimated("soft", "estimated part", "dwr-guess", 0.1)),
        "acceptance-gate/ac-005",
    );
    let source_verifier = ExactCertificateVerifier(vec![hard_certificate]);
    let signature_verifier = ExactRootSignatureVerifier {
        domain: "acceptance-gate/ac-005",
    };
    let capabilities = fs_checker::VerificationCapabilities::deny_all()
        .with_source_certificates(&source_verifier)
        .with_signatures(&signature_verifier);
    let breakdown = pkg
        .color_breakdown_with(&capabilities)
        .expect("the exact source certificate authenticates");
    assert!(
        breakdown.verified == 1 && breakdown.estimated == 1,
        "the package cannot blur colors: {breakdown:?}"
    );
    verdict(
        "ac-005",
        "estimated inputs never launder to verified — enforced by the compose algebra \
         and visible in the package breakdown",
    );
}
