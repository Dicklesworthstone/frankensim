//! EPISTEMIC-ENGINE ACCEPTANCE (bead xpck.8): the TOP-LEVEL runnable
//! integration gate for the whole addendum — a declarative query produces a
//! colored, priced, auditable package and crosses the solver-free checker's
//! typed capability boundary. Scientific promotion consumes an immutable
//! `fs-verify` receipt and independently replays its exact problem/candidate;
//! only the package-root signature policy remains a deterministic fixture, not
//! vendor-independent cryptographic authentication. The retained certificate
//! covers the exact 1-D manufactured proxy; the physical wedge QoI stays an
//! explicit Estimated no-claim until its own upstream certifier exists.
//!
//! The path: admission (typed, teaching refusals) → flywheel discharge
//! (planner + cache) → anytime colored answer (+ the VoI-priced hint)
//! → authentic-receipt-backed evidence package → SOLVER-FREE policy re-check →
//! G5 whole-path replay → the laundering invariant at every hop.
#![cfg(feature = "flywheel-e2e")]

use fs_evidence::{Color, IntervalOp, compose};
use fs_ir::planner::{
    AnswerCache, CachedAnswer, CostTable, MemCache, PlanError, PlanOutcome, ProblemFamily, plan,
};
use fs_ir::{admission, sexpr};
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_verify::estimator::{
    AdmittedVerifierReceipt, PresentedVerifierReceipt, VERIFIER_RECEIPT_QOI,
    VERIFIER_RECEIPT_UNITS, VerifierProducerSourceIdentity, VerifierReceipt,
    admit_verifier_receipt,
};
use fs_verify::fem1d::Poly;

const QOI_ID: &str = VERIFIER_RECEIPT_QOI;
const QOI_UNITS: &str = VERIFIER_RECEIPT_UNITS;
const PHYSICAL_QOI_ID: &str = "cht-wedge-perturbation-growth-min";
const PHYSICAL_QOI_UNITS: &str = "UNRESOLVED:no-authoritative-unit-schema";
const RECEIPT_POLICY: &str = "fs-flywheel-e2e/authentic-receipt-resolver/v1";

fn push_bytes(bytes: &mut Vec<u8>, value: &[u8]) {
    let len = u64::try_from(value.len()).expect("bounded acceptance field length fits u64");
    bytes.extend_from_slice(&len.to_le_bytes());
    bytes.extend_from_slice(value);
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_bytes(bytes, value.as_bytes());
}

fn push_hash(bytes: &mut Vec<u8>, value: fs_checker::ContentHash) {
    bytes.extend_from_slice(value.as_bytes());
}

fn hash_f64_slice(domain: &str, values: &[f64]) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(&mut bytes, domain);
    bytes.extend_from_slice(
        &u64::try_from(values.len())
            .expect("bounded acceptance vector length fits u64")
            .to_le_bytes(),
    );
    for value in values {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    fs_ledger::hash_bytes(&bytes)
}
#[derive(Debug, Clone)]
struct AuthenticReceiptResolver {
    receipt: PresentedVerifierReceipt,
    family: ProblemFamily,
    theta: f64,
    tolerance: f64,
    planner_budget: f64,
    planner_spent: f64,
    rungs: Vec<usize>,
}

struct ResolvedVerifierReceipt<'a> {
    admitted: AdmittedVerifierReceipt<'a>,
}

fn retain_original_receipt(receipt: VerifierReceipt) -> PresentedVerifierReceipt {
    let bytes = receipt
        .canonical_bytes()
        .expect("bounded original receipt transport bytes");
    let presented = PresentedVerifierReceipt::from_retained_bytes(&bytes, receipt.artifact_root())
        .expect("independently rooted original receipt transport");
    assert_eq!(
        presented.artifact_root(),
        receipt.artifact_root(),
        "retained transport must preserve the original verifier receipt identity"
    );
    presented
}

impl ResolvedVerifierReceipt<'_> {
    fn receipt(&self) -> &AdmittedVerifierReceipt<'_> {
        &self.admitted
    }

    fn claim(&self) -> Claim {
        let receipt = self.receipt();
        Claim::from_certificate(
            receipt.qoi(),
            receipt.statement(),
            receipt.bound_lo(),
            receipt.bound_hi(),
            receipt.producer().label(),
            receipt.artifact_root().to_hex(),
        )
    }
}

impl AuthenticReceiptResolver {
    fn capture(
        family: &ProblemFamily,
        theta: f64,
        tolerance: f64,
        planner_budget: f64,
        rungs: &[usize],
    ) -> Result<Self, PlanError> {
        let mut cache = MemCache::default();
        let mut costs = CostTable::new(200.0)?;
        let outcome = plan(
            family,
            theta,
            tolerance,
            planner_budget,
            rungs,
            &mut cache,
            &mut costs,
        )?;
        let (candidate, mesh, spent, receipt) = outcome_candidate(outcome)
            .expect("an acceptance receipt budget must produce a bounded candidate");
        let receipt = retain_original_receipt(receipt);
        let problem = family.at(theta, mesh)?;
        admit_verifier_receipt(&problem, &candidate, tolerance, &receipt)
            .expect("the retained original planner receipt must replay before use");
        Ok(Self {
            receipt,
            family: family.clone(),
            theta,
            tolerance,
            planner_budget,
            planner_spent: spent,
            rungs: rungs.to_vec(),
        })
    }

    fn from_original_step(
        family: &ProblemFamily,
        theta: f64,
        tolerance: f64,
        planner_budget: f64,
        planner_spent: f64,
        rungs: &[usize],
        receipt: VerifierReceipt,
    ) -> Self {
        Self {
            receipt: retain_original_receipt(receipt),
            family: family.clone(),
            theta,
            tolerance,
            planner_budget,
            planner_spent,
            rungs: rungs.to_vec(),
        }
    }

    fn from_presented_step(
        family: &ProblemFamily,
        theta: f64,
        tolerance: f64,
        planner_budget: f64,
        planner_spent: f64,
        rungs: &[usize],
        receipt: PresentedVerifierReceipt,
    ) -> Self {
        Self {
            receipt,
            family: family.clone(),
            theta,
            tolerance,
            planner_budget,
            planner_spent,
            rungs: rungs.to_vec(),
        }
    }

    fn resolve(&self) -> Result<ResolvedVerifierReceipt<'_>, &'static str> {
        let mut cache = MemCache::default();
        let mut costs = CostTable::new(200.0).map_err(|_| "replay cost table refused")?;
        let outcome = plan(
            &self.family,
            self.theta,
            self.tolerance,
            self.planner_budget,
            &self.rungs,
            &mut cache,
            &mut costs,
        )
        .map_err(|_| "independent planner replay refused")?;
        let (candidate, mesh, spent, replay_receipt) =
            outcome_candidate(outcome).ok_or("independent replay produced no candidate")?;
        if spent.to_bits() != self.planner_spent.to_bits() {
            return Err("independent planner consumption differs from the original step");
        }
        if replay_receipt.artifact_root() != self.receipt.artifact_root() {
            return Err("original planner receipt differs from deterministic planner replay");
        }
        let problem = self
            .family
            .at(self.theta, mesh)
            .map_err(|_| "independent problem replay refused")?;
        let admitted = admit_verifier_receipt(&problem, &candidate, self.tolerance, &self.receipt)
            .map_err(|_| "original verifier receipt failed exact production replay")?;
        Ok(ResolvedVerifierReceipt { admitted })
    }
}

fn outcome_candidate(outcome: PlanOutcome) -> Option<(Vec<f64>, Vec<f64>, f64, VerifierReceipt)> {
    match outcome {
        PlanOutcome::Discharged {
            nodal,
            mesh,
            certificate,
            cost,
            ..
        } => Some((nodal, mesh, cost, certificate.receipt().clone())),
        PlanOutcome::RefusedWithBest {
            best_nodal,
            best_mesh,
            best_certificate,
            cost,
            ..
        } => Some((
            best_nodal,
            best_mesh,
            cost,
            best_certificate.receipt().clone(),
        )),
        PlanOutcome::RefusedWithoutAnswer { .. } => None,
    }
}

enum ReceiptPromotion<'a> {
    Verified(ResolvedVerifierReceipt<'a>),
    Gated {
        color: Color,
        no_claim: &'static str,
    },
}

fn resolve_for_promotion(resolver: Option<&AuthenticReceiptResolver>) -> ReceiptPromotion<'_> {
    let Some(resolver) = resolver else {
        return ReceiptPromotion::Gated {
            color: Color::Estimated {
                estimator: "missing-authentic-verifier-receipt".to_string(),
                dispersion: 1.0,
            },
            no_claim: "NO-CLAIM: upstream verifier receipt authority is unavailable",
        };
    };
    match resolver.resolve() {
        Ok(resolved) if resolved.receipt().accepted() => ReceiptPromotion::Verified(resolved),
        Ok(resolved) => ReceiptPromotion::Gated {
            color: Color::Estimated {
                estimator: "verified-bound-above-requested-tolerance".to_string(),
                dispersion: resolved.receipt().bound_hi(),
            },
            no_claim: "NO-CLAIM: authentic bound did not discharge the requested tolerance",
        },
        Err(_) => ReceiptPromotion::Gated {
            color: Color::Estimated {
                estimator: "authentic-verifier-replay-refused".to_string(),
                dispersion: 1.0,
            },
            no_claim: "NO-CLAIM: retained verifier receipt failed independent replay",
        },
    }
}

fn physical_qoi_no_claim(detail: &str, dispersion: f64) -> Claim {
    Claim::estimated(
        PHYSICAL_QOI_ID,
        format!(
            "NO-CLAIM for {PHYSICAL_QOI_ID} ({PHYSICAL_QOI_UNITS}): {detail}; the authentic \
             receipt covers only the named 1-D manufactured proxy"
        ),
        "gated:physical-wedge-certifier-unavailable",
        dispersion,
    )
}

struct ReceiptCertificateVerifier<'a> {
    resolved: ResolvedVerifierReceipt<'a>,
    provenance: Provenance,
    claim_index: usize,
    expected_claim_subject_hash: fs_checker::ContentHash,
}

impl<'a> ReceiptCertificateVerifier<'a> {
    fn from_resolver(
        resolver: &'a AuthenticReceiptResolver,
        provenance: Provenance,
        claim_index: usize,
    ) -> Result<Self, &'static str> {
        let resolved = resolver.resolve()?;
        let expected_claim_subject_hash = resolved
            .claim()
            .declared_source_certificate_subject_hash_unverified();
        Ok(Self {
            resolved,
            provenance,
            claim_index,
            expected_claim_subject_hash,
        })
    }

    fn policy_fingerprint(&self) -> fs_checker::ContentHash {
        let mut bytes = Vec::new();
        push_text(&mut bytes, RECEIPT_POLICY);
        push_hash(
            &mut bytes,
            self.resolved.receipt().artifact_root().content_hash(),
        );
        push_hash(&mut bytes, self.expected_claim_subject_hash);
        push_text(&mut bytes, &self.provenance.code_version);
        push_text(&mut bytes, &self.provenance.constellation_lock);
        bytes.extend_from_slice(
            &u64::try_from(self.claim_index)
                .expect("bounded claim index fits u64")
                .to_le_bytes(),
        );
        fs_ledger::hash_bytes(&bytes)
    }
}

impl fs_checker::SourceCertificateVerifier for ReceiptCertificateVerifier<'_> {
    fn verify(
        &self,
        request: &fs_checker::SourceCertificateRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        let fingerprint = self.policy_fingerprint();
        let receipt = self.resolved.receipt();
        if receipt.accepted()
            && request.package_provenance == &self.provenance
            && request.claim_index == self.claim_index
            && request.claim_subject_hash == self.expected_claim_subject_hash
            && request.claim_id == receipt.qoi()
            && request.statement == receipt.statement()
            && request.lo.to_bits() == receipt.bound_lo().to_bits()
            && request.hi.to_bits() == receipt.bound_hi().to_bits()
            && request.producer == receipt.producer().label()
            && request.certificate_hash == receipt.artifact_root().content_hash()
        {
            fs_checker::VerificationDecision::accept(fingerprint)
        } else {
            fs_checker::VerificationDecision::reject(fingerprint)
        }
    }
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
    let polynomial = Poly::new(c).expect("wedge acceptance polynomial fixture must be admissible");
    ProblemFamily::new(polynomial, "cht-wedge-acceptance")
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
    use fs_plan::voi::{
        Cx, DecisionBudget, LiveDecision, MAX_VOI_EVALUATIONS, MAX_VOI_WORK_UNITS, Probe,
        ProbeKind, UncertaintyNode, rank_purchases,
    };
    // Discharge the query, wrap the answer, and exercise the STANDALONE
    // checker's receipt-policy, composition, and content-address paths. The
    // real verifier runs before the checker receives its sealed capability.
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
    let ranked = rank_purchases(
        &Cx::for_testing(),
        &decision,
        &nodes,
        &menu,
        32,
        DecisionBudget::new(MAX_VOI_EVALUATIONS, MAX_VOI_WORK_UNITS)
            .expect("valid VoI computation budget"),
        "fs-flywheel-acceptance-v1",
        "ac-003-snapshot-v1",
    )
    .expect("valid bounded VoI request");
    let hint = fs_plan::voi::hint_for_query(&ranked);
    let hint_text = hint.render_text();
    // The package: authentic proxy receipt, explicit physical-QoI no-claim,
    // and an Estimated VoI hint, all Merkle-rooted.
    let Color::Verified { lo, hi } = &last.color else {
        panic!("the wedge trajectory ends verified");
    };
    let resolver = AuthenticReceiptResolver::from_original_step(
        &family,
        1.0,
        6e-3,
        last.budget,
        last.spent,
        &RUNGS,
        last.verifier_receipt.clone(),
    );
    let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(Some(&resolver)) else {
        panic!("the authentic final receipt must resolve before package construction");
    };
    let admitted_receipt = resolved.receipt();
    assert_eq!(
        admitted_receipt.bound_lo().to_bits(),
        lo.to_bits(),
        "the package interval is the authentic verifier interval"
    );
    assert_eq!(
        admitted_receipt.bound_hi().to_bits(),
        hi.to_bits(),
        "the package interval is the authentic verifier interval"
    );
    assert_eq!(
        admitted_receipt.flux_hash(),
        last.flux_hash,
        "the package retains the planner's authentic flux reconstruction"
    );
    assert_eq!(
        admitted_receipt.verifier_family(),
        last.verifier_family.as_str(),
        "the package cannot relabel the verifier family"
    );
    assert_eq!(
        resolver.receipt.artifact_root(),
        last.verifier_receipt.artifact_root(),
        "package construction consumes the original anytime receipt, not a fresh mint"
    );
    let qoi_claim = resolved.claim();
    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    let pkg = signed_fixture(
        EvidencePackage::new(provenance.clone())
            .with_claim(qoi_claim)
            .with_claim(physical_qoi_no_claim(
                "no authentic verifier receipt exists for the physical perturbation-growth QoI",
                bound.max(f64::MIN_POSITIVE),
            ))
            .with_claim(Claim::estimated("voi-hint", hint_text, "voi-myopic", 1.0)),
        "acceptance-gate",
    );
    // SOLVER-FREE PACKAGE RE-CHECK: the checker receives only the sealed token
    // produced by the independent replay above. Its injected capability can
    // match the exact request without importing or invoking a solver.
    let source_verifier = ReceiptCertificateVerifier::from_resolver(&resolver, provenance, 0)
        .expect("receipt resolves before the solver-free package check");
    let signature_verifier = ExactRootSignatureVerifier {
        domain: "acceptance-gate",
    };
    let capabilities =
        fs_checker::VerificationCapabilities::deny_all().with_source_certificates(&source_verifier);
    let check =
        fs_checker::check_with_capabilities(&pkg, None, Some(&signature_verifier), &capabilities);
    assert!(
        check.passed(),
        "the package passes the solver-free receipt policy"
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
    match resolve_for_promotion(None) {
        ReceiptPromotion::Gated { color, no_claim } => {
            assert!(matches!(color, Color::Estimated { .. }));
            assert!(no_claim.contains("NO-CLAIM") && no_claim.contains("unavailable"));
        }
        ReceiptPromotion::Verified(_) => {
            panic!("missing upstream authority must never mint a verified claim")
        }
    }
    let pie = check.render_pie();
    println!(
        "{{\"schema_version\":1,\"suite\":\"fs-flywheel-e2e/acceptance\",\
         \"case\":\"ac-003\",\"seed\":\"0x5EED0008\",\"units\":\"{QOI_UNITS}\",\
         \"budget_cells_bits\":\"{:016x}\",\"spent_cells_bits\":\"{:016x}\",\
         \"cancellation\":\"completed-drained-finalized-published\",\
         \"artifact_root\":\"{root}\",\"verifier_family\":\"{}\",\
         \"verifier_receipt\":\"{}\",\"physical_promotion\":\"estimated-no-claim\",\
         \"hint\":{},\"pie\":\"{}\"}}",
        resolver.planner_budget.to_bits(),
        resolver.planner_spent.to_bits(),
        resolver.receipt.verifier_family(),
        resolver.receipt.artifact_root(),
        hint.to_json(),
        pie.replace('"', "'").replace('\n', " | ")
    );
    verdict(
        "ac-003",
        "an immutable fs-verify receipt is independently replayed before its exact 1-D proxy \
         result enters the package; verifier family, flux reconstruction, full receipt address, \
         and VoI hint remain bound, while the physical wedge QoI and missing authority remain \
         explicit Estimated no-claims",
    );
    Ok(())
}

fn package_accepts_claim(
    resolver: &AuthenticReceiptResolver,
    provenance: &Provenance,
    claim: Claim,
) -> bool {
    let package = signed_fixture(
        EvidencePackage::new(provenance.clone()).with_claim(claim),
        "acceptance-gate/ac-003-mutation",
    );
    let source_verifier =
        ReceiptCertificateVerifier::from_resolver(resolver, provenance.clone(), 0)
            .expect("mutation baseline receipt resolves before package checking");
    let signature_verifier = ExactRootSignatureVerifier {
        domain: "acceptance-gate/ac-003-mutation",
    };
    let capabilities =
        fs_checker::VerificationCapabilities::deny_all().with_source_certificates(&source_verifier);
    fs_checker::check_with_capabilities(&package, None, Some(&signature_verifier), &capabilities)
        .passed()
}

#[test]
fn ac_003_production_receipt_and_claim_subject_substitutions_fail_closed() -> Result<(), PlanError>
{
    let family = steep_family()?;
    let resolver = AuthenticReceiptResolver::capture(&family, 1.0, 6e-3, 400.0, &RUNGS)?;
    let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(Some(&resolver)) else {
        panic!("baseline production receipt must promote");
    };
    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    assert!(package_accepts_claim(
        &resolver,
        &provenance,
        resolved.claim()
    ));

    let receipt = resolved.receipt();
    let fake_hash = Claim::from_certificate(
        receipt.qoi(),
        receipt.statement(),
        receipt.bound_lo(),
        receipt.bound_hi(),
        receipt.producer().label(),
        "00".repeat(32),
    );
    assert!(
        !package_accepts_claim(&resolver, &provenance, fake_hash),
        "a fixed/fake certificate address cannot stand in for the production receipt"
    );
    let relabeled = Claim::from_certificate(
        receipt.qoi(),
        receipt.statement(),
        receipt.bound_lo(),
        receipt.bound_hi(),
        "fs-wedge/dwr-certifier",
        receipt.artifact_root().to_hex(),
    );
    assert!(
        !package_accepts_claim(&resolver, &provenance, relabeled),
        "acceptance cannot relabel the production source identity"
    );
    let altered_endpoint = Claim::from_certificate(
        receipt.qoi(),
        receipt.statement(),
        receipt.bound_lo(),
        receipt.bound_hi() + f64::EPSILON,
        receipt.producer().label(),
        receipt.artifact_root().to_hex(),
    );
    assert!(
        !package_accepts_claim(&resolver, &provenance, altered_endpoint),
        "altering an interval endpoint fails exact checker request binding"
    );
    let cross_qoi = Claim::from_certificate(
        format!("{}-foreign", receipt.qoi()),
        receipt.statement(),
        receipt.bound_lo(),
        receipt.bound_hi(),
        receipt.producer().label(),
        receipt.artifact_root().to_hex(),
    );
    assert!(
        !package_accepts_claim(&resolver, &provenance, cross_qoi),
        "the same receipt address cannot authenticate a different QoI subject"
    );

    let cross_problem = AuthenticReceiptResolver::from_presented_step(
        &family,
        0.75,
        resolver.tolerance,
        resolver.planner_budget,
        resolver.planner_spent,
        &resolver.rungs,
        resolver.receipt.clone(),
    );
    assert!(
        cross_problem.resolve().is_err(),
        "an authentic receipt cannot replay against a substituted problem parameter"
    );
    let cross_tolerance = AuthenticReceiptResolver::from_presented_step(
        &family,
        resolver.theta,
        7e-3,
        resolver.planner_budget,
        resolver.planner_spent,
        &resolver.rungs,
        resolver.receipt.clone(),
    );
    assert!(
        cross_tolerance.resolve().is_err(),
        "an authentic receipt cannot replay against a substituted tolerance"
    );

    // Exercise the claim-subject comparison directly: every visible request
    // field remains exact and only the address-free subject digest changes.
    let verifier = ReceiptCertificateVerifier::from_resolver(&resolver, provenance.clone(), 0)
        .expect("baseline production receipt resolves");
    let claim = resolved.claim();
    let statement = receipt.statement();
    let producer = receipt.producer().label();
    let request = fs_checker::SourceCertificateRequest {
        package_provenance: &provenance,
        package_root: fs_checker::ContentHash([7; 32]),
        claim_index: 0,
        claim_id: receipt.qoi(),
        statement: &statement,
        claim_subject_hash: claim.declared_source_certificate_subject_hash_unverified(),
        lo: receipt.bound_lo(),
        hi: receipt.bound_hi(),
        producer: &producer,
        certificate_hash: receipt.artifact_root().content_hash(),
        semantic_witness: None,
    };
    assert!(
        fs_checker::SourceCertificateVerifier::verify(&verifier, &request).accepted(),
        "the exact lower-owned claim subject is recognized"
    );
    let mut changed_subject = *request.claim_subject_hash.as_bytes();
    changed_subject[0] ^= 1;
    let substituted = fs_checker::SourceCertificateRequest {
        claim_subject_hash: fs_checker::ContentHash(changed_subject),
        ..request
    };
    assert!(
        !fs_checker::SourceCertificateVerifier::verify(&verifier, &substituted).accepted(),
        "a subject-only substitution must fail closed"
    );
    Ok(())
}

#[derive(Default)]
struct RecordingCache {
    inner: MemCache,
    committed_mutations: Vec<fs_checker::ContentHash>,
}

impl AnswerCache for RecordingCache {
    fn lookup(&self, key: &str, tolerance: f64) -> Option<CachedAnswer> {
        self.inner.lookup(key, tolerance)
    }

    fn insert(&mut self, key: &str, answer: CachedAnswer) {
        let mut bytes = Vec::new();
        push_text(&mut bytes, "fs-flywheel-e2e/cache-mutation/v1");
        push_text(&mut bytes, key);
        bytes.extend_from_slice(&answer.bound().to_bits().to_le_bytes());
        push_hash(
            &mut bytes,
            hash_f64_slice("cache-candidate/v1", answer.nodal()),
        );
        push_hash(&mut bytes, hash_f64_slice("cache-mesh/v1", answer.mesh()));
        self.committed_mutations.push(fs_ledger::hash_bytes(&bytes));
        self.inner.insert(key, answer);
    }
}

impl RecordingCache {
    fn initial_snapshot_root() -> fs_checker::ContentHash {
        fs_ledger::hash_bytes(b"fs-flywheel-e2e/cache-snapshot/v1:empty")
    }

    fn committed_mutations_root(&self) -> fs_checker::ContentHash {
        let mut bytes = Vec::new();
        push_text(&mut bytes, "fs-flywheel-e2e/cache-mutation-sequence/v1");
        bytes.extend_from_slice(
            &u64::try_from(self.committed_mutations.len())
                .expect("bounded cache mutation count fits u64")
                .to_le_bytes(),
        );
        for mutation in &self.committed_mutations {
            push_hash(&mut bytes, *mutation);
        }
        fs_ledger::hash_bytes(&bytes)
    }
}

fn interval_sequence_root(trajectory: &[fs_ir::anytime::IntervalStep]) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(&mut bytes, "fs-flywheel-e2e/interval-sequence/v1");
    bytes.extend_from_slice(
        &u64::try_from(trajectory.len())
            .expect("bounded interval count fits u64")
            .to_le_bytes(),
    );
    for step in trajectory {
        bytes.extend_from_slice(&step.budget.to_bits().to_le_bytes());
        bytes.extend_from_slice(&step.spent.to_bits().to_le_bytes());
        bytes.extend_from_slice(&step.bound.to_bits().to_le_bytes());
        let Color::Verified { lo, hi } = &step.color else {
            panic!("the acceptance trajectory may retain only verified enclosures");
        };
        bytes.extend_from_slice(&lo.to_bits().to_le_bytes());
        bytes.extend_from_slice(&hi.to_bits().to_le_bytes());
        push_text(&mut bytes, &step.verifier_family);
        bytes.extend_from_slice(&step.flux_hash.to_le_bytes());
        let receipt = &step.verifier_receipt;
        assert_eq!(receipt.bound_hi().to_bits(), step.bound.to_bits());
        assert_eq!(receipt.verifier_family(), step.verifier_family.as_str());
        assert_eq!(receipt.flux_hash(), step.flux_hash);
        for root in [
            receipt.artifact_root().content_hash(),
            receipt.problem_root(),
            receipt.candidate_root(),
            receipt.mesh_root(),
            receipt.operator_root(),
            receipt.coefficient_root(),
            receipt.query_root(),
        ] {
            push_hash(&mut bytes, root);
        }
        push_text(&mut bytes, receipt.qoi());
        push_text(&mut bytes, receipt.units());
        push_text(&mut bytes, &step.hint);
        bytes.push(u8::from(step.discharged));
    }
    fs_ledger::hash_bytes(&bytes)
}

fn ladder_root(budgets: &[f64], rungs: &[usize]) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(
        &mut bytes,
        "fs-flywheel-e2e/ladder-and-escalation-policy/v1",
    );
    bytes.extend_from_slice(
        &u64::try_from(budgets.len())
            .expect("bounded budget ladder fits u64")
            .to_le_bytes(),
    );
    for budget in budgets {
        bytes.extend_from_slice(&budget.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(
        &u64::try_from(rungs.len())
            .expect("bounded rung ladder fits u64")
            .to_le_bytes(),
    );
    for rung in rungs {
        bytes.extend_from_slice(
            &u64::try_from(*rung)
                .expect("bounded rung size fits u64")
                .to_le_bytes(),
        );
    }
    fs_ledger::hash_bytes(&bytes)
}

fn escalation_root(trajectory: &[fs_ir::anytime::IntervalStep]) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(&mut bytes, "fs-flywheel-e2e/escalation-decisions/v1");
    for step in trajectory {
        push_text(&mut bytes, &step.hint);
        bytes.push(u8::from(step.discharged));
        bytes.extend_from_slice(&step.budget.to_bits().to_le_bytes());
    }
    fs_ledger::hash_bytes(&bytes)
}

fn admission_ir_root(report: &admission::AdmissionReport) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(&mut bytes, "fs-flywheel-e2e/admission-ir/v1");
    bytes.extend_from_slice(&report.lowering.ir_version().to_le_bytes());
    push_text(&mut bytes, report.lowering.raw_canonical());
    match report.lowering.lowered_canonical() {
        Some(lowered) => {
            bytes.push(1);
            push_text(&mut bytes, lowered);
        }
        None => bytes.push(0),
    }
    fs_ledger::hash_bytes(&bytes)
}

fn admission_decision_root(
    report: &admission::AdmissionReport,
    cx: &admission::AdmissionContext<'_>,
) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(&mut bytes, "fs-flywheel-e2e/admission-decision/v1");
    push_hash(&mut bytes, admission_ir_root(report));
    push_text(&mut bytes, &report.study);
    bytes.push(u8::from(report.admitted));
    push_text(&mut bytes, &report.diagnosis());
    let capability = cx
        .capability
        .as_ref()
        .expect("acceptance admission has an explicit session capability");
    bytes.extend_from_slice(&capability.cores.to_le_bytes());
    bytes.extend_from_slice(&capability.mem_bytes.to_le_bytes());
    bytes.extend_from_slice(&capability.wall_s.to_bits().to_le_bytes());
    for op in &capability.ops {
        push_text(&mut bytes, op);
    }
    fs_ledger::hash_bytes(&bytes)
}

fn voi_manifest_roots(bound: f64) -> (fs_checker::ContentHash, fs_checker::ContentHash) {
    use fs_plan::voi::{
        Cx, DecisionBudget, LiveDecision, MAX_VOI_EVALUATIONS, MAX_VOI_WORK_UNITS, Probe,
        ProbeKind, UncertaintyNode, rank_purchases,
    };
    let margin = |value: &[f64]| value[0] - 5e-3;
    let decision = LiveDecision {
        margin: &margin,
        arity: 1,
    };
    let nodes = [UncertaintyNode {
        name: "qoi-bound".to_string(),
        lo: 0.0,
        hi: bound.max(1e-6) * 2.0,
        nominal: bound,
    }];
    let menu = [Probe {
        name: "climb-final-rung".to_string(),
        target: "qoi-bound".to_string(),
        cost: 12.0,
        shrink: 0.25,
        kind: ProbeKind::Computational,
    }];
    let ranked = rank_purchases(
        &Cx::for_testing(),
        &decision,
        &nodes,
        &menu,
        32,
        DecisionBudget::new(MAX_VOI_EVALUATIONS, MAX_VOI_WORK_UNITS)
            .expect("valid replay VoI budget"),
        "fs-flywheel-acceptance-v1",
        "ac-004-snapshot-v1",
    )
    .expect("valid replay VoI request");
    let mut inputs = Vec::new();
    push_text(&mut inputs, "fs-flywheel-e2e/voi-inputs/v1");
    inputs.extend_from_slice(&ranked.source_identity_version().to_le_bytes());
    push_text(&mut inputs, "margin/v1:value[0]-threshold");
    inputs.extend_from_slice(&5e-3_f64.to_bits().to_le_bytes());
    push_hash(&mut inputs, ranked.source_context_id());

    let mut outputs = Vec::new();
    push_text(&mut outputs, "fs-flywheel-e2e/voi-outputs/v1");
    outputs.extend_from_slice(&ranked.identity_version().to_le_bytes());
    push_hash(&mut outputs, ranked.context_id());
    (
        fs_ledger::hash_bytes(&inputs),
        fs_ledger::hash_bytes(&outputs),
    )
}

fn whole_path_source_cone_root(
    producer: &VerifierProducerSourceIdentity,
) -> fs_checker::ContentHash {
    let mut bytes = Vec::new();
    push_text(&mut bytes, "fs-flywheel-e2e/whole-path-source-cone/v1");
    push_hash(&mut bytes, producer.producer_source_root());
    push_hash(&mut bytes, producer.dependency_source_root());
    push_hash(&mut bytes, producer.workspace_manifest_root());
    for source in [
        &include_bytes!("acceptance.rs")[..],
        &include_bytes!("../Cargo.toml")[..],
        &include_bytes!("../src/lib.rs")[..],
        &include_bytes!("../../../Cargo.toml")[..],
        &include_bytes!("../../fs-ir/Cargo.toml")[..],
        &include_bytes!("../../fs-ir/src/lib.rs")[..],
        &include_bytes!("../../fs-ir/src/ast.rs")[..],
        &include_bytes!("../../fs-ir/src/sexpr.rs")[..],
        &include_bytes!("../../fs-ir/src/lower.rs")[..],
        &include_bytes!("../../fs-ir/src/study.rs")[..],
        &include_bytes!("../../fs-ir/src/admission.rs")[..],
        &include_bytes!("../../fs-ir/src/planner.rs")[..],
        &include_bytes!("../../fs-ir/src/anytime.rs")[..],
        &include_bytes!("../../fs-plan/Cargo.toml")[..],
        &include_bytes!("../../fs-plan/src/lib.rs")[..],
        &include_bytes!("../../fs-plan/src/voi.rs")[..],
        &include_bytes!("../../fs-package/Cargo.toml")[..],
        &include_bytes!("../../fs-package/src/lib.rs")[..],
        &include_bytes!("../../fs-package/src/origin.rs")[..],
        &include_bytes!("../../fs-checker/Cargo.toml")[..],
        &include_bytes!("../../fs-checker/src/lib.rs")[..],
        &include_bytes!("../../fs-checker/src/semantic.rs")[..],
        &include_bytes!("../../fs-evidence/Cargo.toml")[..],
        &include_bytes!("../../fs-evidence/src/lib.rs")[..],
        &include_bytes!("../../fs-evidence/src/color.rs")[..],
        &include_bytes!("../../fs-evidence/src/admitted.rs")[..],
        &include_bytes!("../../fs-ledger/Cargo.toml")[..],
        &include_bytes!("../../fs-ledger/src/lib.rs")[..],
        &include_bytes!("../../fs-ledger/src/hash.rs")[..],
        &include_bytes!("../../fs-blake3/Cargo.toml")[..],
        &include_bytes!("../../fs-blake3/src/lib.rs")[..],
        &include_bytes!("../../fs-blake3/src/identity.rs")[..],
    ] {
        push_bytes(&mut bytes, source);
    }
    fs_ledger::hash_bytes(&bytes)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct WholePathReplayManifest {
    schema_version: u32,
    admission_ir_root: fs_checker::ContentHash,
    capability_decision_root: fs_checker::ContentHash,
    query_root: fs_checker::ContentHash,
    qoi_root: fs_checker::ContentHash,
    units: String,
    interval_sequence_root: fs_checker::ContentHash,
    verifier_receipts_root: fs_checker::ContentHash,
    ladder_root: fs_checker::ContentHash,
    escalation_root: fs_checker::ContentHash,
    cache_snapshot_root: fs_checker::ContentHash,
    cache_mutations_root: fs_checker::ContentHash,
    voi_inputs_root: fs_checker::ContentHash,
    voi_outputs_root: fs_checker::ContentHash,
    package_root: fs_checker::ContentHash,
    checker_decision_root: fs_checker::ContentHash,
    checker_receipt_root: fs_checker::ContentHash,
    budget_root: fs_checker::ContentHash,
    tolerance_bits: u64,
    planner_budget_bits: u64,
    planner_spent_bits: u64,
    verifier_work: [u128; 6],
    cancellation: String,
    source_cone_root: fs_checker::ContentHash,
    lock_root: fs_checker::ContentHash,
    toolchain_root: fs_checker::ContentHash,
    replay_verifier: String,
}

impl WholePathReplayManifest {
    fn root(&self) -> fs_checker::ContentHash {
        let mut bytes = b"fs-flywheel-e2e/whole-path-manifest/v1".to_vec();
        bytes.extend_from_slice(&self.schema_version.to_le_bytes());
        for root in [
            self.admission_ir_root,
            self.capability_decision_root,
            self.query_root,
            self.qoi_root,
        ] {
            push_hash(&mut bytes, root);
        }
        push_text(&mut bytes, &self.units);
        for root in [
            self.interval_sequence_root,
            self.verifier_receipts_root,
            self.ladder_root,
            self.escalation_root,
            self.cache_snapshot_root,
            self.cache_mutations_root,
            self.voi_inputs_root,
            self.voi_outputs_root,
            self.package_root,
            self.checker_decision_root,
            self.checker_receipt_root,
            self.budget_root,
        ] {
            push_hash(&mut bytes, root);
        }
        bytes.extend_from_slice(&self.tolerance_bits.to_le_bytes());
        bytes.extend_from_slice(&self.planner_budget_bits.to_le_bytes());
        bytes.extend_from_slice(&self.planner_spent_bits.to_le_bytes());
        for work in self.verifier_work {
            bytes.extend_from_slice(&work.to_le_bytes());
        }
        push_text(&mut bytes, &self.cancellation);
        push_hash(&mut bytes, self.source_cone_root);
        push_hash(&mut bytes, self.lock_root);
        push_hash(&mut bytes, self.toolchain_root);
        push_text(&mut bytes, &self.replay_verifier);
        fs_ledger::hash_bytes(&bytes)
    }

    fn verifies(&self, expected: fs_checker::ContentHash) -> bool {
        self.root() == expected
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn ac_004_g5_whole_path_replay() -> Result<(), PlanError> {
    use fs_ir::anytime::run_anytime;
    let run = || -> Result<WholePathReplayManifest, PlanError> {
        let cx = admission_cx();
        let node = sexpr::parse(WEDGE).expect("whole-path replay IR parses");
        let admission = admission::admit(&node, &cx);
        assert!(admission.admitted, "whole-path replay begins admitted");

        let family = steep_family()?;
        let budgets = [30.0, 400.0];
        let mut cache = RecordingCache::default();
        let mut costs = CostTable::new(200.0)?;
        let report = run_anytime(&family, 1.0, 6e-3, &budgets, &RUNGS, &mut cache, &mut costs)?;
        assert!(
            !report.trajectory.is_empty(),
            "replay retains every interval step"
        );

        let mut receipt_sequence = Vec::new();
        push_text(
            &mut receipt_sequence,
            "fs-flywheel-e2e/verifier-receipt-sequence/v1",
        );
        receipt_sequence.extend_from_slice(
            &u64::try_from(report.trajectory.len())
                .expect("bounded receipt sequence fits u64")
                .to_le_bytes(),
        );
        for step in &report.trajectory {
            let resolver = AuthenticReceiptResolver::from_original_step(
                &family,
                1.0,
                6e-3,
                step.budget,
                step.spent,
                &RUNGS,
                step.verifier_receipt.clone(),
            );
            let resolved = resolver
                .resolve()
                .expect("every original trajectory receipt independently replays");
            let admitted_receipt = resolved.receipt();
            assert_eq!(
                admitted_receipt.bound_hi().to_bits(),
                step.bound.to_bits(),
                "receipt sequence binds every interval endpoint"
            );
            assert_eq!(admitted_receipt.flux_hash(), step.flux_hash);
            assert_eq!(
                admitted_receipt.verifier_family(),
                step.verifier_family.as_str()
            );
            assert_eq!(
                resolver.receipt.artifact_root(),
                step.verifier_receipt.artifact_root()
            );
            for root in [
                resolver.receipt.artifact_root().content_hash(),
                resolver.receipt.problem_root(),
                resolver.receipt.candidate_root(),
                resolver.receipt.mesh_root(),
                resolver.receipt.operator_root(),
                resolver.receipt.coefficient_root(),
                resolver.receipt.query_root(),
            ] {
                push_hash(&mut receipt_sequence, root);
            }
        }

        let final_step = report.trajectory.last().expect("final original receipt");
        let resolver = AuthenticReceiptResolver::from_original_step(
            &family,
            1.0,
            6e-3,
            final_step.budget,
            final_step.spent,
            &RUNGS,
            final_step.verifier_receipt.clone(),
        );
        let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(Some(&resolver)) else {
            panic!("the final replay receipt must independently resolve");
        };
        let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
        let pkg = signed_fixture(
            EvidencePackage::new(provenance.clone())
                .with_claim(resolved.claim())
                .with_claim(physical_qoi_no_claim(
                    "whole-path replay has no authentic physical-QoI certifier",
                    report
                        .trajectory
                        .last()
                        .expect("final interval")
                        .bound
                        .max(f64::MIN_POSITIVE),
                )),
            "acceptance-gate/ac-004",
        );
        let package_root = pkg.try_merkle_root().expect("bounded fixture root");
        let source_verifier = ReceiptCertificateVerifier::from_resolver(&resolver, provenance, 0)
            .expect("whole-path receipt resolves before the solver-free package check");
        let signature_verifier = ExactRootSignatureVerifier {
            domain: "acceptance-gate/ac-004",
        };
        let capabilities = fs_checker::VerificationCapabilities::deny_all()
            .with_source_certificates(&source_verifier);
        let check = fs_checker::check_with_capabilities(
            &pkg,
            Some(package_root),
            Some(&signature_verifier),
            &capabilities,
        );
        assert!(
            check.passed(),
            "whole-path checker receipt must be positive"
        );
        assert!(check.validate_decision_hash());
        let checker_receipt = check
            .receipt()
            .expect("positive whole-path check retains its policy receipt");
        assert!(checker_receipt.validate_hash());

        let final_bound = report.trajectory.last().expect("final interval").bound;
        let (voi_inputs_root, voi_outputs_root) = voi_manifest_roots(final_bound);
        let producer = resolver.receipt.producer();
        let mut qoi = Vec::new();
        push_text(&mut qoi, PHYSICAL_QOI_ID);
        push_text(&mut qoi, PHYSICAL_QOI_UNITS);
        push_text(&mut qoi, QOI_ID);
        push_text(&mut qoi, QOI_UNITS);
        push_hash(&mut qoi, resolver.receipt.query_root());
        let mut budget = Vec::new();
        push_text(&mut budget, "fs-flywheel-e2e/whole-path-budgets/v1");
        budget.extend_from_slice(&6e-3_f64.to_bits().to_le_bytes());
        budget.extend_from_slice(&resolver.planner_budget.to_bits().to_le_bytes());
        budget.extend_from_slice(&resolver.planner_spent.to_bits().to_le_bytes());
        for work in resolver.receipt.work_plan() {
            budget.extend_from_slice(&work.to_le_bytes());
        }

        Ok(WholePathReplayManifest {
            schema_version: 1,
            admission_ir_root: admission_ir_root(&admission),
            capability_decision_root: admission_decision_root(&admission, &cx),
            query_root: fs_ledger::hash_bytes(WEDGE.as_bytes()),
            qoi_root: fs_ledger::hash_bytes(&qoi),
            units: format!(
                "physical={PHYSICAL_QOI_UNITS};proxy={QOI_UNITS};physical-promotion=gated"
            ),
            interval_sequence_root: interval_sequence_root(&report.trajectory),
            verifier_receipts_root: fs_ledger::hash_bytes(&receipt_sequence),
            ladder_root: ladder_root(&budgets, &RUNGS),
            escalation_root: escalation_root(&report.trajectory),
            cache_snapshot_root: RecordingCache::initial_snapshot_root(),
            cache_mutations_root: cache.committed_mutations_root(),
            voi_inputs_root,
            voi_outputs_root,
            package_root,
            checker_decision_root: check.decision_hash(),
            checker_receipt_root: checker_receipt.receipt_hash(),
            budget_root: fs_ledger::hash_bytes(&budget),
            tolerance_bits: 6e-3_f64.to_bits(),
            planner_budget_bits: resolver.planner_budget.to_bits(),
            planner_spent_bits: resolver.planner_spent.to_bits(),
            verifier_work: resolver.receipt.work_plan(),
            cancellation: if report.stopped_by_observer() {
                "observer-stopped/partial/no-promotion".to_string()
            } else {
                "completed-drained-finalized-published".to_string()
            },
            source_cone_root: whole_path_source_cone_root(producer),
            lock_root: producer.workspace_lock_root(),
            toolchain_root: producer.toolchain_root(),
            replay_verifier: "fs-flywheel-e2e/whole-path-replay-verifier/v1".to_string(),
        })
    };
    let manifest_a = run()?;
    let manifest_b = run()?;
    assert_eq!(
        manifest_a, manifest_b,
        "every whole-path semantic component replays bit-exact"
    );
    let expected = manifest_a.root();
    assert!(manifest_b.verifies(expected));

    macro_rules! manifest_mutation {
        ($field:literal, $mutation:expr) => {{
            let mut changed = manifest_a.clone();
            $mutation(&mut changed);
            assert!(
                !changed.verifies(expected),
                "mutating or deleting {} must fail the whole-path manifest",
                $field
            );
        }};
    }
    manifest_mutation!("schema version", |m: &mut WholePathReplayManifest| m
        .schema_version +=
        1);
    manifest_mutation!("admission IR", |m: &mut WholePathReplayManifest| m
        .admission_ir_root =
        flip(m.admission_ir_root));
    manifest_mutation!("capability decision", |m: &mut WholePathReplayManifest| m
        .capability_decision_root =
        flip(m.capability_decision_root));
    manifest_mutation!("query", |m: &mut WholePathReplayManifest| m.query_root =
        flip(m.query_root));
    manifest_mutation!("QoI", |m: &mut WholePathReplayManifest| m.qoi_root =
        flip(m.qoi_root));
    manifest_mutation!("units", |m: &mut WholePathReplayManifest| m
        .units
        .push_str("-foreign"));
    manifest_mutation!("interval sequence", |m: &mut WholePathReplayManifest| m
        .interval_sequence_root =
        flip(m.interval_sequence_root));
    manifest_mutation!(
        "authentic receipt sequence",
        |m: &mut WholePathReplayManifest| m.verifier_receipts_root =
            fs_ledger::hash_bytes(b"deleted receipts")
    );
    manifest_mutation!("ladder", |m: &mut WholePathReplayManifest| m.ladder_root =
        flip(m.ladder_root));
    manifest_mutation!("escalation decisions", |m: &mut WholePathReplayManifest| {
        m.escalation_root = flip(m.escalation_root)
    });
    manifest_mutation!("cache snapshot", |m: &mut WholePathReplayManifest| m
        .cache_snapshot_root =
        flip(m.cache_snapshot_root));
    manifest_mutation!("cache mutations", |m: &mut WholePathReplayManifest| m
        .cache_mutations_root =
        fs_ledger::hash_bytes(b"deleted cache mutations"));
    manifest_mutation!("VoI inputs", |m: &mut WholePathReplayManifest| m
        .voi_inputs_root =
        flip(m.voi_inputs_root));
    manifest_mutation!("VoI outputs", |m: &mut WholePathReplayManifest| m
        .voi_outputs_root =
        flip(m.voi_outputs_root));
    manifest_mutation!("package root", |m: &mut WholePathReplayManifest| m
        .package_root =
        flip(m.package_root));
    manifest_mutation!("checker decision", |m: &mut WholePathReplayManifest| m
        .checker_decision_root =
        flip(m.checker_decision_root));
    manifest_mutation!("checker receipt", |m: &mut WholePathReplayManifest| m
        .checker_receipt_root =
        flip(m.checker_receipt_root));
    manifest_mutation!(
        "budgets and consumption",
        |m: &mut WholePathReplayManifest| m.budget_root = flip(m.budget_root)
    );
    manifest_mutation!("tolerance", |m: &mut WholePathReplayManifest| m
        .tolerance_bits ^=
        1);
    manifest_mutation!("planner budget", |m: &mut WholePathReplayManifest| m
        .planner_budget_bits ^=
        1);
    manifest_mutation!("planner consumption", |m: &mut WholePathReplayManifest| {
        m.planner_spent_bits ^= 1
    });
    manifest_mutation!("verifier work", |m: &mut WholePathReplayManifest| m
        .verifier_work[0] +=
        1);
    manifest_mutation!(
        "cancellation disposition",
        |m: &mut WholePathReplayManifest| m.cancellation.push_str("-partial")
    );
    manifest_mutation!("producer source cone", |m: &mut WholePathReplayManifest| {
        m.source_cone_root = flip(m.source_cone_root)
    });
    manifest_mutation!("workspace lock", |m: &mut WholePathReplayManifest| m
        .lock_root =
        flip(m.lock_root));
    manifest_mutation!("toolchain", |m: &mut WholePathReplayManifest| m
        .toolchain_root =
        flip(m.toolchain_root));
    manifest_mutation!("replay verifier", |m: &mut WholePathReplayManifest| m
        .replay_verifier
        .push_str("-substituted"));
    println!(
        "{{\"schema_version\":1,\"suite\":\"fs-flywheel-e2e/acceptance\",\
         \"case\":\"ac-004\",\"seed\":\"0x5EED0008\",\"units\":\"{}\",\
         \"budget_root\":\"{}\",\"tolerance_bits\":\"{:016x}\",\
         \"planner_budget_bits\":\"{:016x}\",\"planner_spent_bits\":\"{:016x}\",\
         \"verifier_work\":{:?},\"cancellation\":\"{}\",\
         \"manifest_root\":\"{expected}\",\"verifier_receipts_root\":\"{}\",\
         \"checker_receipt_root\":\"{}\",\"replay_verifier\":\"{}\",\
         \"physical_promotion\":\"estimated-no-claim\"}}",
        manifest_a.units,
        manifest_a.budget_root,
        manifest_a.tolerance_bits,
        manifest_a.planner_budget_bits,
        manifest_a.planner_spent_bits,
        manifest_a.verifier_work,
        manifest_a.cancellation,
        manifest_a.verifier_receipts_root,
        manifest_a.checker_receipt_root,
        manifest_a.replay_verifier,
    );
    verdict(
        "ac-004",
        "the whole path binds admission, capabilities, query/QoI/units, every interval and \
         authentic verifier receipt, ladder/escalation, cache state/mutations, VoI, package and \
         checker receipts, budgets/cancellation, and source-cone/lock/toolchain inputs; binary \
         attestation remains a no-claim, every bound field replays bit-exact, and every mutation \
         fails closed (G5)",
    );
    Ok(())
}

#[test]
fn ac_005_laundering_invariant_across_the_path() -> Result<(), PlanError> {
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
    let family = steep_family()?;
    let resolver = AuthenticReceiptResolver::capture(&family, 1.0, 6e-3, 400.0, &RUNGS)?;
    let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(Some(&resolver)) else {
        panic!("the authentic hard component must resolve");
    };
    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    let hard_claim = resolved.claim();
    let pkg = signed_fixture(
        EvidencePackage::new(provenance.clone())
            .with_claim(hard_claim)
            .with_claim(Claim::estimated("soft", "estimated part", "dwr-guess", 0.1)),
        "acceptance-gate/ac-005",
    );
    let source_verifier = ReceiptCertificateVerifier::from_resolver(&resolver, provenance, 0)
        .expect("hard-component receipt resolves before package checking");
    let signature_verifier = ExactRootSignatureVerifier {
        domain: "acceptance-gate/ac-005",
    };
    let capabilities = fs_checker::VerificationCapabilities::deny_all()
        .with_source_certificates(&source_verifier)
        .with_signatures(&signature_verifier);
    let breakdown = pkg
        .color_breakdown_with(&capabilities)
        .expect("the authentic verifier receipt authorizes the proxy claim");
    assert!(
        breakdown.verified == 1 && breakdown.estimated == 1,
        "the package cannot blur colors: {breakdown:?}"
    );
    verdict(
        "ac-005",
        "estimated inputs never launder to verified — enforced by the compose algebra \
         and visible in the package breakdown",
    );
    Ok(())
}
