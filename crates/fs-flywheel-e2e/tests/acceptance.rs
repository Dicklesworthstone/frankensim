//! EPISTEMIC-ENGINE ACCEPTANCE (bead xpck.8): the TOP-LEVEL runnable
//! integration gate for the whole addendum — a declarative query produces a
//! colored, priced, auditable package and crosses the solver-free checker's
//! typed capability boundary. Scientific promotion consumes an immutable
//! `fs-verify` receipt and independently replays its exact problem/candidate;
//! scientific promotion additionally consumes injected producer-attestation,
//! executable-observation, and promotion-root capabilities before minting an
//! opaque witness over the stable current-process executable and exact receipt
//! scope.
//! The detached sidecar and package-root signature policies remain
//! deterministic fixtures, not vendor-independent cryptographic
//! authentication. The retained certificate covers the exact 1-D manufactured
//! proxy; the physical wedge QoI stays an explicit Estimated no-claim until its
//! own upstream certifier exists.
//!
//! The path: admission (typed, teaching refusals) → flywheel discharge
//! (planner + cache) → anytime colored answer (+ the VoI-priced hint)
//! → authentic-receipt-backed evidence package → SOLVER-FREE policy re-check →
//! G5 whole-path replay → the laundering invariant at every hop.
#![cfg(feature = "flywheel-e2e")]

use fs_blake3::identity::{
    Admitted, AuthorityAdmitter, AuthorityRef, AuthorityVerifier, CanonicalEncoder,
    CanonicalLimits, CanonicalSchema, ContentId, ExternalAnchorRef, Field, FieldSpec,
    IdentityReceipt, KeyPolicyId, NeverCancel, ObservedIdentity, Presented, PromotionTrustRoot,
    PromotionWitness, SemanticId, StrongIdentity, Verified, VerifierId, WireType,
};
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
const RECEIPT_POLICY: &str = "fs-flywheel-e2e/authentic-receipt-resolver/v2";
const PRODUCER_EXECUTABLE_IDENTITY_VERSION: u32 = 1;
const PRODUCER_EXECUTABLE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-flywheel-e2e.producer-executable.v1";
const PRODUCER_ATTESTATION_CONTEXT: &str = "fs-flywheel-e2e/producer-promotion/v1";
const MAX_PRODUCER_EXECUTABLE_BYTES: u64 = 1 << 30;
const EXECUTABLE_HASH_CHUNK_BYTES: usize = 64 * 1024;
const MAX_PRODUCER_ATTESTATION_SIDECAR_BYTES: usize = 4 * 1024;
const PRODUCER_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(64 * 1024, 16 * 1024, 32, 64, 1024);

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ProducerExecutableIdentity {
    schema_version: u32,
    byte_len: u64,
    raw_hash: fs_checker::ContentHash,
}

impl ProducerExecutableIdentity {
    fn root(self) -> fs_checker::ContentHash {
        let mut bytes = Vec::new();
        push_text(&mut bytes, PRODUCER_EXECUTABLE_IDENTITY_DOMAIN);
        bytes.extend_from_slice(&self.schema_version.to_le_bytes());
        bytes.extend_from_slice(&self.byte_len.to_le_bytes());
        push_hash(&mut bytes, self.raw_hash);
        fs_ledger::hash_bytes(&bytes)
    }
}

fn producer_executable_identity_from_reader(
    reader: &mut impl std::io::Read,
) -> Result<ProducerExecutableIdentity, String> {
    producer_executable_identity_from_reader_with_limit(reader, MAX_PRODUCER_EXECUTABLE_BYTES)
}

fn producer_executable_identity_from_reader_with_limit(
    reader: &mut impl std::io::Read,
    max_bytes: u64,
) -> Result<ProducerExecutableIdentity, String> {
    let mut hasher = fs_ledger::Blake3::new();
    let mut total = 0_u64;
    let mut chunk = Vec::new();
    chunk
        .try_reserve_exact(EXECUTABLE_HASH_CHUNK_BYTES)
        .map_err(|error| format!("cannot reserve executable hash buffer: {error}"))?;
    chunk.resize(EXECUTABLE_HASH_CHUNK_BYTES, 0_u8);
    loop {
        let read = reader
            .read(&mut chunk)
            .map_err(|error| format!("cannot read producer-process executable: {error}"))?;
        if read == 0 {
            break;
        }
        total = total
            .checked_add(u64::try_from(read).expect("read chunk length fits u64"))
            .ok_or_else(|| "producer-process executable length exceeds u64".to_string())?;
        if total > max_bytes {
            return Err(format!(
                "producer-process executable exceeds {max_bytes} bytes"
            ));
        }
        hasher.update(&chunk[..read]);
    }
    Ok(ProducerExecutableIdentity {
        schema_version: PRODUCER_EXECUTABLE_IDENTITY_VERSION,
        byte_len: total,
        raw_hash: hasher.finalize(),
    })
}

fn current_producer_process_executable_identity() -> Result<ProducerExecutableIdentity, String> {
    let path = std::env::current_exe()
        .map_err(|error| format!("cannot resolve producer-process executable: {error}"))?;
    let mut file = std::fs::File::open(&path).map_err(|error| {
        format!(
            "cannot open producer-process executable {}: {error}",
            path.display()
        )
    })?;
    producer_executable_identity_from_reader(&mut file)
}

struct ProducerAttestationSubjectSchemaV1;

impl CanonicalSchema for ProducerAttestationSubjectSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-flywheel-e2e.producer-attestation-subject.v1";
    const NAME: &'static str = "fs-flywheel-producer-attestation-subject";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "Producer executable, receipt scope, and build-input binding";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("purpose", WireType::Utf8),
        FieldSpec::required("evidence_scope_root", WireType::Bytes),
        FieldSpec::required("verifier_receipt_root", WireType::Bytes),
        FieldSpec::required("executable_identity_version", WireType::U64),
        FieldSpec::required("executable_byte_len", WireType::U64),
        FieldSpec::required("executable_raw_hash", WireType::Bytes),
        FieldSpec::required("producer_crate", WireType::Utf8),
        FieldSpec::required("producer_version", WireType::Utf8),
        FieldSpec::required("producer_features", WireType::Utf8),
        FieldSpec::required("producer_source_root", WireType::Bytes),
        FieldSpec::required("dependency_source_root", WireType::Bytes),
        FieldSpec::required("workspace_manifest_root", WireType::Bytes),
        FieldSpec::required("workspace_lock_root", WireType::Bytes),
        FieldSpec::required("toolchain_root", WireType::Bytes),
    ];
}

struct ProducerAttestationVerifierSchemaV1;

impl CanonicalSchema for ProducerAttestationVerifierSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-flywheel-e2e.producer-attestation-verifier.v1";
    const NAME: &'static str = "fs-flywheel-producer-attestation-verifier";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "Pinned detached-sidecar verifier identity";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("implementation", WireType::Utf8),
        FieldSpec::required("sidecar_schema", WireType::U64),
    ];
}

struct ProducerAttestationPolicySchemaV1;

impl CanonicalSchema for ProducerAttestationPolicySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-flywheel-e2e.producer-attestation-policy.v1";
    const NAME: &'static str = "fs-flywheel-producer-attestation-policy";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "Pinned producer-promotion key policy identity";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("policy", WireType::Utf8),
        FieldSpec::required("minimum_sidecar_schema", WireType::U64),
    ];
}

type ProducerAttestationSubject = SemanticId<ProducerAttestationSubjectSchemaV1>;
type PresentedProducerAuthority = AuthorityRef<
    ProducerAttestationSubject,
    ProducerAttestationVerifierSchemaV1,
    ProducerAttestationPolicySchemaV1,
    Presented,
>;
type VerifiedProducerAuthority = AuthorityRef<
    ProducerAttestationSubject,
    ProducerAttestationVerifierSchemaV1,
    ProducerAttestationPolicySchemaV1,
    Verified,
>;
type AdmittedProducerAuthority = AuthorityRef<
    ProducerAttestationSubject,
    ProducerAttestationVerifierSchemaV1,
    ProducerAttestationPolicySchemaV1,
    Admitted,
>;
type ProducerPromotionRoot =
    PromotionTrustRoot<ProducerAttestationVerifierSchemaV1, ProducerAttestationPolicySchemaV1>;
type ProducerPromotionWitness = PromotionWitness<
    ProducerAttestationSubject,
    ProducerAttestationVerifierSchemaV1,
    ProducerAttestationPolicySchemaV1,
>;

#[derive(Debug, Clone, Copy)]
struct ProducerAttestationScope {
    purpose: &'static str,
    evidence_root: fs_checker::ContentHash,
}

fn producer_attestation_verifier_identity(
    implementation: &str,
    sidecar_schema: u64,
) -> IdentityReceipt<VerifierId<ProducerAttestationVerifierSchemaV1>> {
    CanonicalEncoder::<VerifierId<ProducerAttestationVerifierSchemaV1>, _>::new(
        PRODUCER_IDENTITY_LIMITS,
        NeverCancel,
    )
    .expect("valid producer-attestation verifier schema")
    .utf8(Field::new(0, "implementation"), implementation)
    .expect("verifier implementation field")
    .u64(Field::new(1, "sidecar_schema"), sidecar_schema)
    .expect("verifier schema field")
    .finish()
    .expect("producer-attestation verifier identity")
}

fn trusted_producer_attestation_verifier()
-> IdentityReceipt<VerifierId<ProducerAttestationVerifierSchemaV1>> {
    producer_attestation_verifier_identity("fs-flywheel-e2e/exact-detached-sidecar-fixture/v1", 1)
}

fn producer_attestation_policy_identity(
    policy: &str,
    minimum_sidecar_schema: u64,
) -> IdentityReceipt<KeyPolicyId<ProducerAttestationPolicySchemaV1>> {
    CanonicalEncoder::<KeyPolicyId<ProducerAttestationPolicySchemaV1>, _>::new(
        PRODUCER_IDENTITY_LIMITS,
        NeverCancel,
    )
    .expect("valid producer-attestation policy schema")
    .utf8(Field::new(0, "policy"), policy)
    .expect("policy identity field")
    .u64(
        Field::new(1, "minimum_sidecar_schema"),
        minimum_sidecar_schema,
    )
    .expect("policy schema field")
    .finish()
    .expect("producer-attestation policy identity")
}

fn trusted_producer_attestation_policy()
-> IdentityReceipt<KeyPolicyId<ProducerAttestationPolicySchemaV1>> {
    producer_attestation_policy_identity("fs-flywheel-e2e/executable-and-receipt-required/v1", 1)
}

#[derive(Clone, Copy)]
struct ConfiguredProducerPromotion {
    trust_root: ProducerPromotionRoot,
    verifier_observation: fs_blake3::identity::ByteObservation,
    policy_observation: fs_blake3::identity::ByteObservation,
}

fn configured_producer_promotion() -> ConfiguredProducerPromotion {
    let verifier = trusted_producer_attestation_verifier();
    let policy = trusted_producer_attestation_policy();
    ConfiguredProducerPromotion {
        trust_root: ProducerPromotionRoot::configure(
            ObservedIdentity::from_receipt(verifier),
            ObservedIdentity::from_receipt(policy),
            PRODUCER_ATTESTATION_CONTEXT,
        )
        .expect("non-empty pinned producer-promotion context"),
        verifier_observation: ObservedIdentity::from_receipt(verifier).bytes(),
        policy_observation: ObservedIdentity::from_receipt(policy).bytes(),
    }
}

trait ProducerAttestationAuthority {
    fn verify_and_admit(
        &self,
        expected_subject: IdentityReceipt<ProducerAttestationSubject>,
    ) -> Result<AdmittedProducerAuthority, &'static str>;
}

#[derive(Clone, Copy)]
struct ProducerAttestationCapability<'a> {
    authority: &'a dyn ProducerAttestationAuthority,
    observed_executable: ProducerExecutableIdentity,
}

#[derive(Debug, Clone)]
struct DetachedFixtureProducerAuthority {
    subject: IdentityReceipt<ProducerAttestationSubject>,
    anchor_bytes: Vec<u8>,
    verifier: IdentityReceipt<VerifierId<ProducerAttestationVerifierSchemaV1>>,
    policy: IdentityReceipt<KeyPolicyId<ProducerAttestationPolicySchemaV1>>,
}

fn fixture_attestation_sidecar_bytes(
    subject: IdentityReceipt<ProducerAttestationSubject>,
    verifier: IdentityReceipt<VerifierId<ProducerAttestationVerifierSchemaV1>>,
    policy: IdentityReceipt<KeyPolicyId<ProducerAttestationPolicySchemaV1>>,
) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_text(
        &mut bytes,
        "fs-flywheel-e2e/detached-producer-attestation-fixture/v1",
    );
    push_bytes(&mut bytes, subject.id().as_bytes());
    push_bytes(&mut bytes, subject.canonical_preimage().as_bytes());
    bytes.extend_from_slice(&subject.canonical_bytes().to_le_bytes());
    push_bytes(&mut bytes, verifier.id().as_bytes());
    push_bytes(&mut bytes, policy.id().as_bytes());
    push_text(
        &mut bytes,
        "DETERMINISTIC-FIXTURE-NOT-A-PRODUCTION-SIGNATURE",
    );
    assert!(
        bytes.len() <= MAX_PRODUCER_ATTESTATION_SIDECAR_BYTES,
        "bounded producer-attestation sidecar fixture"
    );
    bytes
}

impl DetachedFixtureProducerAuthority {
    fn for_subject(subject: IdentityReceipt<ProducerAttestationSubject>) -> Self {
        let verifier = trusted_producer_attestation_verifier();
        let policy = trusted_producer_attestation_policy();
        let anchor_bytes = fixture_attestation_sidecar_bytes(subject, verifier, policy);
        Self {
            subject,
            anchor_bytes,
            verifier,
            policy,
        }
    }

    fn presented_anchor(&self) -> ExternalAnchorRef {
        ExternalAnchorRef::presented(ContentId::of_bytes(&self.anchor_bytes))
    }

    fn exact_sidecar_matches(&self) -> bool {
        self.anchor_bytes
            == fixture_attestation_sidecar_bytes(self.subject, self.verifier, self.policy)
    }
}

impl ProducerAttestationAuthority for DetachedFixtureProducerAuthority {
    fn verify_and_admit(
        &self,
        expected_subject: IdentityReceipt<ProducerAttestationSubject>,
    ) -> Result<AdmittedProducerAuthority, &'static str> {
        if expected_subject != self.subject {
            return Err("detached producer attestation belongs to another subject");
        }
        AuthorityRef::present(
            self.subject,
            self.presented_anchor(),
            self.verifier.id(),
            self.policy.id(),
        )
        .verify(self)?
        .admit(self)
    }
}

impl
    AuthorityVerifier<
        ProducerAttestationSubject,
        ProducerAttestationVerifierSchemaV1,
        ProducerAttestationPolicySchemaV1,
    > for DetachedFixtureProducerAuthority
{
    type Error = &'static str;

    fn verify(&self, presented: &PresentedProducerAuthority) -> Result<(), Self::Error> {
        if presented.receipt() == self.subject
            && presented.anchor() == self.presented_anchor()
            && presented.verifier() == self.verifier.id()
            && presented.key_policy() == self.policy.id()
            && self.exact_sidecar_matches()
        {
            Ok(())
        } else {
            Err("detached producer attestation failed exact verification")
        }
    }
}

impl
    AuthorityAdmitter<
        ProducerAttestationSubject,
        ProducerAttestationVerifierSchemaV1,
        ProducerAttestationPolicySchemaV1,
    > for DetachedFixtureProducerAuthority
{
    type Error = &'static str;

    fn admit(&self, verified: &VerifiedProducerAuthority) -> Result<(), Self::Error> {
        if verified.receipt() == self.subject
            && verified.anchor() == self.presented_anchor()
            && verified.verifier() == self.verifier.id()
            && verified.key_policy() == self.policy.id()
            && self.exact_sidecar_matches()
        {
            Ok(())
        } else {
            Err("producer attestation failed exact fixture policy admission")
        }
    }
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
    producer_executable: ProducerExecutableIdentity,
}

struct ReplayedVerifierReceipt<'a> {
    admitted: AdmittedVerifierReceipt<'a>,
}

struct PromotedVerifierReceipt<'a> {
    scientific: ReplayedVerifierReceipt<'a>,
    producer_executable: ProducerExecutableIdentity,
    promotion: ProducerPromotionWitness,
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

impl ReplayedVerifierReceipt<'_> {
    fn receipt(&self) -> &AdmittedVerifierReceipt<'_> {
        &self.admitted
    }
}

impl PromotedVerifierReceipt<'_> {
    fn receipt(&self) -> &AdmittedVerifierReceipt<'_> {
        self.scientific.receipt()
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

    fn producer_executable_root(&self) -> fs_checker::ContentHash {
        self.producer_executable.root()
    }

    fn producer_attestation_anchor_root(&self) -> fs_checker::ContentHash {
        fs_checker::ContentHash(*self.promotion.anchor().content_id().as_bytes())
    }

    fn producer_attestation_policy_root(&self) -> fs_checker::ContentHash {
        fs_checker::ContentHash(*self.promotion.key_policy().id().as_bytes())
    }

    fn producer_attestation_root(&self) -> fs_checker::ContentHash {
        let audit = self.promotion.audit();
        let mut bytes = Vec::new();
        push_text(
            &mut bytes,
            "fs-flywheel-e2e/admitted-producer-attestation/v1",
        );
        push_bytes(&mut bytes, self.promotion.subject().id().as_bytes());
        push_bytes(
            &mut bytes,
            self.promotion.subject().canonical_preimage().as_bytes(),
        );
        bytes.extend_from_slice(&self.promotion.subject().canonical_bytes().to_le_bytes());
        push_bytes(&mut bytes, self.promotion.anchor().content_id().as_bytes());
        push_text(&mut bytes, audit.verifier_domain);
        push_bytes(
            &mut bytes,
            audit.verifier_observation.content_id().as_bytes(),
        );
        bytes.extend_from_slice(&audit.verifier_observation.length().to_le_bytes());
        push_text(&mut bytes, audit.key_policy_domain);
        push_bytes(
            &mut bytes,
            audit.key_policy_observation.content_id().as_bytes(),
        );
        bytes.extend_from_slice(&audit.key_policy_observation.length().to_le_bytes());
        push_text(&mut bytes, audit.context);
        fs_ledger::hash_bytes(&bytes)
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
        let producer_executable = current_producer_process_executable_identity()
            .expect("producer executable must be readable before verifier work");
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
            producer_executable,
        })
    }

    fn from_original_step(
        family: &ProblemFamily,
        theta: f64,
        tolerance: f64,
        planner_budget: f64,
        planner_spent: f64,
        rungs: &[usize],
        producer_executable: ProducerExecutableIdentity,
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
            producer_executable,
        }
    }

    fn resolve(&self) -> Result<ReplayedVerifierReceipt<'_>, &'static str> {
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
        Ok(ReplayedVerifierReceipt { admitted })
    }

    fn producer_attestation_subject(
        &self,
        scope: ProducerAttestationScope,
        executable: ProducerExecutableIdentity,
    ) -> Result<IdentityReceipt<ProducerAttestationSubject>, &'static str> {
        let producer = self.receipt.producer();
        CanonicalEncoder::<ProducerAttestationSubject, _>::new(
            PRODUCER_IDENTITY_LIMITS,
            NeverCancel,
        )
        .map_err(|_| "producer-attestation subject schema refused")?
        .utf8(Field::new(0, "purpose"), scope.purpose)
        .map_err(|_| "producer-attestation purpose refused")?
        .bytes(
            Field::new(1, "evidence_scope_root"),
            scope.evidence_root.as_bytes(),
        )
        .map_err(|_| "producer-attestation evidence scope refused")?
        .bytes(
            Field::new(2, "verifier_receipt_root"),
            self.receipt.artifact_root().content_hash().as_bytes(),
        )
        .map_err(|_| "producer-attestation verifier receipt refused")?
        .u64(
            Field::new(3, "executable_identity_version"),
            u64::from(executable.schema_version),
        )
        .map_err(|_| "producer-attestation executable schema refused")?
        .u64(Field::new(4, "executable_byte_len"), executable.byte_len)
        .map_err(|_| "producer-attestation executable length refused")?
        .bytes(
            Field::new(5, "executable_raw_hash"),
            executable.raw_hash.as_bytes(),
        )
        .map_err(|_| "producer-attestation executable hash refused")?
        .utf8(Field::new(6, "producer_crate"), producer.crate_name())
        .map_err(|_| "producer-attestation crate refused")?
        .utf8(Field::new(7, "producer_version"), producer.crate_version())
        .map_err(|_| "producer-attestation version refused")?
        .utf8(Field::new(8, "producer_features"), producer.features())
        .map_err(|_| "producer-attestation features refused")?
        .bytes(
            Field::new(9, "producer_source_root"),
            producer.producer_source_root().as_bytes(),
        )
        .map_err(|_| "producer-attestation source root refused")?
        .bytes(
            Field::new(10, "dependency_source_root"),
            producer.dependency_source_root().as_bytes(),
        )
        .map_err(|_| "producer-attestation dependency root refused")?
        .bytes(
            Field::new(11, "workspace_manifest_root"),
            producer.workspace_manifest_root().as_bytes(),
        )
        .map_err(|_| "producer-attestation workspace manifest refused")?
        .bytes(
            Field::new(12, "workspace_lock_root"),
            producer.workspace_lock_root().as_bytes(),
        )
        .map_err(|_| "producer-attestation workspace lock refused")?
        .bytes(
            Field::new(13, "toolchain_root"),
            producer.toolchain_root().as_bytes(),
        )
        .map_err(|_| "producer-attestation toolchain refused")?
        .finish()
        .map_err(|_| "producer-attestation subject publication refused")
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
    Verified(PromotedVerifierReceipt<'a>),
    Gated {
        color: Color,
        no_claim: &'static str,
    },
}

fn producer_attestation_fixture(
    resolver: &AuthenticReceiptResolver,
    scope: ProducerAttestationScope,
) -> DetachedFixtureProducerAuthority {
    let subject = resolver
        .producer_attestation_subject(scope, resolver.producer_executable)
        .expect("bounded fixture producer-attestation subject");
    DetachedFixtureProducerAuthority::for_subject(subject)
}

fn fixture_producer_attestation_capability(
    authority: &DetachedFixtureProducerAuthority,
) -> ProducerAttestationCapability<'_> {
    ProducerAttestationCapability {
        authority,
        observed_executable: current_producer_process_executable_identity()
            .expect("fixture producer executable must remain readable before promotion"),
    }
}

fn producer_attestation_gated<'a>(
    estimator: &'static str,
    no_claim: &'static str,
) -> ReceiptPromotion<'a> {
    ReceiptPromotion::Gated {
        color: Color::Estimated {
            estimator: estimator.to_string(),
            dispersion: 1.0,
        },
        no_claim,
    }
}

fn resolve_for_promotion<'a>(
    resolver: Option<&'a AuthenticReceiptResolver>,
    scope: ProducerAttestationScope,
    attestation: Option<ProducerAttestationCapability<'_>>,
    configured_promotion: Option<ConfiguredProducerPromotion>,
) -> ReceiptPromotion<'a> {
    let Some(resolver) = resolver else {
        return ReceiptPromotion::Gated {
            color: Color::Estimated {
                estimator: "missing-authentic-verifier-receipt".to_string(),
                dispersion: 1.0,
            },
            no_claim: "NO-CLAIM: upstream verifier receipt authority is unavailable",
        };
    };
    let Some(attestation) = attestation else {
        return producer_attestation_gated(
            "missing-producer-attestation-capability",
            "NO-CLAIM: producer-attestation capability is unavailable",
        );
    };
    let Some(configured_promotion) = configured_promotion else {
        return producer_attestation_gated(
            "missing-producer-promotion-trust-root",
            "NO-CLAIM: independently configured producer-promotion trust root is unavailable",
        );
    };
    let observed_executable = attestation.observed_executable;
    if observed_executable != resolver.producer_executable {
        return producer_attestation_gated(
            "producer-process-executable-drift",
            "NO-CLAIM: retained producer-process executable observation differs from verifier publication",
        );
    }
    match resolver.resolve() {
        Ok(scientific) if scientific.receipt().accepted() => {
            let Ok(expected_subject) =
                resolver.producer_attestation_subject(scope, observed_executable)
            else {
                return producer_attestation_gated(
                    "producer-attestation-subject-refused",
                    "NO-CLAIM: producer attestation subject could not be canonicalized",
                );
            };
            let Ok(policy_relative) = attestation.authority.verify_and_admit(expected_subject)
            else {
                return producer_attestation_gated(
                    "producer-attestation-capability-refusal",
                    "NO-CLAIM: detached producer-attestation capability was refused",
                );
            };
            let Ok(promotion) = configured_promotion.trust_root.admit_for_promotion(
                &policy_relative,
                configured_promotion.verifier_observation,
                configured_promotion.policy_observation,
            ) else {
                return producer_attestation_gated(
                    "producer-attestation-trust-root-refusal",
                    "NO-CLAIM: producer attestation did not match the pinned promotion trust root",
                );
            };
            if promotion.subject() != expected_subject
                || promotion.context() != PRODUCER_ATTESTATION_CONTEXT
            {
                return producer_attestation_gated(
                    "producer-attestation-witness-mismatch",
                    "NO-CLAIM: producer promotion witness did not bind the expected subject",
                );
            }
            ReceiptPromotion::Verified(PromotedVerifierReceipt {
                scientific,
                producer_executable: observed_executable,
                promotion,
            })
        }
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

fn assert_producer_attestation_gated(
    promotion: ReceiptPromotion<'_>,
    expected_estimator: &str,
    scenario: &str,
) {
    match promotion {
        ReceiptPromotion::Gated { color, no_claim } => {
            let Color::Estimated { estimator, .. } = color else {
                panic!("{scenario} must remain Estimated");
            };
            assert_eq!(
                estimator, expected_estimator,
                "{scenario} must fail at the intended attestation boundary"
            );
            assert!(
                no_claim.contains("NO-CLAIM"),
                "{scenario} must retain an explicit no-claim"
            );
        }
        ReceiptPromotion::Verified(_) => panic!("{scenario} must not mint promotion authority"),
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

struct ReceiptCertificateVerifier<'receipt, 'promotion> {
    promoted: &'promotion PromotedVerifierReceipt<'receipt>,
    provenance: Provenance,
    claim_index: usize,
    expected_claim_subject_hash: fs_checker::ContentHash,
}

impl<'receipt, 'promotion> ReceiptCertificateVerifier<'receipt, 'promotion> {
    fn from_promoted(
        promoted: &'promotion PromotedVerifierReceipt<'receipt>,
        provenance: Provenance,
        claim_index: usize,
    ) -> Self {
        let expected_claim_subject_hash = promoted
            .claim()
            .declared_source_certificate_subject_hash_unverified();
        Self {
            promoted,
            provenance,
            claim_index,
            expected_claim_subject_hash,
        }
    }

    fn policy_fingerprint(&self) -> fs_checker::ContentHash {
        let mut bytes = Vec::new();
        push_text(&mut bytes, RECEIPT_POLICY);
        push_hash(
            &mut bytes,
            self.promoted.receipt().artifact_root().content_hash(),
        );
        push_hash(&mut bytes, self.promoted.producer_executable_root());
        push_hash(&mut bytes, self.promoted.producer_attestation_root());
        push_hash(&mut bytes, self.promoted.producer_attestation_policy_root());
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

impl fs_checker::SourceCertificateVerifier for ReceiptCertificateVerifier<'_, '_> {
    fn verify(
        &self,
        request: &fs_checker::SourceCertificateRequest<'_>,
    ) -> fs_checker::VerificationDecision {
        let fingerprint = self.policy_fingerprint();
        let receipt = self.promoted.receipt();
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
        cost_freshness: None,
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
    let producer_executable = current_producer_process_executable_identity()
        .expect("capture producer executable before ac-003 verifier work");
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
    assert!(
        report.refusal.is_none(),
        "the final package rung must discharge without a teaching refusal"
    );
    let last = report.trajectory.last().expect("answer");
    assert!(
        last.discharged,
        "the retained package receipt must be the discharged final answer"
    );
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
        producer_executable,
        last.verifier_receipt.clone(),
    );
    let attestation_scope = ProducerAttestationScope {
        purpose: "fs-flywheel-e2e/ac-003-final-receipt/v1",
        evidence_root: resolver.receipt.artifact_root().content_hash(),
    };
    let producer_authority = producer_attestation_fixture(&resolver, attestation_scope);
    let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(
        Some(&resolver),
        attestation_scope,
        Some(fixture_producer_attestation_capability(&producer_authority)),
        Some(configured_producer_promotion()),
    ) else {
        panic!("the authentic final receipt and producer attestation must resolve");
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
    // SOLVER-FREE PACKAGE RE-CHECK: the checker receives only the promoted
    // receipt carrying both independent scientific replay and the opaque
    // producer witness. Its injected capability can match the exact request
    // without importing or invoking a solver.
    let source_verifier = ReceiptCertificateVerifier::from_promoted(&resolved, provenance, 0);
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
    match resolve_for_promotion(None, attestation_scope, None, None) {
        ReceiptPromotion::Gated { color, no_claim } => {
            assert!(matches!(color, Color::Estimated { .. }));
            assert!(no_claim.contains("NO-CLAIM") && no_claim.contains("unavailable"));
        }
        ReceiptPromotion::Verified(_) => {
            panic!("missing upstream authority must never mint a verified claim")
        }
    }
    match resolve_for_promotion(
        Some(&resolver),
        attestation_scope,
        None,
        Some(configured_producer_promotion()),
    ) {
        ReceiptPromotion::Gated { color, no_claim } => {
            assert!(matches!(color, Color::Estimated { .. }));
            assert!(no_claim.contains("NO-CLAIM") && no_claim.contains("attestation"));
        }
        ReceiptPromotion::Verified(_) => {
            panic!("a scientific receipt without producer attestation must stay gated")
        }
    }
    let pie = check.render_pie();
    println!(
        "{{\"schema_version\":2,\"suite\":\"fs-flywheel-e2e/acceptance\",\
         \"case\":\"ac-003\",\"seed\":\"0x5EED0008\",\"units\":\"{QOI_UNITS}\",\
         \"budget_cells_bits\":\"{:016x}\",\"spent_cells_bits\":\"{:016x}\",\
         \"cancellation\":\"completed-drained-finalized-published\",\
         \"artifact_root\":\"{root}\",\"verifier_family\":\"{}\",\
         \"verifier_receipt\":\"{}\",\"producer_executable_root\":\"{}\",\
         \"producer_attestation_root\":\"{}\",\"physical_promotion\":\"estimated-no-claim\",\
         \"hint\":{},\"pie\":\"{}\"}}",
        resolver.planner_budget.to_bits(),
        resolver.planner_spent.to_bits(),
        resolver.receipt.verifier_family(),
        resolver.receipt.artifact_root(),
        resolved.producer_executable_root(),
        resolved.producer_attestation_root(),
        hint.to_json(),
        pie.replace('"', "'").replace('\n', " | ")
    );
    verdict(
        "ac-003",
        "an immutable fs-verify receipt is independently replayed and bound to the stable \
         producer-process executable through a pinned, fixture-backed opaque promotion witness \
         before its exact 1-D proxy result enters the package; the physical wedge QoI and missing \
         authority remain explicit Estimated no-claims",
    );
    Ok(())
}

fn package_accepts_claim(
    promoted: &PromotedVerifierReceipt<'_>,
    provenance: &Provenance,
    claim: Claim,
) -> bool {
    let package = signed_fixture(
        EvidencePackage::new(provenance.clone()).with_claim(claim),
        "acceptance-gate/ac-003-mutation",
    );
    let source_verifier =
        ReceiptCertificateVerifier::from_promoted(promoted, provenance.clone(), 0);
    let signature_verifier = ExactRootSignatureVerifier {
        domain: "acceptance-gate/ac-003-mutation",
    };
    let capabilities =
        fs_checker::VerificationCapabilities::deny_all().with_source_certificates(&source_verifier);
    fs_checker::check_with_capabilities(&package, None, Some(&signature_verifier), &capabilities)
        .passed()
}

#[test]
fn ac_003_producer_executable_identity_is_content_exact_and_bounded() {
    struct FragmentedReader<'a> {
        bytes: &'a [u8],
        offset: usize,
        max_chunk: usize,
    }

    impl std::io::Read for FragmentedReader<'_> {
        fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
            let remaining = &self.bytes[self.offset..];
            let take = remaining.len().min(output.len()).min(self.max_chunk);
            output[..take].copy_from_slice(&remaining[..take]);
            self.offset += take;
            Ok(take)
        }
    }

    let payload = b"producer-executable-content-fixture".repeat(4_096);
    let mut contiguous = std::io::Cursor::new(payload.as_slice());
    let exact = producer_executable_identity_from_reader(&mut contiguous)
        .expect("bounded contiguous executable fixture");
    assert_eq!(exact.schema_version, PRODUCER_EXECUTABLE_IDENTITY_VERSION);
    assert_eq!(
        exact.byte_len,
        u64::try_from(payload.len()).expect("fixture length fits u64")
    );
    assert_eq!(exact.raw_hash, fs_ledger::hash_bytes(&payload));
    let mut fragmented = FragmentedReader {
        bytes: &payload,
        offset: 0,
        max_chunk: 37,
    };
    let rechunked = producer_executable_identity_from_reader(&mut fragmented)
        .expect("bounded fragmented executable fixture");
    assert_eq!(
        exact, rechunked,
        "read partitioning and path-free reader choice are nonsemantic"
    );

    let mut over_cap = std::io::Cursor::new(payload.as_slice());
    assert!(
        producer_executable_identity_from_reader_with_limit(
            &mut over_cap,
            u64::try_from(payload.len() - 1).expect("fixture length fits u64"),
        )
        .is_err(),
        "the executable byte cap refuses before identity publication"
    );

    let mut stale_schema = exact;
    stale_schema.schema_version += 1;
    assert_ne!(stale_schema.root(), exact.root());
    let mut changed_length = exact;
    changed_length.byte_len += 1;
    assert_ne!(changed_length.root(), exact.root());
    let mut changed_content = exact;
    changed_content.raw_hash = flip(changed_content.raw_hash);
    assert_ne!(changed_content.root(), exact.root());
}

#[test]
#[allow(clippy::too_many_lines)] // one hostile receipt/attestation substitution battery
fn ac_003_production_receipt_and_claim_subject_substitutions_fail_closed() -> Result<(), PlanError>
{
    let family = steep_family()?;
    let resolver = AuthenticReceiptResolver::capture(&family, 1.0, 6e-3, 400.0, &RUNGS)?;
    let attestation_scope = ProducerAttestationScope {
        purpose: "fs-flywheel-e2e/ac-003-mutation-receipt/v1",
        evidence_root: resolver.receipt.artifact_root().content_hash(),
    };
    let producer_authority = producer_attestation_fixture(&resolver, attestation_scope);
    let producer_capability = fixture_producer_attestation_capability(&producer_authority);
    let configured_promotion = configured_producer_promotion();
    let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(
        Some(&resolver),
        attestation_scope,
        Some(producer_capability),
        Some(configured_promotion),
    ) else {
        panic!("baseline production receipt must promote");
    };
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            None,
            Some(configured_promotion),
        ),
        "missing-producer-attestation-capability",
        "missing detached attestation authority",
    );
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(producer_capability),
            None,
        ),
        "missing-producer-promotion-trust-root",
        "an attestation authority cannot synthesize its own promotion trust root",
    );
    let mut changed_verifier_observation = configured_promotion;
    changed_verifier_observation.verifier_observation = fs_blake3::identity::ByteObservation::new(
        configured_promotion.verifier_observation.content_id(),
        configured_promotion.verifier_observation.length() + 1,
    );
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(producer_capability),
            Some(changed_verifier_observation),
        ),
        "producer-attestation-trust-root-refusal",
        "changed retained verifier bytes",
    );
    let mut changed_policy_observation = configured_promotion;
    changed_policy_observation.policy_observation = fs_blake3::identity::ByteObservation::new(
        configured_promotion.policy_observation.content_id(),
        configured_promotion.policy_observation.length() + 1,
    );
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(producer_capability),
            Some(changed_policy_observation),
        ),
        "producer-attestation-trust-root-refusal",
        "changed retained key-policy bytes",
    );

    let mut stale_hash = resolver.clone();
    stale_hash.producer_executable.raw_hash = flip(stale_hash.producer_executable.raw_hash);
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&stale_hash),
            attestation_scope,
            Some(producer_capability),
            Some(configured_promotion),
        ),
        "producer-process-executable-drift",
        "foreign executable content hash",
    );

    let mut foreign_executable = resolver.producer_executable;
    foreign_executable.raw_hash = flip(foreign_executable.raw_hash);
    let foreign_executable_subject = resolver
        .producer_attestation_subject(attestation_scope, foreign_executable)
        .expect("bounded foreign-executable attestation subject");
    let foreign_executable_authority =
        DetachedFixtureProducerAuthority::for_subject(foreign_executable_subject);
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(fixture_producer_attestation_capability(
                &foreign_executable_authority,
            )),
            Some(configured_promotion),
        ),
        "producer-attestation-capability-refusal",
        "valid sidecar transplanted from a foreign executable",
    );

    let mut truncated_sidecar = producer_authority.clone();
    truncated_sidecar
        .anchor_bytes
        .pop()
        .expect("fixture sidecar is non-empty");
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(fixture_producer_attestation_capability(&truncated_sidecar)),
            Some(configured_promotion),
        ),
        "producer-attestation-capability-refusal",
        "truncated detached attestation sidecar",
    );

    let foreign_scope = ProducerAttestationScope {
        purpose: attestation_scope.purpose,
        evidence_root: flip(attestation_scope.evidence_root),
    };
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            foreign_scope,
            Some(producer_capability),
            Some(configured_promotion),
        ),
        "producer-attestation-capability-refusal",
        "attestation transplanted to another receipt scope",
    );

    let mut foreign_verifier = producer_authority.clone();
    foreign_verifier.verifier =
        producer_attestation_verifier_identity("foreign/permit-all-verifier/v1", 1);
    foreign_verifier.anchor_bytes = fixture_attestation_sidecar_bytes(
        foreign_verifier.subject,
        foreign_verifier.verifier,
        foreign_verifier.policy,
    );
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(fixture_producer_attestation_capability(&foreign_verifier)),
            Some(configured_promotion),
        ),
        "producer-attestation-trust-root-refusal",
        "foreign verifier identity",
    );

    let mut foreign_policy = producer_authority.clone();
    foreign_policy.policy = producer_attestation_policy_identity("foreign/policy/v1", 1);
    foreign_policy.anchor_bytes = fixture_attestation_sidecar_bytes(
        foreign_policy.subject,
        foreign_policy.verifier,
        foreign_policy.policy,
    );
    assert_producer_attestation_gated(
        resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(fixture_producer_attestation_capability(&foreign_policy)),
            Some(configured_promotion),
        ),
        "producer-attestation-trust-root-refusal",
        "foreign producer-promotion key policy",
    );

    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    assert!(package_accepts_claim(
        &resolved,
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
        !package_accepts_claim(&resolved, &provenance, fake_hash),
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
        !package_accepts_claim(&resolved, &provenance, relabeled),
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
        !package_accepts_claim(&resolved, &provenance, altered_endpoint),
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
        !package_accepts_claim(&resolved, &provenance, cross_qoi),
        "the same receipt address cannot authenticate a different QoI subject"
    );

    let mut cross_problem = resolver.clone();
    cross_problem.theta = 0.75;
    assert!(
        cross_problem.resolve().is_err(),
        "an authentic receipt cannot replay against a substituted problem parameter"
    );
    let mut cross_tolerance = resolver.clone();
    cross_tolerance.tolerance = 7e-3;
    assert!(
        cross_tolerance.resolve().is_err(),
        "an authentic receipt cannot replay against a substituted tolerance"
    );

    // Exercise the claim-subject comparison directly: every visible request
    // field remains exact and only the address-free subject digest changes.
    let verifier = ReceiptCertificateVerifier::from_promoted(&resolved, provenance.clone(), 0);
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
    producer_executable_root: fs_checker::ContentHash,
    producer_attestation_root: fs_checker::ContentHash,
    producer_attestation_anchor_root: fs_checker::ContentHash,
    producer_attestation_policy_root: fs_checker::ContentHash,
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
        let mut bytes = b"fs-flywheel-e2e/whole-path-manifest/v2".to_vec();
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
            self.producer_executable_root,
            self.producer_attestation_root,
            self.producer_attestation_anchor_root,
            self.producer_attestation_policy_root,
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
        let producer_executable = current_producer_process_executable_identity()
            .expect("capture producer executable before ac-004 whole-path work");
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
                producer_executable,
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
        let verifier_receipts_root = fs_ledger::hash_bytes(&receipt_sequence);

        let final_step = report.trajectory.last().expect("final original receipt");
        let resolver = AuthenticReceiptResolver::from_original_step(
            &family,
            1.0,
            6e-3,
            final_step.budget,
            final_step.spent,
            &RUNGS,
            producer_executable,
            final_step.verifier_receipt.clone(),
        );
        let attestation_scope = ProducerAttestationScope {
            purpose: "fs-flywheel-e2e/ac-004-receipt-sequence/v1",
            evidence_root: verifier_receipts_root,
        };
        let producer_authority = producer_attestation_fixture(&resolver, attestation_scope);
        let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(
            Some(&resolver),
            attestation_scope,
            Some(fixture_producer_attestation_capability(&producer_authority)),
            Some(configured_producer_promotion()),
        ) else {
            panic!("the final replay receipt sequence and producer attestation must resolve");
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
        let source_verifier = ReceiptCertificateVerifier::from_promoted(&resolved, provenance, 0);
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
            schema_version: 2,
            admission_ir_root: admission_ir_root(&admission),
            capability_decision_root: admission_decision_root(&admission, &cx),
            query_root: fs_ledger::hash_bytes(WEDGE.as_bytes()),
            qoi_root: fs_ledger::hash_bytes(&qoi),
            units: format!(
                "physical={PHYSICAL_QOI_UNITS};proxy={QOI_UNITS};physical-promotion=gated"
            ),
            interval_sequence_root: interval_sequence_root(&report.trajectory),
            verifier_receipts_root,
            ladder_root: ladder_root(&budgets, &RUNGS),
            escalation_root: escalation_root(&report.trajectory),
            cache_snapshot_root: RecordingCache::initial_snapshot_root(),
            cache_mutations_root: cache.committed_mutations_root(),
            voi_inputs_root,
            voi_outputs_root,
            package_root,
            checker_decision_root: check.decision_hash(),
            checker_receipt_root: checker_receipt.receipt_hash(),
            producer_executable_root: resolved.producer_executable_root(),
            producer_attestation_root: resolved.producer_attestation_root(),
            producer_attestation_anchor_root: resolved.producer_attestation_anchor_root(),
            producer_attestation_policy_root: resolved.producer_attestation_policy_root(),
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
            replay_verifier: "fs-flywheel-e2e/whole-path-replay-verifier/v2".to_string(),
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
    manifest_mutation!("producer executable", |m: &mut WholePathReplayManifest| m
        .producer_executable_root =
        flip(m.producer_executable_root));
    manifest_mutation!("producer attestation", |m: &mut WholePathReplayManifest| {
        m.producer_attestation_root = flip(m.producer_attestation_root)
    });
    manifest_mutation!(
        "producer attestation anchor",
        |m: &mut WholePathReplayManifest| m.producer_attestation_anchor_root =
            flip(m.producer_attestation_anchor_root)
    );
    manifest_mutation!(
        "producer attestation policy",
        |m: &mut WholePathReplayManifest| m.producer_attestation_policy_root =
            flip(m.producer_attestation_policy_root)
    );
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
        "{{\"schema_version\":2,\"suite\":\"fs-flywheel-e2e/acceptance\",\
         \"case\":\"ac-004\",\"seed\":\"0x5EED0008\",\"units\":\"{}\",\
         \"budget_root\":\"{}\",\"tolerance_bits\":\"{:016x}\",\
         \"planner_budget_bits\":\"{:016x}\",\"planner_spent_bits\":\"{:016x}\",\
         \"verifier_work\":{:?},\"cancellation\":\"{}\",\
         \"manifest_root\":\"{expected}\",\"verifier_receipts_root\":\"{}\",\
         \"checker_receipt_root\":\"{}\",\"producer_executable_root\":\"{}\",\
         \"producer_attestation_root\":\"{}\",\"producer_attestation_anchor_root\":\"{}\",\
         \"producer_attestation_policy_root\":\"{}\",\"replay_verifier\":\"{}\",\
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
        manifest_a.producer_executable_root,
        manifest_a.producer_attestation_root,
        manifest_a.producer_attestation_anchor_root,
        manifest_a.producer_attestation_policy_root,
        manifest_a.replay_verifier,
    );
    verdict(
        "ac-004",
        "the whole path binds admission, capabilities, query/QoI/units, every interval and \
         authentic verifier receipt, the full receipt-sequence-scoped fixture producer \
         executable attestation, ladder/escalation, cache state/mutations, VoI, package/checker \
         receipts, budgets/cancellation, and source-cone/lock/toolchain inputs; every bound field \
         replays bit-exact and every mutation fails closed (G5)",
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
    let attestation_scope = ProducerAttestationScope {
        purpose: "fs-flywheel-e2e/ac-005-hard-component/v1",
        evidence_root: resolver.receipt.artifact_root().content_hash(),
    };
    let producer_authority = producer_attestation_fixture(&resolver, attestation_scope);
    let ReceiptPromotion::Verified(resolved) = resolve_for_promotion(
        Some(&resolver),
        attestation_scope,
        Some(fixture_producer_attestation_capability(&producer_authority)),
        Some(configured_producer_promotion()),
    ) else {
        panic!("the authentic and producer-attested hard component must resolve");
    };
    let provenance = Provenance::new("acceptance-e2e", "Cargo.lock");
    let hard_claim = resolved.claim();
    let pkg = signed_fixture(
        EvidencePackage::new(provenance.clone())
            .with_claim(hard_claim)
            .with_claim(Claim::estimated("soft", "estimated part", "dwr-guess", 0.1)),
        "acceptance-gate/ac-005",
    );
    let source_verifier = ReceiptCertificateVerifier::from_promoted(&resolved, provenance, 0);
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
