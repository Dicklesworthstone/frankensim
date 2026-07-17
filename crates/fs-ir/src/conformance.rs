//! Conformance suites speak IR (bead frankensim-epic-gauntlet-6nb.8,
//! slice 1): each conformance case exposes an equivalent FrankenScript IR
//! program, and running the suite routes every case through the REAL
//! parse → lower → admit entry path — admission findings, capability
//! decisions, and the exact canonical program identities land in the
//! structured fs-casebook record, so conformance exercises the same door
//! a production study walks through, never a test-only side entrance.
//!
//! An [`IrCase`] binds the stable case id, the IR program source, the
//! declared tolerance, and the expected artifact as a content address.
//! [`run_ir_suite`] admits every program through the supplied
//! [`AdmissionContext`]; a refusal becomes a failing record carrying the
//! full deterministic diagnosis, and an admitted case executes its kernel
//! (which receives the [`AdmissionReport`] so budgets/provenance are in
//! scope) and compares the produced artifact's domain-separated content
//! hash against the expectation. Both canonical identities from the
//! [`LoweringReceipt`] ride the record as evidence pointers — the
//! cross-agent negotiation anchor: two agents agree on the exact
//! canonical program, not on prose.
//!
//! No-claims: fs-ir has no general study executor — execution stays in
//! the case kernel (each domain crate runs its own physics); this slice
//! makes ADMISSION and IDENTITY real, not simulation. Golden-ledger
//! unification via ledger-timetravel and IR-level cross-crate contract
//! checking are follow-on slices tracked on the bead.

use fs_blake3::{ContentHash, hash_domain};
use fs_casebook::{CaseOutcome, Suite, SuiteReport, ToleranceSpec, fnv1a64};

use crate::VersionedProgram;
use crate::admission::{AdmissionContext, AdmissionReport, admit_versioned};

const ARTIFACT_DOMAIN: &str = "org.frankensim.fs-ir.conformance-artifact.v1";

/// Compute the domain-separated content address of one conformance
/// artifact's exact bytes.
#[must_use]
pub fn artifact_hash(bytes: &[u8]) -> ContentHash {
    hash_domain(ARTIFACT_DOMAIN, bytes)
}

/// One conformance case as an IR program plus its expectation.
pub struct IrCase {
    /// Stable case id (the cross-run identity).
    pub id: &'static str,
    /// The FrankenScript program source (s-expression form). Inputs are
    /// IR values inside the program; the Five Explicits are checked by
    /// the real admission path, not restated here.
    pub program_sexpr: String,
    /// The declared tolerance model, recorded verbatim in the record.
    pub tolerance: ToleranceSpec,
    /// The expected artifact content address.
    pub expected_artifact: ContentHash,
    /// The execution kernel: runs ONLY after real admission passes, and
    /// receives the admission report so budgets and provenance are in
    /// scope. Returns the artifact's exact bytes.
    pub kernel: Box<dyn FnOnce(&AdmissionReport) -> Vec<u8>>,
}

/// Run a conformance suite of IR cases through the REAL admission path,
/// producing the standard fs-casebook structured report.
///
/// Per case: parse → version-bind → lower → admit under `cx`. A parse or
/// admission refusal is a failing record carrying the deterministic
/// diagnosis. An admitted case executes its kernel and compares the
/// artifact hash; the raw and lowered canonical identities ride the
/// record as evidence pointers.
#[must_use]
pub fn run_ir_suite(
    suite_name: &'static str,
    cases: Vec<IrCase>,
    cx: &AdmissionContext<'_>,
) -> SuiteReport {
    let mut suite = Suite::new(suite_name);
    for case in cases {
        // Admission runs eagerly (the context borrows locals); the
        // casebook closure then owns the finished outcome.
        let inputs_digest = fnv1a64(case.program_sexpr.as_bytes());
        let outcome = evaluate_case(case, cx);
        let (id, tolerance) = (outcome.0, outcome.1);
        let case_outcome = outcome.2;
        suite = suite.case(id, inputs_digest, tolerance, move || case_outcome);
    }
    suite.run()
}

fn evaluate_case(
    case: IrCase,
    cx: &AdmissionContext<'_>,
) -> (&'static str, ToleranceSpec, CaseOutcome) {
    let id = case.id;
    let tolerance = case.tolerance;
    // Fresh source: syntax parse, then version-bind (persisted replays
    // would enter through the version envelope instead).
    let node = match crate::sexpr::parse(&case.program_sexpr) {
        Ok(node) => node,
        Err(err) => {
            return (
                id,
                tolerance,
                CaseOutcome::fail(format!("IR program refused before admission: {err}")),
            );
        }
    };
    let program = match VersionedProgram::try_current(node) {
        Ok(program) => program,
        Err(err) => {
            return (
                id,
                tolerance,
                CaseOutcome::fail(format!("IR program refused before admission: {err}")),
            );
        }
    };
    let report = admit_versioned(&program, cx);
    let raw_identity = format!(
        "ir-raw:{:016x}",
        fnv1a64(report.lowering.raw_canonical().as_bytes())
    );
    if !report.admitted {
        return (
            id,
            tolerance,
            CaseOutcome::fail(format!(
                "REAL admission refused the conformance program:\n{}",
                report.diagnosis()
            ))
            .with_evidence(raw_identity),
        );
    }
    let lowered_identity = report
        .lowering
        .lowered_canonical()
        .map(|c| format!("ir-lowered:{:016x}", fnv1a64(c.as_bytes())));
    let bytes = (case.kernel)(&report);
    let produced = artifact_hash(&bytes);
    let mut outcome = if produced == case.expected_artifact {
        CaseOutcome::pass(format!(
            "admitted through the real path ({} findings) and artifact matches {}",
            report.findings.len(),
            produced.to_hex()
        ))
    } else {
        CaseOutcome::fail(format!(
            "artifact drifted: produced {} expected {}",
            produced.to_hex(),
            case.expected_artifact.to_hex()
        ))
    };
    outcome = outcome.with_evidence(raw_identity);
    if let Some(lowered) = lowered_identity {
        outcome = outcome.with_evidence(lowered);
    }
    (id, tolerance, outcome)
}
