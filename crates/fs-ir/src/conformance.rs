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
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome};

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

/// One ledgered conformance run: the standard structured report plus the
/// ledger op ids that now constitute the golden.
pub struct LedgeredRun {
    /// The structured casebook report.
    pub report: SuiteReport,
    /// One ledger op per case, in case order.
    pub op_ids: Vec<i64>,
}

/// Run a conformance suite through the real admission path AND record it
/// as a golden ledger: one finished op per case (frozen IR = the exact
/// canonical program, the suite's Five Explicits, outcome = the case
/// verdict, produced artifacts content-addressed and linked). CI then
/// replays conformance the same way it replays features:
/// `fs_ledger::travel` `replay_verdict` over two runs of the same suite
/// is the one replay/compare mechanism (golden-ledger unification,
/// bead 6nb.8 slice 2).
///
/// Timestamps are caller-supplied logical nanoseconds (`t0_ns + index`),
/// never a clock read; the travel replay contract excludes run-envelope
/// timing from semantic comparison.
///
/// # Errors
/// [`LedgerError`] on any ledger refusal (the run is not silently
/// half-recorded: the first refusal aborts).
pub fn run_ir_suite_ledgered(
    suite_name: &'static str,
    cases: Vec<IrCase>,
    cx: &AdmissionContext<'_>,
    ledger: &Ledger,
    explicits: &FiveExplicits<'_>,
    t0_ns: i64,
) -> Result<LedgeredRun, LedgerError> {
    let mut suite = Suite::new(suite_name);
    let mut op_ids = Vec::new();
    for (index, case) in cases.into_iter().enumerate() {
        let inputs_digest = fnv1a64(case.program_sexpr.as_bytes());
        let (id, tolerance, outcome, admitted, produced) = evaluate_case_ledgered(case, cx);
        let t_start = t0_ns + 2 * index as i64;
        // The frozen op IR is the exact canonical program identity the
        // admission path bound (raw form; the lowered identity rides the
        // record evidence). Refused cases record as Error ops so the
        // golden retains the refusal.
        let ir_json = format!(
            "{{\"conformance\":\"{suite_name}\",\"case\":\"{id}\",\"program_fnv\":\"{inputs_digest:016x}\"}}"
        );
        let op = ledger.begin_op(None, &ir_json, explicits, t_start)?;
        if let Some(bytes) = &produced {
            let artifact = ledger.put_artifact("conformance-artifact", bytes, None)?;
            ledger.link(op, &artifact.hash, EdgeRole::Out)?;
        }
        let op_outcome = if admitted && outcome.pass {
            OpOutcome::Ok
        } else {
            OpOutcome::Error
        };
        let diag = if outcome.pass {
            None
        } else {
            // The ledger requires JSON diagnostics: wrap and escape.
            let mut escaped = String::with_capacity(outcome.details.len() + 16);
            for c in outcome.details.chars() {
                match c {
                    '"' => escaped.push_str("\\\""),
                    '\\' => escaped.push_str("\\\\"),
                    '\n' => escaped.push_str("\\n"),
                    '\r' => escaped.push_str("\\r"),
                    '\t' => escaped.push_str("\\t"),
                    c if (c as u32) < 0x20 => {
                        escaped.push_str(&format!("\\u{:04x}", c as u32));
                    }
                    c => escaped.push(c),
                }
            }
            Some(format!("{{\"detail\":\"{escaped}\"}}"))
        };
        ledger.finish_op(op, op_outcome, diag.as_deref(), t_start + 1)?;
        op_ids.push(op);
        suite = suite.case(id, inputs_digest, tolerance, move || outcome);
    }
    Ok(LedgeredRun {
        report: suite.run(),
        op_ids,
    })
}

fn evaluate_case_ledgered(
    case: IrCase,
    cx: &AdmissionContext<'_>,
) -> (
    &'static str,
    ToleranceSpec,
    CaseOutcome,
    bool,
    Option<Vec<u8>>,
) {
    let expected = case.expected_artifact;
    let id = case.id;
    let tolerance = case.tolerance;
    let node = match crate::sexpr::parse(&case.program_sexpr) {
        Ok(node) => node,
        Err(err) => {
            return (
                id,
                tolerance,
                CaseOutcome::fail(format!("IR program refused before admission: {err}")),
                false,
                None,
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
                false,
                None,
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
            false,
            None,
        );
    }
    let lowered_identity = report
        .lowering
        .lowered_canonical()
        .map(|c| format!("ir-lowered:{:016x}", fnv1a64(c.as_bytes())));
    let bytes = (case.kernel)(&report);
    let produced = artifact_hash(&bytes);
    let mut outcome = if produced == expected {
        CaseOutcome::pass(format!(
            "admitted through the real path ({} findings) and artifact matches {}",
            report.findings.len(),
            produced.to_hex()
        ))
    } else {
        CaseOutcome::fail(format!(
            "artifact drifted: produced {} expected {}",
            produced.to_hex(),
            expected.to_hex()
        ))
    };
    outcome = outcome.with_evidence(raw_identity);
    if let Some(lowered) = lowered_identity {
        outcome = outcome.with_evidence(lowered);
    }
    (id, tolerance, outcome, true, Some(bytes))
}

fn evaluate_case(
    case: IrCase,
    cx: &AdmissionContext<'_>,
) -> (&'static str, ToleranceSpec, CaseOutcome) {
    let (id, tolerance, outcome, _, _) = evaluate_case_ledgered(case, cx);
    (id, tolerance, outcome)
}
