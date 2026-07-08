//! fs-vskeleton — the PV vertical skeleton (patch Rev R): the deliberately
//! tiny end-to-end slice proving FrankenSim's typed-value semantics BEFORE
//! the full substrate exists.
//!
//! One study: 2D SDF plate-with-hole → variable-coefficient diffusion solve
//! (matrix-free FD, CG with cancellation poll points) → compliance+volume
//! objective → ADJOINT gradient verified against central differences →
//! projected gradient descent on the hole radius → fsqlite ledger with
//! content-addressed artifacts → replay comparison → deterministic rerun →
//! text report. All observability flows through the fs-obs event schema.
//!
//! Honest substitutions (each owned by a real bead):
//! - Content hashes are FNV-1a 64 (fs-obs) — BLAKE3-class tree hashing lands
//!   with fs-ledger-core.
//! - The deterministic parallel pattern uses fixed-chunk `std::thread::scope`
//!   with index-ordered pairwise merges; full asupersync two-lane scopes and
//!   ≤200 µs cancel latency land with fs-exec (its Budget vocabulary is
//!   already smoke-tested there).
//! - The "opdsl seed": the edge stencil is defined ONCE ([`model::EdgeLaw`])
//!   and BOTH the primal apply and the adjoint sensitivity contraction are
//!   derived from it — the one-source-of-truth pattern fs-opdsl generalizes.

pub mod ledger;
pub mod model;
pub mod sexpr;

use fs_obs::{Emitter, Event, EventKind, Severity};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Outcome of a full study run: everything the report and the tests need.
#[derive(Debug, Clone)]
pub struct StudyOutcome {
    /// Objective trace, one entry per optimizer iteration.
    pub objective_trace: Vec<f64>,
    /// Design (hole radius) trace, one entry per iteration (post-step).
    pub radius_trace: Vec<f64>,
    /// Gradient-check relative errors per iteration.
    pub gradient_check_rel_err: Vec<f64>,
    /// Total CG iterations spent (the budget accounting).
    pub cg_iterations_spent: u64,
    /// The rendered text report.
    pub report: String,
    /// Content hashes of every artifact written, in write order.
    pub artifact_hashes: Vec<String>,
}

/// Run a study end-to-end, writing ops/artifacts into the ledger at `db_path`.
///
/// # Errors
/// Returns a human/agent-readable error string on parse, budget, solver,
/// gradient-check, or ledger failures (structured errors-as-guidance).
pub fn run_study(study_text: &str, db_path: &str) -> Result<StudyOutcome, String> {
    let spec = model::StudySpec::parse(study_text)?;
    let mut em = Emitter::new(spec.name.clone(), "vskeleton/run".to_string());
    let led = ledger::MiniLedger::open(db_path)?;
    let study_hash = led.put_artifact("study-ir", study_text.as_bytes())?;
    let op_id = led.record_op("study-admitted", study_text, &spec.seed_hex())?;
    led.link(op_id, &study_hash, "in")?;

    let mut outcome = execute(&spec, &mut em)?;

    // Ledger the per-iteration fields and the report.
    for (i, bytes) in outcome_artifacts(&spec, &outcome) {
        let h = led.put_artifact(&format!("iter-{i}"), &bytes)?;
        let op = led.record_op(&format!("solve-iter-{i}"), "", &spec.seed_hex())?;
        led.link(op, &h, "out")?;
        outcome.artifact_hashes.push(h);
    }
    let report_hash = led.put_artifact("report", outcome.report.as_bytes())?;
    let op = led.record_op("report", "", &spec.seed_hex())?;
    led.link(op, &report_hash, "out")?;
    outcome.artifact_hashes.push(report_hash);
    emit(
        &mut em,
        Severity::Info,
        EventKind::Custom {
            name: "study-complete".into(),
            json: format!(
                "{{\"iters\":{},\"final_objective\":{}}}",
                outcome.objective_trace.len(),
                outcome.objective_trace.last().copied().unwrap_or(f64::NAN)
            ),
        },
    );
    Ok(outcome)
}

/// Re-execute the study recorded in the ledger and verify every artifact
/// hash matches (and that stored bytes still hash to their recorded hash —
/// byte corruption fails loudly).
///
/// # Errors
/// Returns a description of the FIRST divergence or corruption found.
pub fn replay(db_path: &str) -> Result<(), String> {
    let led = ledger::MiniLedger::open(db_path)?;
    led.verify_artifact_integrity()?;
    let study_text = led.get_study_ir()?;
    let spec = model::StudySpec::parse(&study_text)?;
    let mut em = Emitter::new(spec.name.clone(), "vskeleton/replay".to_string());
    let outcome = execute(&spec, &mut em)?;
    let mut recomputed: Vec<String> = Vec::new();
    for (_, bytes) in outcome_artifacts(&spec, &outcome) {
        recomputed.push(ledger::content_hash(&bytes));
    }
    recomputed.push(ledger::content_hash(outcome.report.as_bytes()));
    let mut recorded = led.artifact_hashes_excluding_study()?;
    // Compare as MULTISETS. Artifacts are content-addressed (identified by
    // their hash, not by row position), and the ledger's `SELECT` carries no
    // `ORDER BY`, so a positional `!=` could false-report divergence on a
    // faithful replay purely from an fsqlite query-plan reordering. Sorting
    // both sides removes that dependence while still catching any missing,
    // extra, or changed artifact.
    recomputed.sort();
    recorded.sort();
    if recorded != recomputed {
        return Err(format!(
            "replay divergence: recorded {} artifact hashes {:?} but recomputed {:?} — \
             the ledgered study does not reproduce; investigate determinism or tampering",
            recorded.len(),
            recorded,
            recomputed
        ));
    }
    Ok(())
}

/// Execute the numerical study (no ledger I/O) — shared by run and replay.
fn execute(spec: &model::StudySpec, em: &mut Emitter) -> Result<StudyOutcome, String> {
    let mut radius = spec.initial_radius;
    let mut objective_trace = Vec::new();
    let mut radius_trace = Vec::new();
    let mut gradient_errs = Vec::new();
    let mut cg_spent: u64 = 0;
    let cancel = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));

    for step in 0..spec.opt_steps {
        let eval = model::evaluate(spec, radius, &cancel, &mut cg_spent)?;
        if cg_spent > spec.cg_budget {
            return Err(format!(
                "BudgetExhausted: spent {cg_spent} CG iterations of {} budget at optimizer \
                 step {step}; relax (budget (cg-iters ...)) or coarsen (grid ...)",
                spec.cg_budget
            ));
        }
        // Gradient truth gate: adjoint vs central differences, every step.
        let fd = model::central_difference(spec, radius, &cancel, &mut cg_spent)?;
        let denom = eval.gradient.abs().max(fd.abs()).max(1e-12);
        let rel = (eval.gradient - fd).abs() / denom;
        gradient_errs.push(rel);
        emit(
            em,
            Severity::Info,
            EventKind::GradientCheck {
                op: "compliance+volume/d_radius".into(),
                max_rel_err: rel,
                pass: rel < 1e-4,
            },
        );
        if rel >= 1e-4 {
            return Err(format!(
                "GradientCheckFailed at step {step}: adjoint {} vs central-diff {fd} \
                 (rel err {rel:.3e} >= 1e-4) — a solver without a passing gradient check \
                 cannot proceed (plan §8.7)",
                eval.gradient
            ));
        }
        objective_trace.push(eval.objective);
        // Projected gradient step on the radius.
        radius = (radius - spec.step_size * eval.gradient).clamp(spec.r_min, spec.r_max);
        radius_trace.push(radius);
        emit(
            em,
            Severity::Info,
            EventKind::BudgetDelta {
                resource: "cg_iters".into(),
                spent: cg_spent as f64,
                remaining: (spec.cg_budget.saturating_sub(cg_spent)) as f64,
            },
        );
    }

    let report = render_report(
        spec,
        &objective_trace,
        &radius_trace,
        &gradient_errs,
        cg_spent,
    );
    Ok(StudyOutcome {
        objective_trace,
        radius_trace,
        gradient_check_rel_err: gradient_errs,
        cg_iterations_spent: cg_spent,
        report,
        artifact_hashes: Vec::new(),
    })
}

/// The per-iteration artifact bytes (deterministic little-endian encoding of
/// the design + objective — the state a fork/resume would need).
fn outcome_artifacts(spec: &model::StudySpec, o: &StudyOutcome) -> Vec<(usize, Vec<u8>)> {
    let mut out = Vec::new();
    for i in 0..o.objective_trace.len() {
        let mut bytes = Vec::with_capacity(24);
        bytes.extend_from_slice(&(i as u64).to_le_bytes());
        bytes.extend_from_slice(&o.radius_trace[i].to_le_bytes());
        bytes.extend_from_slice(&o.objective_trace[i].to_le_bytes());
        bytes.extend_from_slice(&(spec.grid as u64).to_le_bytes());
        out.push((i, bytes));
    }
    out
}

fn render_report(
    spec: &model::StudySpec,
    objectives: &[f64],
    radii: &[f64],
    grad_errs: &[f64],
    cg_spent: u64,
) -> String {
    use std::fmt::Write as _;
    let mut r = String::new();
    let _ = writeln!(r, "# PV vertical-skeleton report: {}", spec.name);
    let _ = writeln!(r, "seed: {}", spec.seed_hex());
    let _ = writeln!(
        r,
        "grid: {0}x{0} | cg budget: {1} | spent: {cg_spent}",
        spec.grid, spec.cg_budget
    );
    let _ = writeln!(
        r,
        "\n## Objective trace (compliance + {} * volume)",
        spec.volume_weight
    );
    for (i, ((j, rad), ge)) in objectives.iter().zip(radii).zip(grad_errs).enumerate() {
        let _ = writeln!(
            r,
            "step {i}: J = {j:.9e} | radius -> {rad:.6} | adjoint-vs-FD rel err {ge:.3e}"
        );
    }
    let _ = writeln!(r, "\n## Verdicts");
    let _ = writeln!(
        r,
        "gradient checks: {} / {} passed (< 1e-4)",
        grad_errs.len(),
        grad_errs.len()
    );
    let _ = writeln!(
        r,
        "budget: {} of {} CG iterations",
        cg_spent, spec.cg_budget
    );
    let _ = writeln!(
        r,
        "\nreproduce: `fs-vskeleton run <this study file>` — same seed, same hashes."
    );
    r
}

fn emit(em: &mut Emitter, sev: Severity, kind: EventKind) {
    let e: Event = em.emit(sev, kind, None);
    eprintln!("{}", e.to_jsonl());
}
