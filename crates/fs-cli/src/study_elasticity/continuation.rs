//! Native study continuation over the existing fs-topols checkpoint and ledger.
//! No second optimizer or checkpoint file format: the design already contains
//! exact phi bits and the final iteration row already contains the multiplier.
//! A new receipt binds direct continuation to this executable. Legacy receipts
//! are reconstructed once by replay, never once per subsequent update.
use super::*;
use fs_blake3::DomainHasher;
use fs_topols::OptimizeCheckpoint;
use fs_topols::checkpoint::CheckpointStage;
use std::ops::ControlFlow;

// Bounds additional CG iterations, not wall-clock latency of assembly,
// smoothing, advection or an individual vector operation.
const CG_POLL_ITERS: usize = 32;

#[path = "continuation/projected.rs"]
mod projected;
pub(super) use projected::{Controls as ProjectedControls, PROJECTED_SCOPE, parse_controls};

fn stop_status(cancelled: bool, consumed_wall: f64, wall_limit: f64) -> Option<&'static str> {
    if cancelled { Some("cancelled") }
    else if consumed_wall >= wall_limit { Some("budget-exhausted") }
    else { None }
}

pub(super) struct Evidence {
    producer: ContentHash,
    updates: usize,
    legacy_replayed: usize,
    projected: Option<projected::ConstraintEvidence>,
}

impl Evidence {
    pub(super) fn json(&self) -> String {
        let constraints = self.constraint_fields();
        format!("{{\"version\":1,\"producer\":\"{}\",\"updates_this_invocation\":{},\"legacy_prefix_updates_replayed\":{},\"mode\":\"accepted-state-continuation\"{constraints}}}",
            self.producer.to_hex(), self.updates, self.legacy_replayed)
    }
}

impl Evidence {
    pub(super) fn constraint_fields(&self) -> String {
        self.projected.as_ref().map_or_else(String::new,
            |state| format!(",\"constraints\":{}", state.json()))
    }

    pub(super) fn projected_current(&self) -> Option<&fs_topols::SampledStressEvaluation> {
        self.projected.as_ref().map(projected::ConstraintEvidence::current)
    }

    pub(super) fn constraint_html(&self) -> String {
        self.projected.as_ref().map_or_else(String::new, |state| state.html())
    }
}

fn producer_identity() -> Result<ContentHash> {
    #[cfg(target_os = "linux")]
    let path = std::path::PathBuf::from("/proc/self/exe");
    #[cfg(not(target_os = "linux"))]
    let path = std::env::current_exe()
        .map_err(|error| fail("cli-study-elasticity-producer", error.to_string()))?;
    let mut file = std::fs::File::open(path)
        .map_err(|error| fail("cli-study-elasticity-producer", error.to_string()))?;
    let mut hash = DomainHasher::new("org.frankensim.elasticity-study.executable.v1");
    let mut buffer = [0_u8; 65_536];
    loop {
        let n = file.read(&mut buffer)
            .map_err(|error| fail("cli-study-elasticity-producer", error.to_string()))?;
        if n == 0 { break; }
        hash.update(&buffer[..n]);
    }
    Ok(hash.finalize())
}

fn malformed(what: &str) -> Failure { fail("cli-study-elasticity-continuation", what) }

fn document(bytes: &[u8]) -> Result<JsonValue> {
    let text = std::str::from_utf8(bytes).map_err(|_| malformed("retained state is not UTF-8"))?;
    JsonValue::parse(text).map_err(|error| malformed(&error.to_string()))
}

fn hexadecimal(value: &str, prefix: bool) -> Result<u64> {
    let digits = if prefix { value.strip_prefix("0x") } else { Some(value) }
        .ok_or_else(|| malformed("retained snapshot has no hexadecimal prefix"))?;
    if digits.len() != 16 || !digits.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(malformed("retained state requires exactly sixteen hexadecimal digits"));
    }
    u64::from_str_radix(digits, 16).map_err(|_| malformed("invalid retained state bits"))
}

fn number<'a>(row: &'a JsonValue, key: &str) -> Result<(f64, &'a str)> {
    let value = row.get(key).ok_or_else(|| malformed("missing iteration value"))?;
    let number = value.as_f64().filter(|v| v.is_finite())
        .ok_or_else(|| malformed("non-finite or nonnumeric iteration value"))?;
    let raw = value.number_raw().ok_or_else(|| malformed("iteration value has no numeric spelling"))?;
    Ok((number, raw))
}

/// Decode only the bounded, already linked artifacts. Preserve each numeric
/// spelling when rebuilding row bytes, so the existing trace hash remains an
/// EXACT check, not a tolerance comparison or a re-rounded JSON digest.
fn decode(spec: &ElasticitySpec, receipt: &JsonValue, design: &JsonValue,
    iterations: &JsonValue) -> Result<(GridSdf, OptimizeReport)> {
    let count = integer(receipt, "iterations_completed")?;
    if count > spec.steps || integer(receipt, "target_iterations")? != spec.steps
        || (receipt.str_field("status") == Some("completed") && count != spec.steps) {
        return Err(malformed("retained iteration count disagrees with the study"));
    }
    if iterations.str_field("schema") != Some("elasticity-study-iterations-v1")
        || iterations.str_field("study_id") != Some(spec.id.to_hex().as_str()) {
        return Err(malformed("iteration artifact belongs to another study or schema"));
    }
    let n = integer(design, "n")?;
    let expected_n = 1usize << spec.base.physics.as_ref().expect("admitted physics").mesh_level;
    if n != expected_n { return Err(malformed("retained lattice differs from the admitted mesh level")); }
    let bits = design.get("phi_bits").and_then(JsonValue::as_array)
        .ok_or_else(|| malformed("missing exact retained level-set bits"))?;
    if bits.len() != (n + 1) * (n + 1) {
        return Err(malformed("retained level-set node count is wrong"));
    }
    let mut phi = GridSdf::from_fn(n, &|_, _| 0.0);
    for (slot, bits) in phi.nodes_mut().iter_mut().zip(bits) {
        *slot = f64::from_bits(hexadecimal(bits.as_str()
            .ok_or_else(|| malformed("level-set bits must be strings"))?, false)?);
        if !slot.is_finite() { return Err(malformed("retained level set is non-finite")); }
    }
    let actual_snapshot = snapshot(&phi);
    let expected_snapshot = hexadecimal(design.str_field("snapshot")
        .ok_or_else(|| malformed("missing design snapshot"))?, true)?;
    if expected_snapshot != actual_snapshot { return Err(malformed("design snapshot does not match phi bits")); }
    let rows = iterations.get("iterations").and_then(JsonValue::as_array)
        .ok_or_else(|| malformed("missing retained iteration rows"))?;
    if rows.len() != count { return Err(malformed("receipt and iteration artifact lengths disagree")); }
    let mut report = OptimizeReport::default();
    for (index, row) in rows.iter().enumerate() {
        if integer(row, "iter")? != index { return Err(malformed("iteration ordinals are not contiguous")); }
        let (compliance, c) = number(row, "compliance")?;
        let (volume, v) = number(row, "volume")?;
        let (ell, e) = number(row, "ell")?;
        let (drift, d) = number(row, "drift_h")?;
        if compliance < 0.0 || volume <= 0.0 || ell < 0.0 || drift < 0.0 {
            return Err(malformed("invalid physical or optimizer iteration value"));
        }
        let pad = integer(row, "load_pad_nodes")?;
        let snap = hexadecimal(row.str_field("snapshot")
            .ok_or_else(|| malformed("missing iteration snapshot"))?, true)?;
        report.rows.push(format!("{{\"iter\":{index},\"compliance\":{c},\"volume\":{v},\"ell\":{e},\"drift_h\":{d},\"load_pad_nodes\":{pad},\"snapshot\":\"{snap:#018x}\"}}"));
        report.compliance.push(compliance);
        report.volume.push(volume);
        report.ell.push(ell);
        report.snapshots.push(snap);
        report.load_pad_nodes.push(pad);
    }
    // Detailed redistance/event objects were never serialized by this CLI.
    // Their retained row projection is preserved verbatim; do not fabricate
    // replacement audit/event objects while loading the reporting arrays.
    if receipt.str_field("trace_hash") != Some(trace_hash(&report.rows).to_hex().as_str()) {
        return Err(malformed("retained iteration bytes disagree with the trace hash"));
    }
    if report.snapshots.last().is_some_and(|&last| last != actual_snapshot) {
        return Err(malformed("final iteration and retained geometry are different designs"));
    }
    if count == 0 && spec.projected.is_none() && phi.nodes().iter().zip(initial_phi(spec).nodes())
        .any(|(a, b)| a.to_bits() != b.to_bits()) {
        return Err(malformed("zero-update geometry differs from the declared initial design"));
    }
    Ok((phi, report))
}

fn append(target: &mut OptimizeReport, mut step: OptimizeReport) {
    target.rows.append(&mut step.rows);
    target.compliance.append(&mut step.compliance);
    target.volume.append(&mut step.volume);
    target.ell.append(&mut step.ell);
    target.snapshots.append(&mut step.snapshots);
    target.load_pad_nodes.append(&mut step.load_pad_nodes);
    target.audits.append(&mut step.audits);
    target.events.append(&mut step.events);
}

fn retained_error(mut error: Failure, last: Option<&Outcome>) -> Failure {
    if let Some(last) = last {
        error.message.push_str(&format!("; last accepted state is {}", last.pointer));
    }
    error
}

pub(super) fn drive(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>) -> Result<Outcome> {
    if spec.projected.is_some() {
        projected::drive(spec, ledger, cap, gate, prior)
    } else {
        drive_observed(spec, ledger, cap, gate, prior, |_, _| {})
    }
}

// The observer allows deterministic request injection at real kernel
// boundaries. It cannot supply fields, objectives, stop statuses or receipts.
fn drive_observed(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>,
    mut observe: impl FnMut(usize, CheckpointStage)) -> Result<Outcome> {
    if ledger.in_transaction() {
        return Err(fail("cli-study-elasticity-transaction", "study requires its own ledger transaction"));
    }
    let start = Instant::now();
    let mut evidence = Evidence { producer: producer_identity()?, updates: 0, legacy_replayed: 0, projected: None };
    let mut predecessor = prior.map(|loaded| loaded.hash);
    let mut retained_wall = 0.0;
    let mut last = None;
    let (phi, mut report) = match prior {
        None => (initial_phi(spec), OptimizeReport::default()),
        Some(loaded) => {
            if loaded.value.str_field("study_id") != Some(spec.id.to_hex().as_str()) {
                return Err(malformed("retained study identity changed"));
            }
            retained_wall = loaded.value.f64_field("consumed_wall_s")
                .filter(|v| v.is_finite() && *v >= 0.0)
                .ok_or_else(|| malformed("invalid retained wall charge"))?;
            let design = document(&linked(ledger, &loaded.value, "design", "study-design")?)?;
            let iterations = document(&linked(ledger, &loaded.value, "iterations", "study-iterations")?)?;
            let retained = decode(spec, &loaded.value, &design, &iterations)?;
            if let Some(binding) = loaded.value.get("continuation") {
                if integer(binding, "version")? != 1
                    || binding.str_field("producer") != Some(evidence.producer.to_hex().as_str()) {
                    return Err(malformed("direct continuation requires the identical executable; prior artifacts are unchanged"));
                }
            } else if !retained.1.rows.is_empty() && retained.1.rows.len() < spec.steps {
                // Older receipts do not bind executable bytes. Preserve their
                // previous replay validation ONCE before permitting extension.
                if gate.is_requested() || retained_wall + start.elapsed().as_secs_f64() >= spec.wall_s {
                    return Err(Failure { code: "cli-study-elasticity-resume-budget",
                        message: format!("legacy replay not started; retained state is study-{}", loaded.hash.to_hex()),
                        exit: if gate.is_requested() { exit::CANCELLED } else { exit::BUDGET } });
                }
                let replay = run_prefix(spec, retained.1.rows.len())?;
                if replay.1.rows != retained.1.rows || snapshot(&replay.0) != snapshot(&retained.0)
                    || replay.0.nodes().iter().zip(retained.0.nodes()).any(|(a,b)| a.to_bits() != b.to_bits()) {
                    return Err(malformed("legacy replay does not reproduce the retained design and trace"));
                }
                evidence.legacy_replayed = retained.1.rows.len();
            }
            let status = match loaded.value.str_field("status") {
                Some("completed") => "completed",
                Some("cancelled") => "cancelled",
                Some("budget-exhausted") => "budget-exhausted",
                Some("running") => "running",
                _ => return Err(malformed("unknown retained study status")),
            };
            last = Some(Outcome { pointer: format!("study-{}", loaded.hash.to_hex()),
                receipt: loaded.bytes.clone(), status });
            retained
        }
    };
    let completed = report.rows.len();
    if completed == spec.steps && last.as_ref().is_some_and(|out| out.status == "completed") {
        return last.ok_or_else(|| malformed("complete study has no retained receipt"));
    }
    if retained_wall >= spec.wall_s {
        return Err(retained_error(Failure { code: "cli-study-elasticity-resume-budget",
            message: "retained wall charge exhausts the declared budget".into(), exit: exit::BUDGET }, last.as_ref()));
    }
    let ell = report.ell.last().copied().unwrap_or(spec.ell0);
    let mut state = OptimizeCheckpoint::restore(phi, fixture(spec), settings(spec, spec.steps), completed, ell)
        .map_err(|error| fail("cli-study-elasticity-solve", format!("{error:?}")))?;
    let target = spec.steps.min(completed.saturating_add(cap.unwrap_or(spec.steps - completed)));
    // Retain the seed too: a refusal on the FIRST solve still has an explicit
    // recoverable design, rather than silently dropping the entire study.
    if last.is_none() {
        let initial = persist(spec, ledger, state.geometry(), &report, "running",
            retained_wall + start.elapsed().as_secs_f64(), predecessor, &evidence)?;
        predecessor = initial.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        last = Some(initial);
    }
    loop {
        let status = if gate.is_requested() { "cancelled" }
            else if retained_wall + start.elapsed().as_secs_f64() >= spec.wall_s { "budget-exhausted" }
            else if state.next_iteration() == target {
                if state.is_complete() { "completed" } else { "budget-exhausted" }
            } else { "running" };
        if status != "running" {
            return persist(spec, ledger, state.geometry(), &report, status,
                retained_wall + start.elapsed().as_secs_f64(), predecessor, &evidence)
                .map_err(|error| retained_error(error, last.as_ref()));
        }
        let ordinal = state.next_iteration();
        let step = state.advance_one_controlled(CG_POLL_ITERS, |stage| {
            observe(ordinal, stage);
            match stop_status(gate.is_requested(),
                retained_wall + start.elapsed().as_secs_f64(), spec.wall_s) {
                Some(status) => ControlFlow::Break(status),
                None => ControlFlow::Continue(()),
            }
        }).map_err(|error| retained_error(
            fail("cli-study-elasticity-solve", format!("{error:?}")), last.as_ref()))?;
        let step = match step {
            ControlFlow::Continue(Some(step)) => step,
            ControlFlow::Continue(None) => return Err(malformed("optimizer completed before the declared target")),
            ControlFlow::Break(status) => {
                // The checkpoint owns rollback of the UNPUBLISHED candidate.
                // Persist only the unchanged accepted geometry and prior rows;
                // work spent on the discarded attempt still consumes wall time.
                return persist(spec, ledger, state.geometry(), &report, status,
                    retained_wall + start.elapsed().as_secs_f64(), predecessor, &evidence)
                    .map_err(|error| retained_error(error, last.as_ref()));
            }
        };
        append(&mut report, step);
        evidence.updates += 1;
        let accepted = persist(spec, ledger, state.geometry(), &report, "running",
            retained_wall + start.elapsed().as_secs_f64(), predecessor, &evidence)
            .map_err(|error| retained_error(error, last.as_ref()))?;
        predecessor = accepted.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        last = Some(accepted);
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod interruption;
