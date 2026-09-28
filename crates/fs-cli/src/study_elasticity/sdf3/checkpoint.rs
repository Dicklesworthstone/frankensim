//! Durable stage boundaries with exact, budget-charged prefix replay.
//! Retained JSON is evidence to reproduce, never a source of solver vectors.
use super::*;
use fs_blake3::DomainHasher;

pub(super) const MODE: &str = "verified-stage-replay-v1";
pub(super) const KIND: &str = "study-sdf3-checkpoint";
pub(super) const SCHEMA: &str = "sdf3-stage-state.v1";

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Spent {
    pub(super) wall_s: f64,
    pub(super) linear: SolveWork,
    pub(super) geometry: QuadratureWork3,
}

fn malformed(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-checkpoint", message)
}

impl Spent {
    pub(super) fn validate(&self, spec: &Spec) -> Result<()> {
        if !self.wall_s.is_finite() || self.wall_s < 0.0
            || self.linear.linear_iterations > spec.linear
            || self.geometry.boxes > spec.boxes || self.geometry.points > spec.points
        {
            return Err(malformed("retained consumption is invalid or exceeds a discrete work limit"));
        }
        Ok(())
    }

    pub(super) fn add(&self, linear: SolveWork, geometry: QuadratureWork3, wall_s: f64) -> Result<Self> {
        let add = |a: usize, b: usize| a.checked_add(b)
            .ok_or_else(|| malformed("cumulative work counter overflow"));
        let total_wall = self.wall_s + wall_s;
        if !wall_s.is_finite() || wall_s < 0.0 || !total_wall.is_finite() {
            return Err(malformed("invalid cumulative wall consumption"));
        }
        Ok(Self {
            wall_s: total_wall,
            linear: SolveWork {
                linear_iterations: add(self.linear.linear_iterations, linear.linear_iterations)?,
                linear_solves: add(self.linear.linear_solves, linear.linear_solves)?,
                preconditioner_operator_applications: add(self.linear.preconditioner_operator_applications,
                    linear.preconditioner_operator_applications)?,
                preconditioner_galerkin_products: add(self.linear.preconditioner_galerkin_products,
                    linear.preconditioner_galerkin_products)?,
            },
            geometry: QuadratureWork3 {
                boxes: add(self.geometry.boxes, geometry.boxes)?,
                points: add(self.geometry.points, geometry.points)?,
                field_evaluations: add(self.geometry.field_evaluations, geometry.field_evaluations)?,
            },
        })
    }

    fn read(receipt: &JsonValue) -> Result<Self> {
        let work = receipt.get("work").ok_or_else(|| malformed("missing cumulative work"))?;
        Ok(Self {
            wall_s: receipt.get("consumed_wall_s").and_then(JsonValue::as_f64)
                .ok_or_else(|| malformed("missing wall consumption"))?,
            linear: SolveWork {
                linear_iterations: integer(work, "linear_iterations")?,
                linear_solves: integer(work, "linear_solves")?,
                preconditioner_operator_applications: integer(work, "preconditioner_operator_applications")?,
                preconditioner_galerkin_products: integer(work, "preconditioner_galerkin_products")?,
            },
            geometry: QuadratureWork3 {
                boxes: integer(work, "geometry_boxes")?,
                points: integer(work, "geometry_points")?,
                field_evaluations: integer(work, "geometry_field_evaluations")?,
            },
        })
    }

    pub(super) fn json(&self) -> String {
        format!(
            "{{\"linear_iterations\":{},\"linear_solves\":{},\"preconditioner_operator_applications\":{},\"preconditioner_galerkin_products\":{},\"geometry_boxes\":{},\"geometry_points\":{},\"geometry_field_evaluations\":{}}}",
            self.linear.linear_iterations, self.linear.linear_solves,
            self.linear.preconditioner_operator_applications, self.linear.preconditioner_galerkin_products,
            self.geometry.boxes, self.geometry.points, self.geometry.field_evaluations,
        )
    }

    fn require_remaining(&self, spec: &Spec) -> Result<()> {
        self.validate(spec)?;
        if self.wall_s >= spec.wall_s || self.linear.linear_iterations >= spec.linear
            || self.geometry.boxes >= spec.boxes || self.geometry.points >= spec.points
        {
            return Err(Failure {
                code: "cli-study-sdf3-resume-budget",
                message: "original study allowance is exhausted; resume does not renew work budgets".into(),
                exit: exit::BUDGET,
            });
        }
        Ok(())
    }
}

pub(super) struct Retention {
    pub(super) producer: ContentHash,
    pub(super) predecessor: Option<ContentHash>,
    pub(super) replayed_stages: usize,
}

fn cancelled(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(Failure { code: "cli-study-sdf3-cancelled", message: "study cancelled before publication".into(), exit: exit::CANCELLED })
    } else { Ok(()) }
}

/// Bind all linked code and numerical dependencies, not only this adapter's
/// version string. The running executable is read in bounded, cancellable tiles.
fn producer_identity(gate: &CancelGate) -> Result<ContentHash> {
    #[cfg(target_os = "linux")]
    let path = std::path::PathBuf::from("/proc/self/exe");
    #[cfg(not(target_os = "linux"))]
    let path = std::env::current_exe().map_err(|error| malformed(error.to_string()))?;
    let mut file = std::fs::File::open(path).map_err(|error| malformed(error.to_string()))?;
    let mut hasher = DomainHasher::new("org.frankensim.sdf3-study.executable.v1");
    let mut buffer = [0_u8; 65_536];
    let mut total = 0_u64;
    loop {
        cancelled(gate)?;
        let count = file.read(&mut buffer).map_err(|error| malformed(error.to_string()))?;
        if count == 0 { break; }
        total += count as u64;
        if total > 4 * 1024 * 1024 * 1024 {
            return Err(malformed("executable exceeds the 4 GiB identity envelope"));
        }
        hasher.update(&buffer[..count]);
    }
    cancelled(gate)?;
    Ok(hasher.finalize())
}

struct Expected {
    stages: usize,
    state: ContentHash,
    spent: Spent,
}

fn expected(spec: &Spec, ledger: &Ledger, old: &Loaded, producer: ContentHash) -> Result<Expected> {
    let value = &old.value;
    if value.str_field("resume_mode") != Some(MODE)
        || value.get("resume_supported") != Some(&JsonValue::Bool(true))
    {
        return Err(fail("cli-study-sdf3-resume-unsupported",
            "this receipt has no complete replayable 3-D stage; retain its report and use a completed-stage checkpoint"));
    }
    if value.str_field("driver") != Some(SDF3_DRIVER)
        || value.str_field("study_id") != Some(spec.id.to_hex().as_str())
        || value.str_field("producer") != Some(producer.to_hex().as_str())
        || integer(value, "target_stages")? != spec.schedule.len()
        || integer(value, "target_iterations")? != spec.updates * spec.schedule.len()
    {
        return Err(malformed("checkpoint source, executable or original target changed"));
    }
    let stages = integer(value, "stages_completed")?;
    if stages == 0 || stages > spec.schedule.len()
        || integer(value, "iterations_completed")? > stages * spec.updates
        || ((value.str_field("status") == Some("completed")) != (stages == spec.schedule.len()))
        || !matches!(value.str_field("status"),
            Some("completed" | "checkpointed" | "budget-exhausted" | "cancelled" | "numerical-failure"))
    {
        return Err(malformed("checkpoint stage count or terminal status disagrees with its original schedule"));
    }
    let state = value.str_field("checkpoint").and_then(ContentHash::from_hex)
        .ok_or_else(|| malformed("missing exact stage-state artifact"))?;
    let op = ledger.artifact_output_seal(&old.hash)?
        .ok_or_else(|| malformed("unsealed checkpoint receipt"))?;
    if !ledger.edge_exists(op, &state, EdgeRole::Out)? {
        return Err(malformed("stage-state artifact is not an output of the sealed receipt operation"));
    }
    // Read-integrity and kind admission happen before replay. The blob's entire
    // canonical spelling is then independently reproduced, not trusted as state.
    let bytes = artifact(ledger, state, KIND)?;
    let doc = JsonValue::parse(std::str::from_utf8(&bytes)
        .map_err(|error| malformed(error.to_string()))?)
        .map_err(|error| malformed(error.to_string()))?;
    if doc.str_field("schema") != Some(SCHEMA) || integer(&doc, "stages_completed")? != stages
        || doc.str_field("design") != value.str_field("design")
        || doc.str_field("iterations") != value.str_field("iterations")
    {
        return Err(malformed("stage-state artifact does not bind the retained fields and history"));
    }
    if !matches!(value.get("predecessor"), Some(JsonValue::Null)) && value.str_field("predecessor").is_none() {
        return Err(malformed("checkpoint predecessor must be a hash or null"));
    }
    if let Some(previous) = value.str_field("predecessor") {
        let previous = ContentHash::from_hex(previous).ok_or_else(|| malformed("invalid predecessor"))?;
        if previous == old.hash || !ledger.edge_exists(op, &previous, EdgeRole::In)? {
            return Err(malformed("checkpoint predecessor is not linked input evidence"));
        }
        artifact(ledger, previous, RECEIPT_KIND)?;
    }
    let spent = Spent::read(value)?;
    spent.validate(spec)?;
    Ok(Expected { stages, state, spent })
}

/// Progress follows commit, so an interrupted process has a usable receipt ID.
/// A closed diagnostic stream does not undo an already committed checkpoint.
pub(super) fn announce(out: &Outcome, stages: usize) -> Result<()> {
    use std::io::Write as _;
    let _ = writeln!(std::io::stderr().lock(),
        "{{\"schema\":\"frankensim.cli.sdf3-progress.v1\",\"run_id\":{},\"stages_completed\":{stages}}}",
        quoted(&out.pointer));
    Ok(())
}

pub(super) fn drive(
    spec: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
    old: Option<&Loaded>,
    mut published: impl FnMut(&Outcome, usize) -> Result<()>,
) -> Result<Outcome> {
    let admission = Instant::now();
    cancelled(gate)?;
    let producer = producer_identity(gate)?;
    let expected = old.map(|old| expected(spec, ledger, old, producer)).transpose()?;
    if let (Some(old), Some(expected)) = (old, &expected) {
        if old.value.str_field("status") == Some("completed") && expected.stages == spec.schedule.len() {
            return Ok(Outcome { pointer: format!("study-{}", old.hash.to_hex()), receipt: old.bytes.clone(), status: "completed" });
        }
        expected.spent.require_remaining(spec)?;
    }
    let before = expected.as_ref().map_or(0, |old| old.stages);
    let prior = expected.as_ref().map_or(Spent::default(), |old| old.spent)
        .add(SolveWork::default(), QuadratureWork3::default(), admission.elapsed().as_secs_f64())?;
    let requested = before.saturating_add(cap.unwrap_or(spec.schedule.len())).min(spec.schedule.len());
    let mut verified = expected.is_none();
    let mut retention = Retention { producer, predecessor: old.map(|old| old.hash), replayed_stages: before };
    let mut latest: Option<Outcome> = None;
    let result = compute_observed(spec, gate, prior, requested, |view| {
        cancelled(gate)?;
        if view.completed_stages <= before {
            if view.completed_stages == before {
                let manifest = output::state(view)?;
                if Some(hash_bytes(manifest.as_bytes())) != expected.as_ref().map(|old| old.state) {
                    return Err(fail("cli-study-sdf3-replay-mismatch",
                        "replayed mesh, native fields, densities, history or installed refinement evidence differs from the checkpoint"));
                }
                verified = true;
            }
            return Ok(());
        }
        if !verified { return Err(malformed("cannot publish before verifying the retained stage prefix")); }
        let out = output::persist(spec, ledger, view, &retention)?;
        retention.predecessor = out.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        latest = Some(out);
        published(latest.as_ref().expect("just persisted"), view.completed_stages)
    });
    let result = result.map_err(|mut error| {
        if let Some(out) = &latest {
            error.message.push_str(&format!("; last durable checkpoint: {}", out.pointer));
        } else if let Some(old) = old {
            error.message.push_str(&format!("; retained checkpoint unchanged: study-{}", old.hash.to_hex()));
        }
        error
    })?;
    if !verified {
        return Err(Failure {
            code: "cli-study-sdf3-replay-incomplete",
            message: format!("retained prefix could not be reproduced within the remaining allowance; no checkpoint was replaced ({}); retained run study-{}",
                result.status, old.expect("unverified resume").hash.to_hex()),
            exit: if gate.is_requested() { exit::CANCELLED }
                else if result.status == "budget-exhausted" { exit::BUDGET } else { exit::REFUSED },
        });
    }
    if result.report.continuation.termination == ContinuationTermination::ScheduleComplete {
        return latest.ok_or_else(|| malformed("completed execution produced no durable stage"));
    }
    // Terminal failures retain the previously accepted model, plus all spent
    // work and any rejected refinement evidence. A failed initial stage has no
    // resumable stage but still exports its honest partial numerical result.
    output::persist(spec, ledger, &result.view(), &retention)
}

#[cfg(test)]
#[path = "checkpoint_tests.rs"]
mod tests;
