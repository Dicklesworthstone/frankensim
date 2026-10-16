//! Execute exactly the native model admitted during program binding.
//!
//! This owns the existing Model, including its immutable geometry and card
//! bytes. Neither execution nor pin checking reopens the source paths. The
//! existing sampler, physical driver, receipt identity and recovery are reused.

use std::path::Path;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes};
use fs_exec::CancelGate;
use fs_ledger::Ledger;
use fs_project::{DecodedProject, uncertainty::BoundStudy};

use super::{Model, invalid};
use super::super::{Execution, Result, budget, drive, fail, output_error, render};
use crate::{CommandOutput, OutputMode};

/// Constraints checked against the objects that will actually execute.
#[derive(Debug, Clone, Copy)]
pub(crate) struct StudyPins {
    pub(crate) project: ContentHash,
    pub(crate) source: Option<ContentHash>,
    pub(crate) wall_seconds: Option<f64>,
}

impl StudyPins {
    fn check(self, bound: &BoundStudy, base: &DecodedProject) -> Result<()> {
        if base.hash() != self.project {
            return Err(invalid("native study's physical project differs from the cooling.project binding"));
        }
        let source = hash_bytes(bound.study().canonical().as_bytes());
        if self.source.is_some_and(|pin| pin != source) {
            return Err(invalid(format!("native study hashes to {}, not its pinned :hash", source.to_hex())));
        }
        if self.wall_seconds.is_some_and(|wall|
            !wall.is_finite() || wall <= 0.0 || wall < bound.study().wall_seconds())
        {
            return Err(invalid("native study exceeds the program's wall allowance; :budget limits evaluations, not wall time"));
        }
        Ok(())
    }
}

/// An admitted, bounded model snapshot. No input path survives into execution.
pub(crate) struct PreparedStudy {
    model: Model,
    input_bytes: u64,
}

impl std::fmt::Debug for PreparedStudy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PreparedStudy")
            .field("identity", &self.model.identity().to_hex())
            .field("input_bytes", &self.input_bytes)
            .finish_non_exhaustive()
    }
}

impl PreparedStudy {
    /// The caller's remaining input-storage allowance also bounds asset reads.
    /// Pin/support checks run before geometry or card resources are consumed.
    pub(crate) fn load(path: &Path, pins: StudyPins, input_limit: u64) -> std::result::Result<Self, String> {
        let model = Model::load_checked(path, input_limit, |bound, base| pins.check(bound, base))
            .map_err(|error| error.to_string())?;
        // Includes the joint law's numerical PSD and executor-specific gates,
        // which syntax-only study binding does not establish.
        Execution::new(&model).map_err(|error| error.to_string())?;
        let input_bytes = model.input_bytes().map_err(|error| error.to_string())?;
        Ok(Self { model, input_bytes })
    }

    pub(crate) fn input_bytes(&self) -> u64 { self.input_bytes }

    /// Reusing a shared snapshot still checks every clause's own pin/grant.
    pub(crate) fn check(&self, pins: StudyPins) -> std::result::Result<(), String> {
        pins.check(&self.model.bound, &self.model.base).map_err(|error| error.to_string())
    }

    pub(crate) fn run(&self, ledger_path: &Path, override_text: Option<&str>, mode: OutputMode) -> CommandOutput {
        let result = (|| {
            let cap = budget(override_text)?;
            let ledger = Ledger::open(ledger_path.to_str()
                .ok_or_else(|| fail("cli-uncertainty-ledger", "ledger path is not UTF-8"))?)?;
            drive(&self.model, &ledger, cap, &CancelGate::new_clock_free(), Instant::now(), None)
        })();
        match result {
            Ok(out) => render(out, mode),
            Err(error) => output_error("study", mode, error),
        }
    }
}

#[cfg(test)]
#[path = "prepared/tests.rs"]
mod tests;
