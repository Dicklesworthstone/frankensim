//! Explicit fail-closed seam for the report producer.

use std::path::Path;

use crate::{CommandOutput, OutputMode, unavailable};

// 2026-08-25: fabricated report values were removed; do not restore them.
#[must_use]
pub fn report_path(run_id: &str, ledger_path: Option<&Path>, mode: OutputMode) -> CommandOutput {
    unavailable(
        mode,
        "report",
        run_id,
        "frankensim-rc-root-q61wp.12",
        ledger_path,
    )
}
