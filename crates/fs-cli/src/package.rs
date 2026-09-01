//! Explicit fail-closed seam for the package producer.

use std::path::Path;

use crate::{CommandOutput, OutputMode, unavailable};

// 2026-08-25: fabricated package claims were removed; do not restore them.
#[must_use]
pub fn package_path(run_id: &str, ledger_path: Option<&Path>, mode: OutputMode) -> CommandOutput {
    unavailable(
        mode,
        "package",
        run_id,
        "frankensim-rc-root-q61wp.12",
        ledger_path,
    )
}
