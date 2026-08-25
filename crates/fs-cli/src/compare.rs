//! CLI `compare` command — evidence-aware semantic run differences.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.14.1`
//!
//! The retained-run loader and semantic comparator required by the owning Bead
//! have not landed. This command therefore fails closed instead of emitting
//! fixture values that could be mistaken for evidence from the named runs.

use std::path::Path;

use crate::{CommandOutput, OutputMode, unavailable};

/// Execute the `compare` command comparing two runs.
#[must_use]
pub fn compare_path(
    left_run: &str,
    right_run: &str,
    ledger_override: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    let _ = ledger_override;
    let subject = format!("{left_run}..{right_run}");
    unavailable(
        mode,
        "compare",
        &subject,
        "frankensim-extreal-program-f85xj.6.14.1",
    )
}
