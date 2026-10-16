//! Read-only projection of the original native-UQ receipts. Their stored
//! checkpoints predate resumable model admission and are never reinterpreted.

use std::path::Path;

use fs_blake3::ContentHash;
use fs_ledger::EdgeRole;

use super::{Result, artifact, fail, hash_field, integer, output_error, parse, quoted};
use crate::{CommandOutput, OutputMode, exit};

const PREFIX: &str = "study-uq-";
const DRIVER: &str = "native-cooling-uncertainty-v1";

pub(super) fn owns(pointer: &str) -> bool {
    pointer.starts_with(PREFIX)
}

pub(super) fn resume(mode: OutputMode) -> CommandOutput {
    output_error(
        "study",
        mode,
        fail(
            "cli-uncertainty-legacy-resume",
            "this legacy native-UQ revision retained checkpoints but did not admit resumption; its original receipt remains available through report/package; start a new study to use native recovery",
        ),
    )
}

pub(super) fn export(
    command: &'static str,
    pointer: &str,
    path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    let (ledger, _) = match crate::report::open_export_ledger(command, pointer, path, mode) {
        Ok(opened) => opened,
        Err(output) => return output,
    };
    let result: Result<CommandOutput> = (|| {
        let hash = pointer
            .strip_prefix(PREFIX)
            .and_then(ContentHash::from_hex)
            .ok_or_else(|| {
                fail(
                    "cli-uncertainty-legacy",
                    "invalid legacy native study pointer",
                )
            })?;
        let bytes = artifact(
            &ledger,
            hash,
            "native-uncertainty-receipt",
            16 * 1024 * 1024,
        )?;
        let value = parse(&bytes)?;
        let op_id = ledger.artifact_output_seal(&hash)?.ok_or_else(|| {
            fail(
                "cli-uncertainty-legacy",
                "legacy native receipt is not sealed",
            )
        })?;
        let op = ledger.op(op_id)?.ok_or_else(|| {
            fail(
                "cli-uncertainty-legacy",
                "legacy receipt producer is missing",
            )
        })?;
        let model = hash_field(&value, "model_identity")?;
        let n = integer(&value, "samples_accepted")?;
        let expected_ir = format!(
            "{{\"driver\":{DRIVER:?},\"model\":{},\"accepted\":{n}}}",
            quoted(&model.to_hex())
        );
        if value.str_field("schema") != Some(super::STUDY_RUN_RECEIPT_SCHEMA)
            || value.str_field("driver") != Some(DRIVER)
            || n > fs_project::uncertainty::MAX_SAMPLES
            || op.session.as_deref() != Some(model.as_bytes().as_slice())
            || op.ir != expected_ir
            || op.outcome.as_deref() != Some("ok")
            || !ledger.edge_exists(op_id, &hash, EdgeRole::Out)?
        {
            return Err(fail(
                "cli-uncertainty-legacy",
                "legacy receipt differs from its completed sealed producer",
            ));
        }
        let status = value
            .str_field("status")
            .ok_or_else(|| fail("cli-uncertainty-legacy", "legacy receipt has no status"))?;
        let exit_code = match status {
            "complete" => exit::SUCCESS,
            "cancelled" => exit::CANCELLED,
            "refused" => exit::REFUSED,
            "budget-truncated" => exit::BUDGET,
            _ => {
                return Err(fail(
                    "cli-uncertainty-legacy",
                    "unknown legacy receipt status",
                ));
            }
        };
        let receipt = std::str::from_utf8(&bytes)
            .map_err(|error| fail("cli-uncertainty-legacy", error.to_string()))?;
        // Both legacy verbs return the original receipt. Do not invent an
        // HTML report or evidence package the old producer never retained.
        let stdout = match mode {
            OutputMode::Json => format!(
                "{{\"command\":{},\"status\":{},\"run\":{},\"run_id\":{},\"receipt\":{receipt}}}\n",
                quoted(command),
                quoted(status),
                quoted(pointer),
                quoted(pointer)
            ),
            OutputMode::Text => format!(
                "command={command}\nstatus={status}\nrun={pointer}\nauthority=Estimated\nreceipt={receipt}\n"
            ),
        };
        Ok(CommandOutput {
            exit_code,
            stdout,
            stderr: String::new(),
        })
    })();
    match result {
        Ok(output) => output,
        Err(error) => output_error(command, mode, error),
    }
}
