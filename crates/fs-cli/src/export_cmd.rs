//! Read-only mesh and tabular artifact export verbs (`bead frankensim-vm286`):
//! - `frankensim [--json] mesh export <promoted-mesh-hash> <ledger.db> --format glb|3mf|vtk --out <path>`
//! - `frankensim [--json] table export <solution-hash> <ledger.db> --format csv|arrow --out <path>`
//!
//! Both verbs open existing ledgers through [`crate::report::open_export_ledger`]
//! (never creating a database file on a missing ledger) and write output through
//! [`crate::report::write_retained`] (idempotent on identical bytes, refusing on
//! any differing pre-existing file). Neither verb mutates any ledger transaction.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use fs_blake3::{ContentHash, hash_bytes};
use fs_io::tabular_export::{
    TabularColumn, TabularDataset, export_arrow_ipc_stream, export_tabular_csv,
};
use fs_io::{export_3mf, export_glb, export_vtk};

use crate::json_read::JsonValue;
use crate::report::{open_export_ledger, refuse, write_retained};
use crate::{
    CommandOutput, Diagnostic, OutputMode, RESULT_SCHEMA, escape_text, exit, is_operand,
    push_json_string,
};

const MESH_ARTIFACT_KIND: &str = "geometry-mesh-ply";
const SOLUTION_ARTIFACT_KIND: &str = "solve-conduction-solution";
const SOLUTION_SCHEMA: &str = "frankensim.cli.solve-conduction-solution.v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MeshExportFormat {
    Glb,
    ThreeMf,
    Vtk,
}

impl MeshExportFormat {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "glb" => Some(Self::Glb),
            "3mf" => Some(Self::ThreeMf),
            "vtk" => Some(Self::Vtk),
            _ => None,
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::Glb => "glb",
            Self::ThreeMf => "3mf",
            Self::Vtk => "vtk",
        }
    }

    const fn precision(self) -> &'static str {
        match self {
            Self::Glb => "f32-preview",
            Self::ThreeMf | Self::Vtk => "f64",
        }
    }

    const fn no_claim(self) -> &'static str {
        match self {
            Self::Glb => {
                "GLB stores IEEE-754 f32 preview coordinates; the lossless canonical f64 PLY remains in the ledger"
            }
            Self::ThreeMf | Self::Vtk => {
                "read-only export projection of the retained canonical geometry-mesh-ply artifact"
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TableExportFormat {
    Csv,
    Arrow,
}

impl TableExportFormat {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "csv" => Some(Self::Csv),
            "arrow" => Some(Self::Arrow),
            _ => None,
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::Csv => "csv",
            Self::Arrow => "arrow",
        }
    }

    const fn no_claim(self) -> &'static str {
        match self {
            Self::Csv => {
                "fs_io::export_tabular_csv formats Float64 columns at fixed {:.10e} scientific precision (~1e-10 rel), not bit-exact IEEE-754 round-trip; use --format arrow or the retained solve-conduction-solution artifact for bit-exact f64 values"
            }
            Self::Arrow => {
                "fs_io::export_arrow_ipc_stream writes fs-io ARROW1 binary columnar buffers carrying bit-exact little-endian IEEE-754 f64 values from the retained solve-conduction-solution artifact"
            }
        }
    }
}

fn push_kv(out: &mut String, key: &str, val: &str) {
    let _ = write!(out, ",\"{key}\":");
    push_json_string(out, val);
}

fn parse_format_and_out<F>(
    command: &'static str,
    code: &'static str,
    usage_hint: &'static str,
    args: &[String],
    parse_fmt: impl Fn(&str) -> Option<F>,
) -> Result<(F, PathBuf), Box<Diagnostic>> {
    let err = |msg: String| Box::new(Diagnostic::new(command, code, msg, usage_hint));
    let (mut format, mut out, mut i) = (None, None, 0);
    while i < args.len() {
        match args[i].as_str() {
            "--format" if format.is_none() => {
                let Some(raw) = args.get(i + 1).filter(|v| is_operand(v)) else {
                    return Err(err("`--format` requires a format name".into()));
                };
                let Some(parsed) = parse_fmt(raw) else {
                    return Err(err(format!("unknown export format `{raw}`")));
                };
                format = Some(parsed);
            }
            "--out" if out.is_none() => {
                let Some(path) = args.get(i + 1).filter(|v| is_operand(v)) else {
                    return Err(err("`--out` requires a destination path".into()));
                };
                out = Some(PathBuf::from(path));
            }
            other => return Err(err(format!("unexpected or repeated flag `{other}`"))),
        }
        i += 2;
    }
    let Some(format) = format else {
        return Err(err("missing required `--format` option".into()));
    };
    let Some(out) = out else {
        return Err(err("missing required `--out` option".into()));
    };
    Ok((format, out))
}

/// Parse `--format glb|3mf|vtk --out <path>` for `mesh export`.
pub(crate) fn parse_mesh_export_args(
    args: &[String],
) -> Result<(MeshExportFormat, PathBuf), Box<Diagnostic>> {
    parse_format_and_out(
        "mesh-export",
        "cli-mesh-export-usage",
        "use `frankensim [--json] mesh export <promoted-mesh-hash> <ledger.db> --format glb|3mf|vtk --out <path>`",
        args,
        MeshExportFormat::parse,
    )
}

/// Parse `--format csv|arrow --out <path>` for `table export`.
pub(crate) fn parse_table_export_args(
    args: &[String],
) -> Result<(TableExportFormat, PathBuf), Box<Diagnostic>> {
    parse_format_and_out(
        "table-export",
        "cli-table-export-usage",
        "use `frankensim [--json] table export <solution-hash> <ledger.db> --format csv|arrow --out <path>`",
        args,
        TableExportFormat::parse,
    )
}

fn load_typed_artifact(
    command: &'static str,
    artifact_hex: &str,
    ledger_path: &Path,
    expected_kind: &str,
    mode: OutputMode,
) -> Result<Vec<u8>, Box<CommandOutput>> {
    let fail = |code: &str, msg: String, fix: String| {
        Box::new(refuse(
            mode,
            command,
            exit::REFUSED,
            code,
            artifact_hex,
            msg,
            fix,
        ))
    };
    let Some(hash) = ContentHash::from_hex(artifact_hex) else {
        return Err(Box::new(refuse(
            mode,
            command,
            exit::USAGE,
            "cli-export-artifact-hash",
            artifact_hex,
            format!("invalid artifact hash `{artifact_hex}`"),
            "pass a 64-character lowercase hex BLAKE3 artifact hash",
        )));
    };
    let (ledger, _) =
        open_export_ledger(command, artifact_hex, Some(ledger_path), mode).map_err(Box::new)?;
    let info = match ledger.artifact_info(&hash) {
        Ok(Some(info)) => info,
        Ok(None) => {
            return Err(fail(
                "cli-export-unknown-artifact",
                format!("artifact `{artifact_hex}` is not present in this ledger"),
                "pass an artifact hash retained in the ledger".into(),
            ));
        }
        Err(err) => {
            return Err(fail(
                "cli-export-ledger-read",
                format!("cannot inspect artifact `{artifact_hex}`: {err}"),
                "verify ledger integrity before exporting".into(),
            ));
        }
    };
    if info.kind != expected_kind {
        return Err(fail(
            "cli-export-artifact-kind",
            format!(
                "artifact `{artifact_hex}` has kind `{}`, expected `{expected_kind}`",
                info.kind
            ),
            format!("pass a retained `{expected_kind}` artifact hash"),
        ));
    }
    match ledger.get_artifact(&hash) {
        Ok(Some(bytes)) => Ok(bytes),
        Ok(None) => Err(fail(
            "cli-export-unknown-artifact",
            format!("artifact `{artifact_hex}` bytes are missing from the ledger"),
            "verify ledger integrity before exporting".into(),
        )),
        Err(err) => Err(fail(
            "cli-export-ledger-read",
            format!("cannot read artifact `{artifact_hex}`: {err}"),
            "verify ledger integrity before exporting".into(),
        )),
    }
}

/// Execute `frankensim [--json] mesh export <promoted-mesh-hash> <ledger.db> --format glb|3mf|vtk --out <path>`.
pub(crate) fn mesh_export_path(
    artifact_hex: &str,
    ledger_path: &Path,
    format: MeshExportFormat,
    out_path: &Path,
    mode: OutputMode,
) -> CommandOutput {
    let ply_bytes = match load_typed_artifact(
        "mesh-export",
        artifact_hex,
        ledger_path,
        MESH_ARTIFACT_KIND,
        mode,
    ) {
        Ok(b) => b,
        Err(out) => return *out,
    };
    let soup = match fs_io::ply::read_ply(&ply_bytes) {
        Ok(s) => s,
        Err(err) => {
            return refuse(
                mode,
                "mesh-export",
                exit::REFUSED,
                "cli-mesh-export-decode",
                artifact_hex,
                format!("cannot decode retained PLY artifact `{artifact_hex}`: {err}"),
                "verify the promoted geometry-mesh-ply artifact in the ledger",
            );
        }
    };
    let bytes = match format {
        MeshExportFormat::Glb => export_glb(&soup),
        MeshExportFormat::ThreeMf => export_3mf(&soup),
        MeshExportFormat::Vtk => export_vtk(&soup, None).into_bytes(),
    };
    if let Err(why) = write_retained(out_path, &bytes) {
        return refuse(
            mode,
            "mesh-export",
            exit::REFUSED,
            "cli-export-output-conflict",
            artifact_hex,
            why,
            "choose a new `--out` path or remove the differing file; exports never overwrite a differing artifact",
        );
    }
    let content_hash = hash_bytes(&bytes).to_hex();
    let out_str = out_path.to_string_lossy();
    let vertices = soup.positions.len();
    let triangles = soup.triangles.len();
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", "mesh-export");
            push_kv(&mut out, "status", "ok");
            push_kv(&mut out, "subject", artifact_hex);
            push_kv(&mut out, "artifact", artifact_hex);
            push_kv(&mut out, "format", format.as_str());
            push_kv(&mut out, "precision", format.precision());
            push_kv(&mut out, "out", &out_str);
            let _ = write!(
                out,
                ",\"bytes\":{},\"vertices\":{vertices},\"triangles\":{triangles}",
                bytes.len()
            );
            push_kv(&mut out, "content_hash", &content_hash);
            push_kv(&mut out, "authority", "geometry-export-projection");
            push_kv(&mut out, "no_claim", format.no_claim());
            out.push_str("}\n");
            out
        }
        OutputMode::Text => format!(
            "status=ok\ncommand=mesh-export\nsubject={artifact_hex}\nartifact={artifact_hex}\nformat={}\nprecision={}\nout={}\nbytes={}\nvertices={vertices}\ntriangles={triangles}\ncontent_hash={content_hash}\nauthority=geometry-export-projection\nno_claim={}\n",
            format.as_str(),
            format.precision(),
            escape_text(&out_str),
            bytes.len(),
            escape_text(format.no_claim()),
        ),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

fn parse_solution_dataset(bytes: &[u8]) -> Result<TabularDataset, String> {
    let text = std::str::from_utf8(bytes).map_err(|e| format!("solution is not UTF-8: {e}"))?;
    let val = JsonValue::parse(text).map_err(|e| format!("solution JSON parse failed: {e}"))?;
    if val.str_field("schema") != Some(SOLUTION_SCHEMA) {
        return Err(format!(
            "unexpected solution schema {:?}, expected `{SOLUTION_SCHEMA}`",
            val.str_field("schema")
        ));
    }
    let unit = val.str_field("temperature_unit").unwrap_or("K");
    let raw_temps = val
        .get("temperature")
        .and_then(JsonValue::as_array)
        .ok_or_else(|| "solution JSON missing numeric `temperature` array".to_string())?;
    let mut indices = Vec::with_capacity(raw_temps.len());
    let mut temps = Vec::with_capacity(raw_temps.len());
    for (idx, item) in raw_temps.iter().enumerate() {
        let Some(t) = item.as_f64().filter(|v| v.is_finite()) else {
            return Err(format!(
                "non-finite or non-numeric temperature at index {idx}"
            ));
        };
        let node_idx = i64::try_from(idx).map_err(|_| "node index overflow".to_string())?;
        indices.push(node_idx);
        temps.push(t);
    }
    TabularDataset::new(
        "nodal_temperature",
        vec![
            TabularColumn::new_i64("node_index", "1", indices),
            TabularColumn::new_f64("temperature", unit, temps),
        ],
    )
}

/// Execute `frankensim [--json] table export <solution-hash> <ledger.db> --format csv|arrow --out <path>`.
pub(crate) fn table_export_path(
    artifact_hex: &str,
    ledger_path: &Path,
    format: TableExportFormat,
    out_path: &Path,
    mode: OutputMode,
) -> CommandOutput {
    let sol_bytes = match load_typed_artifact(
        "table-export",
        artifact_hex,
        ledger_path,
        SOLUTION_ARTIFACT_KIND,
        mode,
    ) {
        Ok(b) => b,
        Err(out) => return *out,
    };
    let dataset = match parse_solution_dataset(&sol_bytes) {
        Ok(ds) => ds,
        Err(msg) => {
            return refuse(
                mode,
                "table-export",
                exit::REFUSED,
                "cli-table-export-decode",
                artifact_hex,
                format!("cannot decode retained solution `{artifact_hex}`: {msg}"),
                "pass a valid retained `solve-conduction-solution` artifact",
            );
        }
    };
    let (bytes, receipt) = match format {
        TableExportFormat::Csv => {
            let (csv, receipt) = export_tabular_csv(&dataset);
            (csv.into_bytes(), receipt)
        }
        TableExportFormat::Arrow => export_arrow_ipc_stream(&dataset),
    };
    if let Err(why) = write_retained(out_path, &bytes) {
        return refuse(
            mode,
            "table-export",
            exit::REFUSED,
            "cli-export-output-conflict",
            artifact_hex,
            why,
            "choose a new `--out` path or remove the differing file; exports never overwrite a differing artifact",
        );
    }
    let content_hash = hash_bytes(&bytes).to_hex();
    let receipt_hash = receipt.content_hash.to_hex();
    let out_str = out_path.to_string_lossy();
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", "table-export");
            push_kv(&mut out, "status", "ok");
            push_kv(&mut out, "subject", artifact_hex);
            push_kv(&mut out, "artifact", artifact_hex);
            push_kv(&mut out, "format", format.as_str());
            push_kv(&mut out, "out", &out_str);
            push_kv(&mut out, "table_name", &receipt.table_name);
            let _ = write!(
                out,
                ",\"bytes\":{},\"column_count\":{},\"row_count\":{}",
                bytes.len(),
                receipt.column_count,
                receipt.row_count
            );
            push_kv(&mut out, "content_hash", &content_hash);
            push_kv(&mut out, "receipt_hash", &receipt_hash);
            push_kv(&mut out, "authority", "tabular-export-projection");
            push_kv(&mut out, "no_claim", format.no_claim());
            out.push_str("}\n");
            out
        }
        OutputMode::Text => format!(
            "status=ok\ncommand=table-export\nsubject={artifact_hex}\nartifact={artifact_hex}\nformat={}\nout={}\ntable_name={}\nbytes={}\ncolumn_count={}\nrow_count={}\ncontent_hash={content_hash}\nreceipt_hash={receipt_hash}\nauthority=tabular-export-projection\nno_claim={}\n",
            format.as_str(),
            escape_text(&out_str),
            escape_text(&receipt.table_name),
            bytes.len(),
            receipt.column_count,
            receipt.row_count,
            escape_text(format.no_claim()),
        ),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}
