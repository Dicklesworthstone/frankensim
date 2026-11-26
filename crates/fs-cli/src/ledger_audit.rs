//! Read-only design-ledger audit verbs (`bead frankensim-dhqf7`):
//! - `frankensim [--json] lineage <artifact-hash> <ledger.db> [--depth <N>]`
//! - `frankensim [--json] replay-audit <a.db> <b.db> [--branch-a <b1>] [--branch-b <b2>]`
//!
//! Both verbs open existing ledgers through [`crate::report::open_export_ledger`]
//! and refuse missing files with `exit::INPUT` before touching SQLite, so
//! neither verb ever creates a database file.

use std::fmt::Write as _;
use std::path::Path;

use fs_blake3::ContentHash;
use fs_ledger::{BranchDiff, BranchInfo, Ledger, MAIN_BRANCH, ReplayMismatch, ReplayVerdict};

use crate::report::{open_export_ledger, refuse};
use crate::{
    CommandOutput, Diagnostic, OutputMode, RESULT_SCHEMA, escape_text, exit, format_diagnostic,
    is_operand, push_json_string,
};

pub(crate) const DEFAULT_LINEAGE_DEPTH: usize = 64;
const MAX_LINEAGE_DEPTH: usize = 512;

fn push_kv(out: &mut String, key: &str, val: &str) {
    let _ = write!(out, ",\"{key}\":");
    push_json_string(out, val);
}

/// Parse optional `[--depth <N>]` flags for `lineage`.
pub(crate) fn parse_lineage_args(args: &[String]) -> Result<usize, Box<Diagnostic>> {
    let err = |m: String| {
        Box::new(Diagnostic::new(
            "lineage",
            "cli-lineage-usage",
            m,
            "use `frankensim [--json] lineage <hash> <ledger.db> [--depth <N>]` (1..=512)",
        ))
    };
    let (mut depth, mut i) = (None, 0);
    while i < args.len() {
        if args[i] != "--depth" || depth.is_some() {
            return Err(err(format!("unexpected or repeated flag `{}`", args[i])));
        }
        let Some(n) = args
            .get(i + 1)
            .filter(|v| is_operand(v))
            .and_then(|raw| raw.parse::<usize>().ok())
        else {
            return Err(err("`--depth` requires an integer in 1..=512".into()));
        };
        if !(1..=MAX_LINEAGE_DEPTH).contains(&n) {
            return Err(err(format!(
                "`--depth` must be in 1..={MAX_LINEAGE_DEPTH}, got {n}"
            )));
        }
        depth = Some(n);
        i += 2;
    }
    Ok(depth.unwrap_or(DEFAULT_LINEAGE_DEPTH))
}

/// Parse optional `[--branch-a <b1>] [--branch-b <b2>]` flags for `replay-audit`.
pub(crate) fn parse_replay_audit_args(
    args: &[String],
) -> Result<(Option<String>, Option<String>), Box<Diagnostic>> {
    let err = |m: String| {
        Box::new(Diagnostic::new(
            "replay-audit",
            "cli-replay-audit-usage",
            m,
            "use `frankensim [--json] replay-audit <a.db> <b.db> [--branch-a <b1>] [--branch-b <b2>]`",
        ))
    };
    let (mut a, mut b, mut i) = (None, None, 0);
    while i < args.len() {
        let slot = match args[i].as_str() {
            "--branch-a" if a.is_none() => &mut a,
            "--branch-b" if b.is_none() => &mut b,
            other => return Err(err(format!("unexpected or repeated flag `{other}`"))),
        };
        let Some(val) = args.get(i + 1).filter(|v| is_operand(v)) else {
            return Err(err(format!(
                "`{}` requires a branch name or numeric id",
                args[i]
            )));
        };
        *slot = Some(val.clone());
        i += 2;
    }
    Ok((a, b))
}

/// Execute `frankensim [--json] lineage <artifact-hash> <ledger.db> [--depth <N>]`.
pub(crate) fn lineage_path(
    artifact_hex: &str,
    ledger_path: &Path,
    max_depth: usize,
    mode: OutputMode,
) -> CommandOutput {
    let Some(hash) = ContentHash::from_hex(artifact_hex) else {
        return refuse(
            mode,
            "lineage",
            exit::USAGE,
            "cli-lineage-hash",
            artifact_hex,
            format!("invalid artifact hash `{artifact_hex}`"),
            "pass a 64-character lowercase hex BLAKE3 artifact hash",
        );
    };
    let (ledger, _) = match open_export_ledger("lineage", artifact_hex, Some(ledger_path), mode) {
        Ok(opened) => opened,
        Err(out) => return out,
    };
    let node = match ledger.explain(&hash, max_depth) {
        Ok(Some(node)) => node,
        Ok(None) => {
            return refuse(
                mode,
                "lineage",
                exit::REFUSED,
                "cli-lineage-unknown-artifact",
                artifact_hex,
                format!("artifact `{artifact_hex}` is not present in this ledger"),
                "pass an artifact hash recorded in the ledger",
            );
        }
        Err(err) => {
            return refuse(
                mode,
                "lineage",
                exit::REFUSED,
                "cli-lineage-error",
                artifact_hex,
                format!("cannot walk artifact lineage: {err}"),
                "verify ledger integrity before retrying lineage queries",
            );
        }
    };
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", "lineage");
            push_kv(&mut out, "status", "ok");
            push_kv(&mut out, "subject", artifact_hex);
            push_kv(&mut out, "artifact", &node.hash_hex);
            push_kv(&mut out, "kind", &node.kind);
            let _ = writeln!(
                out,
                ",\"max_depth\":{max_depth},\"truncated\":{},\"producer_count\":{},\"lineage\":{}}}",
                node.truncated,
                node.produced_by.len(),
                node.to_json()
            );
            out
        }
        OutputMode::Text => format!(
            "status=ok\ncommand=lineage\nsubject={}\nmax_depth={max_depth}\n{}",
            escape_text(artifact_hex),
            node.render_text()
        ),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

fn resolve_branch(
    ledger: &Ledger,
    spec: Option<&str>,
    side: &str,
) -> Result<BranchInfo, (&'static str, String)> {
    let lookup = match spec {
        None => ledger.branch(MAIN_BRANCH),
        Some(raw) => {
            if let Ok(id) = raw.parse::<i64>() {
                ledger.branch(id)
            } else {
                ledger.branch_by_name(raw)
            }
        }
    };
    match lookup {
        Ok(Some(info)) => Ok(info),
        Ok(None) => Err((
            "cli-replay-audit-branch",
            format!(
                "branch `{}` not found on ledger {side}",
                spec.unwrap_or("main")
            ),
        )),
        Err(err) => Err((
            "cli-replay-audit-error",
            format!("cannot inspect branch on ledger {side}: {err}"),
        )),
    }
}

fn push_mismatches(out: &mut String, key: &str, items: &[ReplayMismatch]) {
    let push_strs = |dst: &mut String, xs: &[String]| {
        dst.push('[');
        for (j, h) in xs.iter().enumerate() {
            if j > 0 {
                dst.push(',');
            }
            push_json_string(dst, h);
        }
        dst.push(']');
    };
    let _ = write!(out, ",\"{key}\":[");
    for (i, m) in items.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"position\":{},\"op_a\":{},\"op_b\":{},\"only_a\":",
            m.position, m.op_a, m.op_b
        );
        push_strs(out, &m.only_a);
        out.push_str(",\"only_b\":");
        push_strs(out, &m.only_b);
        out.push('}');
    }
    out.push(']');
}

fn push_i64s(out: &mut String, items: &[i64]) {
    out.push('[');
    for (i, v) in items.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "{v}");
    }
    out.push(']');
}

fn render_replay_audit(
    subject: &str,
    ba: &BranchInfo,
    bb: &BranchInfo,
    verdict: &ReplayVerdict,
    diff: Option<&BranchDiff>,
    mode: OutputMode,
) -> CommandOutput {
    let clean = verdict.is_replay_clean();
    let (status, exit_code) = if clean {
        ("clean", exit::SUCCESS)
    } else {
        ("diverged", exit::REFUSED)
    };
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", "replay-audit");
            push_kv(&mut out, "status", status);
            push_kv(&mut out, "subject", subject);
            let _ = write!(out, ",\"branch_a_id\":{}", ba.id);
            push_kv(&mut out, "branch_a", &ba.name);
            let _ = write!(out, ",\"branch_b_id\":{}", bb.id);
            push_kv(&mut out, "branch_b", &bb.name);
            let _ = write!(out, ",\"clean\":{clean},\"compared\":{}", verdict.compared);
            match &verdict.structure_mismatch {
                Some(msg) => push_kv(&mut out, "structure_mismatch", msg),
                None => out.push_str(",\"structure_mismatch\":null"),
            }
            push_mismatches(
                &mut out,
                "deterministic_mismatches",
                &verdict.deterministic_mismatches,
            );
            push_mismatches(&mut out, "fast_divergences", &verdict.fast_divergences);
            if let Some(d) = diff {
                let _ = write!(
                    out,
                    ",\"branch_diff\":{{\"shared\":{},\"only_a\":",
                    d.shared
                );
                push_i64s(&mut out, &d.only_a);
                out.push_str(",\"only_b\":");
                push_i64s(&mut out, &d.only_b);
                out.push('}');
            }
            out.push_str("}\n");
            out
        }
        OutputMode::Text => {
            let mut out = format!(
                "status={status}\ncommand=replay-audit\nsubject={}\nbranch_a={} ({})\nbranch_b={} ({})\nclean={clean}\ncompared={}\ndeterministic_mismatches={}\nfast_divergences={}\n",
                escape_text(subject),
                escape_text(&ba.name),
                ba.id,
                escape_text(&bb.name),
                bb.id,
                verdict.compared,
                verdict.deterministic_mismatches.len(),
                verdict.fast_divergences.len()
            );
            if let Some(msg) = &verdict.structure_mismatch {
                let _ = writeln!(out, "structure_mismatch={}", escape_text(msg));
            }
            if let Some(d) = diff {
                let _ = writeln!(
                    out,
                    "branch_diff_shared={}\nbranch_diff_only_a={}\nbranch_diff_only_b={}",
                    d.shared,
                    d.only_a.len(),
                    d.only_b.len()
                );
            }
            out
        }
    };
    let stderr = if clean {
        String::new()
    } else {
        let detail = verdict.structure_mismatch.clone().unwrap_or_else(|| {
            format!(
                "{} deterministic output mismatch(es) across {} compared op(s)",
                verdict.deterministic_mismatches.len(),
                verdict.compared
            )
        });
        let diag = Diagnostic::new(
            "replay-audit",
            "cli-replay-audit-mismatch",
            detail,
            "inspect mismatched op positions and inputs; deterministic ops must reproduce identical hashes",
        )
        .with_subject(subject);
        format_diagnostic(mode, &diag)
    };
    CommandOutput {
        exit_code,
        stdout,
        stderr,
    }
}

/// Execute `frankensim [--json] replay-audit <a.db> <b.db> [--branch-a <b1>] [--branch-b <b2>]`.
pub(crate) fn replay_audit_path(
    left_path: &Path,
    right_path: &Path,
    branch_a: Option<&str>,
    branch_b: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let subject = format!("{} vs {}", left_path.display(), right_path.display());
    let (left, _) = match open_export_ledger("replay-audit", &subject, Some(left_path), mode) {
        Ok(opened) => opened,
        Err(out) => return out,
    };
    let (right, _) = match open_export_ledger("replay-audit", &subject, Some(right_path), mode) {
        Ok(opened) => opened,
        Err(out) => return out,
    };
    let ba = match resolve_branch(&left, branch_a, "a") {
        Ok(b) => b,
        Err((code, msg)) => {
            return refuse(
                mode,
                "replay-audit",
                exit::REFUSED,
                code,
                &subject,
                msg,
                "pass an existing branch name or id",
            );
        }
    };
    let bb = match resolve_branch(&right, branch_b, "b") {
        Ok(b) => b,
        Err((code, msg)) => {
            return refuse(
                mode,
                "replay-audit",
                exit::REFUSED,
                code,
                &subject,
                msg,
                "pass an existing branch name or id",
            );
        }
    };
    let verdict = match left.replay_verdict(ba.id, &right, bb.id) {
        Ok(v) => v,
        Err(err) => {
            return refuse(
                mode,
                "replay-audit",
                exit::REFUSED,
                "cli-replay-audit-error",
                &subject,
                format!("replay audit failed: {err}"),
                "verify both ledgers before re-running replay-audit",
            );
        }
    };
    let same_file = matches!(
        (std::fs::canonicalize(left_path), std::fs::canonicalize(right_path)),
        (Ok(a), Ok(b)) if a == b
    );
    let diff = if same_file && ba.id != bb.id {
        match left.branch_diff(ba.id, bb.id) {
            Ok(d) => Some(d),
            Err(err) => {
                return refuse(
                    mode,
                    "replay-audit",
                    exit::REFUSED,
                    "cli-replay-audit-error",
                    &subject,
                    format!("branch diff failed: {err}"),
                    "verify both branches exist in the ledger",
                );
            }
        }
    } else {
        None
    };
    render_replay_audit(&subject, &ba, &bb, &verdict, diff.as_ref(), mode)
}
