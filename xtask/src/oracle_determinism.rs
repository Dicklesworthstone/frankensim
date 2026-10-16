//! Cross-environment oracle determinism matrix (bead `frankensim-uf7cw`),
//! companion gate to [`crate::constellation_drift`] in the constellation
//! family.
//!
//! Why this gate exists. The b2can incident (2026-08-26) showed fs-query
//! bore golden outcomes flipping GREEN->RED across runs of identical
//! committed source while sibling constellations drifted between agent
//! syncs — and the ablation receipts proved the attempted heuristic patch
//! was causally inert while its author still saw green elsewhere. A repo
//! whose determinism law is "bit-stable across runs on the same ISA"
//! cannot tolerate result classes that depend on invisible
//! build-environment composition. PASS/FAIL text cannot diagnose that;
//! bit-level numeric receipts can.
//!
//! What it charges. Every leg captures one machine-readable payload set:
//! per-case pass verdicts plus order-sensitive digests over station areas,
//! equivalent radii, arc lengths, the raw medial-pole cloud, and selected
//! scalar bit patterns — inside `git archive` frozen trees where ONLY the
//! sibling state varies. Red violations:
//!
//! - two legs of the SAME architecture class disagreeing on any verdict,
//!   count, closure class, or digest;
//! - cross-architecture legs disagreeing on verdict or structure while
//!   their declared inputs agree (bit textures may differ across ISA by
//!   law; flips may not);
//! - a leg script failure, probe refusal record, or unparseable payload;
//! - a golden run red INSIDE a frozen tree (the flip reproduced under
//!   controlled inputs — labeled as such, never as environment noise);
//! - matrix incompleteness for any lane left unproduced (fail-closed: an
//!   absent witness never counts as agreement).
//!
//! What it deliberately does not charge. Cross-ISA numeric-texture drift
//! under agreeing verdicts is rendered as an ok note — expected freedom,
//! not drift. The instrument measures COMMITTED source only, so it is
//! structurally blind to uncommitted shared-tree churn by construction
//! (the flip class came from sibling states, which it fully fingerprints).
//!
//! Leg production lives in `scripts/ci/oracle_determinism_leg.sh`: local
//! lanes are driven directly by this check; remote x86 lanes run through
//! RCH job mode and are ingested here via FS_ORACLE_REMOTE_DIRS
//! ("label=path,..."). Remote job template (one call per worker+mode):
//!   RCH_WORKER=<worker> rch exec --job --result-dir <stamp-<label>> \
//!     -- env MODE=pinned OUT_DIR=<stamp-<label>> scripts/ci/oracle_determinism_leg.sh
//! Gauntlet tier: G5 determinism audit.

use std::path::{Path, PathBuf};

use crate::PolicyNote;

pub(crate) const CHECK: &str = "oracle-determinism";

/// Environment variable listing pre-captured remote leg directories as
/// `label=path` pairs (comma-separated). Labels must be unique across the
/// whole matrix; collisions are charged, never silently deduplicated.
pub(crate) const REMOTE_DIRS_ENV: &str = "FS_ORACLE_REMOTE_DIRS";

struct LegPayload {
    label: String,
    mode: String,
    arch: String,
    goldens_pass: bool,
    cases: Vec<CasePayload>,
    refusals: Vec<String>,
}

struct CasePayload {
    case: String,
    pass: bool,
    stations: usize,
    boundary_digest: String,
    areas_digest: String,
    radius_digest: String,
    arc_digest: String,
    pole_count: usize,
    poles_digest: String,
    closure: String,
}

#[derive(Debug)]
pub(crate) struct Report {
    pub(crate) violations: Vec<crate::Violation>,
    pub(crate) notes: Vec<PolicyNote>,
}

impl Report {
    fn new() -> Self {
        Self {
            violations: Vec::new(),
            notes: Vec::new(),
        }
    }
    fn violation(&mut self, detail: impl Into<String>) {
        self.violations.push(crate::Violation {
            check: CHECK,
            crate_name: "xtask".to_string(),
            detail: detail.into(),
        });
    }
    fn note(&mut self, verdict: &'static str, detail: impl Into<String>) {
        self.notes.push(PolicyNote {
            check: CHECK,
            crate_name: "xtask".to_string(),
            verdict,
            detail: detail.into(),
        });
    }
}

/// One flat top-level string/bool/number field out of our strict probe
/// JSONL lines. Values contain no escapes by construction (hex strings,
/// fixed enum words, integers); anything else fails extraction and is
/// reported rather than guessed at.
fn flat_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\":");
    let start = line.find(&needle)? + needle.len();
    let rest = &line[start..];
    if rest.starts_with('"') {
        let end = rest[1..].find('"')? + 1;
        Some(&rest[1..end])
    } else {
        let end = rest
            .find(|c: char| c == ',' || c == '}')
            .unwrap_or(rest.len());
        Some(rest[..end].trim())
    }
}

fn bool_field(line: &str, key: &str) -> Option<bool> {
    match flat_field(line, key)? {
        "true" => Some(true),
        "false" => Some(false),
        _ => None,
    }
}

fn usize_field(line: &str, key: &str) -> Option<usize> {
    flat_field(line, key)?.parse().ok()
}

fn parse_case_line(line: &str) -> Option<CasePayload> {
    Some(CasePayload {
        case: flat_field(line, "case")?.to_string(),
        pass: bool_field(line, "pass")?,
        stations: usize_field(line, "stations")?,
        boundary_digest: flat_field(line, "boundary_digest")?.to_string(),
        areas_digest: flat_field(line, "areas_digest")?.to_string(),
        radius_digest: flat_field(line, "radius_digest")?.to_string(),
        arc_digest: flat_field(line, "arc_digest")?.to_string(),
        pole_count: usize_field(line, "pole_count")?,
        poles_digest: flat_field(line, "poles_digest")?.to_string(),
        closure: flat_field(line, "closure")?.to_string(),
    })
}

fn read_nonempty(path: &Path) -> Result<String, String> {
    std::fs::read_to_string(path)
        .map_err(|error| format!("{} unreadable: {error}", path.display()))
        .and_then(|text| {
            if text.trim().is_empty() {
                Err(format!("{} is empty", path.display()))
            } else {
                Ok(text)
            }
        })
}

fn ingest_leg(label: &str, dir: &Path) -> Result<LegPayload, String> {
    let payload_text = read_nonempty(&dir.join("probe.jsonl"))?;
    let receipt_text = read_nonempty(&dir.join("leg_receipt.json"))?;
    let receipt_flat = receipt_text.replace(['\n', '\r'], "");
    let mut cases = Vec::new();
    let mut refusals = Vec::new();
    for line in payload_text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        match flat_field(line, "record") {
            Some("header") | None => {}
            Some("case") => match parse_case_line(line) {
                Some(payload) => cases.push(payload),
                None => refusals.push(format!("unparseable case line: {line}")),
            },
            Some("error") => refusals.push(format!(
                "probe error (case {}): {}",
                flat_field(line, "case").unwrap_or("?"),
                flat_field(line, "reason").unwrap_or("?")
            )),
            Some(other) => refusals.push(format!("unknown record kind {other}")),
        }
    }
    if cases.is_empty() && refusals.is_empty() {
        return Err(format!("{label}: no case records in probe.jsonl"));
    }
    Ok(LegPayload {
        label: label.to_string(),
        mode: flat_field(&receipt_flat, "mode").unwrap_or("?").to_string(),
        arch: flat_field(&receipt_flat, "arch").unwrap_or("?").to_string(),
        goldens_pass: bool_field(&receipt_flat, "goldens_pass").unwrap_or(false),
        cases,
        refusals,
    })
}

fn same_arch_compare(report: &mut Report, mode: &str, left: &LegPayload, right: &LegPayload) {
    for lcase in &left.cases {
        let Some(rcase) = right.cases.iter().find(|c| c.case == lcase.case) else {
            report.violation(format!(
                "[{mode}] case {} present in {} but missing from {}",
                lcase.case, left.label, right.label
            ));
            continue;
        };
        if lcase.pass != rcase.pass {
            report.violation(format!(
                "[{mode}] verdict FLIP on case {}: {}={} vs {}={}",
                lcase.case, left.label, lcase.pass, right.label, rcase.pass
            ));
        }
        let mut mismatches = Vec::new();
        let digests = [
            (
                "boundary_digest",
                (&lcase.boundary_digest, &rcase.boundary_digest),
            ),
            ("areas_digest", (&lcase.areas_digest, &rcase.areas_digest)),
            (
                "radius_digest",
                (&lcase.radius_digest, &rcase.radius_digest),
            ),
            ("arc_digest", (&lcase.arc_digest, &rcase.arc_digest)),
            ("poles_digest", (&lcase.poles_digest, &rcase.poles_digest)),
        ];
        for (field, (mine, theirs)) in digests {
            if mine != theirs {
                mismatches.push(format!("{field}: {} vs {}", mine, theirs));
            }
        }
        if lcase.stations != rcase.stations {
            mismatches.push(format!(
                "stations: {} vs {}",
                lcase.stations, rcase.stations
            ));
        }
        if lcase.pole_count != rcase.pole_count {
            mismatches.push(format!(
                "pole_count: {} vs {}",
                lcase.pole_count, rcase.pole_count
            ));
        }
        if lcase.closure != rcase.closure {
            mismatches.push(format!("closure: {} vs {}", lcase.closure, rcase.closure));
        }
        if !mismatches.is_empty() {
            report.violation(format!(
                "[{mode}] numeric determinism divergence on case {} between same-arch \
                 legs `{}` and `{}`: {}; committed source is identical, so the culprit \
                 is dependency/environment composition — bisect sibling checkouts \
                 against constellation.lock pins next",
                lcase.case,
                left.label,
                right.label,
                mismatches.join("; ")
            ));
        } else {
            report.note(
                "ok",
                format!(
                    "[{mode}] case {} bit-stable across legs `{}` == `{}`",
                    lcase.case, left.label, right.label
                ),
            );
        }
    }
    for rcase in &right.cases {
        if !left.cases.iter().any(|c| c.case == rcase.case) {
            report.violation(format!(
                "[{mode}] case {} present in {} but missing from {}",
                rcase.case, right.label, left.label
            ));
        }
    }
}

fn cross_arch_compare(report: &mut Report, mode: &str, left: &LegPayload, right: &LegPayload) {
    let names: Vec<_> = left.cases.iter().map(|c| c.case.as_str()).collect();
    for rcase in &right.cases {
        let Some(lcase) = left.cases.iter().find(|c| c.case == rcase.case) else {
            report.violation(format!(
                "[{mode}] case {} present in {} but missing from {}",
                rcase.case, right.label, left.label
            ));
            continue;
        };
        if lcase.pass != rcase.pass
            || lcase.closure != rcase.closure
            || lcase.stations != rcase.stations
        {
            report.violation(format!(
                "[{mode}] structural divergence across architectures on case {}: \
                 {}/{}=pass {},{} verdict/closure/stations {}/{}/{} vs {}/{}/{}",
                rcase.case,
                left.label,
                right.label,
                lcase.pass,
                rcase.pass,
                lcase.pass,
                lcase.closure,
                lcase.stations,
                rcase.pass,
                rcase.closure,
                rcase.stations
            ));
        } else {
            report.note(
                "ok",
                format!(
                    "[{mode}] case {} agrees across architectures {} <-> {} \
                     (bit textures lawfully free to differ)",
                    rcase.case, left.label, right.label
                ),
            );
        }
    }
    let right_names: Vec<_> = right.cases.iter().map(|c| c.case.as_str()).collect();
    for name in &names {
        if !right_names.contains(name) {
            report.violation(format!(
                "[{mode}] case {} present in {} but missing from {}",
                name, left.label, right.label
            ));
        }
    }
}

fn compare_pair(report: &mut Report, mode: &str, left: &LegPayload, right: &LegPayload) {
    if left.arch == right.arch {
        same_arch_compare(report, mode, left, right);
    } else {
        cross_arch_compare(report, mode, left, right);
    }
}

/// The standing gate: produce both local legs fresh, ingest any remote
/// lanes named in the environment, compare within modes, and refuse to
/// bless anything when evidence is missing.
pub(crate) fn check(root: &Path) -> Report {
    let mut report = Report::new();
    // Frozen-tree builds are full dep-graph compiles; park every byte on
    // RCH_TARGET_BASE (external NVMe class storage) exactly like the rch
    // lanes do, never on the internal boot volume under disk pressure.
    let base = match std::env::var("RCH_TARGET_BASE") {
        Ok(dir) if !dir.trim().is_empty() => PathBuf::from(dir).join("frankensim-oracle-det"),
        _ => root.join("target/oracle-determinism"),
    };
    if std::fs::create_dir_all(&base).is_err() {
        report.violation(format!("cannot create scratch base {}", base.display()));
        return report;
    }

    // Local lanes always regenerate into run-stamped directories so stale
    // payloads can never masquerade as fresh evidence.
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or_default();
    let mut legs: Vec<LegPayload> = Vec::new();
    for mode in ["pinned", "live"] {
        let dir = base.join(format!("{stamp}-{mode}"));
        let script = root.join("scripts/ci/oracle_determinism_leg.sh");
        let status = std::process::Command::new(&script)
            .env("MODE", mode)
            .env("OUT_DIR", &dir)
            .current_dir(root)
            .status();
        match status {
            Ok(status) if status.success() => match ingest_leg(mode, &dir) {
                Ok(payload) => legs.push(payload),
                Err(error) => report.violation(error),
            },
            Ok(status) => report.violation(format!(
                "local {mode} leg exited {status}; receipts retained in the leg \
                 directory; rerun MODE={mode} OUT_DIR={} scripts/ci/oracle_determinism_leg.sh \
                 directly",
                dir.display()
            )),
            Err(error) => report.violation(format!("local {mode} leg spawn failed: {error}")),
        }
    }

    if let Ok(remote_spec) = std::env::var(REMOTE_DIRS_ENV) {
        for entry in remote_spec.split(',') {
            let entry = entry.trim();
            if entry.is_empty() {
                continue;
            }
            let Some((label, path)) = entry.split_once('=') else {
                report.violation(format!(
                    "{REMOTE_DIRS_ENV} entry `{entry}` is not label=path"
                ));
                continue;
            };
            match ingest_leg(label, Path::new(path)) {
                Ok(payload) => {
                    if legs.iter().any(|existing| existing.label == payload.label) {
                        report.violation(format!("duplicate matrix label `{}`", payload.label));
                    } else {
                        legs.push(payload);
                    }
                }
                Err(error) => report.violation(format!("remote leg {label}: {error}")),
            }
        }
    }

    // Fail-closed completeness: fewer than two distinct environments (or a
    // missing lane entirely) means the matrix cannot support any claim, so
    // it charges loudly instead of blessing partial agreement.
    let mut labels: Vec<&str> = legs.iter().map(|leg| leg.label.as_str()).collect();
    labels.sort_unstable();
    labels.dedup();
    if labels.len() < 2 {
        report.violation(format!(
            "oracle-determinism matrix incomplete: {} leg(s) [{}]. Produce the remote \
             x86 lanes via RCH job mode (scripts/ci/oracle_determinism_leg.sh with \
             MODE=pinned|live on workers ovh-a and ovh-b), then point \
             {REMOTE_DIRS_ENV}=\"<label>=<result-dir>,...\" at their outputs",
            legs.len(),
            labels.join(", ")
        ));
    }

    for group in ["pinned", "live"] {
        let group_legs: Vec<&LegPayload> = legs.iter().filter(|leg| leg.mode == group).collect();
        for window in group_legs.windows(2) {
            compare_pair(&mut report, group, window[0], window[1]);
        }
    }
    for leg in &legs {
        if !leg.refusals.is_empty() {
            report.violation(format!(
                "leg `{}` reported {} refusal record(s): {}",
                leg.label,
                leg.refusals.len(),
                leg.refusals.join(" | ")
            ));
        }
        if !leg.goldens_pass {
            report.violation(format!(
                "leg `{}` ran the gb_003/gb_004 goldens RED inside its frozen tree: \
                 the b2can-class flip reproduced under controlled inputs",
                leg.label
            ));
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe_case(case: &str, pass: bool, digest: &str) -> CasePayload {
        CasePayload {
            case: case.to_string(),
            pass,
            stations: 33,
            boundary_digest: "aaaa".to_string(),
            areas_digest: digest.to_string(),
            radius_digest: digest.to_string(),
            arc_digest: digest.to_string(),
            pole_count: 42,
            poles_digest: digest.to_string(),
            closure: "Skipped".to_string(),
        }
    }

    fn leg(label: &str, mode: &str, arch: &str, cases: Vec<CasePayload>) -> LegPayload {
        LegPayload {
            label: label.to_string(),
            mode: mode.to_string(),
            arch: arch.to_string(),
            goldens_pass: true,
            cases,
            refusals: Vec::new(),
        }
    }

    #[test]
    fn same_arch_bit_divergence_is_charged_with_both_sides_named() {
        let mut report = Report::new();
        let a = leg(
            "mac-pinned-a",
            "pinned",
            "arm64",
            vec![probe_case("gb-003", true, "dead")],
        );
        let b = leg(
            "ovh-pinned-b",
            "pinned",
            "x86_64",
            vec![probe_case("gb-003", true, "c0de")],
        );
        // Same digest value but different arch => cross-arch lane; verify
        // structural agreement stays ok.
        cross_arch_compare(&mut report, "pinned", &a, &b);
        assert!(report.violations.is_empty());
        assert_eq!(report.notes.len(), 1);

        // Same arch, different bits => charged with expected/observed sides.
        let mut report2 = Report::new();
        let b2 = leg(
            "ovh2-pinned",
            "pinned",
            "arm64",
            vec![probe_case("gb-003", true, "beef")],
        );
        same_arch_compare(&mut report2, "pinned", &a, &b2);
        assert_eq!(report2.violations.len(), 1);
        assert!(report2.violations[0].detail.contains("areas_digest"));
        assert!(report2.violations[0].detail.contains("beef"));
    }

    #[test]
    fn cross_arch_verdict_flip_is_charged_but_texture_is_not() {
        let mut report = Report::new();
        let a = leg(
            "arm-a",
            "live",
            "arm64",
            vec![probe_case("gb-004", true, "1111")],
        );
        let b = leg(
            "x86-b",
            "live",
            "x86_64",
            vec![probe_case("gb-004", false, "2222")],
        );
        cross_arch_compare(&mut report, "live", &a, &b);
        assert_eq!(report.violations.len(), 1);
        assert!(
            report.violations[0]
                .detail
                .contains("structural divergence")
        );
        assert!(report.notes.is_empty());
    }

    #[test]
    fn missing_case_on_either_side_is_charged() {
        let mut report = Report::new();
        let a = leg(
            "l1",
            "pinned",
            "arm64",
            vec![probe_case("gb-003", true, "aa")],
        );
        let b = leg(
            "l2",
            "pinned",
            "arm64",
            vec![probe_case("gb-004", true, "aa")],
        );
        same_arch_compare(&mut report, "pinned", &a, &b);
        assert_eq!(report.violations.len(), 2);
        for item in &report.violations {
            assert!(item.detail.contains("missing from"), "{}", item.detail);
        }
    }

    #[test]
    fn flat_field_extracts_strings_bools_and_counts() {
        let line = "{\"record\":\"case\",\"case\":\"gb-003\",\"pass\":true,\
                    \"stations\":33,\"closure\":\"Skipped\",\"poles\":7}";
        assert_eq!(flat_field(line, "record"), Some("case"));
        assert_eq!(flat_field(line, "case"), Some("gb-003"));
        assert_eq!(bool_field(line, "pass"), Some(true));
        assert_eq!(usize_field(line, "stations"), Some(33));
        assert_eq!(flat_field(line, "missing"), None);
    }

    #[test]
    fn probe_error_lines_become_refusal_records() {
        let dir = std::env::temp_dir().join(format!("fs-oracle-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("mkdir");
        std::fs::write(dir.join("probe.jsonl"),
            "{\"record\":\"header\",\"schema\":\"x\",\"env\":\"e\",\"arch\":\"arm64\",\"os\":\"m\"}\n\
             {\"record\":\"error\",\"case\":\"gb-003\",\"reason\":\"medial_poles refused\"}\n")
            .expect("write payload");
        std::fs::write(dir.join("leg_receipt.json"),
            "{\"record\":\"leg_receipt\",\"schema\":\"frankensim-oracle-leg-v1\",\"mode\":\"pinned\",\n\
             \"arch\":\"arm64\",\"goldens_pass\":true}")
            .expect("write receipt");
        let leg = ingest_leg("test-leg", &dir).expect("ingest");
        assert!(leg.goldens_pass);
        assert_eq!(leg.refusals.len(), 1);
        assert!(leg.refusals[0].contains("medial_poles refused"));
        // Scratch left for the OS temp sweeper; sessions never delete.
    }
}
