//! Product-spine ratchet (bead frankensim-o5et9).
//!
//! The solve pipeline executes a leading PREFIX of its stages; the rest are
//! typed gaps naming the bead that owns each. This check makes that prefix a
//! monotone, reviewed fact instead of something anyone has to go and measure.
//!
//! Why a ratchet rather than a threshold: an absolute threshold has to be
//! hand-edited on every advance, which is how gates rot. A ratchet records the
//! high-water mark and refuses BOTH directions — a retreat is a regression,
//! and an advance is a stale artifact that must be committed deliberately. The
//! spine can only move forward on the record, and it can never move without
//! someone seeing it.
//!
//! The stage table is read out of `fs-cli`'s own source rather than by
//! depending on the crate: `xtask` is policy tooling and must not pull the
//! whole product stack (and therefore FrankenSQLite) into every policy run.
//! Source-text derivation is the established idiom here — `check-powi`,
//! `check-libm`, `check-casual-print` and `check-terminology` all work this
//! way. To keep that safe, this module derives the SAME table three
//! independent ways (`ALL`, `name()`, `gap_dependency()`) and refuses unless
//! all three agree, so a source shape it no longer understands fails loudly
//! instead of silently reporting a smaller pipeline.

use std::collections::BTreeMap;
use std::path::Path;

use crate::Violation;

pub(crate) const CHECK: &str = "spine-ratchet";
pub(crate) const SOLVE_SOURCE: &str = "crates/fs-cli/src/solve.rs";
const RATCHET_PATH: &str = "spine-ratchet.json";
const SCHEMA: &str = "frankensim-spine-ratchet-v1";

/// One solve stage as derived from the product source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Stage {
    /// Rust variant name, e.g. `FlowNetwork`.
    pub variant: String,
    /// Stable kebab name, e.g. `flow-network`.
    pub name: String,
    /// Owning bead when the stage is still a typed gap.
    pub gap_owner: Option<String>,
}

fn violation(detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: RATCHET_PATH.to_string(),
        detail: detail.into(),
    }
}

/// Take the balanced `{...}` or `[...]` body that follows `marker`.
fn block_after(source: &str, marker: &str, open: char, close: char) -> Option<String> {
    let start = source.find(marker)? + marker.len();
    let rest = &source[start..];
    let open_at = rest.find(open)?;
    let bytes = rest.as_bytes();
    let mut depth = 0usize;
    for (index, byte) in bytes.iter().enumerate().skip(open_at) {
        if *byte == open as u8 {
            depth += 1;
        } else if *byte == close as u8 {
            depth -= 1;
            if depth == 0 {
                return Some(rest[open_at + 1..index].to_string());
            }
        }
    }
    None
}

/// Variant names in `SolveStage::ALL`, in execution order.
fn parse_all(source: &str) -> Result<Vec<String>, String> {
    let body = block_after(source, "const ALL: [SolveStage;", '[', ']')
        .ok_or_else(|| "cannot locate the SolveStage::ALL array".to_string())?;
    let stages: Vec<String> = body
        .split(',')
        .filter_map(|item| item.trim().strip_prefix("SolveStage::"))
        .map(|item| item.trim().to_string())
        .filter(|item| !item.is_empty())
        .collect();
    if stages.is_empty() {
        return Err("SolveStage::ALL parsed to zero stages".to_string());
    }
    Ok(stages)
}

/// `variant => "literal"` arms from a `match self` body.
fn parse_string_match(source: &str, marker: &str) -> Result<BTreeMap<String, String>, String> {
    let body =
        block_after(source, marker, '{', '}').ok_or_else(|| format!("cannot locate {marker}"))?;
    let mut map = BTreeMap::new();
    for line in body.lines() {
        let line = line.trim();
        let Some((left, right)) = line.split_once("=>") else {
            continue;
        };
        let Some(variant) = left.trim().strip_prefix("SolveStage::") else {
            continue;
        };
        let right = right.trim();
        let Some(open) = right.find('"') else {
            continue;
        };
        let Some(close) = right[open + 1..].find('"') else {
            continue;
        };
        map.insert(
            variant.trim().to_string(),
            right[open + 1..open + 1 + close].to_string(),
        );
    }
    if map.is_empty() {
        return Err(format!("{marker} parsed to zero arms"));
    }
    Ok(map)
}

/// `gap_dependency` arms: `None` for executing, `Some("bead")` for a gap.
/// The first arm may list several variants separated by `|`.
fn parse_gaps(source: &str) -> Result<BTreeMap<String, Option<String>>, String> {
    let body = block_after(source, "fn gap_dependency", '{', '}')
        .ok_or_else(|| "cannot locate fn gap_dependency".to_string())?;
    let mut map = BTreeMap::new();
    for line in body.lines() {
        let line = line.trim();
        let Some((left, right)) = line.split_once("=>") else {
            continue;
        };
        if !left.contains("SolveStage::") {
            continue;
        }
        let right = right.trim();
        let owner = if right.starts_with("None") {
            None
        } else if right.starts_with("Some") {
            let open = right
                .find('"')
                .ok_or_else(|| format!("gap arm has no bead literal: {line}"))?;
            let close = right[open + 1..]
                .find('"')
                .ok_or_else(|| format!("gap arm has an unterminated bead literal: {line}"))?;
            Some(right[open + 1..open + 1 + close].to_string())
        } else {
            continue;
        };
        for part in left.split('|') {
            if let Some(variant) = part.trim().strip_prefix("SolveStage::") {
                map.insert(variant.trim().to_string(), owner.clone());
            }
        }
    }
    if map.is_empty() {
        return Err("fn gap_dependency parsed to zero arms".to_string());
    }
    Ok(map)
}

/// Derive the stage table from product source, cross-validating three
/// independent spellings of it.
pub(crate) fn derive_stages(source: &str) -> Result<Vec<Stage>, String> {
    let order = parse_all(source)?;
    let names = parse_string_match(source, "pub const fn name(self) -> &'static str")?;
    let gaps = parse_gaps(source)?;

    if names.len() != order.len() {
        return Err(format!(
            "SolveStage::ALL lists {} stages but name() has {} arms; the source shape changed and \
             this check will not guess",
            order.len(),
            names.len()
        ));
    }
    if gaps.len() != order.len() {
        return Err(format!(
            "SolveStage::ALL lists {} stages but gap_dependency() covers {}; the source shape \
             changed and this check will not guess",
            order.len(),
            gaps.len()
        ));
    }

    let mut stages = Vec::with_capacity(order.len());
    for variant in order {
        let name = names
            .get(&variant)
            .ok_or_else(|| format!("SolveStage::{variant} has no name() arm"))?
            .clone();
        let gap_owner = gaps
            .get(&variant)
            .ok_or_else(|| format!("SolveStage::{variant} has no gap_dependency() arm"))?
            .clone();
        stages.push(Stage {
            variant,
            name,
            gap_owner,
        });
    }
    Ok(stages)
}

/// The executing prefix: stages before the first typed gap.
pub(crate) fn executing_prefix(stages: &[Stage]) -> Vec<&Stage> {
    stages
        .iter()
        .take_while(|stage| stage.gap_owner.is_none())
        .collect()
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Ratchet {
    pub executing: usize,
    pub total: usize,
    pub prefix: Vec<String>,
}

fn parse_ratchet(text: &str) -> Result<Ratchet, String> {
    if text.trim().is_empty() {
        return Err("the tracked ratchet artifact is empty".to_string());
    }
    if !text.contains(SCHEMA) {
        return Err(format!(
            "the tracked ratchet artifact does not declare schema {SCHEMA}"
        ));
    }
    let number = |key: &str| -> Result<usize, String> {
        let tag = format!("\"{key}\":");
        let start = text
            .find(&tag)
            .ok_or_else(|| format!("the tracked ratchet artifact has no `{key}` field"))?
            + tag.len();
        let rest = text[start..].trim_start();
        let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
        digits
            .parse::<usize>()
            .map_err(|_| format!("`{key}` is not a non-negative integer"))
    };
    let prefix_body = crate::spine_ratchet::block_after(text, "\"executing_prefix\"", '[', ']')
        .ok_or_else(|| {
            "the tracked ratchet artifact has no `executing_prefix` array".to_string()
        })?;
    let prefix: Vec<String> = prefix_body
        .split(',')
        .filter_map(|item| {
            let item = item.trim();
            let open = item.find('"')?;
            let close = item[open + 1..].find('"')?;
            Some(item[open + 1..open + 1 + close].to_string())
        })
        .collect();
    Ok(Ratchet {
        executing: number("high_water_stages_executing")?,
        total: number("stages_total")?,
        prefix,
    })
}

/// Render the artifact that describes the current source.
pub(crate) fn render(stages: &[Stage]) -> String {
    let prefix = executing_prefix(stages);
    let first_gap = stages.iter().find(|stage| stage.gap_owner.is_some());
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("  \"schema\": \"{SCHEMA}\",\n"));
    out.push_str("  \"bead\": \"frankensim-o5et9\",\n");
    out.push_str(&format!(
        "  \"high_water_stages_executing\": {},\n",
        prefix.len()
    ));
    out.push_str(&format!("  \"stages_total\": {},\n", stages.len()));
    out.push_str("  \"executing_prefix\": [");
    for (index, stage) in prefix.iter().enumerate() {
        if index > 0 {
            out.push_str(", ");
        }
        out.push_str(&format!("\"{}\"", stage.name));
    }
    out.push_str("],\n");
    match first_gap {
        Some(stage) => {
            out.push_str(&format!("  \"first_gap\": \"{}\",\n", stage.name));
            out.push_str(&format!(
                "  \"first_gap_owner\": \"{}\",\n",
                stage.gap_owner.as_deref().unwrap_or("")
            ));
        }
        None => {
            out.push_str("  \"first_gap\": null,\n");
            out.push_str("  \"first_gap_owner\": null,\n");
        }
    }
    out.push_str(
        "  \"no_claim\": \"counts the executing stage PREFIX declared in fs-cli source; it does \
         not prove those stages are correct, only that the product still admits them\"\n",
    );
    out.push_str("}\n");
    out
}

pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let source = std::fs::read_to_string(root.join(SOLVE_SOURCE))
        .map_err(|error| format!("cannot read {SOLVE_SOURCE}: {error}"))?;
    let stages = derive_stages(&source)?;
    std::fs::write(root.join(RATCHET_PATH), render(&stages))
        .map_err(|error| format!("cannot write {RATCHET_PATH}: {error}"))
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    let source = match std::fs::read_to_string(root.join(SOLVE_SOURCE)) {
        Ok(text) => text,
        Err(error) => return vec![violation(format!("cannot read {SOLVE_SOURCE}: {error}"))],
    };
    let stages = match derive_stages(&source) {
        Ok(stages) => stages,
        // A parser that no longer understands the source REFUSES. It must
        // never fall through to "zero gaps found", which would read as a
        // complete pipeline.
        Err(error) => {
            return vec![violation(format!(
                "cannot derive the solve stage table from {SOLVE_SOURCE}: {error}"
            ))];
        }
    };

    let tracked = match std::fs::read_to_string(root.join(RATCHET_PATH)) {
        Ok(text) => text,
        Err(error) => {
            return vec![violation(format!(
                "the ratchet artifact is missing or unreadable ({error}); a ratchet that cannot be \
                 read is not a ratchet — run `cargo run -p xtask -- generate-spine-ratchet`"
            ))];
        }
    };
    let recorded = match parse_ratchet(&tracked) {
        Ok(recorded) => recorded,
        Err(error) => {
            return vec![violation(format!(
                "{error}; refusing to pass on a ratchet this check cannot read"
            ))];
        }
    };

    let prefix: Vec<String> = executing_prefix(&stages)
        .into_iter()
        .map(|stage| stage.name.clone())
        .collect();
    let mut violations = Vec::new();

    if stages.len() != recorded.total {
        violations.push(violation(format!(
            "the solve pipeline now declares {} stages but the ratchet records {}; regenerate the \
             ratchet deliberately",
            stages.len(),
            recorded.total
        )));
    }

    if prefix.len() < recorded.executing {
        let lost: Vec<&String> = recorded
            .prefix
            .iter()
            .filter(|name| !prefix.contains(name))
            .collect();
        violations.push(violation(format!(
            "SPINE REGRESSION: {} of {} stages execute but {} did before; stage(s) {:?} stopped \
             executing. A stage that once executed may not silently become a gap again",
            prefix.len(),
            stages.len(),
            recorded.executing,
            lost
        )));
    } else if prefix.len() > recorded.executing {
        violations.push(violation(format!(
            "the spine ADVANCED: {} of {} stages now execute against a recorded high-water of {}. \
             Run `cargo run -p xtask -- generate-spine-ratchet` and commit it, so the advance is on \
             the record",
            prefix.len(),
            stages.len(),
            recorded.executing
        )));
    } else if prefix != recorded.prefix {
        // Same count, different stages. A count alone cannot see a swap.
        violations.push(violation(format!(
            "the executing prefix changed without changing length: recorded {:?}, derived {:?}",
            recorded.prefix, prefix
        )));
    }

    violations
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = r#"
impl SolveStage {
    pub const ALL: [SolveStage; 4] = [
        SolveStage::A,
        SolveStage::B,
        SolveStage::C,
        SolveStage::D,
    ];
    pub const fn name(self) -> &'static str {
        match self {
            SolveStage::A => "a",
            SolveStage::B => "b",
            SolveStage::C => "c",
            SolveStage::D => "d",
        }
    }
    pub const fn gap_dependency(self) -> Option<&'static str> {
        match self {
            SolveStage::A | SolveStage::B => None,
            SolveStage::C => Some("bead-c"),
            SolveStage::D => Some("bead-d"),
        }
    }
}
"#;

    fn ratchet_json(executing: usize, total: usize, prefix: &[&str]) -> String {
        let names: Vec<String> = prefix.iter().map(|name| format!("\"{name}\"")).collect();
        format!(
            "{{\"schema\":\"{SCHEMA}\",\"high_water_stages_executing\":{executing},\
             \"stages_total\":{total},\"executing_prefix\":[{}]}}",
            names.join(",")
        )
    }

    #[test]
    fn g0_the_stage_table_is_derived_from_three_agreeing_spellings() {
        let stages = derive_stages(FIXTURE).expect("fixture parses");
        assert_eq!(stages.len(), 4);
        assert_eq!(stages[0].name, "a");
        assert_eq!(stages[0].gap_owner, None);
        assert_eq!(stages[2].gap_owner.as_deref(), Some("bead-c"));
        // The multi-variant `A | B => None` arm must expand, not collapse.
        assert_eq!(stages[1].gap_owner, None);
        assert_eq!(
            executing_prefix(&stages)
                .into_iter()
                .map(|stage| stage.name.as_str())
                .collect::<Vec<_>>(),
            vec!["a", "b"]
        );
    }

    #[test]
    fn g0_a_source_shape_this_parser_cannot_read_refuses_rather_than_under_counting() {
        // gap_dependency covers fewer stages than ALL declares. Silently
        // treating the missing ones as executing would report a LONGER
        // pipeline than exists, which is the one failure this gate must
        // never have.
        let truncated = FIXTURE.replace("            SolveStage::D => Some(\"bead-d\"),\n", "");
        let error = derive_stages(&truncated).expect_err("must refuse");
        assert!(error.contains("gap_dependency"), "{error}");

        let no_all = FIXTURE.replace("const ALL: [SolveStage;", "const NOPE: [SolveStage;");
        assert!(derive_stages(&no_all).is_err());
    }

    #[test]
    fn g0_the_ratchet_refuses_a_regression_and_names_the_lost_stage() {
        let stages = derive_stages(FIXTURE).expect("parses");
        let prefix: Vec<String> = executing_prefix(&stages)
            .into_iter()
            .map(|stage| stage.name.clone())
            .collect();
        assert_eq!(prefix, vec!["a", "b"]);

        // Recorded high-water is 3 but only 2 execute: a real retreat.
        let recorded = parse_ratchet(&ratchet_json(3, 4, &["a", "b", "c"])).expect("parses");
        assert!(prefix.len() < recorded.executing);
        let lost: Vec<&String> = recorded
            .prefix
            .iter()
            .filter(|name| !prefix.contains(name))
            .collect();
        assert_eq!(lost, vec!["c"]);
    }

    #[test]
    fn g0_an_advance_is_also_refused_so_it_has_to_be_committed() {
        let recorded = parse_ratchet(&ratchet_json(1, 4, &["a"])).expect("parses");
        let stages = derive_stages(FIXTURE).expect("parses");
        assert!(executing_prefix(&stages).len() > recorded.executing);
    }

    #[test]
    fn g0_a_same_length_prefix_swap_is_still_caught() {
        // A count alone is blind to a permutation: two stages execute in
        // both, but not the same two.
        let recorded = parse_ratchet(&ratchet_json(2, 4, &["a", "c"])).expect("parses");
        let stages = derive_stages(FIXTURE).expect("parses");
        let prefix: Vec<String> = executing_prefix(&stages)
            .into_iter()
            .map(|stage| stage.name.clone())
            .collect();
        assert_eq!(prefix.len(), recorded.executing);
        assert_ne!(prefix, recorded.prefix);
    }

    #[test]
    fn g0_an_empty_or_foreign_ratchet_refuses_instead_of_reporting_zero_violations() {
        // The specific hazard this repo has already suffered once: an
        // unreadable gate input that reports zero violations reads exactly
        // like a clean run.
        assert!(parse_ratchet("").is_err());
        assert!(parse_ratchet("   \n ").is_err());
        assert!(parse_ratchet("{\"schema\":\"something-else\"}").is_err());
        assert!(parse_ratchet(&format!("{{\"schema\":\"{SCHEMA}\"}}")).is_err());
    }

    #[test]
    fn g0_the_live_solve_source_still_parses() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let source = std::fs::read_to_string(root.join(SOLVE_SOURCE)).expect("live source");
        let stages = derive_stages(&source).expect("the live stage table must remain derivable");
        assert_eq!(stages.len(), 6, "solve declares six stages today");
        let prefix = executing_prefix(&stages);
        assert!(
            !prefix.is_empty(),
            "at least the import-verify prefix executes"
        );
        assert!(
            prefix.len() < stages.len(),
            "if this fires the spine is complete; update the ratchet and celebrate"
        );
    }
}
