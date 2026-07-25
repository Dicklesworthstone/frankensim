//! Cross-repository compatibility report: which sibling moved, and what to run.
//!
//! `check-constellation` compares one FNV-1a-64 aggregate over
//! `lib=version@git_head` lines, so a drift verdict can say only that *a*
//! constellation repo moved. On 2026-07-24 four agents lost hours to a drift
//! that gate had already detected but could not name. This report answers the
//! two questions the aggregate cannot: WHICH sibling moved, and WHICH
//! FrankenSim tests must go green before the new pin may be accepted.
//!
//! It is a report, not a `check-all` gate — matching the existing decision to
//! keep `check-constellation` out of `check-all`, since a pin drift is a
//! release-train event rather than a source-tree defect.

use std::path::Path;

use fs_govern::compatibility::{PinRow, SURFACES, moved, pin_delta, render_registry};

pub(crate) const CHECK: &str = "compatibility-report";

/// Print the registry, the per-sibling pin delta, and the suite selectors.
///
/// With `candidate = None` the comparison is against the LIVE checkouts, which
/// answers "what has drifted?". With a candidate lock path it is against a
/// PROPOSED pin set, which answers "what would this bump change, and what must
/// go green first?" — the question a release train actually asks. A candidate
/// lock is parsed by the same canonical reader as the tracked one, so a
/// malformed proposal is refused rather than silently half-read.
pub(crate) fn report(root: &Path, candidate: Option<&Path>) -> Result<(), String> {
    let lock_text = super::read_constellation_lock(&root.join("constellation.lock"))?;
    let (_lock_hash, rows) = super::parse_lock_rows(&lock_text)?;

    let recorded: Vec<PinRow> = rows
        .iter()
        .map(|row| PinRow::new(&row.lib, &row.version, &row.git_head))
        .collect();

    let (against, basis) = match candidate {
        Some(path) => {
            let text = super::read_constellation_lock(path).map_err(|error| {
                format!("candidate lock {} is unreadable: {error}", path.display())
            })?;
            let (_hash, candidate_rows) = super::parse_lock_rows(&text).map_err(|error| {
                format!(
                    "candidate lock {} is not a canonical constellation lock: {error}",
                    path.display()
                )
            })?;
            (
                candidate_rows
                    .iter()
                    .map(|row| PinRow::new(&row.lib, &row.version, &row.git_head))
                    .collect::<Vec<PinRow>>(),
                format!("candidate lock {}", path.display()),
            )
        }
        None => (
            super::constellation_entries(root)?
                .iter()
                .map(|entry| PinRow::new(&entry.lib, &entry.version, &entry.git_head))
                .collect::<Vec<PinRow>>(),
            "live checkouts".to_string(),
        ),
    };

    let deltas = pin_delta(&recorded, &against);
    let movers = moved(&deltas);

    println!("# Constellation compatibility report\n");
    println!("## Pin delta (recorded lock versus {basis})\n");
    for delta in &deltas {
        println!("- {}", delta.describe());
    }
    println!();

    if movers.is_empty() {
        println!("All siblings are at their recorded pins.\n");
    } else {
        println!(
            "{} sibling(s) are OFF-PIN. A pin change is a tested, recorded event: run each \
             affected surface below against the candidate pins, then adjudicate the bump.\n",
            movers.len()
        );
        for delta in &movers {
            let Some(surface) = fs_govern::compatibility::surface(&delta.lib) else {
                println!(
                    "- {}: NO REGISTERED SURFACE — this sibling cannot be adjudicated",
                    delta.lib
                );
                continue;
            };
            match surface.selector() {
                Some(selector) => println!(
                    "- {} ({}): `{selector}`",
                    delta.lib,
                    surface.priority.slug()
                ),
                None => println!(
                    "- {} ({}): NO COVERAGE — {}",
                    delta.lib,
                    surface.priority.slug(),
                    surface.no_test_reason.unwrap_or("no reason recorded")
                ),
            }
        }
        println!();
    }

    println!("## Registered compatibility surfaces\n");
    print!("{}", render_registry());
    println!();

    let uncovered: Vec<&str> = SURFACES
        .iter()
        .filter(|surface| surface.tests.is_empty())
        .map(|surface| surface.lib)
        .collect();
    if uncovered.is_empty() {
        println!("Every registered sibling carries compatibility coverage.");
    } else {
        println!(
            "Uncovered surfaces (a bump moving any of these is refused): {}",
            uncovered.join(", ")
        );
    }

    Ok(())
}
