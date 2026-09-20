//! Explicit checkpoint-to-finer-grid study fork, using the existing runner.
use super::*;
use fs_topols::refinement::{ProjectedRefinementReport, prolongate_level_set, refine_projected_study};

pub(super) struct Handoff {
    report: ProjectedRefinementReport,
    source_digest: String,
}

impl Handoff {
    pub(super) fn json(&self) -> String {
        let r = &self.report;
        format!(
            "{{\"schema\":\"projected-grid-refinement-v1\",\"operation\":\"new_fine_grid_study\",\"source_checkpoint_hash\":{},\"coarse_level\":{},\"fine_level\":{},\"source_updates\":{},\"source_solves_started\":{},\"coarse_endpoint\":{},\"transferred_area\":{:.17e},\"max_projection_field_change\":{:.17e},\"fine_baseline\":{},\"cross_grid_improvement_claimed\":false,\"discretization_error_bound\":false}}",
            quoted(&self.source_digest), r.coarse_level, r.fine_level, r.source_updates,
            r.source_solves_started, state_json(&r.coarse_endpoint), r.transferred_area,
            r.max_projection_change, state_json(&r.fine_baseline),
        )
    }
}

struct Request {
    updates: usize,
    max_solves: usize,
    recovery_solves: usize,
    pause_after: Option<usize>,
}

fn request(args: &[String]) -> Result<Request, Box<dyn Error>> {
    if args.len() != 8 && args.len() != 10 {
        return Err("usage: --projected --refine SOURCE.fscp NEW_OUTPUT --updates N --max-solves N --recovery-solves N [--pause-after N]".into());
    }
    let (mut updates, mut max_solves, mut recovery, mut pause) = (None, None, None, None);
    let mut i = 2;
    while i < args.len() {
        let value = args[i + 1].parse::<usize>()?;
        match args[i].as_str() {
            "--updates" if updates.is_none() && (1..=200).contains(&value) => updates = Some(value),
            "--max-solves" if max_solves.is_none() && value > 0 => max_solves = Some(value),
            "--recovery-solves" if recovery.is_none() && (1..=128).contains(&value) => recovery = Some(value),
            "--pause-after" if pause.is_none() && value <= 200 => pause = Some(value),
            _ => return Err("unknown, repeated or out-of-range refinement option; physical policies are inherited".into()),
        }
        i += 2;
    }
    Ok(Request {
        updates: updates.ok_or("refinement requires explicit --updates")?,
        max_solves: max_solves.ok_or("refinement requires explicit --max-solves")?,
        recovery_solves: recovery.ok_or("refinement requires explicit --recovery-solves")?,
        pause_after: pause,
    })
}

pub(super) fn run(args: &[String]) -> Result<u8, Box<dyn Error>> {
    let request = request(args)?;
    let output = Path::new(&args[1]);
    if output.try_exists()? { return Err("output directory already exists; refusing to overwrite it".into()); }
    let bytes = checkpoint::load(Path::new(&args[0]))?;
    let source_digest = fs_ledger::hash_bytes(&bytes).to_string();
    let mut remaining = request.recovery_solves;
    let coarse = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut remaining)
        .map_err(|error| format!("{error}; recovery_solves_started={}", request.recovery_solves - remaining))?;
    let spent = request.recovery_solves - remaining;
    if !(2..=6).contains(&coarse.settings().level)
        || !(1..=16).contains(&coarse.controls().max_candidates)
    {
        return Err(format!("refinement exceeds this executable's coarse-level 2..=6 or candidate 1..=16 bounds; recovery_solves_started={spent}").into());
    }
    let (fine, report) = refine_projected_study(&coarse, request.updates, request.max_solves)
        .map_err(|error| format!("{error}; recovery_solves_started={spent}; fine initialization not completed"))?;
    // Export the actual prolonged INPUT, not the repaired fine baseline twice.
    let input = prolongate_level_set(coarse.geometry(), coarse.fixed_nodes())?.geometry;
    let handoff = Handoff { report, source_digest };
    run_optimizer(output, fine, &input, checkpoint::Options {
        enabled: true, pause_after: request.pause_after, recovery_solves: spent,
    }, Some(&handoff))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refinement_requires_new_work_and_rejects_policy_replacement_before_io() {
        let valid = ["source", "out", "--max-solves", "64", "--updates", "2", "--recovery-solves", "4"]
            .map(str::to_string);
        let parsed = request(&valid).unwrap();
        assert_eq!((parsed.updates, parsed.max_solves, parsed.recovery_solves), (2, 64, 4));
        for (index, value) in [(3, "0"), (5, "0"), (5, "201"), (7, "129"),
                               (2, "--updates"), (2, "--stress-limit"), (2, "--area")] {
            let mut bad = valid.clone();
            bad[index] = value.into();
            assert!(request(&bad).is_err());
        }
        assert!(request(&valid[..6]).is_err());
        let mut paused = valid.to_vec();
        paused.extend(["--pause-after".into(), "0".into()]);
        assert_eq!(request(&paused).unwrap().pause_after, Some(0));
    }
}
