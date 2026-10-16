//! Bind native studies to the exact bounded model snapshots that execute.
//! Physical and sampling seeds stay distinct; the ordinary native producer
//! owns calibration, numerical evaluation, retained results and recovery.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use fs_blake3::ContentHash;
use crate::study::{PreparedStudy, StudyPins};
use super::{Binder, Node, ProjectBinding, Refusals, count_of, string_of};
use crate::{CommandOutput, OutputMode};

#[derive(Debug)]
pub(super) struct StudyStep {
    model: Rc<PreparedStudy>,
    budget: Option<String>,
}

impl StudyStep {
    pub(super) fn run(&self, ledger: &Path, mode: OutputMode) -> CommandOutput {
        self.model.run(ledger, self.budget.as_deref(), mode)
    }
}

/// Program-local input snapshots. The same declared source path always means
/// the first admitted snapshot, even if the filesystem changes between clauses.
/// Distinct paths share one aggregate input-byte allowance (not an RSS claim).
pub(super) struct SnapshotCache {
    limit: u64,
    used: u64,
    models: BTreeMap<PathBuf, Rc<PreparedStudy>>,
}

impl SnapshotCache {
    pub(super) fn new(memory: Option<u64>) -> Self {
        // Match the native model's input-storage allocation; leave room for
        // decoded state and numerical work. No implicit unlimited grant.
        Self { limit: memory.map_or(0, |bytes| (bytes / 4).min(256 * 1024 * 1024)),
            used: 0, models: BTreeMap::new() }
    }

    fn prepare(&mut self, source: PathBuf, pins: StudyPins) -> Result<Rc<PreparedStudy>, String> {
        if let Some(model) = self.models.get(&source) {
            model.check(pins)?;
            return Ok(Rc::clone(model));
        }
        let remaining = self.limit.checked_sub(self.used)
            .filter(|remaining| *remaining > 0)
            .ok_or("native study snapshots exceed the program's declared input-storage allowance; declare a sufficient (budget (mem ...))")?;
        let model = PreparedStudy::load(&source, pins, remaining)?;
        let used = self.used.checked_add(model.input_bytes())
            .filter(|used| *used <= self.limit)
            .ok_or("native study snapshot input-byte accounting exceeded the program allowance")?;
        let model = Rc::new(model);
        // A refused load/pin never publishes a cache entry or spends storage.
        self.models.insert(source, Rc::clone(&model));
        self.used = used;
        Ok(model)
    }
}

#[derive(Debug, PartialEq)]
struct Options<'a> {
    source: &'a str,
    hash: Option<&'a str>,
    budget: Option<usize>,
}

fn options<'a>(named: &[(&str, &'a Node)]) -> Result<Options<'a>, String> {
    let (mut source, mut hash, mut budget) = (None, None, None);
    for &(key, value) in named {
        match key {
            "source" if source.is_none() => {
                source = Some(string_of(value).filter(|s| !s.is_empty()
                    && !s.chars().any(char::is_control))
                    .ok_or("cooling.study :source needs a nonempty path string")?);
            }
            "hash" if hash.is_none() => {
                let pin = string_of(value)
                    .ok_or("cooling.study :hash needs a 64-hex canonical-study hash string")?;
                if ContentHash::from_hex(pin).is_none() {
                    return Err("cooling.study :hash is not a 64-hex content hash".into());
                }
                hash = Some(pin);
            }
            "budget" if budget.is_none() => {
                budget = Some(count_of(value).and_then(|v| usize::try_from(v).ok())
                    .filter(|v| *v <= fs_project::uncertainty::MAX_SAMPLES)
                    .ok_or("cooling.study :budget needs an exact evaluation count in 0..=256")?);
            }
            _ => return Err(format!("cooling.study has an unknown or repeated keyword :{key}")),
        }
    }
    Ok(Options {
        source: source.ok_or("cooling.study requires :source \"study.fsim\"")?,
        hash,
        budget,
    })
}

fn plan(binder: &mut Binder<'_>, project: &ProjectBinding, options: Options<'_>) -> Result<StudyStep, String> {
    let source = binder.resolve(options.source);
    let pins = StudyPins {
        project: ContentHash::from_hex(&project.hash).ok_or("invalid bound project identity")?,
        source: options.hash.map(|hash| ContentHash::from_hex(hash)
            .ok_or("invalid canonical-study pin")).transpose()?,
        wall_seconds: binder.wall_seconds,
    };
    let model = binder.native_studies.prepare(source, pins)?;
    Ok(StudyStep { model, budget: options.budget.map(|value| value.to_string()) })
}

pub(super) fn bind(
    binder: &mut Binder<'_>,
    project: &ProjectBinding,
    named: &[(&str, &Node)],
    refusals: &mut Refusals,
) -> Option<StudyStep> {
    match options(named).and_then(|options| plan(binder, project, options)) {
        Ok(step) => Some(step),
        Err(error) => {
            refusals.push("frankenscript-native-study", error,
                "use (cooling.study project :source \"study.fsim\" [:hash \"canonical-study-hash\"] [:budget N]); bind the same physical project and provide readable assets within the program's wall and memory allowances");
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &str) -> Node {
        fs_ir::sexpr::parse(text).unwrap()
    }

    fn arguments(node: &Node) -> Vec<(&str, &Node)> {
        node.items().unwrap()[1..].chunks_exact(2).map(|pair| {
            let fs_ir::ast::NodeKind::Keyword(key) = &pair[0].kind else { panic!("keyword") };
            (key.as_str(), &pair[1])
        }).collect()
    }

    #[test]
    fn evaluation_budget_keeps_zero_and_does_not_override_the_sampling_plan() {
        let node = parse("(options :source \"nested/study.fsim\" :budget 0)");
        let opts = options(&arguments(&node)).unwrap();
        assert_eq!(opts, Options { source: "nested/study.fsim", hash: None, budget: Some(0) });
        let node = parse("(options :source \"study.fsim\")");
        assert_eq!(options(&arguments(&node)).unwrap().budget, None);
        let node = parse("(options :source \"study.fsim\" :budget 17)");
        assert_eq!(options(&arguments(&node)).unwrap().budget, Some(17));
    }

    #[test]
    fn malformed_pins_budgets_and_attempted_inline_model_changes_refuse() {
        for args in [
            ":source \"study.fsim\" :hash 7", ":source \"study.fsim\" :hash \"bad\"",
            ":source \"study.fsim\" :budget -1", ":source \"study.fsim\" :budget 1.5",
            ":source \"study.fsim\" :budget 257",
            ":source \"study.fsim\" :budget 1s", ":source \"study.fsim\" :budget \"1\"",
            ":source \"study.fsim\" :seed 29", ":source \"study.fsim\" :materials ()",
            ":source \"study.fsim\" :method monte-carlo", ":source \"\"", ":budget 0",
            ":source \"a.fsim\" :source \"b.fsim\"", ":source \"a.fsim\" :budget 0 :budget 1",
        ] {
            let node = parse(&format!("(options {args})"));
            assert!(options(&arguments(&node)).is_err(), "{args}");
        }
    }
}

#[cfg(test)]
#[path = "native_study/snapshot_tests.rs"]
mod snapshot_tests;
