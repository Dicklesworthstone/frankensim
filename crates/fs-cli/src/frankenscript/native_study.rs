//! Bind native cooling probability/sensitivity studies without another sampler.
//!
//! The program's root seed belongs to the physical project; the study source
//! explicitly owns a DIFFERENT sampling seed. Neither is replaced here. The
//! native driver owns assets, calibration, evaluation accounting and recovery.

use std::io::Read as _;
use std::path::{Path, PathBuf};

use fs_project::uncertainty::UncertaintyStudy;

use super::{Binder, Node, ProjectBinding, Refusals, count_of, string_of};
use crate::{CommandOutput, OutputMode, read_project_for_solve};

#[derive(Debug)]
pub(super) struct StudyStep {
    source: PathBuf,
    budget: Option<String>,
}

impl StudyStep {
    pub(super) fn run(&self, ledger: &Path, mode: OutputMode) -> CommandOutput {
        // Keep the ordinary native study's complete receipt as the step result,
        // including its resumable run ID on budget exhaustion. Program execution
        // stops on any non-success; later clauses cannot consume a partial study.
        crate::study::study_path(&self.source, ledger, self.budget.as_deref(), mode)
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
                if fs_blake3::ContentHash::from_hex(pin).is_none() {
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

fn read_source(path: &Path) -> Result<UncertaintyStudy, String> {
    let cap = fs_project::uncertainty::MAX_SOURCE_BYTES as u64;
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let metadata = file.metadata().map_err(|e| e.to_string())?;
    if !metadata.is_file() || metadata.len() > cap {
        return Err(format!("native study must be a regular UTF-8 file within {cap} bytes"));
    }
    let mut text = String::new();
    file.take(cap + 1).read_to_string(&mut text).map_err(|e| e.to_string())?;
    if text.len() as u64 > cap {
        return Err("native study exceeded its source limit while being read".into());
    }
    UncertaintyStudy::parse(&text).map_err(|e| format!("{}: {}", e.code, e.detail))
}

fn plan(
    binder: &Binder<'_>,
    project: &ProjectBinding,
    options: Options<'_>,
    mode: OutputMode,
) -> Result<StudyStep, String> {
    let source = binder.resolve(options.source);
    let study = read_source(&source)?;
    let hash = fs_blake3::hash_bytes(study.canonical().as_bytes()).to_hex();
    if options.hash.is_some_and(|pin| pin != hash) {
        return Err(format!("native study hashes to {hash}, not its pinned :hash"));
    }
    if binder.wall_seconds.is_some_and(|wall| wall < study.wall_seconds()) {
        return Err(format!(
            "native study declares {} s, above the program's wall allowance; :budget limits evaluations, not the declared wall allowance",
            study.wall_seconds(),
        ));
    }
    // Paths inside the native source belong to THAT source's directory, not
    // the program's directory. Equal canonical project content is the binding,
    // not path spelling: relocated copies retain the same physical identity.
    let relative = Path::new(study.project_path());
    if relative.is_absolute() {
        return Err("native study project paths must be relative to the study file".into());
    }
    let base_path = source.parent().unwrap_or_else(|| Path::new(".")).join(relative);
    let base = read_project_for_solve(&base_path, mode).map_err(|out| out.stderr)?;
    if base.hash().to_hex() != project.hash {
        return Err("native study's physical project differs from the cooling.project binding".into());
    }
    // Reuse the owner's complete target, support, dependence and policy gates.
    // In particular, a malformed late study clause cannot allow an earlier
    // cooling.import/solve clause to create a ledger first.
    study.bind(&base.spec).map_err(|e| format!("{}: {}", e.code, e.detail))?;
    Ok(StudyStep { source, budget: options.budget.map(|v| v.to_string()) })
}

pub(super) fn bind(
    binder: &Binder<'_>,
    project: &ProjectBinding,
    named: &[(&str, &Node)],
    refusals: &mut Refusals,
) -> Option<StudyStep> {
    let result = options(named).and_then(|options| plan(binder, project, options, refusals.mode));
    match result {
        Ok(step) => Some(step),
        Err(error) => {
            refusals.push("frankenscript-native-study", error,
                "use (cooling.study project :source \"study.fsim\" [:hash \"canonical-study-hash\"] [:budget N]); bind the same physical project and cover the native wall allowance");
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
