//! Native studies through the program executor, not a sampler/solver stub.
#[path = "../src/json_read.rs"]
mod json_read;
use json_read::JsonValue as J;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const SOURCE: &str = r#"(fsim-uncertainty-study
 :version 1 :project "cooling-reference.fsim" :samples 4 :seed 29
 :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
 :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
 :materials ("aa6061.fsmcdpk") :interfaces ()
 :parameters (
  (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
  (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

fn invoke(args: &[&str], expected: u8) -> (J, String) {
    let out = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(out.exit_code, expected, "{args:?}\n{}\n{}", out.stdout, out.stderr);
    (J::parse(&out.stdout).unwrap(), out.stderr)
}
fn program(body: &str) -> String {
    format!(r#"(study "native-study-program"
 (seed 0x7) (versions (constellation :lock "2026-07"))
 (budget (wall 120s) (mem 64MiB))
 (capability :cores 1 :mem 64MiB :wall 120s :ops (cooling.*))
 (let project (cooling.project "assets/cooling-reference.fsim"))
 {body})"#)
}
fn source(qmc: bool, copula: bool, control: bool) -> String {
    let mut s = SOURCE.to_string();
    if qmc { s = s.replace("monte-carlo", "quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)"); }
    if copula { s = s.replace("independent", "(gaussian-copula :latent-correlation ((1 0.5) (0.5 1)))"); }
    if control { s = s.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)"); }
    s
}
struct Fixture { root: PathBuf, assets: PathBuf }
impl Fixture {
    fn new(source: &str) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let root = loop {
            let p = std::env::temp_dir().join(format!("fs-script-study-{}-{}",
                std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            match std::fs::create_dir(&p) {
                Ok(()) => break p,
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => panic!("scratch: {e}"),
            }
        };
        let assets = root.join("assets"); std::fs::create_dir(&assets).unwrap();
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        for name in ["cooling-reference.fsim", "plate.stl", "aa6061.fsmcdpk"] {
            std::fs::copy(reference.join(name), assets.join(name)).unwrap();
        }
        std::fs::write(assets.join("study.fsim"), source).unwrap();
        Self { root, assets }
    }
    fn run(&self, text: &str, tag: &str, expected: u8) -> (J, String) {
        let path = self.root.join(format!("{tag}.fs"));
        std::fs::write(&path, text).unwrap();
        invoke(&["--json", "run", path.to_str().unwrap(), self.ledger(tag).to_str().unwrap()], expected)
    }
    fn direct(&self) -> J {
        invoke(&["--json", "study", self.assets.join("study.fsim").to_str().unwrap(),
            self.ledger("direct").to_str().unwrap()], fs_cli::exit::SUCCESS).0
    }
    fn ledger(&self, tag: &str) -> PathBuf { self.root.join(format!("{tag}.db")) }
    fn report(&self, tag: &str, native: &J) -> Vec<u8> {
        let hash = native.get("receipt").unwrap().str_field("report_json").unwrap();
        fs_ledger::Ledger::open(self.ledger(tag).to_str().unwrap()).unwrap()
            .get_artifact(&fs_blake3::ContentHash::from_hex(hash).unwrap()).unwrap().unwrap()
    }
}
fn child(program: &J) -> &J {
    let steps = program.get("steps").unwrap().as_array().unwrap();
    assert_eq!(steps.len(), 1);
    assert_eq!(steps[0].str_field("verb"), Some("cooling.study"));
    steps[0].get("result").unwrap()
}

#[test]
fn programs_preserve_native_study_identities_samples_and_estimates() {
    let mut cases = vec![source(false,false,false), source(true,false,false),
        source(false,true,false), source(true,true,false), source(true,true,true)];
    cases.push(SOURCE.replace(":samples 4", ":samples 8").replace("monte-carlo", "sobol-sensitivity"));
    for text in cases {
        let f = Fixture::new(&text);
        let direct = f.direct();
        let declaration = fs_project::uncertainty::UncertaintyStudy::parse(&text).unwrap();
        let pin = fs_blake3::hash_bytes(declaration.canonical().as_bytes()).to_hex();
        let (result, _) = f.run(&program(&format!(
            "(cooling.study project :source \"assets/study.fsim\" :hash \"{pin}\")")),
            "script", fs_cli::exit::SUCCESS);
        assert_eq!(result.str_field("status"), Some("completed"));
        let native = child(&result);
        // Receipt pointers include elapsed wall charges; physical model and
        // report identities, not separately timed invocation pointers, match.
        for field in ["study_id", "model"] {
            assert_eq!(native.get("receipt").unwrap().str_field(field),
                direct.get("receipt").unwrap().str_field(field));
        }
        let actual = f.report("script", native);
        assert_eq!(actual, f.report("direct", &direct), "programs must execute the same native experiment");
        let report = J::parse(std::str::from_utf8(&actual).unwrap()).unwrap();
        assert_eq!(report.str_field("seed"), Some("29"), "do not replace sampling seed with physical seed 7");
    }
}

#[test]
fn interrupted_calibration_stops_later_steps_and_resumes_without_program_or_assets() {
    let f = Fixture::new(&source(true,true,true));
    let direct = f.direct();
    let text = program("(cooling.study project :source \"assets/study.fsim\" :budget 1)\n (cooling.run project :materials (\"assets/aa6061.fsmcdpk\"))");
    let (prefix, _) = f.run(&text, "split", fs_cli::exit::BUDGET);
    assert_eq!(prefix.str_field("status"), Some("stopped"));
    assert_eq!(prefix.f64_field("steps_planned"), Some(2.0));
    let native = child(&prefix); // no second step after the calibration-only budget
    let report = f.report("split", native);
    let report = J::parse(std::str::from_utf8(&report).unwrap()).unwrap();
    assert!(report.get("observations").unwrap().as_array().unwrap().is_empty());
    std::fs::rename(&f.assets, f.root.join("moved-assets")).unwrap();
    std::fs::rename(f.root.join("split.fs"), f.root.join("moved-program.fs")).unwrap();
    let (complete, _) = invoke(&["--json", "study", "--resume", native.str_field("run").unwrap(),
        f.ledger("split").to_str().unwrap()], fs_cli::exit::SUCCESS);
    assert_eq!(f.report("split", &complete), f.report("direct", &direct));
    let (again, _) = invoke(&["--json", "study", "--resume", complete.str_field("run").unwrap(),
        f.ledger("split").to_str().unwrap(), "--budget", "0"], fs_cli::exit::SUCCESS);
    assert_eq!(again, complete);
}

#[test]
fn all_study_clauses_bind_before_an_earlier_import_can_write_a_ledger() {
    let f = Fixture::new(SOURCE);
    let prefix = "(cooling.import project :sources (\"assets/plate.stl\") :unit \"m\" :max-hole-edges 0)\n";
    let step = "(cooling.study project :source \"assets/study.fsim\" :budget 0)";
    for (i, tail) in [
        step.replace(":budget 0", ":budget 257"),
        step.replace(":budget 0", ":hash 7"),
        step.replace(":budget 0", &format!(":hash \"{}\"", "0".repeat(64))),
        step.replace("project :source", "unbound :source"),
    ].iter().enumerate() {
        f.run(&program(&format!("{prefix}{tail}")), &format!("bad-{i}"), fs_cli::exit::REFUSED);
        assert!(!f.ledger(&format!("bad-{i}")).exists());
    }
    let too_small = program(&format!("{prefix}{step}")).replace("wall 120s", "wall 60s");
    let (_, error) = f.run(&too_small, "wall", fs_cli::exit::REFUSED);
    assert!(error.contains("wall allowance")); assert!(!f.ledger("wall").exists());
    std::fs::write(f.assets.join("study.fsim"), SOURCE.replace(":entity \"air\"", ":entity \"missing\"" )).unwrap();
    f.run(&program(&format!("{prefix}{step}")), "target", fs_cli::exit::REFUSED);
    assert!(!f.ledger("target").exists());
    let mut other = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
        f.assets.join("cooling-reference.fsim")).unwrap()).unwrap().decoded.spec;
    other.power.as_mut().unwrap()[0].watts.value = 6.0;
    std::fs::write(f.assets.join("other.fsim"), fs_project::print_sexpr(&other).unwrap()).unwrap();
    std::fs::write(f.assets.join("study.fsim"), SOURCE.replace("cooling-reference.fsim", "other.fsim")).unwrap();
    let (_, error) = f.run(&program(&format!("{prefix}{step}")), "project", fs_cli::exit::REFUSED);
    assert!(error.contains("physical project differs")); assert!(!f.ledger("project").exists());
}

#[test]
fn physical_refusal_remains_a_failed_native_step_not_a_completed_program() {
    let f = Fixture::new(&SOURCE.replace(":materials (\"aa6061.fsmcdpk\")", ":materials ()"));
    let (result, _) = f.run(&program("(cooling.study project :source \"assets/study.fsim\")\n (cooling.run project)"),
        "refusal", fs_cli::exit::REFUSED);
    assert_eq!(result.str_field("status"), Some("stopped"));
    let native = child(&result);
    assert_eq!(native.str_field("status"), Some("refused"));
    let report = f.report("refusal", native);
    let report = J::parse(std::str::from_utf8(&report).unwrap()).unwrap();
    assert!(report.get("observations").unwrap().as_array().unwrap().is_empty());
    assert_eq!(report.f64_field("evaluations_attempted"), Some(1.0));
}

#[test]
fn missing_or_invalid_late_study_assets_refuse_before_any_program_step() {
    let prefix = "(cooling.import project :sources (\"assets/plate.stl\") :unit \"m\" :max-hole-edges 0)\n";
    let step = "(cooling.study project :source \"assets/study.fsim\" :budget 0)";
    for missing in ["missing.stl", "missing.fsmcdpk"] {
        let source = if missing.ends_with(".stl") {
            SOURCE.replace("plate.stl", missing)
        } else { SOURCE.replace("aa6061.fsmcdpk", missing) };
        let f = Fixture::new(&source);
        let (_, error) = f.run(&program(&format!("{prefix}{step}")), "missing", fs_cli::exit::REFUSED);
        assert!(error.contains(missing), "{error}");
        assert!(!f.ledger("missing").exists(), "preparation must precede even an otherwise valid import");
    }
    let f = Fixture::new(SOURCE);
    std::fs::write(f.assets.join("aa6061.fsmcdpk"), b"invalid normalized card pack").unwrap();
    let (_, error) = f.run(&program(&format!("{prefix}{step}")), "bad-card", fs_cli::exit::REFUSED);
    assert!(error.contains("frankenscript-native-study"), "{error}");
    assert!(!f.ledger("bad-card").exists());
}
