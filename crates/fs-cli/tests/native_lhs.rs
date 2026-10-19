//! Native LHS product regressions: imported geometry, material cards, retained
//! child solves and resumption use the ordinary product path, not a mock model.

#[path = "../src/json_read.rs"]
mod json_read;

use json_read::JsonValue as J;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_cli::exit;
use fs_uq::{CorrelationModel, LhsExecution, LhsPlan, ParameterUncertainty};

const STUDY: &str = r#"(fsim-uncertainty-study
    :version 1 :project "cooling-reference.fsim" :samples 6 :seed 29
    :wall-time 120s :method latin-hypercube
    :lhs (replicated-latin-hypercube :replicates 2 :samples-per-replicate 3)
    :correlation independent :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces ()
    :parameters (
        (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
        (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;
static NEXT: AtomicU64 = AtomicU64::new(0);

fn command(args: &[&str], expected: u8) -> J {
    let out = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(out.exit_code, expected, "{args:?}\n{}\n{}", out.stdout, out.stderr);
    J::parse(&out.stdout).unwrap()
}

struct Fixture { dir: PathBuf, sources: PathBuf, ledger: PathBuf }
impl Fixture {
    fn new() -> Self {
        let dir = loop {
            let dir = std::env::temp_dir().join(format!("fs-cli-native-lhs-{}-{}",
                std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            match std::fs::create_dir(&dir) {
                Ok(()) => break dir,
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(e) => panic!("scratch: {e}"),
            }
        };
        let sources = dir.join("sources");
        std::fs::create_dir(&sources).unwrap();
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        for file in ["cooling-reference.fsim", "plate.stl", "aa6061.fsmcdpk"] {
            std::fs::copy(reference.join(file), sources.join(file)).unwrap();
        }
        std::fs::write(sources.join("study.fsim"), STUDY).unwrap();
        Self { ledger: dir.join("study.db"), dir, sources }
    }
    fn study(&self, cap: Option<&str>, expected: u8) -> J {
        let source = self.sources.join("study.fsim");
        let mut args = vec!["--json", "study", source.to_str().unwrap(), self.ledger.to_str().unwrap()];
        if let Some(cap) = cap { args.extend(["--budget", cap]); }
        command(&args, expected)
    }
    fn resume(&self, run: &str, cap: Option<&str>, expected: u8) -> J {
        let mut args = vec!["--json", "study", "--resume", run, self.ledger.to_str().unwrap()];
        if let Some(cap) = cap { args.extend(["--budget", cap]); }
        command(&args, expected)
    }
    fn artifact(&self, hash: &str) -> Vec<u8> {
        let hash = fs_blake3::ContentHash::from_hex(hash).unwrap();
        fs_ledger::Ledger::open(self.ledger.to_str().unwrap()).unwrap()
            .get_artifact(&hash).unwrap().unwrap()
    }
    fn report(&self, result: &J) -> J {
        let receipt = result.get("receipt").unwrap();
        assert_eq!(receipt.f64_field("samples_planned"), Some(6.0));
        J::parse(std::str::from_utf8(&self.artifact(receipt.str_field("report_json").unwrap())).unwrap()).unwrap()
    }
}

fn rows(report: &J) -> &[J] { report.get("observations").unwrap().as_array().unwrap() }
fn point(row: &J) -> Vec<f64> {
    row.get("parameters").unwrap().as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}
fn signature(report: &J) -> Vec<(Vec<u64>, u64)> {
    rows(report).iter().map(|row| (point(row).iter().map(|x| x.to_bits()).collect(),
        row.f64_field("value_k").unwrap().to_bits())).collect()
}

#[test]
fn g1_native_lhs_retains_actual_solves_and_uses_the_declared_sampler() {
    let f = Fixture::new();
    let result = f.study(None, exit::SUCCESS);
    let report = f.report(&result);
    assert_eq!(report.str_field("method"), Some("latin-hypercube"));
    assert_eq!(report.str_field("authority"), Some("Estimated"));
    assert_eq!(report.get("statistics"), Some(&J::Null));
    assert!(report.get("qmc").is_none());
    assert!(report.get("compliance").is_none());
    let lhs = report.get("lhs").unwrap();
    assert_eq!(lhs.str_field("sampler"), Some("philox-latin-hypercube"));
    assert_eq!(lhs.f64_field("samples_in_estimate"), Some(6.0));
    assert_eq!(lhs.f64_field("completed_replicates"), Some(2.0));
    let observations = rows(&report);
    assert_eq!(observations.len(), 6);
    let plan = LhsPlan {
        target_qoi: "temperature-max".into(),
        compliance_threshold: report.f64_field("temperature_limit_k"),
        parameters: vec![ParameterUncertainty::uniform("power", 4.0, 6.0, "W"),
            ParameterUncertainty::uniform("ambient", 294.0, 300.0, "K")],
        correlation: CorrelationModel::Independent,
        replicates: 2, samples_per_replicate: 3, seed: 29,
    };
    let mut ordinal = 0;
    let expected = LhsExecution::new(&plan).unwrap().advance(6, || false, |parameters| {
        let row = &observations[ordinal];
        ordinal += 1;
        assert_eq!(parameters.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            point(row).iter().map(|x| x.to_bits()).collect::<Vec<_>>());
        let value = row.f64_field("value_k").unwrap();
        let child = J::parse(std::str::from_utf8(&f.artifact(row.str_field("qoi_receipt").unwrap())).unwrap()).unwrap();
        let child_qoi = child.get("qoi").unwrap().as_array().unwrap().iter()
            .find(|qoi| qoi.str_field("name") == Some("temperature-max")).unwrap();
        assert_eq!(child_qoi.f64_field("value").unwrap().to_bits(), value.to_bits());
        Ok::<f64, &'static str>(value)
    });
    let mean = lhs.get("estimate").unwrap();
    let expected_mean = expected.estimate.unwrap();
    assert_eq!(mean.f64_field("mean_k"), Some(expected_mean.mean));
    assert_eq!(mean.f64_field("sampling_standard_error_k"), expected_mean.standard_error);
    assert_eq!(lhs.get("compliance_estimate").unwrap().f64_field("probability_of_compliance"),
        expected.compliance.map(|c| c.mean));
    for replicate in observations.chunks_exact(3) {
        for (dimension, low, width) in [(0, 4.0, 2.0), (1, 294.0, 6.0)] {
            let mut bins = replicate.iter().map(|row| ((point(row)[dimension] - low) / width * 3.0).floor() as usize).collect::<Vec<_>>();
            bins.sort();
            assert_eq!(bins, [0, 1, 2]);
        }
    }
    let receipt = result.get("receipt").unwrap();
    let package = String::from_utf8(f.artifact(receipt.str_field("package").unwrap())).unwrap();
    let parsed = fs_package::EvidencePackage::from_json(&package).unwrap();
    assert!(fs_checker::check(&parsed).passed());
    assert!(package.contains("fixed-count-native-latin-hypercube-between-replicate-standard-error"));
    let html = String::from_utf8(f.artifact(receipt.str_field("report_html").unwrap())).unwrap();
    assert!(html.contains("Replicated Latin-hypercube sampling"));
    let exported = command(&["--json", "report", result.str_field("run").unwrap(), f.ledger.to_str().unwrap()], exit::SUCCESS);
    assert_eq!(exported.str_field("verification"), Some("sealed-evidence"));
}

#[test]
fn g5_partial_hypercubes_resume_from_retained_inputs_without_repeating_paid_work() {
    let whole = Fixture::new();
    let complete = whole.report(&whole.study(None, exit::SUCCESS));
    let split = Fixture::new();
    let first = split.study(Some("2"), exit::BUDGET);
    let partial = split.report(&first);
    let lhs = partial.get("lhs").unwrap();
    assert_eq!(lhs.f64_field("completed_replicates"), Some(0.0));
    assert_eq!(lhs.f64_field("partial_replicate_samples"), Some(2.0));
    assert_eq!(lhs.get("estimate"), Some(&J::Null));
    std::fs::rename(&split.sources, split.dir.join("parked-sources")).unwrap();
    let second = split.resume(first.str_field("run").unwrap(), Some("1"), exit::BUDGET);
    let one = split.report(&second);
    assert_eq!(&signature(&one)[..2], signature(&partial).as_slice());
    let lhs = one.get("lhs").unwrap();
    assert_eq!(lhs.f64_field("completed_replicates"), Some(1.0));
    assert_eq!(lhs.get("estimate").unwrap().get("sampling_standard_error_k"), Some(&J::Null));
    let last = split.resume(second.str_field("run").unwrap(), None, exit::SUCCESS);
    let resumed = split.report(&last);
    assert_eq!(signature(&resumed), signature(&complete));
    assert_eq!(resumed.get("lhs"), complete.get("lhs"));
    let terminal = split.resume(last.str_field("run").unwrap(), None, exit::SUCCESS);
    assert_eq!(terminal, last);
}

#[test]
fn g0_invalid_lhs_design_refuses_before_creating_a_native_ledger() {
    let f = Fixture::new();
    std::fs::write(f.sources.join("study.fsim"), STUDY.replace(":samples 6", ":samples 7")).unwrap();
    let result = f.study(None, exit::REFUSED);
    assert_eq!(result.str_field("status"), Some("refused"));
    assert!(!f.ledger.exists());
}
