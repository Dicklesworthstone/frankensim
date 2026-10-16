//! Prescribed-wall probability studies through actual contact/material solves.
use super::*;

fn fixture(nonlinear: bool) -> Fixture {
    let mut f = Fixture::new(nonlinear, false);
    fixed_value(&mut f, "cold", 285.0);
    f.project.power.as_mut().unwrap()[0].watts.value = 5.0;
    let (bytes, card) = contact_pack(0.13);
    let binding = &mut f.project.interface_cards.as_mut().unwrap()[0];
    binding.card = card; binding.claim = None;
    let sources = f.dir.join("sources");
    std::fs::create_dir(&sources).unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    for name in ["cold-body.stl", "hot-body.stl"] {
        std::fs::copy(root.join("examples/contact-pair").join(name), sources.join(name)).unwrap();
    }
    std::fs::copy(f.dir.join("solid.fsmcdpk"), sources.join("solid.fsmcdpk")).unwrap();
    std::fs::write(sources.join("contact.fsintpk"), bytes).unwrap();
    std::fs::write(sources.join("contact-pair.fsim"), fs_project::print_sexpr(&f.project).unwrap()).unwrap();
    f
}

fn source(copula: bool, qmc: bool, controlled: bool) -> String {
    let version = if controlled {3} else {1};
    let control = if controlled {":mean-control (nominal-adjoint :max-solves 1)"} else {""};
    let correlation = if copula {"(gaussian-copula :latent-correlation ((1.0 0.4) (0.4 1.0)))"}
        else {"independent"};
    let method = if qmc {"quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)"}
        else {"monte-carlo"};
    format!(r#"(fsim-uncertainty-study
 :version {version} :project "contact-pair.fsim" :samples 4 :seed 59
 :wall-time 120s :method {method} :correlation {correlation} :qoi "temperature-max"
 :geometry ((mesh :role "cold-body" :path "cold-body.stl" :unit "m" :max-hole-edges 0)
            (mesh :role "hot-body" :path "hot-body.stl" :unit "m" :max-hole-edges 0))
 :materials ("solid.fsmcdpk") :interfaces ("contact.fsintpk") {control}
 :parameters (
  (uniform :name "wall" :target fixed-temperature :entity "cold" :low 280K :high 290K)
  (uniform :name "power" :target power :entity "hot" :low 4W :high 6W)))"#)
}

fn study(f: &Fixture, resume: Option<&str>, budget: Option<&str>, expected: u8) -> J {
    let source = f.dir.join("sources/study.fsim");
    let ledger = f.dir.join("study.db");
    let mut args = vec!["--json", "study"];
    if let Some(run) = resume { args.extend(["--resume", run]); }
    else { args.push(source.to_str().unwrap()); }
    args.push(ledger.to_str().unwrap());
    if let Some(budget) = budget { args.extend(["--budget", budget]); }
    let output = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(output.exit_code, expected, "{}\n{}", output.stdout, output.stderr);
    J::parse(&output.stdout).unwrap()
}
fn report(f: &Fixture, run: &J) -> J {
    artifact(&f.dir.join("study.db"), run.get("receipt").unwrap().str_field("report_json").unwrap())
}
fn parameters(report: &J) -> &[J] {
    report.get("observations").unwrap().as_array().unwrap()
}

#[test]
fn prescribed_studies_keep_samples_and_use_the_complete_lift_for_calibration() {
    for nonlinear in [false, true] {
        for (copula, qmc) in [(false,false), (true,false), (false,true), (true,true)] {
            let raw = fixture(nonlinear);
            std::fs::write(raw.dir.join("sources/study.fsim"), source(copula,qmc,false)).unwrap();
            let raw_result = report(&raw, &study(&raw,None,None,fs_cli::exit::SUCCESS));
            let mut f = fixture(nonlinear);
            std::fs::write(f.dir.join("sources/study.fsim"), source(copula,qmc,true)).unwrap();
            let result = report(&f, &study(&f,None,None,fs_cli::exit::SUCCESS));
            assert_eq!(parameters(&result).len(),4);
            for key in ["observations", "statistics", "qmc"] {
                assert_eq!(result.get(key),raw_result.get(key),"calibration must not alter {key}");
            }
            let control = result.get("mean_control").unwrap();
            assert_eq!(control.str_field("status"),Some("frozen"));
            assert_eq!(control.f64_field("probe_solves_completed"),Some(1.0));
            let gradient = control.get("gradient").unwrap().as_array().unwrap();
            assert_eq!(gradient.len(),2);
            assert!(gradient.iter().all(|v| v.as_f64().unwrap()>0.0));
            if !copula && !qmc {
                // Independent physical input binding, not a second call to the
                // uncertainty adapter whose mapping is under test.
                for (i,sample) in parameters(&result).iter().enumerate() {
                    let values = sample.get("parameters").unwrap().as_array().unwrap();
                    fixed_value(&mut f,"cold",values[0].as_f64().unwrap());
                    f.project.power.as_mut().unwrap()[0].watts.value = values[1].as_f64().unwrap();
                    let (physical,_,_) = f.solve(0.13,None,100+i);
                    let value = physical.get("temperature").unwrap().f64_field("max").unwrap();
                    assert_eq!(value.to_bits(),sample.f64_field("value_k").unwrap().to_bits());
                }
                fixed_value(&mut f,"cold",285.0);
                f.project.power.as_mut().unwrap()[0].watts.value = 5.0;
                let (nominal,_,_) = f.solve(0.13,Some(BOUNDARY_OUTPUT),110);
                let adjoint = nominal.get("nominal_adjoint").unwrap();
                let vertex = adjoint.f64_field("selected_vertex").unwrap() as usize;
                let rows = adjoint.get("parameters").unwrap().as_array().unwrap();
                for (i,target) in ["fixed-temperature","power"].into_iter().enumerate() {
                    let row = rows.iter().find(|r|r.str_field("target")==Some(target)).unwrap();
                    assert_eq!(gradient[i].as_f64().unwrap().to_bits(),row.f64_field("derivative").unwrap().to_bits());
                }
                let mut values = Vec::new();
                for (i,sign) in [-1.0,1.0].into_iter().enumerate() {
                    fixed_value(&mut f,"cold",285.0+sign*0.002);
                    let (_,field,_) = f.solve(0.13,None,111+i);
                    values.push(field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap());
                }
                let expected = (values[1]-values[0])/0.004;
                assert!((gradient[0].as_f64().unwrap()-expected).abs()<3e-4*expected.abs().max(0.01));
            }
        }
    }
}

#[test]
fn prescribed_calibration_resumes_without_the_original_source_directory() {
    let complete = fixture(true);
    let text = source(true,true,true);
    std::fs::write(complete.dir.join("sources/study.fsim"),&text).unwrap();
    let expected = report(&complete,&study(&complete,None,None,fs_cli::exit::SUCCESS));
    let split = fixture(true);
    std::fs::write(split.dir.join("sources/study.fsim"),text).unwrap();
    let paused = study(&split,None,Some("1"),fs_cli::exit::BUDGET);
    let before = report(&split,&paused);
    assert!(parameters(&before).is_empty(),"calibration is not an observation");
    let frozen = before.get("mean_control").unwrap().get("calibration").unwrap().clone();
    std::fs::rename(split.dir.join("sources"),split.dir.join("relocated-sources")).unwrap();
    let resumed = report(&split,&study(&split,paused.str_field("run"),None,fs_cli::exit::SUCCESS));
    for key in ["observations", "statistics", "qmc", "mean_control"] {
        assert_eq!(resumed.get(key),expected.get(key),"{key}");
    }
    assert_eq!(resumed.get("mean_control").unwrap().get("calibration"),Some(&frozen));
}
