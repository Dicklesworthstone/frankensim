//! The existing native study executor consumes both physical boundary controls.
use super::*;

fn fixture(natural: bool) -> Fixture {
    let mut f = Fixture::new();
    f.project = project(natural);
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
    std::fs::copy(reference.join("gray-surface.fsmcdpk"), f.sources.join("gray-surface.fsmcdpk")).unwrap();
    std::fs::write(f.sources.join("cooling-reference.fsim"), fs_project::print_sexpr(&f.project).unwrap()).unwrap();
    f
}

fn source(f: &Fixture, natural: bool, copula: bool, qmc: bool, controlled: bool) -> String {
    let surface = &f.project.cooling.as_ref().unwrap().conduction.as_ref().unwrap()
        .radiation.as_ref().unwrap().surfaces[0].name;
    let ambient = if natural { "natural-convection-ambient" } else { "convection-temperature" };
    let text = study_text(copula, qmc, controlled);
    let power = "(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)";
    let old_ambient = "(uniform :name \"ambient\" :target convection-temperature :entity \"air\" :low 294K :high 300K)";
    assert!(text.contains(power) && text.contains(old_ambient));
    text.replace(power, &format!("(uniform :name \"fluid\" :target {ambient} :entity \"air\" :low 294K :high 300K)"))
        .replace(old_ambient, &format!("(uniform :name \"reservoir\" :target radiation-reservoir-temperature :entity \"{surface}\" :low 305K :high 315K)"))
        .replace(":materials (\"aa6061.fsmcdpk\")", ":materials (\"aa6061.fsmcdpk\" \"gray-surface.fsmcdpk\")")
}

fn physical_project(f: &Fixture, fluid: f64, reservoir: f64) -> fs_project::ProjectSpec {
    let mut project = f.project.clone();
    let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    match &mut setup.boundaries[0].condition {
        fs_project::ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } => ambient_temperature.value = fluid,
        fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. } => reference_temperature.value = fluid,
        _ => panic!("fixture boundary"),
    }
    setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value = reservoir;
    project
}

#[test]
fn native_boundary_temperature_calibration_keeps_all_four_sampling_laws_unchanged() {
    for natural in [false, true] {
        for (copula, qmc) in [(false, false), (true, false), (false, true), (true, true)] {
            let raw = fixture(natural);
            std::fs::write(raw.sources.join("study.fsim"), source(&raw, natural, copula, qmc, false)).unwrap();
            let result = raw.study(None, fs_cli::exit::SUCCESS);
            let (_, baseline) = raw.report(&result);
            let f = fixture(natural);
            std::fs::write(f.sources.join("study.fsim"), source(&f, natural, copula, qmc, true)).unwrap();
            let result = f.study(None, fs_cli::exit::SUCCESS);
            let (_, report) = f.report(&result);
            assert_eq!(rows(&report).len(), 4);
            assert_eq!(rows(&report), rows(&baseline), "calibration must not move a probability child");
            for key in ["statistics", "qmc"] { assert_eq!(report.get(key), baseline.get(key)); }
            let control = report.get("mean_control").unwrap();
            assert_eq!(control.str_field("status"), Some("frozen"));
            assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
            let gradient = control.get("gradient").unwrap().as_array().unwrap();
            assert_eq!(gradient.len(), 2);
            assert!(gradient.iter().all(|v| v.as_f64().unwrap() > 0.0), "both independent reservoirs affect this model");
            if !copula && !qmc {
                // Bypass parameter application and hand-bind each actual native
                // project, rather than comparing two callers of the same map.
                for (i, sample) in rows(&report).iter().enumerate() {
                    let p = sample.get("parameters").unwrap().as_array().unwrap();
                    let project = physical_project(&f, p[0].as_f64().unwrap(), p[1].as_f64().unwrap());
                    let solved = run(&f, &project, 900+i);
                    let actual = solved.get("temperature").unwrap().f64_field("max").unwrap();
                    assert_eq!(actual.to_bits(), sample.f64_field("value_k").unwrap().to_bits());
                }
                let mut nominal = physical_project(&f, 297.0, 310.0);
                request(&mut nominal);
                let nominal = run(&f, &nominal, 910);
                let ambient = if natural { "natural-convection-ambient" } else { "convection-temperature" };
                for (i, target) in [ambient, "radiation-reservoir-temperature"].into_iter().enumerate() {
                    assert_eq!(gradient[i].as_f64().unwrap().to_bits(), coefficient(&nominal, target).to_bits());
                }
            }
        }
    }
}

#[test]
fn native_boundary_temperature_resume_retains_calibration_after_sources_move() {
    let full = fixture(true);
    let text = source(&full, true, true, true, true);
    std::fs::write(full.sources.join("study.fsim"), &text).unwrap();
    let complete = full.study(None, fs_cli::exit::SUCCESS);
    let (expected, _) = full.report(&complete);
    let split = fixture(true);
    std::fs::write(split.sources.join("study.fsim"), text).unwrap();
    let partial = split.study(Some("1"), fs_cli::exit::BUDGET);
    let (_, before) = split.report(&partial);
    assert!(rows(&before).is_empty(), "calibration is not a random observation");
    let frozen = before.get("mean_control").unwrap().get("calibration").unwrap().clone();
    std::fs::rename(&split.sources, split.dir.join("relocated-boundary-sources")).unwrap();
    let result = split.resume(partial.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let (actual, report) = split.report(&result);
    assert_eq!(actual, expected);
    assert_eq!(report.get("mean_control").unwrap().get("calibration"), Some(&frozen));
}
