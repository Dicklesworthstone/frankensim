//! G1/G3/G4: a hot initial state must enter uncertainty and probability events
//! even when the same ordinary native trajectory ends below the thermal limit.
use super::*;
use fs_project::{ConductionTransient, TransientInitialTemperature, TransientRegionCapacity, spec::dims};
use fs_qty::QtyAny;

fn fixture() -> Fixture {
    let mut fixture = Fixture::new();
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 321.0;
    requirement.margin.value = 1.0;
    let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient = Some(ConductionTransient {
        initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(330.0, dims::TEMPERATURE)),
        horizon: QtyAny::new(2.0, dims::TIME),
        max_step: QtyAny::new(0.5, dims::TIME),
        max_steps: 12,
        energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
        capacities: setup.regions.iter().map(|row| TransientRegionCapacity {
            region: row.region.clone(),
            volumetric_heat_capacity: QtyAny::new(100.0, dims::VOLUMETRIC_HEAT_CAPACITY),
            source: "synthetic heat capacity for transient-peak probability regression".into(),
        }).collect(),
        enthalpy: None, power_schedules: Vec::new(),
    });
    std::fs::write(fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap()).unwrap();
    fixture
}

fn source() -> String {
    STUDY.replace("temperature-max", "temperature-peak")
        .replace(":low 4W :high 6W", ":low 5W :high 5W")
        .replace(":low 294K :high 300K", ":low 297K :high 297K")
}

#[test]
fn g1_g3_mc_qmc_and_lhs_observe_the_peak_of_the_same_native_child_solve() {
    for method in [
        "monte-carlo",
        "quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)",
        "latin-hypercube :lhs (replicated-latin-hypercube :replicates 2 :samples-per-replicate 2)",
    ] {
        let fixture = fixture();
        let peak_source = source().replace("monte-carlo", method);
        std::fs::write(fixture.sources.join("study.fsim"),
            peak_source.replace("temperature-peak", "temperature-max")).unwrap();
        let final_result = fixture.study(None, fs_cli::exit::SUCCESS);
        let (_, final_report) = fixture.report(&final_result);
        std::fs::write(fixture.sources.join("study.fsim"), peak_source).unwrap();
        let result = fixture.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = fixture.report(&result);
        assert_eq!(report.str_field("qoi"), Some("temperature-peak"));
        assert_eq!(report.str_field("temperature_observation"),
            Some("fine-grid-sampled-peak-including-initial"));
        assert!(report.str_field("no_claim").unwrap().contains("not a continuous-time"));
        for (peak, final_row) in rows(&report).iter().zip(rows(&final_report)) {
            assert_eq!(peak.f64_field("value_k"), Some(330.0));
            assert!(final_row.f64_field("value_k").unwrap() < 320.0);
            assert_eq!(peak.str_field("run"), final_row.str_field("run"));
            assert_eq!(peak.str_field("qoi_receipt"), final_row.str_field("qoi_receipt"));
            assert_eq!(peak.get("parameters"), final_row.get("parameters"));
        }
        let fraction = if method == "monte-carlo" {
            report.path(&["statistics", "empirical_probability_of_compliance"])
        } else {
            report.path(&[if method.starts_with("quasi") { "qmc" } else { "lhs" },
                "compliance_estimate", "probability_of_compliance"])
        };
        assert_eq!(fraction.and_then(JsonValue::as_f64), Some(0.0));
    }
}

#[test]
fn g1_peak_compliance_can_reject_a_trajectory_whose_final_temperature_passes() {
    let fixture = fixture();
    let source = source().replace(":version 1",
        ":version 2 :compliance (bernoulli-mixture :required-probability 0.99 :alpha 0.05 :min-samples 2)");
    std::fs::write(fixture.sources.join("study.fsim"), source).unwrap();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&result);
    assert_eq!(report.str_field("status"), Some("decision-reached"));
    assert_eq!(report.path(&["compliance", "decision"]).and_then(JsonValue::as_str),
        Some("below-probability-target"));
    assert_eq!(report.path(&["compliance", "empirical_probability_of_compliance"])
        .and_then(JsonValue::as_f64), Some(0.0));
    assert_eq!(rows(&report).len(), 2);
    assert_eq!(report.f64_field("temperature_limit_k"), Some(320.0));
}

#[test]
fn g4_g5_peak_recovery_uses_retained_trajectories_without_reopening_sources() {
    let full = fixture();
    std::fs::write(full.sources.join("study.fsim"), source()).unwrap();
    let completed = full.study(None, fs_cli::exit::SUCCESS);
    let expected = full.report(&completed).0;
    let split = fixture();
    std::fs::write(split.sources.join("study.fsim"), source()).unwrap();
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    assert_eq!(rows(&split.report(&prefix).1)[0].f64_field("value_k"), Some(330.0));
    std::fs::rename(&split.sources, split.dir.join("original-inputs-relocated")).unwrap();
    let result = split.resume(prefix.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(split.report(&result).0, expected);
}
