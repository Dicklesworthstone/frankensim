//! G1/G3/G4: a hot initial state must enter uncertainty and probability events
//! even when the same ordinary native trajectory ends below the thermal limit.
mod schedule;
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
        boundary_schedules: Vec::new(),
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

#[test]
fn g1_boundary_cycles_reach_native_probability_peaks_and_history_aware_mean_calibration() {
    let mut fixture = fixture();
    let card = fixture.project.materials.as_ref().unwrap()[0].card.clone();
    let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = fs_project::ThermalBoundaryCondition::FixedTemperature {
        temperature: QtyAny::new(300.0, dims::TEMPERATURE),
    };
    let time = setup.transient.as_mut().unwrap();
    time.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    time.capacities.clear();
    time.enthalpy = Some(fs_project::TransientEnthalpy {
        max_newton_iterations: 32, max_backtracks: 24,
        materials: vec![fs_project::TransientRegionEnthalpy {
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&card).unwrap(),
            source: "synthetic sensible chart for native boundary-cycle sampling".into(),
            reference_density: QtyAny::new(10.0, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            knots: [(-500.0, 243.15), (1000.0, 393.15)].into_iter()
                .map(|(h, temperature)| fs_project::TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                    temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                    liquid_mass_fraction: 0.0,
                }).collect(),
            phase_conductivity: None,
        }],
    });
    time.boundary_schedules = vec![fs_project::TransientBoundarySchedule {
        target: setup.boundaries[0].target.clone(),
        kind: fs_project::TransientBoundaryTarget::FixedTemperature,
        source: "manufactured thermal cycle for native sampled-peak consumers".into(),
        steps: [(0.5, 300.0), (1.0, 340.0), (2.0, 300.0)].into_iter()
            .map(|(until, value)| fs_project::TransientBoundaryStep {
                until: QtyAny::new(until, dims::TIME),
                value: QtyAny::new(value, dims::TEMPERATURE),
            }).collect(),
    }];
    std::fs::write(fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap()).unwrap();
    let source = STUDY.replace("temperature-max", "temperature-peak")
        .replace(":target convection-temperature :entity \"air\"", ":target initial-specific-enthalpy :entity \"air\"")
        .replace(":low 294K :high 300K", ":low -1.5J/kg :high 18.5J/kg");
    std::fs::write(fixture.sources.join("study.fsim"), &source).unwrap();
    let completed = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&completed);
    assert_eq!(rows(&report).len(), 4);
    for row in rows(&report) {
        assert!((row.f64_field("value_k").unwrap() - 340.0).abs() < 1e-10,
            "the native child must retain the hot intermediate boundary state");
    }
    let controlled = source.replace(":version 1",
        ":version 3 :mean-control (nominal-adjoint :max-solves 1)");
    std::fs::write(fixture.sources.join("study.fsim"), controlled).unwrap();
    let calibrated = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, calibrated_report) = fixture.report(&calibrated);
    assert_eq!(rows(&calibrated_report), rows(&report),
        "calibration must preserve every original boundary-driven probability solve");
    assert_eq!(calibrated_report.get("statistics"), report.get("statistics"));
    let control = calibrated_report.get("mean_control").unwrap();
    assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
    for coefficient in control.get("gradient").unwrap().as_array().unwrap() {
        assert!(coefficient.as_f64().unwrap().abs() < 1e-10,
            "a fully prescribed intermediate peak does not depend on prior heat or initial enthalpy");
    }
}
