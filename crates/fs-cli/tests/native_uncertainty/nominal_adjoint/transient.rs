//! G1/G3/G4/G5: actual native transient solves, frozen mean controls and recovery.
use super::*;
use fs_project::{
    ConductionTransient, TransientInitialTemperature, TransientPowerStep, TransientRegionCapacity,
    TransientRegionPower, spec::dims,
};
use fs_qty::QtyAny;

fn fixture() -> Fixture {
    let mut fixture = Fixture::new();
    fixture.project.power.as_mut().unwrap()[0].duty = 0.37;
    let setup = fixture
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    setup.transient = Some(ConductionTransient {
        initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(
            293.15,
            dims::TEMPERATURE,
        )),
        horizon: QtyAny::new(2.0, dims::TIME),
        max_step: QtyAny::new(0.5, dims::TIME),
        max_steps: 12,
        energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
        capacities: setup
            .regions
            .iter()
            .map(|row| TransientRegionCapacity {
                region: row.region.clone(),
                volumetric_heat_capacity: QtyAny::new(100.0, dims::VOLUMETRIC_HEAT_CAPACITY),
                source: "synthetic capacity for native transient mean-control integration".into(),
            })
            .collect(),
        enthalpy: None,
        power_schedules: Vec::new(),
    });
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    write_project(&fixture);
    fixture
}

fn write_project(fixture: &Fixture) {
    std::fs::write(
        fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap(),
    )
    .unwrap();
}

fn study(fixture: &Fixture, controlled: bool, qmc: bool) -> String {
    let mut text = study_text(false, qmc, controlled);
    // Include pressure alongside temperature to exercise the same projection
    // on the report-only calibration project and every raw trajectory.
    let mut project = fixture.project.clone();
    project.outputs.as_mut().unwrap().insert(
        0,
        fs_project::spec::OutputRequest {
            name: "pressure-drop".into(),
            kind: "scalar".into(),
            region: None,
        },
    );
    std::fs::write(
        fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&project).unwrap(),
    )
    .unwrap();
    // The same declared sampling law is used with and without calibration.
    text = text.replace(":wall-time 120s", ":wall-time 240s");
    std::fs::write(fixture.sources.join("study.fsim"), &text).unwrap();
    text
}

fn gradient(report: &JsonValue) -> Vec<f64> {
    report
        .get("mean_control")
        .unwrap()
        .get("gradient")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|row| row.as_f64().unwrap())
        .collect()
}

#[test]
fn g1_g3_transient_calibration_uses_final_history_and_keeps_raw_probability_samples() {
    for qmc in [false, true] {
        let raw = fixture();
        study(&raw, false, qmc);
        let raw_result = raw.study(None, fs_cli::exit::SUCCESS);
        let (_, raw_report) = raw.report(&raw_result);
        let controlled = fixture();
        study(&controlled, true, qmc);
        let result = controlled.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = controlled.report(&result);
        assert_eq!(
            rows(&report),
            rows(&raw_report),
            "calibration cannot change sampled projects or final fields"
        );
        assert_eq!(report.get("statistics"), raw_report.get("statistics"));
        assert_eq!(report.get("qmc"), raw_report.get("qmc"));
        let control = report.get("mean_control").unwrap();
        assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
        assert!(
            control
                .get("estimate")
                .unwrap()
                .f64_field("variance_ratio")
                .unwrap()
                < 1e-4,
            "fixed-grid linear transient response is affine in these physical controls: {control:?}"
        );

        let coefficients = gradient(&report);
        assert!(coefficients.iter().all(|value| *value > 0.0));
        let h = 0.01;
        for coordinate in 0..2 {
            let mut low = [5.0, 297.0];
            let mut high = low;
            low[coordinate] -= h;
            high[coordinate] += h;
            let a = controlled.independent_solve(100 + coordinate * 2, &low).2;
            let b = controlled.independent_solve(101 + coordinate * 2, &high).2;
            let difference = (b - a) / (2.0 * h);
            assert!(
                (coefficients[coordinate] - difference).abs() <= 1e-5 * difference.abs().max(1.0),
                "trajectory coordinate {coordinate}: adjoint={} FD={difference}",
                coefficients[coordinate]
            );
        }
    }
}

#[test]
fn g4_g5_transient_calibration_resumes_from_sealed_trajectory_without_source_paths() {
    let full = fixture();
    study(&full, true, true);
    let result = full.study(None, fs_cli::exit::SUCCESS);
    let (expected, _) = full.report(&result);
    let split = fixture();
    study(&split, true, true);
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    let (_, report) = split.report(&prefix);
    assert!(
        rows(&report).is_empty(),
        "the adjoint trajectory is calibration, not a probability sample"
    );
    let frozen = report
        .get("mean_control")
        .unwrap()
        .get("calibration")
        .unwrap()
        .clone();
    std::fs::rename(&split.sources, split.dir.join("retained-inputs-only")).unwrap();
    let one = split.resume(
        prefix.str_field("run").unwrap(),
        Some("1"),
        fs_cli::exit::BUDGET,
    );
    let (_, report) = split.report(&one);
    assert_eq!(rows(&report).len(), 1);
    assert_eq!(
        report.get("mean_control").unwrap().get("calibration"),
        Some(&frozen)
    );
    let finished = split.resume(one.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(
        split.report(&finished).0,
        expected,
        "restored controls and sampling reproduce the uninterrupted report bytes"
    );
}

#[test]
fn g3_scheduled_watts_override_the_random_static_power_coordinate() {
    let mut fixture = fixture();
    let time = fixture
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap();
    time.power_schedules.push(TransientRegionPower {
        region: "air".into(),
        source: "absolute delivered-power history overrides static power".into(),
        steps: vec![TransientPowerStep {
            until: QtyAny::new(2.0, dims::TIME),
            watts: QtyAny::new(2.0, dims::POWER),
        }],
    });
    study(&fixture, true, false);
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&result);
    assert_eq!(
        gradient(&report)[0],
        0.0,
        "varying a fully overridden declaration has no thermal effect"
    );
    let low = fixture.independent_solve(110, &[4.0, 297.0]).2;
    let high = fixture.independent_solve(111, &[6.0, 297.0]).2;
    assert_eq!(low.to_bits(), high.to_bits());
}

fn enthalpy_fixture() -> Fixture {
    use fs_project::{TransientEnthalpy, TransientEnthalpyKnot, TransientRegionEnthalpy};
    let mut fixture = fixture();
    let card = fixture.project.materials.as_ref().unwrap()[0].card.clone();
    let time = fixture
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap();
    time.capacities.clear();
    time.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&card).unwrap(),
            source: "synthetic sensible enthalpy chart for native UQ integration".into(),
            reference_density: QtyAny::new(10.0, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            knots: [(-500.0, 243.15), (1000.0, 393.15)]
                .into_iter()
                .map(|(h, t)| TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                    temperature: QtyAny::new(t, dims::TEMPERATURE),
                    liquid_mass_fraction: 0.0,
                })
                .collect(),
            phase_conductivity: None,
        }],
    });
    write_project(&fixture);
    fixture
}

#[test]
fn g1_enthalpy_power_mean_control_consumes_the_existing_enthalpy_adjoint() {
    let fixture = enthalpy_fixture();
    // A constant reference is a zero centered input and needs no calibration
    // coefficient even though the enthalpy producer supplies its derivative.
    let text = study_text(false, false, true).replace(":high 300K", ":high 294K");
    std::fs::write(fixture.sources.join("study.fsim"), text).unwrap();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&result);
    let coefficients = gradient(&report);
    assert!(coefficients[0] > 0.0);
    assert_eq!(coefficients[1], 0.0);
    let h = 0.01;
    let low = fixture.independent_solve(120, &[5.0 - h, 294.0]).2;
    let high = fixture.independent_solve(121, &[5.0 + h, 294.0]).2;
    let difference = (high - low) / (2.0 * h);
    assert!((coefficients[0] - difference).abs() <= 1e-5 * difference.abs().max(1.0));
}

#[test]
fn g1_g5_enthalpy_bath_mean_control_matches_physics_and_survives_resume() {
    let bath_fixture = |controlled| {
        let mut fixture = enthalpy_fixture();
        fixture
            .project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .boundaries[0]
            .condition = fs_project::ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(297.0, dims::TEMPERATURE),
        };
        write_project(&fixture);
        let text = study_text(false, false, controlled)
            .replace(":target convection-temperature", ":target fixed-temperature")
            .replace(":wall-time 120s", ":wall-time 240s");
        std::fs::write(fixture.sources.join("study.fsim"), text).unwrap();
        fixture
    };

    let raw = bath_fixture(false);
    let result = raw.study(None, fs_cli::exit::SUCCESS);
    let (_, raw_report) = raw.report(&result);
    let full = bath_fixture(true);
    let result = full.study(None, fs_cli::exit::SUCCESS);
    let (expected, report) = full.report(&result);
    let coefficients = gradient(&report);
    assert_eq!(coefficients[0], 0.0, "fixed nodes reject power influence");
    assert!(
        (coefficients[1] - 1.0).abs() < 1e-12,
        "a fully prescribed body's final temperature follows its bath: {coefficients:?}"
    );
    assert_eq!(
        rows(&report),
        rows(&raw_report),
        "bath calibration must preserve every raw probability trajectory"
    );
    assert_eq!(report.get("statistics"), raw_report.get("statistics"));
    for row in rows(&report) {
        let bath = row.get("parameters").unwrap().as_array().unwrap()[1]
            .as_f64()
            .unwrap();
        assert!((row.f64_field("value_k").unwrap() - bath).abs() < 1e-10);
    }

    // Bind the ordinary native project directly, independently of the
    // uncertainty adapter and its retained derivative report.
    let h = 0.01;
    let mut temperatures = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let mut project = full.project.clone();
        let fs_project::ThermalBoundaryCondition::FixedTemperature { temperature } =
            &mut project
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .boundaries[0]
                .condition
        else {
            panic!("prescribed enthalpy bath fixture");
        };
        temperature.value += sign * h;
        temperatures.push(full.solve_project(130 + side, &project).2);
    }
    let difference = (temperatures[1] - temperatures[0]) / (2.0 * h);
    assert!((coefficients[1] - difference).abs() < 1e-8);

    let split = bath_fixture(true);
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    let (_, prefix_report) = split.report(&prefix);
    assert!(rows(&prefix_report).is_empty());
    assert_eq!(gradient(&prefix_report), coefficients);
    let frozen = prefix_report
        .get("mean_control")
        .unwrap()
        .get("calibration")
        .unwrap()
        .clone();
    std::fs::rename(&split.sources, split.dir.join("retained-bath-inputs")).unwrap();
    let finished = split.resume(
        prefix.str_field("run").unwrap(),
        None,
        fs_cli::exit::SUCCESS,
    );
    let (actual, restored) = split.report(&finished);
    assert_eq!(
        restored.get("mean_control").unwrap().get("calibration"),
        Some(&frozen)
    );
    assert_eq!(
        actual, expected,
        "sealed bath calibration resumes to the uninterrupted native report"
    );
}

#[path = "transient/convection.rs"]
mod convection;


#[path = "transient/peak.rs"]
mod peak;
