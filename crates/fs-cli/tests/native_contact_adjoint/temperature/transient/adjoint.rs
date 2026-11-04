//! Complete native histories and independent physical probes of live R(T).
use super::*;
use fs_project::{TransientBoundarySchedule, TransientBoundaryStep, TransientBoundaryTarget};

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";

fn contact_pack(surface: bool) -> (Vec<u8>, String, String) {
    if !surface {
        return temperature_pack_from_knots(vec![(280.0, 0.02), (400.0, 1.22)]);
    }
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: vec![
            ([280.0, 0.8e6], 0.02),
            ([400.0, 0.8e6], 1.22),
            ([280.0, 1.2e6], 0.02),
            ([400.0, 1.2e6], 1.22),
        ],
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

fn time(f: &mut Fixture) -> &mut ConductionTransient {
    f.project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
}

fn solve(f: &Fixture, ordinal: usize, output: Option<&str>, surface: bool) -> (J, J) {
    let (source, card, ledger) = prepare_with_pack(f, ordinal, output, contact_pack(surface));
    let result = run(f, &source, &card, &ledger);
    assert_eq!(
        result.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        result.stdout,
        result.stderr
    );
    let result = J::parse(&result.stdout).unwrap();
    let run_receipt = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = run_receipt
        .get("stages")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("stage") == Some("conduction"))
        .unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn report(receipt: &J) -> &J {
    receipt.get("transient").unwrap().get("adjoint").unwrap()
}

fn qoi(receipt: &J, peak: bool) -> f64 {
    let transient = receipt.get("transient").unwrap();
    if peak {
        n(transient.get("sampled_peak").unwrap(), "value_k")
    } else {
        n(transient.get("temporal_error").unwrap(), "fine_final_k")
    }
}

fn agrees(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-5 * expected.abs().max(0.01),
        "complete contact history derivative {actual:e} versus independent resolve {expected:e}"
    );
}

fn regional(report: &J, region: &str) -> f64 {
    let row = report
        .get("regional_power")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("region") == Some(region))
        .unwrap();
    n(row, "dtemperature_duniform_added_watts_k_per_w")
}

fn unchanged(plain: &J, plain_field: &J, differentiated: &J, field: &J) {
    assert_eq!(plain_field.get("temperature"), field.get("temperature"));
    for key in ["energy", "interfaces", "solver"] {
        assert_eq!(
            plain.get(key),
            differentiated.get(key),
            "forward {key} changed"
        );
    }
    let plain = plain.get("transient").unwrap();
    let differentiated = differentiated.get("transient").unwrap();
    for key in [
        "coarse",
        "fine",
        "energy",
        "nonlinear",
        "temperature_contact",
        "temporal_error",
        "sampled_peak",
    ] {
        assert_eq!(
            plain.get(key),
            differentiated.get(key),
            "forward trajectory {key} changed"
        );
    }
    let derivative = differentiated.get("adjoint").unwrap();
    assert_eq!(n(derivative, "reconstructed_solid_endpoints"), 0.0);
    assert_eq!(n(derivative, "reconstruction_krylov_iterations"), 0.0);
    assert_eq!(n(derivative, "max_reconstruction_krylov_iterations"), 0.0);
    assert_eq!(
        n(derivative, "revalidated_solid_endpoints"),
        n(derivative, "state_index")
    );
    assert!(n(derivative, "maximum_dual_relative_residual") < 1e-12);
    let contact = derivative.get("temperature_contact").unwrap();
    assert_eq!(
        contact.str_field("method"),
        Some("complete-capacity-contact-tangent")
    );
    assert_eq!(
        contact.str_field("residual_target"),
        Some("original-accepted-forward-physical-tolerance")
    );
    assert_eq!(n(contact, "control_steps"), n(derivative, "state_index"));
    assert!(n(contact, "workspace_estimate_bytes") > 0.0);
}

#[test]
fn live_contact_capacity_adjoint_preserves_the_forward_field_and_matches_full_histories() {
    for (nonlinear, surface) in [(false, false), (true, false), (true, true)] {
        let mut f = fixture(nonlinear, 0.5);
        f.project.power.as_mut().unwrap()[0].watts.value = 24.0;
        f.project.power.as_mut().unwrap()[0].duty = 1.0;
        time(&mut f).power_schedules.clear();
        let (plain, plain_field) = solve(&f, 90, None, surface);
        let (receipt, field) = solve(&f, 91, Some(FINAL), surface);
        unchanged(&plain, &plain_field, &receipt, &field);
        let derivative = report(&receipt);
        assert_eq!(n(derivative, "value_k"), qoi(&plain, false));
        let epsilon = 0.002;
        for (coordinate, expected) in [
            regional(derivative, "hot"),
            n(derivative, "dtemperature_dcapacity_multiplier_k"),
            n(derivative, "dtemperature_duniform_initial_k_per_k"),
        ]
        .into_iter()
        .enumerate()
        {
            let original = f.project.clone();
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                f.project = original.clone();
                match coordinate {
                    0 => f.project.power.as_mut().unwrap()[0].watts.value += sign * epsilon,
                    1 => {
                        for row in &mut time(&mut f).capacities {
                            row.volumetric_heat_capacity.value *= 1.0 + sign * epsilon;
                        }
                    }
                    _ => {
                        time(&mut f).initial_temperature = TransientInitialTemperature::Uniform(
                            QtyAny::new(293.15 + sign * epsilon, dims::TEMPERATURE),
                        )
                    }
                }
                values.push(qoi(
                    &solve(&f, 92 + 2 * coordinate + side, None, surface).0,
                    false,
                ));
            }
            f.project = original;
            agrees(expected, (values[1] - values[0]) / (2.0 * epsilon));
        }
        let repeated = solve(&f, 99, Some(FINAL), surface);
        assert_eq!(report(&repeated.0), derivative);
        assert_eq!(repeated.1.get("temperature"), field.get("temperature"));
    }
}

#[test]
fn live_contact_capacity_prescribed_bath_lift_contains_the_complete_contact_transpose() {
    let mut f = fixture(true, 0.5);
    let (receipt, _) = solve(&f, 100, Some(FINAL), false);
    let bath = &report(&receipt)
        .get("prescribed_temperature")
        .unwrap()
        .get("boundaries")
        .unwrap()
        .as_array()
        .unwrap()[0];
    assert_eq!(bath.str_field("target"), Some("cold"));
    let actual = n(bath, "dtemperature_dprescribed_temperature_k_per_k");
    assert!(
        actual > 0.0,
        "the hot body must respond to the cold bath across its live joint"
    );
    let mut values = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let ThermalBoundaryCondition::FixedTemperature { temperature } = &mut f
            .project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .boundaries[0]
            .condition
        else {
            panic!("prescribed cold fixture");
        };
        temperature.value = 293.15 + sign * 0.002;
        values.push(qoi(&solve(&f, 101 + side, None, false).0, false));
    }
    agrees(actual, (values[1] - values[0]) / 0.004);
}

#[test]
fn live_contact_sampled_peak_adjoint_uses_scheduled_boundaries_and_zero_future_controls() {
    let mut f = fixture(true, 0.5);
    f.project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = ThermalBoundaryCondition::Convection {
        coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
    };
    let policy = time(&mut f);
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(300.0, dims::TEMPERATURE));
    policy.power_schedules[0].steps[0].watts.value = 80.0;
    policy.power_schedules[0].steps[1].watts.value = 0.0;
    policy.boundary_schedules = vec![TransientBoundarySchedule {
        target: "cold".into(),
        kind: TransientBoundaryTarget::ConvectionTemperature,
        source: "synthetic held reservoir history for contact trajectory derivatives".into(),
        steps: [(0.75, 290.15), (2.0, 295.15)]
            .into_iter()
            .map(|(until, value)| TransientBoundaryStep {
                until: QtyAny::new(until, dims::TIME),
                value: QtyAny::new(value, dims::TEMPERATURE),
            })
            .collect(),
    }];
    let (plain, plain_field) = solve(&f, 110, None, false);
    let (receipt, field) = solve(&f, 111, Some(PEAK), false);
    unchanged(&plain, &plain_field, &receipt, &field);
    let derivative = report(&receipt);
    assert_eq!(n(derivative, "value_k"), qoi(&plain, true));
    assert!(qoi(&plain, true) > qoi(&plain, false) + 1e-3);
    assert!(n(derivative, "time_s") > 0.0 && n(derivative, "time_s") < 2.0);
    let selected = n(derivative, "state_index") as usize;
    let intervals = derivative.get("intervals").unwrap().as_array().unwrap();
    for interval in &intervals[selected..] {
        assert_eq!(n(interval, "dtemperature_dsource_multiplier_k"), 0.0);
        for region in interval.get("regional_power").unwrap().as_array().unwrap() {
            assert_eq!(n(region, "dtemperature_ddelivered_watts_k_per_w"), 0.0);
        }
    }
    let pulse: f64 = intervals
        .iter()
        .filter(|row| n(row, "start_s") < 0.75)
        .map(|row| {
            let hot = row
                .get("regional_power")
                .unwrap()
                .as_array()
                .unwrap()
                .iter()
                .find(|row| row.str_field("region") == Some("hot"))
                .unwrap();
            n(hot, "dtemperature_ddelivered_watts_k_per_w")
        })
        .sum();
    let robin = &derivative.get("robin").unwrap().as_array().unwrap()[0];
    assert_eq!(
        robin.get("reference_history_held_fixed"),
        Some(&J::Bool(true))
    );
    assert_eq!(n(robin, "dtemperature_dreference_k_per_k"), 0.0);
    let original = f.project.clone();
    for (coordinate, expected) in [pulse, n(robin, "dtemperature_dlog_htc_k")]
        .into_iter()
        .enumerate()
    {
        let mut values = Vec::new();
        let epsilon: f64 = 0.002;
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            f.project = original.clone();
            if coordinate == 0 {
                time(&mut f).power_schedules[0].steps[0].watts.value += sign * epsilon;
            } else {
                let ThermalBoundaryCondition::Convection { coefficient, .. } = &mut f
                    .project
                    .cooling
                    .as_mut()
                    .unwrap()
                    .conduction
                    .as_mut()
                    .unwrap()
                    .boundaries[0]
                    .condition
                else {
                    panic!("prescribed convection fixture");
                };
                coefficient.value *= (sign * epsilon).exp();
            }
            values.push(qoi(
                &solve(&f, 112 + 2 * coordinate + side, None, false).0,
                true,
            ));
        }
        agrees(expected, (values[1] - values[0]) / (2.0 * epsilon));
    }
}

#[test]
fn live_contact_initial_peak_does_not_reconstruct_or_differentiate_future_steps() {
    let mut f = fixture(false, 0.5);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    let policy = time(&mut f);
    policy.power_schedules.clear();
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(320.0, dims::TEMPERATURE));
    let (plain, plain_field) = solve(&f, 120, None, false);
    let (receipt, field) = solve(&f, 121, Some(PEAK), false);
    unchanged(&plain, &plain_field, &receipt, &field);
    let derivative = report(&receipt);
    assert_eq!(n(derivative, "state_index"), 0.0);
    assert_eq!(n(derivative, "value_k"), 320.0);
    assert_eq!(n(derivative, "adjoint_krylov_iterations"), 0.0);
    assert_eq!(n(derivative, "dtemperature_duniform_initial_k_per_k"), 1.0);
    assert_eq!(n(derivative, "dtemperature_dcapacity_multiplier_k"), 0.0);
    assert_eq!(regional(derivative, "hot"), 0.0);
    assert_eq!(regional(derivative, "cold"), 0.0);
}

#[test]
fn live_contact_capacity_derivatives_refuse_source_corners_without_a_completed_report() {
    let mut f = fixture(false, 0.5);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    f.project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = ThermalBoundaryCondition::FixedTemperature {
        temperature: QtyAny::new(300.0, dims::TEMPERATURE),
    };
    let policy = time(&mut f);
    policy.power_schedules.clear();
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(300.0, dims::TEMPERATURE));
    let kink = || temperature_pack_from_knots(vec![(280.0, 0.1), (300.0, 0.2), (400.0, 0.8)]);
    for (ordinal, output) in [None, Some(FINAL)].into_iter().enumerate() {
        let (source, card, ledger) = prepare_with_pack(&f, 130 + ordinal, output, kink());
        let result = run(&f, &source, &card, &ledger);
        if output.is_none() {
            assert_eq!(
                result.exit_code,
                fs_cli::exit::SUCCESS,
                "{}{}",
                result.stdout,
                result.stderr
            );
        } else {
            assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
            let diagnostic = format!("{}{}", result.stdout, result.stderr);
            assert!(
                diagnostic.contains("derivative") || diagnostic.contains("slope"),
                "{diagnostic}"
            );
            assert!(
                !diagnostic.contains("run_receipt"),
                "refused derivative published a completed run"
            );
        }
    }
}
