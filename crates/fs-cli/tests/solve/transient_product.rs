//! G1/G3 native project -> imported volume -> real time steps -> report.

use super::*;
use fs_project::{ConductionTransient, TransientRegionCapacity};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

fn declaration(spec: &mut ProjectSpec, capacity: f64, horizon: f64) {
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient = Some(ConductionTransient {
        initial_temperature: QtyAny::new(293.15, fs_project::spec::dims::TEMPERATURE),
        horizon: QtyAny::new(horizon, fs_project::spec::dims::TIME),
        max_step: QtyAny::new(0.5, fs_project::spec::dims::TIME),
        max_steps: 600,
        energy_tolerance: QtyAny::new(1e-6, fs_project::spec::dims::ENERGY),
        capacities: setup
            .regions
            .iter()
            .map(|row| TransientRegionCapacity {
                region: row.region.clone(),
                volumetric_heat_capacity: QtyAny::new(
                    capacity,
                    fs_project::spec::dims::VOLUMETRIC_HEAT_CAPACITY,
                ),
                source: "synthetic constant capacity for the G1 thermal-storage oracle".to_string(),
            })
            .collect(),
    });
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    let envelope = spec.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 293.15;
    envelope.ambient_hi.value = 293.15;
}

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

#[test]
fn g1_native_transient_shipped_example_completes_with_retained_time_evidence() {
    let decoded = fs_project::parse_sexpr(include_str!(
        "../../../../data/reference-project/cooling-transient.fsim"
    ))
    .unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(
        &ledger,
        &decoded.spec,
        include_bytes!("../../../../data/reference-project/plate.stl").to_vec(),
    );
    let (_, qoi, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &fixture_cards());
    let parsed = J::parse(&conduction).unwrap();
    let time = parsed.get("transient").unwrap();
    close(n(time, "total_steps"), 12.0, 0.0);
    close(n(time, "final_time_s"), 2.0, 0.0);
    assert!(receipt_number_field(&qoi, "value") > 293.15);
    assert_eq!(
        time.get("temporal_error").unwrap().str_field("status"),
        Some("estimated-order-one")
    );
}

#[test]
fn g1_native_transient_uniform_heating_respects_storage_time_and_capacity() {
    let mut runs = Vec::new();
    for (capacity, horizon) in [(1000.0, 2.0), (1000.0, 4.0), (2000.0, 4.0), (1000.0, 2.0)] {
        let bytes = tetra_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        declaration(&mut spec, capacity, horizon);
        spec.cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .boundaries[0]
            .condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, fs_project::spec::dims::HEAT_FLUX),
        };
        let canonical = print_sexpr(&spec).unwrap();
        assert!(canonical.contains(":transient (transient"));
        let decoded = fs_project::parse_sexpr(&canonical).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (run, qoi, conduction, solution) =
            run_conjugate_to_completion(&ledger, &decoded, &fixture_cards());
        let receipt = J::parse(&conduction).unwrap();
        let transient = receipt
            .get("transient")
            .expect("native time evidence retained");
        assert_eq!(transient.str_field("status"), Some("completed"));
        close(n(transient, "final_time_s"), horizon, 0.0);
        close(
            n(transient, "total_steps"),
            3.0 * (horizon / 0.5).ceil(),
            0.0,
        );
        let expected = 293.15 + 5.0 * horizon / (capacity / 6.0);
        close(receipt_number_field(&qoi, "value"), expected, 2e-10);
        let field = J::parse(&solution).unwrap();
        for value in field.get("temperature").unwrap().as_array().unwrap() {
            close(value.as_f64().unwrap(), expected, 2e-10);
        }
        let energy = transient.get("energy").unwrap();
        close(n(energy, "stored_change_j"), 5.0 * horizon, 2e-7);
        close(n(energy, "integrated_net_input_j"), 5.0 * horizon, 2e-7);
        assert!(n(energy, "maximum_step_residual_j") <= 1e-6);
        let endpoint_energy = receipt.get("energy").unwrap();
        close(n(endpoint_energy, "storage_w"), 5.0, 2e-7);
        close(n(endpoint_energy, "closure_w"), 0.0, 2e-7);
        close(
            n(endpoint_energy, "source_w") + n(endpoint_energy, "dirichlet_in_w")
                - n(endpoint_energy, "neumann_out_w")
                - n(endpoint_energy, "robin_out_w")
                - n(endpoint_energy, "storage_w"),
            n(endpoint_energy, "closure_w"),
            2e-12,
        );
        let qoi_record = J::parse(&qoi).unwrap();
        let time_scope = qoi_record.get("thermal_time_scope").unwrap();
        assert_eq!(time_scope.str_field("kind"), Some("final"));
        close(n(time_scope, "time_s"), horizon, 0.0);
        let terms = J::parse(&qoi)
            .unwrap()
            .get("budget")
            .unwrap()
            .as_array()
            .unwrap()[0]
            .get("terms")
            .unwrap()
            .as_array()
            .unwrap()
            .clone();
        for kind in ["discretization", "solver-algebraic", "roundoff"] {
            assert!(
                terms.iter().any(|term| term.str_field("kind") == Some(kind)
                    && term.str_field("state") == Some("no-data")),
                "steady or temporal-only evidence must not certify {kind}: {qoi}"
            );
        }
        let report_receipts = stage_receipt_hashes(&ledger, &run);
        let report = String::from_utf8(artifact_bytes(&ledger, &report_receipts[6])).unwrap();
        let html = String::from_utf8(artifact_bytes(
            &ledger,
            &receipt_str_field(&report, "report_html"),
        ))
        .unwrap();
        assert!(
            html.contains("backward") || html.contains("transient"),
            "{html}"
        );
        runs.push((run, qoi, conduction, solution));
    }
    assert_eq!(
        runs[0], runs[3],
        "same physical inputs reproduce run and retained bytes"
    );
    assert_ne!(runs[1].0, runs[2].0, "capacity is part of the run identity");
    close(
        receipt_number_field(&runs[0].1, "value"),
        receipt_number_field(&runs[2].1, "value"),
        2e-10,
    );
}

#[test]
fn g1_native_transient_robin_warming_reports_actual_nested_time_difference() {
    let bytes = tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    declaration(&mut spec, 100.0, 2.0);
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes.clone());
    let (_, qoi, conduction, _) =
        run_conjugate_to_completion(&ledger, &decode(&spec), &fixture_cards());
    let parsed = J::parse(&conduction).unwrap();
    let transient = parsed.get("transient").unwrap();
    let time = transient.get("temporal_error").unwrap();
    assert_eq!(time.str_field("status"), Some("estimated-order-one"));
    let difference = (n(time, "fine_final_k") - n(time, "coarse_final_k")).abs();
    assert!(
        difference > 1e-6,
        "the time comparison must resolve real truncation error"
    );
    close(n(time, "absolute_difference_k"), difference, 1e-13);
    close(n(time, "estimated_half_width_k"), 1.25 * difference, 1e-13);
    close(
        receipt_number_field(&qoi, "value"),
        n(time, "fine_final_k"),
        0.0,
    );
    let steady = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    steady.transient = None;
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes);
    let (_, steady_qoi, _, _) =
        run_conjugate_to_completion(&ledger, &decode(&spec), &fixture_cards());
    let final_transient = receipt_number_field(&qoi, "value");
    assert!(final_transient > 293.15);
    assert!(final_transient < receipt_number_field(&steady_qoi, "value"));
}

#[test]
fn g0_native_transient_refuses_steady_only_combinations_before_a_field_is_published() {
    for variant in ["airflow", "ladder", "adjoint"] {
        let bytes = tetra_stl();
        let mut spec = if variant == "airflow" {
            conjugate_fixture_project(7, &bytes, "convection.gnielinski")
        } else {
            conduction_fixture_project(7, &bytes)
        };
        declaration(&mut spec, 1000.0, 2.0);
        if variant == "ladder" {
            spec.solver.as_mut().unwrap().fidelity = "ladder".to_string();
        }
        if variant == "adjoint" {
            spec.outputs.as_mut().unwrap().push(OutputRequest {
                name: "temperature-max-adjoint".to_string(),
                kind: "report".to_string(),
                region: None,
            });
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let error = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decode(&spec),
            &fixture_cards(),
            &mut Vec::new(),
        )
        .unwrap_err();
        assert_eq!(
            error.code, "cli-solve-conduction-transient",
            "{variant}: {error:?}"
        );
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4
        );
    }
}
