//! G1/G3 native project -> imported volume -> real time steps -> report.

use super::*;
use fs_project::{ConductionTransient, TransientRegionCapacity};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

pub(super) fn declaration(spec: &mut ProjectSpec, capacity: f64, horizon: f64) {
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient = Some(ConductionTransient {
        initial_temperature: QtyAny::new(293.15, fs_project::spec::dims::TEMPERATURE),
        horizon: QtyAny::new(horizon, fs_project::spec::dims::TIME),
        max_step: QtyAny::new(0.5, fs_project::spec::dims::TIME),
        max_steps: 600,
        energy_tolerance: QtyAny::new(1e-6, fs_project::spec::dims::ENERGY),
        power_schedules: Vec::new(),
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

pub(super) fn driven_face_project(bytes: &[u8]) -> ProjectSpec {
    let mut spec = conduction_fixture_project(7, bytes);
    declaration(&mut spec, 15.0, 0.5);
    spec.power.as_mut().unwrap()[0].watts.value = 0.0;
    let binding = &mut spec.materials.as_mut().unwrap()[0];
    binding.temp_lo.value = 280.0;
    binding.temp_hi.value = 360.0;
    spec.envelope.as_mut().unwrap().ambient_lo.value = 300.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 300.0;
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
        parent: "enclosure".to_string(),
        name: "driven-face".to_string(),
        display: "Prescribed thermal reservoir".to_string(),
        expect_id: None,
    });
    spec.assignments.as_mut().unwrap().push(GeometryAssignment {
        artifact: "enclosure".to_string(),
        target: "driven-face".to_string(),
        length_unit: "m".to_string(),
        selector: MeshSelector::HalfSpace {
            normal: [1.0, 1.0, 1.0],
            offset: 1.0,
            side: HalfSpaceSide::AtLeast,
            tolerance: 0.0,
        },
        allow_overlap: true,
    });
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient.as_mut().unwrap().initial_temperature.value = 300.0;
    setup.adiabatic_remainder = true;
    setup.boundaries = vec![ThermalBoundary {
        target: "driven-face".to_string(),
        condition: ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(330.0, fs_project::spec::dims::TEMPERATURE),
        },
    }];
    spec
}

pub(super) fn transient_conductivity_cards(kind: &str) -> CardPackSet {
    let (property, interpolation) = match kind {
        "nonlinear" | "constant-curve" => (
            fs_matdb::PropertyValue::Curve {
                abscissa: "T".to_string(),
                abscissa_dims: fs_project::spec::dims::TEMPERATURE,
                knots: if kind == "nonlinear" {
                    vec![(280.0, 1.6), (360.0, 3.2)]
                } else {
                    vec![(280.0, 2.0), (360.0, 2.0)]
                },
                dims: CONDUCTIVITY_DIMS,
            },
            fs_matdb::InterpolationPolicy::LinearInside,
        ),
        "constant-scalar" => (
            fs_matdb::PropertyValue::Scalar {
                value: 2.0,
                dims: CONDUCTIVITY_DIMS,
            },
            fs_matdb::InterpolationPolicy::ConstantWithinValidity,
        ),
        _ => panic!("unknown conductivity fixture"),
    };
    CardPackSet::admit(vec![raw_pack(
        CardPackKind::Material,
        "fixtures/transient-conductivity.fsmcdpk",
        material_pack_bytes_with_property("AA6061", kind, 280.0, 360.0, property, interpolation),
    )])
    .unwrap()
}

#[test]
fn g1_native_transient_nonlinear_driven_face_matches_the_discrete_physical_balance() {
    let mut fields = Vec::new();
    for kind in [
        "nonlinear",
        "constant-curve",
        "constant-scalar",
        "nonlinear",
    ] {
        let bytes = tetra_stl();
        let mut spec = driven_face_project(&bytes);
        let cards = transient_conductivity_cards(kind);
        rebind_to(&mut spec, &cards);
        let decoded = decode(&spec);
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (run, qoi, conduction, solution) =
            run_conjugate_to_completion(&ledger, &decoded, &cards);
        let field = J::parse(&solution).unwrap();
        let mut values: Vec<_> = field
            .get("temperature")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_f64().unwrap())
            .collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(
            values.len(),
            4,
            "the analytic oracle requires the imported unit tetra"
        );

        // Unit tetra: each row-sum capacity is C=15/24=0.625 J/K.
        // Three driven vertices are 330 K. The remaining scalar equation is
        // C*(T-new - T-old) + dt*k((T-new+3*330)/4)/2*(T-new-330)=0.
        // This reference uses no production assembly, matrix or nonlinear solve.
        let mut expected = 300.0;
        for _ in 0..2 {
            let (mut low, mut high) = (expected, 330.0);
            for _ in 0..80 {
                let value = 0.5 * (low + high);
                let k = if kind == "nonlinear" {
                    0.02 * ((value + 3.0 * 330.0) / 4.0) - 4.0
                } else {
                    2.0
                };
                let residual = 0.625 * (value - expected) + 0.25 * 0.5 * k * (value - 330.0);
                if residual > 0.0 {
                    high = value;
                } else {
                    low = value;
                }
            }
            expected = 0.5 * (low + high);
        }
        close(values[0], expected, 2e-7);
        for &fixed in &values[1..] {
            close(fixed, 330.0, 0.0);
        }
        close(receipt_number_field(&qoi, "value"), 330.0, 0.0);
        let receipt = J::parse(&conduction).unwrap();
        let time = receipt.get("transient").unwrap();
        close(n(time, "total_steps"), 3.0, 0.0);
        let energy = time.get("energy").unwrap();
        close(
            n(energy, "stored_change_j"),
            0.625 * (expected - 300.0 + 90.0),
            2e-7,
        );
        assert!(n(energy, "maximum_step_residual_j") <= 1e-6);
        if kind == "nonlinear" {
            let controls = time.get("nonlinear").unwrap();
            assert_eq!(controls.str_field("method"), Some("newton-fgmres"));
            assert!(n(controls, "total_updates") > 3.0);
            let mut updates = 0.0;
            for grid in ["coarse", "fine"] {
                for step in time.get(grid).unwrap().as_array().unwrap() {
                    let nonlinear = step.get("nonlinear").unwrap();
                    assert!(n(nonlinear, "residual_j") <= n(nonlinear, "threshold_j"));
                    assert!(n(nonlinear, "updates") <= n(controls, "max_updates_per_step"));
                    assert!(
                        n(step, "krylov_iterations")
                            <= n(controls, "max_krylov_iterations_per_step")
                    );
                    updates += n(nonlinear, "updates");
                }
            }
            close(n(controls, "total_updates"), updates, 0.0);
            assert!(conduction.contains("fgmres-backward-euler-newton"));
        } else {
            assert!(matches!(time.get("nonlinear"), Some(J::Null)));
            assert!(conduction.contains("pcg-backward-euler"));
        }
        fields.push((run, conduction, solution, values));
    }
    assert_eq!(
        fields[0], fields[3],
        "nonlinear state and work replay exactly"
    );
    assert!(
        (fields[0].3[0] - fields[1].3[0]).abs() > 1.0,
        "freezing conductivity at the initial temperature changes the physical answer"
    );
    for (&curve, &scalar) in fields[1].3.iter().zip(&fields[2].3) {
        close(curve, scalar, 2e-10);
    }
}

#[test]
fn g0_native_transient_nonlinear_domain_failure_publishes_no_endpoint() {
    let bytes = tetra_stl();
    let mut spec = driven_face_project(&bytes);
    spec.power.as_mut().unwrap()[0].watts.value = 1e6;
    let cards = transient_conductivity_cards("nonlinear");
    rebind_to(&mut spec, &cards);
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes);
    let error = run_solve(
        &ledger,
        &CancelGate::new_clock_free(),
        &mut benign_clock(),
        &decode(&spec),
        &cards,
        &mut Vec::new(),
    )
    .unwrap_err();
    assert_eq!(error.code, "cli-solve-conduction-transient", "{error:?}");
    assert_eq!(error.stage, Some("conduction"));
    assert_eq!(
        stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
        4
    );
}

#[test]
fn g0_native_transient_refuses_steady_only_combinations_before_a_field_is_published() {
    for variant in ["natural-convection", "ladder", "adjoint"] {
        let bytes = tetra_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        declaration(&mut spec, 1000.0, 2.0);
        if variant == "natural-convection" {
            spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0]
                .condition = ThermalBoundaryCondition::NaturalConvection {
                characteristic_length: QtyAny::new(0.1, fs_project::spec::dims::LENGTH),
                ambient_temperature: QtyAny::new(293.15, fs_project::spec::dims::TEMPERATURE),
                correlation: "convection.churchill-chu-vertical-plate".to_string(),
            };
        }
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
