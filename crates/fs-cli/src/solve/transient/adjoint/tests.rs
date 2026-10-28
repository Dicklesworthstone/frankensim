//! Tests of the native time producer after parsing a real .fsim project.
//! Geometry is a manufactured unit cube supplied at the existing lowering seam.
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::box_grid;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_project::{OutputRequest, TransientPowerStep, TransientRegionPower};
use fs_qty::QtyAny;
use crate::json_read::JsonValue as J;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 71, kernel_id: 825, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}

#[derive(Clone, Copy)]
struct Inputs { power: f64, capacity: f64, initial: f64, htc: f64, pulse: bool, nonlinear: bool }
impl Default for Inputs {
    fn default() -> Self {
        Self { power: 5.0, capacity: 1000.0, initial: 293.15, htc: 0.0, pulse: false, nonlinear: false }
    }
}

fn project(input: Inputs, output: Option<&str>) -> ProjectSpec {
    let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-transient.fsim"))).unwrap().decoded.spec;
    spec.power.as_mut().unwrap()[0].watts.value = input.power;
    spec.power.as_mut().unwrap()[0].duty = 1.0;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = if input.htc == 0.0 {
        ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, fs_project::spec::dims::HEAT_FLUX),
        }
    } else {
        ThermalBoundaryCondition::Convection {
            coefficient: QtyAny::new(input.htc, fs_project::spec::dims::HEAT_TRANSFER_COEFFICIENT),
            reference_temperature: QtyAny::new(293.15, fs_project::spec::dims::TEMPERATURE),
        }
    };
    let time = setup.transient.as_mut().unwrap();
    time.initial_temperature = fs_project::TransientInitialTemperature::Uniform(
        QtyAny::new(input.initial, fs_project::spec::dims::TEMPERATURE));
    time.horizon.value = 4.0;
    time.max_step.value = 0.5;
    time.max_steps = 1000;
    time.energy_tolerance.value = 1e-6;
    for row in &mut time.capacities { row.volumetric_heat_capacity.value = input.capacity; }
    if input.pulse {
        time.power_schedules = vec![TransientRegionPower {
            region: "solid".into(), source: "manufactured pulse".into(),
            steps: vec![
                TransientPowerStep {
                    until: QtyAny::new(1.0, fs_project::spec::dims::TIME),
                    watts: QtyAny::new(input.power, fs_project::spec::dims::POWER),
                },
                TransientPowerStep {
                    until: QtyAny::new(4.0, fs_project::spec::dims::TIME),
                    watts: QtyAny::new(0.0, fs_project::spec::dims::POWER),
                },
            ],
        }];
    }
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let canonical = fs_project::print_sexpr(&spec).unwrap();
    let decoded = fs_project::parse_sexpr(&canonical).unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    decoded.spec
}

fn run(cx: &Cx<'_>, input: Inputs, output: Option<&str>) -> (String, Vec<f64>) {
    let spec = project(input, output);
    let (complex, positions) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = if input.htc == 0.0 {
        ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap()
    } else {
        ThermalBoundaryBuilder::new(&mesh).region("solid", |_| true,
            ThermalBc::robin(input.htc, 293.15).unwrap()).unwrap().finish().unwrap()
    };
    let material = if input.nonlinear {
        ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280.0, 5.0), (340.0, 20.0)]).unwrap())
    } else { ConductivityModel::isotropic_declared(10.0).unwrap() };
    let source = ScalarField::Uniform(input.power);
    let ids = std::collections::BTreeMap::from([("solid".to_string(), 1)]);
    let labels = vec![1; mesh.element_count()];
    let result = super::super::solve(cx, &spec, ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source,
    }, None, &labels, &ids, None, None, &[], 101325.0,
        LinearConfig { tolerance: 1e-10, max_iterations: 2000, restart: 20 }, None).unwrap();
    (result.receipt, result.solution.conduction.temperature)
}
fn number(row: &J, key: &str) -> f64 { row.f64_field(key).unwrap() }
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:e} vs {expected:e}");
}
fn value(receipt: &str) -> f64 {
    number(J::parse(receipt).unwrap().get("adjoint").unwrap(), "value_k")
}

#[test]
fn native_transient_adjoint_uniform_heating_and_zero_baseline_are_analytic() {
    with_cx(|cx| for power in [0.0, 5.0] {
        let input = Inputs { power, ..Inputs::default() };
        let (receipt, field) = run(cx, input, Some("temperature-final-adjoint"));
        let parsed = J::parse(&receipt).unwrap();
        let g = parsed.get("adjoint").unwrap();
        near(number(g, "value_k"), input.initial + power * 4.0 / 1000.0, 1e-8);
        near(number(g, "dtemperature_duniform_initial_k_per_k"), 1.0, 1e-8);
        near(number(g, "dtemperature_dcapacity_multiplier_k"), -power * 4.0 / 1000.0, 1e-8);
        let region = &g.get("regional_power").unwrap().as_array().unwrap()[0];
        near(number(region, "dtemperature_duniform_added_watts_k_per_w"), 4.0 / 1000.0, 1e-8);
        let (plain, plain_field) = run(cx, input, None);
        assert_eq!(field, plain_field);
        assert_eq!(plain, format!("{}}}", receipt.split_once(",\"adjoint\":").unwrap().0));
        let (repeated, _) = run(cx, input, Some("temperature-final-adjoint"));
        assert_eq!(repeated, receipt);
    });
}

#[test]
fn native_transient_adjoint_pulse_peak_precedes_final_and_future_controls_are_zero() {
    with_cx(|cx| {
        let input = Inputs { power: 20.0, capacity: 50.0, htc: 8.0, pulse: true,
            nonlinear: true, ..Inputs::default() };
        let (peak, _) = run(cx, input, Some("temperature-sampled-peak-adjoint"));
        let (final_time, _) = run(cx, input, Some("temperature-final-adjoint"));
        assert!(value(&peak) > value(&final_time) + 1e-3);
        let parsed = J::parse(&peak).unwrap();
        let g = parsed.get("adjoint").unwrap();
        assert!(number(g, "time_s") < 4.0);
        let state = number(g, "state_index") as usize;
        let intervals = g.get("intervals").unwrap().as_array().unwrap();
        for row in &intervals[state..] {
            near(number(row, "dtemperature_dsource_multiplier_k"), 0.0, 0.0);
            for region in row.get("regional_power").unwrap().as_array().unwrap() {
                near(number(region, "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
            }
        }
        let pulse_derivative: f64 = intervals.iter().filter(|row| number(row, "start_s") < 1.0)
            .map(|row| number(&row.get("regional_power").unwrap().as_array().unwrap()[0],
                "dtemperature_ddelivered_watts_k_per_w")).sum();
        let h = 1e-3;
        let (plus, _) = run(cx, Inputs { power: input.power + h, ..input }, Some("temperature-sampled-peak-adjoint"));
        let (minus, _) = run(cx, Inputs { power: input.power - h, ..input }, Some("temperature-sampled-peak-adjoint"));
        near(pulse_derivative, (value(&plus) - value(&minus)) / (2.0*h), 2e-4);
        let (plus, _) = run(cx, Inputs { capacity: input.capacity*(1.0+h), ..input }, Some("temperature-final-adjoint"));
        let (minus, _) = run(cx, Inputs { capacity: input.capacity*(1.0-h), ..input }, Some("temperature-final-adjoint"));
        let end = J::parse(&final_time).unwrap();
        near(number(end.get("adjoint").unwrap(), "dtemperature_dcapacity_multiplier_k"),
            (value(&plus) - value(&minus)) / (2.0*h), 2e-4);
        let (plus, _) = run(cx, Inputs { htc: input.htc * h.exp(), ..input }, Some("temperature-final-adjoint"));
        let (minus, _) = run(cx, Inputs { htc: input.htc * (-h).exp(), ..input }, Some("temperature-final-adjoint"));
        let robin = &end.get("adjoint").unwrap().get("robin").unwrap().as_array().unwrap()[0];
        near(number(robin, "dtemperature_dlog_htc_k"), (value(&plus)-value(&minus))/(2.0*h), 2e-4);
    });
}

#[test]
fn native_transient_adjoint_initial_peak_requires_no_reverse_steps() {
    with_cx(|cx| {
        let input = Inputs { initial: 320.0, power: 0.0, htc: 8.0, capacity: 50.0, ..Inputs::default() };
        let (receipt, _) = run(cx, input, Some("temperature-sampled-peak-adjoint"));
        let parsed = J::parse(&receipt).unwrap();
        let g = parsed.get("adjoint").unwrap();
        near(number(g, "time_s"), 0.0, 0.0);
        near(number(g, "value_k"), 320.0, 0.0);
        near(number(g, "dtemperature_duniform_initial_k_per_k"), 1.0, 0.0);
        near(number(g, "reconstructed_solid_endpoints"), 0.0, 0.0);
        near(number(g, "dtemperature_dcapacity_multiplier_k"), 0.0, 0.0);
        assert!(number(g, "active_ties") > 1.0);
    });
}

#[test]
fn native_transient_adjoint_rejects_ignored_requests_and_incomplete_tapes() {
    let mut spec = project(Inputs::default(), Some("temperature-final-adjoint"));
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient = None;
    assert!(requested(&spec).is_err());
    with_cx(|cx| {
        let mut tape = Tape::new(cx, &[300.0, 300.0], 1, None).unwrap();
        assert!(select(cx, &tape, Objective::Final, &[0], None).is_err());
        assert!(tape.record(cx, &[300.0], None).is_err());
        assert!(tape.record(cx, &[f64::NAN, 300.0], None).is_err());
        assert_eq!(tape.states.len(), 1);
        tape.record(cx, &[301.0, 301.0], None).unwrap();
        let selected = select(cx, &tape, Objective::SampledPeak, &[0, 1], None).unwrap();
        assert_eq!((selected.state, selected.vertex, selected.ties), (1, 0, 2));
        assert!(tape.record(cx, &[302.0, 302.0], None).is_err());
    });
    assert!(storage_bytes(usize::MAX, 1, 1).is_err());
}
