//! G1/G3 scheduled regional power through the native imported-volume workflow.
//! The energy oracles use the declared watts and physical window, independently
//! of the numerical grid and the producer's workload-summary calculation.

use super::*;
use fs_project::{ConductionTransient, TransientPowerStep, TransientRegionPower, spec::dims};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const INITIAL_K: f64 = 293.15;
const TETRA_VOLUME_M3: f64 = 1.0 / 6.0;

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
}

fn project(capacity: f64, horizon: f64, max_step: f64, insulated: bool) -> ProjectSpec {
    let mut spec = conduction_fixture_project(7, &tetra_stl());
    transient_product::declaration(&mut spec, capacity, horizon);
    time(&mut spec).max_step = QtyAny::new(max_step, dims::TIME);
    if insulated {
        spec.cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .boundaries[0]
            .condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
    }
    spec
}

fn schedule(spec: &mut ProjectSpec, steps: &[(f64, f64)]) {
    time(spec).power_schedules = vec![TransientRegionPower {
        region: "air".into(),
        source: "synthetic delivered regional watts for the thermal history oracle".into(),
        steps: steps
            .iter()
            .map(|&(until, watts)| TransientPowerStep {
                until: QtyAny::new(until, dims::TIME),
                watts: QtyAny::new(watts, dims::POWER),
            })
            .collect(),
    }];
}

struct History {
    run: String,
    qoi: J,
    conduction: J,
    temperatures: Vec<f64>,
}

impl History {
    fn transient(&self) -> &J {
        self.conduction.get("transient").unwrap()
    }

    fn grid(&self, name: &str) -> &[J] {
        self.transient().get(name).unwrap().as_array().unwrap()
    }

    fn final_k(&self) -> f64 {
        n(
            &self.qoi.get("qoi").unwrap().as_array().unwrap()[0],
            "value",
        )
    }
}

fn run(spec: &ProjectSpec) -> History {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, tetra_stl());
    let (run, qoi, conduction, field) =
        run_conjugate_to_completion(&ledger, &decode(spec), &fixture_cards());
    let field = J::parse(&field).unwrap();
    History {
        run,
        qoi: J::parse(&qoi).unwrap(),
        conduction: J::parse(&conduction).unwrap(),
        temperatures: field
            .get("temperature")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect(),
    }
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

fn delivered_source_j(rows: &[J]) -> f64 {
    let mut previous_time = 0.0;
    rows.iter()
        .map(|row| {
            let end = n(row, "time_s");
            assert!(end > previous_time);
            let energy = (end - previous_time) * n(row, "source_w");
            previous_time = end;
            energy
        })
        .sum()
}

#[test]
fn g1_native_workload_delivers_absolute_pulse_energy_without_static_power_or_duty_scaling() {
    let capacity = 1000.0;
    let mut spec = project(capacity, 2.0, 0.5, true);
    // Scheduled watts replace this 1.25 W baseline; they are not multiplied by
    // its duty factor, nor added to its static generation.
    spec.power.as_mut().unwrap()[0].duty = 0.25;
    schedule(&mut spec, &[(0.5, 10.0), (1.5, 0.0), (2.0, 20.0)]);
    let history = run(&spec);
    let expected_j = 10.0 * 0.5 + 20.0 * 0.5;
    let expected_k = INITIAL_K + expected_j / (capacity * TETRA_VOLUME_M3);
    close(history.final_k(), expected_k, 2e-10);
    for &temperature in &history.temperatures {
        close(temperature, expected_k, 2e-10);
    }
    for grid in ["coarse", "fine"] {
        close(delivered_source_j(history.grid(grid)), expected_j, 2e-10);
        for row in history.grid(grid) {
            let end = n(row, "time_s");
            let watts = if end <= 0.5 {
                10.0
            } else if end <= 1.5 {
                0.0
            } else {
                20.0
            };
            close(n(row, "source_w"), watts, 2e-10);
        }
    }
    let workload = history.transient().get("workload").unwrap();
    assert_eq!(
        workload.str_field("mode"),
        Some("piecewise-constant-regional-power")
    );
    close(
        n(workload, "integrated_scheduled_input_j"),
        expected_j,
        2e-12,
    );
    close(
        n(
            history.transient().get("energy").unwrap(),
            "stored_change_j",
        ),
        expected_j,
        2e-7,
    );
    let scope = history.qoi.get("thermal_time_scope").unwrap();
    assert_eq!(scope.str_field("kind"), Some("final"));
    assert_eq!(
        scope.get("continuous_time_peak_measured"),
        Some(&J::Bool(false))
    );
}

#[test]
fn g1_native_workload_resolves_a_pulse_missed_by_both_original_uniform_grids() {
    let capacity = 1000.0;
    let mut spec = project(capacity, 1.0, 0.5, true);
    // Neither the original 0.5 s grid nor its 0.25 s half grid samples this
    // pulse. Its two switch times must become endpoints on BOTH real grids.
    schedule(&mut spec, &[(0.125, 0.0), (0.1875, 16.0), (1.0, 0.0)]);
    let history = run(&spec);
    close(n(history.transient(), "coarse_steps"), 4.0, 0.0);
    close(n(history.transient(), "fine_steps"), 8.0, 0.0);
    close(n(history.transient(), "total_steps"), 12.0, 0.0);
    for grid in ["coarse", "fine"] {
        let rows = history.grid(grid);
        for switch in [0.125, 0.1875, 1.0] {
            assert!(rows.iter().any(|row| n(row, "time_s") == switch));
        }
        close(delivered_source_j(rows), 1.0, 2e-10);
    }
    close(
        history.final_k(),
        INITIAL_K + 1.0 / (capacity * TETRA_VOLUME_M3),
        2e-10,
    );
    close(
        n(history.conduction.get("energy").unwrap(), "source_w"),
        0.0,
        2e-10,
    );
}

#[test]
fn g3_native_workload_equal_energy_early_and_late_heating_have_different_final_temperatures() {
    let mut early = project(100.0, 2.0, 0.25, false);
    schedule(&mut early, &[(1.0, 10.0), (2.0, 0.0)]);
    let mut late = early.clone();
    schedule(&mut late, &[(1.0, 0.0), (2.0, 10.0)]);
    let early = run(&early);
    let late = run(&late);
    for history in [&early, &late] {
        close(delivered_source_j(history.grid("fine")), 10.0, 2e-10);
        assert!(history.final_k() > INITIAL_K);
        assert!(
            n(
                history.transient().get("energy").unwrap(),
                "maximum_step_residual_j"
            ) <= 1e-6
        );
    }
    assert!(
        late.final_k() > early.final_k() + 1e-3,
        "a Robin-cooled body must retain more of the later equal-energy input"
    );
    assert_ne!(
        early.run, late.run,
        "the full physical schedule belongs to run identity"
    );
}

#[test]
fn g3_native_constant_workload_reproduces_static_delivered_power() {
    let mut baseline = project(1000.0, 2.0, 0.5, true);
    baseline.power.as_mut().unwrap()[0].watts.value = 20.0;
    baseline.power.as_mut().unwrap()[0].duty = 0.25;
    let mut explicit = baseline.clone();
    schedule(&mut explicit, &[(2.0, 5.0)]);
    let baseline = run(&baseline);
    let explicit = run(&explicit);
    assert_eq!(baseline.temperatures.len(), explicit.temperatures.len());
    for (&a, &b) in baseline.temperatures.iter().zip(&explicit.temperatures) {
        close(a, b, 2e-10);
    }
    close(baseline.final_k(), explicit.final_k(), 2e-10);
    for grid in ["coarse", "fine"] {
        assert_eq!(baseline.grid(grid).len(), explicit.grid(grid).len());
        close(delivered_source_j(explicit.grid(grid)), 10.0, 2e-10);
    }
}

#[test]
fn g0_native_workload_step_cap_accounts_for_switches_before_numerical_work() {
    let mut spec = project(1000.0, 1.0, 0.5, true);
    schedule(&mut spec, &[(0.125, 0.0), (0.1875, 16.0), (1.0, 0.0)]);
    // Nine steps admit the old uniform 2+4 path, but not the required 4+8 path.
    time(&mut spec).max_steps = 9;
    assert!(
        spec.validate()
            .iter()
            .any(|v| v.code == "project-conduction-transient-steps")
    );
    time(&mut spec).max_steps = 12;
    assert!(spec.validate().is_empty(), "{:?}", spec.validate());
}
