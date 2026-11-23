//! The Q(V) charge fixtures as declared graphs (`battery_motor/graph.rs`).

use fs_system::commutated_drive::DeviceStorage;

use super::super::graph::{Case, Row, compare, kernel_refusal, kernel_run, leak};
use super::*;

fn charge_case<'s>(
    f: &'static ChargeFixture,
    steady: &'s SteadyDrive,
    start: [f64; 11],
    dt_s: f64,
    substeps: usize,
    load_cell: Option<CellId>,
) -> Case<'s> {
    Case {
        fixture: &f.base,
        context: leak(context(300.0)),
        storage_keys: &f.keys,
        storage: DeviceStorage::Charge {
            curves: curves(1.0),
            maximum_cells_per_chord: CHORD_BUDGET,
        },
        steady,
        start,
        dt_s,
        substeps,
        load_cell,
    }
}

/// The module's `steps` steps of the case from its bound drive, each
/// advancing the cells by `advance_s`; its first refusal, as text.
fn module_run(
    f: &ChargeFixture,
    case: &Case<'_>,
    advance_s: f64,
    steps: usize,
) -> Result<Vec<Row>, String> {
    let context = context(300.0);
    let cells = f.base.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let charges = f.bindings();
    let network = coupled(&cells, &resistors, &devices, &charges);
    let mut thermal = f.base.thermal();
    let authored = curves(1.0);
    let native = &case.steady.native;
    let error = |error| format!("{error:?}");
    let mut drive = network
        .bind(&thermal, native, &authored, case.start, || false)
        .map_err(error)?;
    let intervals = [CommutationInterval {
        gates: case.steady.gates,
        duration_s: case.dt_s,
    }];
    let mut rows = Vec::new();
    for _ in 0..steps {
        let next = network
            .step_candidate(
                &thermal,
                thermal.time_s() + advance_s,
                &drive,
                &intervals,
                &[case.steady.command],
                case.substeps,
                case.load_cell,
                power_controls(),
                cycle_controls(),
                drive_controls(7),
                &[0.0; 7],
                || false,
            )
            .map_err(error)?;
        rows.push(Row {
            state: *next.drive().state(),
            kelvin: IDS.map(|id| next.thermal().temperature_k(id).unwrap()),
        });
        (thermal, drive, _) = next.into_parts();
    }
    Ok(rows)
}

#[test]
fn charge_curves_as_a_graph_keep_the_steady_motoring_and_regenerative_drive() {
    for regeneration in [false, true] {
        let f = leak(charge_fixture(1.0, 0.0, 350.0, false, 1.0));
        let steady = steady_drive(regeneration);
        let start = initial(&steady, [0.0; 2]);
        let case = charge_case(f, &steady, start, 1.0 / 64.0, 4, Some(IDS[3]));
        let module = module_run(f, &case, case.dt_s, 3).unwrap();
        // On the export of 1e611fdf33 the two agree bit for bit. The bounds
        // allow for rounding on another platform.
        compare(&module, &kernel_run(&case, 3), [1e-14, 1e-11]);
    }
}

#[test]
fn charge_curves_as_a_graph_follow_the_hot_diode_threshold() {
    // The module's dynamic fixture: cells of 0.01 J/K, a diode threshold
    // rising 20 mV/K, the nodes started off their steady values.
    let f = leak(charge_fixture(0.01, 0.02, 350.0, false, 1.0));
    let steady = steady_drive(false);
    let start = initial(&steady, [-0.45, 0.4]);
    let case = charge_case(f, &steady, start, 1.0 / 1024.0, 1, None);
    let module = module_run(f, &case, case.dt_s, 3).unwrap();
    // On the export of 1e611fdf33 the two agree within 9.3e-15 relative and
    // 4.9e-11 K over three steps; the bounds sit 21 and 20 times above.
    compare(&module, &kernel_run(&case, 3), [2e-13, 1e-9]);
    // The feedback is real: the output high package, whose diode conducts,
    // warms past 300.75 K, so its threshold rises by more than 15 mV.
    assert!(module[2].kelvin[3] > 300.75, "{:?}", module[2].kelvin);
}

#[test]
fn charge_curves_as_a_graph_refuse_what_the_module_refuses() {
    let steady = steady_drive(false);
    let start = initial(&steady, [0.0; 2]);
    let dt = 1.0 / 64.0;
    // Charge curves qualified only up to 300.000001 K, which the step passes.
    let f = leak(charge_fixture(1.0, 0.0, 300.000001, false, 1.0));
    let case = charge_case(f, &steady, start, dt, 4, None);
    let module = module_run(f, &case, dt, 1).unwrap_err();
    let graph = kernel_refusal(&case, dt);
    for text in [&module, &graph] {
        assert!(text.contains("NoClaimInDomain"), "{text}");
        assert!(text.contains("device_output_charge_0"), "{text}");
    }
    assert!(
        graph.contains(r#"Law { port: "heat_input_high""#),
        "{graph}"
    );
    // A step that is not one gate program.
    let f = leak(charge_fixture(1.0, 0.0, 350.0, false, 1.0));
    let case = charge_case(f, &steady, start, dt, 4, None);
    let module = module_run(f, &case, 2.0 * dt, 1).unwrap_err();
    assert!(
        module.contains("same represented thermal endpoint"),
        "{module}"
    );
    let graph = kernel_refusal(&case, 2.0 * dt);
    assert!(graph.contains("one period of its program"), "{graph}");
}
