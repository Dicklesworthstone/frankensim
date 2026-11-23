//! The nine-state fixtures of `battery_motor`, whose switches store no
//! output charge, as declared graphs (`battery_motor/graph.rs`).

use fs_system::commutated_drive::NoDeviceStorage;

use super::graph::{
    Case, Row, compare, kernel_refusal, kernel_run, leak, owner_check, steady_check,
};
use super::*;

fn nine_state_case<'s>(
    f: &'static Fixture,
    steady: &'s SteadyDrive,
    dt_s: f64,
    substeps: usize,
    load_cell: Option<CellId>,
) -> Case<'s, NoDeviceStorage> {
    Case {
        fixture: f,
        context: leak(common::point(300.0)),
        storage_keys: &[],
        storage: NoDeviceStorage,
        steady,
        start: steady.state,
        dt_s,
        substeps,
        load_cell,
        schedule: Vec::new(),
        linked: false,
    }
}

/// The module's `steps` steps of the case from its bound drive, each
/// advancing the cells by `advance_s`; its first refusal, as text.
fn module_run(
    case: &Case<'_, NoDeviceStorage>,
    advance_s: f64,
    steps: usize,
) -> Result<Vec<Row<[f64; 9]>>, String> {
    let f = case.fixture;
    let cells = f.cells(case.context);
    let resistors = drive_resistors(f);
    let devices = f.devices();
    let network = CommutatedBatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
        bearing_cell: IDS[0],
    };
    let mut thermal = f.thermal();
    let error = |error| format!("{error:?}");
    let mut drive = network
        .bind(&thermal, &case.steady.native, case.start, || false)
        .map_err(error)?;
    let mut rows = Vec::new();
    for k in 0..steps {
        let intervals = [CommutationInterval {
            gates: case.gates(k),
            duration_s: case.dt_s,
        }];
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
fn the_nine_state_drive_as_a_graph_keeps_the_steady_motoring_and_regenerative_drive() {
    for regeneration in [false, true] {
        let f = leak(fixture(1.0, 0.0, 350.0, 100.0));
        let steady = steady_drive(regeneration);
        let case = nine_state_case(f, &steady, 1.0 / 64.0, 4, Some(IDS[3]));
        let module = module_run(&case, case.dt_s, 3).unwrap();
        // On 4d192a789c the two agree bit for bit. The bounds allow for
        // rounding on another platform.
        let graph = kernel_run(&case, 3);
        compare(&module, &graph, [1e-14, 1e-11]);
        // Without the module: within 1.1e-15 relative and 2.7e-14 K of the
        // closed form.
        steady_check(&case, &graph, [2e-14, 1e-12]);
    }
}

#[test]
fn the_nine_state_drive_as_a_graph_follows_the_hot_diode_threshold() {
    // The module's feedback fixture: cells of 0.001 J/K and a diode
    // threshold rising 20 mV/K.
    let f = leak(fixture(0.001, 0.02, 350.0, 100.0));
    let steady = steady_drive(false);
    let case = nine_state_case(f, &steady, 1.0 / 1024.0, 1, None);
    let module = module_run(&case, case.dt_s, 3).unwrap();
    // On 4d192a789c the two agree within 4.4e-14 relative and 1.2e-9 K
    // (1.2e-12 J in these cells) over three steps; the bounds sit 11 and 8
    // times above.
    let graph = kernel_run(&case, 3);
    compare(&module, &graph, [5e-13, 1e-8]);
    // Without the module: the owner by hand bit for bit, each cell's
    // balance within 5.3e-17 J.
    owner_check(&case, &graph, 0.02, [1e-14, 1e-14]);
    // The output high package warms past 300.1 K in the first step, and
    // the converter's flux moves.
    assert!(graph[0].kelvin[3] > 300.1, "{:?}", graph[0].kelvin);
    assert!((graph[0].state[3] - steady.state[3]).abs() > 1e-8);
}

#[test]
fn the_nine_state_drive_as_a_graph_refuses_what_the_module_refuses() {
    let steady = steady_drive(false);
    let dt = 1.0 / 1024.0;
    let case = |upper, limit| {
        nine_state_case(
            leak(fixture(0.001, 0.02, upper, limit)),
            &steady,
            dt,
            4,
            None,
        )
    };
    // The card, cells and devices alike, qualified only up to 300.0001 K,
    // which the step passes. The module refuses the cells' state. In the
    // graph the cells' solve rejects the step, the timeline asks for half of
    // it, and the drive refuses a step that is not one period.
    let hot = case(300.0001, 100.0);
    let module = module_run(&hot, dt, 1).unwrap_err();
    assert!(module.contains("MaterialDomain"), "{module}");
    let graph = kernel_refusal(&hot, dt);
    assert!(graph.contains("one period of its program"), "{graph}");
    // The owner's channel limit: the same refusal on both sides.
    let limited = case(350.0, 0.01);
    let module = module_run(&limited, dt, 1).unwrap_err();
    let graph = kernel_refusal(&limited, dt);
    for text in [&module, &graph] {
        assert!(text.contains("leaves authored support"), "{text}");
    }
    assert!(graph.contains("Step("), "{graph}");
    // Qualified to 350 K, the graph takes the step; a step that is not one
    // gate program is refused.
    let case = case(350.0, 100.0);
    assert_eq!(kernel_run(&case, 1).len(), 1);
    let module = module_run(&case, 2.0 * dt, 1).unwrap_err();
    assert!(
        module.contains("same represented thermal endpoint"),
        "{module}"
    );
    let graph = kernel_refusal(&case, 2.0 * dt);
    assert!(graph.contains("one period of its program"), "{graph}");
}
