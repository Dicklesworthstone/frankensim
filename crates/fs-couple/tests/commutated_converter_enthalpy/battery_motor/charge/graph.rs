//! The Q(V) charge fixtures as declared graphs (`battery_motor/graph.rs`).

use fs_system::commutated_drive::DeviceStorage;

use super::super::graph::{
    Case, Row, cells_start, closed_form_step, compare, heat_assembly, kernel_refusal, kernel_run,
    leak,
};
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
        schedule: Vec::new(),
        linked: false,
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

/// The module's `steps` steps of the case through its heat assembly, which
/// carries the case's link (`FixedCooling`).
fn module_linked_run(f: &ChargeFixture, case: &Case<'_>, steps: usize) -> Vec<Row> {
    let start = &cells_start(case);
    heat_assembly(case, start, |cx, transport| {
        let (resistors, devices) = (drive_resistors(&f.base), f.base.devices());
        let drive = CommutatedBatteryMotorHeatNetwork {
            assembly: transport,
            resistors: &resistors,
            devices: &devices,
            bearing_cell: IDS[0],
        };
        let charges = f.bindings();
        let network = ChargeBatteryMotorHeatNetwork {
            drive,
            charges: &charges,
            maximum_cells_per_chord: CHORD_BUDGET,
        };
        let (native, authored) = (&case.steady.native, curves(1.0));
        let mut drive = network
            .bind(cx, start, native, &authored, case.start)
            .unwrap();
        let mut thermal = start.clone();
        let mut rows = Vec::new();
        for k in 0..steps {
            let intervals = [CommutationInterval {
                gates: case.gates(k),
                duration_s: case.dt_s,
            }];
            let next = network
                .step_candidate(
                    cx,
                    &thermal,
                    thermal.time_s() + case.dt_s,
                    &drive,
                    &intervals,
                    &[case.steady.command],
                    case.substeps,
                    case.load_cell,
                    power_controls(),
                    cycle_controls(),
                    drive_controls(7),
                    &[0.0; 7],
                )
                .unwrap();
            rows.push(Row {
                state: *next.drive().state(),
                kelvin: IDS.map(|id| next.thermal().temperature_k(id).unwrap()),
            });
            (thermal, drive, _) = next.into_parts();
        }
        rows
    })
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
fn charge_curves_as_a_graph_follow_a_gate_program_that_changes() {
    // The module's dynamic fixture continues with only the output high
    // switch on; a third step returns to the first gates.
    let f = leak(charge_fixture(0.01, 0.02, 350.0, false, 1.0));
    let steady = steady_drive(false);
    let start = initial(&steady, [-0.45, 0.4]);
    let other = DeviceGateState {
        on: [false, false, true, false],
    };
    let dt = 1.0 / 1024.0;
    let case = Case {
        schedule: vec![steady.gates, other, steady.gates],
        ..charge_case(f, &steady, start, dt, 1, None)
    };
    let module = module_run(f, &case, dt, 3).unwrap();
    // On a8d5a7f02c the two agree within 5.5e-15 relative and 2.7e-11 K
    // over three steps; the bounds sit 18 and 37 times above.
    compare(&module, &kernel_run(&case, 3), [1e-13, 1e-9]);
    // The schedule matters: with the first gates kept, the second step ends
    // elsewhere.
    let kept = module_run(f, &charge_case(f, &steady, start, dt, 1, None), dt, 3).unwrap();
    assert_ne!(
        kept[1].state.map(f64::to_bits),
        module[1].state.map(f64::to_bits)
    );
}

#[test]
fn charge_curves_as_a_graph_carry_the_link_between_cells() {
    // The module's heat-assembly fixture: steady motoring with a link of
    // 0.25 W/K between the inductor's cell and the motor's. Then the hot
    // dynamic fixture with the same link.
    let steady = steady_drive(false);
    let fixtures = [
        (
            charge_fixture(1.0, 0.0, 350.0, false, 1.0),
            [0.0; 2],
            1.0 / 64.0,
            4,
            [1e-14, 1e-11],
        ),
        (
            charge_fixture(0.01, 0.02, 350.0, false, 1.0),
            [-0.45, 0.4],
            1.0 / 1024.0,
            1,
            [1e-13, 1e-9],
        ),
    ];
    // On 72250ec7fc the steady fixture agrees bit for bit in the drive and
    // within 1.7e-13 K in the cells, the hot one within 7.8e-16 relative and
    // 1.6e-11 K. The bounds are those of the same fixtures without the link,
    // 1e-13 in place of 2e-13 for the hot drive.
    for (f, offset, dt, substeps, bounds) in fixtures {
        let f = leak(f);
        let start = initial(&steady, offset);
        let case = Case {
            linked: true,
            ..charge_case(f, &steady, start, dt, substeps, None)
        };
        let module = module_linked_run(f, &case, 3);
        compare(&module, &kernel_run(&case, 3), bounds);
        // The steady fixture's first step also meets the module's own
        // closed form on the graph: within 5.7e-14 K.
        assert!(substeps != 4 || closed_form_step(&case) < 1e-12);
    }
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
