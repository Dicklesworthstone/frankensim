//! The constant-capacitance fixtures as declared graphs (`battery_motor/graph.rs`).

use fs_system::commutated_drive::DeviceStorage;

use super::super::graph::{
    Case, Row, cells_start, closed_form_step, compare, heat_assembly, kernel_refusal, kernel_run,
    leak, owner_check, steady_check,
};
use super::*;

fn capacitive_case<'s>(
    f: &'static CapacitiveFixture,
    steady: &'s SteadyDrive,
    start: [f64; 11],
    dt_s: f64,
    substeps: usize,
    load_cell: Option<CellId>,
) -> Case<'s> {
    Case {
        fixture: &f.base,
        context: leak(common::point(300.0)),
        storage_keys: &f.keys,
        storage: DeviceStorage::Capacitance(COSS),
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
    f: &CapacitiveFixture,
    case: &Case<'_>,
    advance_s: f64,
    steps: usize,
) -> Result<Vec<Row>, String> {
    let context = common::point(300.0);
    let cells = f.base.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let caps = f.bindings();
    let network = coupled(&cells, &resistors, &devices, &caps);
    let mut thermal = f.base.thermal();
    let native = &case.steady.native;
    let error = |error| format!("{error:?}");
    let mut drive = network
        .bind(&thermal, native, COSS, case.start, || false)
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
fn module_linked_run(f: &CapacitiveFixture, case: &Case<'_>, steps: usize) -> Vec<Row> {
    let start = &cells_start(case);
    heat_assembly(case, start, |cx, transport| {
        let (resistors, devices) = (drive_resistors(&f.base), f.base.devices());
        let drive = CommutatedBatteryMotorHeatNetwork {
            assembly: transport,
            resistors: &resistors,
            devices: &devices,
            bearing_cell: IDS[0],
        };
        let caps = f.bindings();
        let network = CapacitiveBatteryMotorHeatNetwork {
            drive,
            capacitances: &caps,
        };
        let native = &case.steady.native;
        let mut drive = network.bind(cx, start, native, COSS, case.start).unwrap();
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
fn constant_capacitors_as_a_graph_keep_the_steady_motoring_and_regenerative_drive() {
    for regeneration in [false, true] {
        let f = leak(capacitive_fixture(1.0, 0.0, 350.0, false, 1.0));
        let steady = steady_drive(regeneration);
        let start = initial(&steady, [0.0; 2]);
        let case = capacitive_case(f, &steady, start, 1.0 / 64.0, 4, Some(IDS[3]));
        let module = module_run(f, &case, case.dt_s, 3).unwrap();
        // On 38d60e3722 the two agree bit for bit. The bounds
        // allow for rounding on another platform.
        let graph = kernel_run(&case, 3);
        compare(&module, &graph, [1e-14, 1e-11]);
        // Without the module: within 1.1e-15 relative and 2.7e-14 K of the
        // closed form.
        steady_check(&case, &graph, [2e-14, 1e-12]);
    }
}

#[test]
fn constant_capacitors_as_a_graph_follow_the_hot_diode_threshold() {
    // The module's dynamic fixture: cells of 0.01 J/K, a diode threshold
    // rising 20 mV/K, the nodes started off their steady values.
    let f = leak(capacitive_fixture(0.01, 0.02, 350.0, false, 1.0));
    let steady = steady_drive(false);
    let start = initial(&steady, [-0.3, 0.4]);
    let case = capacitive_case(f, &steady, start, 1.0 / 1024.0, 1, None);
    let module = module_run(f, &case, case.dt_s, 3).unwrap();
    // On 38d60e3722 the two agree within 4.2e-15 relative and
    // 2.2e-11 K over three steps; the bounds sit 24 and 45 times above.
    let graph = kernel_run(&case, 3);
    compare(&module, &graph, [1e-13, 1e-9]);
    // Without the module: the owner by hand bit for bit, each cell's
    // balance within 5.5e-16 J.
    owner_check(&case, &graph, 0.02, [1e-14, 1e-14]);
    // The feedback is real: the output high package, whose diode conducts,
    // warms past 300.8 K, so its threshold rises by more than 16 mV.
    assert!(module[2].kelvin[3] > 300.8, "{:?}", module[2].kelvin);
}

#[test]
fn constant_capacitors_as_a_graph_follow_a_gate_program_that_changes() {
    // The module's dynamic fixture continues with only the output high
    // switch on; a third step returns to the first gates.
    let f = leak(capacitive_fixture(0.01, 0.02, 350.0, false, 1.0));
    let steady = steady_drive(false);
    let start = initial(&steady, [-0.3, 0.4]);
    let other = DeviceGateState {
        on: [false, false, true, false],
    };
    let dt = 1.0 / 1024.0;
    let case = Case {
        schedule: vec![steady.gates, other, steady.gates],
        ..capacitive_case(f, &steady, start, dt, 1, None)
    };
    let module = module_run(f, &case, dt, 3).unwrap();
    // On a8d5a7f02c the two agree within 2.8e-15 relative and 1.4e-11 K
    // over three steps; the bounds sit 36 times above.
    let graph = kernel_run(&case, 3);
    compare(&module, &graph, [1e-13, 5e-10]);
    // Without the module: the owner by hand bit for bit, each cell's
    // balance within 5.5e-16 J.
    owner_check(&case, &graph, 0.02, [1e-14, 1e-14]);
    // The schedule matters: with the first gates kept, the second step ends
    // elsewhere.
    let kept = module_run(f, &capacitive_case(f, &steady, start, dt, 1, None), dt, 3).unwrap();
    assert_ne!(
        kept[1].state.map(f64::to_bits),
        module[1].state.map(f64::to_bits)
    );
}

#[test]
fn constant_capacitors_as_a_graph_carry_the_link_between_cells() {
    // The module's heat-assembly fixture: steady motoring with a link of
    // 0.25 W/K between the inductor's cell and the motor's. Then the hot
    // dynamic fixture with the same link.
    let steady = steady_drive(false);
    let fixtures = [
        ((1.0, 0.0), [0.0; 2], 1.0 / 64.0, 4, [1e-14, 1e-11]),
        ((0.01, 0.02), [-0.3, 0.4], 1.0 / 1024.0, 1, [1e-13, 1e-9]),
    ];
    // On 72250ec7fc the steady fixture agrees bit for bit in the drive and
    // within 1.7e-13 K in the cells, the hot one bit for bit. The bounds are
    // those of the same fixtures without the link.
    for ((cp, slope), offset, dt, substeps, bounds) in fixtures {
        let f = leak(capacitive_fixture(cp, slope, 350.0, false, 1.0));
        let start = initial(&steady, offset);
        let case = Case {
            linked: true,
            ..capacitive_case(f, &steady, start, dt, substeps, None)
        };
        let module = module_linked_run(f, &case, 3);
        let graph = kernel_run(&case, 3);
        compare(&module, &graph, bounds);
        // Without the module: the owner by hand bit for bit, each cell's
        // balance within 4.4e-14 J.
        owner_check(&case, &graph, slope, [1e-14, 1e-12]);
        // The steady fixture's first step also meets the module's own
        // closed form on the graph: within 5.7e-14 K.
        assert!(substeps != 4 || closed_form_step(&case) < 1e-12);
    }
}

#[test]
fn constant_capacitors_as_a_graph_refuse_what_the_module_refuses() {
    let steady = steady_drive(false);
    let start = initial(&steady, [0.0; 2]);
    let dt = 1.0 / 64.0;
    // Capacitors qualified only up to 300.000001 K, which the step passes.
    let f = leak(capacitive_fixture(1.0, 0.0, 300.000001, false, 1.0));
    let case = capacitive_case(f, &steady, start, dt, 4, None);
    let module = module_run(f, &case, dt, 1).unwrap_err();
    let graph = kernel_refusal(&case, dt);
    for text in [&module, &graph] {
        assert!(text.contains("NoClaimInDomain"), "{text}");
        assert!(text.contains("device_output_capacitance_0"), "{text}");
    }
    assert!(
        graph.contains(r#"Law { port: "heat_input_high""#),
        "{graph}"
    );
    // A step that is not one gate program.
    let f = leak(capacitive_fixture(1.0, 0.0, 350.0, false, 1.0));
    let case = capacitive_case(f, &steady, start, dt, 4, None);
    let module = module_run(f, &case, 2.0 * dt, 1).unwrap_err();
    assert!(
        module.contains("same represented thermal endpoint"),
        "{module}"
    );
    let graph = kernel_refusal(&case, 2.0 * dt);
    assert!(graph.contains("one period of its program"), "{graph}");
}
