//! The nine-state fixtures of `battery_motor`, whose switches store no
//! output charge, as declared graphs (`battery_motor/graph.rs`).

use fs_system::commutated_drive::NoDeviceStorage;

use super::graph::{
    Case, Coolant, Row, cells_start, compare, heat_assembly, kernel_refusal, kernel_run, leak,
    owner_check, steady_check,
};
use super::*;

/// The coolant's cell in the module's heat-assembly fixture.
const COOLANT: CellId = CellId(80);

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
        coolant: None,
        drivers: None,
    }
}

/// The module's coolant: the fixture's claims as a liquid, 2 kg at 299.5 K
/// (`commutated_battery_drive_and_finite_liquid_coolant_share_the_original_heat_assembly`).
pub(super) fn coolant(f: &Fixture) -> Coolant {
    let mut claims = ClaimSet::new();
    for (_, claim) in f.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let id = MaterialStateId {
        chemistry: "synthetic-coolant".into(),
        phase: "liquid".into(),
        process: "authored finite heat-capacity fixture".into(),
        revision: 0,
    };
    let card = MaterialCard::assemble(id, claims, vec![]).unwrap();
    let query = HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS),
        temperature_axis: "T",
        lower: &common::point(299.0),
        upper: &common::point(350.0),
        selection: SelectionPolicy::SingleClaimOnly,
    };
    let caloric =
        ResolvedHeatCapacityCurve::try_new(&card, query, SolidLiquidPhase::Liquid, 0.0, 1e-6)
            .unwrap();
    let cell = DepositedCell {
        id: COOLANT,
        segment_id: 7,
        body_id: 7,
        material_id: MaterialId(1),
        mass_kg: 2.0,
        specific_h_j_kg: 0.5,
        geometry: CellGeometry {
            start_m: [20.0, 0.0, 0.0],
            end_m: [22.0, 0.0, 0.0],
            volume_m3: 2.0,
            surface_area_m2: 10.0,
        },
        boundary: common::boundary(),
    };
    Coolant {
        card,
        caloric,
        context: common::point(299.5),
        cell,
        partner: (IDS[0], super::graph::LINK.2),
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

/// The module's `steps` steps of the case through its heat assembly, which
/// carries the coolant and its link (`FixedCooling`).
fn module_linked_run(
    case: &Case<'_, NoDeviceStorage>,
    steps: usize,
) -> Vec<Row<[f64; 9], [f64; 8]>> {
    let start = &cells_start(case);
    let ids: [CellId; 8] = core::array::from_fn(|i| IDS.get(i).copied().unwrap_or(COOLANT));
    heat_assembly(case, start, |cx, transport| {
        let (resistors, devices) = (drive_resistors(case.fixture), case.fixture.devices());
        let network = CommutatedBatteryMotorHeatNetwork {
            assembly: transport,
            resistors: &resistors,
            devices: &devices,
            bearing_cell: IDS[0],
        };
        let native = &case.steady.native;
        let mut drive = network.bind(cx, start, native, case.start).unwrap();
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
                    drive_controls(8),
                    &[0.0; 8],
                )
                .unwrap();
            rows.push(Row {
                state: *next.drive().state(),
                kelvin: ids.map(|id| next.thermal().temperature_k(id).unwrap()),
            });
            (thermal, drive, _) = next.into_parts();
        }
        rows
    })
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
fn the_nine_state_drive_as_a_graph_shares_a_finite_liquid_coolant() {
    let f = leak(fixture(1.0, 0.0, 350.0, 100.0));
    let steady = steady_drive(false);
    let case = Case {
        coolant: Some(leak(coolant(f))),
        ..nine_state_case(f, &steady, 1.0 / 64.0, 4, None)
    };
    let module = module_linked_run(&case, 3);
    // On 4d192a789c the drive agrees bit for bit and the cells within
    // 1.1e-13 K over three steps.
    let graph = kernel_run(&case, 3);
    compare(&module, &graph, [1e-14, 1e-11]);
    // Without the module: the owner by hand bit for bit, each cell's
    // balance, the coolant's by its own caloric curve, within 8.1e-14 J.
    owner_check(&case, &graph, 0.0, [1e-14, 1e-12]);
    // The module's closed form for the first step, backward Euler between
    // the inductor's cell of 1 J/K and the coolant's 2 J/K from 299.5 K: met
    // exactly on 4d192a789c.
    let (h, g) = (case.dt_s, 0.25);
    let q = per_body(steady.heat_w.map(|p| h * p), 0.0);
    let x = h * g * (0.5 + q[0]) / (1.0 + h * g * 1.5);
    let mut expected = q.map(|q| 300.0 + q).to_vec();
    expected[0] -= x;
    expected.push(299.5 + x / 2.0);
    let worst = graph[0]
        .kelvin
        .iter()
        .zip(&expected)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(worst <= 1e-12, "{worst:e}");
    // The coolant warms through the link.
    assert!(graph[2].kelvin[7] > graph[0].kelvin[7] && graph[0].kelvin[7] > 299.5);
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
    assert_eq!(kernel_run::<_, [f64; 7]>(&case, 1).len(), 1);
    let module = module_run(&case, 2.0 * dt, 1).unwrap_err();
    assert!(
        module.contains("same represented thermal endpoint"),
        "{module}"
    );
    let graph = kernel_refusal(&case, 2.0 * dt);
    assert!(graph.contains("one period of its program"), "{graph}");
}
