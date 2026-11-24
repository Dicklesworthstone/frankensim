//! The module's voltage-reversal trial as a declared graph on the
//! composition kernel.
//!
//! `WindingThermalNetwork::step_hysteretic_candidate` couples fs-em's
//! play-memory winding circuit to the enthalpy cell by hand: the circuit's
//! step with every conductor at the cell's END temperature, Newton on the
//! cell energy. The graph says the same with no code for the pair:
//!
//! - `fs_em::winding::component::PlanarWinding::at_port_temperatures` on
//!   the play-memory owner: one heat port, every triangle's conductor, copper
//!   and core alike ([`CellConductors`]), following the cell's card at the
//!   temperature the port receives; the authored 0.3-ohm graph resistor keeps
//!   its ohms;
//! - the shared kernel's backward-Euler cell (`winding_graph/kernel.rs`);
//! - the source port following the step's start time: 20, 0, -30 and 0 V.
//!
//! Held, in this order, to the ledger, to `reference()` (the scalar P1/play
//! balance) from the kernel's own start of interval, and to the module, step
//! for step.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_em::winding::component::{ConductorLaw, HeatReceiver, PlanarWinding, PlanarWindingState};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_system::{NodeId, System};
use fs_thermal::component::TimelineState;

use super::common::MATERIAL;
use super::kernel::{CellConductors, conduction, graph, ledger_closes, timeline, with_cx};
use super::*;

const STEPS: usize = 4;
const DRIVES: [f64; STEPS] = [20.0, 0.0, -30.0, 0.0];
/// The cell at 300 K with no current, field or retained polarization.
const START: History = History {
    t: 300.0,
    i: 0.0,
    h: 0.0,
    z: 0.0,
    b: 0.0,
};
const TOLERANCE_J: f64 = 1e-9;

type PlayState = PlanarWindingState<'static, HystereticWindingCircuit<'static>>;

/// The source's volts over the step that starts at `t0`.
fn drive(t0: f64) -> f64 {
    DRIVES[(t0 / DT).round() as usize]
}

/// One accepted interval: the cell's temperature \[K\], the winding current
/// \[A\], the free vertices' potentials \[Wb/m\], each triangle's field
/// \[A/m\] and retained polarization \[T\], and the interval's heat, stored
/// change, numerical dissipation and supplied work \[J\].
#[derive(Debug, Clone)]
struct Row {
    kelvin: f64,
    current_a: f64,
    potential: [f64; 2],
    memory: Vec<[f64; 4]>,
    energies: [f64; 4],
}

fn row(kelvin: f64, circuit: &HystereticWindingCircuit<'_>, energies: [f64; 4]) -> Row {
    let potential = circuit.field().potential_wb_m();
    Row {
        kelvin,
        current_a: circuit.currents_a()[0],
        potential: [potential[6], potential[7]],
        memory: circuit
            .field()
            .magnetic_state()
            .cell_states()
            .iter()
            .map(|cell| {
                let (h, z) = (cell.field_a_m(), cell.polarizations_t()[0]);
                [h[0], h[1], z[0], z[1]]
            })
            .collect(),
        energies,
    }
}

/// `reference()` from `from` under `volts`, as a row.
fn reference_row(from: History, volts: f64) -> Row {
    let expected = reference(from, volts);
    let e = expected.state;
    Row {
        kelvin: e.t,
        current_a: e.i,
        potential: [0.5 * e.b, -0.5 * e.b],
        memory: (0..8)
            .map(|cell| {
                let sign = if cell < 4 { 1.0 } else { -1.0 };
                let [x, y] = DIRECTIONS[cell % 4].map(|d| sign * d);
                [x * e.h, y * e.h, x * e.z, y * e.z]
            })
            .collect(),
        energies: [
            expected.heat(),
            expected.stored_after - expected.stored_before,
            expected.numerical,
            expected.supplied,
        ],
    }
}

/// The history `reference()` steps from: triangle 0's field and retained
/// polarization along its direction (+x), and the flux density they give.
fn history(row: &Row) -> History {
    let (h, z) = (row.memory[0][0], row.memory[0][2]);
    History {
        t: row.kelvin,
        i: row.current_a,
        h,
        z,
        b: MU * h + z,
    }
}

fn kernel_run(system: &mut System, winding: NodeId, cells: NodeId) -> Vec<Row> {
    let mut rows = Vec::new();
    let mut before = [0.0; 4];
    with_cx(|cx| {
        for _ in 0..STEPS {
            let report = system.step(DT, cx).unwrap();
            assert_eq!(
                report.dt_s.to_bits(),
                DT.to_bits(),
                "the kernel took the whole interval"
            );
            let ledger = system.ledger();
            let node = &ledger.nodes[winding.index()];
            let now = [
                node.dissipated_j.value(),
                node.stored_change_j.value(),
                node.numerical_loss_j.value(),
                ledger
                    .boundaries
                    .iter()
                    .map(|b| b.energy_in_j.value())
                    .sum(),
            ];
            let kelvin = system
                .state::<TimelineState<'static>>(cells)
                .unwrap()
                .simulation()
                .temperature_k(CELL)
                .unwrap();
            let circuit = system.state::<PlayState>(winding).unwrap().circuit();
            let energies = [0, 1, 2, 3].map(|k| now[k] - before[k]);
            rows.push(row(kelvin, circuit, energies));
            before = now;
        }
    });
    rows
}

/// The module's run: `step_hysteretic_candidate` from `start`, every
/// conductor it reports held to `law` at its own end temperature, bit for
/// bit.
fn module_run(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
    start: &HystereticWindingCircuit<'_>,
    law: &CellConductors,
) -> Vec<Row> {
    let context = point(300.0);
    let cells = cells(card, curve, &context);
    let regions = routes(8);
    let network = layout(&cells, &regions);
    let (mut heat, mut state) = (thermal(curve, 2.0 * DEPTH), start.clone());
    let mut rows = Vec::new();
    for (index, volts) in DRIVES.into_iter().enumerate() {
        let next = network
            .step_hysteretic_candidate(
                &heat,
                (index + 1) as f64 * DT,
                &state,
                &[volts],
                controls(),
                winding_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let kelvin = next.thermal().temperature_k(CELL).unwrap();
        let actual = &next.report().electrical;
        let copper = actual.cells.iter().map(|loss| &loss.source);
        let core = actual.field.cells.iter().map(|loss| &loss.source);
        let sources: Vec<_> = copper
            .map(|s| (s.triangle, &s.conductor))
            .chain(core.map(|s| (s.triangle, &s.conductor)))
            .collect();
        assert_eq!(sources.len(), 8);
        for (triangle, conductor) in sources {
            assert_eq!(
                &law.conductor(triangle, kelvin).unwrap(),
                conductor,
                "step {index}, triangle {triangle}"
            );
        }
        let energies = [
            actual.heat_j,
            actual.stored_after_j - actual.stored_before_j,
            actual.numerical_dissipation_j,
            actual.supplied_work_j,
        ];
        rows.push(row(kelvin, next.circuit(), energies));
        (heat, state, _) = next.into_parts();
    }
    rows
}

/// The largest differences over the intervals: temperature \[K\], current
/// \[A\], potential \[Wb/m\], field \[A/m\], polarization \[T\] and energy
/// \[J\].
fn largest(ours: &[Row], theirs: &[Row]) -> [f64; 6] {
    assert_eq!((ours.len(), theirs.len()), (STEPS, STEPS));
    let mut most = [0.0_f64; 6];
    for (a, b) in ours.iter().zip(theirs) {
        let mut take = |slot: usize, x: f64, y: f64| most[slot] = most[slot].max((x - y).abs());
        take(0, a.kelvin, b.kelvin);
        take(1, a.current_a, b.current_a);
        for (x, y) in a.potential.iter().zip(b.potential) {
            take(2, *x, y);
        }
        assert_eq!(a.memory.len(), b.memory.len());
        for (x, y) in a.memory.iter().zip(&b.memory) {
            for (k, (p, q)) in x.iter().zip(y).enumerate() {
                take(3 + k / 2, *p, *q);
            }
        }
        for (x, y) in a.energies.iter().zip(b.energies) {
            take(5, *x, y);
        }
    }
    most
}

const UNITS: [&str; 6] = ["K", "A", "Wb/m", "A/m", "T", "J"];

fn hold(most: [f64; 6], bounds: [f64; 6], what: &str) {
    let shown: Vec<_> = most
        .iter()
        .zip(UNITS)
        .map(|(value, unit)| format!("{value:.3e} {unit}"))
        .collect();
    println!("{what}: {}", shown.join(", "));
    for ((value, bound), unit) in most.iter().zip(bounds).zip(UNITS) {
        assert!(*value <= bound, "{what}: {value:e} {unit} above {bound:e}");
    }
}

/// Bounds of the kernel against `reference()` and against the module, in
/// [`largest`]'s order. On the export of 187374db3b with the shared kernel
/// file, the kernel agrees with `reference()` within 1.408e-10 K, 1.355e-11
/// A, 5.049e-13 Wb/m, 4.809e-12 A/m, 9.618e-13 T and 7.449e-12 J, and with
/// the module within 2.171e-10 K and the same others to rounding: the
/// interface closes to a relative 1e-10. The bounds sit about fifty times
/// above.
const ORACLE: [f64; 6] = [7e-9, 7e-10, 2.5e-11, 2.5e-10, 5e-11, 4e-10];
const MODULE: [f64; 6] = [1.1e-8, 7e-10, 2.5e-11, 2.5e-10, 5e-11, 4e-10];

/// The fixture of
/// `voltage_reversal_solves_current_core_diffusion_memory_and_temperature_together`:
/// one cell holds every triangle and the graph resistor; the drive reverses.
#[test]
fn the_declared_graph_reproduces_the_voltage_reversal_and_its_memory() {
    // `reference()`'s own first interval is resolved, as the module's test
    // asserts: every heat and the numerical dissipation above 0.01 J, the
    // cell warming by more than 1 K, and the current at the start
    // temperature above the coupled one.
    {
        let (old, drive) = (START, DRIVES[0]);
        let expected = reference(old, drive);
        let e = expected.state;
        assert!(expected.core > 0.01 && expected.copper > 0.01 && expected.pinning > 0.01);
        assert!(expected.numerical > 0.01);
        assert!(e.t > old.t + 1.0);
        assert!(electrical(old, old.t, drive).i > e.i + 1e-3);
    }

    let card: &'static Arc<MaterialCard> = Box::leak(Box::new(card(Some(500.0))));
    let curve: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(caloric(card)));
    let owner: &'static PlanarHystereticMagnetostatics =
        Box::leak(Box::new(problem(card, paired_squares())));
    let start = circuit(card, owner);
    let law = Arc::new(CellConductors::new(
        card,
        point(300.0),
        start.field().conductors(),
        start.conductors(),
    ));
    let winding = PlanarWinding::at_port_temperatures(
        &start,
        vec![law.clone()],
        BTreeMap::new(),
        vec![HeatReceiver {
            triangles: (0..8).collect(),
            resistors: vec![2],
        }],
        winding_controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let cell = thermal(curve, 2.0 * DEPTH);
    let (mut system, electrical, cells) = graph(
        winding,
        PlanarWindingState::new(start.clone()),
        timeline(&cell, conduction(&cell, &[(MATERIAL, curve)], &[])),
        &[CELL],
        Some(drive),
    );
    let kernel = kernel_run(&mut system, electrical, cells);
    ledger_closes(&system, electrical, cells, 0.0);

    // Each interval against `reference()` from the kernel's own start of
    // interval.
    let mut from = START;
    let expected: Vec<Row> = kernel
        .iter()
        .zip(DRIVES)
        .map(|(row, volts)| {
            let expected = reference_row(from, volts);
            from = history(row);
            expected
        })
        .collect();
    let oracle = largest(&kernel, &expected);
    hold(oracle, ORACLE, "kernel against reference");

    // The module, step for step, its conductors held to the law first.
    let module = module_run(card, curve, &start, &law);
    let most = largest(&kernel, &module);
    hold(most, MODULE, "kernel against module");

    // The reversal is real: the cell warms by more than 40 K and the
    // retained polarization ends reversed.
    let end = &kernel[STEPS - 1];
    assert!(end.kelvin > 340.0 && end.memory[0][2] < -0.05, "{end:?}");
}
