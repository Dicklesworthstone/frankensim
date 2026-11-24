//! The solid winding of solid.rs's first test as a declared graph on the
//! composition kernel.
//!
//! `WindingThermalNetwork::step_nonlinear_candidate` steps solid.rs's solid
//! winding with every conductor at its cell's END temperature, Newton on the
//! four cells' energies. The graph says the same with no code for the pair:
//!
//! - `PlanarWinding::at_port_temperatures` on solid.rs's circuit: one heat
//!   port per cell, receiving that cell's four conducting triangles, and one
//!   [`CellConductors`] law, so that each port's triangles follow the card
//!   at its own temperature and the solid regions' profiles follow their
//!   conductivity;
//! - the shared kernel's backward-Euler timeline over the four cells; cell
//!   0's authored 6 W heater enters through the timeline's environment port;
//! - the source held at 1 V.
//!
//! Held, in this order, to the ledger with the heater's heat, to solid.rs's
//! closed form on the kernel's own states, and to the module, step for step.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_em::winding::component::{ConductorLaw, HeatReceiver, PlanarWinding, PlanarWindingState};
use fs_em::winding::nonlinear::NonlinearWindingCircuit;
use fs_system::{NodeId, System};
use fs_thermal::component::TimelineState;

use super::graph::TOLERANCE_J;
use super::kernel::{CellConductors, graph, ledger_closes, with_cx};
use super::solid::{
    CELL_IDS, HEAT_CAPACITY, HEATER_W, Parts, SOLID_MATERIAL, STEP, VOLTAGE, circuit, closed,
    config, holds, parts, receiver, temperatures,
};
use super::*;

const STEPS: usize = 2;

/// An accepted state: the winding circuit and the four cells' temperatures
/// \[K\].
type State = (NonlinearWindingCircuit<'static>, [f64; 4]);

/// One receiver per cell, with the cell's conducting triangles.
fn receivers(circuit: &NonlinearWindingCircuit<'_>) -> Vec<HeatReceiver> {
    let mesh = circuit.field().problem().mesh();
    (0..CELL_IDS.len())
        .map(|cell| HeatReceiver {
            triangles: (0..mesh.faces().len())
                .filter(|&triangle| receiver(mesh, triangle) == Some(cell))
                .collect(),
            resistors: vec![],
        })
        .collect()
}

/// The kernel's accepted states: the start, then the end of each interval.
fn kernel_run(system: &mut System, winding: NodeId, cells: NodeId) -> Vec<State> {
    let read = |system: &System| -> State {
        let circuit = system
            .state::<PlanarWindingState<'static>>(winding)
            .unwrap()
            .circuit()
            .clone();
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        (circuit, temperatures(timeline.simulation()))
    };
    let mut states = vec![read(system)];
    with_cx(|cx| {
        for _ in 0..STEPS {
            let report = system.step(STEP, cx).unwrap();
            assert_eq!(
                report.dt_s.to_bits(),
                STEP.to_bits(),
                "the kernel took the whole interval"
            );
            states.push(read(system));
        }
    });
    states
}

/// The module's accepted states from `parts`, every conductor it reports
/// held to `law` at its cell's end temperature, bit for bit.
fn module_run(parts: &Parts, law: &CellConductors) -> Vec<State> {
    let layout = WindingThermalNetwork {
        cells: parts.cells,
        links: &parts.links,
        magnetic_regions: &parts.regions,
        graph_resistors: &[],
        volume_tolerance_m3: 1e-12,
    };
    let mesh = parts.problem.mesh();
    let mut thermal = parts.thermal.clone();
    let mut circuit = circuit(parts.problem);
    let mut states = vec![(circuit.clone(), temperatures(&thermal))];
    for interval in 1..=STEPS {
        let next = layout
            .step_nonlinear_candidate(
                &thermal,
                interval as f64 * STEP,
                &circuit,
                &[VOLTAGE],
                config(),
                field::controls(),
                &[0.0; 4],
                || false,
            )
            .unwrap();
        let t = temperatures(next.thermal());
        for loss in &next.report().electrical.field.cells {
            let triangle = loss.source.triangle;
            let cell = receiver(mesh, triangle).unwrap();
            assert_eq!(
                law.conductor(triangle, t[cell]).unwrap(),
                loss.source.conductor,
                "interval {interval}, triangle {triangle}"
            );
        }
        (thermal, circuit, _) = next.into_parts();
        states.push((circuit.clone(), t));
    }
    states
}

/// Each interval of `states` against solid.rs's closed form: its identities
/// on the accepted states, and the largest departure of a cell's stored
/// change from its closed-form heat, with the heater's on cell 0 \[J\].
fn closed_form(states: &[State]) -> f64 {
    let mut most = 0.0_f64;
    for ((old, before), (next, after)) in states.iter().zip(&states[1..]) {
        let closed = closed(old, next, *after);
        holds(&closed, old, next);
        for cell in 0..CELL_IDS.len() {
            let heater = if cell == 0 { STEP * HEATER_W } else { 0.0 };
            let stored = HEAT_CAPACITY * (after[cell] - before[cell]);
            most = most.max((stored - closed.cell_heat[cell] - heater).abs());
        }
        // The heater is real: its cell runs ahead of its neighbour.
        assert!(after[0] > after[1] + 0.5, "{after:?}");
    }
    most
}

/// The largest differences over the intervals' ends: temperature \[K\],
/// current \[A\], potential \[Wb/m\] and profile \[1/m²\].
fn largest(ours: &[State], theirs: &[State]) -> [f64; 4] {
    assert_eq!((ours.len(), theirs.len()), (STEPS + 1, STEPS + 1));
    let mut most = [0.0_f64; 4];
    for ((a, ta), (b, tb)) in ours.iter().zip(theirs).skip(1) {
        let mut take = |slot: usize, x: &[f64], y: &[f64]| {
            assert_eq!(x.len(), y.len());
            for (p, q) in x.iter().zip(y) {
                most[slot] = most[slot].max((p - q).abs());
            }
        };
        take(0, ta, tb);
        take(1, a.currents_a(), b.currents_a());
        take(2, a.field().potential_wb_m(), b.field().potential_wb_m());
        take(3, &a.profiles_per_m2()[0], &b.profiles_per_m2()[0]);
    }
    most
}

/// Bounds of [`closed_form`] and [`largest`], about fifty times the
/// departures measured.
const CLOSED_J: f64 = 2.5e-12;
const MODULE: [f64; 4] = [4e-11, 2.5e-11, 2.7e-12, 4e-13];

/// The fixture of `spatial_solid_heat_and_profiles_follow_the_same_finite_temperature_endpoint`:
/// four cells, no link, the heater on cell 0.
#[test]
fn the_declared_graph_reproduces_the_solid_winding_and_its_heater() {
    let parts = parts(0.0, false);
    let start = circuit(parts.problem);
    let law = Arc::new(CellConductors::new(
        parts.card,
        parts.context.clone(),
        start.field().conductors(),
        start.conductors(),
    ));
    let shared: Arc<dyn ConductorLaw> = law.clone();
    let winding = PlanarWinding::at_port_temperatures(
        &start,
        vec![shared; CELL_IDS.len()],
        BTreeMap::new(),
        receivers(&start),
        field::controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let (mut system, electrical, cells) = graph(
        winding,
        PlanarWindingState::new(start),
        (&parts.thermal, &[(SOLID_MATERIAL, parts.curve)], &[]),
        &CELL_IDS,
        Some(|_| VOLTAGE),
    );
    let kernel = kernel_run(&mut system, electrical, cells);
    // The heater's heat is its authored source over the run, not a reading.
    ledger_closes(&system, electrical, cells, HEATER_W * STEP * STEPS as f64);

    let departure = closed_form(&kernel);
    println!("graph S1: kernel against the closed form: {departure:.3e} J");
    assert!(departure <= CLOSED_J, "{departure:e} J above {CLOSED_J:e}");

    let module = module_run(&parts, &law);
    let most = largest(&kernel, &module);
    println!(
        "graph S1: kernel against module: {:.3e} K, {:.3e} A, {:.3e} Wb/m, {:.3e} 1/m2",
        most[0], most[1], most[2], most[3]
    );
    for (value, bound) in most.iter().zip(MODULE) {
        assert!(*value <= bound, "{most:?} against {MODULE:?}");
    }
}
