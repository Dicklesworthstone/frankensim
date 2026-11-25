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
//!
//! The second test declares solid.rs's second, `MaterialHeatNetwork`'s
//! assembly with sixteen magnetic laws, on the magnetic lane:
//! `PlanarWinding::magnetic_at_port_temperatures` with the same laws, and
//! link 17 in the timeline. Its ports return `-T1 (s1 - s0)` beside the
//! Joule heat where the module deposits `-P - (T1 s1 - T0 s0)`; the
//! Component books the difference, `-(P + (T1 - T0) s0)`, as numerical
//! loss, as in graph_caloric.rs.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::eddy::thermal::MagneticThermalLaw;
use fs_couple::circuit_enthalpy::material::MaterialNetwork;
use fs_em::winding::component::{
    ConductorLaw, HeatReceiver, MagneticLaw, MagneticLaws, PlanarWinding, PlanarWindingState,
};
use fs_em::winding::nonlinear::{NonlinearWindingCircuit, NonlinearWindingStep};
use fs_rep_mesh::TriComplex2;
use fs_system::{NodeId, System};
use fs_thermal::component::TimelineState;
use fs_thermal::enthalpy::conduction::ConductiveLink;

use super::graph::TOLERANCE_J;
use super::kernel::{CellConductors, conduction, graph, ledger_closes, timeline, with_cx};
use super::solid::{
    CELL_IDS, END_H, HEAT_CAPACITY, HEATER_W, Parts, SOLID_MATERIAL, STEP, VOLTAGE, check_receipt,
    circuit, closed, config, holds, linear, magnetic_laws, parts, receiver, temperatures,
};
use super::*;

const STEPS: usize = 2;

/// An accepted state: the winding circuit and the four cells' temperatures
/// \[K\].
type State = (NonlinearWindingCircuit<'static>, [f64; 4]);

/// An interval's step as the kernel's state keeps it, and the winding's
/// booked stored change and numerical loss over the interval \[J\].
type Kept = (NonlinearWindingStep, [f64; 2]);

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
fn kernel_run(system: &mut System, winding: NodeId, cells: NodeId) -> (Vec<State>, Vec<Kept>) {
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
    let (mut kept, mut before) = (Vec::new(), [0.0; 2]);
    with_cx(|cx| {
        for _ in 0..STEPS {
            let report = system.step(STEP, cx).unwrap();
            assert_eq!(
                report.dt_s.to_bits(),
                STEP.to_bits(),
                "the kernel took the whole interval"
            );
            states.push(read(system));
            let node = &system.ledger().nodes[winding.index()];
            let now = [node.stored_change_j.value(), node.numerical_loss_j.value()];
            let state = system.state::<PlanarWindingState<'static>>(winding);
            let step = state.unwrap().step().unwrap().clone();
            kept.push((step, [now[0] - before[0], now[1] - before[1]]));
            before = now;
        }
    });
    (states, kept)
}

/// The kernel's kept steps against solid.rs's closed forms, interval by
/// interval, with `check_receipt`'s record checks: each triangle's
/// resistivity, card, query temperature, field, current density, volume,
/// heat and power; the region currents; the voltages; and the supplied
/// work's split. The owner's numerical dissipation stays above 1e-5 J.
/// Then the field's internal energy `psi + T s` in closed form from the
/// step's flux densities, with the end turns': its change against the
/// winding's booked stored change, and in the balance of the supplied work
/// and the heater against it, the cells' heat and the booked numerical
/// loss \[J\].
fn records(card: &MaterialCard, kernel: &[State], kept: &[Kept], slope: f64, what: &str) {
    let mut most = [0.0_f64; 2];
    let intervals = kernel.iter().zip(&kernel[1..]).zip(kept);
    for (((old, before), (next, after)), (step, [stored, numerical])) in intervals {
        check_receipt(card, old, next, *after, step);
        assert!(step.numerical_dissipation_j > 1e-5, "{what}");
        let mesh = next.field().problem().mesh();
        let mut internal = [0.0; 2];
        for magnetic in &step.field.magnetic_cells {
            let owner = receiver(mesh, magnetic.triangle);
            let t = owner.map_or([300.0; 2], |cell| [before[cell], after[cell]]);
            let local = if owner.is_some() { slope } else { 0.0 };
            let volume = mesh.face_measure(magnetic.triangle).unwrap();
            for (i, b) in magnetic.flux_density_t.iter().enumerate() {
                let (mu, b2) = (1.0 + local * (t[i] - 300.0), b[0] * b[0] + b[1] * b[1]);
                let entropy = 0.5 * volume * b2 * local / (mu * mu);
                internal[i] += 0.5 * volume * b2 / mu + t[i] * entropy;
            }
        }
        let end = 0.5 * END_H * (next.currents_a()[0].powi(2) - old.currents_a()[0].powi(2));
        let change = internal[1] - internal[0] + end;
        let cells: f64 = (0..CELL_IDS.len())
            .map(|cell| HEAT_CAPACITY * (after[cell] - before[cell]))
            .sum();
        let work = STEP * (VOLTAGE * next.currents_a()[0] + HEATER_W);
        most[0] = most[0].max((stored - change).abs());
        most[1] = most[1].max((change + cells + numerical - work).abs());
    }
    let shown = most.map(|x| format!("{x:.3e}")).join(" and ");
    println!("graph {what}: internal energy against the booked change and the balance: {shown} J");
    for (value, bound) in most.iter().zip(RECORDS_J) {
        assert!(*value <= bound, "{what}: {most:?} against {RECORDS_J:?}");
    }
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

/// A triangle's flux density \[T\] from a P1 potential \[Wb/m\]: the
/// potential's gradient turned a quarter.
pub(super) fn flux_density(mesh: &TriComplex2, potential: &[f64], triangle: usize) -> [f64; 2] {
    let face = mesh.faces()[triangle].map(|v| v as usize);
    let [p, q, r] = face.map(|v| mesh.vertices()[v]);
    let [a, b, c] = face.map(|v| potential[v]);
    let twice_area = (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]);
    let dx = (a * (q[1] - r[1]) + b * (r[1] - p[1]) + c * (p[1] - q[1])) / twice_area;
    let dy = (a * (r[0] - q[0]) + b * (p[0] - r[0]) + c * (q[0] - p[0])) / twice_area;
    [dy, -dx]
}

/// Each interval of `states` against solid.rs's closed form: its identities
/// on the accepted states, and the largest departure of a cell's stored
/// change from its closed-form heat \[J\]. Beside the Joule heat that is the
/// heater's on cell 0, link 17's conduction between cells 0 and 1 when
/// `linked`, and each magnetic triangle's port heat `-T1 (s1 - s0)` at
/// `slope` per kelvin, `s = V B^2 slope / (2 mu^2)`. Also, over the run, the
/// terms `-(P + (T1 - T0) s0)` that the Component books, `P` the owner's
/// parameter work at the retained B, and the owner's backward-Euler
/// dissipation \[J\].
fn closed_form(states: &[State], slope: f64, linked: bool) -> [f64; 3] {
    let (mut most, mut terms, mut dissipation) = (0.0_f64, 0.0, 0.0);
    for ((old, before), (next, after)) in states.iter().zip(&states[1..]) {
        let closed = closed(old, next, *after);
        holds(&closed, old, next);
        let mut other = [STEP * HEATER_W, 0.0, 0.0, 0.0];
        if linked {
            let transfer = STEP * 0.5 * (after[0] - after[1]);
            other[0] -= transfer;
            other[1] += transfer;
        }
        dissipation += 0.5 * END_H * (next.currents_a()[0] - old.currents_a()[0]).powi(2);
        let mesh = next.field().problem().mesh();
        for triangle in 0..mesh.faces().len() {
            let b = [old, next].map(|c| flux_density(mesh, c.field().potential_wb_m(), triangle));
            let volume = mesh.face_measure(triangle).unwrap();
            let owner = receiver(mesh, triangle);
            let local = if owner.is_some() { slope } else { 0.0 };
            let t = owner.map_or([300.0; 2], |cell| [before[cell], after[cell]]);
            let mu = t.map(|t| 1.0 + local * (t - 300.0));
            let b2 = b.map(|b| b[0] * b[0] + b[1] * b[1]);
            let entropy = [0, 1].map(|i| 0.5 * volume * b2[i] * local / mu[i].powi(2));
            let change = (b[1][0] - b[0][0]).powi(2) + (b[1][1] - b[0][1]).powi(2);
            dissipation += 0.5 * volume * change / mu[1];
            if let Some(cell) = owner {
                let parameter = 0.5 * volume * b2[0] * (1.0 / mu[1] - 1.0 / mu[0]);
                other[cell] -= t[1] * (entropy[1] - entropy[0]);
                terms -= parameter + (t[1] - t[0]) * entropy[0];
            }
        }
        for cell in 0..CELL_IDS.len() {
            let stored = HEAT_CAPACITY * (after[cell] - before[cell]);
            most = most.max((stored - closed.cell_heat[cell] - other[cell]).abs());
        }
        // The heater is real: its cell runs ahead of its neighbour.
        assert!(after[0] > after[1] + 0.5, "{after:?}");
    }
    [most, terms, dissipation]
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
/// Bounds of [`records`]: the internal energy's change against the booked
/// stored change, and the balance \[J\]. Over S1 and S2's three runs the
/// kernel departs by at most 3.5e-18 J and 1.3e-13 J; the bounds sit about
/// fifty times above.
const RECORDS_J: [f64; 2] = [1.8e-16, 6.5e-12];
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
        timeline(
            &parts.thermal,
            conduction(&parts.thermal, &[(SOLID_MATERIAL, parts.curve)], &[]),
        ),
        &CELL_IDS,
        Some(|_| VOLTAGE),
    );
    let (kernel, kept) = kernel_run(&mut system, electrical, cells);
    // The heater's heat is its authored source over the run, not a reading.
    ledger_closes(&system, electrical, cells, HEATER_W * STEP * STEPS as f64);
    records(parts.card, &kernel, &kept, 0.0, "S1");

    let [departure, ..] = closed_form(&kernel, 0.0, false);
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

/// The module's accepted states from `parts` through `MaterialHeatNetwork`'s
/// assembly with `laws`, every conductor it reports held to `law` at its
/// cell's end temperature, bit for bit.
fn heat_module_run(parts: &Parts, law: &CellConductors, laws: &[MagneticThermalLaw]) -> Vec<State> {
    let bulk = MaterialNetwork {
        cells: parts.cells,
        resistors: &[],
        links: &parts.links,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let mesh = parts.problem.mesh();
    let mut thermal = parts.thermal.clone();
    let mut circuit = circuit(parts.problem);
    let mut states = vec![(circuit.clone(), temperatures(&thermal))];
    with_cx(|cx| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        for interval in 1..=STEPS {
            let next = prepared
                .step_nonlinear_winding_candidate(
                    cx,
                    &thermal,
                    interval as f64 * STEP,
                    &circuit,
                    &[VOLTAGE],
                    &parts.regions,
                    &[],
                    1e-12,
                    Some((laws, linear::limits())),
                    config(),
                    field::controls(),
                    &[0.0; 4],
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
    });
    states
}

/// Bounds of the second test, about fifty times the rounding measured over
/// the three variants: [`closed_form`]'s departure and the booked loss, less
/// the owner's dissipation, against the terms \[J\].
const MAGNETIC_J: [f64; 2] = [2.5e-12, 9e-12];

/// At slope 0, the same: [`largest`] and the cells' energy gap \[J\].
const CONTROL: [f64; 5] = [5e-10, 2.5e-11, 2.7e-12, 3.7e-12, 4.3e-10];

/// At slope 0.004, [`largest`] and the gap's departure from the terms as a
/// fraction of them. They measure 1.274e-8 K, 5.285e-12 A, 9.837e-13 Wb/m,
/// 1.644e-10 1/m2 and 1.1e-3, all from the cooler cells' feedback through
/// their resistivity and permeability, as in graph_caloric.rs; the bounds
/// sit about a third above.
const FEEDBACK: [f64; 5] = [1.6e-8, 7e-12, 1.3e-12, 2.1e-10, 1.5e-3];

/// The fixture of `full_heat_assembly_preserves_solid_loss_and_magnetic_caloric_owners`
/// in its three variants: no link at slope 0, link 17 at slope 0, and link
/// 17 at slope 0.004 per kelvin.
#[test]
fn the_declared_graph_and_the_heat_assembly_differ_by_the_booked_term() {
    let mut ends = Vec::new();
    for (linked, slope) in [(false, 0.0), (true, 0.0), (true, 0.004)] {
        let parts = parts(slope, linked);
        let start = circuit(parts.problem);
        let law = Arc::new(CellConductors::new(
            parts.card,
            parts.context.clone(),
            start.field().conductors(),
            start.conductors(),
        ));
        let laws = magnetic_laws(parts.card, parts.problem);
        let magnetic = MagneticLaws {
            triangles: laws
                .iter()
                .map(|binding| MagneticLaw {
                    triangle: binding.triangle,
                    law: binding.law.clone(),
                })
                .collect(),
            assembly: linear::limits(),
        };
        let shared: Arc<dyn ConductorLaw> = law.clone();
        let winding = PlanarWinding::magnetic_at_port_temperatures(
            &start,
            vec![shared; CELL_IDS.len()],
            BTreeMap::new(),
            magnetic,
            receivers(&start),
            field::controls(),
            TOLERANCE_J,
        )
        .unwrap();
        // Link 17 conducts as the module's `ConductionLink` does: the card's
        // 1 W/(m K) over 0.25 m to each cell's centre.
        let links: Vec<_> = parts
            .links
            .iter()
            .map(|link| ConductiveLink {
                link_id: link.link,
                half_distances_m: [link.distance_a_m, link.distance_b_m],
                contact_resistance_m2_k_w: 0.0,
            })
            .collect();
        let (mut system, electrical, cells) = graph(
            winding,
            PlanarWindingState::at_temperatures(start, vec![300.0; CELL_IDS.len()]),
            timeline(
                &parts.thermal,
                conduction(&parts.thermal, &[(SOLID_MATERIAL, parts.curve)], &links),
            ),
            &CELL_IDS,
            Some(|_| VOLTAGE),
        );
        let (kernel, kept) = kernel_run(&mut system, electrical, cells);
        ledger_closes(&system, electrical, cells, HEATER_W * STEP * STEPS as f64);
        let node = &system.ledger().nodes[electrical.index()];
        assert_eq!(node.entropy.as_ref().unwrap().audited_steps, STEPS as u64);
        let what = format!("S2 (link {linked}, slope {slope})");
        records(parts.card, &kernel, &kept, slope, &what);

        let [departure, terms, dissipation] = closed_form(&kernel, slope, linked);
        let booked = (node.numerical_loss_j.value() - dissipation - terms).abs();
        println!(
            "graph S2 (link {linked}, slope {slope}): closed form {departure:.3e} J; booked against the terms \
             {booked:.3e} J; terms {terms:.4e} J, dissipation {dissipation:.4e} J"
        );
        for (value, bound) in [departure, booked].into_iter().zip(MAGNETIC_J) {
            assert!(value <= bound, "{value:e} above {bound:e}");
        }

        ends.push(kernel[STEPS].0.currents_a()[0]);
        let module = heat_module_run(&parts, &law, &laws);
        let most = largest(&kernel, &module);
        let gap: f64 = (0..CELL_IDS.len())
            .map(|cell| HEAT_CAPACITY * (module[STEPS].1[cell] - kernel[STEPS].1[cell]))
            .sum();
        println!(
            "graph S2 (link {linked}, slope {slope}): kernel against module: {:.3e} K, {:.3e} A, {:.3e} Wb/m, \
             {:.3e} 1/m2; cells' gap {gap:.4e} J",
            most[0], most[1], most[2], most[3]
        );
        if slope == 0.0 {
            assert_eq!(terms.to_bits(), 0.0_f64.to_bits());
            let values = [most[0], most[1], most[2], most[3], gap.abs()];
            for (value, bound) in values.into_iter().zip(CONTROL) {
                assert!(value <= bound, "slope 0: {values:?} against {CONTROL:?}");
            }
        } else {
            // The term is resolved, two orders above the bound of its
            // booking; the module's cells end cooler by it, less the feedback
            // of their resistivity and permeability.
            assert!(-terms > 100.0 * MAGNETIC_J[1], "{terms:e}");
            let fraction = (gap - terms).abs() / terms.abs();
            println!(
                "graph S2 (slope {slope}): the cells' gap departs from the terms by {fraction:.3e} of them"
            );
            let values = [most[0], most[1], most[2], most[3], fraction];
            for (value, bound) in values.into_iter().zip(FEEDBACK) {
                assert!(
                    value <= bound,
                    "slope {slope}: {values:?} against {FEEDBACK:?}"
                );
            }
        }
    }
    // Link 17 and the slope each move the current, below `largest`'s
    // bounds, as S2 asserts of the module.
    assert!((ends[1] - ends[0]).abs() > 1e-9, "link: {ends:?}");
    assert!((ends[2] - ends[1]).abs() > 1e-7, "slope: {ends:?}");
}
