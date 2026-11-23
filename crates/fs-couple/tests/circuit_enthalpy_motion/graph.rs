//! The moving field's coupling as a declared graph on the composition kernel.
//!
//! `circuit_enthalpy::motion::step_candidate` advances, by hand, a circuit
//! whose capacitor or inductor is a gap of variable width, the body that
//! moves the gap and the cells the circuit's resistors heat: the circuit,
//! the gap and the body in one discrete-gradient descriptor step
//! (`MaterialFieldMotion`), each resistor at the material resistance of its
//! cell's END temperature, and Newton on the cell energies. The graph says
//! the same with no code for the pair:
//!
//! - `fs_system::circuit::KirchhoffCircuit`: the module's netlist with the
//!   gap's branch replaced by a source port, each resistor following its
//!   card through [`CardResistance`] at the temperature its port receives;
//! - `fs_system::elements::GapField`: the gap's field, joined to that port;
//! - `fs_system::elements::Mass`, the module's body, and `Force`, its held
//!   external force, in one common-flow junction with the gap: the body
//!   sets the velocity at which the gap opens;
//! - the cells of [`timeline`] and one heat route per resistor.
//!
//! The module steps the circuit, gap and body as one descriptor; the graph's
//! components each take their own discrete-gradient step and the kernel's
//! interface solve joins them, so the two agree only within the solves'
//! tolerances.

use std::sync::Arc;

use fs_couple::circuit_enthalpy::motion::{self, FieldMotionEnthalpyInput};
use fs_couple::material_circuit::MaterialCircuit;
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor,
};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_system::circuit::KirchhoffCircuit;
use fs_system::elements::{Force, GapField, GapFieldState, Mass, Motion, SourceAccount, Waveform};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::TimelineState;

use crate::card::CardResistance;
use crate::common::*;
use crate::kernel::{CELLS, timeline, with_cx};

/// The gap's admitted range [m], the module's and the element's.
const LIMITS: [f64; 2] = [0.5, 1.5];
/// The cells' link conductance [W/K].
const CONDUCTANCE: f64 = 1.0;

/// Which storage the gap holds: both fixtures' cards have unit coefficient,
/// area, turns and length, so the coupling scale is 1 and the gap starts at
/// 1 m.
#[derive(Clone, Copy, Debug)]
enum Field {
    Capacitor,
    Inductor,
}

/// A moving-field device: the module's netlist, whose branch 0 is the gap's
/// capacitor or inductor and whose resistances are placeholders the cards
/// replace, the gap's charge \[C\] or flux linkage \[Wb\], the held source
/// ports \[V or A\], the body's mass \[kg\] and velocity \[m/s\], the held
/// force \[N\] toward a wider gap, and the run.
struct Fixture {
    netlist: CircuitGraph,
    field: Field,
    coordinate: f64,
    ports: Vec<f64>,
    mass_kg: f64,
    velocity_m_s: f64,
    force_n: f64,
    dt_s: f64,
    steps: usize,
}

/// One accepted interval: the gap's coordinate, its width \[m\], the body's
/// velocity \[m/s\], the cells' temperatures \[K\] and the Joule heat \[J\].
#[derive(Debug, Clone)]
struct Row {
    coordinate: f64,
    gap_m: f64,
    velocity_m_s: f64,
    kelvin: [f64; 2],
    joule_j: f64,
}

fn capacitor_storage() -> MaterialStorage {
    let card = crate::source::scalar_card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        1.0,
    );
    let capacitor = ResolvedCapacitor::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    MaterialStorage {
        branch: 0,
        element: MaterialStorageElement::Capacitor(capacitor),
    }
}

/// The module's run: `steps` intervals with `motion::step_candidate`, the
/// resistors resolved by the fixtures' `coefficients`.
fn module_run(
    fixture: &Fixture,
    cards: &[MaterialCard; 2],
    curves: &[ResolvedHeatCapacityCurve; 2],
) -> Vec<Row> {
    let mut thermal = initial(curves, CONDUCTANCE, None);
    let resistors = coefficients(&thermal, cards, curves).unwrap().resistors;
    let (fluxes, charges, storage) = match fixture.field {
        Field::Capacitor => (vec![], vec![fixture.coordinate], capacitor_storage()),
        Field::Inductor => (vec![fixture.coordinate], vec![], crate::gap_storage()),
    };
    let circuit = MaterialCircuit::new_with_storage(
        fixture.netlist.clone(),
        resistors,
        vec![storage],
        fluxes,
        charges,
        &fixture.ports,
    )
    .unwrap();
    let mut model = MaterialFieldMotion::new(
        circuit,
        0,
        crate::body(fixture.mass_kg),
        fixture.velocity_m_s,
        LIMITS,
    )
    .unwrap();
    assert_eq!(
        model.gap_m().to_bits(),
        1.0_f64.to_bits(),
        "the graph starts at 1 m"
    );
    let mut rows = Vec::new();
    for k in 1..=fixture.steps {
        let input = FieldMotionEnthalpyInput {
            motion: &model,
            ports: &fixture.ports,
            external_force_n: fixture.force_n,
            destinations: &DESTINATIONS,
        };
        let next = motion::step_candidate(
            &thermal,
            k as f64 * fixture.dt_s,
            input,
            config(),
            &[0.0; 2],
            |state| coefficients(state, cards, curves),
            || false,
        )
        .unwrap();
        let circuit = next.motion().circuit();
        rows.push(Row {
            coordinate: match fixture.field {
                Field::Capacitor => circuit.charges_c()[0],
                Field::Inductor => circuit.fluxes_wb()[0],
            },
            gap_m: next.motion().gap_m(),
            velocity_m_s: next.motion().velocity_m_s(),
            kelvin: CELLS.map(|id| next.thermal().temperature_k(id).unwrap()),
            joule_j: next.report().deposited_joule_j,
        });
        let (t, m, _) = next.into_parts();
        thermal = t;
        model = m;
    }
    rows
}

/// The graph's nodes: the circuit, the gap, the body, the force and the
/// cells.
#[derive(Clone, Copy)]
struct Nodes {
    circuit: NodeId,
    gap: NodeId,
    body: NodeId,
    push: NodeId,
    cells: NodeId,
}

/// The module's netlist as a `KirchhoffCircuit` whose gap branch is a
/// source port after the held ones, each resistor following its card, and
/// the gap's field.
fn kernel_circuit(
    fixture: &Fixture,
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
) -> (KirchhoffCircuit, GapField) {
    let gap_port = fixture.ports.len();
    let start = sourced(curves, CONDUCTANCE, None, [0.0; 2]);
    // The gap's port follows the held ports. A capacitor receives the
    // current the circuit draws through a voltage source; an inductor the
    // voltage across a current source.
    let mut netlist = fixture.netlist.clone();
    let (element, branch) = match fixture.field {
        Field::Capacitor => (
            GapField::capacitor(1.0, LIMITS),
            Branch::VoltageSource { port: gap_port },
        ),
        Field::Inductor => (
            GapField::inductor(1.0, LIMITS),
            Branch::CurrentSource { port: gap_port },
        ),
    };
    netlist.branches[0].2 = branch;
    // The module's placeholders are NaN; the kernel admits only resistances
    // it could step with, which the cards then replace.
    for (_, _, element) in &mut netlist.branches {
        if let Branch::Resistor { ohms } = element {
            *ohms = 1.0;
        }
    }
    let mut circuit = KirchhoffCircuit::new(netlist, vec![], 1e-8).unwrap();
    for (tag, (destination, (card, curve))) in
        (0_u64..).zip(DESTINATIONS.iter().zip(cards.iter().zip(curves.iter())))
    {
        circuit = circuit
            .following(
                destination.branch,
                Arc::new(CardResistance {
                    card,
                    curve,
                    volume_m3: start.cells()[&destination.cell].geometry.volume_m3,
                    tag,
                }),
            )
            .unwrap();
    }
    (circuit, element.unwrap())
}

/// The same device as a graph: the netlist with the gap's branch a source
/// port joined to a [`GapField`], the body and the force at a common
/// velocity with the gap's opening, the cells of [`timeline`] and one heat
/// route per resistor into its destination cell.
fn graph(
    fixture: &Fixture,
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
) -> (System, Nodes) {
    let (timeline, cells_start) = timeline(curves, CONDUCTANCE, [0.0; 2], fixture.dt_s);
    let gap_port = fixture.ports.len();
    let (circuit, element) = kernel_circuit(fixture, cards, curves);
    let state = circuit.state(&[], &[]).unwrap();

    let mut builder = SystemBuilder::new();
    let nodes = Nodes {
        circuit: builder.add("circuit", circuit, state).unwrap(),
        gap: builder
            .add_with_tangent(
                "gap",
                element,
                GapFieldState {
                    coordinate: fixture.coordinate,
                    gap_m: 1.0,
                },
            )
            .unwrap(),
        body: builder
            .add_with_tangent(
                "body",
                Mass::new(fixture.mass_kg).unwrap(),
                Motion {
                    position_m: 1.0,
                    velocity_m_s: fixture.velocity_m_s,
                },
            )
            .unwrap(),
        push: builder
            .add_with_tangent(
                "push",
                Force::new(Waveform::constant(fixture.force_n).unwrap()),
                SourceAccount::default(),
            )
            .unwrap(),
        cells: builder.add("cells", timeline, cells_start).unwrap(),
    };
    builder
        .bond(
            builder
                .port(nodes.circuit, &format!("source{gap_port}"))
                .unwrap(),
            builder.port(nodes.gap, "field").unwrap(),
        )
        .unwrap();
    builder
        .common_flow_ratios(
            builder.port(nodes.body, "p").unwrap(),
            &[
                (builder.port(nodes.gap, "gap").unwrap(), -1.0),
                (builder.port(nodes.push, "p").unwrap(), 1.0),
            ],
        )
        .unwrap();
    for (port, &value) in fixture.ports.iter().enumerate() {
        builder
            .boundary(
                builder
                    .port(nodes.circuit, &format!("source{port}"))
                    .unwrap(),
                Box::new(move |_, _, values| values[0] = value),
            )
            .unwrap();
    }
    heat_routes(&mut builder, nodes);
    let mut system = builder.seal().unwrap();
    for destination in DESTINATIONS {
        let port = system
            .port(nodes.circuit, &format!("heat{}", destination.branch))
            .unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    (system, nodes)
}

/// One heat route per resistor into its destination cell, and the cells'
/// environment to a sink at 300 K.
fn heat_routes(builder: &mut SystemBuilder, nodes: Nodes) {
    for destination in DESTINATIONS {
        builder
            .heat_route(
                builder
                    .port(nodes.circuit, &format!("heat{}", destination.branch))
                    .unwrap(),
                &[(
                    builder
                        .port(nodes.cells, &format!("cell{}", destination.cell.0))
                        .unwrap(),
                    1.0,
                )],
                None,
                RouteTemperature::SinkEquivalent,
            )
            .unwrap();
    }
    builder
        .heat_route(
            builder.port(nodes.cells, "environment").unwrap(),
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
}

fn graph_run(system: &mut System, nodes: Nodes, dt_s: f64, steps: usize) -> Vec<Row> {
    let mut rows = Vec::new();
    let mut joule_before = 0.0;
    with_cx(|cx| {
        for _ in 0..steps {
            let report = system.step(dt_s, cx).unwrap();
            assert_eq!(
                report.dt_s.to_bits(),
                dt_s.to_bits(),
                "the kernel took the whole interval"
            );
            let joule = system.ledger().nodes[nodes.circuit.index()]
                .dissipated_j
                .value();
            let gap = system.state::<GapFieldState>(nodes.gap).unwrap();
            let body = system.state::<Motion>(nodes.body).unwrap();
            let timeline = system.state::<TimelineState<'static>>(nodes.cells).unwrap();
            rows.push(Row {
                coordinate: gap.coordinate,
                gap_m: gap.gap_m,
                velocity_m_s: body.velocity_m_s,
                kelvin: CELLS.map(|id| timeline.simulation().temperature_k(id).unwrap()),
                joule_j: joule - joule_before,
            });
            joule_before = joule;
        }
    });
    rows
}

/// Bounds of one comparison: temperature \[K\], coordinate \[C or Wb\], gap
/// \[m\], velocity \[m/s\] and Joule heat per interval \[J\].
struct Bounds {
    kelvin: f64,
    coordinate: f64,
    gap_m: f64,
    velocity_m_s: f64,
    joule_j: f64,
}

/// Runs the fixture both ways and holds every interval to `bounds`; checks
/// the graph's energy account: the force's work and the sources' supply are
/// the gap's, the body's and the cells' energy change. Returns the module's
/// rows and the graph's.
fn compare(
    fixture: &Fixture,
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    bounds: &Bounds,
) -> [Vec<Row>; 2] {
    let module = module_run(fixture, cards, curves);
    let (mut system, nodes) = graph(fixture, cards, curves);
    let kernel = graph_run(&mut system, nodes, fixture.dt_s, fixture.steps);
    assert_eq!((kernel.len(), module.len()), (fixture.steps, fixture.steps));
    for (k, (ours, theirs)) in kernel.iter().zip(&module).enumerate() {
        for cell in 0..2 {
            let (x, y) = (ours.kelvin[cell], theirs.kelvin[cell]);
            assert!(
                (x - y).abs() <= bounds.kelvin,
                "interval {k}, cell {cell}: {x} vs {y}"
            );
        }
        for (what, x, y, bound) in [
            (
                "coordinate",
                ours.coordinate,
                theirs.coordinate,
                bounds.coordinate,
            ),
            ("gap", ours.gap_m, theirs.gap_m, bounds.gap_m),
            (
                "velocity",
                ours.velocity_m_s,
                theirs.velocity_m_s,
                bounds.velocity_m_s,
            ),
            ("Joule heat", ours.joule_j, theirs.joule_j, bounds.joule_j),
        ] {
            assert!((x - y).abs() <= bound, "interval {k}, {what}: {x} vs {y}");
        }
    }
    let ledger = system.ledger();
    assert!(
        ledger.energy_residual_j().abs() <= 1e-8,
        "{}",
        ledger.energy_residual_j()
    );
    let account = |node: NodeId| ledger.nodes[node.index()].stored_change_j.value();
    let pushed = -system
        .state::<SourceAccount>(nodes.push)
        .unwrap()
        .energy_in_j;
    let supplied: f64 = ledger
        .boundaries
        .iter()
        .map(|boundary| boundary.energy_in_j.value())
        .sum();
    let stored = account(nodes.gap) + account(nodes.body) + account(nodes.cells);
    assert!(
        (pushed + supplied - stored).abs() <= 1e-8,
        "{pushed} + {supplied} vs {stored}"
    );
    [module, kernel]
}

fn cards(
    first: (f64, f64, f64),
    second: (f64, f64, f64),
) -> (
    &'static [MaterialCard; 2],
    &'static [ResolvedHeatCapacityCurve; 2],
) {
    let (a, ca) = material(first.0, first.1, first.2);
    let (b, cb) = material(second.0, second.1, second.2);
    (Box::leak(Box::new([a, b])), Box::leak(Box::new([ca, cb])))
}

/// The gap and its body's velocity at the last interval, and the warmer
/// cell's rise \[K\].
fn reached(rows: &[Row]) -> (f64, f64, f64) {
    let last = rows.last().unwrap();
    (last.gap_m, last.velocity_m_s, last.kelvin[0] - 300.0)
}

#[test]
fn the_declared_graph_reproduces_a_gap_inductor_discharging_through_heated_resistors() {
    // The module's tests' finite coil: the gap's winding shorted through
    // both resistors, cell 30's resistance rising 0.5 % per kelvin.
    let (cards, curves) = cards((1.0, 1.0, 0.005), (3.0, 2.0, 0.0));
    let fixture = Fixture {
        netlist: CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Inductor { henries: f64::NAN }),
                (1, 2, Branch::Resistor { ohms: f64::NAN }),
                (2, 0, Branch::Resistor { ohms: f64::NAN }),
            ],
            transformers: vec![],
        },
        field: Field::Inductor,
        coordinate: 10.0,
        ports: vec![],
        mass_kg: 100.0,
        velocity_m_s: 0.0,
        force_n: 0.0,
        dt_s: 0.1,
        steps: 8,
    };
    // Measured: 1.2e-9 K, 1.9e-12 Wb, 3.4e-15 m, 1.0e-14 m/s, 1.9e-12 J.
    let bounds = Bounds {
        kelvin: 1e-8,
        coordinate: 2e-11,
        gap_m: 1e-13,
        velocity_m_s: 1e-13,
        joule_j: 2e-11,
    };
    let [module, _] = compare(&fixture, cards, curves, &bounds);
    let (gap, velocity, rise) = reached(&module);
    assert!(
        gap < 0.99 && velocity < -0.01 && rise > 5.0,
        "{gap} {velocity} {rise}"
    );
}

#[test]
fn the_declared_graph_reproduces_a_gap_capacitor_pushed_open_while_it_discharges() {
    let (cards, curves) = cards((1.0, 1.0, 0.005), (3.0, 2.0, 0.0));
    let fixture = Fixture {
        netlist: CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 2, Branch::Resistor { ohms: f64::NAN }),
                (2, 0, Branch::Resistor { ohms: f64::NAN }),
            ],
            transformers: vec![],
        },
        field: Field::Capacitor,
        coordinate: 4.0,
        ports: vec![],
        mass_kg: 20.0,
        velocity_m_s: 0.05,
        force_n: 10.0,
        dt_s: 0.1,
        steps: 8,
    };
    // Measured: 2.7e-11 K, 1.4e-14 C, 1.3e-14 m, 1.8e-14 m/s, 1.1e-13 J.
    let bounds = Bounds {
        kelvin: 5e-10,
        coordinate: 2e-13,
        gap_m: 2e-13,
        velocity_m_s: 2e-13,
        joule_j: 2e-12,
    };
    let [module, _] = compare(&fixture, cards, curves, &bounds);
    let (gap, velocity, rise) = reached(&module);
    assert!(
        gap > 1.0 && velocity > 0.05 && rise > 0.1,
        "{gap} {velocity} {rise}"
    );
}

#[test]
fn the_declared_graph_reproduces_a_held_linkage_beside_a_driven_branch() {
    // The module's closed-form fixture: the winding shorted by a 0 V port
    // keeps its linkage while a second port heats both resistors; its work
    // does not drive the body.
    let (cards, curves) = cards((2.0, 2.0, 0.005), (4.0, 3.0, 0.0));
    let fixture = Fixture {
        netlist: CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Inductor { henries: f64::NAN }),
                (2, 0, Branch::Resistor { ohms: f64::NAN }),
                (2, 0, Branch::Resistor { ohms: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
                (2, 0, Branch::VoltageSource { port: 1 }),
            ],
            transformers: vec![],
        },
        field: Field::Inductor,
        coordinate: 0.4,
        ports: vec![0.0, 3.0],
        mass_kg: 2.0,
        velocity_m_s: 0.02,
        force_n: 0.03,
        dt_s: 0.25,
        steps: 4,
    };
    // Measured: 1.3e-8 K, 7.0e-11 J; the coordinate, gap and velocity equal.
    let bounds = Bounds {
        kelvin: 2e-7,
        coordinate: 1e-14,
        gap_m: 1e-13,
        velocity_m_s: 1e-13,
        joule_j: 1e-9,
    };
    let [module, kernel] = compare(&fixture, cards, curves, &bounds);
    // The linkage is held and the force less the attraction `0.4^2 / 2`
    // accelerates the body uniformly, which the midpoint rule integrates
    // exactly: both sides within their solves' tolerances.
    let t = 0.25 * 4.0;
    let a = (0.03 - 0.5 * 0.4 * 0.4) / 2.0;
    for (side, rows, bound) in [("module", &module, 1e-10), ("graph", &kernel, 1e-10)] {
        let (gap, velocity, _) = reached(rows);
        let expected = 1.0 + 0.02 * t + 0.5 * a * t * t;
        assert!(
            (gap - expected).abs() <= bound,
            "{side}: {gap} vs {expected}"
        );
        assert!(
            (velocity - (0.02 + a * t)).abs() <= bound,
            "{side}: {velocity}"
        );
        for row in rows {
            assert!(
                (row.coordinate - 0.4).abs() <= 1e-14,
                "{side}: {}",
                row.coordinate
            );
        }
    }
}
