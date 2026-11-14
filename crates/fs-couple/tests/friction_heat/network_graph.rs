//! The friction-heat network case: two contacts that share a cell, declared
//! as a graph.
//!
//! `fs_couple::friction_heat::network` steps several friction interfaces on
//! one finite enthalpy network. Each contact has its own massless plane, and
//! all of their powers enter one backward-Euler residual, so contacts may
//! share a cell. Here the network of its own shared-cell test
//! (`rigid_friction_heat_network.rs`,
//! `g1_zero_friction_shared_cell_conduction_and_duplicate_route_refusal`) is
//! declared from parts: three cells at 330, 300 and 280 K, joined by contact
//! links 10 (cells 0 and 1) and 20 (cells 1 and 2), so that the two contacts
//! share cell 1. Nothing is written for the pair:
//!
//! - fs-thermal's `ThermalTimeline` of the three cells, stepped by its
//!   owner's backward-Euler solve, with both links named as carried outside
//!   and each cell's port returning its end-of-step temperature;
//! - two of the kernel's `HeatNode`s, one plane per contact, each joined to
//!   its two cells by its route's two half conductances;
//! - two of fs-tribo's `SlidingInterface`s under prescribed slip, one on
//!   each plane.
//!
//! One kernel step holds both interfaces. Cell 1's port is one junction with
//! a body port of each plane, so the timeline's one solve takes both planes'
//! heat into cell 1 together, and the timeline, both planes and both
//! interfaces are one loop whose wires the kernel's Newton closes together.
//!
//! The judge is the module. Every kernel step's endpoint goes to
//! `FrictionHeatNetworkStep::evaluate`: the three cells' energy increments,
//! the temperature each interface was given, and each interface's dissipated
//! work as the kernel's ledger booked it from the interface's receipt. The
//! step must then be accepted at the module tests' tolerances (1e-10 J). The
//! module forms every cell's residual and each interface's balance with its
//! own route conductances.
//!
//! The graph has no equivalent of: conductances that follow the materials at
//! the endpoint (the planes' are constants), a contact resistance, the
//! module's entropy audit, or other sources on the network.

use fs_couple::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use fs_exec::{CancelGate, Cx};
use fs_system::elements::HeatNode;
use fs_system::{NodeId, PortRef, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductivityLaw};
use fs_tribo::tangential::component::{PrescribedSlip, SlidingInterface};
use fs_tribo::{FrictionLaw, TriboError};

use super::graph::{Conducting, LawAt, coulomb, mu_falling, mu_three_tenths, with_cx};
use super::*;

/// The shared-cell test's cells: heat capacities [J/K] (unit specific heat
/// on the fixture's `T = 300 K + h`), conductivities [W/(m K)], the half
/// length of a cell [m] and the contact area [m2].
const CAP: [f64; 3] = [0.02, 0.03, 0.04];
const K: [f64; 3] = [40.0, 60.0, 100.0];
const A: f64 = 0.01;
const AREA: f64 = core::f64::consts::PI * 0.02 * 0.02;
/// The contact links: 10 joins cells 0 and 1, 20 joins cells 1 and 2.
const LINKS: [u64; 2] = [10, 20];
const START_K: [f64; 3] = [330.0, 300.0, 280.0];
const DT: f64 = 0.01;

/// The shared-cell test's three cells at `kelvin`, with its two contact
/// links and no film.
fn three_cells(kelvin: [f64; 3]) -> DepositionSimulation<'static> {
    let law: &'static dyn ThermalLaw = &SENSIBLE;
    let mut events = Vec::new();
    for (i, (&mass_kg, &t)) in CAP.iter().zip(&kelvin).enumerate() {
        let id = CellId(i as u64);
        let z = 2.0 * A * i as f64;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: id.0,
                body_id: id.0,
                material_id: MaterialId(1),
                mass_kg,
                specific_h_j_kg: t - 300.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, z - 2.0 * A],
                    end_m: [0.0, 0.0, z],
                    volume_m3: AREA * 2.0 * A,
                    surface_area_m2: 2.0 * AREA + 2.0 * core::f64::consts::PI * 0.02 * 2.0 * A,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for (i, &id) in LINKS.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                area_m2: AREA,
                conductance_w_k: series_conductance_w_k(AREA, [A; 2], [K[i], K[i + 1]]).unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

/// The timeline of the three cells on the contract: both links named as
/// carried outside, so the owner conducts nothing between the cells itself;
/// a heat port on every cell, returning its end-of-step temperature. With
/// no link conducted inside, the conductivity bound here is never asked for
/// a flow.
fn timeline(thermal: &DepositionSimulation<'static>) -> ThermalTimeline<'static> {
    let law: &'static Conducting = Box::leak(Box::new(Conducting { w_m_k: K[1] }));
    let conduction = ConductionObserver::new(
        thermal,
        &[(MaterialId(1), law as &dyn ConductivityLaw)],
        &[],
        &[],
        || Ok(()),
    )
    .expect("conduction")
    .with_links_carried_outside(thermal, &LINKS, || Ok(()))
    .expect("both contacts carried outside");
    let solver = ImplicitConductionGraph::new(
        thermal,
        &conduction,
        &[(MaterialId(1), law as &dyn ImplicitConductionLaw)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .expect("solver");
    let controls = ImplicitConductionControls {
        maximum_dt_s: DT,
        cell_tolerance_j: 1e-13,
        total_tolerance_j: 1e-12,
        ..ImplicitConductionControls::default()
    };
    ThermalTimeline::implicit(
        (0..3).map(CellId).collect(),
        solver,
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .expect("timeline")
    .returning_end_temperatures()
    .expect("end temperatures")
}

fn frictionless(_: f64) -> Result<FrictionLaw, TriboError> {
    coulomb(0.0)
}

/// The declared network and what the module needs to judge its steps.
struct Network {
    system: System,
    cells: NodeId,
    /// The plane each interface heats: two planes, or one plane twice.
    planes: [NodeId; 2],
    interfaces: [NodeId; 2],
    routes: [FrictionHeatRoute; 2],
}

/// What the module made of one kernel step.
struct Judged {
    report: FrictionHeatNetworkReport,
    /// The largest of the module's residuals at the graph's endpoint [J].
    residual_j: f64,
    /// Each interface's dissipated work over the step, from the ledger [J].
    dissipated_j: [f64; 2],
    /// The temperature each interface was given [K].
    given_k: [f64; 2],
    /// The kernel's Newton iterations in the step.
    iterations: u32,
}

fn port(builder: &SystemBuilder, node: NodeId, name: &str) -> PortRef {
    builder.port(node, name).expect("port")
}

fn cell(builder: &SystemBuilder, cells: NodeId, i: usize) -> PortRef {
    port(builder, cells, &format!("cell{i}"))
}

/// The module's rule: each contact heats its own plane, joined to its two
/// cells by its route's halves `halves[i]`. Cell 1 is one junction with a
/// body port of each plane, so the timeline's solve takes both planes'
/// heat into it in the same step.
fn a_plane_each(
    builder: &mut SystemBuilder,
    cells: NodeId,
    interfaces: [NodeId; 2],
    halves: [[f64; 2]; 2],
) -> [NodeId; 2] {
    let planes = [0, 1].map(|i| {
        let plane = HeatNode::new(&halves[i]).expect("plane");
        builder
            .add_with_tangent(&format!("plane{}", LINKS[i]), plane, ())
            .expect("plane")
    });
    for i in 0..2 {
        builder
            .common_effort(
                port(builder, planes[i], "source"),
                &[port(builder, interfaces[i], "heat")],
            )
            .expect("an interface on its plane");
    }
    builder
        .common_effort(
            cell(builder, cells, 0),
            &[port(builder, planes[0], "body0")],
        )
        .expect("cell 0 on plane 10");
    builder
        .common_effort(
            cell(builder, cells, 1),
            &[
                port(builder, planes[0], "body1"),
                port(builder, planes[1], "body0"),
            ],
        )
        .expect("cell 1 on both planes");
    builder
        .common_effort(
            cell(builder, cells, 2),
            &[port(builder, planes[1], "body1")],
        )
        .expect("cell 2 on plane 20");
    planes
}

/// Not the module's rule: both contacts heat ONE plane, joined to cell 0,
/// to cell 1 by both contacts' halves together, and to cell 2.
fn one_plane_for_both(
    builder: &mut SystemBuilder,
    cells: NodeId,
    interfaces: [NodeId; 2],
    halves: [[f64; 2]; 2],
) -> [NodeId; 2] {
    let plane =
        HeatNode::new(&[halves[0][0], halves[0][1] + halves[1][0], halves[1][1]]).expect("plane");
    let plane = builder.add_with_tangent("plane", plane, ()).expect("plane");
    builder
        .common_effort(
            port(builder, plane, "source"),
            &[
                port(builder, interfaces[0], "heat"),
                port(builder, interfaces[1], "heat"),
            ],
        )
        .expect("both interfaces on the plane");
    for i in 0..3 {
        builder
            .common_effort(
                cell(builder, cells, i),
                &[port(builder, plane, &format!("body{i}"))],
            )
            .expect("a cell on the plane");
    }
    [plane; 2]
}

/// The network with interface laws `laws` sliding at `slip_m_s` under the
/// tests' 10 N. With `one_plane`, both contacts heat ONE plane, joined to
/// cell 0, cell 1 (both contacts' halves) and cell 2; that is not the
/// module's rule, which gives each contact its own plane.
fn declare(laws: [LawAt; 2], slip_m_s: [f64; 2], one_plane: bool) -> Network {
    let thermal = three_cells(START_K);
    let routes = [0, 1].map(|i| {
        FrictionHeatRoute::from_link(&thermal, LINKS[i], [K[i], K[i + 1]]).expect("route")
    });
    // Each route's two half conductances, A (k / d), in its own arithmetic.
    let halves = [0, 1].map(|i| {
        let (area, d) = (routes[i].contact_area_m2(), routes[i].distances_m());
        [area * (K[i] / d[0]), area * (K[i + 1] / d[1])]
    });
    let timeline = timeline(&thermal);
    let start = timeline.start(&thermal).expect("start");

    let mut builder = SystemBuilder::new();
    let cells = builder.add("cells", timeline, start).expect("cells");
    let interfaces = [0, 1].map(|i| {
        let interface = SlidingInterface::new(
            laws[i],
            PrescribedSlip {
                normal_n: input(DT).normal_n,
                slip_rate_m_s: [slip_m_s[i], 0.0],
                nominal_pressure_pa: input(DT).nominal_pressure_pa,
                solve: input(DT).solve,
                work_tolerance_j: 1e-12,
            },
        )
        .expect("interface");
        let unloaded = interface
            .start(0.5 * (START_K[i] + START_K[i + 1]))
            .expect("unloaded");
        builder
            .add(&format!("interface{}", LINKS[i]), interface, unloaded)
            .expect("interface")
    });
    let planes = if one_plane {
        one_plane_for_both(&mut builder, cells, interfaces, halves)
    } else {
        a_plane_each(&mut builder, cells, interfaces, halves)
    };
    builder
        .heat_route(
            port(&builder, cells, "environment"),
            &[],
            Some(heat_sink("room", 300.0).expect("sink")),
            RouteTemperature::Fixed(300.0),
        )
        .expect("environment");
    let mut system = builder.seal().expect("sealed");
    // The temperatures in the loops need a start: each plane's body ports at
    // their cells, each interface between its two cells.
    let mut guess = |node: NodeId, name: &str, kelvin: f64| {
        let at = system.port(node, name).expect("port");
        assert!(system.guess(at, &[kelvin]));
    };
    if one_plane {
        for (i, kelvin) in START_K.into_iter().enumerate() {
            guess(planes[0], &format!("body{i}"), kelvin);
        }
    } else {
        for i in 0..2 {
            guess(planes[i], "body0", START_K[i]);
            guess(planes[i], "body1", START_K[i + 1]);
        }
    }
    for i in 0..2 {
        guess(interfaces[i], "heat", 0.5 * (START_K[i] + START_K[i + 1]));
    }
    Network {
        system,
        cells,
        planes,
        interfaces,
        routes,
    }
}

impl Network {
    /// The timeline's state.
    fn thermal(&self) -> &DepositionSimulation<'static> {
        self.system
            .state::<TimelineState<'static>>(self.cells)
            .expect("cells")
            .simulation()
    }

    /// The temperature each interface was given at the last step [K].
    fn given_k(&self) -> [f64; 2] {
        self.interfaces.map(|node| {
            self.system
                .last_input(self.system.port(node, "heat").expect("port"))[0]
        })
    }

    /// The temperature each interface's plane returned at the last step [K].
    fn plane_k(&self) -> [f64; 2] {
        self.planes.map(|node| {
            self.system
                .last_output(self.system.port(node, "source").expect("port"))[0]
        })
    }

    /// Each interface's dissipated work since the start, as the ledger
    /// booked it from the interface's receipts [J].
    fn dissipated_j(&self) -> [f64; 2] {
        let ledger = self.system.ledger();
        LINKS.map(|id| {
            let name = format!("interface{id}");
            ledger
                .nodes
                .iter()
                .find(|node| node.name == name)
                .expect("node")
                .dissipated_j
                .value()
        })
    }

    /// One kernel step of `DT`, then the module's verdict on its endpoint.
    fn step(&mut self, cx: &Cx<'_>) -> Result<Judged, DepositionError> {
        let old = self.thermal().clone();
        let before = self.dissipated_j();
        let taken = self.system.step(DT, cx).expect("kernel step");
        assert_eq!(taken.dt_s.to_bits(), DT.to_bits());
        let after = self.dissipated_j();
        let judged = self.judge(&old, [after[0] - before[0], after[1] - before[1]])?;
        Ok(Judged {
            iterations: taken.iterations,
            ..judged
        })
    }

    /// Hand the graph's present state, reached from `old` in one step, to
    /// the module: its own step on `old` with both routes, evaluated at the
    /// graph's energy increments, at the temperature each interface was
    /// given and at the work each dissipated, and accepted at the tests'
    /// tolerances.
    fn judge(
        &self,
        old: &DepositionSimulation<'static>,
        dissipated_j: [f64; 2],
    ) -> Result<Judged, DepositionError> {
        let new = self.thermal();
        let module = FrictionHeatNetworkStep::new(
            old,
            old.time_s() + DT,
            &[&self.routes[0], &self.routes[1]],
            2,
            || Ok(()),
        )
        .expect("module step");
        assert_eq!(module.link_ids(), &LINKS);
        let energy: Vec<f64> = module
            .cell_ids()
            .iter()
            .map(|id| {
                let (before, after) = (&old.cells()[id], &new.cells()[id]);
                before.mass_kg * (after.specific_h_j_kg - before.specific_h_j_kg)
            })
            .collect();
        let given_k = self.given_k();
        let trial = module.evaluate(
            &energy,
            &given_k,
            |_, _| {
                Ok::<_, DepositionError>(FrictionHeatNetworkInput {
                    interfaces: (0..2)
                        .map(|i| InterfaceHeatInput {
                            dissipated_j: dissipated_j[i],
                            conductivity_w_m_k: [K[i], K[i + 1]],
                        })
                        .collect(),
                    transport: StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    },
                })
            },
            || Ok(()),
        )?;
        let residual_j = trial
            .residual_j()
            .iter()
            .fold(0.0_f64, |worst, r| worst.max(r.abs()));
        let refused = {
            let report = trial.report();
            let [r10, r20] = [0, 1].map(|i| report.interfaces[i].interface_residual_j);
            format!(
                "interface residuals {r10:e} J (link 10), {r20:e} J (link 20); total {:e} J",
                report.total_energy_residual_j
            )
        };
        let (next, (), report) = trial
            .accept(tolerances(), || Ok(()))
            .inspect_err(|error| println!("the module refuses: {error:?}; {refused}"))?;
        // The evidence is the acceptance above. The cells the module
        // publishes are the increments handed to it added back, so they
        // return the graph's to rounding by construction; this only checks
        // the hand-over.
        assert_eq!(next.time_s().to_bits(), new.time_s().to_bits());
        for (id, cell) in next.cells() {
            let held = new.cells()[id].specific_h_j_kg;
            let scale = held.abs() + old.cells()[id].specific_h_j_kg.abs();
            assert!(
                (cell.specific_h_j_kg - held).abs() <= 4.0 * f64::EPSILON * scale,
                "cell {}: {} vs {held}",
                id.0,
                cell.specific_h_j_kg
            );
        }
        Ok(Judged {
            report,
            residual_j,
            dissipated_j,
            given_k,
            iterations: 0,
        })
    }
}

/// The independent backward Euler of the three cells joined by the two
/// contacts' series conductances, with no friction: one step of `DT` from
/// `kelvin`, by elimination of the tridiagonal system written here.
fn tridiagonal_step(kelvin: [f64; 3]) -> [f64; 3] {
    let g = [0, 1].map(|i| DT * series_conductance_w_k(AREA, [A; 2], [K[i], K[i + 1]]).unwrap());
    let left = CAP[0] + g[0];
    let right = CAP[2] + g[1];
    let middle = CAP[1] + g[0] + g[1] - g[0] * g[0] / left - g[1] * g[1] / right;
    let t1 =
        (CAP[1] * kelvin[1] + g[0] * CAP[0] * kelvin[0] / left + g[1] * CAP[2] * kelvin[2] / right)
            / middle;
    [
        (CAP[0] * kelvin[0] + g[0] * t1) / left,
        t1,
        (CAP[2] * kelvin[2] + g[1] * t1) / right,
    ]
}

/// The cells' temperatures [K].
fn kelvin(thermal: &DepositionSimulation<'static>) -> [f64; 3] {
    [0, 1, 2].map(|i| thermal.temperature_k(CellId(i)).expect("T"))
}

#[test]
fn two_contacts_that_share_a_cell_are_one_step_of_the_network() {
    // Interface 10 with the third case's falling coefficient at 0.1 m/s,
    // interface 20 at 0.3 and 0.2 m/s, both under 10 N, 40 steps of 10 ms.
    let gate = CancelGate::new_clock_free();
    let mut network = declare([mu_falling, mu_three_tenths], [0.1, 0.2], false);
    let start = kelvin(network.thermal());
    for step in 0..40 {
        let judged = with_cx(&gate, |cx| network.step(cx)).unwrap_or_else(|error| {
            panic!("step {step}: the module refuses the graph's step: {error:?}")
        });
        let report = &judged.report;
        let [r10, r20] = [0, 1].map(|i| report.interfaces[i].interface_residual_j);
        let open_k = [0, 1].map(|i| network.plane_k()[i] - judged.given_k[i]);
        println!(
            "network step {step}: accepted; largest residual {:.3e} J; interface residuals {r10:.3e}, \
             {r20:.3e} J; total {:.3e} J; dissipated {:?} J; interfaces at {:?} K; planes' wires open by \
             {open_k:?} K; {} kernel iterations",
            judged.residual_j,
            report.total_energy_residual_j,
            judged.dissipated_j,
            judged.given_k,
            judged.iterations
        );
        // Both contacts rub, and both planes heat the shared cell in this
        // one step.
        assert!(judged.dissipated_j.iter().all(|&d| d > 0.0));
        assert!(report.interfaces[0].body_heat_j[1] != 0.0);
        assert!(report.interfaces[1].body_heat_j[0] != 0.0);
        assert_eq!(report.interfaces[0].cells[1], CellId(1));
        assert_eq!(report.interfaces[1].cells[0], CellId(1));
    }
    let end = kelvin(network.thermal());
    println!("network: cells from {start:?} K to {end:?} K after 40 steps");
    // The friction heat stays: the network's energy rose by what was
    // dissipated, which the module checked step by step.
    let gained_j: f64 = (0..3).map(|i| CAP[i] * (end[i] - start[i])).sum();
    println!("network: the cells gained {gained_j:.6e} J");
    assert!(gained_j > 0.0);
}

#[test]
fn without_friction_the_three_cells_follow_the_tridiagonal_backward_euler() {
    // No friction: the planes only conduct, two massless nodes in series
    // between the cells, and the cells are the independent tridiagonal
    // backward Euler above, step by step from the graph's own cells. Each
    // cell is held to it within the module tests' cell tolerance (1e-10 J)
    // in energy.
    let gate = CancelGate::new_clock_free();
    let mut network = declare([frictionless, frictionless], [0.1, 0.2], false);
    let mut worst_j = 0.0_f64;
    for step in 0..20 {
        let before = kelvin(network.thermal());
        let judged = with_cx(&gate, |cx| network.step(cx)).unwrap_or_else(|error| {
            panic!("step {step}: the module refuses the graph's step: {error:?}")
        });
        let after = kelvin(network.thermal());
        let expected = tridiagonal_step(before);
        let off_j: [f64; 3] = [0, 1, 2].map(|i| CAP[i] * (after[i] - expected[i]));
        worst_j = off_j
            .iter()
            .fold(worst_j, |worst, off| worst.max(off.abs()));
        println!(
            "frictionless step {step}: cells {after:?} K; off the tridiagonal step by {off_j:?} J; \
             dissipated {:?} J; module residual {:.3e} J; {} kernel iterations",
            judged.dissipated_j, judged.residual_j, judged.iterations
        );
        assert!(judged.dissipated_j.iter().all(|&d| d == 0.0));
        for off in off_j {
            assert!(off.abs() <= tolerances().cell_j, "{off}");
        }
    }
    println!("frictionless: worst cell off the tridiagonal step {worst_j:.3e} J over 20 steps");
    // Heat ran from the hot end to the cold one through the shared cell.
    let end = kelvin(network.thermal());
    assert!(end[0] < START_K[0] && end[2] > START_K[2]);
}

#[test]
fn both_contacts_through_one_plane_are_not_the_networks_rule() {
    // The same two contacts heating ONE plane joined to all three cells. The
    // module gives each contact its own plane and judges each interface's
    // balance at its own temperature through its own route; here both are
    // given the one plane's temperature.
    let gate = CancelGate::new_clock_free();
    let mut network = declare([mu_falling, mu_three_tenths], [0.1, 0.2], true);
    let verdict = with_cx(&gate, |cx| network.step(cx));
    match &verdict {
        Err(error) => println!("one plane: the module refuses the first step: {error:?}"),
        Ok(judged) => println!(
            "one plane: the module ACCEPTS the first step: largest residual {:.3e} J",
            judged.residual_j
        ),
    }
    assert!(
        matches!(verdict, Err(DepositionError::EnergyBalance(_))),
        "one plane for both contacts must be refused"
    );
}
