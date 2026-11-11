//! The friction-heat coupling's six base cases, declared as a graph.
//!
//! `fs_couple::friction_heat` solves the temperature of a massless plane
//! between two cells inside one backward-Euler thermal step, with the
//! friction work generated on the plane. Here the same device is declared
//! from parts, and nothing is written for the pair:
//!
//! - the thermal timeline of the module's own test network, stepped by its
//!   owner's backward-Euler solve, with the contact link named as carried
//!   outside and the two contact cells' ports returning end-of-step
//!   temperatures (`fs_thermal::component::ThermalTimeline`);
//! - the kernel's `HeatNode`, the massless plane, joined to the two cells
//!   by the route's two half conductances;
//! - fs-tribo's `SlidingInterface` under the tests' prescribed slip, its
//!   heat port on the plane.
//!
//! The module is the reference and stays. The criterion is the module's own
//! verdict: the endpoint of every kernel step (the cells' energy increments
//! and the plane's temperature) is handed to `FrictionHeatStep::evaluate`
//! with the module tests' own law callback, and must be accepted at the
//! tests' tolerances of 1e-10 J; what the module then publishes (cells,
//! tangential state, traction) must be what the graph holds. The tests' own
//! independent references (their LU endpoint in temperature coordinates,
//! their closed forms) are checked on the graph as well.
//!
//! The graph has no equivalent of: conductances that follow the materials'
//! conductivities at the endpoint (the plane's are constants, as the
//! conductivities of the six cases are), a contact resistance split about
//! the plane, the module's entropy audit, or its area, moving-patch,
//! third-body, rigid and other variants.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::elements::HeatNode;
use fs_system::{
    NodeId, RouteTemperature, SolveConfig, System, SystemBuilder, SystemError, heat_sink,
};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionError, ImplicitConductionGraph,
    ImplicitConductionLaw, ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_tribo::tangential::component::{PrescribedSlip, SlidingInterface, SlidingState};
use fs_tribo::{DryTangentialState, TangentialStep, TriboError};

use super::*;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// The fixture's caloric law with one side's constant conductivity, as the
/// implicit conduction graph binds a material.
struct Conducting {
    w_m_k: f64,
}

impl ConductivityLaw for Conducting {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.w_m_k)
    }
}

impl ImplicitConductionLaw for Conducting {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        &SENSIBLE
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    /// 50 K to 5300 K on the fixture's `T = 300 K + h`.
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [-250.0, 5000.0]
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: SENSIBLE.temperature_k(specific_h_j_kg)?,
            temperature_slope_kg_k_j: 1.0,
            conductivity_w_m_k: self.w_m_k,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// The module tests' network with each side in its own material, so that
/// the conduction graph can give the sides their own conductivities. The
/// cells, links, films and numbers are those of `network`.
fn sided(
    n: usize,
    lengths: [f64; 2],
    k: [f64; 2],
    capacity_per_volume: [f64; 2],
    cooling: [f64; 2],
    initial_rise: [f64; 2],
) -> DepositionSimulation<'static> {
    network_of(
        n,
        lengths,
        k,
        capacity_per_volume,
        cooling,
        initial_rise,
        &[(MaterialId(1), &SENSIBLE), (MaterialId(2), &SENSIBLE)],
    )
}

/// The friction law at an interface temperature.
type LawAt = fn(f64) -> Result<FrictionLaw, TriboError>;

fn coulomb(mu: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Coulomb {
        static_mu: mu,
        kinetic_mu: mu,
    })
}

fn mu_one(_: f64) -> Result<FrictionLaw, TriboError> {
    coulomb(1.0)
}

/// The third case's declared law: `0.8 - 0.01 (T - 300 K)`.
fn mu_falling(kelvin: f64) -> Result<FrictionLaw, TriboError> {
    coulomb(0.8 - 0.01 * (kelvin - 300.0))
}

fn mu_three_tenths(_: f64) -> Result<FrictionLaw, TriboError> {
    coulomb(0.3)
}

const LAW_REFUSAL: TriboError = TriboError::InvalidInput {
    field: "the test law refuses",
};

/// A law that is there at 300 K, where the interface starts, and refused a
/// quarter of a kelvin above, where the first step takes the plane.
fn refusing(kelvin: f64) -> Result<FrictionLaw, TriboError> {
    if kelvin > 300.25 {
        Err(LAW_REFUSAL)
    } else {
        coulomb(1.0)
    }
}

/// The declared graph and what the module needs to judge its steps.
struct Graph {
    system: System,
    cells: NodeId,
    plane: NodeId,
    interface: NodeId,
    route: FrictionHeatRoute,
    k: [f64; 2],
    law_at: LawAt,
    slip_rate_m_s: f64,
}

/// What the module made of one kernel step.
struct Verdict {
    report: FrictionHeatReport,
    friction: TangentialStep<DryTangentialState>,
    /// Largest of the module's residuals at the graph's endpoint [J].
    residual_j: f64,
    /// What the kernel's loop left open on the wire from the plane to the
    /// interface: the temperature the plane returned less the one the
    /// interface was last given [K].
    wire_k: f64,
}

/// The timeline of `thermal`, a network of `sided`, on the contract: every
/// link but the contact conducted by the owner's implicit solve, the
/// contact named as carried outside, a heat port on each of `contact`'s two
/// cells, thermal intervals up to `dt`. With `end`, the cell ports return
/// end-of-step temperatures; without, the default means of the two ends.
fn timeline(
    thermal: &DepositionSimulation<'static>,
    lengths: [f64; 2],
    k: [f64; 2],
    contact: [CellId; 2],
    dt: f64,
    end: bool,
) -> ThermalTimeline<'static> {
    let n = thermal.cells().len() / 2;
    let laws: [&'static Conducting; 2] = k.map(|w_m_k| &*Box::leak(Box::new(Conducting { w_m_k })));
    // Every link but the contact is the timeline's own, with the half
    // distances its authored conductance was written from.
    let inside: Vec<ConductiveLink> = thermal
        .links()
        .values()
        .filter(|link| link.id != CONTACT)
        .map(|link| {
            let side = usize::from(link.a.0 as usize > n);
            ConductiveLink {
                link_id: link.id,
                half_distances_m: [0.5 * (lengths[side] / n as f64); 2],
                contact_resistance_m2_k_w: 0.0,
            }
        })
        .collect();
    let conduction = ConductionObserver::new(
        thermal,
        &[
            (MaterialId(1), laws[0] as &dyn ConductivityLaw),
            (MaterialId(2), laws[1] as &dyn ConductivityLaw),
        ],
        &inside,
        &[],
        || Ok(()),
    )
    .expect("conduction")
    .with_links_carried_outside(thermal, &[CONTACT], || Ok(()))
    .expect("the contact carried outside");
    let solver = ImplicitConductionGraph::new(
        thermal,
        &conduction,
        &[
            (MaterialId(1), laws[0] as &dyn ImplicitConductionLaw),
            (MaterialId(2), laws[1] as &dyn ImplicitConductionLaw),
        ],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .expect("solver");
    let controls = ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-13,
        total_tolerance_j: 1e-12,
        ..ImplicitConductionControls::default()
    };
    let timeline =
        ThermalTimeline::implicit(contact.to_vec(), solver, controls, ImplicitScheme::BackwardEuler)
            .expect("timeline");
    if end {
        timeline.returning_end_temperatures().expect("end temperatures")
    } else {
        timeline
    }
}

/// The graph on `thermal`, a network of `sided`: conductivities `k`, the
/// law `law_at` sliding at `slip_rate_m_s` under the tests' 10 N, steps up
/// to `dt`, and the cell ports of [`timeline`] as `end` says.
fn declare(
    thermal: &DepositionSimulation<'static>,
    lengths: [f64; 2],
    k: [f64; 2],
    law_at: LawAt,
    slip_rate_m_s: f64,
    dt: f64,
    end: bool,
) -> Graph {
    let route = FrictionHeatRoute::from_link(thermal, CONTACT, k).expect("route");
    let [a, b] = route.cells();
    let timeline = timeline(thermal, lengths, k, [a, b], dt, end);
    let start = timeline.start(thermal).expect("start");
    // The plane's two conductances in the route's own arithmetic, A (k / d).
    let (area, distances) = (route.contact_area_m2(), route.distances_m());
    let plane = HeatNode::new(&[
        area * (k[0] / distances[0]),
        area * (k[1] / distances[1]),
    ])
    .expect("plane");
    let interface = SlidingInterface::new(
        law_at,
        PrescribedSlip {
            normal_n: input(dt).normal_n,
            slip_rate_m_s: [slip_rate_m_s, 0.0],
            nominal_pressure_pa: input(dt).nominal_pressure_pa,
            solve: input(dt).solve,
            work_tolerance_j: 1e-12,
        },
    )
    .expect("interface");
    let unloaded = interface.start(300.0).expect("unloaded");

    let mut builder = SystemBuilder::new();
    let cells = builder.add("cells", timeline, start).expect("cells");
    let plane = builder.add_with_tangent("plane", plane, ()).expect("plane");
    let interface = builder.add("interface", interface, unloaded).expect("interface");
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).expect("port");
    for (cell, body) in [(a, "body0"), (b, "body1")] {
        let side = port(&builder, cells, &format!("cell{}", cell.0));
        builder
            .common_effort(side, &[port(&builder, plane, body)])
            .expect("a cell on the plane");
    }
    builder
        .common_effort(port(&builder, plane, "source"), &[port(&builder, interface, "heat")])
        .expect("the interface on the plane");
    builder
        .heat_route(
            port(&builder, cells, "environment"),
            &[],
            Some(heat_sink("room", 300.0).expect("sink")),
            RouteTemperature::Fixed(300.0),
        )
        .expect("environment");
    let mut system = builder.seal().expect("sealed");
    // The module's tests accept 1e-10 J. The kernel's default closes a loop
    // to 1e-12 of its largest unknown of a kind, which for temperatures
    // near 300 K is 3e-10 K, and through these conductances over a step
    // that is as many joules: at the default the module refused the third
    // case's ninth step with 2.4e-10 J (observed 2026-11-11). The loop is
    // closed two orders tighter here.
    system.set_config(SolveConfig {
        relative_tolerance: 1.0e-14,
        ..SolveConfig::default()
    });
    // The temperatures in the loop need a start.
    let kelvin = [a, b].map(|cell| thermal.temperature_k(cell).expect("T"));
    for (body, start_k) in ["body0", "body1"].into_iter().zip(kelvin) {
        let at = system.port(plane, body).expect("port");
        assert!(system.guess(at, &[start_k]));
    }
    let at = system.port(interface, "heat").expect("port");
    assert!(system.guess(at, &[0.5 * (kelvin[0] + kelvin[1])]));
    Graph {
        system,
        cells,
        plane,
        interface,
        route,
        k,
        law_at,
        slip_rate_m_s,
    }
}

impl Graph {
    /// The timeline's state.
    fn thermal(&self) -> &DepositionSimulation<'static> {
        self.system
            .state::<TimelineState<'static>>(self.cells)
            .expect("cells")
            .simulation()
    }

    /// The interface's state.
    fn friction(&self) -> &SlidingState<DryTangentialState> {
        self.system.state(self.interface).expect("interface")
    }

    /// The temperature the plane returned at the last step [K].
    fn plane_k(&self) -> f64 {
        self.system.last_output(self.system.port(self.plane, "source").expect("port"))[0]
    }

    /// The temperature the interface was given at the last step [K]: the
    /// one its law was evaluated at and its state holds. The kernel solves
    /// the loop these two are on by iteration, so it equals the plane's
    /// only to what that solve leaves (an ulp or so here).
    fn interface_k(&self) -> f64 {
        self.system.last_input(self.system.port(self.interface, "heat").expect("port"))[0]
    }

    /// The traction at the last step [N].
    fn traction_n(&self) -> [f64; 2] {
        let wire = self.system.last_output(self.system.port(self.interface, "traction").expect("port"));
        [wire[0], wire[1]]
    }

    /// Heat that entered each contact cell from the plane in the last step
    /// of `dt` [J], read off the two wires.
    fn body_heat_j(&self, dt: f64) -> [f64; 2] {
        self.route.cells().map(|cell| {
            let at = self.system.port(self.cells, &format!("cell{}", cell.0)).expect("port");
            dt * self.system.last_input(at)[0]
        })
    }

    /// The slip the tests' input prescribes over `dt`, at this graph's rate.
    fn slip(&self, dt: f64) -> TangentialStepInput {
        TangentialStepInput {
            slip_increment_m: [self.slip_rate_m_s * dt, 0.0],
            ..input(dt)
        }
    }

    /// One kernel step of `dt`, then the module's verdict on its endpoint.
    fn step(&mut self, dt: f64, cx: &Cx<'_>) -> Result<Verdict, DepositionError> {
        let old = self.thermal().clone();
        let history = self.friction().tangential().clone();
        let taken = self.system.step(dt, cx).expect("kernel step");
        assert_eq!(taken.dt_s, dt);
        self.verdict(&old, &history, dt)
    }

    /// Hand the graph's present state, reached from `old` and `history` in
    /// one step of `dt`, to the module: its own step on `old`, evaluated at
    /// the graph's energy increments and at the temperature the graph's
    /// interface was given, with the module tests' law callback, and
    /// accepted at the tests' tolerances.
    fn verdict(
        &self,
        old: &DepositionSimulation<'static>,
        history: &DryTangentialState,
        dt: f64,
    ) -> Result<Verdict, DepositionError> {
        let new = self.thermal();
        let module =
            FrictionHeatStep::new(old, old.time_s() + dt, &self.route, || Ok(())).expect("module step");
        let energy: Vec<f64> = module
            .cell_ids()
            .iter()
            .map(|id| {
                let (before, after) = (&old.cells()[id], &new.cells()[id]);
                before.mass_kg * (after.specific_h_j_kg - before.specific_h_j_kg)
            })
            .collect();
        let (law_at, k, slip) = (self.law_at, self.k, self.slip(dt));
        let trial = module.evaluate(
            &energy,
            self.interface_k(),
            |_, temperature| {
                let law = law_at(temperature).expect("law");
                let friction =
                    tangential_heat_candidate(&law, history, slip, temperature, 1e-12).expect("friction");
                Ok::<_, DepositionError>(source(friction.dissipated_j, k, friction))
            },
            || Ok(()),
        )?;
        let residual_j = trial.residual_j().iter().fold(0.0_f64, |worst, r| worst.max(r.abs()));
        let (next, friction, report) = trial.accept(tolerances(), || Ok(()))?;
        // What the module publishes is what the graph holds: the cells to
        // the rounding of forming an increment and adding it back, the
        // tangential state and the traction exactly.
        assert_eq!(next.time_s(), new.time_s());
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
        assert_eq!(&friction.state, self.friction().tangential());
        assert_eq!(friction.traction_n, self.traction_n());
        assert_eq!(report.interface_temperature_k, self.interface_k());
        // The loop the plane and the interface are on is closed to the
        // tolerance set in `declare` and the rounding the kernel allows a
        // temperature (16 units).
        let wire_k = self.plane_k() - self.interface_k();
        assert!(wire_k.abs() <= (1e-14 + 32.0 * f64::EPSILON) * self.plane_k(), "{wire_k}");
        Ok(Verdict {
            report,
            friction,
            residual_j,
            wire_k,
        })
    }
}

#[test]
fn friction_work_enters_the_network_once_through_the_plane() {
    // The module's first case: one cell a side, half cells of 0.5 m and
    // 1.5 m, conductivities 2 and 9, 10 W of friction for 0.2 s.
    let gate = CancelGate::new_clock_free();
    let (lengths, k, dt) = ([1.0, 3.0], [2.0, 9.0], 0.2);
    let old = sided(1, lengths, k, [2.0, 2.0], [0.0; 2], [0.0; 2]);
    let mut graph = declare(&old, lengths, k, mu_one, 1.0, dt, true);
    assert_eq!(graph.route.distances_m(), [0.5, 1.5]);
    let verdict = with_cx(&gate, |cx| graph.step(dt, cx)).expect("the module accepts the step");
    eprintln!(
        "first case: the module's largest residual at the graph's endpoint {:.3e} J (it accepts \
         1e-10 J); the plane's wire open by {:.3e} K",
        verdict.residual_j, verdict.wire_k
    );

    // The tests' own reference, an LU solve in temperature coordinates.
    let mut lu = factor(&old, dt, [4.0, 6.0], 0.0);
    let (energy, t) = endpoint(&old, &mut lu, dt, 10.0);
    close(graph.plane_k(), t, 1e-12);
    for (cell, by_lu) in graph.thermal().cells().values().zip(&energy) {
        close(cell.mass_kg * cell.specific_h_j_kg, *by_lu, 1e-12);
    }
    // Effective endpoint admittances include each finite heat capacity.
    let a = 4.0 * 2.0 / (2.0 + dt * 4.0);
    let b = 6.0 * 6.0 / (6.0 + dt * 6.0);
    let on_the_wires = graph.body_heat_j(dt);
    close(verdict.report.body_heat_j[0], 2.0 * a / (a + b), 1e-12);
    close(on_the_wires[0], 2.0 * a / (a + b), 1e-12);
    close(on_the_wires.iter().sum(), verdict.friction.dissipated_j, 1e-12);
    close(verdict.friction.dissipated_j, 2.0, 1e-12);
    close(graph.thermal().energy_j() - old.energy_j(), 2.0, 1e-12);
    assert_eq!(graph.thermal().contact_patch(CONTACT), old.contact_patch(CONTACT));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(graph.friction().tangential().contact_temperature_k(), Some(graph.interface_k()));
    assert_eq!(graph.thermal().cells()[&CellId(1)].boundary.source_w, 0.0);

    // The kernel's account of the same step: 2 J dissipated by the
    // interface, all of it stored in the cells, nothing in the room.
    let ledger = graph.system.ledger();
    let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).expect("node");
    close(account("interface").dissipated_j.value(), 2.0, 1e-12);
    close(account("cells").stored_change_j.value(), 2.0, 1e-12);
    close(account("interface").stored_change_j.value(), -2.0, 1e-12);
    assert_eq!(ledger.sink_heat_j(), 0.0);
    assert!(ledger.energy_residual_j().abs() <= 1e-12, "{}", ledger.energy_residual_j());
    // The work became heat at the plane's temperature, and the plane then
    // conducted it down to the cells.
    close(account("interface").entropy_generated_j_per_k.value(), 2.0 / graph.plane_k(), 1e-15);
    assert!(account("plane").entropy_generated_j_per_k.value() > 0.0);
}

#[test]
fn without_slip_the_plane_only_conducts() {
    // The module's second case: no slip, the first cell 10 K warmer.
    let gate = CancelGate::new_clock_free();
    let (lengths, k, dt) = ([1.0, 3.0], [2.0, 9.0], 0.2);
    let old = sided(1, lengths, k, [2.0, 2.0], [0.0; 2], [10.0, 0.0]);
    let mut graph = declare(&old, lengths, k, mu_one, 0.0, dt, true);
    let verdict = with_cx(&gate, |cx| graph.step(dt, cx)).expect("the module accepts the step");
    let exchange = 2.4 * 10.0 * dt / (1.0 + dt * 2.4 * (0.5 + 1.0 / 6.0));
    let on_the_wires = graph.body_heat_j(dt);
    for (side, sign) in [(0, -1.0), (1, 1.0)] {
        close(verdict.report.body_heat_j[side], sign * exchange, 1e-12);
        close(on_the_wires[side], sign * exchange, 1e-12);
    }
    close(graph.thermal().energy_j(), old.energy_j(), 1e-12);
    assert_eq!(verdict.report.dissipated_j, 0.0);
    assert_eq!(graph.traction_n(), [0.0; 2]);
    let ledger = graph.system.ledger();
    let interface = ledger.nodes.iter().find(|node| node.name == "interface").expect("node");
    assert_eq!(interface.dissipated_j.value(), 0.0);
    eprintln!(
        "second case: the module's largest residual at the graph's endpoint {:.3e} J",
        verdict.residual_j
    );
}

#[test]
fn friction_that_falls_with_temperature_settles_where_the_module_says() {
    // The module's third case: mu = 0.8 - 0.01 (T - 300 K), films of 3 and
    // 5 W/K to a 300 K room, 80 steps of 0.25 s. Beside the graph, the
    // module's own run from its tests' LU endpoints.
    let gate = CancelGate::new_clock_free();
    let (lengths, k, dt) = ([1.0; 2], [1.0, 2.0], 0.25);
    let start = sided(1, lengths, k, [1.0; 2], [3.0, 5.0], [0.0; 2]);
    let mut graph = declare(&start, lengths, k, mu_falling, 1.0, dt, true);
    let route = FrictionHeatRoute::from_link(&start, CONTACT, k).expect("route");
    let mut lu = factor(&start, dt, [2.0, 4.0], 0.1);
    let mut module = start.clone();
    let mut history = mu_falling(300.0).expect("law").initial_state();
    let (mut worst_residual_j, mut worst_apart_k, mut to_room_j) = (0.0_f64, 0.0_f64, 0.0_f64);
    let (mut last_force, mut worst_wire_k) = (0.0, 0.0_f64);
    with_cx(&gate, |cx| {
        for step in 0..80 {
            let verdict = graph
                .step(dt, cx)
                .unwrap_or_else(|refusal| panic!("step {step}: the module refuses: {refusal:?}"));
            worst_residual_j = worst_residual_j.max(verdict.residual_j);
            worst_wire_k = worst_wire_k.max(verdict.wire_k.abs());
            to_room_j -= verdict.report.external_heat_j;
            last_force = verdict.friction.traction_n[0];

            // The module's own step, as its test takes it.
            let (energy, t) = endpoint(&module, &mut lu, dt, 8.0);
            let own = FrictionHeatStep::new(&module, module.time_s() + dt, &route, || Ok(()))
                .expect("module step");
            let trial = own
                .evaluate(
                    &energy,
                    t,
                    |_, temperature| {
                        let law = mu_falling(temperature).expect("law");
                        let friction =
                            tangential_heat_candidate(&law, &history, input(dt), temperature, 1e-12)
                                .expect("friction");
                        Ok::<_, DepositionError>(source(friction.dissipated_j, k, friction))
                    },
                    || Ok(()),
                )
                .expect("trial");
            let (next, friction, report) = trial.accept(tolerances(), || Ok(())).expect("accepted");
            worst_apart_k =
                worst_apart_k.max((graph.plane_k() - report.interface_temperature_k).abs());
            history = friction.state;
            module = next;
        }
    });
    let effective = 2.0 * 3.0 / (2.0 + 3.0) + 4.0 * 5.0 / (4.0 + 5.0);
    let expected = 300.0 + 8.0 / (effective + 0.1);
    eprintln!(
        "third case, 80 steps: the module's largest residual at a graph endpoint \
         {worst_residual_j:.3e} J; the plane's wire open by {worst_wire_k:.3e} K at worst; the \
         plane {worst_apart_k:.3e} K at worst from the module's own run; {:.3e} K from the \
         analytic steady state; {to_room_j:.6} J to the room",
        (graph.plane_k() - expected).abs()
    );
    close(graph.plane_k(), expected, 1e-9);
    assert!(last_force < 8.0);
    close(last_force, 8.0 - 0.1 * (expected - 300.0), 1e-9);
    close(graph.traction_n()[0], last_force, 0.0);
    assert_eq!(graph.friction().tangential().contact_temperature_k(), Some(graph.interface_k()));
    assert!(worst_apart_k <= 1e-10, "{worst_apart_k}");
    // The films' heat is in the room, as the module's reports have it.
    let ledger = graph.system.ledger();
    assert!(to_room_j > 1.0, "{to_room_j}");
    close(ledger.sink_heat_j(), to_room_j, 1e-9 * to_room_j);
    assert!(ledger.energy_residual_j().abs() <= 1e-9 * to_room_j, "{}", ledger.energy_residual_j());
}

#[test]
fn a_refusal_publishes_nothing_and_mean_temperature_ports_are_not_the_modules_rule() {
    // The module's fourth case, for the graph: a refusal inside a step
    // leaves every state as it was.
    let (lengths, k, dt) = ([1.0, 3.0], [2.0, 9.0], 0.2);
    let old = sided(1, lengths, k, [2.0, 2.0], [0.0; 2], [0.0; 2]);
    let mut graph = declare(&old, lengths, k, refusing, 1.0, dt, true);
    let unloaded = graph.friction().clone();
    let gate = CancelGate::new_clock_free();
    let refused = with_cx(&gate, |cx| graph.system.step(dt, cx));
    assert!(
        matches!(&refused, Err(SystemError::Component { node, .. }) if node == "interface"),
        "{refused:?}"
    );
    assert_eq!(graph.thermal().cells(), old.cells());
    assert_eq!((graph.system.time_s(), graph.system.steps()), (0.0, 0));
    assert_eq!(graph.friction(), &unloaded);

    let mut graph = declare(&old, lengths, k, mu_one, 1.0, dt, true);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let stopped = with_cx(&cancelled, |cx| graph.system.step(dt, cx));
    assert!(matches!(stopped, Err(SystemError::Cancelled)), "{stopped:?}");
    assert_eq!(graph.thermal().cells(), old.cells());
    assert_eq!((graph.system.time_s(), graph.system.steps()), (0.0, 0));

    // The same graph with the default ports, which return the mean of the
    // two ends, is another discrete rule: the kernel takes its step, and
    // the module refuses the endpoint.
    let mut mean = declare(&old, lengths, k, mu_one, 1.0, dt, false);
    let refusal = with_cx(&gate, |cx| mean.step(dt, cx)).err();
    eprintln!("mean-temperature ports: the module refuses the graph's endpoint with {refusal:?}");
    let Some(DepositionError::EnergyBalance(defect_j)) = refusal else {
        panic!("{refusal:?}");
    };
    assert!(defect_j.abs() > 1e-3, "{defect_j}");
}

#[test]
fn the_contact_has_one_owner_in_the_graph_too() {
    // The module's fifth case refuses a second owner of the interface
    // conduction. In the graph the conduction binding is where that is
    // decided: the contact cannot be both conducted and carried outside,
    // and it cannot be left unsaid.
    let (lengths, k) = ([1.0; 2], [1.0; 2]);
    let old = sided(1, lengths, k, [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, k).expect("route");
    let law: &'static Conducting = Box::leak(Box::new(Conducting { w_m_k: 1.0 }));
    let materials = [
        (MaterialId(1), law as &dyn ConductivityLaw),
        (MaterialId(2), law as &dyn ConductivityLaw),
    ];
    let conducted = [ConductiveLink {
        link_id: CONTACT,
        half_distances_m: route.distances_m(),
        contact_resistance_m2_k_w: 0.0,
    }];
    let both = ConductionObserver::new(&old, &materials, &conducted, &[], || Ok(()))
        .expect("conducted")
        .with_links_carried_outside(&old, &[CONTACT], || Ok(()));
    assert_eq!(
        both.err(),
        Some(DepositionError::InvalidInput(
            "a link carried outside is also bound to a conductance"
        ))
    );
    let unsaid = ConductionObserver::new(&old, &materials, &[], &[], || Ok(())).expect("unsaid");
    let solver = ImplicitConductionGraph::new(
        &old,
        &unsaid,
        &[
            (MaterialId(1), law as &dyn ImplicitConductionLaw),
            (MaterialId(2), law as &dyn ImplicitConductionLaw),
        ],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    );
    assert!(matches!(
        solver.err(),
        Some(ImplicitConductionError::Invalid(
            "thermal graph bindings changed or direct-chain topology is invalid"
        ))
    ));
    // The module's own refusals stand beside these: a route whose link is
    // not the series conductance of the conductivities given.
    assert!(FrictionHeatRoute::from_link(&old, CONTACT, [2.0; 2]).is_err());
}

#[test]
fn the_partition_and_the_plane_temperature_converge_on_three_grids() {
    // The module's sixth case: constant surface generation between two
    // initially isothermal half spaces, share_A = e_A / (e_A + e_B) and
    // dT = 2 P sqrt(t) / (sqrt(pi) (e_A + e_B)), with e = [1, 2] and
    // P = 3 W; 512 steps to t = 0.25 s on 16, 32 and 64 cells a side.
    let gate = CancelGate::new_clock_free();
    let (lengths, k) = ([4.0; 2], [1.0, 4.0]);
    let dt = 0.25 / 512.0;
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let start = sided(n, lengths, k, [1.0; 2], [0.0; 2], [0.0; 2]);
        let mut graph = declare(&start, lengths, k, mu_three_tenths, 1.0, dt, true);
        let (mut received, mut worst_residual_j, mut worst_wire_k) = ([0.0; 2], 0.0_f64, 0.0_f64);
        with_cx(&gate, |cx| {
            for step in 0..512 {
                let verdict = graph.step(dt, cx).unwrap_or_else(|refusal| {
                    panic!("n {n}, step {step}: the module refuses: {refusal:?}")
                });
                worst_residual_j = worst_residual_j.max(verdict.residual_j);
                worst_wire_k = worst_wire_k.max(verdict.wire_k.abs());
                let on_the_wires = graph.body_heat_j(dt);
                for side in 0..2 {
                    received[side] += on_the_wires[side];
                }
            }
        });
        close(received.iter().sum(), 0.75, 1e-10);
        close(graph.thermal().energy_j(), 0.75, 1e-10);
        let partition_error = (received[0] / 0.75 - 1.0 / 3.0).abs();
        let temperature_error =
            ((graph.plane_k() - 300.0) - 1.0 / core::f64::consts::PI.sqrt()).abs();
        eprintln!(
            "sixth case, n={n}: partition_error={partition_error:e} \
             temperature_error={temperature_error:e}; the module's largest residual at a graph \
             endpoint {worst_residual_j:.3e} J over 512 steps; the plane's wire open by \
             {worst_wire_k:.3e} K at worst"
        );
        errors.push((partition_error, temperature_error));
    }
    for pair in errors.windows(2) {
        assert!(pair[1].0 < pair[0].0);
        assert!(pair[1].1 < pair[0].1);
    }
    assert!(errors[2].0 / (1.0 / 3.0) < 0.01);
    assert!(errors[2].1 / (1.0 / core::f64::consts::PI.sqrt()) < 0.02);
}
