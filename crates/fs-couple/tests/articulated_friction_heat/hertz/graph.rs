//! The Hertz case of the friction-heat family, declared as a graph.
//!
//! `fs_couple::friction_heat::hertz::propose_rigid_hertz_heat` takes one step
//! of two rigid bodies rubbing on one Hertz patch: the joint slip and spin
//! friction of fs-tribo's `HertzWrenchContact`, the bodies' recoil, and the
//! heat into a finite cell network through a massless face, all at the
//! interface temperature its own nonlinear solve finds. Here the same device
//! is declared from parts, and nothing is written for the pair:
//!
//! - fs-mbd's `RigidBody`, twice, with the module tests' masses, inertias,
//!   poses and momenta;
//! - the lever: one `power_map_members` connection from the two bodies'
//!   `force` and `torque` ports onto the patch. Its rows are the contact
//!   point's slip, `v + w x r` at each body's contact arm, and the patch
//!   radius times the relative spin; the bodies receive the adjoint, the
//!   force at the arm and its moment;
//! - fs-tribo's `HertzPatchInterface` on the kernel's load line. Its law at
//!   an interface temperature is the tests' own card refreshed there
//!   (`SourceThermalCoulomb::hertz_at`), as the module's `refresh_source`
//!   asks it;
//! - the kernel's `HeatNode` as the face, with the route's two half
//!   conductances, and fs-thermal's `ThermalTimeline` of the tests' three
//!   cells, stepped by its owner's backward-Euler solve with end-of-step
//!   temperatures: link 20, the face, carried outside; link 10 its own.
//!
//! The module is the reference and stays. The criterion is what its own
//! tests assert of a base case, at their tolerance (`close`, 3e-9 relative):
//! the interface temperature, the heat, the three cells, the force or spin
//! moment, the bodies' momenta (the arms' moments included), the slip or
//! spin after the step and the first law. The graph is held to the module's
//! candidate on the same inputs and to the tests' independent scalar root.
//!
//! The Hertz case is the module's: the normal load (10 N) and the patch
//! radius are given, and "Hertz" is the patch's pressure profile. No normal
//! contact law is stepped, here or in the module.
//!
//! Limit: the lever is a constant map, fixed when the graph is admitted, so
//! its arms are those of the start pose. A lever whose arms turn with its
//! bodies keeps the ratios of that pose. The module forms its pair's
//! compliance and impulse at the start pose as well. In the spin case the
//! arms lie on the spin axis; in the sliding case the bodies turn about y by
//! about 4e-4 rad within the step.
//!
//! In the spin case body a spins at 100 rad/s, a radian within the step, and
//! a torque across the spin axis would precess it. The line the kernel
//! differences from the bodies carries that gyroscopic coupling as an
//! antisymmetric slip block, which the module's start-pose compliance leaves
//! out. `HertzPatchInterface` steps on the full line (its "The symmetric
//! part"). The patch's force is zero in a pure spin, so the block moves the
//! step only through the line's rounding and the kernel's differencing; each
//! call of the patch is recorded (`Recorded`), and the test prints its owner
//! steps, its work balance against the rounding allowance, and the columns
//! the kernel differences.
//!
//! The graph has no equivalent of: the module's entropy audit, its receipts
//! of the card's pinned claims (the law's coefficients are compared here,
//! not its receipts), its check of the face's geometry against the contact
//! point, or its trial count.

use std::cell::RefCell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, Cx, ExecMode, StreamKey};
use fs_mbd::component::RigidBody;
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_ports::component::{Class, Component, LoadLine, PortSpec, Proposal};
use fs_system::elements::HeatNode;
use fs_system::{NodeId, PortRef, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, ThermalLink, TimedEvent, series_conductance_w_k,
};
use fs_tribo::TriboError;
use fs_tribo::sliding_patch::component::{HertzPatchInterface, HertzPatchLoad, HertzPatchState};
use fs_tribo::sliding_patch::contact::HertzWrenchContact;

use super::super::{CALORIC, LINK};
use super::*;

/// The normal load `Pair::contact` gives the module [N].
const NORMAL_N: f64 = 10.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// The fixture's caloric law, `T = 300 K + h`, with one cell's constant
/// conductivity, as the implicit conduction graph binds a material.
struct Conducting {
    w_m_k: f64,
}

static CONDUCTING: [Conducting; 3] = [
    Conducting { w_m_k: K[0] },
    Conducting { w_m_k: K[1] },
    Conducting { w_m_k: K[2] },
];

impl ConductivityLaw for Conducting {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.w_m_k)
    }
}

impl ImplicitConductionLaw for Conducting {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        &CALORIC
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
            temperature_k: CALORIC.temperature_k(specific_h_j_kg)?,
            temperature_slope_kg_k_j: 1.0,
            conductivity_w_m_k: self.w_m_k,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// The tests' network, `thermal([300.0; 3], So3::identity())`, with each cell
/// in its own material, so that the conduction graph can give each cell its
/// conductivity. The cells, links and numbers are `thermal`'s.
fn own_materials() -> DepositionSimulation<'static> {
    let endpoints = [[-A, A], [-3.0 * A, -A], [A, 3.0 * A]];
    let mut events: Vec<_> = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: CAP[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, endpoints[i][0]],
                    end_m: [0.0, 0.0, endpoints[i][1]],
                    volume_m3: AREA * 2.0 * A,
                    surface_area_m2: 2.0 * AREA + 8.0 * AREA.sqrt() * A,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINK[i],
                a: CellId(0),
                b: CellId(i as u64 + 1),
                area_m2: AREA,
                conductance_w_k: series_conductance_w_k(AREA, [A; 2], [K[0], K[i + 1]]).unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &CALORIC as &dyn ThermalLaw),
            (MaterialId(1), &CALORIC as &dyn ThermalLaw),
            (MaterialId(2), &CALORIC as &dyn ThermalLaw),
        ],
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

/// The three cells on the contract: link 10 conducted by the owner's
/// backward-Euler solve with the half distances its authored conductance
/// was written from, link 20 (the face) carried outside, heat ports on the
/// face's two cells returning end-of-step temperatures.
fn cells(network: &DepositionSimulation<'static>) -> ThermalTimeline<'static> {
    let conductivity: Vec<(MaterialId, &'static dyn ConductivityLaw)> = (0..3)
        .map(|i| {
            (
                MaterialId(i as u64),
                &CONDUCTING[i] as &'static dyn ConductivityLaw,
            )
        })
        .collect();
    let inside = [ConductiveLink {
        link_id: LINK[1],
        half_distances_m: [A, A],
        contact_resistance_m2_k_w: 0.0,
    }];
    let conduction = ConductionObserver::new(network, &conductivity, &inside, &[], || Ok(()))
        .expect("conduction")
        .with_links_carried_outside(network, &[LINK[0]], || Ok(()))
        .expect("the face carried outside");
    let implicit: Vec<(MaterialId, &'static dyn ImplicitConductionLaw)> = (0..3)
        .map(|i| {
            (
                MaterialId(i as u64),
                &CONDUCTING[i] as &'static dyn ImplicitConductionLaw,
            )
        })
        .collect();
    let solver = ImplicitConductionGraph::new(
        network,
        &conduction,
        &implicit,
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
        vec![CellId(0), CellId(1)],
        solver,
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .expect("timeline")
    .returning_end_temperatures()
    .expect("end temperatures")
}

/// The patch's law at an interface temperature: the tests' card refreshed
/// there through its pinned claims, as the module's `refresh_source` asks
/// it (no nominal pressure, the tests' normal load, the face's radius).
fn on_card(owner: &SourceThermalCoulomb, kelvin: f64) -> Result<HertzWrenchContact, TriboError> {
    owner
        .hertz_at(kelvin, None, NORMAL_N, radius())
        .map(|source| *source.law())
        .map_err(|_| TriboError::InvalidInput {
            field: "the card's Hertz law",
        })
}

/// One call of the patch on the line that took its step: what the kernel
/// handed it and what it gave.
#[derive(Clone, Copy, Debug)]
struct Pass {
    inputs: [f64; 4],
    line: [f64; 9],
    outputs: [f64; 4],
    /// The owner's steps the call took (`HertzPatchState::owner_steps`).
    owner_steps: usize,
    /// The port's work less the heat, unsigned [J].
    residual_j: f64,
    /// The work's rounding allowance, `balance_rounding_j` [J].
    rounding_j: f64,
}

/// The patch, its calls on the line recorded. Every call is handed on
/// unchanged; the record is the test's, not the component's.
struct Recorded<F> {
    patch: HertzPatchInterface<F>,
    passes: Rc<RefCell<Vec<Pass>>>,
}

impl<F> Component for Recorded<F>
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    type State = HertzPatchState;
    type Error = TriboError;

    fn ports(&self) -> &[PortSpec] {
        self.patch.ports()
    }

    fn class(&self) -> Class {
        self.patch.class()
    }

    fn stored_energy_j(&self, state: &HertzPatchState) -> Result<f64, TriboError> {
        self.patch.stored_energy_j(state)
    }

    fn propose(
        &self,
        old: &HertzPatchState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut HertzPatchState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        self.patch.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

impl<F> LoadLine for Recorded<F>
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    fn line_ports(&self) -> &[usize] {
        self.patch.line_ports()
    }

    fn propose_on_line(
        &self,
        old: &HertzPatchState,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut HertzPatchState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        let proposal = self
            .patch
            .propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)?;
        if let Proposal::Step(receipt) = &proposal {
            self.passes.borrow_mut().push(Pass {
                inputs: inputs.try_into().expect("the patch's four inputs"),
                line: compliance.try_into().expect("the patch's 3 x 3 line"),
                outputs: (&*outputs).try_into().expect("the patch's four outputs"),
                owner_steps: new.owner_steps(),
                residual_j: receipt.residual,
                rounding_j: receipt.energy_tolerance_j,
            });
        }
        Ok(proposal)
    }
}

/// Prints every recorded pass of a case: its owner steps and its work
/// balance against the rounding allowance; then each column the kernel
/// differenced, from a pass on the same line whose inputs differ in one
/// coefficient by a differencing step; then the owner steps per pass and
/// the largest.
fn print_passes(case: &str, passes: &[Pass]) {
    for (j, pass) in passes.iter().enumerate() {
        println!(
            "Hertz {case} pass {j}: inputs {:?}, outputs {:?}, owner steps {}, work residual {:e} J, rounding {:e} J, {}",
            pass.inputs,
            pass.outputs,
            pass.owner_steps,
            pass.residual_j,
            pass.rounding_j,
            if pass.residual_j <= pass.rounding_j {
                "within"
            } else {
                "EXCEEDS"
            }
        );
    }
    for (j, pass) in passes.iter().enumerate() {
        let same_line = |base: &Pass| {
            base.line
                .iter()
                .zip(&pass.line)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        };
        let moved = |base: &Pass| {
            let differ: Vec<usize> = (0..4)
                .filter(|&k| base.inputs[k].to_bits() != pass.inputs[k].to_bits())
                .collect();
            match differ[..] {
                [k] if (pass.inputs[k] - base.inputs[k]).abs()
                    <= 1e-6 * base.inputs[k].abs().max(1.0) =>
                {
                    Some(k)
                }
                _ => None,
            }
        };
        let found = (0..j).rev().find_map(|i| {
            moved(&passes[i])
                .filter(|_| same_line(&passes[i]))
                .map(|k| (i, k))
        });
        if let Some((i, k)) = found {
            let base = &passes[i];
            let step = pass.inputs[k] - base.inputs[k];
            let column: [f64; 4] =
                core::array::from_fn(|r| (pass.outputs[r] - base.outputs[r]) / step);
            println!(
                "Hertz {case} column d/d input {k} at pass {j} from pass {i}, step {step:e}: {column:?}"
            );
        }
    }
    let steps: Vec<usize> = passes.iter().map(|pass| pass.owner_steps).collect();
    println!(
        "Hertz {case}: {} passes, owner steps per pass {steps:?}, largest {}",
        passes.len(),
        steps.iter().max().copied().unwrap_or(0)
    );
}

/// The declared graph.
struct Graph {
    system: System,
    bodies: [NodeId; 2],
    patch: NodeId,
    plane: NodeId,
    cells: NodeId,
    /// The cells before the step.
    start: DepositionSimulation<'static>,
    /// The patch's calls on the line.
    passes: Rc<RefCell<Vec<Pass>>>,
}

/// What a graph's step gave, in the module tests' terms.
struct Outcome {
    /// The temperature the patch was given [K].
    interface_k: f64,
    /// The temperature the face returned [K].
    plane_k: f64,
    /// The patch's heat over the step [J].
    heat_j: f64,
    /// The magnitude of the friction force [N], or of the spin moment
    /// [N m] in the spin case.
    effort: f64,
    bodies: [RigidBodyState; 2],
    cells_k: [f64; 3],
    /// The slip of the contact point along x after the step [m/s], or the
    /// relative spin [rad/s] in the spin case.
    after: f64,
    /// Kinetic energy and cell energy gained, less heat from outside [J].
    first_law_j: f64,
}

/// The lever's members: each body's `force` and `torque` port with its map
/// onto the patch's `[v_x, v_y, a w]`, row-major 3 x 3.
///
/// Body a (0) is on the slip's positive side and body b (1) on its
/// negative; the connection gives the patch `-sum P_k f_k`, so a's maps
/// carry -1 and b's +1. A force port's map takes the tangent components of
/// the centre's velocity. A torque port's takes `w x r` at the arm
/// `r = (0, 0, arm)` of the start pose, whose tangent components are
/// `(arm w_y, -arm w_x)`, and the patch radius times `w_z`.
fn lever(builder: &SystemBuilder, bodies: [NodeId; 2], arms: [f64; 2]) -> Vec<(PortRef, Vec<f64>)> {
    let r = radius();
    let mut members = Vec::new();
    for (i, sign) in [(0, -1.0), (1, 1.0)] {
        let arm = arms[i];
        let port = |name: &str| builder.port(bodies[i], name).expect("port");
        members.push((
            port("force"),
            vec![sign, 0.0, 0.0, 0.0, sign, 0.0, 0.0, 0.0, 0.0],
        ));
        members.push((
            port("torque"),
            vec![
                0.0,
                sign * arm,
                0.0,
                -sign * arm,
                0.0,
                0.0,
                0.0,
                0.0,
                sign * r,
            ],
        ));
    }
    members
}

/// The device of the module's base case (`pair`, sliding or spinning), the
/// patch's law `law_at`, and the lever's contact arms `arms`.
fn declare<F>(pair: &Pair, spin: bool, law_at: F, arms: [f64; 2]) -> Graph
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError> + 'static,
{
    let network = own_materials();
    let route = FrictionHeatRoute::from_link(&network, LINK[0], [K[0], K[1]]).expect("route");
    let timeline = cells(&network);
    let cells_start = timeline.start(&network).expect("cells start");
    // The face's two conductances in the route's own arithmetic, A (k / d).
    let (area, distances) = (route.contact_area_m2(), route.distances_m());
    let face =
        HeatNode::new(&[area * (K[0] / distances[0]), area * (K[1] / distances[1])]).expect("face");
    let patch = HertzPatchInterface::new(
        law_at,
        HertzPatchLoad {
            normal_n: NORMAL_N,
            radius_m: radius(),
            solve: controls().mechanical.contact,
        },
    )
    .expect("patch");
    let moving = if spin {
        [0.0, 0.0, radius() * SPIN]
    } else {
        [SPEED, 0.0, 0.0]
    };
    let rubbing = patch.start(moving).expect("patch start");
    let passes = Rc::new(RefCell::new(Vec::new()));
    let patch = Recorded {
        patch,
        passes: Rc::clone(&passes),
    };

    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).expect("port");
    let bodies = [0, 1].map(|i| {
        let body = RigidBody::new(RigidBodyIntegrator::new(Gravity::ZERO), pair.properties[i])
            .expect("body");
        builder
            .add(&format!("body{i}"), body, pair.bodies[i])
            .expect("body")
    });
    let patch = builder.add_on_line("patch", patch, rubbing).expect("patch");
    let plane = builder.add_with_tangent("plane", face, ()).expect("plane");
    let cells = builder.add("cells", timeline, cells_start).expect("cells");
    let members = lever(&builder, bodies, arms);
    builder
        .power_map_members(members, port(&builder, patch, "wrench"))
        .expect("the lever");
    builder
        .common_effort(
            port(&builder, plane, "source"),
            &[port(&builder, patch, "heat")],
        )
        .expect("the patch on the face");
    for (cell, side) in [("cell0", "body0"), ("cell1", "body1")] {
        builder
            .common_effort(port(&builder, cells, cell), &[port(&builder, plane, side)])
            .expect("a cell on the face");
    }
    builder
        .heat_route(
            port(&builder, cells, "environment"),
            &[],
            Some(heat_sink("room", 300.0).expect("sink")),
            RouteTemperature::Fixed(300.0),
        )
        .expect("environment");
    let mut system = builder.seal().expect("sealed");
    // The temperatures in the loop need a start.
    for (node, name) in [(plane, "body0"), (plane, "body1"), (patch, "heat")] {
        let at = system.port(node, name).expect("port");
        assert!(system.guess(at, &[300.0]));
    }
    Graph {
        system,
        bodies,
        patch,
        plane,
        cells,
        start: network,
        passes,
    }
}

impl Graph {
    /// One kernel step of the module's `DT`; what it gave.
    fn step(&mut self, spin: bool, before: &[RigidBodyState; 2]) -> Outcome {
        let report =
            with_cx(|cx| self.system.step(DT, cx)).unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(report.dt_s.to_bits(), DT.to_bits());
        assert_eq!(report.shortenings, 0);
        let wire = |node, name: &str| self.system.port(node, name).expect("port");
        let (heat, wrench) = (wire(self.patch, "heat"), wire(self.patch, "wrench"));
        let bodies = self
            .bodies
            .map(|node| *self.system.state::<RigidBodyState>(node).expect("body"));
        let thermal = self
            .system
            .state::<TimelineState<'static>>(self.cells)
            .expect("cells")
            .simulation();
        let cells_k = [0, 1, 2].map(|i| thermal.temperature_k(CellId(i)).expect("T"));
        let effort = if spin {
            radius() * self.system.last_output(wrench)[2]
        } else {
            self.system.last_output(wrench)[0]
        };
        let after = if spin {
            bodies[0].angular_momentum_body().z / I[0][2]
                - bodies[1].angular_momentum_body().z / I[1][2]
        } else {
            let at_contact = |i: usize| {
                bodies[i].linear_momentum_world().x / MASS[i]
                    + ARMS[i] * bodies[i].angular_momentum_body().y / I[i][1]
            };
            at_contact(0) - at_contact(1)
        };
        let outside_j = DT * self.system.last_output(wire(self.cells, "environment"))[0];
        Outcome {
            interface_k: self.system.last_input(heat)[0],
            plane_k: self.system.last_output(wire(self.plane, "source"))[0],
            heat_j: -DT * self.system.last_output(heat)[0],
            effort,
            bodies,
            cells_k,
            after,
            first_law_j: kinetic(&bodies) - kinetic(before) + energy(thermal)
                - energy(&self.start)
                - outside_j,
        }
    }
}

/// `a` differs from `b` by more than `close` allows.
fn misses(a: f64, b: f64) -> bool {
    (a - b).abs() > 3e-9 * b.abs().max(1.0)
}

/// The module tests' base case on the graph: held to the module's candidate
/// and to the tests' scalar root.
fn the_graph_is_the_module(spin: bool) {
    let pair = Pair::new(
        if spin { 0.0 } else { SPEED },
        if spin { SPIN } else { 0.0 },
    );
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let module = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let (interface_k, cells_k, heat_j, effort) = reference(spin, DT);

    let card = owner.clone();
    let mut graph = declare(&pair, spin, move |kelvin| on_card(&card, kelvin), ARMS);
    let got = graph.step(spin, &pair.bodies);
    println!(
        "Hertz {}: graph interface {} K, heat {} J, effort {}; module {} K, {} J; root {} K, {} J, {}",
        if spin { "spin" } else { "slide" },
        got.interface_k,
        got.heat_j,
        got.effort,
        module.heat.interfaces[0].interface_temperature_k,
        module.mechanical.dissipated_energy_j,
        interface_k,
        heat_j,
        effort
    );
    print_passes(if spin { "spin" } else { "slide" }, &graph.passes.borrow());

    // The interface temperature, the face's and the heat.
    close(
        got.interface_k,
        module.heat.interfaces[0].interface_temperature_k,
    );
    close(got.interface_k, interface_k);
    close(got.plane_k, interface_k);
    close(got.heat_j, module.mechanical.dissipated_energy_j);
    close(got.heat_j, heat_j);
    // The three cells; cell 2 warms through the ordinary link 10.
    for (id, (&kelvin, &root)) in (0_u64..).zip(got.cells_k.iter().zip(&cells_k)) {
        close(kelvin, module.thermal.temperature_k(CellId(id)).unwrap());
        close(kelvin, root);
    }
    assert!(got.cells_k[2] > 300.0);
    // The force or spin moment, which the warmed face has lowered.
    close(got.effort, effort);
    assert!(got.effort < effort_scale(spin) * MU[0]);
    // The recoil.
    let [a, b] = got.bodies;
    let [module_a, module_b] = [module.bodies_after[0], module.bodies_after[1]];
    if spin {
        close(module.mechanical.contact.spin_moment_nm, -got.effort);
        close(got.after, module.mechanical.spin_after_rad_s);
        close(got.after, SPIN - DT * mobility(true) * effort);
        close(
            a.angular_momentum_body().z,
            module_a.angular_momentum_body().z,
        );
        close(a.angular_momentum_body().z, I[0][2] * SPIN - DT * effort);
        close(
            b.angular_momentum_body().z,
            module_b.angular_momentum_body().z,
        );
        close(b.angular_momentum_body().z, DT * effort);
    } else {
        close(module.mechanical.contact.force_n[0], -got.effort);
        close(got.after, module.mechanical.slip_velocity_after_m_s[0]);
        close(got.after, SPEED - DT * mobility(false) * effort);
        close(
            a.linear_momentum_world().x,
            module_a.linear_momentum_world().x,
        );
        close(a.linear_momentum_world().x, MASS[0] * SPEED - DT * effort);
        close(
            b.linear_momentum_world().x,
            module_b.linear_momentum_world().x,
        );
        close(b.linear_momentum_world().x, DT * effort);
        close(
            a.angular_momentum_body().y,
            module_a.angular_momentum_body().y,
        );
        close(a.angular_momentum_body().y, -ARMS[0] * DT * effort);
        close(
            b.angular_momentum_body().y,
            module_b.angular_momentum_body().y,
        );
        close(b.angular_momentum_body().y, ARMS[1] * DT * effort);
    }
    // The first law over the step, and the kernel's ledger.
    close(got.first_law_j, 0.0);
    close(graph.system.ledger().energy_residual_j(), 0.0);
}

#[test]
fn the_sliding_case_on_the_kernel_is_the_modules_step() {
    the_graph_is_the_module(false);
}

#[test]
fn the_spinning_case_on_the_kernel_is_the_modules_step() {
    the_graph_is_the_module(true);
}

#[test]
fn a_law_that_ignores_the_face_temperature_misses_the_module() {
    // The module needs more than two thermal trials: the coefficient falls
    // as the face warms. A patch asked at 300 K whatever the face does must
    // miss the interface temperature and the heat.
    let pair = Pair::new(SPEED, 0.0);
    let (interface_k, _, heat_j, _) = reference(false, DT);
    let card = law(0);
    let mut graph = declare(&pair, false, move |_| on_card(&card, 300.0), ARMS);
    let got = graph.step(false, &pair.bodies);
    println!(
        "Hertz slide, law held at 300 K: interface {} K against {interface_k}, heat {} J against {heat_j}",
        got.interface_k, got.heat_j
    );
    assert!(misses(got.interface_k, interface_k), "{}", got.interface_k);
    assert!(misses(got.heat_j, heat_j), "{}", got.heat_j);
    assert!(got.heat_j > heat_j && got.interface_k > interface_k);
}

#[test]
fn a_lever_without_its_arms_misses_the_bodies_turn() {
    // The module turns body a by -ARMS[0] DT effort about y and body b by
    // ARMS[1] DT effort. A lever whose arms are zero puts the force on the
    // centres only, and the bodies do not turn.
    let pair = Pair::new(SPEED, 0.0);
    let (.., effort) = reference(false, DT);
    let card = law(0);
    let mut graph = declare(&pair, false, move |kelvin| on_card(&card, kelvin), [0.0; 2]);
    let got = graph.step(false, &pair.bodies);
    for (i, turn) in [(0, -ARMS[0] * DT * effort), (1, ARMS[1] * DT * effort)] {
        let held = got.bodies[i].angular_momentum_body().y;
        println!("Hertz slide, no arms: body {i} turns {held} N m s against {turn}");
        assert!(held.abs() <= f64::MIN_POSITIVE, "{held}");
        assert!(misses(held, turn));
    }
}
