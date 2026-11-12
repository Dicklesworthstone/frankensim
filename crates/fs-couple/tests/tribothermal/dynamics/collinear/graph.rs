//! The collinear stopping and reversal cases, declared as a graph.
//!
//! `fs_couple::tribothermal::dynamics::collinear` advances two translating
//! bodies under held forces through a stop: each sliding segment on the
//! midpoint relation, a segment that would reverse cut at zero relative
//! velocity, and from rest the friction law's own static return map. Here
//! the same device is declared from parts, and nothing is written for the
//! pair:
//!
//! - the kernel's `Mass`, one for each body along the tangent axis of the
//!   case and one across it, under the tests' held forces;
//! - fs-tribo's `FrictionInterface` between the two masses of the axis,
//!   with the Coulomb coefficients the tests' interface card holds, behind
//!   the kernel's load line;
//! - the two finite bodies as `UniformEnthalpyBody` on the tests' own
//!   caloric charts, with the masses of the tests' thermal cells;
//! - a heat route that divides the friction heat between them by the
//!   tests' partition.
//!
//! The module is the reference and stays. The criterion is its report: the
//! instant of the stop, the velocities and displacements of both bodies,
//! the external work, the kinetic energy, the heat deposited and the two
//! temperatures of its own run of each case, and the first law over the
//! whole run. The interface arrests through the kernel's event path, as the
//! module does: the step that holds the stop is cut once, to end there.
//! `elements::Slider`, which arrests inside its step, is run on the second
//! case beside it and misses the module's numbers by the first-order error
//! of that arrest.
//!
//! The graph has no equivalent of: coefficients that follow a material card
//! at the trial temperatures and speed (the module's third test), the
//! thermal step as a backward-Euler solve of a cell network with other
//! sources, or two tangent axes solved together.

// Exact float comparisons here assert what is fixed exactly: masses that
// are the cases' own, a pair that never slid and took no heat, a step that
// was not cut.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::component::UniformEnthalpyBody;
use fs_system::elements::{
    Force, HeatCapacity, HeatNode, Mass, Motion, Slider, Slip, SourceAccount, Spring, Waveform,
};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder};
use fs_tribo::tangential::component::{FrictionInterface, FrictionState, InterfaceLoad};
use fs_tribo::{DryTangentialState, FrictionLaw, TriboError};

use super::*;

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

/// The normal load of every case [N], and the coefficient its interface
/// card holds, static and kinetic.
const LOAD_N: f64 = 10.0;
const MU: f64 = 0.2;

/// One of the module's own test cases.
struct Case {
    /// The first test resolves its bodies from material regions (2 kg and
    /// 3 kg of 1 m3); the others take the masses of their thermal cells.
    sourced: bool,
    masses_kg: [f64; 2],
    /// Specific heat of each body [J/(kg K)].
    cp_j_kg_k: [f64; 2],
    /// Share of the friction heat each body receives.
    partition: [f64; 2],
    velocities_m_s: [[f64; 2]; 2],
    forces_n: [[f64; 2]; 2],
    /// The tangent axis the relative motion is on.
    axis: usize,
    /// The speed the contact is declared with at the start [m/s].
    declared_m_s: f64,
    until_s: f64,
}

/// `sourced_bodies_stop_at_the_physical_time_and_all_relative_kinetic_energy_heats_them`.
const STOPS: Case = Case {
    sourced: true,
    masses_kg: [2.0, 3.0],
    cp_j_kg_k: [2.0, 3.0],
    partition: [0.4, 0.6],
    velocities_m_s: [[1.0, 0.0], [0.0; 2]],
    forces_n: [[0.0; 2]; 2],
    axis: 0,
    declared_m_s: 1.0,
    until_s: 1.0,
};

/// `one_call_stops_then_reverses_with_exact_force_work_and_unsigned_sliding_heat`.
const REVERSES: Case = Case {
    sourced: false,
    masses_kg: [1.0, 1.0],
    cp_j_kg_k: [2.0, 3.0],
    partition: [0.5, 0.5],
    velocities_m_s: [[1.0, 0.0], [0.0; 2]],
    forces_n: [[-6.0, 0.0], [0.0; 2]],
    axis: 0,
    declared_m_s: 1.0,
    until_s: 0.5,
};

/// `common_transverse_acceleration_and_static_relative_motion_have_correct_work`.
const HOLDS: Case = Case {
    sourced: false,
    masses_kg: [2.0, 3.0],
    cp_j_kg_k: [2.0, 3.0],
    partition: [0.5, 0.5],
    velocities_m_s: [[1.0, 0.0], [1.0, 0.0]],
    forces_n: [[4.0, 1.0], [6.0, 0.0]],
    axis: 1,
    declared_m_s: 0.0,
    until_s: 0.2,
};

/// What a run came to, the module's or a graph's.
#[derive(Debug, Clone, PartialEq)]
struct Outcome {
    /// The instants at which sliding ended at rest [s].
    stops_s: Vec<f64>,
    velocities_m_s: [[f64; 2]; 2],
    displacements_m: [[f64; 2]; 2],
    kinetic_change_j: f64,
    external_work_j: f64,
    /// Heat the two bodies took: the rise of their enthalpy [J].
    heat_j: [f64; 2],
    temperatures_k: [f64; 2],
    /// Kinetic and thermal energy gained, less the external work [J].
    residual_j: f64,
}

/// The caloric charts and thermal cells of a case, and the module's own run
/// of it.
fn module(case: &Case) -> (Outcome, [fs_material::phase::EquilibriumEnthalpyPhaseCurve; 2], [f64; 2]) {
    let point = context();
    let mut bodies = Vec::new();
    let (cards, curves): (Vec<_>, Vec<_>) = (0..2)
        .map(|side| {
            if case.sourced {
                let (card, curve, body) = sourced_body(case.masses_kg[side], case.cp_j_kg_k[side]);
                bodies.push(body);
                (card, curve)
            } else {
                common::material(1.0, case.cp_j_kg_k[side], 0.0)
            }
        })
        .unzip();
    let (Ok(cards), Ok(curves)) = (
        <[MaterialCard; 2]>::try_from(cards),
        <[ResolvedHeatCapacityCurve; 2]>::try_from(curves),
    ) else {
        panic!("two cards and two charts");
    };
    let thermal = if case.sourced {
        inventory_with_geometry(&curves, case.masses_kg, [1.0, 1.0])
    } else {
        inventory(&curves, case.masses_kg)
    };
    let card = interface_with_static(&cards, MU, "speed", SPEED, [(0.0, MU), (4.0, MU)], 2000.0);
    let partition = HeatPartition::new(case.partition[0], case.partition[1], 0.0).unwrap();
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(&cards, &curves, &card, &point, partition),
        LOAD_N,
        case.declared_m_s,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let pair = if case.sourced {
        let Ok(bodies) = <[_; 2]>::try_from(bodies) else {
            panic!("two bodies");
        };
        SlidingPair::from_material_bodies(&thermal, contact, bodies, case.velocities_m_s).unwrap()
    } else {
        SlidingPair::new(&thermal, contact, case.velocities_m_s).unwrap()
    };
    assert_eq!(pair.masses_kg(), case.masses_kg);
    let result = pair
        .step_collinear_candidate(
            &thermal,
            case.until_s,
            MotionDrive {
                normal_load_n: LOAD_N,
                body_forces_n: case.forces_n,
            },
            event_controls(case.axis),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let report = result.report();
    let cell = |side: usize| &thermal.cells()[&CELLS[side]];
    let after = |side: usize| &result.thermal().cells()[&CELLS[side]];
    let heat_j =
        [0, 1].map(|side| cell(side).mass_kg * (after(side).specific_h_j_kg - cell(side).specific_h_j_kg));
    // The module's own accounts agree with its states.
    close(heat_j[0] + heat_j[1], report.deposited_heat_j);
    close(heat_j[0] + heat_j[1], report.thermal_energy_change_j);
    let outcome = Outcome {
        stops_s: report
            .segments
            .iter()
            .filter(|segment| segment.reached_stop)
            .map(|segment| segment.end_time_s)
            .collect(),
        velocities_m_s: result.pair().velocities_mps(),
        displacements_m: report.body_displacements_m,
        kinetic_change_j: report.kinetic_energy_change_j,
        external_work_j: report.external_work_j,
        heat_j,
        temperatures_k: [0, 1].map(|side| result.thermal().temperature_k(CELLS[side]).unwrap()),
        residual_j: report.total_energy_residual_j,
    };
    let charts = [0, 1].map(|side| curves[side].curve().clone());
    let enthalpies = [0, 1].map(|side| cell(side).specific_h_j_kg);
    (outcome, charts, enthalpies)
}

/// What rubs between the two masses of the axis.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Rubbing {
    /// fs-tribo's interface: a stop cuts the step (the event path).
    Interface,
    /// The kernel's ideal slider: a stop is taken inside its step.
    Slider,
}

type LawAt = fn(f64) -> Result<FrictionLaw, TriboError>;
type Rubbed = FrictionState<DryTangentialState>;

// A law source: its type is the interface's. So for the others below.
#[allow(clippy::unnecessary_wraps)]
fn coulomb(_: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Coulomb {
        static_mu: MU,
        kinetic_mu: MU,
    })
}

/// The load of the cases with the tangential budgets of the module's own
/// motion controls and its impulse tolerance.
fn pressed() -> InterfaceLoad {
    InterfaceLoad {
        normal_n: LOAD_N,
        nominal_pressure_pa: None,
        solve: motion_controls().sliding.tangential,
        impulse_tolerance_n_s: event_controls(0).impulse_tolerance_n_s,
        work_tolerance_j: 1e-12,
    }
}

struct Graph {
    system: System,
    /// The two masses along the case's axis, then the two across it.
    masses: [[NodeId; 2]; 2],
    rubbing: NodeId,
    bodies: [NodeId; 2],
    kind: Rubbing,
}

/// The graph of `case` on the module's charts, each body starting at the
/// enthalpy of its thermal cell.
fn declare(
    case: &Case,
    charts: &[fs_material::phase::EquilibriumEnthalpyPhaseCurve; 2],
    enthalpies_j_kg: [f64; 2],
    kind: Rubbing,
) -> Graph {
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let across = 1 - case.axis;
    // Row 0 along the case's axis, row 1 across it; in each, body a then b.
    let (mut masses, mut pushes) = (Vec::new(), Vec::new());
    for (row, axis) in [case.axis, across].into_iter().enumerate() {
        for side in 0..2 {
            let name = format!("{}_{}", ["a", "b"][side], ["along", "across"][row]);
            masses.push(
                builder
                    .add_with_tangent(
                        &name,
                        Mass::new(case.masses_kg[side]).unwrap(),
                        Motion {
                            position_m: 0.0,
                            velocity_m_s: case.velocities_m_s[side][axis],
                        },
                    )
                    .unwrap(),
            );
            pushes.push(
                builder
                    .add_with_tangent(
                        &format!("push_{name}"),
                        Force::new(Waveform::constant(case.forces_n[side][axis]).unwrap()),
                        SourceAccount::default(),
                    )
                    .unwrap(),
            );
        }
    }
    let masses = [[masses[0], masses[1]], [masses[2], masses[3]]];
    let pushes = [[pushes[0], pushes[1]], [pushes[2], pushes[3]]];
    // A member of a junction receives minus the flow of the port that sets
    // it: the slip of the two ports is minus the bodies' relative velocity.
    let slip = -(case.velocities_m_s[0][case.axis] - case.velocities_m_s[1][case.axis]);
    let rubbing = match kind {
        Rubbing::Interface => {
            let interface = FrictionInterface::new(coulomb as LawAt, pressed()).unwrap();
            let start = interface.start(300.0, slip).unwrap();
            builder.add_on_line("rubbing", interface, start).unwrap()
        }
        Rubbing::Slider => {
            let start = if slip == 0.0 { Slip::stuck() } else { Slip::sliding(slip) };
            builder
                .add_on_line("rubbing", Slider::new(LOAD_N, MU, MU).unwrap(), start)
                .unwrap()
        }
    };
    for side in 0..2 {
        builder
            .common_flow(
                port(&builder, masses[0][side], "p"),
                &[
                    port(&builder, rubbing, ["a", "b"][side]),
                    port(&builder, pushes[0][side], "p"),
                ],
            )
            .unwrap();
        builder
            .common_flow(port(&builder, masses[1][side], "p"), &[port(&builder, pushes[1][side], "p")])
            .unwrap();
    }
    let bodies = [0, 1].map(|side| {
        let body = UniformEnthalpyBody::new(charts[side].clone(), case.masses_kg[side]).unwrap();
        let start = charts[side].state_at_specific_enthalpy(enthalpies_j_kg[side]).unwrap();
        builder
            .add_with_tangent(&format!("body_{}", ["a", "b"][side]), body, start)
            .unwrap()
    });
    builder
        .heat_route(
            port(&builder, rubbing, "heat"),
            &[
                (port(&builder, bodies[0], "heat"), case.partition[0]),
                (port(&builder, bodies[1], "heat"), case.partition[1]),
            ],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(system.port(rubbing, "heat").unwrap(), &[300.0]));
    Graph {
        system,
        masses,
        rubbing,
        bodies,
        kind,
    }
}

impl Graph {
    fn at_rest(&self) -> bool {
        match self.kind {
            Rubbing::Interface => self.system.state::<Rubbed>(self.rubbing).unwrap().at_rest(),
            Rubbing::Slider => self.system.state::<Slip>(self.rubbing).unwrap().slip_m_s() == 0.0,
        }
    }

    fn motion(&self, row: usize, side: usize) -> Motion {
        *self.system.state::<Motion>(self.masses[row][side]).unwrap()
    }

    fn body(&self, side: usize) -> EquilibriumPhaseState {
        *self.system.state::<EquilibriumPhaseState>(self.bodies[side]).unwrap()
    }
}

/// How a graph's run went, besides where it ended.
#[derive(Debug, Default)]
struct Steps {
    taken: u32,
    cuts: u32,
    /// The most cuts of one step.
    most_cuts: u32,
    updates: u32,
}

/// Steps of at most `dt` to the case's end.
fn run(case: &Case, graph: &mut Graph, enthalpies_j_kg: [f64; 2], dt: f64) -> (Outcome, Steps) {
    let mut steps = Steps::default();
    let mut stops_s = Vec::new();
    with_cx(|cx| {
        // The last step may be a rounding of the end short.
        while case.until_s - graph.system.time_s() > 1e-12 {
            let was_at_rest = graph.at_rest();
            let ask = dt.min(case.until_s - graph.system.time_s());
            let report = graph
                .system
                .step(ask, cx)
                .unwrap_or_else(|error| panic!("dt {dt}, step {}: {error}", steps.taken));
            steps.taken += 1;
            steps.cuts += report.shortenings;
            steps.most_cuts = steps.most_cuts.max(report.shortenings);
            steps.updates += report.iterations;
            if graph.at_rest() && !was_at_rest {
                stops_s.push(report.time_s);
            }
        }
    });
    let across = 1 - case.axis;
    let mut outcome = Outcome {
        stops_s,
        velocities_m_s: [[0.0; 2]; 2],
        displacements_m: [[0.0; 2]; 2],
        kinetic_change_j: 0.0,
        external_work_j: 0.0,
        heat_j: [0.0; 2],
        temperatures_k: [0.0; 2],
        residual_j: 0.0,
    };
    for (side, before_j_kg) in enthalpies_j_kg.into_iter().enumerate() {
        for (row, axis) in [case.axis, across].into_iter().enumerate() {
            let motion = graph.motion(row, side);
            let before = case.velocities_m_s[side][axis];
            outcome.velocities_m_s[side][axis] = motion.velocity_m_s;
            outcome.displacements_m[side][axis] = motion.position_m;
            // (v - v0)(v + v0) / 2, not a difference of squares.
            outcome.kinetic_change_j += 0.5
                * case.masses_kg[side]
                * (motion.velocity_m_s - before)
                * (motion.velocity_m_s + before);
            outcome.external_work_j += case.forces_n[side][axis] * motion.position_m;
        }
        let body = graph.body(side);
        outcome.heat_j[side] =
            case.masses_kg[side] * (body.specific_enthalpy_j_kg() - before_j_kg);
        outcome.temperatures_k[side] = body.temperature_k();
    }
    outcome.residual_j = outcome.kinetic_change_j + outcome.heat_j[0] + outcome.heat_j[1]
        - outcome.external_work_j;
    (outcome, steps)
}

/// The largest difference between two outcomes, each quantity in its own
/// unit: seconds, m/s, m, J, K.
fn apart(graph: &Outcome, module: &Outcome) -> f64 {
    assert_eq!(graph.stops_s.len(), module.stops_s.len(), "{graph:?} against {module:?}");
    let mut worst = 0.0_f64;
    let mut see = |a: f64, b: f64| worst = worst.max((a - b).abs());
    for (a, b) in graph.stops_s.iter().zip(&module.stops_s) {
        see(*a, *b);
    }
    for side in 0..2 {
        for axis in 0..2 {
            see(graph.velocities_m_s[side][axis], module.velocities_m_s[side][axis]);
            see(graph.displacements_m[side][axis], module.displacements_m[side][axis]);
        }
        see(graph.heat_j[side], module.heat_j[side]);
        see(graph.temperatures_k[side], module.temperatures_k[side]);
    }
    see(graph.kinetic_change_j, module.kinetic_change_j);
    see(graph.external_work_j, module.external_work_j);
    worst
}

/// Whether an instant falls on the ends of steps of `dt` from zero.
fn on_the_grid(instant_s: f64, dt: f64) -> bool {
    let steps = instant_s / dt;
    (steps - steps.round()).abs() < 1e-9
}

#[test]
fn the_graph_gives_the_modules_stop_its_velocities_its_heat_and_its_temperatures() {
    for (name, case) in [("stops", &STOPS), ("reverses", &REVERSES), ("holds", &HOLDS)] {
        let (reference, charts, enthalpies) = module(case);
        // The step lengths put the stop of each case on a step's end and
        // inside a step: 0.6 s and 0.1 s are ends of steps of 0.05 s.
        for dt in [0.05, 0.07, 0.013, 0.2] {
            let mut graph = declare(case, &charts, enthalpies, Rubbing::Interface);
            let (outcome, steps) = run(case, &mut graph, enthalpies, dt);
            // Each quantity in its own unit. The module finds its stop to
            // 1e-10 s and 1e-10 N s, and its thermal solve closes a cell to
            // 1e-7 J; the two agree far inside that.
            let off = apart(&outcome, &reference);
            assert!(off <= 1e-10, "{name}, dt {dt}: {off:e}\n{outcome:?}\n{reference:?}");
            assert!(outcome.residual_j.abs() <= 1e-11, "{name}, dt {dt}: {:e}", outcome.residual_j);
            // A stop inside a step cuts that step once, to end there; a
            // stop at a step's end cuts nothing; nothing else is cut.
            let inside = reference.stops_s.iter().filter(|stop| !on_the_grid(**stop, dt)).count();
            assert_eq!(steps.cuts as usize, inside, "{name}, dt {dt}");
            assert!(steps.most_cuts <= 1, "{name}, dt {dt}");
            // The interface's record is the slip of its two ports: minus
            // the bodies' relative velocity, and exactly none at rest.
            let slip = graph.system.state::<Rubbed>(graph.rubbing).unwrap().slip_m_s();
            let relative =
                outcome.velocities_m_s[0][case.axis] - outcome.velocities_m_s[1][case.axis];
            assert!((slip + relative).abs() <= 1e-12, "{name}, dt {dt}: {slip} {relative}");
            println!(
                "collinear graph, {name}, dt {dt}: {} steps, {} cut, {} updates; stops {:?} \
                 against {:?}; off the module by {off:e} at most; first law off by {:e} J (the \
                 module's {:e} J)",
                steps.taken,
                steps.cuts,
                steps.updates,
                outcome.stops_s,
                reference.stops_s,
                outcome.residual_j,
                reference.residual_j,
            );
        }
    }
    // What the three cases are: a stop that lasts, a stop that is left the
    // other way under the drive, and a pair that never slides.
    let (stops, ..) = module(&STOPS);
    assert_eq!(stops.velocities_m_s[0][0], stops.velocities_m_s[1][0]);
    let (reverses, ..) = module(&REVERSES);
    assert!(reverses.velocities_m_s[0][0] < reverses.velocities_m_s[1][0]);
    let (holds, ..) = module(&HOLDS);
    assert!(holds.stops_s.is_empty() && holds.heat_j == [0.0; 2]);
}

#[test]
fn the_kernels_slider_arrests_inside_its_step_and_misses_the_modules_numbers() {
    // The planted negative for the event path. `elements::Slider` with the
    // same force is put in the same graph. It takes the stop inside the
    // step that holds it, by its own law of a step, and no step is cut.
    // On the second case that is first order in the step: at 0.07 s body
    // `a` ends 0.02 m/s too fast and 0.0166 J of heat too many are booked;
    // at 0.013 s 0.004 m/s too slow.
    // (The numbers are those of an emulation of the two arrests written
    // before either graph was run.) The interface, in the test above, is
    // within 1e-10 at both.
    let (reference, charts, enthalpies) = module(&REVERSES);
    for (dt, velocity_a, heat_j, work_j) in
        [(0.07, -1.420, 0.4366, 1.1130), (0.013, -1.396, 0.413_092, 1.069_908)]
    {
        let mut graph = declare(&REVERSES, &charts, enthalpies, Rubbing::Slider);
        let (outcome, steps) = run(&REVERSES, &mut graph, enthalpies, dt);
        assert_eq!(steps.cuts, 0);
        let heat = outcome.heat_j[0] + outcome.heat_j[1];
        assert!((outcome.velocities_m_s[0][0] - velocity_a).abs() <= 1e-9, "{outcome:?}");
        assert!((heat - heat_j).abs() <= 1e-9, "{heat}");
        assert!((outcome.external_work_j - work_j).abs() <= 1e-9, "{outcome:?}");
        // Its own first law is exact all the same: the error is in where
        // the stop is taken, not in the accounts.
        assert!(outcome.residual_j.abs() <= 1e-11, "{:e}", outcome.residual_j);
        let missed = (outcome.velocities_m_s[0][0] - reference.velocities_m_s[0][0]).abs();
        assert!(missed > 0.25 * dt, "dt {dt}: {missed}");
        println!(
            "collinear graph with the kernel's slider, dt {dt}: body a ends at {} m/s against \
             the module's {}, heat {heat} J against {}",
            outcome.velocities_m_s[0][0],
            reference.velocities_m_s[0][0],
            reference.heat_j[0] + reference.heat_j[1],
        );
    }
}

#[allow(clippy::unnecessary_wraps)]
fn dry(_: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Coulomb {
        static_mu: 0.5,
        kinetic_mu: 0.3,
    })
}

/// The least and the largest of some values.
fn range(values: &[f64]) -> (f64, f64) {
    values.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(least, most), value| {
        (least.min(*value), most.max(*value))
    })
}

/// A block that starts at the spring's rest point, at this velocity.
fn thrown(velocity_m_s: f64) -> Motion {
    Motion {
        position_m: 0.0,
        velocity_m_s,
    }
}

/// A block of 1 kg whose port carries the interface's `a`, the other end
/// `b` on a side that moves as `side_m_s` says (of a step's start), and the
/// interface's heat in a body of 40 J/K. With `spring` the block is on a
/// spring to a wall, of that stiffness [N/m] and cubic term [N/m3]. Under
/// 10 N the interface holds 5 N and slides against 3 N.
struct OnASide {
    system: System,
    block: NodeId,
    rubbing: NodeId,
    body: NodeId,
}

const BODY_J_K: f64 = 40.0;

fn on_a_side(
    start: Motion,
    spring: Option<(f64, f64)>,
    side_m_s: impl Fn(f64) -> f64 + 'static,
) -> OnASide {
    let velocity_m_s = start.velocity_m_s;
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let block = builder.add_with_tangent("block", Mass::new(1.0).unwrap(), start).unwrap();
    let interface = FrictionInterface::new(dry as LawAt, pressed()).unwrap();
    // The ports receive minus the velocities of what they are joined to.
    let slipping = interface.start(300.0, -(velocity_m_s - side_m_s(0.0))).unwrap();
    let rubbing = builder.add_on_line("rubbing", interface, slipping).unwrap();
    let mut joined = vec![port(&builder, rubbing, "a")];
    if let Some((stiffness, cubic)) = spring {
        // Stretched by where the block starts.
        let spring = builder
            .add_with_tangent("spring", Spring::new(stiffness, cubic).unwrap(), start.position_m)
            .unwrap();
        joined.push(port(&builder, spring, "a"));
        builder
            .boundary(port(&builder, spring, "b"), Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
    }
    builder.common_flow(port(&builder, block, "p"), &joined).unwrap();
    builder
        .boundary(
            port(&builder, rubbing, "b"),
            Box::new(move |t0, _, values| values[0] = -side_m_s(t0)),
        )
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    OnASide {
        system,
        block,
        rubbing,
        body,
    }
}

impl OnASide {
    fn at_rest(&self) -> bool {
        self.system.state::<Rubbed>(self.rubbing).unwrap().at_rest()
    }

    fn block(&self) -> Motion {
        *self.system.state::<Motion>(self.block).unwrap()
    }

    /// Heat the body has taken [J], by its temperature: known to what
    /// adding small rises to 300 K leaves of them, half a unit in the last
    /// place of 300 K a step at worst.
    fn heat_j(&self) -> (f64, f64) {
        let rounding_j = BODY_J_K * 300.0 * f64::EPSILON * (self.system.steps() as f64 + 1.0);
        (BODY_J_K * (self.system.state::<f64>(self.body).unwrap() - 300.0), rounding_j)
    }
}

#[test]
fn a_mass_on_a_spring_turns_at_the_closed_forms_points_and_stays_in_the_dead_band() {
    // 1 kg on a spring of 100 N/m, on ground that holds 5 N and slides
    // against 3 N, thrown from the spring's rest point at sqrt(22) m/s. Its
    // first turning point A has k A^2 / 2 + 3 N A = m v0^2 / 2: 0.44 m.
    // In each half cycle it swings about a centre 0.03 m against its
    // motion, so the turning points fall by 0.06 m: -0.38, 0.32, -0.26,
    // 0.20, -0.14, 0.08. At -0.02 m the spring pulls with 2 N, inside the
    // 5 N held, and the mass stays.
    //
    // Every turning point is a stop that the interface cuts a step at, and
    // a reversal the static map then lets go. On the midpoint rule a
    // linear spring keeps its energy exactly and a constant force does
    // exactly its work over a step, so the turning points of the discrete
    // motion are those of the closed form at any step length, to where the
    // stop is located.
    let (stiffness, v0) = (100.0_f64, 22.0_f64.sqrt());
    let turns = [0.44, -0.38, 0.32, -0.26, 0.20, -0.14, 0.08, -0.02];
    for dt in [1.0e-3, 0.7e-3] {
        let mut graph = on_a_side(thrown(v0), Some((stiffness, 0.0)), |_| 0.0);
        let (mut found, mut cuts, mut slivers) = (Vec::new(), 0, 0);
        let mut taken = 0_u32;
        with_cx(|cx| {
            while graph.system.time_s() < 4.0 {
                let was_at_rest = graph.at_rest();
                let report = graph.system.step(dt, cx).unwrap_or_else(|error| panic!("{error}"));
                taken += 1;
                cuts += report.shortenings;
                if report.dt_s < 1e-3 * dt {
                    slivers += 1;
                }
                if graph.at_rest() && !was_at_rest {
                    found.push(graph.block().position_m);
                }
            }
        });
        assert_eq!(found.len(), turns.len(), "{found:?}");
        let worst = found
            .iter()
            .zip(turns)
            .map(|(at, turn)| (at - turn).abs())
            .fold(0.0_f64, f64::max);
        // At rest in the dead band, and all the energy it was thrown with,
        // less what the spring keeps, is heat in the body.
        let ends = graph.block();
        assert!(graph.at_rest() && ends.velocity_m_s.abs() <= 1e-12, "{ends:?}");
        let (heat_j, rounding_j) = graph.heat_j();
        let expected_j = 0.5 * v0 * v0 - 0.5 * stiffness * 0.02 * 0.02;
        println!(
            "a mass on a spring through the interface, dt {dt}: turning points off by {worst:e} m \
             at most; {taken} steps, {cuts} cut, {slivers} under a thousandth of a step; heat \
             {heat_j} J against {expected_j}"
        );
        assert!(worst <= 1e-12, "{found:?}");
        assert!((heat_j - expected_j).abs() <= rounding_j + 1e-12, "{heat_j}");
        // The spring is not a body under a held force, so a stop is not
        // located in one cut (fs-tribo's component, "Where the stop is
        // found"): two cuts a turning point here, and no step after it.
        assert!(cuts <= 3 * turns.len() as u32 && slivers == 0, "{cuts} {slivers}");
    }
}

#[test]
fn a_block_on_a_belt_sticks_and_slips_with_the_closed_forms_period() {
    // Stuck, the block rides the belt until the spring pulls with the
    // static limit, at x_s = 5 N / k. Sliding, it swings about
    // x_k = 3 N / k as a free oscillator, from x_k + d at the belt's speed
    // (d = x_s - x_k) round to x_k - d at the belt's speed again, where
    // the belt catches it. One period is the ride 2 d / V and the swing
    // 2 (pi - atan(d w / V)) / w.
    //
    // Both events are located. The catch is a stop of the slip, and the
    // interface cuts the step there (the spring is not a body under a
    // held force: up to three cuts a catch). The break-away is the
    // spring's force reaching the 5 N held, inside a step: the step is
    // cut there, once, and the next one lets go from its start. So every
    // period is the closed form's to second order in the step (the swing
    // is a midpoint oscillator), 3.4 dt^2 measured, and the same from
    // cycle to cycle.
    //
    // The belt's speeds are chosen so that the ride is not a whole number
    // of steps. (At 0.1 m/s it is 0.4 s and every break-away falls on a
    // step's end at each of these step lengths, located or not. The first
    // version of this test ran there, with a break-away that was not
    // located, and printed 3.4 dt^2 as its result: an accident of the
    // numbers. Not located, at these speeds, the break-away was taken up
    // to half a step early or late and one period was up to 0.54 of a
    // step off, different from cycle to cycle: 9f4fd342b2.)
    //
    // At 0.137 m/s and 0.25 ms the third catch is the regression of the
    // stop's own row: the interface as first landed returned there the
    // force that ends the step at no slip exactly, stiff in the released
    // flow and only inside a band the kernel's differencing does not see,
    // and the step was refused after 9009 steps (`did not converge in 25
    // iterations; worst input block.p has residual -1.4865e-11`). A catch
    // is now taken under the kinetic 3 N, and leaves a slip inside the
    // impulse tolerance that the first step at rest takes out: both are
    // measured below, at every catch.
    let stiffness = 100.0_f64;
    let (omega, d) = (stiffness.sqrt(), (5.0 - 3.0) / stiffness);
    for belt in [0.0937_f64, 0.137] {
        let period =
            2.0 * d / belt + 2.0 * (std::f64::consts::PI - (d * omega / belt).atan()) / omega;
        for halvings in 0..4 {
            let dt = 2.0e-3 / f64::from(1_u32 << halvings);
            let mut graph = on_a_side(thrown(belt), Some((stiffness, 0.0)), move |_| belt);
            // The instant and place of each catch (the start counts as
            // one), and the start of each step that breaks away.
            let (mut caught, mut broke) = (vec![(0.0, 0.0)], Vec::new());
            let (mut cuts, mut most_cuts, mut slivers, mut taken) = (0, 0, 0, 0_u32);
            // Cuts of steps that start at rest: the break-aways'.
            let mut cuts_at_rest = 0;
            let mut worked_j = 0.0;
            // The largest slip a catch leaves [m/s], and the most the first
            // step at rest holds with beside the spring's force [N].
            let (mut left_m_s, mut beside_n) = (0.0_f64, 0.0_f64);
            with_cx(|cx| {
                let (at_a, at_b) = (
                    graph.system.port(graph.rubbing, "a").unwrap(),
                    graph.system.port(graph.rubbing, "b").unwrap(),
                );
                let mut after_a_catch = false;
                while graph.system.time_s() < 6.0 {
                    let was_at_rest = graph.at_rest();
                    let before_m = graph.block().position_m;
                    let report = graph
                        .system
                        .step(dt, cx)
                        .unwrap_or_else(|error| panic!("belt {belt}, dt {dt}: {error}"));
                    let record = graph.system.state::<Rubbed>(graph.rubbing).unwrap().slip_m_s();
                    if after_a_catch {
                        // Riding the belt, the block needs the spring's
                        // force at mid-step. This one step also takes out
                        // what the catch left, and then nothing is left.
                        let spring_n = stiffness * 0.5 * (before_m + graph.block().position_m);
                        let held_n = graph.system.last_output(at_a)[0].abs();
                        beside_n = beside_n.max((held_n - spring_n.abs()).abs());
                        assert_eq!(record, 0.0);
                        after_a_catch = false;
                    }
                    if graph.at_rest() && !was_at_rest {
                        left_m_s = left_m_s.max(record.abs());
                        after_a_catch = true;
                    }
                    taken += 1;
                    cuts += report.shortenings;
                    most_cuts = most_cuts.max(report.shortenings);
                    if was_at_rest {
                        cuts_at_rest += report.shortenings;
                    }
                    if report.dt_s < 1e-3 * dt {
                        slivers += 1;
                    }
                    // What the two ends gave up in the step.
                    let slip = graph.system.last_input(at_a)[0] - graph.system.last_input(at_b)[0];
                    worked_j += report.dt_s * graph.system.last_output(at_a)[0] * slip;
                    if was_at_rest && !graph.at_rest() {
                        broke.push(report.time_s - report.dt_s);
                    }
                    if graph.at_rest() && !was_at_rest {
                        caught.push((report.time_s, graph.block().position_m));
                    }
                }
            });
            assert!(caught.len() >= 7, "{caught:?}");
            // Riding from its catch at the belt's speed, the block is where
            // the spring pulls with 5 N after (0.05 m - x) / V.
            let late: Vec<f64> = broke
                .iter()
                .zip(&caught)
                .map(|(broke, (at_s, at_m))| broke - (at_s + (0.05 - at_m) / belt))
                .collect();
            // Catch to catch, after the first cycle, which starts from
            // another place.
            let off: Vec<f64> =
                caught[1..].windows(2).map(|pair| pair[1].0 - pair[0].0 - period).collect();
            let (heat_j, rounding_j) = graph.heat_j();
            // The spring's force passes the limit at k V, so the impulse
            // tolerance is this much time either way.
            let within_s = (2.0 * pressed().impulse_tolerance_n_s / (stiffness * belt)).sqrt();
            println!(
                "a block on a belt at {belt} m/s, dt {dt}: {} catches in {taken} steps, {cuts} cut \
                 ({most_cuts} at most in a step, {cuts_at_rest} for {} break-aways); the \
                 break-away is taken from {:+.2e} to {:+.2e} s off the instant the limit is \
                 reached (within {within_s:.1e} s); a period is from {:.4} to {:.4} dt^2 off \
                 {period} s, {:.1e} s apart",
                caught.len() - 1,
                broke.len(),
                range(&late).0,
                range(&late).1,
                range(&off).0 / (dt * dt),
                range(&off).1 / (dt * dt),
                range(&off).1 - range(&off).0,
            );
            assert!(late.iter().all(|late| late.abs() <= within_s), "{late:?}");
            assert_eq!(cuts_at_rest as usize, broke.len());
            assert!(off.iter().all(|off| *off > 0.0 && *off < 4.0 * dt * dt), "{off:?}");
            assert!(range(&off).1 - range(&off).0 <= 1e-7 * dt, "{off:?}");
            // What a catch leaves is inside the impulse tolerance of
            // 1e-10 N s on the block of 1 kg, and the one step that takes
            // it out holds with at most that tolerance over the step
            // beside the spring's force.
            println!(
                "   a catch leaves a slip of {left_m_s:e} m/s at most; the first step at rest \
                 holds with {beside_n:e} N at most beside the spring's force (the tolerance \
                 over the step: {:e} N)",
                1e-10 / dt
            );
            assert!(left_m_s <= 1e-10 && beside_n <= 1e-10 / dt, "{left_m_s} {beside_n}");
            // The work of the two ends is the heat in the body, to rounding.
            assert!((worked_j - heat_j).abs() <= rounding_j, "{worked_j} {heat_j}");
            assert!(most_cuts <= 3 && slivers == 0, "{most_cuts} {slivers}");
        }
    }
}

#[test]
fn a_side_that_jumps_drags_the_block_under_the_kinetic_force_until_the_located_catch() {
    // A block at rest on a side at 0.1 m/s. After 0.095 s the side moves
    // at 0.3 m/s, at once: the step after the jump starts with a slip of
    // 0.2 m/s, which the interface reads off the change of the flow its
    // prescribed end is given. 3 N drag the block of 1 kg after the side,
    // and it is caught 0.2 / 3 s after the jump, inside a step: one cut.
    let dt = 0.01_f64;
    let mut graph = on_a_side(thrown(0.1), None, |t0| if t0 > 0.095 { 0.3 } else { 0.1 });
    let (mut caught_at, mut cuts) = (None, 0);
    with_cx(|cx| {
        for step in 0..30 {
            let was_at_rest = graph.at_rest();
            let report =
                graph.system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
            cuts += report.shortenings;
            if graph.at_rest() && !was_at_rest {
                assert!(caught_at.is_none());
                caught_at = Some(report.time_s);
            }
            if step < 10 {
                assert!(graph.at_rest(), "step {step}");
            }
        }
    });
    let caught_at = caught_at.expect("the block is caught");
    assert!((caught_at - (0.1 + 0.2 / 3.0)).abs() <= 1e-12, "{caught_at}");
    assert_eq!(cuts, 1);
    assert!((graph.block().velocity_m_s - 0.3).abs() <= 1e-12);
    // The heat of the drag: 3 N over the slip's mean of 0.1 m/s.
    let (heat_j, rounding_j) = graph.heat_j();
    assert!((heat_j - 3.0 * 0.1 * 0.2 / 3.0).abs() <= rounding_j, "{heat_j}");
}

#[test]
fn a_stop_beside_a_spring_that_stiffens_is_reached_in_a_few_short_steps() {
    // Where the stop is found. The interface locates a stop from the line
    // of the step it is asked for, as if what it is joined to were bodies
    // under held forces. A spring's force changes within the step, so the
    // located instant is only as good as the step it is located from, and
    // the cut step ends near the stop: beyond it, or short of it. Beyond
    // it, it is cut again, from a line that is now close. Short of it, it
    // is taken, sliding, and the stop is left to the next step, which is
    // asked at full length again.
    //
    // A block of 1 kg against 3 N and a spring with a cubic term, in steps
    // of 0.07, 0.01 and 0.001 s, thrown outwards from the spring's rest
    // point and let go 0.3 m out moving inwards: 4, 3 and 2 cuts and 5, 3
    // and 2, all in the step that holds the stop, and no step after it.
    // The block comes to rest at the same place at all three step lengths
    // (to 2e-13 m; printed).
    //
    // With two bodies (the next test) the cut step ends short of the stop.
    for (name, start, spring) in [
        ("thrown outwards", thrown(1.0), (1.0, 200.0)),
        (
            "moving inwards",
            Motion {
                position_m: -0.3,
                velocity_m_s: 0.5,
            },
            (1.0, 60.0),
        ),
    ] {
        let mut stops = Vec::new();
        for dt in [0.07, 0.01, 0.001] {
            let mut graph = on_a_side(start, Some(spring), |_| 0.0);
            let (mut cuts, mut slivers, mut taken, mut stopped_at) = (0, 0, 0_u32, None);
            with_cx(|cx| {
                while stopped_at.is_none() {
                    let report =
                        graph.system.step(dt, cx).unwrap_or_else(|error| panic!("{error}"));
                    taken += 1;
                    cuts += report.shortenings;
                    if report.dt_s < 1e-3 * dt {
                        slivers += 1;
                    }
                    if graph.at_rest() {
                        stopped_at = Some(report.time_s);
                    }
                    assert!(taken < 1000);
                }
                // It holds there: the spring pulls with less than 5 N.
                for _ in 0..5 {
                    let report = graph.system.step(dt, cx).unwrap();
                    assert!(graph.at_rest() && report.shortenings == 0);
                }
            });
            let ends = graph.block();
            println!(
                "a stop beside a stiffening spring, {name}, dt {dt}: at rest at {} s and {} m \
                 after {taken} steps, {cuts} cut, {slivers} of the steps under a thousandth of \
                 a step",
                stopped_at.unwrap(),
                ends.position_m
            );
            assert!(ends.velocity_m_s.abs() <= 1e-9, "{ends:?}");
            assert!(cuts <= 5 && slivers == 0, "{name}, dt {dt}: {cuts} cut, {slivers} slivers");
            stops.push((stopped_at.unwrap(), ends.position_m));
        }
        // The three runs are three discretisations of one motion: the
        // instant of the stop converges, and its place is the same.
        assert!((stops[0].0 - stops[2].0).abs() < 2e-3, "{name}: {stops:?}");
        assert!((stops[1].0 - stops[2].0).abs() < 1e-4, "{name}: {stops:?}");
        assert!((stops[0].1 - stops[2].1).abs() < 1e-12, "{name}: {stops:?}");
        assert!((stops[1].1 - stops[2].1).abs() < 1e-12, "{name}: {stops:?}");
    }
}

#[test]
fn a_stop_of_two_bodies_beside_that_spring_is_left_to_three_short_steps() {
    // The first case of the module (2 kg at 1 m/s on 3 kg, 2 N of
    // friction) with body `a` on a spring of 1 N/m and 40 N/m3, in steps
    // of 0.07 s. The step that holds the stop is cut twice and ends short
    // of it, sliding, with a millionth of a metre per second left. Each
    // step after it is asked at 0.07 s again, is cut once, and ends with
    // about a fortieth of the slip it began with: 3.6e-7, 9.1e-9 and
    // 2.3e-10 s long, the last one inside the impulse tolerance and at
    // rest.
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let blocks = [(2.0, 1.0), (3.0, 0.0)].map(|(mass_kg, velocity_m_s)| {
        builder
            .add_with_tangent(
                &format!("block_{mass_kg}"),
                Mass::new(mass_kg).unwrap(),
                thrown(velocity_m_s),
            )
            .unwrap()
    });
    let interface = FrictionInterface::new(coulomb as LawAt, pressed()).unwrap();
    let slipping = interface.start(300.0, -1.0).unwrap();
    let rubbing = builder.add_on_line("rubbing", interface, slipping).unwrap();
    let spring = builder.add_with_tangent("spring", Spring::new(1.0, 40.0).unwrap(), 0.0).unwrap();
    builder
        .common_flow(
            port(&builder, blocks[0], "p"),
            &[port(&builder, rubbing, "a"), port(&builder, spring, "a")],
        )
        .unwrap();
    builder
        .common_flow(port(&builder, blocks[1], "p"), &[port(&builder, rubbing, "b")])
        .unwrap();
    builder
        .boundary(port(&builder, spring, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    let (mut cuts, mut taken, mut lengths, mut slips) = (0, 0_u32, Vec::new(), Vec::new());
    with_cx(|cx| {
        while !system.state::<Rubbed>(rubbing).unwrap().at_rest() {
            let report = system.step(0.07, cx).unwrap_or_else(|error| panic!("{error}"));
            taken += 1;
            cuts += report.shortenings;
            if report.shortenings > 0 {
                lengths.push(report.dt_s);
                slips.push(system.state::<Rubbed>(rubbing).unwrap().slip_m_s());
            }
            assert!(taken < 1000);
        }
    });
    let velocities = blocks.map(|block| system.state::<Motion>(block).unwrap().velocity_m_s);
    println!(
        "a stop of two bodies beside a stiffening spring, dt 0.07: at rest at {} s after {taken} \
         steps, {cuts} cut; the steps that were cut are {lengths:?} s long and end with slips \
         of {slips:?} m/s",
        system.time_s()
    );
    assert!((velocities[0] - velocities[1]).abs() <= 1e-9, "{velocities:?}");
    // Four cut steps, the first near the full length and the others far
    // shorter; five cuts in all; the slip falls by more than ten from each
    // to the next, and the last is at rest with a slip inside the impulse
    // tolerance of 1e-10 N s on the reduced mass of 1.2 kg.
    assert_eq!((lengths.len(), cuts), (4, 5), "{lengths:?}");
    assert!(lengths[0] > 0.06 && lengths[1..].iter().all(|length| *length < 1e-6), "{lengths:?}");
    assert!(slips.windows(2).all(|pair| pair[1].abs() < 0.1 * pair[0].abs()), "{slips:?}");
    assert!(slips[3].abs() <= 1e-10 / 1.2, "{slips:?}");
}

#[test]
fn a_side_that_accelerates_is_outside_the_claim() {
    // Not covered, and what it costs. A side that speeds up at 1 m/s2
    // under a block of 1 kg, which 1 N of the 5 N held would carry. A
    // prescribed end is held over its step, so the side is a staircase,
    // and each step of it is a jump: the block is dragged under the 3 N
    // of sliding until it is caught, which takes a third of the step
    // before (the jump is the acceleration times that step, and 3 N on
    // 1 kg make up 3 m/s in a second), and the step is cut there. Then
    // the side has moved again, by a third as much. So every step asked
    // for is followed by a train of steps, each a third of the one before,
    // until what is left is inside the impulse tolerance: more than ten
    // steps for each one asked, and for each train 1.125 times
    // m (a dt)^2 / 2 of heat (the squares of the thirds add up to 1.125)
    // where a carried block makes none. The block ends a step behind the
    // side. (As first landed, with the stop's stiff row, such a run ended
    // in a step that did not converge.)
    for dt in [0.01_f64, 0.001] {
        let mut graph = on_a_side(thrown(0.0), None, |t0| t0);
        let mut lengths = Vec::new();
        with_cx(|cx| {
            while graph.system.time_s() < 0.2 {
                let report = graph.system.step(dt, cx).unwrap_or_else(|error| panic!("{error}"));
                lengths.push(report.dt_s);
            }
        });
        for pair in lengths[..8].windows(2) {
            assert!((pair[1] / pair[0] - 1.0 / 3.0).abs() <= 1e-9, "dt {dt}: {lengths:?}");
        }
        let whole = lengths.iter().filter(|length| **length == dt).count();
        let (heat_j, _) = graph.heat_j();
        let per_step_j = 0.5 * dt * dt;
        let (time_s, block_m_s) = (graph.system.time_s(), graph.block().velocity_m_s);
        println!(
            "a side that accelerates, dt {dt}: {} steps reach {time_s} s, {whole} of them whole; \
             the block at {block_m_s} m/s; heat {heat_j:e} J, {:.3} of m (a dt)^2 / 2 for each \
             whole step",
            lengths.len(),
            heat_j / (whole as f64 * per_step_j)
        );
        assert!(lengths.len() as f64 > 10.0 * 0.2 / dt, "dt {dt}: {}", lengths.len());
        assert!(heat_j > whole as f64 * per_step_j && heat_j < 1.2 * whole as f64 * per_step_j);
        assert!((block_m_s - time_s).abs() <= 1.5 * dt, "dt {dt}: {block_m_s} at {time_s}");
    }
}

/// `mu = 0.3 - 0.002 (T - 300 K)`, static and kinetic.
#[allow(clippy::unnecessary_wraps)]
fn weakening(kelvin: f64) -> Result<FrictionLaw, TriboError> {
    let mu = 0.3 - 0.002 * (kelvin - 300.0);
    Ok(FrictionLaw::Coulomb {
        static_mu: mu,
        kinetic_mu: mu,
    })
}

/// Two sides held at a slip of 0.1 m/s under 10 N, the interface's heat on
/// a massless plane between a body of 2 J/K behind 0.5 W/K and one of
/// 5 J/K behind 0.2 W/K. Returns the time mean of the plane's temperature
/// as the interface received it over `until_s`, the two bodies'
/// temperatures at the end, and the largest difference between a step's
/// heat rate and the law's traction at that step's temperature times the
/// slip.
fn warmed_plane(law_at: LawAt, dt: f64, until_s: f64) -> (f64, [f64; 2], f64) {
    let (capacities, conductances) = ([2.0, 5.0], [0.5, 0.2]);
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let interface = FrictionInterface::new(law_at, pressed()).unwrap();
    let start = interface.start(300.0, 0.1).unwrap();
    let rubbing = builder.add_on_line("rubbing", interface, start).unwrap();
    for (name, flow) in [("a", 0.1), ("b", 0.0)] {
        builder
            .boundary(port(&builder, rubbing, name), Box::new(move |_, _, values| values[0] = flow))
            .unwrap();
    }
    let plane = builder
        .add_with_tangent("plane", HeatNode::new(&conductances).unwrap(), ())
        .unwrap();
    let bodies = [0, 1].map(|side| {
        let body = builder
            .add_with_tangent(
                &format!("body{side}"),
                HeatCapacity::new(capacities[side]).unwrap(),
                300.0,
            )
            .unwrap();
        builder
            .common_effort(
                port(&builder, body, "heat"),
                &[port(&builder, plane, &format!("body{side}"))],
            )
            .unwrap();
        body
    });
    builder
        .common_effort(port(&builder, plane, "source"), &[port(&builder, rubbing, "heat")])
        .unwrap();
    let mut system = builder.seal().unwrap();
    let warmth = system.port(rubbing, "heat").unwrap();
    assert!(system.guess(warmth, &[300.0]));
    for side in 0..2 {
        let at = system.port(plane, &format!("body{side}")).unwrap();
        assert!(system.guess(at, &[300.0]));
    }
    let (mut mean_k, mut off_w) = (0.0, 0.0_f64);
    with_cx(|cx| {
        while until_s - system.time_s() > 1e-9 {
            let report = system.step(dt, cx).unwrap_or_else(|error| panic!("{error}"));
            assert_eq!(report.shortenings, 0);
            let kelvin = system.last_input(warmth)[0];
            mean_k += report.dt_s * kelvin;
            let FrictionLaw::Coulomb { kinetic_mu, .. } = law_at(kelvin).unwrap() else {
                unreachable!()
            };
            off_w = off_w.max((system.last_output(warmth)[0] + kinetic_mu * LOAD_N * 0.1).abs());
        }
    });
    let ends = bodies.map(|body| *system.state::<f64>(body).unwrap());
    (mean_k / until_s, ends, off_w)
}

#[test]
fn the_interfaces_mean_temperature_is_the_closed_form_of_two_lumped_bodies() {
    // A coefficient of 0.3 under 10 N at 0.1 m/s: a constant 0.3 W into
    // the plane. With g_i = G_i / (G_a + G_b) and the two conductances in
    // series, G_s, the bodies' temperature difference relaxes as
    // D(t) = D_inf (1 - exp(-t / tau)), tau = 1 / (G_s (1 / C_a + 1 / C_b)),
    // D_inf = P tau (g_a / C_a - g_b / C_b), about their capacity mean
    // T_0 + P t / (C_a + C_b). The plane is at
    // g_a T_a + g_b T_b + P / (G_a + G_b).
    let (c, g, power) = ([2.0_f64, 5.0], [0.5_f64, 0.2], 0.3 * LOAD_N * 0.1);
    let (sum_c, sum_g) = (c[0] + c[1], g[0] + g[1]);
    let series = g[0] * g[1] / sum_g;
    let tau = 1.0 / (series * (1.0 / c[0] + 1.0 / c[1]));
    let d_inf = power * tau * (g[0] / sum_g / c[0] - g[1] / sum_g / c[1]);
    let until = 20.0_f64;
    // The plane's weight on the difference, and the difference's time mean.
    let weight = (g[0] * c[1] - g[1] * c[0]) / (sum_g * sum_c);
    let mean_d = d_inf * (1.0 - tau / until * (1.0 - (-until / tau).exp()));
    let mean_plane = 300.0 + 0.5 * power * until / sum_c + weight * mean_d + power / sum_g;
    let d_end = d_inf * (1.0 - (-until / tau).exp());
    let ends = [
        300.0 + power * until / sum_c + c[1] / sum_c * d_end,
        300.0 + power * until / sum_c - c[0] / sum_c * d_end,
    ];
    let mut errors = Vec::new();
    for dt in [0.2, 0.1, 0.05] {
        let (mean_k, bodies, off_w) = warmed_plane(coulomb_03, dt, until);
        assert!(off_w <= 1e-15, "{off_w}");
        let error = (mean_k - mean_plane).abs();
        let at_ends = (bodies[0] - ends[0]).abs().max((bodies[1] - ends[1]).abs());
        println!(
            "two lumped bodies, dt {dt}: mean interface temperature {mean_k} K against \
             {mean_plane}, off by {error:e}; bodies end at {bodies:?} against {ends:?}, off by \
             {at_ends:e}"
        );
        errors.push((error, at_ends));
    }
    // Second order in the step, both.
    for pair in errors.windows(2) {
        assert!(pair[1].0 < 0.3 * pair[0].0 && pair[1].1 < 0.3 * pair[0].1, "{errors:?}");
    }
    assert!(errors[2].0 <= 1e-5 && errors[2].1 <= 1e-5, "{errors:?}");

    // A law that follows the temperature is asked at the plane's
    // temperature of each step: the heat of a step is the traction at the
    // temperature the interface received in it, times the slip. The plane
    // warms and the heat falls.
    let (mean_k, bodies, off_w) = warmed_plane(weakening, 0.1, until);
    assert!(off_w <= 1e-15, "{off_w}");
    assert!(mean_k > 300.2 && mean_k < mean_plane, "{mean_k}");
    assert!(bodies[0] < ends[0] && bodies[1] < ends[1], "{bodies:?}");
}

#[allow(clippy::unnecessary_wraps)]
fn coulomb_03(_: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Coulomb {
        static_mu: 0.3,
        kinetic_mu: 0.3,
    })
}

/// A block of 1 kg at rest on the ground under a push, with the interface
/// (5 N held, 3 N sliding) or the kernel's slider with the same forces
/// between it and the ground, and the heat left at 300 K.
fn pushed_on_the_ground(push: Waveform, kind: Rubbing) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let block = builder.add_with_tangent("block", Mass::new(1.0).unwrap(), thrown(0.0)).unwrap();
    let pusher = builder
        .add_with_tangent("push", Force::new(push), SourceAccount::default())
        .unwrap();
    let rubbing = match kind {
        Rubbing::Interface => {
            let interface = FrictionInterface::new(dry as LawAt, pressed()).unwrap();
            let start = interface.start(300.0, 0.0).unwrap();
            builder.add_on_line("rubbing", interface, start).unwrap()
        }
        Rubbing::Slider => builder
            .add_on_line("rubbing", Slider::new(LOAD_N, 0.5, 0.3).unwrap(), Slip::stuck())
            .unwrap(),
    };
    builder
        .common_flow(
            port(&builder, block, "p"),
            &[port(&builder, rubbing, "a"), port(&builder, pusher, "p")],
        )
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    (builder.seal().unwrap(), block, rubbing)
}

/// Steps of `dt` until `until_s` is passed: the start of the step that
/// lets go, the cuts, and the block and the time at the end.
fn let_go(push: Waveform, kind: Rubbing, dt: f64, until_s: f64) -> (Option<f64>, u32, Motion, f64) {
    let (mut system, block, rubbing) = pushed_on_the_ground(push, kind);
    let at_rest = |system: &System| match kind {
        Rubbing::Interface => system.state::<Rubbed>(rubbing).unwrap().at_rest(),
        Rubbing::Slider => system.state::<Slip>(rubbing).unwrap().is_stuck(),
    };
    let (mut broke, mut cuts) = (None, 0);
    with_cx(|cx| {
        while system.time_s() < until_s {
            let was_at_rest = at_rest(&system);
            let report = system.step(dt, cx).unwrap_or_else(|error| panic!("dt {dt}: {error}"));
            cuts += report.shortenings;
            if was_at_rest && !at_rest(&system) && broke.is_none() {
                broke = Some(report.time_s - report.dt_s);
            }
        }
    });
    (broke, cuts, *system.state::<Motion>(block).unwrap(), system.time_s())
}

#[test]
fn a_push_that_rises_through_the_limit_lets_go_at_the_instant_it_reaches_it() {
    // The break-away is located. A block of 1 kg at rest under a push
    // that rises at 10 N/s: the 5 N the interface holds are reached at
    // 0.5 s, inside a step at each of these step lengths. The step is cut
    // there, once; the cut step is held; the next one lets go from its
    // start. From then on the block slides against 3 N under a push that
    // goes on rising, so its velocity is 2 tau + 5 tau^2 and its place
    // tau^2 + 5/3 tau^3 after a time tau. The midpoint rule gives the
    // velocity of a linearly rising force exactly, and the place to
    // second order (a quarter of dt^2 here, measured and asserted).
    //
    // The planted negative: the kernel's slider with the same forces in
    // the same graph decides each step by its mean, and lets go at the
    // start of the first step whose middle is beyond the limit: a third
    // to nearly half a step off at these three lengths, with no cut, and
    // the velocity off by about twice that lateness.
    for dt in [0.03_f64, 0.007, 0.0011] {
        let rising = || Waveform::table(vec![(0.0, 0.0), (2.0, 20.0)]).unwrap();
        let (broke, cuts, ends, time_s) = let_go(rising(), Rubbing::Interface, dt, 0.8);
        let broke = broke.expect("it lets go");
        let tau = time_s - 0.5;
        let (v, x) = (2.0 * tau + 5.0 * tau * tau, tau * tau + 5.0 / 3.0 * tau * tau * tau);
        assert!((broke - 0.5).abs() <= 1e-12, "dt {dt}: {broke}");
        assert_eq!(cuts, 1, "dt {dt}");
        assert!((ends.velocity_m_s - v).abs() <= 1e-12, "dt {dt}: {} {v}", ends.velocity_m_s);
        let behind = (ends.position_m - x) / (dt * dt);
        assert!(behind > 0.2 && behind < 0.3, "dt {dt}: {behind}");

        let (late, cuts, slid, time_s) = let_go(rising(), Rubbing::Slider, dt, 0.8);
        let late = (late.expect("it lets go") - 0.5) / dt;
        let tau = time_s - 0.5;
        let off = slid.velocity_m_s - (2.0 * tau + 5.0 * tau * tau);
        println!(
            "a push rising through the limit, dt {dt}: the interface lets go {:+.1e} s off 0.5 s \
             with one cut, its velocity {:+.1e} m/s off the closed form and its place {behind:.3} \
             dt^2; the kernel's slider lets go {late:+.3} of a step off, its velocity {off:+.3e} m/s",
            broke - 0.5,
            ends.velocity_m_s - v
        );
        assert_eq!(cuts, 0, "dt {dt}");
        assert!(late.abs() <= 0.5 && late.abs() > 0.3, "dt {dt}: {late}");
        assert!(off.abs() > 0.5 * dt, "dt {dt}: {off}");
    }
}

#[test]
fn a_push_that_curves_is_let_go_within_the_impulse_tolerance_of_the_instant() {
    // The load at a step's two ends is read off the straight line through
    // the means of the last step and of this one, which is exact only for
    // a load that is straight in time. Under 8 N sin(pi t), which reaches
    // the 5 N held at asin(5 / 8) / pi, the first cut lands near the
    // instant; the cut step is located again from its own, shorter, line,
    // and so on until what letting go there leaves is inside the impulse
    // tolerance: three to five cuts at these step lengths. The load
    // passes the limit at 19.6 N/s, so a tolerance of 1e-10 N s is 3.2e-6 s
    // either way (half the rate times the time squared), and the instant
    // is within that. The kernel's slider is within half a step.
    let instant = (5.0_f64 / 8.0).asin() / std::f64::consts::PI;
    let rate = 8.0 * std::f64::consts::PI * (1.0 - (5.0_f64 / 8.0).powi(2)).sqrt();
    let within_s = (2.0 * pressed().impulse_tolerance_n_s / rate).sqrt();
    for dt in [0.03_f64, 0.015, 0.0075, 0.00375] {
        let curving = || Waveform::sine(8.0, 0.5, 0.0, 0.0).unwrap();
        let (broke, cuts, ..) = let_go(curving(), Rubbing::Interface, dt, 0.4);
        let off_s = broke.expect("it lets go") - instant;
        let (late, none, ..) = let_go(curving(), Rubbing::Slider, dt, 0.4);
        let late = (late.expect("it lets go") - instant) / dt;
        println!(
            "a push that curves, dt {dt}: the interface lets go {off_s:+.3e} s off the instant \
             (within {within_s:.2e} s) after {cuts} cuts; the kernel's slider {late:+.3} of a step off"
        );
        assert!(off_s.abs() <= within_s, "dt {dt}: {off_s}");
        assert!((1..=5).contains(&cuts), "dt {dt}: {cuts}");
        assert_eq!(none, 0);
        assert!(late.abs() <= 0.5 && late.abs() > 0.1, "dt {dt}: {late}");
    }
}
