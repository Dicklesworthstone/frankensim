//! G0: the composition kernel on small lumped systems.
//!
//! Every device here is built from one-port elements and connections. No test
//! contains stepping code for a combination: the reference trajectories are
//! computed by a separate monolithic formula and compared with what the
//! kernel produces from the graph.

// Exact float comparisons here assert exact zeros and unchanged clocks.
#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;
use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_ports::{PortKind, StableId};
use fs_system::{
    AdmissionError, NodeId, Phase, RouteTemperature, SolveConfig, System, SystemBuilder,
    SystemError, heat_sink,
};

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

fn id(text: &str) -> StableId {
    StableId::new(text).expect("stable id")
}

fn mechanical(name: &str, role: PortRole) -> PortSpec {
    PortSpec::power(
        id(name),
        PortKind::MechanicalForceVelocity,
        role,
        NonZeroUsize::MIN,
    )
    .expect("power port")
}

fn heat(name: &str, role: PortRole) -> PortSpec {
    PortSpec::heat(id(name), role).expect("heat port")
}

fn roundoff(scale: f64) -> f64 {
    32.0 * f64::EPSILON * scale
}

/// Inertia: receives force, returns the mean velocity of the step.
struct Inertia {
    mass_kg: f64,
    ports: [PortSpec; 1],
}

impl Inertia {
    fn new(mass_kg: f64) -> Self {
        Self {
            mass_kg,
            ports: [mechanical("p", PortRole::EffortIn)],
        }
    }
}

impl Component for Inertia {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, velocity: &f64) -> Result<f64, Cancelled> {
        Ok(0.5 * self.mass_kg * velocity * velocity)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        *new = old + dt_s * inputs[0] / self.mass_kg;
        let mean = 0.5 * (old + *new);
        outputs[0] = mean;
        let stored_change_j = self.mass_kg * mean * (*new - old);
        let scale = self.stored_energy_j(old)?
            + self.stored_energy_j(new)?
            + (dt_s * inputs[0] * mean).abs();
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            roundoff(scale),
        )))
    }
}

impl PortTangent for Inertia {
    fn jvp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        d_outputs[0] = 0.5 * dt_s / self.mass_kg * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        bar_inputs[0] = 0.5 * dt_s / self.mass_kg * bar_outputs[0];
        Ok(())
    }
}

/// Spring: receives flow, returns the effort of the mean displacement. With a
/// cubic term the effort is the discrete gradient of the potential, so the
/// stored change equals the port work for any step.
struct Spring {
    stiffness_n_m: f64,
    cubic_n_m3: f64,
    ports: [PortSpec; 1],
}

impl Spring {
    fn new(stiffness_n_m: f64, cubic_n_m3: f64) -> Self {
        Self {
            stiffness_n_m,
            cubic_n_m3,
            ports: [mechanical("p", PortRole::FlowIn)],
        }
    }

    fn potential(&self, q: f64) -> f64 {
        0.5 * self.stiffness_n_m * q * q + 0.25 * self.cubic_n_m3 * q * q * q * q
    }
}

impl Component for Spring {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        if self.cubic_n_m3 == 0.0 {
            Class::ExactLinear
        } else {
            Class::Smooth
        }
    }

    fn stored_energy_j(&self, q: &f64) -> Result<f64, Cancelled> {
        Ok(self.potential(*q))
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let (a, b) = (*old, old + dt_s * inputs[0]);
        *new = b;
        // (V(b) - V(a)) / (b - a), written without the division.
        let effort =
            0.5 * self.stiffness_n_m * (a + b) + 0.25 * self.cubic_n_m3 * (a + b) * (a * a + b * b);
        outputs[0] = effort;
        let stored_change_j = effort * (b - a);
        let scale = self.potential(a) + self.potential(b) + stored_change_j.abs();
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            roundoff(scale),
        )))
    }
}

impl PortTangent for Spring {
    fn jvp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        let (a, b) = (*old, old + dt_s * inputs[0]);
        let slope =
            0.5 * self.stiffness_n_m + 0.25 * self.cubic_n_m3 * (a * a + 2.0 * a * b + 3.0 * b * b);
        d_outputs[0] = slope * dt_s * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

/// Viscous damper: receives flow, returns `c f`; the loss leaves as heat.
struct Damper {
    damping_n_s_m: f64,
    ports: [PortSpec; 2],
}

impl Damper {
    fn new(damping_n_s_m: f64) -> Self {
        Self {
            damping_n_s_m,
            ports: [
                mechanical("p", PortRole::FlowIn),
                heat("heat", PortRole::EffortIn),
            ],
        }
    }
}

impl Component for Damper {
    type State = ();
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let (flow, temperature_k) = (inputs[0], inputs[1]);
        let loss_w = self.damping_n_s_m * flow * flow;
        outputs[0] = self.damping_n_s_m * flow;
        outputs[1] = -loss_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: dt_s * loss_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * loss_w / temperature_k,
            energy_tolerance_j: roundoff(dt_s * loss_w),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

/// Lumped thermal mass: receives heat rate, returns its mean temperature.
struct ThermalMass {
    capacity_j_k: f64,
    ports: [PortSpec; 1],
}

impl ThermalMass {
    fn new(capacity_j_k: f64) -> Self {
        Self {
            capacity_j_k,
            ports: [heat("heat", PortRole::FlowIn)],
        }
    }
}

impl Component for ThermalMass {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, temperature_k: &f64) -> Result<f64, Cancelled> {
        Ok(self.capacity_j_k * temperature_k)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let heat_j = dt_s * inputs[0];
        *new = old + heat_j / self.capacity_j_k;
        outputs[0] = 0.5 * (old + *new);
        let scale = self.capacity_j_k * (old.abs() + new.abs());
        Ok(Proposal::Step(Receipt::lossless(heat_j, roundoff(scale))))
    }
}

/// A force that depends on time, with the work it has done as its account.
struct ForceSource {
    amplitude_n: f64,
    omega_rad_s: f64,
    ports: [PortSpec; 1],
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct SourceState {
    time_s: f64,
    energy_j: f64,
}

impl ForceSource {
    fn new(amplitude_n: f64, omega_rad_s: f64) -> Self {
        Self {
            amplitude_n,
            omega_rad_s,
            ports: [mechanical("p", PortRole::FlowIn)],
        }
    }

    fn force_n(&self, time_s: f64) -> f64 {
        self.amplitude_n * (self.omega_rad_s * time_s).cos()
    }
}

impl Component for ForceSource {
    type State = SourceState;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &SourceState) -> Result<f64, Cancelled> {
        Ok(state.energy_j)
    }

    fn propose(
        &self,
        old: &SourceState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SourceState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let force_n = self.force_n(old.time_s + 0.5 * dt_s);
        outputs[0] = force_n;
        let change_j = dt_s * force_n * inputs[0];
        *new = SourceState {
            time_s: old.time_s + dt_s,
            energy_j: old.energy_j + change_j,
        };
        let scale = old.energy_j.abs() + new.energy_j.abs() + change_j.abs();
        Ok(Proposal::Step(Receipt::lossless(change_j, roundoff(scale))))
    }
}

fn step_all(system: &mut System, dt_s: f64, steps: u32, gate: &CancelGate) -> u32 {
    with_cx(gate, |cx| {
        let mut iterations = 0;
        for _ in 0..steps {
            iterations += system.step(dt_s, cx).expect("step").iterations;
        }
        iterations
    })
}

const MASS: f64 = 0.25;
const STIFFNESS: f64 = 400.0;
const DAMPING: f64 = 1.5;

#[test]
fn a_bonded_mass_and_spring_follow_the_monolithic_midpoint_rule() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), -0.3)
        .unwrap();
    // The spring's displacement is minus the mass position: its flow is minus
    // the mass velocity.
    let spring = builder
        .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
        .unwrap();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    let mut system = builder.seal().unwrap();
    assert_eq!(system.plan_shape(), (1, 1));
    system.set_config(SolveConfig {
        relative_tolerance: 1.0e-14,
        ..SolveConfig::default()
    });

    let (dt, steps) = (1.0e-3, 5_000_u32);
    let iterations = step_all(&mut system, dt, steps, &gate);
    // Linear block with exact tangents: one Newton iteration per step, plus
    // an occasional second one when roundoff sits on the 1e-14 tolerance.
    assert!(
        iterations >= steps && iterations <= steps + steps / 50,
        "{iterations}"
    );

    // Monolithic implicit midpoint for m x'' = -k x.
    let (mut x, mut v) = (2.0e-3_f64, -0.3_f64);
    let half = 0.5 * dt / MASS;
    for _ in 0..steps {
        let mean = (v - half * STIFFNESS * x) / (1.0 + 0.5 * dt * half * STIFFNESS);
        x += dt * mean;
        v = 2.0 * mean - v;
    }
    let velocity = *system.state::<f64>(mass).unwrap();
    let position = -*system.state::<f64>(spring).unwrap();
    assert!((velocity - v).abs() <= 1.0e-10 * 0.3, "{velocity} vs {v}");
    assert!(
        (position - x).abs() <= 1.0e-10 * 2.0e-3,
        "{position} vs {x}"
    );

    // A lossless graph: stored energy is constant and the ledger closes.
    let ledger = system.ledger();
    let start_j = 0.5 * MASS * 0.09 + 0.5 * STIFFNESS * 4.0e-6;
    assert!(ledger.stored_change_j().abs() <= 1.0e-11 * start_j);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * start_j);
    assert!(ledger.connection_defect_j().abs() <= 1.0e-11 * start_j);
    assert_eq!(ledger.entropy_generated_j_per_k(), 0.0);
    assert_eq!(ledger.steps, u64::from(steps));
    assert!((system.time_s() - 5.0).abs() <= 1.0e-12);
}

struct DampedGraph {
    system: System,
    mass: NodeId,
    spring: NodeId,
    damper: NodeId,
    hot: NodeId,
    cold: NodeId,
}

const FRACTION_HOT: f64 = 0.6;
const FRACTION_COLD: f64 = 0.3;

fn damped_graph() -> DampedGraph {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), -0.3)
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
        .unwrap();
    let damper = builder.add("damper", Damper::new(DAMPING), ()).unwrap();
    let source = builder
        .add(
            "source",
            ForceSource::new(0.8, 180.0),
            SourceState {
                time_s: 0.0,
                energy_j: 0.0,
            },
        )
        .unwrap();
    let hot = builder.add("hot", ThermalMass::new(2.0), 320.0).unwrap();
    let cold = builder.add("cold", ThermalMass::new(5.0), 290.0).unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let members = [
        port(&builder, spring, "p"),
        port(&builder, damper, "p"),
        port(&builder, source, "p"),
    ];
    builder
        .common_flow(port(&builder, mass, "p"), &members)
        .unwrap();
    builder
        .heat_route(
            port(&builder, damper, "heat"),
            &[
                (port(&builder, hot, "heat"), FRACTION_HOT),
                (port(&builder, cold, "heat"), FRACTION_COLD),
            ],
            Some(heat_sink("ambient", 293.15).unwrap()),
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let damper_heat = system.port(damper, "heat").unwrap();
    system.guess(damper_heat, &[300.0]);
    DampedGraph {
        system,
        mass,
        spring,
        damper,
        hot,
        cold,
    }
}

#[test]
fn a_driven_damped_oscillator_heats_its_sinks_and_the_ledger_closes() {
    let gate = CancelGate::new_clock_free();
    let mut graph = damped_graph();
    // Everything is one loop: the damper reads the sinks' temperatures.
    assert_eq!(graph.system.plan_shape(), (1, 1));
    let (dt, steps) = (2.0e-4, 10_000_u32);
    step_all(&mut graph.system, dt, steps, &gate);

    // Monolithic midpoint rule for m x'' = -k x - c x' + F(t).
    let (mut x, mut v, mut loss_j, mut work_j) = (2.0e-3_f64, -0.3_f64, 0.0_f64, 0.0_f64);
    let half = 0.5 * dt / MASS;
    for step in 0..steps {
        let force = 0.8 * (180.0 * (f64::from(step) + 0.5) * dt).cos();
        let mean = (v + half * (force - STIFFNESS * x))
            / (1.0 + half * DAMPING + 0.5 * dt * half * STIFFNESS);
        x += dt * mean;
        v = 2.0 * mean - v;
        loss_j += dt * DAMPING * mean * mean;
        work_j += dt * force * mean;
    }
    let system = &graph.system;
    let velocity = *system.state::<f64>(graph.mass).unwrap();
    let position = -*system.state::<f64>(graph.spring).unwrap();
    assert!((velocity - v).abs() <= 1.0e-9 * 0.3, "{velocity} vs {v}");
    assert!((position - x).abs() <= 1.0e-9 * 2.0e-3, "{position} vs {x}");

    // Every joule the damper lost arrived somewhere, by the declared split.
    let ledger = system.ledger();
    let lost = ledger.nodes[graph.damper.index()].dissipated_j.value();
    assert!(
        (lost - loss_j).abs() <= 1.0e-9 * loss_j,
        "{lost} vs {loss_j}"
    );
    let gained = |node: NodeId| ledger.nodes[node.index()].stored_change_j.value();
    assert!((gained(graph.hot) - FRACTION_HOT * lost).abs() <= 1.0e-12 * lost);
    assert!((gained(graph.cold) - FRACTION_COLD * lost).abs() <= 1.0e-12 * lost);
    assert_eq!(ledger.sinks[0].name, "ambient");
    assert!((ledger.sink_heat_j() - 0.1 * lost).abs() <= 1.0e-12 * lost);
    let hot_k = *system.state::<f64>(graph.hot).unwrap();
    assert!((hot_k - 320.0 - FRACTION_HOT * lost / 2.0).abs() <= 1.0e-10);

    // First law of the whole window. The source's account is the work it
    // delivered, so the system is closed apart from the ambient sink.
    let scale = work_j.abs() + lost;
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * scale);
    assert!(ledger.connection_defect_j().abs() <= 1.0e-9 * scale);
    assert!((ledger.stored_change_j() + ledger.sink_heat_j()).abs() <= 1.0e-9 * scale);
    assert_eq!(ledger.numerical_loss_j(), 0.0);

    // Second law: the damper generates entropy; the route, fed the
    // entropy-equivalent temperature of its receivers, generates none.
    let generated = ledger.entropy_generated_j_per_k();
    assert!(generated > 0.0);
    let route = ledger.connections[1].entropy_generated_j_per_k.value();
    assert!(
        route.abs() <= 1.0e-9 * generated,
        "route {route} of {generated}"
    );
}

/// A whole oscillator as one component, as an owner would offer it.
struct Oscillator {
    ports: [PortSpec; 2],
}

impl Oscillator {
    fn new() -> Self {
        Self {
            ports: [
                mechanical("drive", PortRole::EffortIn),
                heat("heat", PortRole::EffortIn),
            ],
        }
    }
}

impl Component for Oscillator {
    type State = (f64, f64);
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &(f64, f64)) -> Result<f64, Cancelled> {
        Ok(0.5 * MASS * state.1 * state.1 + 0.5 * STIFFNESS * state.0 * state.0)
    }

    fn propose(
        &self,
        old: &(f64, f64),
        dt_s: f64,
        inputs: &[f64],
        new: &mut (f64, f64),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let half = 0.5 * dt_s / MASS;
        let mean = (old.1 + half * (inputs[0] - STIFFNESS * old.0))
            / (1.0 + half * DAMPING + 0.5 * dt_s * half * STIFFNESS);
        *new = (old.0 + dt_s * mean, 2.0 * mean - old.1);
        let loss_w = DAMPING * mean * mean;
        outputs[0] = mean;
        outputs[1] = -loss_w;
        let stored_change_j =
            MASS * mean * (new.1 - old.1) + STIFFNESS * (old.0 + 0.5 * dt_s * mean) * (dt_s * mean);
        let scale = self.stored_energy_j(old)?
            + self.stored_energy_j(new)?
            + (dt_s * inputs[0] * mean).abs()
            + dt_s * loss_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: dt_s * loss_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * loss_w / inputs[1],
            energy_tolerance_j: roundoff(scale),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

#[test]
fn a_one_way_graph_is_swept_once_and_matches_a_hand_staggered_run_bitwise() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let oscillator = builder
        .add("oscillator", Oscillator::new(), (2.0e-3, -0.3))
        .unwrap();
    let hot = builder.add("hot", ThermalMass::new(2.0), 320.0).unwrap();
    let cold = builder.add("cold", ThermalMass::new(5.0), 290.0).unwrap();
    let drive = builder.port(oscillator, "drive").unwrap();
    builder
        .boundary(
            drive,
            Box::new(|t0, dt, value| value[0] = 0.8 * (180.0 * (t0 + 0.5 * dt)).cos()),
        )
        .unwrap();
    let sinks = [
        (builder.port(hot, "heat").unwrap(), 0.75),
        (builder.port(cold, "heat").unwrap(), 0.25),
    ];
    builder
        .heat_route(
            builder.port(oscillator, "heat").unwrap(),
            &sinks,
            None,
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    // No loop anywhere: the oscillator first, then each thermal mass.
    assert_eq!(system.plan_shape(), (3, 0));
    let (dt, steps) = (2.0e-4, 4_000_u32);
    assert_eq!(step_all(&mut system, dt, steps, &gate), 0);

    // The same components stepped by hand in the staggered order.
    let (owner, hot_mass, cold_mass) = (
        Oscillator::new(),
        ThermalMass::new(2.0),
        ThermalMass::new(5.0),
    );
    let (mut state, mut hot_k, mut cold_k) = ((2.0e-3, -0.3), 320.0, 290.0);
    let mut time = fs_system::ledger::Sum::default();
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            let force = 0.8 * (180.0 * (time.value() + 0.5 * dt)).cos();
            let (mut next, mut out) = (state, [0.0; 2]);
            owner
                .propose(&state, dt, &[force, 300.0], &mut next, &mut out, cx)
                .unwrap();
            state = next;
            let (mut next_k, mut mean) = (0.0, [0.0]);
            hot_mass
                .propose(&hot_k, dt, &[-0.75 * out[1]], &mut next_k, &mut mean, cx)
                .unwrap();
            hot_k = next_k;
            cold_mass
                .propose(&cold_k, dt, &[-0.25 * out[1]], &mut next_k, &mut mean, cx)
                .unwrap();
            cold_k = next_k;
            time.add(dt);
        }
    });
    let got = *system.state::<(f64, f64)>(oscillator).unwrap();
    assert_eq!(got.0.to_bits(), state.0.to_bits());
    assert_eq!(got.1.to_bits(), state.1.to_bits());
    assert_eq!(system.state::<f64>(hot).unwrap().to_bits(), hot_k.to_bits());
    assert_eq!(
        system.state::<f64>(cold).unwrap().to_bits(),
        cold_k.to_bits()
    );

    // The boundary did the work; the ledger says where it went.
    let ledger = system.ledger();
    let work = ledger.boundaries[0].energy_in_j.value();
    assert_eq!(ledger.boundaries[0].port, "oscillator.drive");
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * work.abs().max(1.0e-3));
    // A route multiplies by fractions, so its members cancel to roundoff.
    let routed = ledger.nodes[oscillator.index()].dissipated_j.value();
    assert!(ledger.connections[0].energy_defect_j.value().abs() <= 1.0e-14 * routed);
}

#[test]
fn a_nonlinear_spring_without_a_tangent_converges_and_conserves_energy() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), 0.9)
        .unwrap();
    // `add` hides the spring's tangent, so the driver differences its port map.
    let spring = builder
        .add("spring", Spring::new(STIFFNESS, 6.0e7), 0.0)
        .unwrap();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    let mut system = builder.seal().unwrap();
    let steps = 3_000_u32;
    let iterations = step_all(&mut system, 2.0e-4, steps, &gate);
    assert!(
        iterations > steps,
        "a nonlinear loop needs more than one iteration per step"
    );
    assert!(
        iterations < 6 * steps,
        "Newton should converge in a few iterations: {iterations}"
    );
    let start_j = 0.5 * MASS * 0.81;
    let ledger = system.ledger();
    assert!(ledger.stored_change_j().abs() <= 1.0e-9 * start_j);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * start_j);
    // The cubic term matters at this amplitude.
    let reach = system.state::<f64>(spring).unwrap().abs();
    assert!(reach < 0.9 * (MASS * 0.81 / STIFFNESS).sqrt() || reach < 0.02);
}

/// Refuses in `admit` while the shared flag is set.
struct Gatekeeper {
    refuse: Rc<Cell<bool>>,
}

#[derive(Debug, PartialEq)]
struct Refused(&'static str);

impl Component for Gatekeeper {
    type State = u32;
    type Error = Refused;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &u32) -> Result<f64, Refused> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &u32,
        _dt_s: f64,
        _inputs: &[f64],
        new: &mut u32,
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Refused> {
        *new = old + 1;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }

    fn admit(&self, _old: &u32, _new: &u32, _receipt: &Receipt) -> Result<(), Refused> {
        if self.refuse.get() {
            Err(Refused("not yet"))
        } else {
            Ok(())
        }
    }
}

fn gated_system(refuse: &Rc<Cell<bool>>) -> (System, NodeId, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), -0.3)
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
        .unwrap();
    // Added last: it refuses after every other component has proposed.
    let keeper = builder
        .add(
            "keeper",
            Gatekeeper {
                refuse: Rc::clone(refuse),
            },
            0,
        )
        .unwrap();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    (builder.seal().unwrap(), mass, spring, keeper)
}

#[test]
fn a_late_refusal_publishes_nothing_and_a_retry_is_bitwise_a_first_attempt() {
    let gate = CancelGate::new_clock_free();
    let refuse = Rc::new(Cell::new(true));
    let (mut system, mass, spring, keeper) = gated_system(&refuse);
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("refused");
    let SystemError::Component { node, phase, error } = error else {
        panic!("expected the keeper's refusal");
    };
    assert_eq!((node.as_str(), phase), ("keeper", Phase::Admit));
    // The owner's error arrives unchanged.
    assert_eq!(error.downcast_ref::<Refused>(), Some(&Refused("not yet")));
    assert_eq!(
        system.state::<f64>(mass).unwrap().to_bits(),
        (-0.3_f64).to_bits()
    );
    assert_eq!(
        system.state::<f64>(spring).unwrap().to_bits(),
        (-2.0e-3_f64).to_bits()
    );
    assert_eq!(*system.state::<u32>(keeper).unwrap(), 0);
    assert_eq!(
        (system.steps(), system.time_s(), system.ledger().steps),
        (0, 0.0, 0)
    );

    refuse.set(false);
    with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect("admitted");
    let allowed = Rc::new(Cell::new(false));
    let (mut fresh, fresh_mass, fresh_spring, _) = gated_system(&allowed);
    with_cx(&gate, |cx| fresh.step(1.0e-3, cx)).expect("admitted");
    assert_eq!(
        system.state::<f64>(mass).unwrap().to_bits(),
        fresh.state::<f64>(fresh_mass).unwrap().to_bits()
    );
    assert_eq!(
        system.state::<f64>(spring).unwrap().to_bits(),
        fresh.state::<f64>(fresh_spring).unwrap().to_bits()
    );
    assert_eq!(system.ledger(), fresh.ledger());
}

/// Requests cancellation from inside its own step, then carries on.
struct Canceller;

impl Component for Canceller {
    type State = u32;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &u32) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &u32,
        _dt_s: f64,
        _inputs: &[f64],
        new: &mut u32,
        _outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.cancel_gate().request();
        *new = old + 1;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn cancellation_during_a_step_publishes_nothing() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), -0.3)
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
        .unwrap();
    let canceller = builder.add("canceller", Canceller, 0).unwrap();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    let mut system = builder.seal().unwrap();
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("cancelled");
    assert!(matches!(error, SystemError::Cancelled), "{error}");
    assert_eq!(
        system.state::<f64>(mass).unwrap().to_bits(),
        (-0.3_f64).to_bits()
    );
    assert_eq!(*system.state::<u32>(canceller).unwrap(), 0);
    assert_eq!((system.steps(), system.ledger().steps), (0, 0));
    // A gate that is already requested stops the step before any work.
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("cancelled");
    assert!(matches!(error, SystemError::Cancelled));
}

/// Has an event after `remaining` seconds, announced in one of two ways.
struct Timer {
    announce: bool,
}

impl Component for Timer {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _state: &f64) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn next_event_s(&self, remaining: &f64) -> Option<f64> {
        (self.announce && *remaining > 0.0).then_some(*remaining)
    }

    fn propose(
        &self,
        remaining: &f64,
        dt_s: f64,
        _inputs: &[f64],
        new: &mut f64,
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        if *remaining > 0.0 && dt_s > *remaining {
            return Ok(Proposal::Shorten { dt_s: *remaining });
        }
        *new = remaining - dt_s;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn a_component_bounds_the_step_before_or_during_the_trial() {
    let gate = CancelGate::new_clock_free();
    for announce in [false, true] {
        let mut builder = SystemBuilder::new();
        let mass = builder
            .add_with_tangent("mass", Inertia::new(MASS), -0.3)
            .unwrap();
        let spring = builder
            .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
            .unwrap();
        let timer = builder.add("timer", Timer { announce }, 2.5e-4).unwrap();
        let (a, b) = (
            builder.port(mass, "p").unwrap(),
            builder.port(spring, "p").unwrap(),
        );
        builder.bond(a, b).unwrap();
        let mut system = builder.seal().unwrap();
        let report = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect("step");
        assert_eq!(report.dt_s.to_bits(), 2.5e-4_f64.to_bits());
        // Announced events cost no retry; a shortening during the trial does.
        assert_eq!(report.shortenings, u32::from(!announce));
        assert_eq!(*system.state::<f64>(timer).unwrap(), 0.0);
        // The whole system took the short step, not only the timer.
        let report = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect("step");
        assert_eq!(report.dt_s.to_bits(), 1.0e-3_f64.to_bits());
        assert!((system.time_s() - 1.25e-3).abs() <= 1.0e-18);
    }
}

/// A set-valued element that offers no tangent.
struct Latch {
    ports: [PortSpec; 1],
}

impl Component for Latch {
    type State = ();
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        outputs[0] = inputs[0].signum();
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn a_non_smooth_component_in_a_loop_is_refused_by_name() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), -0.3)
        .unwrap();
    let latch = builder
        .add(
            "latch",
            Latch {
                ports: [mechanical("p", PortRole::FlowIn)],
            },
            (),
        )
        .unwrap();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(latch, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    let mut system = builder.seal().unwrap();
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("unsupported");
    let SystemError::Unsupported { node, .. } = error else {
        panic!("expected a refusal that names the component, got {error}");
    };
    assert_eq!(node, "latch");
    assert_eq!(system.steps(), 0);
}

#[test]
fn a_temperature_inside_a_loop_needs_a_starting_value() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let oscillator = builder
        .add("oscillator", Oscillator::new(), (2.0e-3, -0.3))
        .unwrap();
    let hot = builder.add("hot", ThermalMass::new(2.0), 320.0).unwrap();
    let drive = builder.port(oscillator, "drive").unwrap();
    builder
        .boundary(drive, Box::new(|_, _, value| value[0] = 0.0))
        .unwrap();
    builder
        .heat_route(
            builder.port(oscillator, "heat").unwrap(),
            &[(builder.port(hot, "heat").unwrap(), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("no guess");
    let SystemError::GuessRequired { port } = error else {
        panic!("expected a request for a starting temperature, got {error}");
    };
    assert_eq!(port, "oscillator.heat");
    let heat_port = system.port(oscillator, "heat").unwrap();
    system.guess(heat_port, &[300.0]);
    with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect("step");
    // With one receiver the source sees that receiver's mean temperature.
    let seen = system.last_input(heat_port)[0];
    let mean = system.last_output(system.port(hot, "heat").unwrap())[0];
    assert!((seen - mean).abs() <= 1.0e-9, "{seen} vs {mean}");
}

#[test]
fn descriptions_that_cannot_run_are_refused_at_admission() {
    let two = || {
        let mut builder = SystemBuilder::new();
        let mass = builder.add("mass", Inertia::new(MASS), 0.0).unwrap();
        let spring = builder
            .add("spring", Spring::new(STIFFNESS, 0.0), 0.0)
            .unwrap();
        (builder, mass, spring)
    };

    let (mut builder, _, _) = two();
    assert!(matches!(
        builder.add("mass", Inertia::new(1.0), 0.0),
        Err(AdmissionError::DuplicateNode { .. })
    ));

    let (builder, mass, _) = two();
    assert!(matches!(
        builder.port(mass, "missing"),
        Err(AdmissionError::UnknownPort { .. })
    ));

    // A port that exchanges power may not be left open.
    let (builder, _, _) = two();
    let Err(AdmissionError::Unconnected { port }) = builder.seal() else {
        panic!("an open power port must be refused");
    };
    assert_eq!(port, "mass.p");

    // Two ports that both receive effort cannot be bonded.
    let mut builder = SystemBuilder::new();
    let first = builder.add("first", Inertia::new(1.0), 0.0).unwrap();
    let second = builder.add("second", Inertia::new(1.0), 0.0).unwrap();
    let (a, b) = (
        builder.port(first, "p").unwrap(),
        builder.port(second, "p").unwrap(),
    );
    assert!(matches!(
        builder.bond(a, b),
        Err(AdmissionError::Role { .. })
    ));

    // A mechanical port cannot be bonded to a heat port.
    let mut builder = SystemBuilder::new();
    let mass = builder.add("mass", Inertia::new(1.0), 0.0).unwrap();
    let body = builder.add("body", ThermalMass::new(1.0), 300.0).unwrap();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(body, "heat").unwrap(),
    );
    assert!(matches!(
        builder.bond(a, b),
        Err(AdmissionError::Mismatch { .. })
    ));

    // A port joins one connection only.
    let (mut builder, mass, spring) = two();
    let (a, b) = (
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    assert!(matches!(
        builder.bond(a, b),
        Err(AdmissionError::AlreadyConnected { .. })
    ));

    // Heat that is given to nobody needs a named sink.
    let route = |fractions: [f64; 2], sink: bool| {
        let mut builder = SystemBuilder::new();
        let damper = builder.add("damper", Damper::new(1.0), ()).unwrap();
        let hot = builder.add("hot", ThermalMass::new(1.0), 300.0).unwrap();
        let cold = builder.add("cold", ThermalMass::new(1.0), 300.0).unwrap();
        let sinks = [
            (builder.port(hot, "heat").unwrap(), fractions[0]),
            (builder.port(cold, "heat").unwrap(), fractions[1]),
        ];
        builder.heat_route(
            builder.port(damper, "heat").unwrap(),
            &sinks,
            sink.then(|| heat_sink("ambient", 293.15).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
    };
    assert!(matches!(
        route([0.6, 0.3], false),
        Err(AdmissionError::SinkRequired { .. })
    ));
    assert!(route([0.6, 0.3], true).is_ok());
    assert!(route([0.6, 0.4], false).is_ok());
    assert!(matches!(
        route([0.7, 0.4], true),
        Err(AdmissionError::Fractions { .. })
    ));
    assert!(matches!(
        route([0.7, -0.1], true),
        Err(AdmissionError::Fractions { .. })
    ));
}
