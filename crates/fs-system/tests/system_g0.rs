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
    let (mut iterations, mut reach) = (0_u32, 0.0_f64);
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            iterations += system.step(2.0e-4, cx).expect("step").iterations;
            reach = reach.max(system.state::<f64>(spring).unwrap().abs());
        }
    });
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
    // The turning point follows from energy alone, independently of the
    // integrator: E = k q^2 / 2 + k3 q^4 / 4. The largest displacement seen
    // over several periods must reach it and may not pass it.
    let cubic = 6.0e7;
    let turning =
        ((-STIFFNESS + (STIFFNESS * STIFFNESS + 4.0 * cubic * start_j).sqrt()) / cubic).sqrt();
    assert!(
        reach <= turning * (1.0 + 1.0e-9),
        "{reach} passes {turning}"
    );
    assert!(reach >= 0.995 * turning, "{reach} falls short of {turning}");
    // The cubic term matters at this amplitude: a linear spring would reach
    // more than twice as far.
    assert!(turning < 0.5 * (2.0 * start_j / STIFFNESS).sqrt());
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
        // The whole system took the short step, not only the timer: the
        // oscillator is where a twin without a timer is after 0.25 ms.
        let mut twin = SystemBuilder::new();
        let twin_mass = twin
            .add_with_tangent("mass", Inertia::new(MASS), -0.3)
            .unwrap();
        let twin_spring = twin
            .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
            .unwrap();
        let (a, b) = (
            twin.port(twin_mass, "p").unwrap(),
            twin.port(twin_spring, "p").unwrap(),
        );
        twin.bond(a, b).unwrap();
        let mut twin = twin.seal().unwrap();
        with_cx(&gate, |cx| twin.step(2.5e-4, cx)).expect("twin step");
        assert_eq!(
            system.state::<f64>(mass).unwrap().to_bits(),
            twin.state::<f64>(twin_mass).unwrap().to_bits()
        );
        assert_eq!(
            system.state::<f64>(spring).unwrap().to_bits(),
            twin.state::<f64>(twin_spring).unwrap().to_bits()
        );
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

// ---------------------------------------------------------------------------
// Regressions from the independent review of 2026-11-10.
// ---------------------------------------------------------------------------

/// Mass, spring, optional damper and a constant force on one common-flow
/// junction. The damper's heat goes to an ambient sink.
fn loaded_oscillator(damping: Option<f64>, position_m: f64) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), 0.0)
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -position_m)
        .unwrap();
    let source = builder
        .add(
            "load",
            ForceSource::new(0.7, 0.0),
            SourceState {
                time_s: 0.0,
                energy_j: 0.0,
            },
        )
        .unwrap();
    let mut members = vec![
        builder.port(spring, "p").unwrap(),
        builder.port(source, "p").unwrap(),
    ];
    if let Some(damping) = damping {
        let damper = builder.add("damper", Damper::new(damping), ()).unwrap();
        members.push(builder.port(damper, "p").unwrap());
        builder
            .heat_route(
                builder.port(damper, "heat").unwrap(),
                &[],
                Some(heat_sink("ambient", 293.15).unwrap()),
                RouteTemperature::Fixed(293.15),
            )
            .unwrap();
    }
    builder
        .common_flow(builder.port(mass, "p").unwrap(), &members)
        .unwrap();
    (builder.seal().unwrap(), mass, spring)
}

#[test]
fn a_loaded_oscillator_settles_to_rest_under_the_default_tolerances() {
    // Near rest the force on the mass is a small difference of large efforts.
    // The convergence test must be scaled by those efforts, or the residual
    // can never get under it and the run stops with NotConverged.
    let gate = CancelGate::new_clock_free();
    let (mut system, mass, spring) = loaded_oscillator(Some(DAMPING), 0.0);
    // The envelope decays as exp(-3 t): ten seconds leave 1e-13 of it.
    step_all(&mut system, 1.0e-3, 10_000, &gate);
    let rest_m = 0.7 / STIFFNESS;
    let position = -*system.state::<f64>(spring).unwrap();
    assert!(
        (position - rest_m).abs() <= 1.0e-9 * rest_m,
        "{position} vs {rest_m}"
    );
    assert!(system.state::<f64>(mass).unwrap().abs() <= 1.0e-9);
    // The work of the load ended as spring energy and ambient heat.
    let ledger = system.ledger();
    assert!((ledger.sink_heat_j() - 0.5 * 0.7 * rest_m).abs() <= 1.0e-9 * 0.7 * rest_m);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-12);

    // An undamped oscillator started exactly at its loaded equilibrium stays
    // there, at a small and at a large step.
    for dt in [1.0e-3, 5.0e-2] {
        let (mut system, mass, spring) = loaded_oscillator(None, rest_m);
        step_all(&mut system, dt, 2_000, &gate);
        assert!(system.state::<f64>(mass).unwrap().abs() <= 1.0e-12);
        assert!((-*system.state::<f64>(spring).unwrap() - rest_m).abs() <= 1.0e-12 * rest_m);
    }
}

#[test]
fn a_refused_connection_leaves_the_builder_as_it_was() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Inertia::new(MASS), -0.3)
        .unwrap();
    let first = builder
        .add_with_tangent("first", Spring::new(STIFFNESS, 0.0), 0.0)
        .unwrap();
    let second = builder
        .add_with_tangent("second", Spring::new(STIFFNESS, 0.0), 0.0)
        .unwrap();
    let other = builder
        .add_with_tangent("other", Inertia::new(MASS), 0.1)
        .unwrap();
    let port = |builder: &SystemBuilder, node| builder.port(node, "p").unwrap();
    let (m, p1, p2, o) = (
        port(&builder, mass),
        port(&builder, first),
        port(&builder, second),
        port(&builder, other),
    );
    builder.bond(o, p2).unwrap();
    // `second` is taken, so the junction is refused as a whole...
    assert!(matches!(
        builder.common_flow(m, &[p1, p2]),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    // ...and so is one that names a port twice.
    assert!(matches!(
        builder.common_flow(m, &[p1, p1]),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    // Neither refusal claimed `mass.p` or `first.p`: the valid call succeeds,
    // the system seals and it runs.
    builder.common_flow(m, &[p1]).unwrap();
    let mut system = builder.seal().unwrap();
    step_all(&mut system, 1.0e-3, 10, &gate);
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-15);

    // A refused heat route leaves no sink behind in the ledger.
    let mut builder = SystemBuilder::new();
    let damper = builder.add("damper", Damper::new(1.0), ()).unwrap();
    let hot = builder.add("hot", ThermalMass::new(1.0), 300.0).unwrap();
    let (source, sink) = (
        builder.port(damper, "heat").unwrap(),
        builder.port(hot, "heat").unwrap(),
    );
    assert!(matches!(
        builder.heat_route(
            source,
            &[(sink, 0.4), (sink, 0.4)],
            Some(heat_sink("ambient", 293.15).unwrap()),
            RouteTemperature::Fixed(300.0),
        ),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    builder
        .heat_route(
            source,
            &[(sink, 0.4)],
            Some(heat_sink("ambient", 293.15).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mechanical_port = builder.port(damper, "p").unwrap();
    builder
        .boundary(mechanical_port, Box::new(|_, _, flow| flow[0] = 0.5))
        .unwrap();
    let system = builder.seal().unwrap();
    assert_eq!(system.ledger().sinks.len(), 1);
}

#[test]
fn heat_ports_cannot_share_one_heat_rate_among_several_members() {
    let mut builder = SystemBuilder::new();
    let oscillator = builder
        .add("oscillator", Oscillator::new(), (2.0e-3, -0.3))
        .unwrap();
    let hot = builder.add("hot", ThermalMass::new(2.0), 320.0).unwrap();
    let cold = builder.add("cold", ThermalMass::new(5.0), 290.0).unwrap();
    let setter = builder.port(oscillator, "heat").unwrap();
    let members = [
        builder.port(hot, "heat").unwrap(),
        builder.port(cold, "heat").unwrap(),
    ];
    let Err(AdmissionError::HeatCommonFlow { port }) = builder.common_flow(setter, &members) else {
        panic!("a heat rate handed to two members would be counted twice");
    };
    assert_eq!(port, "oscillator.heat");
    // One member is an ordinary bond and is admitted.
    builder.common_flow(setter, &members[..1]).unwrap();
}

/// Emits a large constant signal and the cosine of the signal it receives.
struct CosineLoop {
    large: f64,
    ports: [PortSpec; 3],
}

impl CosineLoop {
    fn new(large: f64) -> Self {
        Self {
            large,
            ports: [
                PortSpec::signal_in(id("x"), NonZeroUsize::MIN),
                PortSpec::signal_out(id("large"), NonZeroUsize::MIN),
                PortSpec::signal_out(id("cos"), NonZeroUsize::MIN),
            ],
        }
    }
}

impl Component for CosineLoop {
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
        _dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        outputs[0] = self.large;
        outputs[1] = inputs[0].cos();
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

/// Passes a signal through unchanged.
struct Relay {
    ports: [PortSpec; 2],
}

impl Component for Relay {
    type State = ();
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
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
        outputs[0] = inputs[0];
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn a_signal_loop_is_judged_by_its_own_size_not_by_a_large_neighbour() {
    // x = cos(x) has the root 0.7390851332151607. A second signal of 1e12 in
    // the same block must not loosen the test on x.
    let gate = CancelGate::new_clock_free();
    for large in [1.0, 1.0e12] {
        let mut builder = SystemBuilder::new();
        let looped = builder.add("loop", CosineLoop::new(large), ()).unwrap();
        let relay = builder
            .add(
                "relay",
                Relay {
                    ports: [
                        PortSpec::signal_in(id("in"), NonZeroUsize::MIN),
                        PortSpec::signal_out(id("out"), NonZeroUsize::MIN),
                    ],
                },
                (),
            )
            .unwrap();
        // cos -> x closes the loop; large -> relay -> nothing puts the large
        // signal's consumer next to it.
        builder
            .signal(
                builder.port(looped, "cos").unwrap(),
                builder.port(looped, "x").unwrap(),
            )
            .unwrap();
        builder
            .signal(
                builder.port(looped, "large").unwrap(),
                builder.port(relay, "in").unwrap(),
            )
            .unwrap();
        let mut system = builder.seal().unwrap();
        with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect("step");
        let x = system.last_input(system.port(looped, "x").unwrap())[0];
        assert!(
            (x - 0.739_085_133_215_160_7).abs() <= 1.0e-10,
            "large {large}: x = {x}"
        );
        // The relay saw the large signal exactly.
        let relayed = system.last_output(system.port(relay, "out").unwrap())[0];
        assert_eq!(relayed.to_bits(), large.to_bits());
    }
}

#[test]
fn two_masses_on_one_spring_share_its_force() {
    // Common effort with two members: both masses feel the spring force and
    // the spring shortens by the sum of their displacements.
    let gate = CancelGate::new_clock_free();
    let (m1, m2, k) = (0.25, 0.4, STIFFNESS);
    let (v1, v2, q) = (-0.3, 0.2, 1.5e-3);
    let mut builder = SystemBuilder::new();
    let first = builder
        .add_with_tangent("first", Inertia::new(m1), v1)
        .unwrap();
    let second = builder
        .add_with_tangent("second", Inertia::new(m2), v2)
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(k, 0.0), q)
        .unwrap();
    let members = [
        builder.port(first, "p").unwrap(),
        builder.port(second, "p").unwrap(),
    ];
    builder
        .common_effort(builder.port(spring, "p").unwrap(), &members)
        .unwrap();
    let mut system = builder.seal().unwrap();
    let (dt, steps) = (1.0e-3, 3_000_u32);
    step_all(&mut system, dt, steps, &gate);

    // Monolithic midpoint rule: e = k q_mid, q' = -(v1 + v2), m_i v_i' = e.
    let (mut a, mut b, mut c) = (v1, v2, q);
    for _ in 0..steps {
        let effort =
            k * (c - 0.5 * dt * (a + b)) / (1.0 + 0.25 * k * dt * dt * (1.0 / m1 + 1.0 / m2));
        let (a_mean, b_mean) = (a + 0.5 * dt * effort / m1, b + 0.5 * dt * effort / m2);
        c -= dt * (a_mean + b_mean);
        a = 2.0 * a_mean - a;
        b = 2.0 * b_mean - b;
    }
    let got = |node| *system.state::<f64>(node).unwrap();
    assert!((got(first) - a).abs() <= 1.0e-10, "{} vs {a}", got(first));
    assert!((got(second) - b).abs() <= 1.0e-10, "{} vs {b}", got(second));
    assert!((got(spring) - c).abs() <= 1.0e-11, "{} vs {c}", got(spring));
    let start_j = 0.5 * m1 * v1 * v1 + 0.5 * m2 * v2 * v2 + 0.5 * k * q * q;
    assert!(system.ledger().stored_change_j().abs() <= 1.0e-11 * start_j);
}

/// Returns `gain * input + offset` on its one mechanical port; it keeps no
/// accounts, because the steps it is used in are refused before the audit.
struct Affine {
    gain: f64,
    offset: f64,
    ports: [PortSpec; 1],
}

impl Component for Affine {
    type State = ();
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
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
        outputs[0] = self.gain * inputs[0] + self.offset;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

impl PortTangent for Affine {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        d_outputs[0] = self.gain * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        bar_inputs[0] = self.gain * bar_outputs[0];
        Ok(())
    }
}

#[test]
fn solver_failures_are_typed_and_name_what_failed() {
    let gate = CancelGate::new_clock_free();
    let oscillator = || {
        let mut builder = SystemBuilder::new();
        let mass = builder
            .add_with_tangent("mass", Inertia::new(MASS), -0.3)
            .unwrap();
        let spring = builder
            .add_with_tangent("spring", Spring::new(STIFFNESS, 0.0), -2.0e-3)
            .unwrap();
        let (a, b) = (
            builder.port(mass, "p").unwrap(),
            builder.port(spring, "p").unwrap(),
        );
        builder.bond(a, b).unwrap();
        (builder.seal().unwrap(), mass, spring)
    };

    // No iteration allowed: the loop cannot be solved and says where.
    let (mut system, _, _) = oscillator();
    system.set_config(SolveConfig {
        max_iterations: 0,
        ..SolveConfig::default()
    });
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("not converged");
    let SystemError::NotConverged {
        nodes,
        iterations,
        worst,
        residual,
    } = error
    else {
        panic!("expected NotConverged, got {error}");
    };
    assert_eq!((nodes.len(), iterations), (2, 0));
    assert!(worst == "mass.p" || worst == "spring.p", "{worst}");
    assert!(residual != 0.0);
    assert_eq!(system.steps(), 0);

    // A tolerance so loose that a wrong first guess is accepted leaves the
    // two sides of the bond with different energies; the interface check
    // refuses what the solver let through.
    let (mut system, mass, spring) = oscillator();
    system.set_config(SolveConfig {
        relative_tolerance: 10.0,
        ..SolveConfig::default()
    });
    for node in [mass, spring] {
        let port = system.port(node, "p").expect("port");
        assert!(system.guess(port, &[1.0]));
    }
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("interface");
    assert!(
        matches!(error, SystemError::Interface { connection: 0, .. }),
        "{error}"
    );
    assert_eq!(system.steps(), 0);

    // A conductance of 2 against a resistance of -1/2: the interface
    // Jacobian [[1, 1/2], [2, 1]] is singular. Exact tangents make it so
    // exactly; a differenced tangent would only make it nearly so.
    let mut builder = SystemBuilder::new();
    let conductance = builder
        .add_with_tangent(
            "conductance",
            Affine {
                gain: 2.0,
                offset: 0.0,
                ports: [mechanical("p", PortRole::EffortIn)],
            },
            (),
        )
        .unwrap();
    let resistance = builder
        .add_with_tangent(
            "resistance",
            Affine {
                gain: -0.5,
                offset: 1.0,
                ports: [mechanical("p", PortRole::FlowIn)],
            },
            (),
        )
        .unwrap();
    let (a, b) = (
        builder.port(conductance, "p").unwrap(),
        builder.port(resistance, "p").unwrap(),
    );
    builder.bond(a, b).unwrap();
    let mut system = builder.seal().unwrap();
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("singular");
    let SystemError::Singular { nodes } = error else {
        panic!("expected Singular, got {error}");
    };
    assert_eq!(nodes, ["conductance", "resistance"]);
}

/// Answers the step negotiation with fixed values, valid or not.
// The shared `_s` suffix is the unit, seconds.
#[allow(clippy::struct_field_names)]
struct Negotiator {
    limit_s: f64,
    event_s: Option<f64>,
    shorten_to_s: Option<f64>,
}

impl Component for Negotiator {
    type State = ();
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn limit_dt_s(&self, _state: &(), _proposed_s: f64) -> Result<f64, Cancelled> {
        Ok(self.limit_s)
    }

    fn next_event_s(&self, _state: &()) -> Option<f64> {
        self.event_s
    }

    fn propose(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        _new: &mut (),
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        Ok(match self.shorten_to_s {
            Some(dt_s) => Proposal::Shorten { dt_s },
            None => Proposal::Step(Receipt::lossless(0.0, 0.0)),
        })
    }
}

#[test]
fn a_bad_answer_in_the_step_negotiation_is_refused_and_names_the_component() {
    let gate = CancelGate::new_clock_free();
    let attempt = |limit_s: f64, event_s: Option<f64>, shorten_to_s: Option<f64>| {
        let mut builder = SystemBuilder::new();
        builder
            .add(
                "negotiator",
                Negotiator {
                    limit_s,
                    event_s,
                    shorten_to_s,
                },
                (),
            )
            .unwrap();
        let mut system = builder.seal().unwrap();
        with_cx(&gate, |cx| system.step(1.0e-3, cx))
    };
    // A valid limit and event are taken.
    assert_eq!(
        attempt(5.0e-4, Some(2.5e-4), None).unwrap().dt_s.to_bits(),
        2.5e-4_f64.to_bits()
    );
    // A limit that is not a number is refused even beside a valid event, and
    // so is an event that is not a number beside a valid limit.
    for (limit, event) in [
        (f64::NAN, Some(2.5e-4)),
        (5.0e-4, Some(f64::NAN)),
        (0.0, None),
        (2.0e-3, None),
        (5.0e-4, Some(0.0)),
    ] {
        let error = attempt(limit, event, None).expect_err("refused");
        let SystemError::Negotiation { node, .. } = error else {
            panic!("limit {limit}, event {event:?}: expected Negotiation, got {error}");
        };
        assert_eq!(node, "negotiator");
    }
    // A request to shorten to a step that is not shorter names its author.
    for shorter in [1.0e-3, 2.0e-3, 0.0, f64::NAN] {
        let error = attempt(1.0e-3, None, Some(shorter)).expect_err("refused");
        let SystemError::Contract { node, .. } = error else {
            panic!("shorten to {shorter}: expected Contract, got {error}");
        };
        assert_eq!(node, "negotiator");
    }
}

/// Gives off a fixed heat rate from a finite store and claims to generate no
/// entropy doing so.
struct Heater {
    watts: f64,
    ports: [PortSpec; 1],
}

impl Component for Heater {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, store_j: &f64) -> Result<f64, Cancelled> {
        Ok(*store_j)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        *new = old - dt_s * self.watts;
        outputs[0] = -self.watts;
        Ok(Proposal::Step(Receipt::lossless(
            -dt_s * self.watts,
            roundoff(old.abs() + new.abs()),
        )))
    }
}

#[test]
fn a_route_that_would_destroy_entropy_is_refused() {
    // The heater says its heat leaves at 100 K and generates no entropy, and
    // the route delivers that heat to a body at 300 K: entropy Q/100 leaves
    // the heater and only Q/300 arrives. Nothing accounts for the difference.
    let gate = CancelGate::new_clock_free();
    let build = |source_k: f64| {
        let mut builder = SystemBuilder::new();
        let heater = builder
            .add(
                "heater",
                Heater {
                    watts: 5.0,
                    ports: [heat("heat", PortRole::EffortIn)],
                },
                100.0,
            )
            .unwrap();
        let body = builder.add("body", ThermalMass::new(2.0), 300.0).unwrap();
        builder
            .heat_route(
                builder.port(heater, "heat").unwrap(),
                &[(builder.port(body, "heat").unwrap(), 1.0)],
                None,
                RouteTemperature::Fixed(source_k),
            )
            .unwrap();
        (builder.seal().unwrap(), body)
    };
    let (mut system, body) = build(100.0);
    let error = with_cx(&gate, |cx| system.step(1.0e-2, cx)).expect_err("second law");
    let SystemError::SecondLaw { generated_j_per_k } = error else {
        panic!("expected SecondLaw, got {error}");
    };
    assert!(generated_j_per_k < 0.0);
    assert_eq!(
        system.state::<f64>(body).unwrap().to_bits(),
        300.0_f64.to_bits()
    );
    // Heat that leaves at 1000 K and arrives at 300 K generates entropy and
    // is admitted.
    let (mut system, body) = build(1000.0);
    with_cx(&gate, |cx| system.step(1.0e-2, cx)).expect("step");
    assert!(*system.state::<f64>(body).unwrap() > 300.0);
    assert!(system.ledger().entropy_generated_j_per_k() > 0.0);
}

/// Requests cancellation and then notices it itself, returning its own error.
struct SelfCancelling;

impl Component for SelfCancelling {
    type State = ();
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        _new: &mut (),
        _outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.cancel_gate().request();
        cx.checkpoint()?;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn cancellation_noticed_by_a_component_is_reported_as_cancellation() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    builder.add("stopper", SelfCancelling, ()).unwrap();
    let mut system = builder.seal().unwrap();
    let error = with_cx(&gate, |cx| system.step(1.0e-3, cx)).expect_err("cancelled");
    assert!(matches!(error, SystemError::Cancelled), "{error}");
    assert_eq!(system.steps(), 0);
}
