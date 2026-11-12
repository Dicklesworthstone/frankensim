//! G0/G3: stream inventories are independently audited before publication.
//! Energy closure alone cannot detect a component creating gas mass.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};
use fs_system::elements::{GasContents, GasVolume, Restriction};
use fs_system::{SystemBuilder, SystemError};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}, tolerance {tolerance}"
    );
}

fn air() -> GasVolume {
    GasVolume::new(287.0, 1.4, 0).unwrap()
}

#[test]
fn two_vessels_and_a_restriction_conserve_audited_mass_through_flow_reversal() {
    let mut builder = SystemBuilder::new();
    let pipe = builder
        .add_with_tangent("pipe", Restriction::new(2.0e-8, 0).unwrap(), ())
        .unwrap();
    let mut nodes = Vec::new();
    let mut initial_mass = Vec::new();
    for (side, (name, volume)) in [("a", 0.01), ("b", 0.03)].into_iter().enumerate() {
        let gas = air();
        let contents = gas.contents(2.0e5, 300.0, volume, &[]).unwrap();
        initial_mass.push(contents.mass_kg);
        let node = builder.add_with_tangent(name, gas, contents).unwrap();
        builder
            .bond(
                builder.port(node, "stream").unwrap(),
                builder.port(pipe, name).unwrap(),
            )
            .unwrap();
        builder
            .boundary(
                builder.port(node, "heat").unwrap(),
                Box::new(move |t0, dt, values| {
                    let time = t0 + 0.5 * dt;
                    values[0] = if side == 0 && time < 1.0 {
                        800.0
                    } else if side == 1 && (1.0..2.0).contains(&time) {
                        4_000.0
                    } else {
                        0.0
                    };
                }),
            )
            .unwrap();
        nodes.push(node);
    }
    let outlet = builder.port(pipe, "a").unwrap();
    let mut system = builder.seal().unwrap();
    let total_initial: f64 = initial_mass.iter().sum();
    let (mut forward, mut reverse) = (false, false);
    with_cx(|cx| {
        for _ in 0..300 {
            system.step(0.01, cx).unwrap();
            let rate = system.last_output(outlet)[0];
            forward |= rate > 1.0e-8;
            reverse |= rate < -1.0e-8;
            let total: f64 = nodes
                .iter()
                .map(|&node| system.state::<GasContents>(node).unwrap().mass_kg)
                .sum();
            close(total, total_initial, 2.0e-14);
            close(system.ledger().mass_residual_kg().unwrap(), 0.0, 2.0e-14);
        }
    });
    assert!(forward && reverse, "both donors must be exercised");
    let ledger = system.ledger();
    for (&node, &initial) in nodes.iter().zip(&initial_mass) {
        let held = system.state::<GasContents>(node).unwrap().mass_kg;
        let account = ledger.nodes[node.index()].mass.as_ref().unwrap();
        close(account.stored_change_kg.value(), held - initial, 1.0e-15);
        close(account.port_mass_in_kg.value(), held - initial, 2.0e-14);
        close(account.balance_defect_kg.value(), 0.0, 2.0e-14);
        assert_eq!(account.audited_steps, system.steps());
    }
    let transport = ledger.nodes[pipe.index()].mass.as_ref().unwrap();
    assert_eq!(transport.stored_change_kg.value(), 0.0);
    assert_eq!(transport.port_mass_in_kg.value(), 0.0);
    assert_eq!(transport.audited_steps, system.steps());
    for connection in &ledger.connections {
        close(connection.mass_defect_kg.value(), 0.0, 2.0e-14);
    }
    for boundary in &ledger.boundaries {
        assert_eq!(boundary.mass_in_kg.value(), 0.0);
    }
}

#[test]
fn filling_books_boundary_mass_and_checkpoint_replay_restores_every_mass_account() {
    let gas = air();
    let start = gas.contents(1.0e5, 300.0, 0.02, &[]).unwrap();
    let rate = 5.0e-3;
    let enthalpy = 1.4 * gas.heat_capacity_j_kg_k() * 300.0;
    let mut builder = SystemBuilder::new();
    let vessel = builder
        .add_with_tangent("vessel", gas, start.clone())
        .unwrap();
    builder
        .boundary(
            builder.port(vessel, "stream").unwrap(),
            Box::new(move |_, _, values| values.copy_from_slice(&[rate, rate * enthalpy])),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(vessel, "heat").unwrap(),
            Box::new(|_, _, values| values.fill(0.0)),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..100 {
            system.step(0.01, cx).unwrap();
        }
        let checkpoint = system.checkpoint().unwrap();
        let ledger_at_checkpoint = system.ledger().clone();
        for _ in 0..100 {
            system.step(0.01, cx).unwrap();
        }
        let final_state = system.state::<GasContents>(vessel).unwrap().clone();
        let final_ledger = system.ledger().clone();
        let account = final_ledger.nodes[vessel.index()].mass.as_ref().unwrap();
        let added = rate * 2.0;
        close(final_state.mass_kg, start.mass_kg + added, 2.0e-15);
        close(
            final_ledger.boundaries[0].mass_in_kg.value(),
            added,
            1.0e-17,
        );
        close(account.port_mass_in_kg.value(), added, 1.0e-17);
        close(account.stored_change_kg.value(), added, 2.0e-15);
        close(final_ledger.mass_residual_kg().unwrap(), 0.0, 2.0e-15);
        assert_eq!(account.audited_steps, 200);

        system.restore(&checkpoint).unwrap();
        assert_eq!(system.ledger(), &ledger_at_checkpoint);
        assert_eq!(system.steps(), 100);
        for _ in 0..100 {
            system.step(0.01, cx).unwrap();
        }
        assert_eq!(system.state::<GasContents>(vessel).unwrap(), &final_state);
        assert_eq!(system.ledger(), &final_ledger);
    });
}

/// An energy-neutral inventory whose injected mass leak is invisible to the
/// energy audit. The external switch repairs a refused trial without touching
/// the component's accepted state or the system's private trial buffers.
struct Inventory {
    ports: [PortSpec; 1],
    leak_kg: Rc<Cell<f64>>,
    measure: bool,
    fixed_mass_rate_kg_s: Option<f64>,
}

impl Inventory {
    fn new(leak_kg: Rc<Cell<f64>>, measure: bool) -> Self {
        Self {
            ports: [
                PortSpec::stream(StableId::new("stream").unwrap(), PortRole::FlowIn, 0).unwrap(),
            ],
            leak_kg,
            measure,
            fixed_mass_rate_kg_s: None,
        }
    }

    fn donor(rate_kg_s: f64) -> Self {
        let mut owner = Self::new(Rc::new(Cell::new(0.0)), true);
        owner.ports[0] =
            PortSpec::stream(StableId::new("stream").unwrap(), PortRole::EffortIn, 0).unwrap();
        owner.fixed_mass_rate_kg_s = Some(rate_kg_s);
        owner
    }
}

impl Component for Inventory {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _: &f64) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn stored_mass_kg(&self, mass: &f64) -> Result<Option<f64>, Cancelled> {
        Ok(self.measure.then_some(*mass))
    }

    fn propose(
        &self,
        old: &f64,
        dt: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let rate = if let Some(rate) = self.fixed_mass_rate_kg_s {
            outputs.copy_from_slice(&[rate, 0.0]);
            rate
        } else {
            outputs.copy_from_slice(&[1.0e5, 0.0]);
            inputs[0]
        };
        *new = old + dt * rate + self.leak_kg.get();
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn an_energy_balanced_mass_leak_refuses_without_publishing_any_component() {
    let leak = Rc::new(Cell::new(0.0));
    let mut builder = SystemBuilder::new();
    let good = builder
        .add("good", Inventory::new(Rc::new(Cell::new(0.0)), true), 1.0)
        .unwrap();
    // Put the faulty component last so an earlier passing mass audit cannot
    // accidentally publish that component before the later refusal.
    let bad = builder
        .add("leaky", Inventory::new(Rc::clone(&leak), true), 2.0)
        .unwrap();
    for node in [good, bad] {
        builder
            .boundary(
                builder.port(node, "stream").unwrap(),
                Box::new(|_, _, values| values.copy_from_slice(&[0.25, 0.0])),
            )
            .unwrap();
    }
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.1, cx).unwrap();
        let before = (
            *system.state::<f64>(good).unwrap(),
            *system.state::<f64>(bad).unwrap(),
            system.time_s(),
            system.steps(),
            system.ledger().clone(),
        );
        leak.set(0.01);
        assert!(matches!(
            system.step(0.1, cx),
            Err(SystemError::Contract {
                node,
                violation: ContractViolation::MassBalance { .. },
            }) if node == "leaky"
        ));
        assert_eq!(*system.state::<f64>(good).unwrap(), before.0);
        assert_eq!(*system.state::<f64>(bad).unwrap(), before.1);
        assert_eq!(system.time_s(), before.2);
        assert_eq!(system.steps(), before.3);
        assert_eq!(system.ledger(), &before.4);

        leak.set(0.0);
        system.step(0.1, cx).unwrap();
        close(*system.state::<f64>(good).unwrap(), 1.05, 1.0e-15);
        close(*system.state::<f64>(bad).unwrap(), 2.05, 1.0e-15);
        close(system.ledger().mass_residual_kg().unwrap(), 0.0, 1.0e-15);
        for node in &system.ledger().nodes {
            assert_eq!(node.mass.as_ref().unwrap().audited_steps, 2);
        }
    });
}

#[test]
fn an_unmeasured_stream_inventory_has_no_mass_conservation_verdict() {
    let mut builder = SystemBuilder::new();
    let node = builder
        .add(
            "unknown",
            Inventory::new(Rc::new(Cell::new(0.0)), false),
            1.0,
        )
        .unwrap();
    builder
        .boundary(
            builder.port(node, "stream").unwrap(),
            Box::new(|_, _, values| values.copy_from_slice(&[0.01, 0.0])),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| system.step(0.1, cx).map(|_| ()).unwrap());
    assert!(system.ledger().nodes[node.index()].has_streams);
    assert!(system.ledger().nodes[node.index()].mass.is_none());
    close(
        system.ledger().boundaries[0].mass_in_kg.value(),
        0.001,
        1.0e-18,
    );
    assert_eq!(system.ledger().mass_residual_kg(), None);
}

#[test]
fn a_finite_balanced_mass_transfer_does_not_overflow_the_connection_audit() {
    let mut builder = SystemBuilder::new();
    let donor = builder
        .add("donor", Inventory::donor(-1.0e308), 1.0e308)
        .unwrap();
    let receiver = builder
        .add(
            "receiver",
            Inventory::new(Rc::new(Cell::new(0.0)), true),
            0.0,
        )
        .unwrap();
    builder
        .bond(
            builder.port(donor, "stream").unwrap(),
            builder.port(receiver, "stream").unwrap(),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| system.step(1.0, cx).map(|_| ()).unwrap());
    assert_eq!(*system.state::<f64>(donor).unwrap(), 0.0);
    assert_eq!(*system.state::<f64>(receiver).unwrap(), 1.0e308);
    assert_eq!(system.steps(), 1);
    assert_eq!(system.ledger().connections[0].mass_defect_kg.value(), 0.0);
    assert_eq!(system.ledger().mass_residual_kg(), Some(0.0));
}
