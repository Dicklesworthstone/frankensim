//! G0/G3: independently measured transported inventories, including a species
//! conversion that leaves total mass and energy unchanged.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::conformance::{Case, Failure, check_component};
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt, StreamBasis,
    StreamQuantity,
};
use fs_qty::chemistry::SpeciesId;
use fs_qty::{Dims, Mass};
use fs_system::elements::{GasContents, GasVolume, Restriction};
use fs_system::{AdmissionError, SystemBuilder, SystemError};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 39,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    });
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}, tolerance {tolerance}"
    );
}

#[test]
fn different_mixtures_conserve_each_inventory_when_the_donor_reverses() {
    assert_eq!(mixture_trajectory(false), mixture_trajectory(true));
}

fn species(name: &str) -> StreamQuantity {
    StreamQuantity::SpeciesMass(SpeciesId::new(name).unwrap())
}

fn mixture_trajectory(identified: bool) -> (Vec<u64>, fs_system::Ledger) {
    let gas = || {
        if identified {
            GasVolume::identified(287.0, 1.4, vec![species("N2"), species("O2")]).unwrap()
        } else {
            GasVolume::new(287.0, 1.4, 2).unwrap()
        }
    };
    let mut builder = SystemBuilder::new();
    let pipe = builder
        .add_with_tangent("pipe", Restriction::for_stream(2.0e-8, &gas().ports()[0]).unwrap(), ())
        .unwrap();
    let mut nodes = Vec::new();
    let mut initial = [0.0; 2];
    for (side, (name, volume, fractions)) in [("a", 0.01, [0.8, 0.2]), ("b", 0.03, [0.1, 0.9])]
        .into_iter()
        .enumerate()
    {
        let gas = gas();
        let contents = gas.contents(2.0e5, 300.0, volume, &fractions).unwrap();
        for i in 0..2 {
            initial[i] += contents.carried_kg[i];
        }
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
                Box::new(move |t, dt, values| {
                    let time = t + 0.5 * dt;
                    values[0] = if side == 0 && time < 1.0 {
                        800.0
                    } else if side == 1 && (1.0..2.0).contains(&time) {
                        4000.0
                    } else {
                        0.0
                    };
                }),
            )
            .unwrap();
        nodes.push(node);
    }
    let port = builder.port(pipe, "a").unwrap();
    let mut system = builder.seal().unwrap();
    let (mut forward, mut reverse) = (false, false);
    let mut trajectory = Vec::new();
    with_cx(|cx| {
        for _ in 0..300 {
            system.step(0.01, cx).unwrap();
            let rate = system.last_output(port)[0];
            forward |= rate > 1e-8;
            reverse |= rate < -1e-8;
            for &node in &nodes {
                let state = system.state::<GasContents>(node).unwrap();
                trajectory.extend([state.mass_kg, state.energy_j, state.volume_m3]
                    .into_iter().chain(state.carried_kg.iter().copied()).map(f64::to_bits));
            }
            for i in 0..2 {
                let total: f64 = nodes
                    .iter()
                    .map(|&node| system.state::<GasContents>(node).unwrap().carried_kg[i])
                    .sum();
                close(total, initial[i], 2e-14);
                close(system.ledger().carried_residual(i).unwrap(), 0.0, 2e-14);
            }
        }
    });
    assert!(forward && reverse);
    for node in &system.ledger().nodes {
        for account in &node.carried {
            assert_eq!(account.as_ref().unwrap().audited_steps, 300);
            close(account.as_ref().unwrap().balance_defect.value(), 0.0, 2e-14);
        }
    }
    for connection in &system.ledger().connections {
        for amount in &connection.carried_defect {
            close(amount.value(), 0.0, 2e-14);
        }
    }
    assert_eq!(system.ledger().carried_residual(2), None);
    (trajectory, system.ledger().clone())
}

#[test]
fn stream_connections_refuse_wrong_species_units_materials_and_energy_datums() {
    let volume = GasVolume::identified(287.0, 1.4, vec![species("N2"), species("O2")]).unwrap();
    let original = volume.ports()[0].stream_basis().unwrap();
    let make =
        |medium, reference, quantities| StreamBasis::new(medium, reference, quantities).unwrap();
    let cases = [
        None,
        Some(make(
            original.medium().clone(),
            original.enthalpy_reference().clone(),
            vec![species("O2"), species("N2")],
        )),
        Some(make(
            original.medium().clone(),
            original.enthalpy_reference().clone(),
            vec![species("N2"), species("Ar")],
        )),
        Some(make(
            original.medium().clone(),
            StableId::new("other-enthalpy-datum").unwrap(),
            original.quantities().to_vec(),
        )),
        Some(make(
            StableId::new("other-gas").unwrap(),
            original.enthalpy_reference().clone(),
            original.quantities().to_vec(),
        )),
    ];
    for (case, basis) in cases.into_iter().enumerate() {
        let port = match basis {
            Some(basis) => {
                PortSpec::identified_stream(StableId::new("s").unwrap(), PortRole::FlowIn, basis)
                    .unwrap()
            }
            None => PortSpec::stream(StableId::new("s").unwrap(), PortRole::FlowIn, 2).unwrap(),
        };
        for junction in 0..3 {
            let mut builder = SystemBuilder::new();
            let contents = volume.contents(1e5, 300.0, 0.01, &[0.8, 0.2]).unwrap();
            let gas = builder
                .add_with_tangent("gas", volume.clone(), contents)
                .unwrap();
            let pipe = builder
                .add_with_tangent("pipe", Restriction::for_stream(1e-8, &port).unwrap(), ())
                .unwrap();
            let a = builder.port(gas, "stream").unwrap();
            let b = builder.port(pipe, "a").unwrap();
            let result = match junction {
                0 => builder.bond(a, b),
                1 => builder.common_effort(a, &[b]),
                _ => builder.common_flow(b, &[a]),
            };
            assert!(
                matches!(result, Err(AdmissionError::StreamBasis { .. })),
                "case {case}, junction {junction}: {result:?}"
            );
            // A rejected connection must not claim either endpoint.
            builder
                .boundary(a, Box::new(|_, _, values| values.fill(0.0)))
                .unwrap();
            builder
                .boundary(b, Box::new(|_, _, values| values.fill(0.0)))
                .unwrap();
        }
    }
    // Same named coordinate with different SI inventory units is incompatible.
    let quantity = |dimensions| StreamQuantity::Named {
        id: StableId::new("inventory").unwrap(),
        dimensions,
    };
    let declared = |dimensions| {
        PortSpec::identified_stream(
            StableId::new("s").unwrap(),
            PortRole::FlowIn,
            make(
                original.medium().clone(),
                original.enthalpy_reference().clone(),
                vec![quantity(dimensions)],
            ),
        )
        .unwrap()
    };
    let mass = GasVolume::identified(287.0, 1.4, vec![quantity(Mass::DIMS)]).unwrap();
    let mut builder = SystemBuilder::new();
    let contents = mass.contents(1e5, 300.0, 0.01, &[1.0]).unwrap();
    let gas = builder.add_with_tangent("gas", mass, contents).unwrap();
    let pipe = builder
        .add_with_tangent(
            "pipe",
            Restriction::for_stream(1e-8, &declared(Dims::NONE)).unwrap(),
            (),
        )
        .unwrap();
    assert!(matches!(
        builder.bond(
            builder.port(gas, "stream").unwrap(),
            builder.port(pipe, "a").unwrap()
        ),
        Err(AdmissionError::StreamBasis { .. })
    ));
    assert!(GasVolume::identified(287.0, 1.4, vec![quantity(Dims::NONE)]).is_err());
    // Even with no extra coordinates, changing the gas changes its stream basis.
    let a = GasVolume::identified(287.0, 1.4, vec![]).unwrap();
    for b in [
        GasVolume::identified(300.0, 1.4, vec![]).unwrap(),
        GasVolume::identified(287.0, 1.5, vec![]).unwrap(),
    ] {
        assert_ne!(a.ports()[0].stream_basis(), b.ports()[0].stream_basis());
    }
}

#[test]
fn boundary_inventories_and_their_compensated_accounts_replay_exactly() {
    let gas = GasVolume::new(287.0, 1.4, 2).unwrap();
    let start = gas.contents(1.0e5, 300.0, 0.02, &[0.2, 0.8]).unwrap();
    let rates = [0.003, 0.002];
    let enthalpy = 1.4 * gas.heat_capacity_j_kg_k() * 300.0;
    let mut builder = SystemBuilder::new();
    let vessel = builder
        .add_with_tangent("vessel", gas, start.clone())
        .unwrap();
    builder
        .boundary(
            builder.port(vessel, "stream").unwrap(),
            Box::new(move |_, _, values| {
                values.copy_from_slice(&[0.005, 0.005 * enthalpy, rates[0], rates[1]]);
            }),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(vessel, "heat").unwrap(),
            Box::new(|_, _, values| values.fill(0.0)),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert_eq!(system.ledger().carried_residual(0), None);
    with_cx(|cx| {
        for _ in 0..100 {
            system.step(0.01, cx).unwrap();
        }
        let checkpoint = system.checkpoint().unwrap();
        for _ in 0..100 {
            system.step(0.01, cx).unwrap();
        }
        let state = system.state::<GasContents>(vessel).unwrap().clone();
        let ledger = system.ledger().clone();
        for i in 0..2 {
            close(
                state.carried_kg[i],
                start.carried_kg[i] + 2.0 * rates[i],
                2e-15,
            );
            close(
                ledger.boundaries[0].carried_in[i].value(),
                2.0 * rates[i],
                1e-17,
            );
            close(ledger.carried_residual(i).unwrap(), 0.0, 2e-15);
        }
        system.restore(&checkpoint).unwrap();
        for _ in 0..100 {
            system.step(0.01, cx).unwrap();
        }
        assert_eq!(system.state::<GasContents>(vessel).unwrap(), &state);
        assert_eq!(system.ledger(), &ledger);
    });
}

struct Mixture {
    ports: [PortSpec; 1],
    conversion: Rc<Cell<f64>>,
    measure: Rc<Cell<bool>>,
}

#[test]
fn a_trace_constituent_limits_the_step_by_its_own_inventory() {
    let gas = GasVolume::new(287.0, 1.4, 1).unwrap();
    let start = gas.contents(1.0e5, 300.0, 0.02, &[1.0e-30]).unwrap();
    let rate = -2.0 * start.carried_kg[0];
    let enthalpy = 1.4 * gas.heat_capacity_j_kg_k() * 300.0;
    let mut builder = SystemBuilder::new();
    let vessel = builder
        .add_with_tangent("vessel", gas, start.clone())
        .unwrap();
    builder
        .boundary(
            builder.port(vessel, "stream").unwrap(),
            Box::new(move |_, _, values| {
                values.copy_from_slice(&[rate, rate * enthalpy, rate]);
            }),
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
        let report = system.step(1.0, cx).unwrap();
        assert_eq!(report.dt_s, 0.5);
        assert_eq!(report.shortenings, 1);
        assert_eq!(
            system.state::<GasContents>(vessel).unwrap().carried_kg[0],
            0.0
        );
        close(system.ledger().carried_residual(0).unwrap(), 0.0, 1e-45);
        let account = system.ledger().nodes[vessel.index()].carried[0]
            .as_ref()
            .unwrap();
        assert_eq!(account.port_in.value(), -start.carried_kg[0]);
    });
}

impl Component for Mixture {
    type State = [f64; 2];
    type Error = Cancelled;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &Self::State) -> Result<f64, Cancelled> {
        Ok(0.0)
    }
    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Cancelled> {
        Ok(Some(state[0] + state[1]))
    }
    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        carried: usize,
    ) -> Result<Option<f64>, Cancelled> {
        Ok(self.measure.get().then_some(state[carried]))
    }
    fn propose(
        &self,
        old: &Self::State,
        dt: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        *new = [
            old[0] + dt * inputs[2] + self.conversion.get(),
            old[1] + dt * inputs[3] - self.conversion.get(),
        ];
        outputs.copy_from_slice(&[1e5, 0.0, new[0] / 2.0, new[1] / 2.0]);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn unreported_conversion_is_rejected_atomically_and_unknown_is_not_zero() {
    let conversion = Rc::new(Cell::new(0.0));
    let measure = Rc::new(Cell::new(true));
    let mut builder = SystemBuilder::new();
    let mut nodes = Vec::new();
    for name in ["first", "last"] {
        let node = builder
            .add(
                name,
                Mixture {
                    ports: [PortSpec::stream(
                        StableId::new("stream").unwrap(),
                        PortRole::FlowIn,
                        2,
                    )
                    .unwrap()],
                    conversion: if name == "last" {
                        Rc::clone(&conversion)
                    } else {
                        Rc::new(Cell::new(0.0))
                    },
                    measure: Rc::clone(&measure),
                },
                [1.0, 1.0],
            )
            .unwrap();
        builder
            .boundary(
                builder.port(node, "stream").unwrap(),
                Box::new(|_, _, values| values.fill(0.0)),
            )
            .unwrap();
        nodes.push(node);
    }
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.1, cx).unwrap();
        let ledger = system.ledger().clone();
        conversion.set(0.125);
        assert!(matches!(system.step(0.1, cx), Err(SystemError::Contract {
            node, violation: ContractViolation::CarriedBalance { carried: 0, .. },
        }) if node == "last"));
        assert_eq!(system.steps(), 1);
        assert_eq!(system.ledger(), &ledger);
        for &node in &nodes {
            assert_eq!(system.state::<[f64; 2]>(node).unwrap(), &[1.0, 1.0]);
        }
        let component = Mixture {
            ports: [
                PortSpec::stream(StableId::new("stream").unwrap(), PortRole::FlowIn, 2).unwrap(),
            ],
            conversion: Rc::clone(&conversion),
            measure: Rc::clone(&measure),
        };
        let case = Case {
            state: &[1.0, 1.0],
            stale: &[5.0, 9.0],
            dt_s: 0.1,
            inputs: &[0.0; 4],
        };
        assert!(matches!(
            check_component(&component, &case, cx),
            Err(Failure::Contract(ContractViolation::CarriedBalance {
                carried: 0,
                ..
            }))
        ));
        conversion.set(0.0);
        system.step(0.1, cx).unwrap();
        assert_eq!(system.ledger().carried_residual(0), Some(0.0));
        measure.set(false);
        system.step(0.1, cx).unwrap();
        assert_eq!(system.ledger().carried_residual(0), None);
        measure.set(true);
        system.step(0.1, cx).unwrap();
        assert_eq!(
            system.ledger().carried_residual(0),
            None,
            "later data cannot fill an unaudited step"
        );
        assert_eq!(
            system.ledger().nodes[0].carried[0]
                .as_ref()
                .unwrap()
                .audited_steps,
            3
        );
    });
}
