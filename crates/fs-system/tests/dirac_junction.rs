//! G0/G2/G4/G5: admitted Dirac laws, residual tiers, closed storage dynamics and replay.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{
    Case, Observed, check_component, check_port_tangent, check_stage_residual,
};
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};
use fs_ports::{
    CoordinateBinding, JunctionLaw, JunctionSchema, PortKind, PortOrientation, PortSchema,
    PortTimestamp, PortValueShape, PowerPairing, StableId,
};
use fs_system::{
    AnalysisError, BlockSolveMode, BlockSolverPolicy, MonolithicConfig, NodeId, StateMap,
    StationaryConfig, System, SystemBuilder,
};
use fs_system::elements::{DiracError, DiracJunction, ElementError, Storage};

const KIND: PortKind = PortKind::ElectricalVoltageCurrent;

fn id(value: &str) -> StableId {
    StableId::new(value).unwrap()
}

fn ports(n: usize) -> Vec<PortSchema> {
    (0..n)
        .map(|i| {
            KIND.scalar_seed_schema(
                id(&format!("terminal/{i}")),
                CoordinateBinding::new(
                    id("voltage-current"),
                    id("circuit"),
                    PortOrientation::OutwardFromOwner,
                ),
                PortTimestamp::new(id("circuit-clock"), 0),
            )
            .unwrap()
        })
        .collect()
}

fn cycle(g: f64) -> JunctionSchema {
    JunctionSchema::try_new(
        id("circulator"),
        ports(3),
        JunctionLaw::Dirac {
            rows: 3,
            incidence: vec![0.0, -g, g, g, 0.0, -g, -g, g, 0.0],
        },
    )
    .unwrap()
}

fn euclidean_ports(n: usize, shape: PortValueShape) -> Vec<PortSchema> {
    ports(n)
        .iter()
        .map(|port| {
            PortSchema::try_new(
                port.id().clone(),
                KIND,
                port.effort_dimensions(),
                port.flow_dimensions(),
                shape,
                port.coordinates().clone(),
                PowerPairing::EuclideanDot,
                port.timestamp().clone(),
                port.conservation_roles().iter().copied(),
            )
            .unwrap()
        })
        .collect()
}

fn vector_cycle(g: f64, width: usize) -> JunctionSchema {
    JunctionSchema::try_new(
        id("circulator"),
        euclidean_ports(3, PortValueShape::vector(width).unwrap()),
        cycle(g).law().clone(),
    )
    .unwrap()
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let key = StreamKey {
            seed: 17,
            kernel_id: 11,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn admitted_relation_and_exact_tangent_preserve_power_and_duality() {
    let schema = cycle(2.0);
    let junction = DiracJunction::new(schema.clone()).unwrap();
    assert_eq!(junction.schema(), &schema);
    assert_eq!(junction.class(), Class::ExactLinear);
    assert_eq!(junction.unknowns(), 0);
    for (port, declared) in junction.ports().iter().zip(schema.ports()) {
        assert_eq!(port.id(), declared.id());
    }
    with_cx(|cx| {
        for inputs in [[2.0, 3.0, 5.0], [-4.0, 2.0, 7.0], [1.0, 1.0, 1.0]] {
            let case = Case {
                state: &(),
                stale: &(),
                dt_s: 0.125,
                inputs: &inputs,
            };
            let Observed::Stepped {
                outputs, receipt, ..
            } = check_component(&junction, &case, cx).unwrap()
            else {
                panic!("a Dirac law never shortens the step");
            };
            assert_eq!(
                outputs,
                [
                    2.0 * (inputs[2] - inputs[1]),
                    2.0 * (inputs[0] - inputs[2]),
                    2.0 * (inputs[1] - inputs[0])
                ]
            );
            assert_eq!(
                outputs.iter().zip(inputs).map(|(f, e)| f * e).sum::<f64>(),
                0.0
            );
            assert_eq!(receipt.dissipated_j, 0.0);
            assert_eq!(receipt.stored_change_j, 0.0);
            check_port_tangent(&junction, &case, 1.0e-7, 1.0e-8, cx).unwrap();
            check_stage_residual(&junction, &case, &[], 1.0e-8, 1.0e-10, cx).unwrap();
            let mut stage_outputs = [f64::NAN; 3];
            junction
                .residual(&(), 0.125, &inputs, &[], &mut [], &mut stage_outputs, cx)
                .unwrap();
            assert_eq!(&stage_outputs[..], outputs.as_slice());
            assert_eq!(
                junction.accept(&(), 0.125, &inputs, &[], &mut (), cx).unwrap(),
                receipt
            );
            let (direction, covector) = ([0.25, -0.5, 1.25], [-2.0, 3.0, 0.5]);
            let (mut push, mut pull) = ([0.0; 3], [0.0; 3]);
            junction
                .jvp(&(), 0.125, &inputs, &direction, &mut push, cx)
                .unwrap();
            junction
                .vjp(&(), 0.125, &inputs, &covector, &mut pull, cx)
                .unwrap();
            let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
            assert_eq!(dot(&push, &covector), dot(&direction, &pull));
        }
    });
}

fn capacitor_graph(schema: JunctionSchema, initial: &[f64]) -> (System, Vec<NodeId>) {
    let mut builder = SystemBuilder::new();
    let junction = builder
        .add_with_tangent_and_stage_residual("junction", DiracJunction::new(schema).unwrap(), ())
        .unwrap();
    let capacitors: Vec<_> = initial
        .iter()
        .enumerate()
        .map(|(i, &charge)| {
            let node = builder
                .add_with_tangent_and_stage_residual(
                    &format!("capacitor/{i}"),
                    Storage::in_domain(KIND, 1.0).unwrap(),
                    charge,
                )
                .unwrap();
            builder
                .bond(
                    builder.port(node, "p").unwrap(),
                    builder.port(junction, &format!("terminal/{i}")).unwrap(),
                )
                .unwrap();
            node
        })
        .collect();
    (builder.seal().unwrap(), capacitors)
}

/// For equal unit capacitors, q' = -J q rotates about (1,1,1).
/// Midpoint advances its angle by -2 atan(sqrt(3) g dt / 2) per step.
/// This independent rotation oracle checks the direction and all three
/// trajectories as well as charge and energy conservation.
#[test]
fn three_capacitors_exchange_energy_along_the_midpoint_rotation() {
    check_capacitor_rotation(BlockSolverPolicy::InterfaceNewton);
}

#[test]
fn monolithic_capacitors_follow_the_independent_midpoint_rotation() {
    check_capacitor_rotation(BlockSolverPolicy::MonolithicNewtonKrylov(
        MonolithicConfig::default(),
    ));
}

fn check_capacitor_rotation(policy: BlockSolverPolicy) {
    let (g, dt) = (0.75, 0.04);
    let initial = [1.0_f64, -0.5, 2.0];
    let (mut system, capacitors) = capacitor_graph(cycle(g), &initial);
    let monolithic = matches!(&policy, BlockSolverPolicy::MonolithicNewtonKrylov(_));
    system.set_solver_policy(policy);
    if monolithic {
        assert!(
            system
                .execution_plan_for_dt(dt)
                .blocks
                .iter()
                .any(|block| block.solve_mode == BlockSolveMode::MonolithicNewtonKrylov)
        );
    }
    let root_three = 3.0_f64.sqrt();
    let mean = initial.iter().sum::<f64>() / 3.0;
    let cross = [
        (initial[2] - initial[1]) / root_three,
        (initial[0] - initial[2]) / root_three,
        (initial[1] - initial[0]) / root_three,
    ];
    let energy = 0.5 * initial.iter().map(|q| q * q).sum::<f64>();
    let mut maximum_exchange = [0.0_f64; 3];
    with_cx(|cx| {
        for step in 1..=100 {
            assert_eq!(system.step(dt, cx).unwrap().dt_s, dt);
            let angle = -(step as f64) * 2.0 * (root_three * g * dt / 2.0).atan();
            let charges: Vec<_> = capacitors
                .iter()
                .map(|&node| *system.state::<f64>(node).unwrap())
                .collect();
            for i in 0..3 {
                let expected = mean + (initial[i] - mean) * angle.cos() + cross[i] * angle.sin();
                near(charges[i], expected, 2.0e-10);
                maximum_exchange[i] = maximum_exchange[i].max((charges[i] - initial[i]).abs());
            }
            near(charges.iter().sum(), initial.iter().sum(), 2.0e-10);
            near(
                0.5 * charges.iter().map(|q| q * q).sum::<f64>(),
                energy,
                2.0e-10,
            );
            near(system.ledger().energy_residual_j(), 0.0, 2.0e-10);
        }
    });
    assert!(maximum_exchange.iter().all(|&exchange| exchange > 0.5));
    for node in &system.ledger().nodes {
        assert_eq!(node.dissipated_j.value(), 0.0);
        assert_eq!(node.numerical_loss_j.value(), 0.0);
    }
    for connection in &system.ledger().connections {
        near(connection.energy_defect_j.value(), 0.0, 2.0e-10);
    }
}

struct ChargeCoordinates(Vec<NodeId>);

impl StateMap for ChargeCoordinates {
    fn dim(&self) -> usize {
        self.0.len()
    }

    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        for (value, &node) in out.iter_mut().zip(&self.0) {
            *value = *system
                .state::<f64>(node)
                .ok_or(AnalysisError::Invalid { what: "capacitor state" })?;
        }
        Ok(())
    }

    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        for (&node, &value) in self.0.iter().zip(state) {
            system.set_accepted_state(node, value)?;
        }
        Ok(())
    }
}

#[test]
fn stationary_two_port_junction_finds_zero_charge_and_restores_the_trajectory() {
    // The full-rank skew law has only q = 0 as an equilibrium. This is an
    // operating-point test, not a claim that conservative dynamics decay.
    let schema = JunctionSchema::try_new(
        id("two-port"),
        ports(2),
        JunctionLaw::Dirac {
            rows: 2,
            incidence: vec![0.0, -1.0, 1.0, 0.0],
        },
    )
    .unwrap();
    let (mut system, capacitors) = capacitor_graph(schema, &[1.0, -2.0]);
    system.set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(
        MonolithicConfig::default(),
    ));
    let coordinates = ChargeCoordinates(capacitors);
    with_cx(|cx| {
        for _ in 0..5 {
            system.step(0.04, cx).unwrap();
        }
        let before = system.checkpoint().unwrap();
        let before_ledger = system.ledger().clone();
        let before_charge: Vec<_> = coordinates
            .0
            .iter()
            .map(|&node| system.state::<f64>(node).unwrap().to_bits())
            .collect();
        let solution = system
            .solve_stationary(
                &coordinates,
                StationaryConfig {
                    probe_dt_s: 0.1,
                    ..StationaryConfig::default()
                },
                cx,
            )
            .unwrap();
        for &charge in &solution.state {
            near(charge, 0.0, 1.0e-8);
        }
        assert!(solution.residual_norm <= 1.0e-10);
        assert_eq!(system.time_s().to_bits(), before.time_s().to_bits());
        assert_eq!(system.steps(), before.steps());
        assert_eq!(system.ledger(), &before_ledger);
        for (&node, expected) in coordinates.0.iter().zip(before_charge) {
            assert_eq!(system.state::<f64>(node).unwrap().to_bits(), expected);
        }
    });
}

#[test]
fn reuses_schema_admission_and_refuses_other_laws_and_shapes() {
    let malformed = JunctionLaw::Dirac {
        rows: 2,
        incidence: vec![0.0, 1.0, 1.0, 0.0],
    };
    assert!(JunctionSchema::try_new(id("invalid"), ports(2), malformed).is_err());
    for law in [
        JunctionLaw::CommonEffort {
            flow_signs: vec![1; 3],
        },
        JunctionLaw::CommonFlow {
            effort_signs: vec![1; 3],
        },
    ] {
        assert_eq!(
            DiracJunction::new(JunctionSchema::try_new(id("common"), ports(3), law).unwrap())
                .unwrap_err(),
            DiracError::NotDirac
        );
    }
    let shape = PortValueShape::tensor(1, 2).unwrap();
    let tensor =
        JunctionSchema::try_new(id("tensor"), euclidean_ports(3, shape), cycle(1.0).law().clone())
            .unwrap();
    assert_eq!(
        DiracJunction::new(tensor).unwrap_err(),
        DiracError::UnsupportedShape { port: 0, shape }
    );
    let mut mismatched = euclidean_ports(3, PortValueShape::vector(2).unwrap());
    mismatched[2] = euclidean_ports(3, PortValueShape::vector(3).unwrap())[2].clone();
    assert!(
        JunctionSchema::try_new(id("mismatch"), mismatched, cycle(1.0).law().clone()).is_err()
    );
    // The schema is small; no vector allocation is needed to refuse a
    // component whose complete runtime layout would overflow.
    assert_eq!(
        DiracJunction::new(vector_cycle(1.0, usize::MAX)).unwrap_err(),
        DiracError::WidthOverflow {
            ports: 3,
            width: usize::MAX,
        }
    );
    // A zero skew law is a valid open junction, not a singularity to invert.
    let zero = DiracJunction::new(cycle(0.0)).unwrap();
    with_cx(|cx| {
        let mut output = [91.0; 3];
        zero.propose(&(), 1.0, &[2.0, 3.0, 4.0], &mut (), &mut output, cx)
            .unwrap();
        assert_eq!(output, [0.0; 3]);
        check_stage_residual(
            &zero,
            &Case {
                state: &(),
                stale: &(),
                dt_s: 1.0,
                inputs: &[2.0, 3.0, 4.0],
            },
            &[],
            1.0e-8,
            1.0e-10,
            cx,
        )
        .unwrap();
    });
}

#[test]
fn invalid_calls_and_cancellation_refuse_without_panicking() {
    let junction = DiracJunction::new(cycle(2.0)).unwrap();
    with_cx(|cx| {
        let mut output = [91.0; 3];
        for (dt, inputs) in [
            (0.0, &[1.0, 2.0, 3.0][..]),
            (1.0, &[1.0, 2.0][..]),
            (1.0, &[1.0, f64::NAN, 3.0][..]),
        ] {
            assert!(
                junction
                    .propose(&(), dt, inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [91.0; 3]);
            assert!(
                junction
                    .residual(&(), dt, inputs, &[], &mut [], &mut output, cx)
                    .is_err()
            );
            assert!(junction.accept(&(), dt, inputs, &[], &mut (), cx).is_err());
            assert_eq!(output, [91.0; 3]);
        }
        assert!(
            junction
                .residual(&(), 1.0, &[1.0; 3], &[0.0], &mut [], &mut output, cx)
                .is_err()
        );
        assert!(junction.accept(&(), 1.0, &[1.0; 3], &[0.0], &mut (), cx).is_err());
        assert_eq!(output, [91.0; 3]);
        cx.cancel_gate().request();
        assert_eq!(
            junction.propose(&(), 1.0, &[1.0; 3], &mut (), &mut output, cx),
            Err(DiracError::Cancelled)
        );
        assert_eq!(
            junction.residual(&(), 1.0, &[1.0; 3], &[], &mut [], &mut output, cx),
            Err(DiracError::Cancelled)
        );
        assert_eq!(
            junction.jvp_z(&(), 1.0, &[1.0; 3], &[], &[], &mut [], &mut output, cx),
            Err(DiracError::Cancelled)
        );
        assert_eq!(
            junction.vjp_z(&(), 1.0, &[1.0; 3], &[], &[], &[1.0; 3], &mut [], cx),
            Err(DiracError::Cancelled)
        );
        assert_eq!(
            junction.accept(&(), 1.0, &[1.0; 3], &[], &mut (), cx),
            Err(DiracError::Cancelled)
        );
        assert_eq!(output, [91.0; 3]);
    });
}

#[test]
fn vector_relation_keeps_coordinates_separate_and_preserves_power_and_duality() {
    let junction = DiracJunction::new(vector_cycle(2.0, 2)).unwrap();
    assert!(junction.ports().iter().all(|port| port.width() == 2));
    with_cx(|cx| {
        // Port-major pairs: J acts separately on (2, 5, 7) and (3, -4, 11).
        let inputs = [2.0, 3.0, 5.0, -4.0, 7.0, 11.0];
        let case = Case {
            state: &(),
            stale: &(),
            dt_s: 0.125,
            inputs: &inputs,
        };
        let Observed::Stepped {
            outputs, receipt, ..
        } = check_component(&junction, &case, cx).unwrap()
        else {
            panic!("a Dirac law never shortens the step");
        };
        assert_eq!(outputs, [4.0, 30.0, -10.0, -16.0, 6.0, -14.0]);
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        assert_eq!(dot(&inputs, &outputs), 0.0);
        assert_eq!(receipt.stored_change_j, 0.0);
        assert_eq!(receipt.dissipated_j, 0.0);
        check_port_tangent(&junction, &case, 1.0e-7, 1.0e-8, cx).unwrap();
        check_stage_residual(&junction, &case, &[], 1.0e-8, 1.0e-10, cx).unwrap();
        let mut stage_outputs = [f64::NAN; 6];
        junction
            .residual(&(), 0.125, &inputs, &[], &mut [], &mut stage_outputs, cx)
            .unwrap();
        assert_eq!(&stage_outputs[..], outputs.as_slice());
        assert_eq!(
            junction.accept(&(), 0.125, &inputs, &[], &mut (), cx).unwrap(),
            receipt
        );

        let direction = [0.25, 1.5, -0.5, 0.25, 1.25, -2.0];
        let covector = [-2.0, 0.5, 3.0, 1.25, 0.5, -0.75];
        let (mut push, mut pull) = ([0.0; 6], [0.0; 6]);
        junction
            .jvp(&(), 0.125, &inputs, &direction, &mut push, cx)
            .unwrap();
        junction
            .vjp(&(), 0.125, &inputs, &covector, &mut pull, cx)
            .unwrap();
        assert_eq!(push, [3.5, -4.5, -2.0, 7.0, -1.5, -2.5]);
        assert_eq!(dot(&push, &covector), dot(&direction, &pull));
    });
}

/// A test-only vector port around two unchanged native unit capacitors.
/// Every coordinate delegates its storage law, step and tangent to Storage.
struct CapacitorPair {
    scalar: Storage,
    ports: [PortSpec; 1],
}

impl CapacitorPair {
    fn new() -> Self {
        Self {
            scalar: Storage::in_domain(KIND, 1.0).unwrap(),
            ports: [PortSpec::power(
                id("p"),
                KIND,
                PortRole::FlowIn,
                NonZeroUsize::new(2).unwrap(),
            )
            .unwrap()],
        }
    }
}

impl Component for CapacitorPair {
    type State = [f64; 2];
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, ElementError> {
        Ok(self.scalar.stored_energy_j(&state[0])? + self.scalar.stored_energy_j(&state[1])?)
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let mut receipt = Receipt::lossless(0.0, 0.0);
        for i in 0..2 {
            match self.scalar.propose(
                &old[i],
                dt_s,
                &inputs[i..i + 1],
                &mut new[i],
                &mut outputs[i..i + 1],
                cx,
            )? {
                Proposal::Step(step) => {
                    receipt.stored_change_j += step.stored_change_j;
                    receipt.energy_tolerance_j += step.energy_tolerance_j;
                }
                shorter => return Ok(shorter),
            }
        }
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for CapacitorPair {
    fn jvp(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        for i in 0..2 {
            self.scalar.jvp(
                &old[i],
                dt_s,
                &inputs[i..i + 1],
                &d_inputs[i..i + 1],
                &mut d_outputs[i..i + 1],
                cx,
            )?;
        }
        Ok(())
    }

    fn vjp(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // Each native capacitor's scalar tangent is its own transpose.
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

#[test]
fn vector_capacitors_follow_independent_rotations_and_replay_from_checkpoint() {
    let (g, dt) = (0.75, 0.04);
    let initial = [[1.0_f64, -0.25], [-0.5, 1.5], [2.0, 0.75]];
    let mut builder = SystemBuilder::new();
    let junction = builder
        .add_with_tangent(
            "junction",
            DiracJunction::new(vector_cycle(g, 2)).unwrap(),
            (),
        )
        .unwrap();
    let capacitors: Vec<_> = initial
        .iter()
        .enumerate()
        .map(|(i, &charge)| {
            let node = builder
                .add_with_tangent(&format!("capacitor/{i}"), CapacitorPair::new(), charge)
                .unwrap();
            builder
                .bond(
                    builder.port(node, "p").unwrap(),
                    builder.port(junction, &format!("terminal/{i}")).unwrap(),
                )
                .unwrap();
            node
        })
        .collect();
    let mut system = builder.seal().unwrap();
    let root_three = 3.0_f64.sqrt();
    let energy = 0.5 * initial.iter().flatten().map(|q| q * q).sum::<f64>();
    let observe = |system: &System| {
        (
            system.time_s().to_bits(),
            system.steps(),
            capacitors
                .iter()
                .map(|&node| (*system.state::<[f64; 2]>(node).unwrap()).map(f64::to_bits))
                .collect::<Vec<_>>(),
            system.ledger().clone(),
        )
    };
    with_cx(|cx| {
        let mut checkpoint = None;
        for step in 1..=100 {
            assert_eq!(system.step(dt, cx).unwrap().dt_s, dt);
            let angle = -(step as f64) * 2.0 * (root_three * g * dt / 2.0).atan();
            let charges: Vec<_> = capacitors
                .iter()
                .map(|&node| *system.state::<[f64; 2]>(node).unwrap())
                .collect();
            for coordinate in 0..2 {
                let q = initial.map(|pair| pair[coordinate]);
                let mean = q.iter().sum::<f64>() / 3.0;
                let cross = [
                    (q[2] - q[1]) / root_three,
                    (q[0] - q[2]) / root_three,
                    (q[1] - q[0]) / root_three,
                ];
                for i in 0..3 {
                    let expected =
                        mean + (q[i] - mean) * angle.cos() + cross[i] * angle.sin();
                    near(charges[i][coordinate], expected, 2.0e-10);
                }
                near(
                    charges.iter().map(|pair| pair[coordinate]).sum(),
                    q.iter().sum(),
                    2.0e-10,
                );
            }
            near(
                0.5 * charges.iter().flatten().map(|q| q * q).sum::<f64>(),
                energy,
                2.0e-10,
            );
            near(system.ledger().energy_residual_j(), 0.0, 2.0e-10);
            if step == 37 {
                checkpoint = Some(system.checkpoint().unwrap());
            }
        }
        let expected = observe(&system);
        system.restore(checkpoint.as_ref().unwrap()).unwrap();
        for _ in 37..100 {
            system.step(dt, cx).unwrap();
        }
        assert_eq!(observe(&system), expected);
    });
}

#[test]
fn vector_shape_and_incidence_remain_bound_to_checkpoint_identity() {
    let identity = |junction: &DiracJunction| {
        let mut words = Vec::new();
        junction.identity_words(&mut words);
        words
    };
    let scalar = DiracJunction::new(cycle(2.0)).unwrap();
    let vector_one = DiracJunction::new(vector_cycle(2.0, 1)).unwrap();
    let vector_two = DiracJunction::new(vector_cycle(2.0, 2)).unwrap();
    let changed = DiracJunction::new(vector_cycle(2.0_f64.next_up(), 2)).unwrap();
    assert_ne!(identity(&scalar), identity(&vector_one));
    assert_ne!(identity(&vector_one), identity(&vector_two));
    assert_ne!(identity(&vector_two), identity(&changed));
    assert_eq!(
        identity(&vector_two),
        identity(&DiracJunction::new(vector_cycle(2.0, 2)).unwrap())
    );
    let mut state_words = Vec::new();
    assert!(vector_two.encode_state(&(), &mut state_words));
    assert!(state_words.is_empty());
    assert_eq!(vector_two.decode_state(&state_words), Some(()));
    assert_eq!(vector_two.decode_state(&[0]), None);
}

#[test]
fn residual_acceptance_checks_flows_even_when_all_power_terms_are_zero() {
    let junction = DiracJunction::new(
        JunctionSchema::try_new(
            id("star"),
            ports(3),
            JunctionLaw::Dirac {
                rows: 3,
                incidence: vec![0.0, 1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
            },
        )
        .unwrap(),
    )
    .unwrap();
    with_cx(|cx| {
        // The first flow is MAX + MAX, but its own effort is zero.
        // Every term of effort_i * J_ij * effort_j is exactly zero.
        let inputs = [0.0, f64::MAX, f64::MAX];
        let mut output = [91.0; 3];
        let error = junction
            .propose(&(), 1.0, &inputs, &mut (), &mut output, cx)
            .unwrap_err();
        assert_eq!(
            junction
                .residual(&(), 1.0, &inputs, &[], &mut [], &mut output, cx)
                .unwrap_err(),
            error
        );
        assert_eq!(
            junction.accept(&(), 1.0, &inputs, &[], &mut (), cx).unwrap_err(),
            error
        );
        assert_eq!(output, [91.0; 3]);
    });
}

#[test]
fn vector_calls_refuse_wrong_extents_nonfinite_values_and_cancellation() {
    let junction = DiracJunction::new(vector_cycle(2.0, 2)).unwrap();
    with_cx(|cx| {
        let mut output = [91.0; 6];
        for (dt, inputs) in [
            (0.0, &[1.0; 6][..]),
            (1.0, &[1.0; 5][..]),
            (1.0, &[1.0, 2.0, 3.0, 4.0, 5.0, f64::NAN][..]),
        ] {
            assert!(
                junction
                    .propose(&(), dt, inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [91.0; 6]);
        }
        let mut short = [91.0; 5];
        assert!(
            junction
                .propose(&(), 1.0, &[1.0; 6], &mut (), &mut short, cx)
                .is_err()
        );
        assert_eq!(short, [91.0; 5]);
        assert!(
            junction
                .jvp(&(), 1.0, &[1.0; 6], &[1.0; 5], &mut output, cx)
                .is_err()
        );
        assert!(
            junction
                .vjp(&(), 1.0, &[1.0; 6], &[1.0; 6], &mut short, cx)
                .is_err()
        );
        assert_eq!(output, [91.0; 6]);
        assert_eq!(short, [91.0; 5]);
        assert!(
            junction
                .propose(&(), 1.0, &[f64::MAX; 6], &mut (), &mut output, cx)
                .is_err()
        );
        cx.cancel_gate().request();
        assert_eq!(
            junction.propose(&(), 1.0, &[1.0; 6], &mut (), &mut output, cx),
            Err(DiracError::Cancelled)
        );
    });
}
