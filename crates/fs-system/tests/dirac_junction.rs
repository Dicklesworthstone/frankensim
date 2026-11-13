//! Admitted N-way Dirac laws, their tangent, and a closed capacitor network.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortTangent};
use fs_ports::{
    CoordinateBinding, JunctionLaw, JunctionSchema, PortKind, PortOrientation, PortSchema,
    PortTimestamp, PortValueShape, PowerPairing, StableId,
};
use fs_system::SystemBuilder;
use fs_system::elements::{DiracError, DiracJunction, Storage};

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

/// For equal unit capacitors, q' = -J q rotates about (1,1,1).
/// Midpoint advances its angle by -2 atan(sqrt(3) g dt / 2) per step.
/// This independent rotation oracle checks the direction and all three
/// trajectories as well as charge and energy conservation.
#[test]
fn three_capacitors_exchange_energy_along_the_midpoint_rotation() {
    let (g, dt) = (0.75, 0.04);
    let initial = [1.0_f64, -0.5, 2.0];
    let mut builder = SystemBuilder::new();
    let junction = builder
        .add_with_tangent("junction", DiracJunction::new(cycle(g)).unwrap(), ())
        .unwrap();
    let capacitors: Vec<_> = initial
        .iter()
        .enumerate()
        .map(|(i, &charge)| {
            let node = builder
                .add_with_tangent(
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
    let mut system = builder.seal().unwrap();
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
    let vector_ports = ports(3)
        .iter()
        .map(|port| {
            PortSchema::try_new(
                port.id().clone(),
                KIND,
                port.effort_dimensions(),
                port.flow_dimensions(),
                PortValueShape::vector(2).unwrap(),
                port.coordinates().clone(),
                PowerPairing::EuclideanDot,
                port.timestamp().clone(),
                port.conservation_roles().iter().copied(),
            )
            .unwrap()
        })
        .collect();
    let vector =
        JunctionSchema::try_new(id("vector"), vector_ports, cycle(1.0).law().clone()).unwrap();
    assert_eq!(
        DiracJunction::new(vector).unwrap_err(),
        DiracError::NonScalar { port: 0 }
    );
    // A zero skew law is a valid open junction, not a singularity to invert.
    let zero = DiracJunction::new(cycle(0.0)).unwrap();
    with_cx(|cx| {
        let mut output = [91.0; 3];
        zero.propose(&(), 1.0, &[2.0, 3.0, 4.0], &mut (), &mut output, cx)
            .unwrap();
        assert_eq!(output, [0.0; 3]);
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
        }
        cx.cancel_gate().request();
        assert_eq!(
            junction.propose(&(), 1.0, &[1.0; 3], &mut (), &mut output, cx),
            Err(DiracError::Cancelled)
        );
        assert_eq!(output, [91.0; 3]);
    });
}
