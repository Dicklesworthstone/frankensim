//! G0/G1: nonmatching field traces retain their physical duality metrics.
//! The trajectory oracle condenses the Galerkin midpoint equations directly;
//! it does not infer correctness from the kernel's conservation ledger.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::{
    Class, Component, DualityMetric, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};
use fs_system::{AdmissionError, SolveConfig, SystemBuilder, SystemError};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
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

/// Independent bilinear form, deliberately not using the port's energy API.
fn pair(matrix: &[f64], a: &[f64], b: &[f64]) -> f64 {
    let mut value = 0.0;
    for i in 0..a.len() {
        for j in 0..b.len() {
            value += a[i] * matrix[i * b.len() + j] * b[j];
        }
    }
    value
}

fn field_port(role: PortRole, matrix: &[f64], width: usize) -> PortSpec {
    let metric = DualityMetric::new(NonZeroUsize::new(width).unwrap(), matrix.to_vec()).unwrap();
    PortSpec::field_power(
        StableId::new("trace").unwrap(),
        PortKind::MechanicalForceVelocity,
        role,
        metric,
    )
    .unwrap()
}

/// A Galerkin inertia or elastic storage with a constant trace mass matrix.
/// Inertia receives effort and returns mean velocity; elasticity receives
/// flow and returns mean stress. Density or stiffness is uniform.
struct FieldStorage {
    ports: [PortSpec; 1],
    matrix: Vec<f64>,
    coefficient: f64,
}

impl FieldStorage {
    fn new(role: PortRole, matrix: &[f64], width: usize, coefficient: f64) -> Self {
        Self {
            ports: [field_port(role, matrix, width)],
            matrix: matrix.to_vec(),
            coefficient,
        }
    }

    fn inertial(&self) -> bool {
        self.ports[0].role() == PortRole::EffortIn
    }
}

impl Component for FieldStorage {
    type State = Vec<f64>;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, Cancelled> {
        Ok(0.5 * self.coefficient * pair(&self.matrix, state, state))
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt: f64,
        input: &[f64],
        new: &mut Vec<f64>,
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let inertial = self.inertial();
        for i in 0..old.len() {
            new[i] = old[i] + dt * input[i] / if inertial { self.coefficient } else { 1.0 };
            output[i] = 0.5 * (old[i] + new[i]) * if inertial { 1.0 } else { self.coefficient };
        }
        let mut change = 0.0;
        for i in 0..old.len() {
            for j in 0..old.len() {
                change += self.coefficient
                    * 0.5
                    * (old[i] + new[i])
                    * self.matrix[i * old.len() + j]
                    * (new[j] - old[j]);
            }
        }
        let scale = self.stored_energy_j(old)? + self.stored_energy_j(new)? + change.abs();
        Ok(Proposal::Step(Receipt::lossless(
            change,
            256.0 * f64::EPSILON * scale,
        )))
    }
}

impl PortTangent for FieldStorage {
    fn jvp(
        &self,
        _old: &Vec<f64>,
        dt: f64,
        _input: &[f64],
        direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        cx.checkpoint()?;
        let slope = 0.5
            * dt
            * if self.inertial() {
                1.0 / self.coefficient
            } else {
                self.coefficient
            };
        for (out, value) in output.iter_mut().zip(direction) {
            *out = slope * value;
        }
        Ok(())
    }

    fn vjp(
        &self,
        old: &Vec<f64>,
        dt: f64,
        input: &[f64],
        direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        self.jvp(old, dt, input, direction, output, cx)
    }
}

#[test]
fn nonmatching_field_oscillators_follow_the_weighted_monolithic_equations() {
    let metric_a = [2.0, 0.6, 0.6, 1.3];
    let metric_b = [1.8, -0.2, 0.4, -0.2, 2.5, 0.3, 0.4, 0.3, 1.5];
    // A full-rank interpolation and one with an unconstrained inertial mode.
    for map in [
        [[1.0, 2.0], [0.5, -1.0], [-2.0, 0.25]],
        [[1.0, 2.0], [2.0, 4.0], [-0.5, -1.0]],
    ] {
        let (mass, stiffness) = (1.7, 9.3);
        let (mut velocity, mut extension) = ([0.35, -0.7], [0.1, -0.2, 0.4]);
        let energy = |v: &[f64], q: &[f64]| {
            0.5 * mass * pair(&metric_a, v, v) + 0.5 * stiffness * pair(&metric_b, q, q)
        };
        let energy0 = energy(&velocity, &extension);
        let mut builder = SystemBuilder::new();
        let inertias = builder
            .add_with_tangent(
                "inertias",
                FieldStorage::new(PortRole::EffortIn, &metric_a, 2, mass),
                velocity.to_vec(),
            )
            .unwrap();
        let springs = builder
            .add_with_tangent(
                "springs",
                FieldStorage::new(PortRole::FlowIn, &metric_b, 3, stiffness),
                extension.to_vec(),
            )
            .unwrap();
        let a = builder.port(inertias, "trace").unwrap();
        let b = builder.port(springs, "trace").unwrap();
        builder
            .power_map(a, b, map.into_iter().flatten().collect())
            .unwrap();
        let mut system = builder.seal().unwrap();
        system.set_config(SolveConfig {
            max_iterations: 1,
            ..SolveConfig::default()
        });
        let dt = 0.02;
        // Independent Galerkin condensation: (m M_a + h² k/4 Pᵀ M_b P)
        // v_mean = m M_a v_old + h k/2 Pᵀ M_b q_old. The 2x2 inverse is
        // explicit, with no use of the kernel's metric factor or map.
        let mut schur = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                schur[i][j] = mass * metric_a[2 * i + j];
                for r in 0..3 {
                    for s in 0..3 {
                        schur[i][j] += 0.25
                            * dt
                            * dt
                            * stiffness
                            * map[r][i]
                            * metric_b[3 * r + s]
                            * map[s][j];
                    }
                }
            }
        }
        let det = schur[0][0] * schur[1][1] - schur[0][1] * schur[1][0];
        with_cx(|cx| {
            for _ in 0..600 {
                let mut rhs = [0.0; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        rhs[i] += mass * metric_a[2 * i + j] * velocity[j];
                    }
                    for r in 0..3 {
                        for s in 0..3 {
                            rhs[i] += 0.5
                                * dt
                                * stiffness
                                * map[r][i]
                                * metric_b[3 * r + s]
                                * extension[s];
                        }
                    }
                }
                let mean = [
                    (schur[1][1] * rhs[0] - schur[0][1] * rhs[1]) / det,
                    (schur[0][0] * rhs[1] - schur[1][0] * rhs[0]) / det,
                ];
                for i in 0..2 {
                    velocity[i] = 2.0 * mean[i] - velocity[i];
                }
                for r in 0..3 {
                    extension[r] -= dt * (map[r][0] * mean[0] + map[r][1] * mean[1]);
                }
                system.step(dt, cx).unwrap();
                let v = system.state::<Vec<f64>>(inertias).unwrap();
                let q = system.state::<Vec<f64>>(springs).unwrap();
                for i in 0..2 {
                    close(v[i], velocity[i], 2.0e-10);
                }
                for r in 0..3 {
                    close(q[r], extension[r], 2.0e-10);
                }
                close(energy(v, q), energy0, 2.0e-11 * energy0);
                let power_a = pair(&metric_a, system.last_input(a), system.last_output(a));
                let power_b = pair(&metric_b, system.last_output(b), system.last_input(b));
                close(power_a + power_b, 0.0, 2.0e-11 * energy0 / dt);
            }
        });
        close(system.ledger().energy_residual_j(), 0.0, 2.0e-11 * energy0);
        close(
            system.ledger().connections[0].energy_defect_j.value(),
            0.0,
            2.0e-11 * energy0,
        );
    }
}

#[test]
fn incompatible_metrics_refuse_plain_connections_without_claiming_the_ports() {
    let mut builder = SystemBuilder::new();
    let inertia = builder
        .add_with_tangent(
            "inertia",
            FieldStorage::new(PortRole::EffortIn, &[2.0, 0.5, 0.5, 1.5], 2, 1.0),
            vec![0.3, -0.2],
        )
        .unwrap();
    let elasticity = builder
        .add_with_tangent(
            "elasticity",
            FieldStorage::new(PortRole::FlowIn, &[1.0, -0.25, -0.25, 3.0], 2, 2.0),
            vec![0.1, 0.4],
        )
        .unwrap();
    let a = builder.port(inertia, "trace").unwrap();
    let b = builder.port(elasticity, "trace").unwrap();
    for result in [
        builder.bond(a, b),
        builder.common_flow(a, &[b]),
        builder.common_effort(b, &[a]),
        builder.common_flow_ratios(a, &[(b, 2.0)]),
        builder.common_effort_ratios(b, &[(a, 2.0)]),
        builder.transformer(a, b, 2.0),
    ] {
        assert!(matches!(result, Err(AdmissionError::Mismatch { .. })));
    }
    // Identity interpolation remains valid, with a nontrivial effort
    // pullback to respect the different physical metrics.
    builder.power_map(a, b, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.01, cx).unwrap();
    });
    assert_eq!(system.steps(), 1);
    assert_eq!(system.ledger().connections.len(), 1);
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-13);
}

#[test]
fn field_maps_allow_an_ordinary_port_on_either_side_and_identity_bonds() {
    let identity = [1.0, 0.0, 0.0, 1.0];
    let weighted = [2.0, 0.5, 0.5, 1.5];
    for (plain_left, nontrivial_metric) in
        [(false, false), (true, false), (false, true), (true, true)]
    {
        let field_matrix = if nontrivial_metric {
            &weighted
        } else {
            &identity
        };
        let matrix_a = if plain_left { &identity } else { field_matrix };
        let matrix_b = if plain_left { field_matrix } else { &identity };
        let mut inertia = FieldStorage::new(PortRole::EffortIn, matrix_a, 2, 1.0);
        let mut elasticity = FieldStorage::new(PortRole::FlowIn, matrix_b, 2, 2.0);
        let plain = if plain_left {
            &mut inertia
        } else {
            &mut elasticity
        };
        plain.ports[0] = PortSpec::power(
            StableId::new("trace").unwrap(),
            PortKind::MechanicalForceVelocity,
            plain.ports[0].role(),
            NonZeroUsize::new(2).unwrap(),
        )
        .unwrap();
        let energy = |v: &[f64], q: &[f64]| 0.5 * pair(matrix_a, v, v) + pair(matrix_b, q, q);
        let energy0 = energy(&[0.3, -0.2], &[0.1, 0.4]);
        let mut builder = SystemBuilder::new();
        let inertias = builder
            .add_with_tangent("inertia", inertia, vec![0.3, -0.2])
            .unwrap();
        let springs = builder
            .add_with_tangent("elasticity", elasticity, vec![0.1, 0.4])
            .unwrap();
        let a = builder.port(inertias, "trace").unwrap();
        let b = builder.port(springs, "trace").unwrap();
        if nontrivial_metric {
            builder.power_map(a, b, identity.to_vec()).unwrap();
        } else {
            // An explicit identity metric is the ordinary power pairing.
            builder.bond(a, b).unwrap();
        }
        let mut system = builder.seal().unwrap();
        with_cx(|cx| {
            for _ in 0..200 {
                system.step(0.01, cx).unwrap();
            }
        });
        let v = system.state::<Vec<f64>>(inertias).unwrap();
        let q = system.state::<Vec<f64>>(springs).unwrap();
        assert_ne!(v.as_slice(), &[0.3, -0.2]);
        close(energy(v, q), energy0, 1.0e-11 * energy0);
        close(system.ledger().energy_residual_j(), 0.0, 1.0e-11 * energy0);
    }
}

#[test]
fn an_unrepresentable_metric_adjoint_does_not_claim_either_port() {
    let mut builder = SystemBuilder::new();
    let inertia = builder
        .add_with_tangent(
            "inertia",
            FieldStorage::new(PortRole::EffortIn, &[1.0], 1, 1.0),
            vec![0.0],
        )
        .unwrap();
    let elasticity = builder
        .add_with_tangent(
            "elasticity",
            FieldStorage::new(PortRole::FlowIn, &[1.0e300], 1, 1.0),
            vec![0.0],
        )
        .unwrap();
    let a = builder.port(inertia, "trace").unwrap();
    let b = builder.port(elasticity, "trace").unwrap();
    // Finite declared coefficients can still produce an unrepresentable
    // derived P^T M_b. Such a refusal must leave both endpoints available.
    assert!(matches!(
        builder.power_map(a, b, vec![1.0e100]),
        Err(AdmissionError::MapAdjoint { .. })
    ));
    builder.power_map(a, b, vec![1.0e-300]).unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.01, cx).unwrap();
    });
    assert_eq!(system.steps(), 1);
    assert_eq!(system.ledger().connections.len(), 1);
    assert_eq!(system.ledger().energy_residual_j(), 0.0);
}

/// Prescribed field effort or flow, keeping the supplied work as its state.
/// Work is formed as M_ij * (e_i f_j), independently of the production
/// pairing and scale. These examples keep that product representable even
/// where a different multiplication order would overflow.
struct FixedFieldEndpoint {
    ports: [PortSpec; 1],
    matrix: Vec<f64>,
    returned: Vec<f64>,
}

impl FixedFieldEndpoint {
    fn new(role: PortRole, matrix: &[f64], returned: Vec<f64>) -> Self {
        Self {
            ports: [field_port(role, matrix, returned.len())],
            matrix: matrix.to_vec(),
            returned,
        }
    }
}

impl Component for FixedFieldEndpoint {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &f64) -> Result<f64, Cancelled> {
        Ok(*state)
    }

    fn propose(
        &self,
        old: &f64,
        dt: f64,
        input: &[f64],
        new: &mut f64,
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        output.copy_from_slice(&self.returned);
        let (effort, flow) = if self.ports[0].role() == PortRole::EffortIn {
            (input, &self.returned[..])
        } else {
            (&self.returned[..], input)
        };
        let mut power = 0.0;
        for (i, &e) in effort.iter().enumerate() {
            for (j, &f) in flow.iter().enumerate() {
                power += self.matrix[i * flow.len() + j] * (e * f);
            }
        }
        let change = dt * power;
        *new = old + change;
        Ok(Proposal::Step(Receipt::lossless(
            change,
            128.0 * f64::EPSILON * (old.abs() + new.abs() + change.abs()),
        )))
    }
}

impl PortTangent for FixedFieldEndpoint {
    fn jvp(
        &self,
        _old: &f64,
        _dt: f64,
        _input: &[f64],
        _direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        cx.checkpoint()?;
        output.fill(0.0);
        Ok(())
    }

    fn vjp(
        &self,
        old: &f64,
        dt: f64,
        input: &[f64],
        direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        self.jvp(old, dt, input, direction, output, cx)
    }
}

#[test]
fn field_power_scaling_keeps_finite_work_and_refuses_an_infinite_allowance() {
    for (matrix, effort, flow, should_accept) in [
        // M*f is one and e*(M*f) is finite; e*M would overflow first.
        (vec![1.0e200], vec![1.0e200], vec![1.0e-200], true),
        // Signed powers cancel, but their absolute term sum exceeds f64::MAX
        // and cannot be used as an infinite acceptance tolerance.
        (
            vec![1.0e308, 5.0e307, 5.0e307, 1.0e308],
            vec![1.0, -1.0],
            vec![1.0, 1.0],
            false,
        ),
    ] {
        let mut builder = SystemBuilder::new();
        let a = builder
            .add_with_tangent(
                "flow_source",
                FixedFieldEndpoint::new(PortRole::EffortIn, &matrix, flow),
                0.0,
            )
            .unwrap();
        let b = builder
            .add_with_tangent(
                "effort_source",
                FixedFieldEndpoint::new(PortRole::FlowIn, &matrix, effort),
                0.0,
            )
            .unwrap();
        builder
            .bond(
                builder.port(a, "trace").unwrap(),
                builder.port(b, "trace").unwrap(),
            )
            .unwrap();
        let mut system = builder.seal().unwrap();
        with_cx(|cx| {
            let ledger_before = system.ledger().clone();
            let result = system.step(0.1, cx);
            if should_accept {
                result.unwrap();
                close(*system.state::<f64>(a).unwrap(), 1.0e199, 1.0e184);
                close(*system.state::<f64>(b).unwrap(), -1.0e199, 1.0e184);
                assert_eq!(system.steps(), 1);
            } else {
                assert!(matches!(
                    result,
                    Err(SystemError::Interface { tolerance_j, .. }) if !tolerance_j.is_finite()
                ));
                assert_eq!(*system.state::<f64>(a).unwrap(), 0.0);
                assert_eq!(*system.state::<f64>(b).unwrap(), 0.0);
                assert_eq!(system.steps(), 0);
                assert_eq!(system.time_s(), 0.0);
                assert_eq!(system.ledger(), &ledger_before);
            }
        });
    }
}
