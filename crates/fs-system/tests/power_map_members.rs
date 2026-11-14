//! G0/G2: one power map from several effort-receiving ports, as a lever takes
//! the motion of two bodies onto one contact. Trajectories are checked against
//! the midpoint equations of the stacked map, solved here by elimination, not
//! against the kernel's ledger; the scalar case is checked against the
//! ratio junction to the bit.

#![allow(clippy::float_cmp, clippy::needless_range_loop)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::{
    Carrier, Class, Component, DualityMetric, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};
use fs_system::{AdmissionError, PortRef, SolveConfig, SystemBuilder};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
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

fn power(role: PortRole, width: usize) -> PortSpec {
    PortSpec::power(
        StableId::new("p").unwrap(),
        PortKind::MechanicalForceVelocity,
        role,
        NonZeroUsize::new(width).unwrap(),
    )
    .unwrap()
}

fn field_port(role: PortRole, matrix: &[f64], width: usize) -> PortSpec {
    let metric = DualityMetric::new(NonZeroUsize::new(width).unwrap(), matrix.to_vec()).unwrap();
    PortSpec::field_power(
        StableId::new("p").unwrap(),
        PortKind::MechanicalForceVelocity,
        role,
        metric,
    )
    .unwrap()
}

/// `a . G b` for a row-major square `G`.
fn pair(matrix: &[f64], a: &[f64], b: &[f64]) -> f64 {
    let mut value = 0.0;
    for i in 0..a.len() {
        for j in 0..b.len() {
            value += a[i] * matrix[i * b.len() + j] * b[j];
        }
    }
    value
}

/// A storage bank with a constant metric `G` (the identity for an ordinary
/// port) and one coefficient: an inertia when effort is the input, a spring
/// when flow is the input. State is velocity or extension; the energy is
/// `coefficient x . G x / 2`.
struct Storage {
    ports: [PortSpec; 1],
    metric: Vec<f64>,
    coefficient: f64,
}

impl Storage {
    fn plain(role: PortRole, width: usize, coefficient: f64) -> Self {
        let mut metric = vec![0.0; width * width];
        for i in 0..width {
            metric[i * width + i] = 1.0;
        }
        Self {
            ports: [power(role, width)],
            metric,
            coefficient,
        }
    }

    fn field(role: PortRole, metric: &[f64], width: usize, coefficient: f64) -> Self {
        Self {
            ports: [field_port(role, metric, width)],
            metric: metric.to_vec(),
            coefficient,
        }
    }

    fn inertial(&self) -> bool {
        self.ports[0].role() == PortRole::EffortIn
    }
}

impl Component for Storage {
    type State = Vec<f64>;
    type Error = Cancelled;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, Cancelled> {
        Ok(0.5 * self.coefficient * pair(&self.metric, state, state))
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
        let mean: Vec<f64> = old
            .iter()
            .zip(new.iter())
            .map(|(o, n)| 0.5 * (o + n))
            .collect();
        let step: Vec<f64> = old.iter().zip(new.iter()).map(|(o, n)| n - o).collect();
        let change = self.coefficient * pair(&self.metric, &mean, &step);
        let scale = self.stored_energy_j(old)? + self.stored_energy_j(new)? + change.abs();
        Ok(Proposal::Step(Receipt::lossless(
            change,
            256.0 * f64::EPSILON * scale,
        )))
    }
}

impl PortTangent for Storage {
    fn jvp(
        &self,
        _old: &Vec<f64>,
        dt: f64,
        _input: &[f64],
        direction: &[f64],
        output: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
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

/// Gaussian elimination with partial pivoting.
fn solve(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&i, &j| matrix[i][column].abs().total_cmp(&matrix[j][column].abs()))
            .unwrap();
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        for row in column + 1..n {
            let factor = matrix[row][column] / matrix[column][column];
            for k in column..n {
                matrix[row][k] -= factor * matrix[column][k];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = (row + 1..n).map(|k| matrix[row][k] * x[k]).sum();
        x[row] = (rhs[row] - tail) / matrix[row][row];
    }
    x
}

/// One midpoint step of inertias with energy matrix `inertia` (n x n, the
/// members' blocks on its diagonal) joined by the stacked map `map` (rows x
/// n) to springs with energy matrix `stiffness` (rows x rows):
/// `(A + h^2/4 P^T K P) v_mean = A v + h/2 P^T K q`, then `v = 2 v_mean - v`
/// and `q -= h P v_mean`.
fn midpoint_step(
    inertia: &[Vec<f64>],
    stiffness: &[Vec<f64>],
    map: &[Vec<f64>],
    dt: f64,
    velocity: &mut [f64],
    extension: &mut [f64],
) {
    let (n, rows) = (velocity.len(), extension.len());
    // P^T K, n x rows.
    let mut weighted = vec![vec![0.0; rows]; n];
    for i in 0..n {
        for s in 0..rows {
            weighted[i][s] = (0..rows).map(|r| map[r][i] * stiffness[r][s]).sum();
        }
    }
    let mut schur = inertia.to_vec();
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            schur[i][j] +=
                0.25 * dt * dt * (0..rows).map(|s| weighted[i][s] * map[s][j]).sum::<f64>();
            rhs[i] += inertia[i][j] * velocity[j];
        }
        rhs[i] += 0.5
            * dt
            * (0..rows)
                .map(|s| weighted[i][s] * extension[s])
                .sum::<f64>();
    }
    let mean = solve(schur, rhs);
    for i in 0..n {
        velocity[i] = 2.0 * mean[i] - velocity[i];
    }
    for r in 0..rows {
        extension[r] -= dt * (0..n).map(|i| map[r][i] * mean[i]).sum::<f64>();
    }
}

/// `blocks` placed on the diagonal of one square matrix.
fn block_diagonal(blocks: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let widths: Vec<usize> = blocks.iter().map(|b| b.len().isqrt()).collect();
    let n: usize = widths.iter().sum();
    let mut matrix = vec![vec![0.0; n]; n];
    let mut start = 0;
    for (block, width) in blocks.iter().zip(widths) {
        for i in 0..width {
            for j in 0..width {
                matrix[start + i][start + j] = block[i * width + j];
            }
        }
        start += width;
    }
    matrix
}

fn scaled(matrix: &[f64], by: f64) -> Vec<f64> {
    matrix.iter().map(|value| by * value).collect()
}

/// Row-major `rows x columns` member maps placed side by side: `[P_1 ... P_n]`.
fn stacked(maps: &[(&[f64], usize)], rows: usize) -> Vec<Vec<f64>> {
    (0..rows)
        .map(|r| {
            maps.iter()
                .flat_map(|(map, columns)| map[r * columns..(r + 1) * columns].iter().copied())
                .collect()
        })
        .collect()
}

/// The ratio junction and the power map with members `r I` after the same
/// steps: the members' states and the setter's, to the bit.
fn levered_bank(as_map: bool) -> Vec<Vec<f64>> {
    let mut builder = SystemBuilder::new();
    let springs = builder
        .add_with_tangent(
            "springs",
            Storage::plain(PortRole::FlowIn, 2, 7.5),
            vec![0.1, -0.3],
        )
        .unwrap();
    let mut nodes = vec![springs];
    let mut members = Vec::new();
    for (name, mass, velocity, ratio) in [
        ("near", 1.2, [0.4, -0.2], 1.6),
        ("far", 2.1, [-0.1, 0.5], -0.45),
    ] {
        let node = builder
            .add_with_tangent(
                name,
                Storage::plain(PortRole::EffortIn, 2, mass),
                velocity.to_vec(),
            )
            .unwrap();
        members.push((builder.port(node, "p").unwrap(), ratio));
        nodes.push(node);
    }
    let setter = builder.port(springs, "p").unwrap();
    if as_map {
        let members = members
            .iter()
            .map(|&(port, ratio)| (port, vec![ratio, 0.0, 0.0, ratio]))
            .collect();
        builder.power_map_members(members, setter).unwrap();
    } else {
        builder.common_effort_ratios(setter, &members).unwrap();
    }
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..800 {
            system.step(0.02, cx).unwrap();
        }
    });
    nodes
        .iter()
        .map(|&node| system.state::<Vec<f64>>(node).unwrap().clone())
        .collect()
}

#[test]
fn scalar_members_are_the_common_effort_junction_to_the_bit() {
    let by_junction = levered_bank(false);
    let by_map = levered_bank(true);
    assert_ne!(by_junction[0], vec![0.1, -0.3]);
    for (junction, map) in by_junction.iter().zip(&by_map) {
        let bits = |values: &[f64]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(junction), bits(map));
    }
}

#[test]
fn two_members_follow_the_midpoint_equations_of_their_stacked_map() {
    let left_map = [1.0, 2.0, 0.5, -1.0, -2.0, 0.25];
    // A column the stacked map adds to full rank, the sum of the left
    // member's columns (rank two, one unconstrained inertial mode), and a
    // member the map does not reach.
    for right_map in [[0.7, 1.5, -0.4], [3.0, -0.5, -1.75], [0.0; 3]] {
        let (left_mass, right_mass, stiffness) = (1.5, 0.8, 9.0);
        let (mut velocity, mut extension) = (vec![0.35, -0.7, 0.45], vec![0.1, -0.2, 0.4]);
        let inertia = block_diagonal(&[vec![left_mass, 0.0, 0.0, left_mass], vec![right_mass]]);
        let springs_matrix = block_diagonal(&[scaled(
            &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            stiffness,
        )]);
        let map = stacked(&[(&left_map[..], 2), (&right_map[..], 1)], 3);
        let energy = |v: &[f64], q: &[f64]| -> f64 {
            0.5 * left_mass * (v[0] * v[0] + v[1] * v[1])
                + 0.5 * right_mass * v[2] * v[2]
                + 0.5 * stiffness * q.iter().map(|x| x * x).sum::<f64>()
        };
        let energy0 = energy(&velocity, &extension);

        let mut builder = SystemBuilder::new();
        let left = builder
            .add_with_tangent(
                "left",
                Storage::plain(PortRole::EffortIn, 2, left_mass),
                velocity[..2].to_vec(),
            )
            .unwrap();
        let right = builder
            .add_with_tangent(
                "right",
                Storage::plain(PortRole::EffortIn, 1, right_mass),
                velocity[2..].to_vec(),
            )
            .unwrap();
        let springs = builder
            .add_with_tangent(
                "springs",
                Storage::plain(PortRole::FlowIn, 3, stiffness),
                extension.clone(),
            )
            .unwrap();
        let members = vec![
            (builder.port(left, "p").unwrap(), left_map.to_vec()),
            (builder.port(right, "p").unwrap(), right_map.to_vec()),
        ];
        let b = builder.port(springs, "p").unwrap();
        builder.power_map_members(members, b).unwrap();
        let mut system = builder.seal().unwrap();
        // The loop is affine, so one exact Newton update closes it.
        system.set_config(SolveConfig {
            max_iterations: 1,
            ..SolveConfig::default()
        });
        let dt = 0.03;
        with_cx(|cx| {
            for _ in 0..1_000 {
                midpoint_step(
                    &inertia,
                    &springs_matrix,
                    &map,
                    dt,
                    &mut velocity,
                    &mut extension,
                );
                system.step(dt, cx).unwrap();
                let l = system.state::<Vec<f64>>(left).unwrap();
                let r = system.state::<Vec<f64>>(right).unwrap();
                let q = system.state::<Vec<f64>>(springs).unwrap();
                let v = [l[0], l[1], r[0]];
                for i in 0..3 {
                    close(v[i], velocity[i], 2.0e-10);
                    close(q[i], extension[i], 2.0e-10);
                }
                close(energy(&v, q), energy0, 2.0e-11 * energy0);
            }
        });
        if right_map == [0.0; 3] {
            // Nothing reaches the right member, so it keeps its motion.
            assert_eq!(system.state::<Vec<f64>>(right).unwrap(), &vec![0.45]);
        }
        assert_eq!(system.ledger().connections.len(), 1);
        close(
            system.ledger().connections[0].energy_defect_j.value(),
            0.0,
            1.0e-11 * energy0,
        );
        close(system.ledger().energy_residual_j(), 0.0, 1.0e-11 * energy0);
    }
}

#[test]
fn each_member_takes_its_own_metric_adjoint() {
    // A field trace of width two and an ordinary port, mapped onto a field
    // trace of width three. Member k receives M_k^-1 P_k^T M_b e_b, so the
    // power the three ports receive sums to zero.
    let metric_left = [2.0, 0.6, 0.6, 1.3];
    let metric_b = [1.8, -0.2, 0.4, -0.2, 2.5, 0.3, 0.4, 0.3, 1.5];
    let left_map = [1.0, 2.0, 0.5, -1.0, -2.0, 0.25];
    let right_map = [0.7, 1.5, -0.4];
    let (left_mass, right_mass, stiffness) = (1.7, 0.9, 9.3);
    let (mut velocity, mut extension) = (vec![0.35, -0.7, 0.45], vec![0.1, -0.2, 0.4]);
    let inertia = block_diagonal(&[scaled(&metric_left, left_mass), vec![right_mass]]);
    let springs_matrix = block_diagonal(&[scaled(&metric_b, stiffness)]);
    let map = stacked(&[(&left_map[..], 2), (&right_map[..], 1)], 3);
    let energy = |v: &[f64], q: &[f64]| -> f64 {
        0.5 * left_mass * pair(&metric_left, &v[..2], &v[..2])
            + 0.5 * right_mass * v[2] * v[2]
            + 0.5 * stiffness * pair(&metric_b, q, q)
    };
    let energy0 = energy(&velocity, &extension);

    let mut builder = SystemBuilder::new();
    let left = builder
        .add_with_tangent(
            "left",
            Storage::field(PortRole::EffortIn, &metric_left, 2, left_mass),
            velocity[..2].to_vec(),
        )
        .unwrap();
    let right = builder
        .add_with_tangent(
            "right",
            Storage::plain(PortRole::EffortIn, 1, right_mass),
            velocity[2..].to_vec(),
        )
        .unwrap();
    let springs = builder
        .add_with_tangent(
            "springs",
            Storage::field(PortRole::FlowIn, &metric_b, 3, stiffness),
            extension.clone(),
        )
        .unwrap();
    let (a_left, a_right, b) = (
        builder.port(left, "p").unwrap(),
        builder.port(right, "p").unwrap(),
        builder.port(springs, "p").unwrap(),
    );
    builder
        .power_map_members(
            vec![(a_left, left_map.to_vec()), (a_right, right_map.to_vec())],
            b,
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig {
        max_iterations: 1,
        ..SolveConfig::default()
    });
    let dt = 0.02;
    with_cx(|cx| {
        for _ in 0..600 {
            midpoint_step(
                &inertia,
                &springs_matrix,
                &map,
                dt,
                &mut velocity,
                &mut extension,
            );
            system.step(dt, cx).unwrap();
            let l = system.state::<Vec<f64>>(left).unwrap();
            let r = system.state::<Vec<f64>>(right).unwrap();
            let q = system.state::<Vec<f64>>(springs).unwrap();
            let v = [l[0], l[1], r[0]];
            for i in 0..3 {
                close(v[i], velocity[i], 2.0e-10);
                close(q[i], extension[i], 2.0e-10);
            }
            close(energy(&v, q), energy0, 2.0e-11 * energy0);
            let received = pair(
                &metric_left,
                system.last_input(a_left),
                system.last_output(a_left),
            ) + system.last_input(a_right)[0] * system.last_output(a_right)[0]
                + pair(&metric_b, system.last_output(b), system.last_input(b));
            close(received, 0.0, 2.0e-11 * energy0 / dt);
        }
    });
    close(system.ledger().energy_residual_j(), 0.0, 2.0e-11 * energy0);
}

/// Zero exchange, except that a heat receiver reports 300 K. Used only to
/// exercise admission and prove a repaired builder can step.
struct TerminatedPort([PortSpec; 1]);
impl Component for TerminatedPort {
    type State = ();
    type Error = Cancelled;
    fn ports(&self) -> &[PortSpec] {
        &self.0
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, (): &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }
    fn propose(
        &self,
        (): &(),
        _: f64,
        _: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        outputs.fill(
            if self.0[0].carrier() == Carrier::Heat && self.0[0].role() == PortRole::FlowIn {
                300.0
            } else {
                0.0
            },
        );
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

fn terminated(builder: &mut SystemBuilder, name: &str, spec: PortSpec) -> PortRef {
    let node = builder.add(name, TerminatedPort([spec]), ()).unwrap();
    builder.port(node, "p").unwrap()
}

#[test]
fn refused_members_leave_every_port_usable() {
    let mut builder = SystemBuilder::new();
    let left = terminated(&mut builder, "left", power(PortRole::EffortIn, 2));
    let right = terminated(&mut builder, "right", power(PortRole::EffortIn, 1));
    let tip = terminated(&mut builder, "tip", power(PortRole::FlowIn, 3));
    let pushed = terminated(&mut builder, "pushed", power(PortRole::FlowIn, 1));
    let spare_a = terminated(&mut builder, "spare_a", power(PortRole::EffortIn, 1));
    let spare_b = terminated(&mut builder, "spare_b", power(PortRole::FlowIn, 1));
    let hot = terminated(
        &mut builder,
        "hot",
        PortSpec::heat(StableId::new("p").unwrap(), PortRole::EffortIn).unwrap(),
    );
    let left_map = || vec![1.0, 2.0, 0.5, -1.0, -2.0, 0.25];

    assert!(matches!(
        builder.power_map_members(Vec::new(), tip),
        Err(AdmissionError::Empty)
    ));
    // Members receive effort and the common port receives flow; heat is not
    // a power pair. The common port named as a member is refused by its role.
    for (members, b) in [
        (vec![(left, left_map()), (pushed, vec![1.0; 3])], tip),
        (vec![(left, left_map()), (hot, vec![1.0; 3])], tip),
        (vec![(left, left_map()), (tip, vec![1.0; 9])], tip),
        (vec![(spare_a, vec![1.0])], right),
    ] {
        assert!(matches!(
            builder.power_map_members(members, b),
            Err(AdmissionError::Role { .. })
        ));
    }
    // Each member's map has the common port's rows and its own columns.
    assert!(matches!(
        builder.power_map_members(vec![(left, left_map()), (right, vec![1.0; 2])], tip),
        Err(AdmissionError::MapShape {
            rows: 3,
            columns: 1,
            entries: 2
        })
    ));
    // The index of a non-finite coefficient counts through the maps in order.
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            builder.power_map_members(vec![(left, left_map()), (right, vec![1.0, value, 1.0])], tip),
            Err(AdmissionError::MapCoefficient { index: 7, value: rejected }) if rejected.to_bits() == value.to_bits()
        ));
    }
    assert!(matches!(
        builder.power_map_members(vec![(right, vec![1.0; 3]), (right, vec![1.0; 3])], tip),
        Err(AdmissionError::AlreadyConnected { .. })
    ));

    builder
        .power_map_members(vec![(left, left_map()), (right, vec![0.7, 1.5, -0.4])], tip)
        .unwrap();
    // A refusal on a claimed member must not claim the others.
    assert!(matches!(
        builder.power_map_members(vec![(spare_a, vec![1.0]), (right, vec![1.0])], spare_b),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    builder.power_map(spare_a, spare_b, vec![2.0]).unwrap();
    builder
        .boundary(pushed, Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    builder
        .boundary(hot, Box::new(|_, _, values| values.fill(300.0)))
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.01, cx).unwrap();
    });
    assert_eq!(system.steps(), 1);
    assert_eq!(system.ledger().connections.len(), 2);
    assert_eq!(system.ledger().energy_residual_j(), 0.0);
}
