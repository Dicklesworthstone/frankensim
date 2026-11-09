//! G1/G3: native nonlinear material memory, coherent drives and temporal refinement.

use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactTones, PreparedNonlinearTones,
};
use fs_phs::mechanical::LinearMechanicalSystem;
use fs_phs::{
    MappedRelaxationBranch, PortHamiltonian, QuadraticStorage, RelaxationBranch,
    RelaxationStrainMap, ResistanceWitness, StepWorkspace,
};
use std::cell::Cell;
use std::rc::Rc;

#[derive(Debug)]
struct Geometry {
    dimensions: usize,
    linear: [f64; 2],
    hessian: [f64; 3],
    calls: Rc<(Cell<usize>, Cell<usize>)>,
}
impl RelaxationStrainMap for Geometry {
    fn state_dimension(&self) -> usize {
        self.dimensions
    }
    fn strain(&self, x: &[f64]) -> f64 {
        self.linear[0] * x[0]
            + self.linear[1] * x[1]
            + 0.5 * self.hessian[0] * x[0] * x[0]
            + self.hessian[1] * x[0] * x[1]
            + 0.5 * self.hessian[2] * x[1] * x[1]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.calls.0.set(self.calls.0.get() + 1);
        out.fill(0.0);
        out[0] = self.linear[0] + self.hessian[0] * x[0] + self.hessian[1] * x[1];
        out[1] = self.linear[1] + self.hessian[1] * x[0] + self.hessian[2] * x[1];
    }
    fn hessian_vector(&self, _: &[f64], v: &[f64], out: &mut [f64]) {
        self.calls.1.set(self.calls.1.get() + 1);
        out.fill(0.0);
        out[0] = self.hessian[0] * v[0] + self.hessian[1] * v[1];
        out[1] = self.hessian[1] * v[0] + self.hessian[2] * v[1];
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn dense() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}
fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 128,
        max_omega_dt: 0.04,
    }
}

fn scalar(calls: Rc<(Cell<usize>, Cell<usize>)>) -> PortHamiltonian {
    LinearMechanicalSystem::new(&[2.0], &[0.4], &[3.0], &[1.0], 1, 1, vec![])
        .unwrap()
        .into_system()
        .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
            strain: Rc::new(Geometry {
                dimensions: 2,
                linear: [1.0, 0.0],
                hessian: [1.0, 0.0, 0.0],
                calls,
            }),
            stiffness: 4.0,
            relaxation_time_s: 0.2,
        }])
        .unwrap()
}

#[test]
fn component_work_preserves_kinetics_and_the_actual_prestressed_endpoint_tangent() {
    let calls = Rc::new((Cell::new(0), Cell::new(0)));
    let system = scalar(calls.clone());
    let counts = (calls.0.get(), calls.1.get());
    assert!(counts.0 > 0 && counts.1 > 0);
    let (a, b) = ([0.3, 5.0, 0.1], [0.7, 9.0, 0.25]);
    let mid = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(a[i], b[i]));
    let rule = system
        .storage()
        .discrete_gradient_rule()
        .expect("native Maxwell rule");
    let mut effort = [f64::NAN; 3];
    rule.evaluate(&a, &b, &mid, &mut effort);
    let stress = |x: [f64; 3]| 2.0 * (x[0] + 0.5 * x[0] * x[0]) - x[2];
    let mean_stress = 0.5 * (stress(a) + stress(b));
    let strain_gradient = 2.0 * (1.0 + mid[0]);
    close(
        effort[0],
        3.0 * mid[0] + mean_stress * strain_gradient,
        2e-15,
    );
    let mut midpoint_effort = [f64::NAN; 3];
    system.storage().gradient(&mid, &mut midpoint_effort);
    // Preserve the original represented inverse mass, including its rounding.
    assert_eq!(effort[1], midpoint_effort[1]);
    close(effort[1], 0.5 * mid[1], 2e-15);
    close(effort[2], -mean_stress, 2e-15);
    let change = |x: [f64; 3]| 1.5 * x[0] * x[0] + 0.25 * x[1] * x[1] + 0.5 * stress(x).powi(2);
    close(
        effort
            .iter()
            .enumerate()
            .map(|(i, e)| e * (b[i] - a[i]))
            .sum(),
        change(b) - change(a),
        3e-15,
    );
    let direction = [0.4, -0.2, 0.3];
    let mut action = [f64::NAN; 3];
    assert!(rule.end_jacobian_action(&a, &b, &direction, &mut action));
    let de = 0.5 * (2.0 * (1.0 + b[0]) * direction[0] - direction[2]);
    close(
        action[0],
        1.5 * direction[0] + de * strain_gradient + mean_stress * direction[0],
        2e-15,
    );
    close(action[1], 0.25 * direction[1], 2e-15);
    close(action[2], -de, 2e-15);
    assert!(rule.end_jacobian_action(&a, &a, &direction, &mut action));
    let mut continuous = [f64::NAN; 3];
    assert!(
        system
            .storage()
            .continuous_tangent()
            .unwrap()
            .apply(&a, &direction, &mut continuous)
    );
    for (a, c) in action.into_iter().zip(continuous) {
        close(a, 0.5 * c, 2e-15);
    }
    // Large stored prestress must not inject the branch's geometric work into
    // an independent kinetic coordinate, or obscure a pure memory variation.
    let high_a = [0.3, 1e8, 1e12];
    let high_b = [0.7, 1e8 + 1.0, 1e12];
    let high_mid = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(high_a[i], high_b[i]));
    rule.evaluate(&high_a, &high_b, &high_mid, &mut effort);
    system.storage().gradient(&high_mid, &mut midpoint_effort);
    assert_eq!(effort[1], midpoint_effort[1]);
    close(effort[1], 50_000_000.25, 2e-15);
    assert!(rule.end_jacobian_action(&high_a, &high_b, &[0.0, 0.0, 1.0], &mut action));
    close(action[0], -0.5 * strain_gradient, 2e-15);
    assert_eq!(action[1], 0.0);
    assert_eq!(action[2], 0.5);
    // The native force generates its material overtones; this is not an
    // authored output waveform. q=A*cos(theta), z=0 gives third-harmonic
    // elastic reaction k*A^3/8 for epsilon=q+q^2/2.
    let amplitude = 0.4;
    let mut third = 0.0;
    for i in 0..64 {
        let theta = std::f64::consts::TAU * i as f64 / 64.0;
        let q = amplitude * theta.cos();
        let mut gradient = [0.0; 3];
        system.storage().gradient(&[q, 0.0, 0.0], &mut gradient);
        third += (gradient[0] - 3.0 * q) * (3.0 * theta).cos() / 32.0;
    }
    close(third, 4.0 * amplitude.powi(3) / 8.0, 3e-15);
    assert_eq!(counts, (calls.0.get(), calls.1.get()));
}

#[test]
fn prepared_material_step_keeps_original_kinematic_velocity_and_branch_heat() {
    let calls = Rc::new((Cell::new(0), Cell::new(0)));
    let system = scalar(calls.clone());
    let counts = (calls.0.get(), calls.1.get());
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let (initial, dt) = ([0.4, 0.7, -0.1], 0.04);
    let (mut next, mut output) = ([f64::NAN; 3], [f64::NAN; 1]);
    let record = workspace
        .step_into(&system, &initial, &[1.3], dt, &mut next, &mut output)
        .unwrap();
    close(
        (next[0] - initial[0]) / dt,
        0.25 * (next[1] + initial[1]),
        2e-12,
    );
    close(record.supplied, 1.3 * (next[0] - initial[0]), 2e-12);
    let effort = workspace.last_step_discrete_effort().unwrap();
    close(
        record.dissipated,
        dt * (0.4 * effort[1].powi(2) + 5.0 * effort[2].powi(2)),
        2e-13,
    );
    assert!(record.dissipated > 0.0);
    close(record.balance_residual(), 0.0, 2e-12);
    assert_eq!(counts, (calls.0.get(), calls.1.get()));
}

fn coupled() -> PortHamiltonian {
    let base = LinearMechanicalSystem::new(
        &[2.0, 0.2, 0.2, 1.5],
        &[0.4, 0.3, -0.05, 0.5],
        &[5.0, 0.6, 0.6, 8.0],
        &[1.0, 0.3, -0.2, 0.8],
        2,
        2,
        vec![],
    )
    .unwrap()
    .into_system();
    let definitions = [
        ([1.0, -0.2], [0.8, 0.3, 0.0], 3.0, 0.08),
        ([0.5, 1.0], [0.0, 0.0, 0.4], 1.7, 0.23),
    ];
    base.with_mapped_relaxation_branches(
        definitions
            .into_iter()
            .map(
                |(linear, hessian, stiffness, relaxation_time_s)| MappedRelaxationBranch {
                    strain: Rc::new(Geometry {
                        dimensions: 4,
                        linear,
                        hessian,
                        calls: Rc::new((Cell::new(0), Cell::new(0))),
                    }),
                    stiffness,
                    relaxation_time_s,
                },
            )
            .collect(),
    )
    .unwrap()
}

const FREQUENCIES: [f64; 3] = [0.0, 2.1, 3.7];
const C0: [f64; 6] = [0.2, 0.7, -0.1, -0.15, 0.3, 0.2];
const S0: [f64; 6] = [0.0, -0.2, 0.3, 0.0, 0.4, -0.1];
const C1: [f64; 6] = [-0.1, 0.3, 0.4, 0.2, -0.1, 0.5];
const S1: [f64; 6] = [0.0, 0.4, -0.2, 0.0, 0.2, 0.3];
const DURATION: f64 = 0.8;

fn physical_initial() -> [f64; 10] {
    [0.15, -0.08, 0.2, -0.1, 0.01, -0.03, 0.0, 0.0, 0.0, 0.0]
}

// Independent physical [q0,q1,v0,v1,z0,z1,W,Q,Q_arm0,Q_arm1] equations.
// The complete consistent mass, nonproportional/skew damping, strain Hessians
// and two material relaxation times appear explicitly; no production effort,
// energy, flow, quadrature or discrete-gradient callback enters this reference.
fn reference(affine: bool) -> [f64; 10] {
    let rhs = |t: f64, x: [f64; 10]| {
        let mut u = [0.0; 2];
        for p in 0..2 {
            for j in 0..3 {
                let i = 3 * p + j;
                let ratio = if affine { t / DURATION } else { 0.0 };
                let c = C0[i] + ratio * (C1[i] - C0[i]);
                let s = S0[i] + ratio * (S1[i] - S0[i]);
                let (sn, cs) = (FREQUENCIES[j] * t).sin_cos();
                u[p] += c * cs + s * sn;
            }
        }
        let strain0 = x[0] - 0.2 * x[1] + 0.4 * x[0] * x[0] + 0.3 * x[0] * x[1];
        let strain1 = 0.5 * x[0] + x[1] + 0.2 * x[1] * x[1];
        let e0 = 3.0_f64.sqrt() * strain0 - x[4];
        let e1 = 1.7_f64.sqrt() * strain1 - x[5];
        let f = [u[0] + 0.3 * u[1], -0.2 * u[0] + 0.8 * u[1]];
        let p0 = f[0]
            - 5.0 * x[0]
            - 0.6 * x[1]
            - 0.4 * x[2]
            - 0.3 * x[3]
            - 3.0_f64.sqrt() * (1.0 + 0.8 * x[0] + 0.3 * x[1]) * e0
            - 1.7_f64.sqrt() * 0.5 * e1;
        let p1 = f[1] - 0.6 * x[0] - 8.0 * x[1] + 0.05 * x[2]
            - 0.5 * x[3]
            - 3.0_f64.sqrt() * (-0.2 + 0.3 * x[0]) * e0
            - 1.7_f64.sqrt() * (1.0 + 0.4 * x[1]) * e1;
        let q0 = e0 * e0 / 0.08;
        let q1 = e1 * e1 / 0.23;
        [
            x[2],
            x[3],
            (1.5 * p0 - 0.2 * p1) / 2.96,
            (-0.2 * p0 + 2.0 * p1) / 2.96,
            e0 / 0.08,
            e1 / 0.23,
            f[0] * x[2] + f[1] * x[3],
            0.4 * x[2] * x[2] + 0.25 * x[2] * x[3] + 0.5 * x[3] * x[3] + q0 + q1,
            q0,
            q1,
        ]
    };
    let mut x = physical_initial();
    let dt = DURATION / 65536.0;
    let add = |x: [f64; 10], y: [f64; 10], h: f64| core::array::from_fn(|i| x[i] + h * y[i]);
    for i in 0..65536 {
        let t = i as f64 * dt;
        let a = rhs(t, x);
        let b = rhs(t + 0.5 * dt, add(x, a, 0.5 * dt));
        let c = rhs(t + 0.5 * dt, add(x, b, 0.5 * dt));
        let d = rhs(t + dt, add(x, c, dt));
        for j in 0..10 {
            x[j] += dt / 6.0 * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]);
        }
    }
    x
}

fn channels() -> [ResistanceWitness; 2] {
    core::array::from_fn(|i| {
        let rate = if i == 0 { 1.0 / 0.08 } else { 1.0 / 0.23 };
        let mut matrix = vec![0.0; 36];
        let mut row = vec![0.0; 6];
        matrix[(4 + i) * 6 + 4 + i] = rate;
        row[4 + i] = f64::sqrt(rate);
        ResistanceWitness::try_from_source_rows(6, matrix, row).unwrap()
    })
}

fn simulate(affine: bool, steps: usize) -> [f64; 10] {
    let h = DURATION / steps as f64;
    let constructor = if affine {
        PreparedNonlinearTones::new_owned_affine_with_resistance_channels
    } else {
        PreparedNonlinearTones::new_owned_with_resistance_channels
    };
    let mut prepared =
        constructor(coupled(), &FREQUENCIES, h, split(), dense(), &channels()).unwrap();
    assert!(
        prepared
            .system()
            .storage()
            .discrete_gradient_rule()
            .is_some()
    );
    assert!(!prepared.system().storage().quadratic_part().unwrap().exact);
    let r = physical_initial();
    let mut x = [
        r[0],
        r[1],
        2.0 * r[2] + 0.2 * r[3],
        0.2 * r[2] + 1.5 * r[3],
        r[4],
        r[5],
    ];
    let (mut next, mut mean, mut channel) = ([0.0; 6], [0.0; 2], [0.0; 2]);
    let mut integrals = [0.0; 4];
    for step in 0..steps {
        let t = step as f64 * h;
        let mut knots = [[0.0; 6]; 4];
        for p in 0..2 {
            for j in 0..3 {
                let i = 3 * p + j;
                let (sn, cs) = (FREQUENCIES[j] * t).sin_cos();
                for (offset, time) in [(0, t), (2, t + h)] {
                    let ratio = if affine { time / DURATION } else { 0.0 };
                    let c = C0[i] + ratio * (C1[i] - C0[i]);
                    let s = S0[i] + ratio * (S1[i] - S0[i]);
                    knots[offset][i] = c * cs + s * sn;
                    knots[offset + 1][i] = -c * sn + s * cs;
                }
            }
        }
        let record = if affine {
            prepared.step_affine_observed_into(
                &x,
                &knots[0],
                &knots[1],
                &knots[2],
                &knots[3],
                &mut next,
                &mut mean,
                &mut channel,
            )
        } else {
            prepared.step_observed_into(
                &x,
                &knots[0],
                &knots[1],
                &mut next,
                &mut mean,
                &mut channel,
            )
        }
        .unwrap();
        assert!(record.dissipated >= 0.0 && channel.iter().all(|v| *v >= 0.0));
        close(record.balance_residual(), 0.0, 2e-11);
        integrals[0] += record.supplied;
        integrals[1] += record.dissipated;
        integrals[2] += channel[0];
        integrals[3] += channel[1];
        x = next;
    }
    [
        x[0],
        x[1],
        (1.5 * x[2] - 0.2 * x[3]) / 2.96,
        (-0.2 * x[2] + 2.0 * x[3]) / 2.96,
        x[4],
        x[5],
        integrals[0],
        integrals[1],
        integrals[2],
        integrals[3],
    ]
}

#[test]
fn coherent_constant_and_affine_material_memory_converges_with_independent_heat() {
    for affine in [false, true] {
        let reference = reference(affine);
        let mut errors = [0.0_f64; 3];
        for (rung, steps) in [32, 64, 128].into_iter().enumerate() {
            let result = simulate(affine, steps);
            errors[rung] = result
                .iter()
                .zip(reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
        }
        assert!(
            errors[0] > 1e-8,
            "manufactured nonlinear axis must be live: {errors:?}"
        );
        assert!(
            errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1e-4, "{errors:?}");
    }
}

#[test]
fn fixed_geometry_recovers_exact_memory_heat_and_coherent_duplicate_sources() {
    let fixed = PortHamiltonian::new(
        2,
        1,
        vec![0.0; 4],
        vec![0.0; 4],
        vec![0.0; 2],
        Box::new(QuadraticStorage::new(vec![0.0; 4], 2).unwrap()),
    )
    .unwrap()
    .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain: Rc::new(Geometry {
            dimensions: 2,
            linear: [1.0, 0.0],
            hessian: [1.0, 0.0, 0.0],
            calls: Rc::new((Cell::new(0), Cell::new(0))),
        }),
        stiffness: 4.0,
        relaxation_time_s: 0.02,
    }])
    .unwrap();
    let h = 0.4;
    let initial = [0.3, 0.0, -0.1];
    let target = 2.0 * (0.3 + 0.5 * 0.3 * 0.3);
    let mut prepared =
        PreparedNonlinearTones::new_owned(fixed, &[2.1], h, split(), dense()).unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 1]);
    let result = prepared
        .step_into(&initial, &[0.7], &[-0.2], &mut next, &mut mean)
        .unwrap();
    assert_eq!(next[0], initial[0]);
    assert_eq!(next[1], initial[1]);
    close(
        next[2],
        target + (initial[2] - target) * (-h / 0.02).exp(),
        3e-13,
    );
    close(
        result.dissipated,
        0.5 * (target - initial[2]).powi(2) * (1.0 - (-2.0 * h / 0.02).exp()),
        3e-13,
    );
    assert_eq!(result.supplied, 0.0);
    let mut one = PreparedNonlinearTones::new_owned(
        scalar(Rc::new((Cell::new(0), Cell::new(0)))),
        &[2.1],
        0.03,
        split(),
        dense(),
    )
    .unwrap();
    let mut duplicates = PreparedNonlinearTones::new_owned(
        scalar(Rc::new((Cell::new(0), Cell::new(0)))),
        &[2.1, 2.1, 2.1],
        0.03,
        split(),
        dense(),
    )
    .unwrap();
    let (mut one_next, mut dupe_next, mut one_y, mut dupe_y) =
        ([0.0; 3], [0.0; 3], [0.0; 1], [0.0; 1]);
    let a = one
        .step_into(&initial, &[1.0], &[-0.2], &mut one_next, &mut one_y)
        .unwrap();
    let b = duplicates
        .step_into(
            &initial,
            &[1e16, 1.0, -1e16],
            &[0.3, -0.2, -0.3],
            &mut dupe_next,
            &mut dupe_y,
        )
        .unwrap();
    assert_eq!(one_next, dupe_next);
    assert_eq!(one_y, dupe_y);
    assert_eq!(a.supplied, b.supplied);
    assert_eq!(a.dissipated, b.dissipated);
}

#[derive(Debug)]
struct HintedGeometry(Geometry);
impl RelaxationStrainMap for HintedGeometry {
    fn state_dimension(&self) -> usize {
        self.0.state_dimension()
    }
    fn strain(&self, x: &[f64]) -> f64 {
        self.0.strain(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.0.gradient(x, out);
    }
    fn hessian_vector(&self, x: &[f64], direction: &[f64], out: &mut [f64]) {
        self.0.hessian_vector(x, direction, out);
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
    fn zero_constant_linear_strain(&self) -> Option<Vec<f64>> {
        let mut row = vec![0.0; self.0.dimensions];
        row[..2].copy_from_slice(&self.0.linear);
        Some(row)
    }
}

#[test]
fn declared_linear_row_preserves_the_complete_maxwell_dyad_once_under_coherent_tones() {
    let frequencies = [0.0, 2.1, 7.3];
    let (c0, s0) = ([0.2, 0.8, -0.4], [7.0, -0.1, 0.2]);
    let (c1, s1) = ([-0.3, -0.2, 0.7], [-11.0, 0.4, -0.3]);
    let initial = [0.2, -0.3, 0.1];
    // The row hint captures the entire axial/history block. Zero curvature
    // leaves no nonlinear effort to duplicate, even across 35 memory times.
    for affine in [false, true] {
        let mapped = LinearMechanicalSystem::new(&[2.0], &[0.4], &[3.0], &[1.0], 1, 1, vec![])
            .unwrap()
            .into_system()
            .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
                strain: Rc::new(HintedGeometry(Geometry {
                    dimensions: 2,
                    linear: [1.3, 0.0],
                    hessian: [0.0; 3],
                    calls: Rc::new((Cell::new(0), Cell::new(0))),
                })),
                stiffness: 4.0,
                relaxation_time_s: 0.002,
            }])
            .unwrap();
        let linear = LinearMechanicalSystem::new(
            &[2.0],
            &[0.4],
            &[3.0],
            &[1.0],
            1,
            1,
            vec![RelaxationBranch {
                projection: vec![1.3],
                stiffness: 4.0,
                relaxation_time_s: 0.002,
            }],
        )
        .unwrap()
        .into_system();
        let mapped_q = mapped.storage().quadratic_part().unwrap();
        let linear_q = linear.storage().quadratic_part().unwrap();
        assert_eq!(mapped_q.matrix, linear_q.matrix);
        assert!(!mapped_q.exact); // The declaration supplies a row, not an exact-energy claim.
        assert!(mapped.storage().discrete_gradient_rule().is_some());
        let mut mapped = if affine {
            PreparedNonlinearTones::new_owned_affine(mapped, &frequencies, 0.07, split(), dense())
        } else {
            PreparedNonlinearTones::new_owned(mapped, &frequencies, 0.07, split(), dense())
        }
        .unwrap();
        let mut exact = if affine {
            PreparedExactTones::new_owned_affine(linear, &frequencies, 0.07, split(), dense())
        } else {
            PreparedExactTones::new_owned(linear, &frequencies, 0.07, split(), dense())
        }
        .unwrap();
        let (mut actual, mut expected, mut actual_y, mut expected_y) =
            ([0.0; 3], [0.0; 3], [0.0; 1], [0.0; 1]);
        let a = if affine {
            mapped.step_affine_into(&initial, &c0, &s0, &c1, &s1, &mut actual, &mut actual_y)
        } else {
            mapped.step_into(&initial, &c0, &s0, &mut actual, &mut actual_y)
        }
        .unwrap();
        let b = if affine {
            exact.step_affine_into(&initial, &c0, &s0, &c1, &s1, &mut expected, &mut expected_y)
        } else {
            exact.step_into(&initial, &c0, &s0, &mut expected, &mut expected_y)
        }
        .unwrap();
        for (a, b) in actual
            .iter()
            .chain(&actual_y)
            .zip(expected.iter().chain(&expected_y))
        {
            close(*a, *b, 8e-12);
        }
        close(a.supplied, b.supplied, 8e-12);
        close(a.dissipated, b.dissipated, 8e-12);
        close(a.delta_h, b.delta_h, 8e-12);
        close(a.balance_residual(), 0.0, 8e-12);
        assert!(a.dissipated > 0.0);
    }
}

#[derive(Debug)]
struct InvalidGeometry {
    incomplete: bool,
}
impl RelaxationStrainMap for InvalidGeometry {
    fn state_dimension(&self) -> usize {
        2
    }
    fn strain(&self, x: &[f64]) -> f64 {
        x[0] * x[1]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[x[1], x[0]]);
    }
    fn hessian_vector(&self, _: &[f64], v: &[f64], out: &mut [f64]) {
        out[0] = v[1];
        if !self.incomplete {
            out[1] = 2.0 * v[0];
        }
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
}

#[derive(Debug)]
struct PositiveDomain;
impl RelaxationStrainMap for PositiveDomain {
    fn state_dimension(&self) -> usize {
        2
    }
    fn strain(&self, x: &[f64]) -> f64 {
        if x[0] > 0.0 {
            x[0] + 0.5 * x[0] * x[0]
        } else {
            f64::NAN
        }
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if x[0] > 0.0 {
            out.copy_from_slice(&[1.0 + x[0], 0.0]);
        } else {
            out.fill(f64::NAN);
        }
    }
    fn hessian_vector(&self, x: &[f64], v: &[f64], out: &mut [f64]) {
        if x[0] > 0.0 {
            out.copy_from_slice(&[v[0], 0.0]);
        } else {
            out.fill(f64::NAN);
        }
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
}

#[test]
fn unavailable_images_retain_source_domains_and_cancelled_trials_are_atomic() {
    for incomplete in [false, true] {
        let base = LinearMechanicalSystem::new(&[2.0], &[0.4], &[3.0], &[1.0], 1, 1, vec![])
            .unwrap()
            .into_system();
        let uncached = base
            .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
                strain: Rc::new(InvalidGeometry { incomplete }),
                stiffness: 4.0,
                relaxation_time_s: 0.2,
            }])
            .unwrap();
        assert!(uncached.storage().discrete_gradient_rule().is_none());
    }
    // A polynomial declaration does not widen the physical source domain.
    // Missing analytic data at zero declines the cache, not the material.
    let domain = LinearMechanicalSystem::new(&[2.0], &[0.4], &[3.0], &[1.0], 1, 1, vec![])
        .unwrap()
        .into_system()
        .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
            strain: Rc::new(PositiveDomain),
            stiffness: 4.0,
            relaxation_time_s: 0.2,
        }])
        .unwrap();
    assert!(domain.storage().discrete_gradient_rule().is_none());
    let effort = domain.effort(&[0.3, 0.2, 0.1]);
    close(effort[0], 0.9 + 2.6 * (0.69 - 0.1), 2e-15);
    let mut workspace = StepWorkspace::new(&domain).unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 1]);
    let record = workspace
        .step_into(
            &domain,
            &[0.3, 0.2, 0.1],
            &[0.7],
            0.01,
            &mut next,
            &mut mean,
        )
        .unwrap();
    assert!(next[0] > 0.0 && record.dissipated > 0.0);
    let mut prepared = PreparedNonlinearTones::new_owned(
        scalar(Rc::new((Cell::new(0), Cell::new(0)))),
        &[2.1],
        0.03,
        split(),
        dense(),
    )
    .unwrap();
    let initial = [0.3, 0.2, -0.1];
    let (mut next, mut mean) = ([17.0; 3], [19.0; 1]);
    assert!(
        prepared
            .step_into_controlled(&initial, &[0.7], &[-0.2], &mut next, &mut mean, || true)
            .is_err()
    );
    assert_eq!(next, [17.0; 3]);
    assert_eq!(mean, [19.0; 1]);
    let original = initial;
    let accepted = prepared
        .step_into(&initial, &[0.7], &[-0.2], &mut next, &mut mean)
        .unwrap();
    assert_eq!(initial, original);
    assert!(accepted.dissipated > 0.0);
    let saved = (next, mean, accepted);
    prepared
        .step_into(&initial, &[0.7], &[-0.2], &mut next, &mut mean)
        .unwrap();
    assert_eq!(next, saved.0);
    assert_eq!(mean, saved.1);
    // Linear relaxation around a native constitutive child retains the child
    // rule and its material rows instead of discarding that capability.
    let nested = scalar(Rc::new((Cell::new(0), Cell::new(0))))
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![0.5, 0.0, 0.0],
            stiffness: 2.0,
            relaxation_time_s: 0.1,
        }])
        .unwrap();
    let rule = nested.storage().discrete_gradient_rule().unwrap();
    let mut action = [0.0; 4];
    assert!(rule.end_jacobian_action(
        &[0.2, 0.3, 0.1, -0.2],
        &[0.3, 0.4, -0.1, 0.2],
        &[0.3, -0.2, 0.1, 0.4],
        &mut action
    ));
    assert!(action.iter().all(|v| v.is_finite()));
    let linear_first = LinearMechanicalSystem::new(
        &[2.0],
        &[0.4],
        &[3.0],
        &[1.0],
        1,
        1,
        vec![RelaxationBranch {
            projection: vec![0.5],
            stiffness: 2.0,
            relaxation_time_s: 0.1,
        }],
    )
    .unwrap()
    .into_system()
    .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain: Rc::new(Geometry {
            dimensions: 3,
            linear: [1.0, 0.0],
            hessian: [1.0, 0.0, 0.0],
            calls: Rc::new((Cell::new(0), Cell::new(0))),
        }),
        stiffness: 4.0,
        relaxation_time_s: 0.2,
    }])
    .unwrap();
    assert!(
        linear_first
            .storage()
            .discrete_gradient_rule()
            .unwrap()
            .end_jacobian_action(
                &[0.2, 0.3, -0.2, 0.1],
                &[0.3, 0.4, 0.2, -0.1],
                &[0.3, -0.2, 0.4, 0.1],
                &mut action,
            )
    );
}
