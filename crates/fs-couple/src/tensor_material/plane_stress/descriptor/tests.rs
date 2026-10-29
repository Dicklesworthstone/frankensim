use super::*;
use fs_material::visco_tensor::TensorPronyBranch;
use fs_material::{IsotropicElastic, Tangent6};

fn projection_budget() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 16,
        max_scalar_terms: 1_000_000,
    }
}
fn body_budget() -> TensorMaterialBudget {
    TensorMaterialBudget {
        max_state_coordinates: 128,
        max_material_points: 4,
    }
}
// An orthonormal basis with no coordinate-aligned vector. The exact physical
// source declares the number of independent viscous normal channels below.
const FRAME: [[f64; 3]; 3] = [
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
    [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
];
const MASS: [f64; 4] = [1.4, 0.3, 0.3, 0.9];

fn point(rank: usize) -> DescriptorPlaneStressPoint {
    let mut b: [[f64; 6]; 6] = [[0.0; 6]; 6];
    b[0] = [2.0, -0.8, 0.4, 0.7, 0.3, -0.2];
    b[1] = [-0.5, 1.8, 0.7, -0.3, 0.1, 0.6];
    let branches = vec![
        TensorPronyBranch::from_gram_rows(b, 0.04).unwrap(),
        TensorPronyBranch::from_gram_rows(
            core::array::from_fn(|i| core::array::from_fn(|j| 0.7 * b[j][i])),
            0.13,
        )
        .unwrap(),
    ];
    let base = TensorProny::from_tensor_shear(
        IsotropicElastic::new(120.0, 0.24, 0.1).unwrap().stiffness(),
        branches,
        0.1,
    )
    .unwrap();
    let mut factors: Tangent6 = [[0.0; 6]; 6];
    factors[3][0] = 1.2;
    factors[4][1] = 1.4;
    factors[5][3] = 0.9;
    for i in 0..rank {
        factors[i][0] = 0.2 * (i + 1) as f64;
        factors[i][1] = -0.3;
        for k in 0..3 {
            factors[i][[2, 4, 5][k]] = (1.6 + 0.4 * i as f64) * FRAME[i][k]
                / if k == 0 {
                    1.0
                } else {
                    core::f64::consts::SQRT_2
                };
        }
    }
    let viscosity = core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|k| factors[k][i] * factors[k][j]).sum::<f64>())
    });
    let material = TensorProny::new(
        *base.c_eq(),
        Some(viscosity),
        base.branches().to_vec(),
        (0.0, f64::MAX),
    )
    .unwrap();
    DescriptorPlaneStressPoint::new(
        &material,
        &[[1.0, 0.0, 0.2], [0.0, 1.0, -0.3]],
        0.2,
        projection_budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn body(rank: usize) -> DescriptorPlaneStressBody {
    DescriptorPlaneStressBody::from_mass_matrix(
        2,
        MASS.to_vec(),
        vec![point(rank)],
        body_budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn physical_state(h: &DescriptorPlaneStressHistory<'_>) -> Vec<f64> {
    let mut x = h.state()[..h.body.differential_dimension()].to_vec();
    x[2..4].copy_from_slice(&h.velocity_m_s());
    x
}

// Independent continuous equations in physical velocities and all 3 normal
// strains. Only the SOURCE's explicitly constructed nullspace appears here.
// The production descriptor never reads FRAME or the source rank.
struct Reference {
    k: Vec<f64>,
    d: Vec<f64>,
    a: Vec<Vec<f64>>,
    rates: Vec<f64>,
    normal: fs_la::factor::Lu,
    rank: usize,
    instantaneous: Vec<f64>,
}
impl Reference {
    fn new(point: &DescriptorPlaneStressPoint, rank: usize) -> Self {
        let source = point.projection();
        let k = source.equilibrium_stiffness().to_vec();
        let d = source.viscosity_matrix().to_vec();
        let a: Vec<_> = source
            .branches()
            .iter()
            .flat_map(|b| b.rows().iter().cloned())
            .collect();
        let rates = source
            .branches()
            .iter()
            .flat_map(|b| [b.relaxation_time_s().recip(); 6])
            .collect();
        let instantaneous = (0..25)
            .map(|ij| k[ij] + a.iter().map(|row| row[ij / 5] * row[ij % 5]).sum::<f64>())
            .collect::<Vec<_>>();
        let mut matrix = vec![0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                matrix[3 * i + j] = (0..3)
                    .map(|r| {
                        FRAME[i][r]
                            * if i < rank {
                                d[(2 + r) * 5 + 2 + j]
                            } else {
                                instantaneous[(2 + r) * 5 + 2 + j]
                            }
                    })
                    .sum();
            }
        }
        let normal = fs_la::factor::lu(&matrix, 3).unwrap();
        Self {
            k,
            d,
            a,
            rates,
            normal,
            rank,
            instantaneous,
        }
    }
    fn rhs(&self, x: &[f64], force: [f64; 2]) -> Vec<f64> {
        let strain = [x[0], x[1], x[16], x[17], x[18]];
        let mut stress = [0.0; 5];
        for i in 0..5 {
            stress[i] = (0..5).map(|j| self.k[i * 5 + j] * strain[j]).sum::<f64>()
                + self
                    .a
                    .iter()
                    .enumerate()
                    .map(|(j, row)| row[i] * x[4 + j])
                    .sum::<f64>();
        }
        let mut normal_rate = vec![0.0; 3];
        for i in 0..3 {
            normal_rate[i] = (0..3)
                .map(|r| {
                    FRAME[i][r]
                        * if i < self.rank {
                            -stress[2 + r]
                                - self.d[(2 + r) * 5] * x[2]
                                - self.d[(2 + r) * 5 + 1] * x[3]
                        } else {
                            -self.instantaneous[(2 + r) * 5] * x[2]
                                - self.instantaneous[(2 + r) * 5 + 1] * x[3]
                                + self
                                    .a
                                    .iter()
                                    .enumerate()
                                    .map(|(j, row)| row[2 + r] * self.rates[j] * x[4 + j])
                                    .sum::<f64>()
                        }
                })
                .sum();
        }
        self.normal.solve(&mut normal_rate);
        let rates = [x[2], x[3], normal_rate[0], normal_rate[1], normal_rate[2]];
        let mut dx = vec![0.0; 21];
        dx[..2].copy_from_slice(&x[2..4]);
        for i in 0..5 {
            let loss = (0..5).map(|j| self.d[i * 5 + j] * rates[j]).sum::<f64>();
            stress[i] += loss;
            dx[20] += rates[i] * loss;
        }
        let determinant = MASS[0] * MASS[3] - MASS[1] * MASS[2];
        dx[2] = (MASS[3] * (force[0] - stress[0]) - MASS[1] * (force[1] - stress[1])) / determinant;
        dx[3] =
            (-MASS[2] * (force[0] - stress[0]) + MASS[0] * (force[1] - stress[1])) / determinant;
        for i in 0..12 {
            dx[4 + i] = self.a[i].iter().zip(rates).map(|(a, v)| a * v).sum::<f64>()
                - self.rates[i] * x[4 + i];
            dx[20] += self.rates[i] * x[4 + i] * x[4 + i];
        }
        dx[16..19].copy_from_slice(&normal_rate);
        dx[19] = force[0] * x[2] + force[1] * x[3];
        for f in &stress[2..] {
            assert!(f.abs() < 1e-11, "source normal traction {f}");
        }
        dx
    }
}
fn rk4(x: &mut [f64], h: f64, rhs: impl Fn(&[f64]) -> Vec<f64>) {
    let a = rhs(x);
    let stage = |k: &[f64], scale: f64| {
        x.iter()
            .zip(k)
            .map(|(x, k)| x + scale * k)
            .collect::<Vec<_>>()
    };
    let b = rhs(&stage(&a, 0.5 * h));
    let c = rhs(&stage(&b, 0.5 * h));
    let d = rhs(&stage(&c, h));
    for i in 0..x.len() {
        x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
    }
}
fn error(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn oblique_normal_constraints_retain_full_mass_memory_and_physical_second_order_motion() {
    let gate = CancelGate::new_clock_free();
    let force = [0.03, -0.02];
    let duration = 0.04;
    for rank in 0..=3 {
        let body = body(rank);
        assert_eq!(body.state_dimension(), 22);
        assert_eq!(body.differential_dimension(), 19);
        let source = Reference::new(&body.points[0], rank);
        let start = body
            .history_relaxed_at_rest(&[0.001, -0.002], &gate)
            .unwrap();
        let mut reference = physical_state(&start);
        reference.extend([0.0, 0.0]);
        for _ in 0..4096 {
            rk4(&mut reference, duration / 4096.0, |x| source.rhs(x, force));
        }
        let mut errors = Vec::new();
        let mut heat_errors = Vec::new();
        let mut work_errors = Vec::new();
        for steps in [20, 40, 80] {
            let mut history = body
                .history_relaxed_at_rest(&[0.001, -0.002], &gate)
                .unwrap();
            for _ in 0..steps {
                let record = history
                    .advance(&force, duration / f64::from(steps), &gate)
                    .unwrap();
                assert!(record.dissipated >= 0.0);
                for stress in history.plane_stress_normal_stress_pa(0).unwrap() {
                    assert!(stress.abs() < 2e-10);
                }
            }
            assert!(history.balance_residual_j().abs() < 1e-12);
            errors.push(error(&physical_state(&history), &reference[..19]));
            work_errors.push((history.supplied_work_j() - reference[19]).abs());
            heat_errors.push((history.dissipated_heat_j() - reference[20]).abs());
        }
        eprintln!("rank {rank}: state {errors:?}; work {work_errors:?}; heat {heat_errors:?}");
        for values in [&errors, &work_errors, &heat_errors] {
            assert!(
                values[0] / values[1] > 3.8 && values[1] / values[2] > 3.8,
                "rank {rank}: {values:?}"
            );
        }
        assert!(errors[2] < 2e-7);
        assert!(work_errors[2] < 1e-10);
        assert!(heat_errors[2] < 1e-9);
    }
}

#[test]
fn complete_restart_preserves_constraint_state_and_refusals_publish_nothing() {
    let gate = CancelGate::new_clock_free();
    let body = body(1);
    let mut history = body
        .history_relaxed_at_rest(&[0.001, -0.002], &gate)
        .unwrap();
    for _ in 0..5 {
        history.advance(&[0.03, -0.02], 0.001, &gate).unwrap();
    }
    let accepted = history.state().to_vec();
    let work = history.supplied_work_j();
    let heat = history.dissipated_heat_j();
    let mut restored = body.history_from_state(accepted.clone(), &gate).unwrap();
    let mut polls = 0;
    assert!(matches!(
        history.advance_controlled(&[0.03, -0.02], 0.001, || {
            polls += 1;
            polls == 3
        }),
        Err(TensorMaterialError::Cancelled)
    ));
    assert_eq!(polls, 3);
    assert_eq!(history.state(), accepted);
    assert_eq!(
        (history.supplied_work_j(), history.dissipated_heat_j()),
        (work, heat)
    );
    assert!(history.advance(&[f64::NAN, 0.0], 0.001, &gate).is_err());
    assert_eq!(history.state(), accepted);
    let a = history.advance(&[-0.02, 0.01], 0.0005, &gate).unwrap();
    let b = restored.advance(&[-0.02, 0.01], 0.0005, &gate).unwrap();
    assert_eq!(a.x, b.x);
    assert_eq!(a.dissipated, b.dissipated);
    assert_eq!(a.supplied, b.supplied);
    let mut inconsistent = accepted.clone();
    inconsistent[16] += 0.001;
    assert!(body.history_from_state(inconsistent, &gate).is_err());
    let mut moving = accepted;
    moving[2] += 1.0;
    assert!(body.history_from_state(moving, &gate).is_err());
}

#[test]
fn full_coordinate_capacity_and_original_material_domain_are_explicit() {
    let gate = CancelGate::new_clock_free();
    let point = point(2);
    assert!(
        DescriptorPlaneStressBody::from_mass_matrix(
            2,
            MASS.to_vec(),
            vec![point.clone()],
            TensorMaterialBudget {
                max_state_coordinates: 21,
                ..body_budget()
            },
            &gate
        )
        .is_err()
    );
    let body = DescriptorPlaneStressBody::from_mass_matrix(
        2,
        MASS.to_vec(),
        vec![point],
        TensorMaterialBudget {
            max_state_coordinates: 22,
            ..body_budget()
        },
        &gate,
    )
    .unwrap();
    assert_eq!(body.state_dimension(), 22);
    assert_eq!(body.mass_matrix_kg(), MASS);
    let mut x = body
        .history_relaxed_at_rest(&[0.0, 0.0], &gate)
        .unwrap()
        .state()
        .to_vec();
    x[0] = f64::NAN;
    assert!(body.history_from_state(x, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        body.history_relaxed_at_rest(&[0.0, 0.0], &cancelled)
            .is_err()
    );
    let bounded = TensorProny::from_tensor_shear(
        IsotropicElastic::new(120.0, 0.24, 0.1).unwrap().stiffness(),
        vec![],
        0.001,
    )
    .unwrap();
    let bounded_point = DescriptorPlaneStressPoint::new(
        &bounded,
        &[[1.0, 0.0, 0.0]],
        0.2,
        projection_budget(),
        &gate,
    )
    .unwrap();
    let bounded_body =
        DescriptorPlaneStressBody::new(vec![1.0], vec![bounded_point], body_budget(), &gate)
            .unwrap();
    assert!(bounded_body.history_relaxed_at_rest(&[0.1], &gate).is_err());
}
