//! G1/G3: original Maxwell equations, full fields, physical work and scaling.
use super::super::*;

fn matrix(n: usize, m: usize, values: &[f64]) -> Csr {
    let mut a = Coo::new(n, m);
    for i in 0..n {
        for j in 0..m {
            if values[i * m + j] != 0.0 {
                a.push(i, j, values[i * m + j]);
            }
        }
    }
    a.assemble()
}
fn scalar(v: f64) -> Csr {
    matrix(1, 1, &[v])
}
fn simple() -> SparseForms {
    SparseForms {
        mass: [scalar(2.0), scalar(0.3)],
        stiffness: [scalar(4.0), scalar(0.0)],
        damping: [scalar(0.05), scalar(0.02)],
        coupling: scalar(0.7),
    }
}
fn prepare(forms: SparseForms, dt: f64, sources: &[BranchSource<'_>]) -> SparseMidpoint {
    SparseMidpoint::new_relaxing(
        forms,
        dt,
        SparseStructuralTimeLimits::default(),
        sources,
        SparseStructuralRelaxationLimits::default(),
        &|| false,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}: {:e}",
        (a - b).abs()
    );
}
fn state_values(state: &SparseStructuralTimeState<'_>) -> Vec<f64> {
    state
        .position
        .iter()
        .chain(&state.velocity)
        .chain(&state.stress)
        .copied()
        .collect()
}

#[test]
fn dark_original_histories_keep_storage_and_positive_heat_without_spatial_motion() {
    for tau in [0.00001, 0.6, 200.0] {
        let dt = 0.02;
        let k = 3.7;
        let owner = prepare(
            simple(),
            dt,
            &[BranchSource {
                strain: &[0.0, 0.0],
                stiffness: k,
                relaxation_time_s: tau,
            }],
        );
        let mut state = owner
            .state_with_stresses(0.0, &[0.0], &[0.0], &[0.0], &[0.0], &[1.2], &|| false)
            .unwrap();
        near(
            owner.energy(&state, &|| false).unwrap().relaxation_j,
            1.2 * 1.2 / (2.0 * k),
            2e-16,
        );
        let ratio = (tau - 0.5 * dt) / (tau + 0.5 * dt);
        let mut sigma = 1.2;
        let mut heat = 0.0;
        let initial_energy = owner.energy(&state, &|| false).unwrap().total_j();
        for _ in 0..25 {
            let step = owner.step(&state, &[0.0], &[0.0], &|| false).unwrap();
            let expected_next = ratio * sigma;
            let midpoint = 0.5 * (sigma + expected_next);
            let expected_heat = dt * midpoint * midpoint / (k * tau);
            near(step.state.branch_stresses()[0], expected_next, 3e-14);
            near(step.memory[0].dissipation_j, expected_heat, 3e-14);
            assert!(step.memory[0].dissipation_j >= 0.0);
            assert_eq!(step.state.position, [0.0, 0.0]);
            assert_eq!(step.state.velocity, [0.0, 0.0]);
            assert_eq!(step.work.structural_to_relaxation_work_j, 0.0);
            assert_eq!(step.work.acoustic_to_relaxation_work_j, 0.0);
            heat += step.memory[0].dissipation_j;
            sigma = expected_next;
            state = step.state;
        }
        near(
            owner.energy(&state, &|| false).unwrap().total_j() + heat,
            initial_energy,
            2e-13,
        );
    }
}

#[test]
fn driven_original_solid_and_wall_memories_converge_to_independent_continuous_fields() {
    // Entire exact trajectory prescribed independently of any discretized
    // matrix or propagator: q=[0.3,-0.2] exp(lambda*t), sigma=k*b*lambda*q/(lambda+1/tau).
    let lambda = -0.3;
    let duration = 0.7;
    let q = [0.3, -0.2];
    let bs = [0.8, -0.6];
    let ks = [1.3, 0.9];
    let tau = [0.7, 1.9];
    let stresses = [
        ks[0] * bs[0] * lambda * q[0] / (lambda + 1.0 / tau[0]),
        ks[1] * bs[1] * lambda * q[1] / (lambda + 1.0 / tau[1]),
    ];
    let force = (2.0 * lambda * lambda + 0.05 * lambda + 4.0) * q[0] - 0.7 * lambda * q[1]
        + bs[0] * stresses[0];
    let source =
        (0.3 * lambda * lambda + 0.02 * lambda) * q[1] + 0.7 * lambda * q[0] + bs[1] * stresses[1];
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let dt = duration / count as f64;
        let owner = prepare(
            simple(),
            dt,
            &[
                BranchSource {
                    strain: &[bs[0], 0.0],
                    stiffness: ks[0],
                    relaxation_time_s: tau[0],
                },
                BranchSource {
                    strain: &[0.0, bs[1]],
                    stiffness: ks[1],
                    relaxation_time_s: tau[1],
                },
            ],
        );
        let mut state = owner
            .state_with_stresses(
                0.0,
                &q[..1],
                &[lambda * q[0]],
                &q[1..],
                &[lambda * q[1]],
                &stresses,
                &|| false,
            )
            .unwrap();
        let mut heat = 0.0;
        for step in 0..count {
            let t = (step as f64 + 0.5) * dt;
            let scale = (lambda * t).exp();
            let next = owner
                .step(&state, &[force * scale], &[source * scale], &|| false)
                .unwrap();
            for m in &next.memory {
                near(m.balance_residual_j(), 0.0, 2e-13);
                assert!(m.dissipation_j > 0.0);
            }
            assert!(next.work.structural_to_relaxation_work_j > 0.0);
            assert!(next.work.acoustic_to_relaxation_work_j > 0.0);
            near(next.work.structural_balance_residual_j(), 0.0, 2e-13);
            near(next.work.acoustic_balance_residual_j(), 0.0, 2e-13);
            near(next.work.relaxation_balance_residual_j(), 0.0, 2e-13);
            near(next.work.balance_residual_j(), 0.0, 2e-13);
            heat += next.work.relaxation_dissipation_j;
            state = next.state;
        }
        let scale = (lambda * duration).exp();
        let expected = [
            q[0] * scale,
            q[1] * scale,
            lambda * q[0] * scale,
            lambda * q[1] * scale,
            stresses[0] * scale,
            stresses[1] * scale,
        ];
        let error = state_values(&state)
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>();
        errors.push(error);
        let exact_heat = (stresses[0] * stresses[0] / (ks[0] * tau[0])
            + stresses[1] * stresses[1] / (ks[1] * tau[1]))
            * ((2.0 * lambda * duration).exp() - 1.0)
            / (2.0 * lambda);
        if count == 64 {
            near(heat, exact_heat, 1e-6);
        }
    }
    eprintln!("sparse original Maxwell continuous refinement: {errors:?}");
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
    assert!(errors[2] < 1e-5);
}

fn small() -> SparseForms {
    SparseForms {
        mass: [
            matrix(2, 2, &[2.0, 0.3, 0.3, 1.1]),
            matrix(3, 3, &[0.4, 0.03, 0.0, 0.03, 0.3, 0.02, 0.0, 0.02, 0.5]),
        ],
        stiffness: [
            matrix(2, 2, &[6.0, -1.0, -1.0, 3.0]),
            matrix(3, 3, &[2.0, -1.0, -1.0, -1.0, 1.5, -0.5, -1.0, -0.5, 1.5]),
        ],
        damping: [
            matrix(2, 2, &[0.2, 0.03, 0.03, 0.1]),
            matrix(3, 3, &[0.05, 0.01, 0.0, 0.01, 0.04, 0.0, 0.0, 0.0, 0.02]),
        ],
        coupling: matrix(3, 2, &[0.2, 0.07, -0.04, 0.3, 0.15, -0.08]),
    }
}
fn directions() -> [[f64; 5]; 5] {
    [
        [0.7, -0.2, 0.0, 0.0, 0.0],
        [-0.1, 0.9, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.3, -0.5, 0.8],
        [0.0, 0.0, 0.4, 0.7, -0.2],
        [0.0; 5],
    ]
}
fn gaussian(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for i in 0..n {
        let pivot = (i..n)
            .max_by(|&j, &k| a[j * n + i].abs().total_cmp(&a[k * n + i].abs()))
            .unwrap();
        for k in i..n {
            a.swap(i * n + k, pivot * n + k);
        }
        b.swap(i, pivot);
        for j in i + 1..n {
            let s = a[j * n + i] / a[i * n + i];
            for k in i..n {
                a[j * n + k] -= s * a[i * n + k];
            }
            b[j] -= s * b[i];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            b[i] -= a[i * n + j] * b[j];
        }
        b[i] /= a[i * n + i];
    }
    b
}

#[test]
fn sparse_augmented_memories_match_independent_original_endpoint_equations() {
    let h = 0.07;
    let directions = directions();
    let branches = directions
        .iter()
        .enumerate()
        .map(|(i, b)| BranchSource {
            strain: b,
            stiffness: 0.9 + 0.2 * i as f64,
            relaxation_time_s: 0.1 + 0.13 * i as f64,
        })
        .collect::<Vec<_>>();
    let owner = prepare(small(), h, &branches);
    let q = [0.2, -0.1, 0.01, 0.04, -0.02];
    let v = [0.1, 0.2, 0.3, -0.1, 0.07];
    let f = [1.2, -0.2, 0.4, 0.02, -0.03];
    let stress = [0.2, -0.3, 0.4, 0.1, -0.5];
    let initial = owner
        .state_with_stresses(0.0, &q[..2], &v[..2], &q[2..], &v[2..], &stress, &|| false)
        .unwrap();
    let step = owner.step(&initial, &f[..2], &f[2..], &|| false).unwrap();
    // Independent nonsymmetric ORIGINAL equations, in endpoint velocity and
    // physical stress. No signed acoustic row, sqrt(k) transform or LDLT.
    let n = 5;
    let dim = n + branches.len();
    let mut a = vec![0.0; dim * dim];
    let mut rhs = vec![0.0; dim];
    for i in 0..n {
        let d = usize::from(i >= 2);
        let off = if d == 0 { 0 } else { 2 };
        let end = if d == 0 { 2 } else { 5 };
        let row = i - off;
        rhs[i] = f[i];
        for j in off..end {
            let col = j - off;
            let m = owner.forms.mass[d].get(row, col);
            let c = owner.forms.damping[d].get(row, col);
            let k = owner.forms.stiffness[d].get(row, col);
            a[i * dim + j] = m / h + 0.5 * c + 0.25 * h * k;
            rhs[i] += (m / h - 0.5 * c - 0.25 * h * k) * v[j] - k * q[j];
        }
    }
    for i in 0..3 {
        for j in 0..2 {
            let b = owner.forms.coupling.get(i, j);
            a[j * dim + 2 + i] = -0.5 * b;
            rhs[j] += 0.5 * b * v[2 + i];
            a[(2 + i) * dim + j] = 0.5 * b;
            rhs[2 + i] -= 0.5 * b * v[j];
        }
    }
    for (r, b) in branches.iter().enumerate() {
        let index = n + r;
        rhs[index] = (1.0 - 0.5 * h / b.relaxation_time_s) * stress[r];
        a[index * dim + index] = 1.0 + 0.5 * h / b.relaxation_time_s;
        for i in 0..n {
            a[i * dim + index] = 0.5 * b.strain[i];
            rhs[i] -= 0.5 * b.strain[i] * stress[r];
            a[index * dim + i] = -0.5 * h * b.stiffness * b.strain[i];
            rhs[index] += 0.5 * h * b.stiffness * b.strain[i] * v[i];
        }
    }
    let expected = gaussian(a, rhs);
    for i in 0..n {
        near(step.state.velocity[i], expected[i], 2e-13);
        near(
            step.state.position[i],
            q[i] + 0.5 * h * (v[i] + expected[i]),
            2e-13,
        );
    }
    for (j, &sigma) in step.state.branch_stresses().iter().enumerate() {
        near(sigma, expected[n + j], 2e-13);
    }
    assert!(step.work.structural_to_relaxation_work_j.abs() > 1e-6);
    assert!(step.work.acoustic_to_relaxation_work_j.abs() > 1e-6);
    assert_eq!(step.memory[4].structural_work_j, 0.0);
    assert_eq!(step.memory[4].acoustic_work_j, 0.0);
    assert!(step.memory[4].dissipation_j > 0.0);
}

fn chain(n: usize, diagonal: f64, edge: f64) -> Csr {
    let mut a = Coo::new(n, n);
    for i in 0..n {
        a.push(i, i, diagonal);
    }
    for i in 1..n {
        for (r, c, x) in [
            (i, i, edge),
            (i - 1, i - 1, edge),
            (i, i - 1, -edge),
            (i - 1, i, -edge),
        ] {
            a.push(r, c, x);
        }
    }
    a.assemble()
}
#[test]
fn large_full_fields_retain_nonlocal_material_maps_without_dense_rank_updates() {
    let (ns, na) = (48, 64);
    let n = ns + na;
    let mut coupling = Coo::new(na, ns);
    for a in 0..na {
        coupling.push(a, a % ns, 0.07 + 0.0001 * a as f64);
        coupling.push(a, (a + 3) % ns, -0.02);
    }
    let forms = SparseForms {
        mass: [chain(ns, 1.7, 0.13), chain(na, 0.4, 0.02)],
        stiffness: [chain(ns, 0.0, 7.0), chain(na, 0.0, 1.3)],
        damping: [chain(ns, 0.03, 0.013), chain(na, 0.01, 0.0)],
        coupling: coupling.assemble(),
    };
    let solid = (0..n)
        .map(|i| if i < ns { 0.02 + 0.001 * i as f64 } else { 0.0 })
        .collect::<Vec<_>>();
    let wall = (0..n)
        .map(|i| {
            if i >= ns {
                -0.01 - 0.0001 * i as f64
            } else {
                0.0
            }
        })
        .collect::<Vec<_>>();
    // Dense rank updates alone would require 48^2+64^2 entries.
    let mut limits = SparseStructuralTimeLimits::default();
    limits.maximum_matrix_entries = 3000;
    let owner = SparseMidpoint::new_relaxing(
        forms,
        0.002,
        limits,
        &[
            BranchSource {
                strain: &solid,
                stiffness: 1.2,
                relaxation_time_s: 0.07,
            },
            BranchSource {
                strain: &wall,
                stiffness: 0.8,
                relaxation_time_s: 0.11,
            },
        ],
        SparseStructuralRelaxationLimits::default(),
        &|| false,
    )
    .unwrap();
    let mut state = owner
        .state_with_stresses(
            0.0,
            &vec![0.0; ns],
            &vec![0.01; ns],
            &vec![0.0; na],
            &vec![-0.03; na],
            &[0.1, -0.2],
            &|| false,
        )
        .unwrap();
    let initial = owner.energy(&state, &|| false).unwrap().total_j();
    let mut work = 0.0;
    let mut heat = 0.0;
    for _ in 0..30 {
        let step = owner
            .step(&state, &vec![0.03; ns], &vec![0.002; na], &|| false)
            .unwrap();
        assert!(step.memory.iter().all(|m| m.dissipation_j > 0.0));
        near(step.work.interface_balance_residual_j(), 0.0, 1e-14);
        near(step.work.balance_residual_j(), 0.0, 1e-12);
        work += step.work.structural_input_work_j + step.work.acoustic_input_work_j;
        heat += step.work.structural_dissipation_j
            + step.work.wall_dissipation_j
            + step.work.relaxation_dissipation_j;
        state = step.state;
    }
    near(
        owner.energy(&state, &|| false).unwrap().total_j() + heat,
        initial + work,
        3e-12,
    );
}

#[test]
fn original_material_strain_units_preserve_fields_storage_and_work() {
    let q = [0.2, -0.1];
    let v = [0.1, 0.3];
    let h = 0.03;
    let base = prepare(
        simple(),
        h,
        &[
            BranchSource {
                strain: &[0.8, 0.0],
                stiffness: 1.3,
                relaxation_time_s: 0.7,
            },
            BranchSource {
                strain: &[0.0, -0.6],
                stiffness: 0.9,
                relaxation_time_s: 1.9,
            },
        ],
    );
    let state = base
        .state_with_stresses(
            0.0,
            &q[..1],
            &v[..1],
            &q[1..],
            &v[1..],
            &[0.2, -0.4],
            &|| false,
        )
        .unwrap();
    let expected = base.step(&state, &[0.4], &[-0.2], &|| false).unwrap();
    for scale in [1e-8, 1e8] {
        let owner = prepare(
            simple(),
            h,
            &[
                BranchSource {
                    strain: &[0.8 / scale, 0.0],
                    stiffness: 1.3 * scale * scale,
                    relaxation_time_s: 0.7,
                },
                BranchSource {
                    strain: &[0.0, -0.6 / scale],
                    stiffness: 0.9 * scale * scale,
                    relaxation_time_s: 1.9,
                },
            ],
        );
        let state = owner
            .state_with_stresses(
                0.0,
                &q[..1],
                &v[..1],
                &q[1..],
                &v[1..],
                &[0.2 * scale, -0.4 * scale],
                &|| false,
            )
            .unwrap();
        let actual = owner.step(&state, &[0.4], &[-0.2], &|| false).unwrap();
        for (a, b) in actual
            .state
            .position
            .iter()
            .chain(&actual.state.velocity)
            .zip(
                expected
                    .state
                    .position
                    .iter()
                    .chain(&expected.state.velocity),
            )
        {
            near(*a, *b, 3e-14);
        }
        for (a, b) in actual.state.stress.iter().zip(&expected.state.stress) {
            near(*a / scale, *b, 3e-14);
        }
        near(
            actual.work.final_energy.total_j(),
            expected.work.final_energy.total_j(),
            3e-14,
        );
        for (a, b) in actual.memory.iter().zip(&expected.memory) {
            near(a.dissipation_j, b.dissipation_j, 3e-14);
            near(a.structural_work_j, b.structural_work_j, 3e-14);
            near(a.acoustic_work_j, b.acoustic_work_j, 3e-14);
        }
    }
}

#[test]
fn original_wall_memory_retains_passive_frequency_response_and_the_relaxed_gauge() {
    // Isolated uniform pressure, imposed sinusoidal forcing, continuous analytic
    // wall stress sigma=k*b*P/(1/tau+i*omega). Exact initial periodic history
    // removes startup fitting; refining dt checks the physical admittance.
    let omega = 2.0;
    let (k, tau, b, amplitude): (f64, f64, f64, f64) = (1.1, 0.8, 0.7, 0.2);
    let real = k * b * amplitude * (1.0 / tau) / ((1.0 / tau).powi(2) + omega * omega);
    let imag = -k * b * amplitude * omega / ((1.0 / tau).powi(2) + omega * omega);
    let period = 2.0 * std::f64::consts::PI / omega;
    let mut errors = Vec::new();
    for count in [80, 160, 320] {
        let dt = period / count as f64;
        let mut forms = simple();
        forms.coupling = scalar(0.0);
        let owner = prepare(
            forms,
            dt,
            &[BranchSource {
                strain: &[0.0, b],
                stiffness: k,
                relaxation_time_s: tau,
            }],
        );
        let mut state = owner
            .state_with_stresses(0.0, &[0.0], &[0.0], &[3.0], &[amplitude], &[real], &|| {
                false
            })
            .unwrap();
        let mut heat = 0.0;
        for step in 0..count {
            let t = (step as f64 + 0.5) * dt;
            let p = amplitude * (omega * t).cos();
            let p_dot = -amplitude * omega * (omega * t).sin();
            let sigma = real * (omega * t).cos() - imag * (omega * t).sin();
            let drive = 0.3 * p_dot + 0.02 * p + b * sigma;
            let next = owner.step(&state, &[0.0], &[drive], &|| false).unwrap();
            heat += next.work.relaxation_dissipation_j;
            state = next.state;
        }
        let expected_heat = period * (real * real + imag * imag) / (2.0 * k * tau);
        let error = (state.pressure_pa()[0] - amplitude).abs()
            + (state.branch_stresses()[0] - real).abs()
            + (state.pressure_potential_pa_s()[0] - 3.0).abs()
            + (heat - expected_heat).abs();
        errors.push(error);
    }
    eprintln!("sparse causal wall admittance refinement: {errors:?}");
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
    assert!(errors[2] < 1e-4);
}

#[test]
fn memory_source_scope_budgets_and_cancelled_candidates_refuse_without_history_loss() {
    let invalid = [
        BranchSource {
            strain: &[0.2, 0.3],
            stiffness: 1.0,
            relaxation_time_s: 0.2,
        },
        BranchSource {
            strain: &[0.2],
            stiffness: 1.0,
            relaxation_time_s: 0.2,
        },
        BranchSource {
            strain: &[f64::NAN, 0.0],
            stiffness: 1.0,
            relaxation_time_s: 0.2,
        },
        BranchSource {
            strain: &[0.2, 0.0],
            stiffness: 0.0,
            relaxation_time_s: 0.2,
        },
        BranchSource {
            strain: &[0.2, 0.0],
            stiffness: 1.0,
            relaxation_time_s: 0.0,
        },
    ];
    for branch in &invalid {
        assert!(
            SparseMidpoint::new_relaxing(
                simple(),
                0.01,
                SparseStructuralTimeLimits::default(),
                std::slice::from_ref(branch),
                SparseStructuralRelaxationLimits::default(),
                &|| false,
            )
            .is_err()
        );
    }
    let branches = [BranchSource {
        strain: &[0.8, 0.0],
        stiffness: 1.3,
        relaxation_time_s: 0.7,
    }];
    for limits in [
        SparseStructuralRelaxationLimits {
            maximum_branches: 0,
            ..SparseStructuralRelaxationLimits::default()
        },
        SparseStructuralRelaxationLimits {
            maximum_source_entries: 1,
            ..SparseStructuralRelaxationLimits::default()
        },
    ] {
        assert!(matches!(
            SparseMidpoint::new_relaxing(
                simple(),
                0.01,
                SparseStructuralTimeLimits::default(),
                &branches,
                limits,
                &|| false,
            ),
            Err(SparseStructuralTimeError::Limit { .. })
        ));
    }
    let owner = prepare(simple(), 0.01, &branches);
    let foreign = prepare(simple(), 0.01, &branches);
    assert!(
        owner
            .state_with_stresses(0.0, &[0.0], &[0.0], &[0.0], &[0.0], &[], &|| false)
            .is_err()
    );
    assert!(
        owner
            .state_with_stresses(0.0, &[0.0], &[0.0], &[0.0], &[0.0], &[f64::NAN], &|| false)
            .is_err()
    );
    let state = owner
        .state_with_stresses(0.0, &[0.1], &[0.3], &[0.2], &[-0.1], &[0.4], &|| false)
        .unwrap();
    let before = state_values(&state);
    let expected = owner.step(&state, &[0.2], &[0.1], &|| false).unwrap();
    assert!(foreign.step(&state, &[0.2], &[0.1], &|| false).is_err());
    assert!(owner.step(&state, &[f64::NAN], &[0.1], &|| false).is_err());
    let calls = std::cell::Cell::new(0usize);
    let counted = || {
        calls.set(calls.get() + 1);
        false
    };
    owner.step(&state, &[0.2], &[0.1], &counted).unwrap();
    let total = calls.get();
    for stop in [1, total / 2, total] {
        let count = std::cell::Cell::new(0);
        let cancelled = || {
            count.set(count.get() + 1);
            count.get() >= stop
        };
        assert!(matches!(
            owner.step(&state, &[0.2], &[0.1], &cancelled),
            Err(SparseStructuralTimeError::Cancelled)
        ));
        assert_eq!(state_values(&state), before);
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(
            state_values(&owner.step(&state, &[0.2], &[0.1], &|| false).unwrap().state),
            state_values(&expected.state),
        );
    }
}
