//! G1/G3: independent physical fields, relaxation/dispersion and original work.
use super::*;

fn matrix(n: usize, values: &[f64]) -> Csr {
    let mut c = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if values[i * n + j] != 0.0 {
                c.push(i, j, values[i * n + j]);
            }
        }
    }
    c.assemble()
}
fn scalar(x: f64) -> Csr {
    matrix(1, &[x])
}
fn arm(strength: f64, tau: f64) -> BulkMemoryArm {
    BulkMemoryArm {
        delta_compressibility_pa_inv: strength,
        relaxation_time_s: tau,
    }
}
fn near(x: f64, y: f64, tol: f64) {
    assert!((x - y).abs() <= tol, "{x:e} != {y:e}, tolerance {tol:e}");
}
fn one(dt: f64) -> PreparedBulkDynamics {
    PreparedBulkDynamics::new(
        scalar(2.0),
        scalar(0.0),
        scalar(0.0),
        vec![Region {
            dofs: vec![0],
            mass: scalar(0.4),
            arms: vec![arm(0.8, 0.2)],
        }],
        vec![0],
        dt,
        BulkTimeLimits::default(),
        &mut || false,
    )
    .unwrap()
}
#[test]
fn closed_original_compression_relaxes_with_conserved_strain_and_independent_heat() {
    // e=beta_inf*p+delta_beta*q is constant. One physical decay rate follows
    // directly from the two continuous equations, independently of midpoint.
    let mut errors = Vec::new();
    for h in [0.04, 0.02, 0.01] {
        let owner = one(h);
        let mut state = owner
            .state(&[0.0], &[1.0], &[vec![0.0]], &mut || false)
            .unwrap();
        let start = owner.energy(&state, &mut || false).unwrap().total_j();
        let mut heat = 0.0;
        for _ in 0..(0.8 / h).round() as usize {
            let step = owner.step(&state, &[0.0], &mut || false).unwrap();
            heat += step.work.molecular_dissipation_j;
            assert!(step.work.molecular_dissipation_j > 0.0);
            state = step.state;
            near(
                2.0 * state.pressure[0] + 0.32 * state.memories[0][0],
                2.0,
                2e-13,
            );
        }
        let q = 2.0 / 2.32 * (1.0 - det::exp(-1.16 * 0.8 / 0.2));
        let p = (2.0 - 0.32 * q) / 2.0;
        errors.push(
            (state.pressure[0] - p)
                .abs()
                .max((state.memories[0][0] - q).abs()),
        );
        near(
            owner.energy(&state, &mut || false).unwrap().total_j() + heat,
            start,
            3e-13,
        );
    }
    eprintln!("bulk closed relaxation errors {errors:?}");
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
}
fn driven(h: f64) -> PreparedBulkDynamics {
    PreparedBulkDynamics::new(
        matrix(2, &[2.0, 0.2, 0.2, 1.0]),
        matrix(2, &[1.0, -1.0, -1.0, 1.0]),
        matrix(2, &[0.1, 0.02, 0.02, 0.2]),
        vec![Region {
            dofs: vec![0, 1],
            mass: matrix(2, &[0.6, 0.1, 0.1, 0.4]),
            arms: vec![arm(0.3, 0.2), arm(0.5, 0.06)],
        }],
        vec![0, 1],
        h,
        BulkTimeLimits::default(),
        &mut || false,
    )
    .unwrap()
}
fn continuous(t: f64) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
    let rates = [0.7, -0.4];
    let p: Vec<_> = rates.map(|a| det::exp(a * t)).to_vec();
    let psi: Vec<_> = (0..2).map(|i| (p[i] - 1.0) / rates[i]).collect();
    let q: Vec<Vec<f64>> = [0.2, 0.06]
        .iter()
        .map(|&tau| (0..2).map(|i| p[i] / (1.0 + rates[i] * tau)).collect())
        .collect();
    let mut f = vec![0.0; 2];
    let m = [[2.0, 0.2], [0.2, 1.0]];
    let k = [[1.0, -1.0], [-1.0, 1.0]];
    let c = [[0.1, 0.02], [0.02, 0.2]];
    let b = [[0.6, 0.1], [0.1, 0.4]];
    for i in 0..2 {
        for j in 0..2 {
            f[i] += m[i][j] * rates[j] * p[j]
                + k[i][j] * psi[j]
                + c[i][j] * p[j]
                + b[i][j] * (0.3 * rates[j] * q[0][j] + 0.5 * rates[j] * q[1][j]);
        }
    }
    (psi, p, q, f)
}
#[test]
fn complete_non_diagonal_fields_converge_to_independent_material_forcing() {
    let mut errors = Vec::new();
    for h in [0.04, 0.02, 0.01] {
        let owner = driven(h);
        let (psi, p, q, _) = continuous(0.0);
        let mut state = owner.state(&psi, &p, &q, &mut || false).unwrap();
        for step in 0..(0.4 / h).round() as usize {
            let (_, _, _, f) = continuous((step as f64 + 0.5) * h);
            let next = owner.step(&state, &f, &mut || false).unwrap();
            assert!(next.memories.iter().all(|m| m.dissipation_j > 0.0));
            state = next.state;
        }
        let (psi, p, q, _) = continuous(0.4);
        let mut e = 0.0_f64;
        for (a, b) in state
            .pressure
            .iter()
            .zip(&p)
            .chain(state.potential.iter().zip(&psi))
            .chain(state.memories.iter().flatten().zip(q.iter().flatten()))
        {
            e = e.max((a - b).abs());
        }
        errors.push(e);
    }
    eprintln!("bulk manufactured errors {errors:?}");
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
    assert!(errors[2] < 2e-5);
}
#[test]
fn original_harmonic_compressibility_recovers_passive_dispersion() {
    // An independently prescribed harmonic pressure supplies the exact source
    // from beta(omega); q has the positive causal lag at each original tau.
    let omega = 2.3;
    let mut errors = Vec::new();
    for h in [0.04, 0.02, 0.01] {
        let owner = one(h);
        let a = omega * 0.2;
        let initial = vec![vec![1.0 / (1.0 + a * a)]];
        let mut state = owner
            .state(&[0.0], &[1.0], &initial, &mut || false)
            .unwrap();
        for i in 0..(2.0 / h).round() as usize {
            let phase = omega * (i as f64 + 0.5) * h;
            let (s, c) = det::sin_cos(phase);
            let qdot = omega * (-s + a * c) / (1.0 + a * a);
            let f = -2.0 * omega * s + 0.32 * qdot;
            state = owner.step(&state, &[f], &mut || false).unwrap().state;
        }
        let (s, c) = det::sin_cos(omega * 2.0);
        let q = (c + a * s) / (1.0 + a * a);
        errors.push(
            (state.pressure[0] - c)
                .abs()
                .max((state.memories[0][0] - q).abs()),
        );
    }
    eprintln!("bulk harmonic errors {errors:?}");
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
}
#[test]
fn pressure_release_retains_original_volume_histories_and_their_heat() {
    let owner = PreparedBulkDynamics::new(
        matrix(2, &[1.0, 0.1, 0.1, 1.0]),
        matrix(2, &[1.0, -1.0, -1.0, 1.0]),
        matrix(2, &[0.0; 4]),
        vec![Region {
            dofs: vec![0, 1],
            mass: matrix(2, &[1.0, 0.3, 0.3, 1.0]),
            arms: vec![arm(0.2, 0.05)],
        }],
        vec![1],
        0.01,
        BulkTimeLimits::default(),
        &mut || false,
    )
    .unwrap();
    let state = owner
        .state(&[0.0, 0.0], &[0.0, 0.0], &[vec![1.0, 0.0]], &mut || false)
        .unwrap();
    let step = owner.step(&state, &[0.0, 0.0], &mut || false).unwrap();
    assert_eq!(step.state.pressure[0], 0.0);
    assert_eq!(step.state.potential[0], 0.0);
    assert!(step.state.memories[0][0] < 1.0);
    assert!(step.work.molecular_dissipation_j > 0.0);
    assert!(step.state.pressure[1] > 0.0);
    near(step.work.balance_residual_j(), 0.0, 2e-14);
}
#[test]
fn material_partition_and_coordinate_scaling_preserve_full_fields_and_work() {
    let original = driven(0.03);
    let scale = [2.0, 0.5];
    let transform = |a: &Csr| {
        let mut c = Coo::new(2, 2);
        for i in 0..2 {
            let (cols, v) = a.row(i);
            for (&j, &v) in cols.iter().zip(v) {
                c.push(i, j, v * scale[i] * scale[j]);
            }
        }
        c.assemble()
    };
    let transformed = PreparedBulkDynamics::new(
        transform(&original.mass),
        transform(&original.stiffness),
        transform(&original.conductance),
        vec![Region {
            dofs: vec![0, 1],
            mass: transform(&original.regions[0].mass),
            arms: original.regions[0].arms.clone(),
        }],
        vec![0, 1],
        0.03,
        BulkTimeLimits::default(),
        &mut || false,
    )
    .unwrap();
    let (psi, p, q, f) = continuous(0.0);
    let a = original.state(&psi, &p, &q, &mut || false).unwrap();
    let rescale = |x: &[f64]| x.iter().zip(scale).map(|(&x, s)| x / s).collect::<Vec<_>>();
    let b = transformed
        .state(
            &rescale(&psi),
            &rescale(&p),
            &q.iter().map(|q| rescale(q)).collect::<Vec<_>>(),
            &mut || false,
        )
        .unwrap();
    let x = original.step(&a, &f, &mut || false).unwrap();
    let y = transformed
        .step(
            &b,
            &f.iter().zip(scale).map(|(&f, s)| f * s).collect::<Vec<_>>(),
            &mut || false,
        )
        .unwrap();
    for i in 0..2 {
        near(x.state.pressure[i], y.state.pressure[i] * scale[i], 2e-13);
        for j in 0..2 {
            near(
                x.state.memories[j][i],
                y.state.memories[j][i] * scale[i],
                2e-13,
            );
        }
    }
    near(
        x.work.final_energy.total_j(),
        y.work.final_energy.total_j(),
        5e-13,
    );
    near(
        x.work.molecular_dissipation_j,
        y.work.molecular_dissipation_j,
        2e-14,
    );
}
#[test]
fn empty_memories_keep_lossless_gauge_and_large_sparse_coordinate_space() {
    let n = 256;
    let mut k = Coo::new(n, n);
    let mut m = Coo::new(n, n);
    for i in 0..n {
        m.push(i, i, 1.0);
        if i + 1 < n {
            k.push(i, i, 0.4);
            k.push(i + 1, i + 1, 0.4);
            k.push(i, i + 1, -0.4);
            k.push(i + 1, i, -0.4);
        }
    }
    let owner = PreparedBulkDynamics::new(
        m.assemble(),
        k.assemble(),
        Coo::new(n, n).assemble(),
        vec![],
        (0..n).collect(),
        0.1,
        BulkTimeLimits::default(),
        &mut || false,
    )
    .unwrap();
    let state = owner
        .state(&vec![3.0; n], &vec![0.0; n], &[], &mut || false)
        .unwrap();
    let next = owner.step(&state, &vec![0.0; n], &mut || false).unwrap();
    assert_eq!(next.state.potential, vec![3.0; n]);
    assert_eq!(next.state.pressure, vec![0.0; n]);
    assert_eq!(next.work.molecular_dissipation_j, 0.0);
}
#[test]
fn original_source_limits_cancellation_and_foreign_histories_refuse_transactionally() {
    let owner = one(0.01);
    let foreign = one(0.01);
    let state = owner
        .state(&[0.2], &[0.4], &[vec![0.1]], &mut || false)
        .unwrap();
    assert!(foreign.step(&state, &[0.0], &mut || false).is_err());
    assert!(owner.state(&[0.0], &[0.0], &[], &mut || false).is_err());
    assert!(owner.step(&state, &[f64::NAN], &mut || false).is_err());
    for stop in [1, 3, 8, 20] {
        let mut calls = 0;
        let result = owner.step(&state, &[0.8], &mut || {
            calls += 1;
            calls == stop
        });
        assert!(matches!(result, Err(BulkTimeError::Cancelled)));
        assert_eq!(state.pressure, [0.4]);
        assert_eq!(state.memories, [vec![0.1]]);
    }
    let a = owner.step(&state, &[0.8], &mut || false).unwrap();
    let b = owner.step(&state, &[0.8], &mut || false).unwrap();
    assert_eq!(a.state.pressure, b.state.pressure);
    assert_eq!(a.state.memories, b.state.memories);
    let mut limits = BulkTimeLimits::default();
    limits.maximum_memory_dofs = 0;
    assert!(
        PreparedBulkDynamics::new(
            scalar(2.0),
            scalar(0.0),
            scalar(0.0),
            vec![Region {
                dofs: vec![0],
                mass: scalar(0.4),
                arms: vec![arm(0.8, 0.2)]
            }],
            vec![0],
            0.1,
            limits,
            &mut || false
        )
        .is_err()
    );
}
