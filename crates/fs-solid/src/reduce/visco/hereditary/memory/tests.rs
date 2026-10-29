use super::*;

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 2e-12 * scale.max(1e-300), "{a:e} != {b:e}");
}
fn gram(rows: &[Vec<f64>], i: usize, j: usize) -> f64 {
    rows.iter().map(|r| r[i] * r[j]).sum()
}
fn energy_force(rows: &[Vec<f64>], elastic: &[f64]) -> (f64, Vec<f64>) {
    let mut force = vec![0.0; elastic.len()];
    let mut energy = 0.0;
    for row in rows {
        let z = row.iter().zip(elastic).map(|(a, x)| a * x).sum::<f64>();
        energy += 0.5 * z * z;
        for (f, a) in force.iter_mut().zip(row) {
            *f += a * z;
        }
    }
    (energy, force)
}

#[test]
fn equal_source_times_preserve_full_storage_reaction_heat_and_complex_response() {
    let first = vec![vec![2., -1., 0.4], vec![0., 3., -1.], vec![1., 0., 2.]];
    let second = vec![vec![-0.7, 1.2, 0.], vec![3., 0.4, -2.]];
    let third = vec![vec![0.2, 0.9, 1.5], vec![0., 0., 0.]];
    let grouped = compact(
        &[(0.03, &first), (0.03, &second), (0.1, &third)],
        3,
        100_000,
        || false,
    )
    .unwrap();
    assert_eq!(grouped.len(), 2);
    assert!(grouped.iter().all(|b| b.rows().len() <= 3));
    let all = [first, second].concat();
    let sources = [(0.03, all.as_slice()), (0.1, third.as_slice())];
    for ((tau, original), compact) in sources.iter().zip(&grouped) {
        assert_eq!(*tau, compact.relaxation_time_s());
        // Arbitrary elastic strain after a shared modal history, covering both
        // cold initializations and any subsequent displacement driving.
        for elastic in [[0.02, -0.03, 0.01], [-0.1, 0.07, 0.004], [0., 0., 0.]] {
            let (old_energy, old_force) = energy_force(original, &elastic);
            let (new_energy, new_force) = energy_force(compact.rows(), &elastic);
            close(new_energy, old_energy, old_energy);
            close(
                2. * new_energy / tau,
                2. * old_energy / tau,
                2. * old_energy / tau,
            );
            for (a, b) in new_force.iter().zip(old_force) {
                close(*a, b, 1.);
            }
        }
        for omega in [0., 1., 10., 33., 100., 1e6] {
            let x = omega * tau;
            let storage = x * x / (1. + x * x);
            let loss = x / (1. + x * x);
            for i in 0..3 {
                for j in 0..3 {
                    let old = gram(original, i, j);
                    let new = gram(compact.rows(), i, j);
                    let scale = gram(original, i, i).sqrt() * gram(original, j, j).sqrt();
                    close(storage * new, storage * old, storage * scale);
                    close(loss * new, loss * old, loss * scale);
                }
            }
        }
    }
}
#[test]
fn weak_physical_coordinate_survives_without_a_relative_rank_floor() {
    let rows = vec![
        vec![1e60, 0., 1.],
        vec![0., 1e-60, -0.3],
        vec![1e60, 1e-60, 0.7],
    ];
    let image = compact(&[(0.2, &rows)], 3, 10_000, || false).unwrap();
    let kept = image[0].rows();
    assert!(gram(kept, 1, 1) > 0.);
    for i in 0..3 {
        for j in 0..3 {
            close(
                gram(kept, i, j),
                gram(&rows, i, j),
                gram(&rows, i, i).sqrt() * gram(&rows, j, j).sqrt(),
            );
        }
    }
}
#[test]
fn exact_time_identity_and_explicit_work_and_cancellation_are_preserved() {
    let rows = vec![vec![1., -2.], vec![0.3, 0.1]];
    let other = f64::from_bits(0.1_f64.to_bits() + 1);
    assert_eq!(
        compact(&[(0.1, &rows), (other, &rows)], 2, 1000, || false)
            .unwrap()
            .len(),
        2
    );
    assert!(compact(&[(0.1, &rows)], 2, 1, || false).is_err());
    assert_eq!(
        compact(&[(0.1, &rows)], 2, 1000, || true),
        Err(CompressionError::Cancelled)
    );
    assert!(compact(&[(0.0, &rows)], 2, 1000, || false).is_err());
    let zero = vec![vec![0., 0.]; 6];
    assert!(compact(&[(0.1, &zero)], 2, 0, || false).unwrap().is_empty());
}

#[test]
fn compact_history_matches_original_tensor_material_ramps_and_independent_power() {
    use fs_material::visco::tensor::TensorProny;
    use fs_material::visco_tensor::{TensorProjectionBudget, TensorPronyBranch};
    let equilibrium = core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { 30. + 7. * i as f64 } else { 0. })
    });
    let mut viscosity = equilibrium.map(|r| r.map(|x| x * 0.002));
    viscosity[0][2] = 0.015;
    viscosity[2][0] = 0.015;
    let mut a = [[0.; 6]; 6];
    a[0] = [2., -3., 1., 0.4, -0.2, 0.7];
    a[1] = [0., 0.1, -0.7, 1.2, 0.8, -0.4];
    let mut b = [[0.; 6]; 6];
    b[0] = [-0.3, 0.5, 2.1, -1.5, 0., 0.9];
    let law = TensorProny::new(
        equilibrium,
        Some(viscosity),
        vec![
            TensorPronyBranch::from_gram_rows(a, 0.007).unwrap(),
            TensorPronyBranch::from_gram_rows(b, 0.04).unwrap(),
        ],
        (0., 1e6),
    )
    .unwrap();
    let columns = vec![
        [0.7, 0.2, 0.3, 0.1, -0.04, 0.08],
        [-0.1, 0.9, -0.3, 0.13, 0.15, -0.07],
        [0.4, -0.5, 0.6, -0.2, 0.09, 0.11],
    ];
    let other = columns
        .iter()
        .map(|c| c.map(|x| -1.3 * x))
        .collect::<Vec<_>>();
    let budget = TensorProjectionBudget {
        max_coordinates: 3,
        max_scalar_terms: 100_000,
    };
    let points = vec![
        law.project_strain(&columns, 0.023, budget, || false)
            .unwrap(),
        law.project_strain(&other, 0.037, budget, || false).unwrap(),
    ];
    let sources = points
        .iter()
        .flat_map(|p| p.branches())
        .map(|b| (b.relaxation_time_s(), b.rows()))
        .collect::<Vec<_>>();
    let compact = compact(&sources, 3, 1_000_000, || false).unwrap();
    let mut histories = points
        .iter()
        .map(|p| p.material().history())
        .collect::<Vec<_>>();
    let mut memory = compact
        .iter()
        .map(|b| vec![0.; b.rows().len()])
        .collect::<Vec<_>>();
    let mut k = vec![0.; 9];
    let mut c = vec![0.; 9];
    for point in &points {
        for i in 0..9 {
            k[i] += point.equilibrium_stiffness()[i];
            c[i] += point.viscosity_matrix()[i];
        }
    }
    let mut old = [0.; 3];
    for (next, dt) in [
        ([2e-4, -1e-4, 3e-4], 0.003),
        ([-1e-4, 4e-4, 1e-4], 0.002),
        ([3e-4, 2e-4, -2e-4], 0.0008),
        ([3e-4, 2e-4, -2e-4], 0.01),
    ] {
        let velocity: [f64; 3] = core::array::from_fn(|i| (next[i] - old[i]) / dt);
        let mut source_force = [0.; 3];
        let mut source_energy = 0.;
        let mut source_heat = 0.;
        let mut source_work = 0.;
        for (point, history) in points.iter().zip(&mut histories) {
            let report = history
                .advance(point.strain(&next).unwrap(), dt, || false)
                .unwrap();
            source_energy += point.volume_m3() * report.stored_energy;
            source_heat += point.volume_m3() * report.dissipated;
            source_work += point.volume_m3() * report.work;
            for (f, value) in source_force
                .iter_mut()
                .zip(point.conjugate_force(&report.stress).unwrap())
            {
                *f += value;
            }
        }
        let mut energy = 0.;
        let mut heat = 0.;
        let mut work = 0.;
        let mut force = [0.; 3];
        for i in 0..3 {
            for j in 0..3 {
                force[i] += k[i * 3 + j] * next[j] + c[i * 3 + j] * velocity[j];
                energy += 0.5 * next[i] * k[i * 3 + j] * next[j];
                heat += dt * velocity[i] * c[i * 3 + j] * velocity[j];
                work += 0.5 * (old[i] + next[i]) * k[i * 3 + j] * (next[j] - old[j])
                    + dt * velocity[i] * c[i * 3 + j] * velocity[j];
            }
        }
        for (branch, memory) in compact.iter().zip(&mut memory) {
            let tau = branch.relaxation_time_s();
            for (row, z) in branch.rows().iter().zip(memory) {
                let rate = row.iter().zip(velocity).map(|(a, v)| a * v).sum::<f64>();
                let start = *z;
                // Independent scalar stress-coordinate trajectory, sampled for
                // power; production uses source exact-ramp and exact block owners.
                let value = |t: f64| {
                    let decay = (-t / tau).exp();
                    decay * start - tau * (-t / tau).exp_m1() * rate
                };
                let panels = 128;
                for i in 0..=panels {
                    let weight = if i == 0 || i == panels {
                        1.
                    } else if i % 2 == 0 {
                        2.
                    } else {
                        4.
                    };
                    let stress = value(dt * i as f64 / panels as f64);
                    heat += dt / (3 * panels) as f64 * weight * stress * stress / tau;
                    work += dt / (3 * panels) as f64 * weight * stress * rate;
                }
                *z = value(dt);
                energy += 0.5 * (*z) * (*z);
                for (f, a) in force.iter_mut().zip(row) {
                    *f += a * (*z);
                }
            }
        }
        assert!((energy - source_energy).abs() < 2e-11 * source_energy.max(1e-30));
        assert!((heat - source_heat).abs() < 2e-8 * source_heat.max(1e-30));
        assert!((work - source_work).abs() < 2e-8 * (source_work.abs() + source_heat).max(1e-30));
        for (a, b) in force.into_iter().zip(source_force) {
            assert!((a - b).abs() < 2e-11 * (a.abs() + b.abs()).max(1e-30));
        }
        old = next;
    }
}
