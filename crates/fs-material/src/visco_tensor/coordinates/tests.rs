use super::*;
use crate::visco_tensor::TensorPronyError;

fn image(rows: &[Vec<f64>], times: &[f64], n: usize) -> RelaxationCoordinates {
    RelaxationCoordinates::new(
        rows,
        times,
        n,
        RelaxationCoordinateBudget {
            maximum_states: 256,
            maximum_scalar_terms: 1_000_000,
        },
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}
fn target(rows: &[Vec<f64>], q: &[f64]) -> Vec<f64> {
    rows.iter()
        .map(|row| row.iter().zip(q).map(|(a, q)| a * q).sum())
        .collect()
}
fn force(rows: &[Vec<f64>], q: &[f64], z: &[f64]) -> Vec<f64> {
    let goal = target(rows, q);
    (0..q.len())
        .map(|i| {
            rows.iter()
                .enumerate()
                .map(|(j, row)| row[i] * (goal[j] - z[j]))
                .sum()
        })
        .collect()
}

#[test]
fn g1_shared_coordinates_preserve_every_original_history_force_storage_and_heat() {
    let rows = vec![
        vec![1.0, 2.0],
        vec![-3.0, 1.0],
        vec![0.0, 0.0],
        vec![0.3, -0.8],
        vec![2.0, 1.0],
        vec![-0.2, 1.8],
        vec![0.0, -0.4],
    ];
    let times = vec![0.2, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7];
    let shared = image(&rows, &times, 2);
    assert_eq!(shared.coordinates(), 4);
    assert_eq!(shared.original_coordinates(), 7);
    let map = shared.memory_map(0..rows.len()).unwrap();
    for relaxed in [false, true] {
        let initial = [0.3, -0.6];
        let mut full = if relaxed {
            target(&rows, &initial)
        } else {
            vec![0.0; rows.len()]
        };
        let mut compact = if relaxed {
            target(shared.coupling_rows(), &initial)
        } else {
            vec![0.0; shared.coordinates()]
        };
        let mut full_heat = 0.0;
        let mut shared_heat = 0.0;
        for q in [[0.2, -0.1], [-0.7, 0.4], [0.0, 0.0], [1.0, 0.3]] {
            for (z, (goal, time)) in full
                .iter_mut()
                .zip(target(&rows, &q).into_iter().zip(&times))
            {
                let e = fs_math::det::exp(-0.03 / time);
                let mismatch = goal - *z;
                full_heat += 0.5 * mismatch * mismatch * (1.0 - e * e);
                *z = goal - mismatch * e;
            }
            for (z, (goal, time)) in compact.iter_mut().zip(
                target(shared.coupling_rows(), &q)
                    .into_iter()
                    .zip(shared.relaxation_times_s()),
            ) {
                let e = fs_math::det::exp(-0.03 / time);
                let mismatch = goal - *z;
                shared_heat += 0.5 * mismatch * mismatch * (1.0 - e * e);
                *z = goal - mismatch * e;
            }
            for (a, b) in full.iter().zip(map.reconstruct(&compact).unwrap()) {
                close(*a, b);
            }
            for (a, b) in
                force(&rows, &q, &full)
                    .iter()
                    .zip(force(shared.coupling_rows(), &q, &compact))
            {
                close(*a, b);
            }
            let energy = |rows: &[Vec<f64>], z: &[f64]| {
                target(rows, &q)
                    .iter()
                    .zip(z)
                    .map(|(a, z)| 0.5 * (a - z).powi(2))
                    .sum::<f64>()
            };
            close(
                energy(&rows, &full),
                energy(shared.coupling_rows(), &compact),
            );
            close(full_heat, shared_heat);
            let point = shared.memory_map(2..5).unwrap();
            for (a, b) in point.reconstruct(&compact).unwrap().iter().zip(&full[2..5]) {
                close(*a, *b);
            }
        }
    }
}

#[test]
fn g0_shared_coordinates_keep_weak_zero_and_rank_deficient_directions_and_distinct_times() {
    let rows = vec![
        vec![1.0, 1e-200],
        vec![2.0, 2e-200],
        vec![0.0, 0.0],
        vec![-1.0, -1e-200],
    ];
    let shared = image(&rows, &[1.0; 4], 2);
    assert_eq!(shared.coordinates(), 2);
    assert!(shared.coupling_rows().iter().any(|r| r[1] != 0.0));
    let q = [0.4, 0.1];
    let actual = shared
        .memory_map(0..4)
        .unwrap()
        .reconstruct(&target(shared.coupling_rows(), &q))
        .unwrap();
    for (a, b) in actual.into_iter().zip(target(&rows, &q)) {
        close(a, b);
    }
    let nearby = f64::from_bits(1.0_f64.to_bits() + 1);
    let times = [1.0, nearby, 1.0, nearby];
    let separated = image(&rows, &times, 2);
    assert_eq!(separated.coordinates(), 4);
    assert_eq!(separated.relaxation_times_s(), &[1.0, 1.0, nearby, nearby]);
    let zero = image(&vec![vec![0.0; 2]; 6], &[1.0; 6], 2);
    assert_eq!(zero.coordinates(), 2);
    assert_eq!(
        zero.memory_map(0..6)
            .unwrap()
            .reconstruct(&[0.0, 0.0])
            .unwrap(),
        vec![0.0; 6]
    );
}

#[test]
fn g0_shared_coordinates_refuse_capacity_domains_and_cancellation_without_partial_image() {
    let rows = vec![vec![1.0, 2.0]; 4];
    let times = vec![0.5; 4];
    let budget = RelaxationCoordinateBudget {
        maximum_states: 2,
        maximum_scalar_terms: 256,
    };
    assert!(RelaxationCoordinates::new(&rows, &times, 2, budget, || false).is_ok());
    assert!(
        RelaxationCoordinates::new(
            &rows,
            &times,
            2,
            RelaxationCoordinateBudget {
                maximum_states: 1,
                ..budget
            },
            || false
        )
        .is_err()
    );
    assert!(
        RelaxationCoordinates::new(
            &rows,
            &times,
            2,
            RelaxationCoordinateBudget {
                maximum_scalar_terms: 255,
                ..budget
            },
            || false
        )
        .is_err()
    );
    assert!(RelaxationCoordinates::new(&rows, &[0.0; 4], 2, budget, || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        RelaxationCoordinates::new(&rows, &times, 2, budget, || {
            calls += 1;
            calls == 8
        }),
        Err(TensorPronyError::Cancelled)
    ));
    let shared = image(&rows, &times, 2);
    assert!(shared.memory_map(0..5).is_err());
    assert!(shared.memory_map(3..2).is_err());
    let map = shared.memory_map(0..4).unwrap();
    assert!(map.reconstruct(&[0.0]).is_err());
    assert!(map.reconstruct(&[f64::NAN, 0.0]).is_err());
}
