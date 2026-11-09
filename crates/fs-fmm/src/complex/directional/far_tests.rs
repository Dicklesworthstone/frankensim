use super::*;

fn cloud() -> (Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<C64>) {
    let mut points = Vec::new();
    let mut normals = Vec::new();
    let mut charges = Vec::new();
    for center in [-1.0, 1.0] {
        for x in 0..3 {
            for y in 0..3 {
                for z in 0..3 {
                    points.push([
                        0.009 * (x as f64 - 1.0),
                        0.007 * (y as f64 - 1.0),
                        center + 0.008 * (z as f64 - 1.0),
                    ]);
                    normals.push(if (x + y + z) % 2 == 0 {
                        [0.6, 0.0, 0.8]
                    } else {
                        [0.0, 0.8, 0.6]
                    });
                    charges.push(C64::new(
                        1.0 + 0.03 * (x + 2 * y) as f64,
                        0.07 * (z as f64 - 1.0),
                    ));
                }
            }
        }
    }
    (points, normals, charges)
}
fn settings(order: usize) -> DirectionalOptions {
    DirectionalOptions {
        order,
        leaf_capacity: 16,
        ..DirectionalOptions::for_tolerance(1e-2).unwrap()
    }
}
fn engine(layer: usize, k: f64, order: usize, far: bool, max_near_pairs: usize) -> FmmC64<'static> {
    let (points, normals, _) = cloud();
    let options = DirectionalOptions {
        max_near_pairs,
        ..settings(order)
    };
    match (layer, far) {
        (0, true) => {
            FmmC64::new_single_far_field_directional(HelmholtzSingle { k }, points, options)
        }
        (1, true) => FmmC64::new_double_far_field_directional(
            HelmholtzDouble { k },
            points,
            normals,
            options,
        ),
        (_, true) => FmmC64::new_hypersingular_far_field_directional(
            HelmholtzHypersingular { k },
            points,
            normals,
            options,
        ),
        (0, false) => FmmC64::new_single_directional(HelmholtzSingle { k }, &points, options),
        (1, false) => {
            FmmC64::new_double_directional(HelmholtzDouble { k }, points, normals, options)
        }
        (_, false) => FmmC64::new_hypersingular_directional(
            HelmholtzHypersingular { k },
            points,
            normals,
            options,
        ),
    }
    .unwrap()
}
fn direct_parts(op: &FmmC64<'_>, charges: &[C64]) -> (Vec<C64>, Vec<C64>) {
    let blocks = op.near_field_blocks().collect::<Vec<_>>();
    let neighbours = op.near_field_neighbor_leaves();
    let mut leaf = vec![0; charges.len()];
    for (i, points) in blocks.iter().enumerate() {
        for &p in points.iter() {
            leaf[p] = i;
        }
    }
    let (mut far, mut near) = (
        vec![C64::ZERO; charges.len()],
        vec![C64::ZERO; charges.len()],
    );
    for t in 0..charges.len() {
        for s in 0..charges.len() {
            if t == s {
                continue;
            }
            let value = op.pair(t, s) * charges[s];
            if neighbours[leaf[t]].binary_search(&leaf[s]).is_ok() {
                near[t] = near[t] + value;
            } else {
                far[t] = far[t] + value;
            }
        }
    }
    (far, near)
}
fn relative(a: &[C64], b: &[C64]) -> f64 {
    (a.iter()
        .zip(b)
        .map(|(a, b)| (*a - *b).norm_sq())
        .sum::<f64>()
        / b.iter().map(|v| v.norm_sq()).sum::<f64>())
    .sqrt()
}

#[test]
fn directional_far_layers_match_partitioned_direct_under_order_refinement() {
    let (_, _, charges) = cloud();
    for layer in 0..3 {
        let mut errors = Vec::new();
        for order in [4, 6, 8] {
            let op = engine(layer, 34.0, order, true, 1);
            let stats = op.directional_stats().unwrap();
            assert!(stats.high_frequency_pairs > 0);
            assert_eq!(stats.near_pairs, 0);
            let expected = direct_parts(&op, &charges).0;
            let actual = op.far_potentials(&charges).unwrap();
            errors.push(relative(&actual, &expected));
            assert!(op.potentials(&charges).is_err());
        }
        eprintln!("directional far layer={layer}: order4/6/8 errors={errors:?}");
        assert!(errors[1] < errors[0] && errors[2] < errors[1]);
        assert!(errors[2] < 1e-5);
    }
}

#[test]
fn directional_far_partition_reconstructs_full_fields_and_keeps_worker_bits() {
    let (_, _, charges) = cloud();
    for layer in 0..3 {
        let mut far = engine(layer, 34.0, 5, true, 1);
        let full = engine(layer, 34.0, 5, false, 100000);
        assert_eq!(
            far.near_field_neighbor_leaves(),
            full.near_field_neighbor_leaves()
        );
        let values = far.far_potentials(&charges).unwrap();
        let exact_near = direct_parts(&far, &charges).1;
        let combined = values
            .iter()
            .zip(exact_near)
            .map(|(&f, n)| f + n)
            .collect::<Vec<_>>();
        assert!(relative(&combined, &full.potentials(&charges).unwrap()) < 3e-14);
        far.set_workers(3).unwrap();
        assert_eq!(values, far.far_potentials(&charges).unwrap());
        assert!(full.directional_stats().unwrap().near_pairs > 1);
    }
    let (points, _, _) = cloud();
    assert!(
        FmmC64::new_single_directional(
            HelmholtzSingle { k: 34.0 },
            &points,
            DirectionalOptions {
                max_near_pairs: 1,
                ..settings(5)
            }
        )
        .is_err()
    );
}

#[test]
fn directional_far_low_frequency_delegates_to_original_w5a_bits() {
    let (points, normals, charges) = cloud();
    for layer in 0..3 {
        let a = engine(layer, 0.8, 4, true, 1);
        let b = match layer {
            0 => FmmC64::new_single_far_field_within(
                HelmholtzSingle { k: 0.8 },
                points.clone(),
                4,
                16,
                ComplexBudget::V1,
            ),
            1 => FmmC64::new_double_far_field_within(
                HelmholtzDouble { k: 0.8 },
                points.clone(),
                normals.clone(),
                4,
                16,
                ComplexBudget::V1,
            ),
            _ => FmmC64::new_hypersingular_far_field_within(
                HelmholtzHypersingular { k: 0.8 },
                points.clone(),
                normals.clone(),
                4,
                16,
                ComplexBudget::V1,
            ),
        }
        .unwrap();
        assert!(a.directional_stats().is_none());
        assert_eq!(
            a.near_field_neighbor_leaves(),
            b.near_field_neighbor_leaves()
        );
        assert_eq!(
            a.far_potentials(&charges).unwrap(),
            b.far_potentials(&charges).unwrap()
        );
    }
}
