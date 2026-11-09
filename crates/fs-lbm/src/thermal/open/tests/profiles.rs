use super::*;

fn plane_index(cell: usize, dims: [usize; 3], axis: usize) -> usize {
    let x = cell % dims[0];
    let y = cell / dims[0] % dims[1];
    let z = cell / (dims[0] * dims[1]);
    match axis {
        0 => y + dims[1] * z,
        1 => x + dims[0] * z,
        _ => x + dims[0] * y,
    }
}

#[test]
fn spatial_profiles_resolve_each_inflow_and_ignore_outflow_targets_on_every_axis() {
    let dims = [5, 4, 3];
    let cells = dims.into_iter().product();
    for axis in 0..3 {
        let entries = cells / dims[axis];
        let velocity = |cell| {
            let mut u = [0.0; 3];
            u[axis] = if plane_index(cell, dims, axis) % 2 == 0 {
                0.04
            } else {
                -0.03
            };
            u
        };
        let build = || {
            ThermalOpenField3::try_new::<ThermalOpenError3>(
                dims,
                &vec![false; cells],
                &vec![1.0; cells],
                0.8,
                faces(axis, 0.3, 1.7),
                limits(),
                0,
                |cell| Ok(velocity(cell)),
                || Ok(()),
            )
            .unwrap()
        };
        let mut actual = build();
        let mut changed_outflow = build();
        let profiles = [
            (0..entries)
                .map(|k| 0.2 + 0.013 * k as f64)
                .collect::<Vec<_>>(),
            (0..entries)
                .map(|k| 1.2 + 0.017 * k as f64)
                .collect::<Vec<_>>(),
        ];
        let mut expected_in = [0.0; 2];
        let mut expected_out = [0.0; 2];
        for side in 0..2 {
            let mut changed = profiles[side].clone();
            for k in 0..entries {
                let u = if k % 2 == 0 { 0.04 } else { -0.03 };
                let inward = if side == 0 { u } else { -u };
                if inward > 0.0 {
                    expected_in[side] += inward * profiles[side][k];
                } else {
                    // Initial outgoing fluid is exactly T=1. Outflow carries
                    // its own temperature, regardless of the prescribed bath.
                    expected_out[side] += inward;
                    changed[k] += 3.0;
                }
            }
            actual
                .set_reservoir_temperature_profile(2 * axis + side, &profiles[side])
                .unwrap();
            changed_outflow
                .set_reservoir_temperature_profile(2 * axis + side, &changed)
                .unwrap();
        }
        let receipt = step(&mut actual, &vec![0.0; cells], None);
        let other = step(&mut changed_outflow, &vec![0.0; cells], None);
        assert_eq!(bits(&actual), bits(&changed_outflow));
        assert_eq!(receipt, other);
        for side in 0..2 {
            let exchange = receipt.reservoirs[2 * axis + side];
            near(exchange.inflow_heat, expected_in[side], 3e-15);
            near(exchange.outflow_heat, expected_out[side], 3e-15);
            near(
                exchange.heat_in,
                expected_in[side] + expected_out[side],
                3e-15,
            );
        }
        // Independent population zeroth-moment prediction at each original
        // face node: inflow changes T by u_in*(T_res-T_initial); outflow does not.
        for cell in 0..cells {
            let coordinate = match axis {
                0 => cell % dims[0],
                1 => cell / dims[0] % dims[1],
                _ => cell / (dims[0] * dims[1]),
            };
            if coordinate != 0 && coordinate + 1 != dims[axis] {
                continue;
            }
            let side = usize::from(coordinate != 0);
            let inward = if side == 0 {
                velocity(cell)[axis]
            } else {
                -velocity(cell)[axis]
            };
            let k = plane_index(cell, dims, axis);
            let expected = if inward > 0.0 {
                1.0 + inward * (profiles[side][k] - 1.0)
            } else {
                1.0
            };
            near(actual.temperature(cell).unwrap(), expected, 2e-15);
        }
        near(
            receipt.stored_change,
            expected_in.iter().sum::<f64>() + expected_out.iter().sum::<f64>(),
            3e-14,
        );
        near(receipt.combined_defect, 0.0, 3e-14);
    }
}

#[test]
fn profile_updates_are_atomic_bounded_reusable_and_restore_the_uniform_fallback() {
    let dims = [5, 4, 3];
    let cells = dims.into_iter().product();
    let entries = cells / dims[0];
    let probe = field(dims, [0.03, 0.0, 0.0], -0.7, 0.8, faces(0, -0.7, 0.6));
    let mut bounded = limits();
    bounded.maximum_storage_bytes =
        probe.base_storage_bytes + entries * core::mem::size_of::<f64>();
    let mut actual = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &vec![false; cells],
        &vec![-0.7; cells],
        0.8,
        faces(0, -0.7, 0.6),
        bounded,
        0,
        |_| Ok([0.03, 0.0, 0.0]),
        || Ok(()),
    )
    .unwrap();
    let before = bits(&actual);
    let before_history = format!("{:?}", actual.history);
    let uniform = vec![-0.7; entries];
    actual
        .set_reservoir_temperature_profile(0, &uniform)
        .unwrap();
    let pointer = actual.reservoir_temperature_profile(0).unwrap().as_ptr();
    let updated: Vec<_> = (0..entries).map(|k| -0.6 + 0.01 * k as f64).collect();
    actual
        .set_reservoir_temperature_profile(0, &updated)
        .unwrap();
    assert_eq!(
        actual.reservoir_temperature_profile(0).unwrap().as_ptr(),
        pointer
    );
    assert!(
        actual
            .set_reservoir_temperature_profile(1, &uniform)
            .is_err()
    );
    assert!(actual.reservoir_temperature_profile(1).is_none());
    for value in [f64::NAN, f64::INFINITY, 10.01] {
        let mut refused = updated.clone();
        refused[entries - 1] = value;
        assert!(
            actual
                .set_reservoir_temperature_profile(0, &refused)
                .is_err()
        );
        assert_eq!(actual.reservoir_temperature_profile(0).unwrap(), updated);
    }
    assert!(
        actual
            .set_reservoir_temperature_profile(0, &updated[..entries - 1])
            .is_err()
    );
    for face in [2, 6] {
        assert!(
            actual
                .set_reservoir_temperature_profile(face, &updated)
                .is_err()
        );
        assert!(actual.clear_reservoir_temperature_profile(face).is_err());
    }
    assert!(actual.set_reservoir_temperature(0, 0.2).is_err());
    assert!(
        actual
            .set_wall_boundary(0, ThermalWall::Temperature(0.2))
            .is_err()
    );
    assert_eq!(
        actual.boundaries().faces[0],
        ThermalOpenFace3::Reservoir { temperature: -0.7 }
    );
    assert_eq!(bits(&actual), before);
    assert_eq!(format!("{:?}", actual.history), before_history);
    assert_eq!(actual.completed_steps(), 0);
    assert!(actual.last_balance().is_none());
    actual.clear_reservoir_temperature_profile(0).unwrap();
    assert!(actual.reservoir_temperature_profile(0).is_none());
    assert_eq!(
        actual.boundaries().faces[0],
        ThermalOpenFace3::Reservoir { temperature: -0.7 }
    );
    actual
        .set_reservoir_temperature_profile(1, &uniform)
        .unwrap();
    actual.set_reservoir_temperature(0, 0.2).unwrap();
}
