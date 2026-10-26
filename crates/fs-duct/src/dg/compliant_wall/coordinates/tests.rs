use super::*;

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale.max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn product(matrix: &[f64], x: &[f64]) -> Vec<f64> {
    matrix
        .chunks_exact(x.len())
        .map(|row| row.iter().zip(x).map(|(a, b)| a * b).sum())
        .collect()
}

/// Independent sealed-pocket solutions with different bulk moduli. Only the
/// original production coordinate/stamp functions construct the operator;
/// RK4 here is a test integrator, not a replacement for the DG exponential.
#[test]
fn heterogeneous_two_region_dynamics_match_independent_oscillators() {
    let (area, perimeter, sigma, stiffness) = (0.0003_f64, 0.06_f64, 0.2_f64, 1e7_f64);
    let lengths: [f64; 2] = [0.18, 0.39];
    let bulk: [f64; 2] = [1e5, 2.5e5];
    let pressure_map = [
        (bulk[0] / (area * lengths[0])).sqrt(),
        0.0,
        0.0,
        (bulk[1] / (area * lengths[1])).sqrt(),
    ];
    let root = [lengths[0].sqrt(), 0.0, 0.0, lengths[1].sqrt()];
    let inverse = [1.0 / root[0], 0.0, 0.0, 1.0 / root[3]];
    let coordinates = WallCoordinates::new(
        2,
        &pressure_map,
        &root,
        &inverse,
        perimeter,
        sigma,
        &mut || false,
    )
    .unwrap();
    let mut operator = vec![0.0; 36];
    coordinates
        .stamp(&mut operator, 6, 0, 4, &mut || false)
        .unwrap();
    for i in 0..2 {
        operator[(2 + i) * 6 + 4 + i] = (stiffness / sigma).sqrt();
        operator[(4 + i) * 6 + 2 + i] = -(stiffness / sigma).sqrt();
    }
    let p0 = 10.0;
    let mut x = vec![
        p0 / pressure_map[0],
        p0 / pressure_map[3],
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let initial: f64 = 0.5 * x.iter().map(|v| v * v).sum::<f64>();
    let dt = 2e-8;
    for step in 1..=500 {
        let k1 = product(&operator, &x);
        let shifted = |k: &[f64], scale| {
            x.iter()
                .zip(k)
                .map(|(x, k)| x + scale * k)
                .collect::<Vec<_>>()
        };
        let k2 = product(&operator, &shifted(&k1, 0.5 * dt));
        let k3 = product(&operator, &shifted(&k2, 0.5 * dt));
        let k4 = product(&operator, &shifted(&k3, dt));
        for i in 0..6 {
            x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
        if step % 50 != 0 {
            continue;
        }
        let velocity = product(&coordinates.velocity_map, &x[4..]);
        let displacement = product(&coordinates.velocity_map, &x[2..4]);
        let mut physical_energy = 0.0;
        for i in 0..2 {
            let beta = perimeter * bulk[i] / area;
            let omega = ((stiffness + beta) / sigma).sqrt();
            let (sin, cos) = (omega * step as f64 * dt).sin_cos();
            let w = p0 / (stiffness + beta) * (1.0 - cos);
            let v = p0 / (sigma * omega) * sin;
            let p = pressure_map[i * 3] * x[i];
            close(p, p0 - beta * w, p0);
            close(
                displacement[i] * (sigma / stiffness).sqrt(),
                w,
                p0 / (stiffness + beta),
            );
            close(velocity[i], v, p0 / (sigma * omega));
            physical_energy += 0.5
                * lengths[i]
                * (area / bulk[i] * p * p + perimeter * (sigma * v * v + stiffness * w * w));
        }
        close(physical_energy, initial, initial);
        close(0.5 * x.iter().map(|v| v * v).sum::<f64>(), initial, initial);
    }
    // At 10 microseconds the independent solutions differ by about 0.075 Pa.
    assert!((pressure_map[0] * x[0] - pressure_map[3] * x[1]).abs() > 0.01);
}

/// Roots of heterogeneous pressure and geometric masses need not commute.
/// Check physical weak continuity and wall acceleration, not just skewness.
#[test]
fn noncommuting_roots_preserve_weak_force_and_physical_energy() {
    let coordinates = WallCoordinates::new(
        2,
        &[2.0, 0.0, 0.0, 1.0],
        &[2.0, 0.5, 0.5, 1.0],
        &[4.0 / 7.0, -2.0 / 7.0, -2.0 / 7.0, 8.0 / 7.0],
        3.0,
        2.0,
        &mut || false,
    )
    .unwrap();
    let mut operator = vec![0.0; 16];
    coordinates
        .stamp(&mut operator, 4, 0, 2, &mut || false)
        .unwrap();
    assert_ne!(operator[3], operator[6]);
    let x = [0.4, -0.7, 0.3, 0.8];
    let derivative = product(&operator, &x);
    let velocity = product(&coordinates.velocity_map, &x[2..]);
    let acceleration = product(&coordinates.velocity_map, &derivative[2..]);
    let gram = [4.25, 1.5, 1.5, 1.25];
    let swept_flow = product(&gram, &velocity);
    let pressure = [2.0 * x[0], x[1]];
    for i in 0..2 {
        let weak_derivative = derivative[i] * [0.5, 1.0][i];
        close(weak_derivative, -3.0 * swept_flow[i], 10.0);
        close(acceleration[i], pressure[i] / 2.0, 1.0);
    }
    close(
        x.iter().zip(&derivative).map(|(x, d)| x * d).sum(),
        0.0,
        10.0,
    );
    let physical_energy = 0.5
        * (0.25 * pressure[0].powi(2)
            + pressure[1].powi(2)
            + 6.0
                * velocity
                    .iter()
                    .zip(&swept_flow)
                    .map(|(v, gv)| v * gv)
                    .sum::<f64>());
    close(
        physical_energy,
        0.5 * x.iter().map(|v| v * v).sum::<f64>(),
        1.0,
    );
}

#[test]
fn uniform_limit_and_refusals_preserve_every_positive_coupling() {
    let root = [2.0, 0.5, 0.5, 1.0];
    let inverse = [4.0 / 7.0, -2.0 / 7.0, -2.0 / 7.0, 8.0 / 7.0];
    let pressure: Vec<_> = inverse.iter().map(|v| 7.0 * v).collect();
    let coordinates =
        WallCoordinates::new(2, &pressure, &root, &inverse, 3.0, 2.0, &mut || false).unwrap();
    for (i, value) in coordinates.coupling.iter().enumerate() {
        close(
            *value,
            if i == 0 || i == 3 {
                7.0 * (3.0_f64 / 2.0).sqrt()
            } else {
                0.0
            },
            10.0,
        );
    }
    assert!(matches!(
        WallCoordinates::new(2, &pressure, &root, &inverse, 3.0, 2.0, &mut || true),
        Err(CoordinateError::Cancelled)
    ));
    let mut calls = 0;
    assert!(matches!(
        WallCoordinates::new(2, &pressure, &root, &inverse, 3.0, 2.0, &mut || {
            calls += 1;
            calls == 3
        }),
        Err(CoordinateError::Cancelled)
    ));
    for mass in [0.0, -1.0, f64::INFINITY, 1e-320] {
        assert!(
            WallCoordinates::new(2, &pressure, &root, &inverse, 3.0, mass, &mut || false).is_err()
        );
    }
    assert!(
        WallCoordinates::new(2, &pressure[..3], &root, &inverse, 3.0, 2.0, &mut || false).is_err()
    );
    let mut operator = vec![0.0; 16];
    assert_eq!(
        coordinates.stamp(&mut operator, 4, 0, 1, &mut || false),
        Err(CoordinateError::Invalid)
    );
    assert!(operator.iter().all(|v| *v == 0.0));
    assert_eq!(
        coordinates.stamp(&mut operator, 4, 0, 2, &mut || true),
        Err(CoordinateError::Cancelled)
    );
    assert!(operator.iter().all(|v| *v == 0.0));
}
