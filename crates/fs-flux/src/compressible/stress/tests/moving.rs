use super::*;

fn moving_balance<const D: usize>(shape: [usize; D], spacing: [f64; D]) {
    let count: usize = shape.iter().product();
    let velocity: Vec<[f64; D]> = (0..count)
        .map(|i| std::array::from_fn(|a| ((7 * i + 3 * a + 1) as f64).sin()))
        .collect();
    let material: Vec<_> = (0..count)
        .map(|i| Viscosity {
            shear: 10.0_f64.powf(3.0 * ((5 * i + 1) as f64).sin()),
            bulk: 0.2 + ((11 * i + 2) as f64).cos().abs(),
        })
        .collect();
    let walls = std::array::from_fn(|axis| {
        std::array::from_fn(|side| {
            std::array::from_fn(|a| ((axis * 13 + side * 7 + a + 1) as f64).cos())
        })
    });
    for mode in 0..3 {
        let ends = std::array::from_fn(|axis| match mode {
            0 => [End::NoSlip; 2],
            1 => [End::Navier(0.17 * (axis + 1) as f64), End::Slip],
            _ => [End::Slip; 2],
        });
        let (rate, loss, force, power) =
            balance_with_walls(shape, spacing, ends, walls, &velocity, &material, true);
        assert!(loss >= 0.0);
        assert!(power.abs() > 1e-5, "normal wall work must remain active");
        assert!((power - rate - loss).abs() < 8e-14 * (power.abs() + rate.abs() + loss));

        let boost = std::array::from_fn::<_, D, _>(|a| 0.25 * (a + 1) as f64);
        let shifted: Vec<_> = velocity
            .iter()
            .map(|u| std::array::from_fn(|a| u[a] + boost[a]))
            .collect();
        let shifted_walls = std::array::from_fn(|axis| {
            std::array::from_fn(|side| std::array::from_fn(|a| walls[axis][side][a] + boost[a]))
        });
        let (shift_rate, shift_loss, shift_force, shift_power) = balance_with_walls(
            shape,
            spacing,
            ends,
            shifted_walls,
            &shifted,
            &material,
            true,
        );
        for (a, b) in force.iter().zip(&shift_force) {
            for component in 0..D {
                assert!((a[component] - b[component]).abs() < 2e-11 * (1.0 + a[component].abs()));
            }
        }
        let added_power: f64 = force
            .iter()
            .map(|f| (0..D).map(|a| boost[a] * f[a]).sum::<f64>())
            .sum();
        let scale = (power.abs() + rate.abs() + loss).max(1.0);
        assert!((shift_loss - loss).abs() < 8e-14 * scale);
        assert!((shift_power - power - added_power).abs() < 8e-14 * scale);
        assert!((shift_rate - rate - added_power).abs() < 8e-14 * scale);
    }
}

#[test]
fn g0_moving_normal_traces_close_independent_viscous_power_at_corners() {
    moving_balance([3, 2], [0.7, 1.3]);
    moving_balance([3, 2, 2], [0.7, 1.3, 0.9]);
    moving_balance([1, 3, 1], [0.7, 1.3, 0.9]);
    moving_balance([1, 1, 1], [0.7, 1.3, 0.9]);
}

#[test]
fn g1_normal_wall_recoil_uses_longitudinal_viscosity_and_laboratory_work() {
    // One cell, lower wall at rest, upper wall speed 1/2, cell speed 1/8.
    // mu=3/4 and bulk=1/2 give longitudinal viscosity 3/2. The two normal
    // gradients are 1/4 and 3/4, so stresses are 3/8 and 9/8. Their difference
    // gives force 3/4, kinetic power 3/32 and wall power 9/16.
    for axis in 0..3 {
        for condition in [End::NoSlip, End::Navier(0.23), End::Slip] {
            let mut ends = [[End::Periodic; 2]; 3];
            ends[axis] = [condition; 2];
            let mut walls = [[[0.0; 3]; 2]; 3];
            walls[axis][1][axis] = 0.5;
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.125;
            let (rate, loss, force, power) = balance_with_walls(
                [1; 3],
                [1.0; 3],
                ends,
                walls,
                &[velocity],
                &[Viscosity {
                    shear: 0.75,
                    bulk: 0.5,
                }],
                true,
            );
            let mut expected = [0.0; 3];
            expected[axis] = 0.75;
            assert_eq!(force, [expected]);
            assert_eq!(rate, 3.0 / 32.0);
            assert_eq!(power, 9.0 / 16.0);
            assert_eq!(loss, 15.0 / 32.0);
        }
    }
}

#[test]
fn g3_comoving_walls_have_zero_stress_for_every_slip_law() {
    let velocity = [0.5, -0.25, 0.75];
    let viscosity = Viscosity {
        shear: 0.7,
        bulk: 0.2,
    };
    for slip in [Slip::Free, Slip::Length(0.0), Slip::Length(0.3)] {
        let gradient = velocity_gradient([0.3, 0.7, 1.1], velocity, |_, _| GradientSample::Wall {
            velocity,
            slip,
        })
        .unwrap();
        assert_eq!(gradient, [[0.0; 3]; 3]);
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let result = wall_stress(
                    axis,
                    sign,
                    0.2,
                    velocity,
                    velocity,
                    slip,
                    viscosity,
                    |a, j| gradient[a][j],
                )
                .unwrap();
                assert_eq!(result.traction, [0.0; 3]);
                assert_eq!(result.slip_dissipation, 0.0);
            }
        }
    }
}
