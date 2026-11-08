use super::*;
use crate::acoustics::thermoviscous::time::ThermoviscousTimeLimits;

const GAMMA: f64 = 2000.0;
const SIZE: [f64; 3] = [0.002; 3];

fn two_coordinates<'c>(
    mesh: &'c TetComplex,
    pos: &[[f64; 3]],
    signs: [f64; 2],
    rates: [f64; 2],
) -> ThermoviscousFe<'c> {
    let (medium, thermo) = liquid();
    let geometry = HelmholtzFe::new(
        mesh,
        pos,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let amplitude = 1.0 / det::sqrt(medium.density_kg_m3 * SIZE.iter().product::<f64>());
    // Two independent mass-normalized source coordinates share one physical
    // boundary trace. Their differential motion is internal to the structure;
    // their common motion accelerates the actual resolved fluid.
    let mut shapes = vec![vec![[0.0; 3]; SimplexSpace::new(mesh, 2).ndof]; 2];
    for r in 0..2 {
        for face in &faces {
            for &v in face {
                shapes[r][v as usize][0] = signs[r] * amplitude;
            }
        }
    }
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            };
            2
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
    .with_damping_rates(&rates)
    .unwrap();
    ThermoviscousFe::new_with_moving_wall(
        mesh,
        pos,
        medium,
        thermo,
        &ThermoviscousBoundary {
            no_slip_faces: faces.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        MovingWall {
            faces,
            velocity_shapes: shapes,
            structure,
        },
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn relative_channel(signs: [f64; 2]) -> Vec<Vec<f64>> {
    let root = det::sqrt(GAMMA);
    vec![vec![signs[0] * root, -signs[1] * root]]
}

fn complex_close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-25),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn g1_full_relative_damping_keeps_common_fluid_loading_and_undamped_motion() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    let model = two_coordinates(&mesh, &pos, [1.0; 2], [0.0; 2])
        .with_moving_wall_damping_channels(relative_channel([1.0; 2]), &mut || false)
        .unwrap();
    let matrix = model.moving_wall_damping_matrix_per_s().unwrap();
    assert_eq!(matrix[0], -matrix[1]);
    assert_eq!(matrix[0], matrix[3]);
    let frequency = 30.0;
    let omega = std::f64::consts::TAU * frequency;
    let prepared = model.prepare(frequency, &mut || false).unwrap();
    let common = prepared
        .response_with_wall_forces(&[], &[C64::ONE; 2], &mut || false)
        .unwrap();
    let velocities = &common.moving_wall.modal_velocity_m_sqrt_kg_per_s;
    complex_close(velocities[0], velocities[1], 2e-10);
    // Fluid modal inertia is [[1,1],[1,1]]: common acceleration has total
    // modal mass 3, while the relative damper's force and heat are zero.
    complex_close(velocities[0], C64::new(0.0, -1.0 / (3.0 * omega)), 2e-5);
    assert!(common.moving_wall.fluid_force_n_per_sqrt_kg[0].abs() > 0.6);
    let diagonal_heat = GAMMA * velocities.iter().map(|v| v.norm_sq()).sum::<f64>();
    assert!(common.power.structural_damping_w < 1e-16 * diagonal_heat);
    let apparent = velocities.iter().map(|v| 0.5 * v.abs()).sum::<f64>();
    assert!(common.power.balance_residual_w().abs() < 2e-8 * apparent);

    let relative = prepared
        .response_with_wall_forces(&[], &[C64::ONE, -C64::ONE], &mut || false)
        .unwrap();
    let exact = C64::ONE / C64::new(2.0 * GAMMA, omega);
    complex_close(
        relative.moving_wall.modal_velocity_m_sqrt_kg_per_s[0],
        exact,
        2e-9,
    );
    complex_close(
        relative.moving_wall.modal_velocity_m_sqrt_kg_per_s[1],
        -exact,
        2e-9,
    );
    let expected_heat = 2.0 * GAMMA * exact.norm_sq();
    assert!((relative.power.structural_damping_w / expected_heat - 1.0).abs() < 2e-9);
    assert!(
        relative
            .moving_wall
            .fluid_force_n_per_sqrt_kg
            .iter()
            .all(|f| f.abs() < 2e-8)
    );
    assert_work(&relative);
}

#[test]
fn g1_causal_full_damping_matches_the_common_mass_and_differential_dashpot() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    let model = two_coordinates(&mesh, &pos, [1.0; 2], [0.0; 2])
        .with_moving_wall_damping_channels(relative_channel([1.0; 2]), &mut || false)
        .unwrap();
    let dt = 1e-4;
    let prepared = model
        .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();

    // Initialize the exact continuum pressure gradient for common acceleration.
    // P1 pressure and P2 velocity represent this solution without a spatial fit.
    let mut initial = prepared.zero_state();
    let amplitude = 1.0 / det::sqrt(1000.0 * SIZE.iter().product::<f64>());
    for (p, xyz) in initial.pressure_pa.iter_mut().zip(&pos) {
        *p = -1000.0 * amplitude * (2.0 / 3.0) * (xyz[0] - 0.5 * SIZE[0]);
    }
    let original = initial.clone();
    let common = prepared
        .step(&initial, &[], &[1.0; 2], &mut || false)
        .unwrap();
    assert_eq!(initial, original);
    for r in 0..2 {
        assert!((common.state.wall_velocity_m_sqrt_kg_per_s[r] / (dt / 3.0) - 1.0).abs() < 2e-8);
        assert!((common.state.wall_displacement_m_sqrt_kg[r] / (dt * dt / 6.0) - 1.0).abs() < 2e-8);
    }
    assert!(common.work.structural_damping_heat_j.abs() < 1e-14 * common.work.structural_work_j);
    assert!(common.work.balance_residual_j().abs() < 2e-8 * common.work.supplied_j());

    let relative = prepared
        .step(&prepared.zero_state(), &[], &[1.0, -1.0], &mut || false)
        .unwrap();
    // Exact internal mode: v'=1-2*gamma*v; no fluid boundary moves.
    let rate = 2.0 * GAMMA;
    let velocity = -det::expm1(-rate * dt) / rate;
    let displacement = (dt - velocity) / rate;
    let heat = 2.0 * displacement - velocity * velocity;
    for (r, sign) in [1.0, -1.0].into_iter().enumerate() {
        assert!(
            (relative.state.wall_velocity_m_sqrt_kg_per_s[r] / (sign * velocity) - 1.0).abs()
                < 2e-8
        );
        assert!(
            (relative.state.wall_displacement_m_sqrt_kg[r] / (sign * displacement) - 1.0).abs()
                < 2e-8
        );
    }
    assert!((relative.work.structural_damping_heat_j / heat - 1.0).abs() < 2e-8);
    assert!(relative.work.fluid_viscous_heat_j.abs() < 1e-14 * heat);
    assert!(relative.work.balance_residual_j().abs() < 2e-8 * relative.work.supplied_j());
}

#[test]
fn g3_signed_coordinate_changes_preserve_full_damping_and_acoustic_fields() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    let force = [C64::new(0.8, 0.3), C64::new(-0.2, 0.4)];
    let original = two_coordinates(&mesh, &pos, [1.0; 2], [0.0; 2])
        .with_moving_wall_damping_channels(relative_channel([1.0; 2]), &mut || false)
        .unwrap();
    let reversed = two_coordinates(&mesh, &pos, [-1.0, 1.0], [0.0; 2])
        .with_moving_wall_damping_channels(relative_channel([-1.0, 1.0]), &mut || false)
        .unwrap();
    let a = original
        .prepare(250.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &force, &mut || false)
        .unwrap();
    let b = reversed
        .prepare(250.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[-force[0], force[1]], &mut || false)
        .unwrap();
    for (r, sign) in [-1.0, 1.0].into_iter().enumerate() {
        complex_close(
            b.moving_wall.modal_velocity_m_sqrt_kg_per_s[r],
            a.moving_wall.modal_velocity_m_sqrt_kg_per_s[r].scale(sign),
            2e-9,
        );
    }
    let scale = a.pressure_pa.iter().map(|p| p.abs()).fold(0.0, f64::max);
    for (x, y) in a.pressure_pa.iter().zip(&b.pressure_pa) {
        assert!((*x - *y).abs() < 2e-9 * scale);
    }
    assert!((a.power.structural_damping_w / b.power.structural_damping_w - 1.0).abs() < 2e-9);
    assert_work(&a);
    assert_work(&b);
}

#[test]
fn g0_channels_require_one_finite_complete_structural_loss_owner() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    for channels in [
        vec![vec![1.0]],
        vec![vec![1.0, f64::NAN]],
        vec![vec![1e200; 2]],
        vec![vec![1e-200; 2]],
        vec![vec![1.0; 2]; 65],
    ] {
        assert!(
            two_coordinates(&mesh, &pos, [1.0; 2], [0.0; 2])
                .with_moving_wall_damping_channels(channels, &mut || false)
                .is_err()
        );
    }
    assert!(
        two_coordinates(&mesh, &pos, [1.0; 2], [1.0, 0.0])
            .with_moving_wall_damping_channels(relative_channel([1.0; 2]), &mut || false)
            .is_err()
    );
    let zero = two_coordinates(&mesh, &pos, [1.0; 2], [0.0; 2])
        .with_moving_wall_damping_channels(vec![], &mut || false)
        .unwrap();
    assert_eq!(
        zero.moving_wall_damping_matrix_per_s(),
        Some([0.0; 4].as_slice())
    );
    assert!(
        zero.with_moving_wall_damping_channels(vec![], &mut || false)
            .is_err()
    );
    assert!(matches!(
        two_coordinates(&mesh, &pos, [1.0; 2], [0.0; 2])
            .with_moving_wall_damping_channels(relative_channel([1.0; 2]), &mut || true),
        Err(ThermoviscousError::Cancelled),
    ));
}
