use super::*;
use crate::compressible::Error;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}",
    );
}

fn exact(
    omega: [f64; 3],
    acceleration: [f64; 3],
    density: f64,
    momentum: [f64; 3],
    energy: f64,
    duration: f64,
) -> Kick<3> {
    kick::<Error, 3>(
        BodyForce::CoriolisExact {
            angular_velocity: &[omega],
            acceleration: &[acceleration],
        },
        0,
        density,
        momentum,
        energy,
        duration,
    )
    .unwrap()
}

// Independent Rodrigues oracle for the axis (1,2,2)/3. Its exact norm is
// known, so the physical phase does not depend on production normalization.
fn rotated(momentum: [f64; 3], duration: f64) -> [f64; 3] {
    let axis = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let projection = (momentum[0] + 2.0 * momentum[1] + 2.0 * momentum[2]) / 3.0;
    let cross = [
        (2.0 * momentum[2] - 2.0 * momentum[1]) / 3.0,
        (2.0 * momentum[0] - momentum[2]) / 3.0,
        (momentum[1] - 2.0 * momentum[0]) / 3.0,
    ];
    let (sine, cosine) = (-6.0 * duration).sin_cos();
    core::array::from_fn(|a| {
        cosine * momentum[a] + sine * cross[a] + (1.0 - cosine) * projection * axis[a]
    })
}

#[test]
fn exact_coriolis_preserves_physical_phase_speed_and_axial_momentum() {
    let initial: [f64; 3] = [0.23, -0.17, 0.11];
    let initial_norm = initial.into_iter().map(|p| p * p).sum::<f64>();
    for duration in [1e-12, 1e-5, 0.125, std::f64::consts::PI / 6.0, 1.9, 6.3] {
        let after = exact([1.0, 2.0, 2.0], [0.0; 3], 1.7, initial, 4.0, duration);
        let expected = rotated(initial, duration);
        assert_eq!(after.energy.to_bits(), 4.0_f64.to_bits());
        assert_eq!(after.work_density, 0.0);
        close(
            after.momentum.into_iter().map(|p| p * p).sum(),
            initial_norm,
            5e-15,
        );
        close(
            after.momentum[0] + 2.0 * after.momentum[1] + 2.0 * after.momentum[2],
            initial[0] + 2.0 * initial[1] + 2.0 * initial[2],
            5e-15,
        );
        for a in 0..3 {
            close(after.momentum[a], expected[a], 5e-15);
            assert_eq!(after.impulse_density[a], after.momentum[a] - initial[a]);
        }
    }

    // More than thirty physical turns: phase error must remain at roundoff,
    // rather than accumulating the midpoint/Cayley dispersion error.
    let mut momentum = initial;
    for _ in 0..256 {
        momentum = exact([1.0, 2.0, 2.0], [0.0; 3], 1.7, momentum, 4.0, 0.13).momentum;
    }
    let expected = rotated(initial, 256.0 * 0.13);
    for a in 0..3 {
        close(momentum[a], expected[a], 2e-13);
    }
}

#[test]
fn exact_coriolis_integrates_parallel_and_transverse_drive_without_heat() {
    let rho = 1.3;
    let velocity = [0.2, -0.1, 0.3];
    let momentum = velocity.map(|u| rho * u);
    let acceleration = [0.4, -0.6, 0.9];
    let duration = 0.83;
    let frequency = 3.5;
    let initial_internal = 50.0;
    let initial_kinetic = 0.5 * rho * velocity.into_iter().map(|u| u * u).sum::<f64>();
    let initial_energy = initial_internal + initial_kinetic;
    let after = exact(
        [0.0, 0.0, 0.5 * frequency],
        acceleration,
        rho,
        momentum,
        initial_energy,
        duration,
    );
    // Closed form for du_x/dt = f*u_y+a_x,
    // du_y/dt = -f*u_x+a_y and du_z/dt = a_z.
    let (sine, cosine) = (frequency * duration).sin_cos();
    let expected = [
        cosine * velocity[0]
            + sine * velocity[1]
            + (sine * acceleration[0] + (1.0 - cosine) * acceleration[1]) / frequency,
        -sine * velocity[0]
            + cosine * velocity[1]
            + (-(1.0 - cosine) * acceleration[0] + sine * acceleration[1]) / frequency,
        velocity[2] + duration * acceleration[2],
    ];
    let final_kinetic = 0.5 * rho * expected.into_iter().map(|u| u * u).sum::<f64>();
    for a in 0..3 {
        close(after.momentum[a], rho * expected[a], 2e-15);
    }
    close(after.work_density, final_kinetic - initial_kinetic, 3e-15);
    close(after.energy - final_kinetic, initial_internal, 3e-15);

    // The held affine source is a flow map: one full interval and two unequal
    // subintervals must agree, including drive, axial motion and work.
    for omega in [[0.0, 0.0, 1.75], [1.0, -2.0, 0.5], [-1.0, 0.0, 0.0]] {
        let full = exact(omega, acceleration, rho, momentum, initial_energy, duration);
        let first = exact(
            omega,
            acceleration,
            rho,
            momentum,
            initial_energy,
            0.3 * duration,
        );
        let second = exact(
            omega,
            acceleration,
            rho,
            first.momentum,
            first.energy,
            0.7 * duration,
        );
        for a in 0..3 {
            close(second.momentum[a], full.momentum[a], 4e-15);
            close(
                first.impulse_density[a] + second.impulse_density[a],
                full.impulse_density[a],
                4e-15,
            );
        }
        close(second.energy, full.energy, 3e-15);
        close(
            first.work_density + second.work_density,
            full.work_density,
            4e-15,
        );
    }
}

#[test]
fn exact_coriolis_retains_geostrophic_balance_and_zero_spin_arithmetic() {
    let momentum = [1.0, -0.5, 0.25];
    for duration in [1.0 / 64.0, 0.5, 4.0, 128.0] {
        let after = exact(
            [0.0, 0.0, 2.0],
            [1.0, 2.0, 0.0],
            2.0,
            momentum,
            8.0,
            duration,
        );
        assert_eq!(after.momentum, momentum);
        assert_eq!(after.energy.to_bits(), 8.0_f64.to_bits());
        assert_eq!(after.impulse_density, [0.0; 3]);
        assert_eq!(after.work_density, 0.0);
    }
    let acceleration = [0.4, -0.2, 0.1];
    let reference = kick::<Error, 3>(
        BodyForce::UniformAcceleration(acceleration),
        0,
        1.3,
        momentum,
        8.0,
        0.07,
    )
    .unwrap();
    let after = exact([0.0; 3], acceleration, 1.3, momentum, 8.0, 0.07);
    assert_eq!(after.momentum, reference.momentum);
    assert_eq!(after.energy.to_bits(), reference.energy.to_bits());
    assert_eq!(after.impulse_density, reference.impulse_density);
    assert_eq!(
        after.work_density.to_bits(),
        reference.work_density.to_bits()
    );
}

#[test]
fn exact_coriolis_admits_complete_fields_and_planar_geometry() {
    let rotation = [[0.0, 0.0, 1.0], [0.0; 3], [0.0, 0.0, -1.0]];
    let acceleration = [[0.4, -0.2]];
    let field = BodyForce::CoriolisExact {
        angular_velocity: &rotation,
        acceleration: &acceleration,
    };
    assert!(field.admit::<Error>(3, &mut || false).unwrap());
    for cell in 0..3 {
        let rho = 0.5 + cell as f64;
        let momentum = [0.2 * rho, -0.1 * rho];
        let actual = kick::<Error, 2>(field, cell, rho, momentum, 4.0, 0.1).unwrap();
        let uniform = BodyForce::CoriolisExact {
            angular_velocity: &rotation[cell..=cell],
            acceleration: &acceleration,
        };
        let expected = kick::<Error, 2>(uniform, 0, rho, momentum, 4.0, 0.1).unwrap();
        assert_eq!(actual.momentum, expected.momentum);
        assert_eq!(actual.energy, expected.energy);
    }
    for invalid in [
        BodyForce::CoriolisExact {
            angular_velocity: &[],
            acceleration: &acceleration,
        },
        BodyForce::CoriolisExact {
            angular_velocity: &rotation,
            acceleration: &[],
        },
        BodyForce::CoriolisExact {
            angular_velocity: &rotation[..2],
            acceleration: &acceleration,
        },
        BodyForce::CoriolisExact {
            angular_velocity: &[[1.0, 0.0, 1.0]],
            acceleration: &acceleration,
        },
        BodyForce::CoriolisExact {
            angular_velocity: &rotation,
            acceleration: &[[f64::NAN, 0.0]],
        },
    ] {
        assert!(matches!(
            invalid.admit::<Error>(3, &mut || false),
            Err(Error::Invalid(_))
        ));
    }
    assert!(matches!(
        field.admit::<Error>(3, &mut || true),
        Err(Error::Cancelled)
    ));
}
