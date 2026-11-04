//! Physical acceleration and hydrostatics in the original sparse liquid owner.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurface3, FreeSurfaceInitialCell3};

const SIDE: usize = 8;

fn liquid(
    tau: f64,
    density: impl Fn([usize; 3]) -> Option<f64>,
    velocity: [f64; 3],
) -> FreeSurface3 {
    let mut initial = Vec::new();
    for z in 0..SIDE {
        for y in 0..SIDE {
            for x in 0..SIDE {
                let coordinates = [x, y, z];
                if let Some(density) = density(coordinates) {
                    initial.push(FreeSurfaceInitialCell3 {
                        coordinates,
                        fill_fraction: 1.0,
                        density,
                        velocity,
                    });
                }
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xACC3,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        FreeSurface3::from_cells(
            [SIDE; 3],
            tau,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            &[],
            &initial,
            &cx,
        )
        .unwrap()
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tolerance {tolerance:e})"
    );
}

#[test]
fn equal_acceleration_at_unequal_densities_survives_load_switches() {
    let gate = CancelGate::new_clock_free();
    let initial_velocity = [0.02, -0.01, 0.015];
    let acceleration = [0.0007, -0.0012, 0.0009];
    for density in [0.5, 0.8, 1.3, 2.4] {
        let mut flow = liquid(0.8, |_| Some(density), initial_velocity);
        assert_eq!(flow.body_acceleration(), None);
        let mut expected_velocity = initial_velocity;
        // Three steps cannot carry a wall reflection to this interior cell.
        for acceleration in [acceleration, acceleration.map(|a| -0.5 * a), [0.0; 3]] {
            flow.rebind_body_acceleration(acceleration, &gate).unwrap();
            assert_eq!(flow.body_acceleration(), Some(acceleration));
            let rebound = flow.sample_at(4, 4, 4).unwrap();
            for axis in 0..3 {
                // Replacing a held load has no elapsed time and no impulse.
                near(rebound.velocity[axis], expected_velocity[axis], 3e-14);
            }
            flow.step_with_gate(&gate).unwrap();
            let after = flow.sample_at(4, 4, 4).unwrap();
            let force = flow.force_density_at(4, 4, 4).unwrap();
            near(after.density, density, 8e-15);
            for axis in 0..3 {
                expected_velocity[axis] += acceleration[axis];
                near(after.velocity[axis], expected_velocity[axis], 3e-14);
                near(force[axis], after.density * acceleration[axis], 1e-17);
            }
        }
    }
}

#[test]
fn affine_density_uses_the_endpoint_force_without_a_fictitious_reframe() {
    let gate = CancelGate::new_clock_free();
    let density = 1.1;
    let gradient = [0.01, -0.006, 0.008];
    let acceleration = [0.004, -0.003, 0.002];
    let density_change = -0.5
        * gradient
            .iter()
            .zip(acceleration)
            .map(|(g, a)| g * a)
            .sum::<f64>();
    let endpoint_density = density + density_change;

    // At u=0, the initialized Guo populations collide to
    // w_q*rho*(1 + 3 e_q.a/2), independent of relaxation time. For affine
    // rho, the lattice's second moment then gives rho_1=rho_0-grad(rho).a/2
    // and j_raw,1=-grad(rho)/3+rho_0*a/2. Adding the ENDPOINT half force
    // produces the trapezoidal physical impulse (rho_0+rho_1)*a/2.
    // Freezing the initial force or reframing the streamed populations loses
    // that density-dependent half impulse; this oracle uses neither kernel.
    for tau in [0.63, 0.9, 1.25] {
        let mut flow = liquid(
            tau,
            |p| {
                Some(
                    density
                        + (0..3)
                            .map(|axis| gradient[axis] * (p[axis] as f64 - 4.0))
                            .sum::<f64>(),
                )
            },
            [0.0; 3],
        );
        flow.rebind_body_acceleration(acceleration, &gate).unwrap();
        flow.step_with_gate(&gate).unwrap();
        let after = flow.sample_at(4, 4, 4).unwrap();
        let force = flow.force_density_at(4, 4, 4).unwrap();
        assert!((after.density - density).abs() > 3e-5);
        near(after.density, endpoint_density, 3e-15);
        for axis in 0..3 {
            let momentum =
                -gradient[axis] / 3.0 + 0.5 * (density + endpoint_density) * acceleration[axis];
            near(after.density * after.velocity[axis], momentum, 3e-15);
            near(force[axis], after.density * acceleration[axis], 1e-17);
        }
    }
}

#[test]
fn density_proportional_hydrostatics_balance_atmosphere_and_walls() {
    let gate = CancelGate::new_clock_free();
    let gravity = 0.003;
    // Exact discrete hydrostatics, not a constant-force linear density ramp:
    // opposing streamed populations balance when adjacent densities have this
    // ratio. At the free surface, rho_surface*(1-3g/2)=rho_atmosphere=1.
    let ratio: f64 = (1.0 + 1.5 * gravity) / (1.0 - 1.5 * gravity);
    let surface_density = 1.0 / (1.0 - 1.5 * gravity);
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let density = |p: [usize; 3]| {
                let depth = if sign < 0.0 {
                    3usize.checked_sub(p[axis])?
                } else {
                    p[axis].checked_sub(4)?
                };
                Some(surface_density * ratio.powi(depth as i32))
            };
            let mut acceleration = [0.0; 3];
            acceleration[axis] = sign * gravity;
            let mut flow = liquid(0.8, density, [0.0; 3]);
            flow.rebind_body_acceleration(acceleration, &gate).unwrap();
            let mass = flow.ledger_mass();
            for _ in 0..16 {
                flow.step_with_gate(&gate).unwrap();
            }
            for z in 0..SIDE {
                for y in 0..SIDE {
                    for x in 0..SIDE {
                        let Some(expected_density) = density([x, y, z]) else {
                            continue;
                        };
                        let sample = flow.sample_at(x, y, z).unwrap();
                        near(sample.density, expected_density, 2e-12);
                        near(sample.fill_fraction, 1.0, 2e-12);
                        for velocity in sample.velocity {
                            near(velocity, 0.0, 2e-12);
                        }
                    }
                }
            }
            near(flow.ledger_mass(), mass, 2e-12 * mass);
        }
    }
}
