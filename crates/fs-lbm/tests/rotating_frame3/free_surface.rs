//! The original sparse VOF state consumes the same rigid-frame source law.

use super::*;
use fs_lbm::d3q19::freesurface3::{
    ContactModel3, ConversionRedistribution3, FreeSurface3, FreeSurfaceError3,
    FreeSurfaceInitialCell3,
};

fn seed(
    cells: [usize; 3],
    gravity: [f64; 3],
    initial: Vec<FreeSurfaceInitialCell3>,
) -> FreeSurface3 {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0xacc3,
                kernel_id: 0,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        FreeSurface3::from_cells(
            cells,
            0.9,
            gravity,
            0.0,
            ContactModel3::Neutral,
            &[],
            &initial,
            &cx,
        )
        .unwrap()
    })
}

fn uniform(cells: [usize; 3], rho: f64, velocity: [f64; 3]) -> Vec<FreeSurfaceInitialCell3> {
    let mut initial = Vec::new();
    for z in 0..cells[2] {
        for y in 0..cells[1] {
            for x in 0..cells[0] {
                initial.push(FreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    density: rho,
                    velocity,
                });
            }
        }
    }
    initial
}

// One tile's documented q-major state snapshot allows an independent physical
// first/second-moment check, including nonequilibrium after a real native step.
fn single_tile_populations(surface: &FreeSurface3, [x, y, z]: [usize; 3]) -> [f64; Q3] {
    assert_eq!(surface.domain_cells(), [4; 3]);
    let bits = surface.grid().state_bits();
    let lane = x + 4 * y + 16 * z;
    core::array::from_fn(|q| f64::from_bits(bits[q * 64 + lane]))
}

#[test]
fn sparse_coriolis_source_retains_physical_moments_and_has_no_midpoint_work() {
    let gate = CancelGate::new_clock_free();
    let mut initial = uniform([4; 3], 1.3, [0.03, -0.015, 0.021]);
    initial
        .iter_mut()
        .find(|cell| cell.coordinates == [1, 2, 1])
        .unwrap()
        .fill_fraction = 0.6;
    let mut surface = seed([4; 3], [0.0002, -0.0001, 0.0003], initial);
    let cell = [1, 2, 1];
    let original = surface.sample_at(cell[0], cell[1], cell[2]).unwrap();
    let old_force = surface.force_density_at(cell[0], cell[1], cell[2]).unwrap();
    let old_moments = moments(single_tile_populations(&surface, cell), old_force);
    let mass = surface.ledger_mass();
    let omega = [0.04, -0.03, 0.02];
    RigidFrame3 {
        origin: [1.5, 2.5, 1.5],
        angular_velocity: omega,
        ..Default::default()
    }
    .prepare_free_surface_step(&mut surface, [0.0; 3], &gate)
    .unwrap();
    let next = surface.sample_at(cell[0], cell[1], cell[2]).unwrap();
    assert_eq!(next.density, original.density);
    assert_eq!(next.fill_fraction, original.fill_fraction);
    assert_eq!(surface.ledger_mass(), mass);
    let force = surface.force_density_at(cell[0], cell[1], cell[2]).unwrap();
    let preserved = moments(single_tile_populations(&surface, cell), force);
    // Read the actual shared collision result, not only the force descriptor.
    let post = collide_cell3(
        single_tile_populations(&surface, cell),
        surface.collision_model(),
        force,
    )
    .unwrap();
    let after_collision = moments(post, force);
    let velocity = after_collision.1.map(|p| p / after_collision.0);
    let magnitude = dot(omega, omega).sqrt();
    let axis = omega.map(|w| w / magnitude);
    let angle = 2.0 * magnitude.atan();
    let parallel = dot(axis, original.velocity);
    let turn = cross(axis, original.velocity);
    for a in 0..3 {
        near(next.velocity[a], original.velocity[a], 3e-15);
        near(preserved.1[a], old_moments.1[a], 4e-15);
        for b in 0..3 {
            near(preserved.2[a][b], old_moments.2[a][b], 4e-15);
        }
        let expected = axis[a] * parallel
            + angle.cos() * (original.velocity[a] - axis[a] * parallel)
            - angle.sin() * turn[a];
        near(velocity[a], expected, 4e-15);
    }
    near(
        dot(velocity, velocity),
        dot(original.velocity, original.velocity),
        4e-16,
    );
    let midpoint = core::array::from_fn(|d| 0.5 * (velocity[d] + original.velocity[d]));
    near(dot(force, midpoint), 0.0, 4e-16);
}

#[test]
fn sparse_force_history_produces_density_independent_acceleration_and_free_fall() {
    let gate = CancelGate::new_clock_free();
    let body = [0.0002, -0.0003, 0.0001];
    for rho in [0.7, 1.4] {
        let original = seed([12; 3], body, uniform([12; 3], rho, [0.0; 3]));
        for origin in [body, [-0.0001, 0.0001, 0.0]] {
            let mut surface = original.clone();
            let frame = RigidFrame3 {
                origin_acceleration: origin,
                ..Default::default()
            };
            for step in 1..=3 {
                let before = surface.sample_at(6, 6, 6).unwrap();
                let mass = surface.ledger_mass();
                frame
                    .prepare_free_surface_step(&mut surface, body, &gate)
                    .unwrap();
                let prepared = surface.sample_at(6, 6, 6).unwrap();
                assert_eq!(surface.ledger_mass(), mass);
                assert_eq!(prepared.density, before.density);
                for d in 0..3 {
                    near(prepared.velocity[d], before.velocity[d], 4e-15);
                    near(
                        surface.force_density_at(6, 6, 6).unwrap()[d],
                        before.density * (body[d] - origin[d]),
                        2e-17,
                    );
                }
                surface.step_with_gate(&gate).unwrap();
                let after = surface.sample_at(6, 6, 6).unwrap();
                for d in 0..3 {
                    near(
                        after.velocity[d],
                        f64::from(step) * (body[d] - origin[d]),
                        2e-14,
                    );
                }
            }
        }
    }
}

#[test]
fn sparse_euler_and_centrifugal_terms_use_half_cell_positions_and_current_density() {
    let gate = CancelGate::new_clock_free();
    let frame = RigidFrame3 {
        origin: [1.5, 1.5, 2.5],
        angular_velocity: [0.01, -0.02, 0.03],
        angular_acceleration: [0.002, 0.003, -0.001],
        origin_acceleration: [0.001, -0.002, 0.005],
    };
    // At (2.5,3.5,1.5), centrifugal=(.002,.0016,.0004), Euler=(.001,-.001,-.001).
    let expected = [0.006, -0.0034, -0.0026];
    let velocity = expected.map(|a| -0.5 * a); // zero midpoint, hence zero Coriolis.
    let mut surface = seed([4; 3], [0.0; 3], uniform([4; 3], 1.25, velocity));
    frame
        .prepare_free_surface_step(&mut surface, [0.004, -0.006, 0.003], &gate)
        .unwrap();
    let actual = surface.force_density_at(2, 3, 1).unwrap();
    let sample = surface.sample_at(2, 3, 1).unwrap();
    for d in 0..3 {
        near(actual[d], sample.density * expected[d], 2e-16);
    }
}

#[test]
fn conversion_carries_the_force_frame_and_refreshes_new_wet_cells() {
    let gate = CancelGate::new_clock_free();
    let initial = uniform([4; 3], 1.0, [0.04, 0.0, 0.0])
        .into_iter()
        .map(|mut cell| {
            for d in 0..3 {
                cell.coordinates[d] += 4;
            }
            cell
        })
        .collect();
    let mut surface = seed([12; 3], [0.0; 3], initial);
    let previous = surface.clone();
    let frame = RigidFrame3::default();
    let acceleration = [0.0002, -0.0001, 0.00005];
    frame
        .prepare_free_surface_step(&mut surface, acceleration, &gate)
        .unwrap();
    surface.step_with_gate(&gate).unwrap();
    assert!(surface.conversions().gas_to_interface > 0);
    let mut newly_wet = 0;
    for z in 0..12 {
        for y in 0..12 {
            for x in 0..12 {
                if previous.sample_at(x, y, z).is_none() && surface.sample_at(x, y, z).is_some() {
                    newly_wet += 1;
                    assert!(
                        surface
                            .force_density_at(x, y, z)
                            .unwrap()
                            .into_iter()
                            .all(f64::is_finite)
                    );
                }
            }
        }
    }
    assert!(newly_wet > 0);
    let mass = surface.ledger_mass();
    frame
        .prepare_free_surface_step(&mut surface, acceleration, &gate)
        .unwrap();
    assert_eq!(surface.ledger_mass(), mass);
    for z in 0..12 {
        for y in 0..12 {
            for x in 0..12 {
                if let Some(sample) = surface.sample_at(x, y, z) {
                    let force = surface.force_density_at(x, y, z).unwrap();
                    for d in 0..3 {
                        near(force[d], sample.density * acceleration[d], 1e-18);
                    }
                }
            }
        }
    }
}

#[test]
fn invalid_frame_sources_and_cancellation_leave_sparse_state_unchanged() {
    let gate = CancelGate::new_clock_free();
    let mut surface = seed(
        [4; 3],
        [0.001, -0.002, 0.003],
        uniform([4; 3], 1.0, [0.01; 3]),
    );
    let before = surface.grid().state_bits();
    let mass = surface.ledger_mass();
    let force = surface.force_density_at(1, 2, 3);
    for (frame, acceleration) in [
        (RigidFrame3::default(), [f64::NAN, 0.0, 0.0]),
        (
            RigidFrame3 {
                angular_velocity: [f64::MAX, 0.0, 0.0],
                ..Default::default()
            },
            [0.0; 3],
        ),
        (
            RigidFrame3 {
                origin: [f64::MAX; 3],
                angular_acceleration: [f64::MAX; 3],
                ..Default::default()
            },
            [0.0; 3],
        ),
    ] {
        assert!(
            frame
                .prepare_free_surface_step(&mut surface, acceleration, &gate)
                .is_err()
        );
        assert_eq!(surface.grid().state_bits(), before);
        assert_eq!(surface.ledger_mass(), mass);
        assert_eq!(surface.force_density_at(1, 2, 3), force);
    }
    let mut historical = surface
        .clone()
        .with_conversion_redistribution(ConversionRedistribution3::LegacyGlobal);
    assert!(
        RigidFrame3::default()
            .prepare_free_surface_step(&mut historical, [0.0; 3], &gate)
            .is_err()
    );
    assert_eq!(historical.grid().state_bits(), before);
    for acceleration in [[0.0; 3], [0.0003, -0.0002, 0.0001]] {
        let mut independent = surface.clone();
        independent
            .rebind_body_acceleration(acceleration, &gate)
            .unwrap();
        let populations = independent.grid().state_bits();
        let mass = independent.ledger_mass();
        let force = independent.force_density_at(1, 2, 3);
        assert_eq!(
            RigidFrame3::default().prepare_free_surface_step(&mut independent, [0.0; 3], &gate),
            Err(FreeSurfaceError3::Transport(
                "rigid-frame forcing cannot replace an independent body-acceleration owner",
            ))
        );
        assert_eq!(independent.grid().state_bits(), populations);
        assert_eq!(independent.ledger_mass(), mass);
        assert_eq!(independent.force_density_at(1, 2, 3), force);
        assert_eq!(independent.body_acceleration(), Some(acceleration));
    }
    gate.request();
    assert_eq!(
        RigidFrame3::default().prepare_free_surface_step(&mut surface, [0.0; 3], &gate),
        Err(FreeSurfaceError3::Cancelled)
    );
    assert_eq!(surface.grid().state_bits(), before);
    assert_eq!(surface.force_density_at(1, 2, 3), force);
}
