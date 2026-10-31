use super::*;

fn tank(wall: PlanarWall3) -> FreeSurface3 {
    FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, y, z| {
            let p = [x, y, z][wall.axis()];
            if wall.upper() { p >= 4 } else { p < 4 }
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn liquid_momentum(surface: &FreeSurface3) -> [f64; 3] {
    let mut result = [0.0; 3];
    for (&key, tile) in &surface.tiles {
        let slot = surface.grid.slot_of(key).unwrap();
        for (lane, &cell) in tile.cells.iter().enumerate() {
            if cell == Cell3::Gas {
                continue;
            }
            let p = FreeSurface3::global_of(key, lane);
            let state = surface
                .sample_at(p.0 as usize, p.1 as usize, p.2 as usize)
                .unwrap();
            let mass = if cell == Cell3::Fluid {
                surface.grid.populations(slot, lane).iter().sum()
            } else {
                tile.mass[lane]
            };
            for axis in 0..3 {
                result[axis] += mass * state.velocity[axis];
            }
        }
    }
    result
}

#[test]
fn g1_six_flooded_moving_planes_exchange_actual_populations_mass_impulse_and_work() {
    for wall in [
        PlanarWall3::XMin,
        PlanarWall3::XMax,
        PlanarWall3::YMin,
        PlanarWall3::YMax,
        PlanarWall3::ZMin,
        PlanarWall3::ZMax,
    ] {
        let mut surface = tank(wall);
        let before = liquid_momentum(&surface);
        let tangent = (wall.axis() + 1) % 3;
        let mut velocity = [0.0; 3];
        velocity[tangent] = 0.01;
        let report = surface
            .step_with_planar_wall_motion(wall, velocity, 1e-12, &CancelGate::new_clock_free())
            .unwrap();
        let after = liquid_momentum(&surface);
        // Independently sum the two moving tangential D3Q19 diagonals:
        // 2*(2/36)*rho*u/cs² = rho*u/3 per fully wetted wall cell.
        let impulse = 64.0 * velocity[tangent] / 3.0;
        close(report.fluid_impulse[tangent], impulse);
        close(after[tangent] - before[tangent], impulse);
        close(report.fluid_work, velocity[tangent] * impulse);
        close(report.wall_work, -report.fluid_work);
        close(report.wetted_area, 64.0);
        let inward = if wall.upper() { -1.0 } else { 1.0 };
        close(report.fluid_impulse[wall.axis()], inward * 64.0 / 3.0);
        for axis in 0..3 {
            assert_eq!(report.wall_impulse[axis], -report.fluid_impulse[axis]);
        }
        close(report.mass_after, report.mass_before);
        assert!(report.relative_mass_defect < 1e-12);
        let mut p = [4; 3];
        p[wall.axis()] = if wall.upper() { 7 } else { 0 };
        close(
            surface.sample_at(p[0], p[1], p[2]).unwrap().velocity[tangent],
            velocity[tangent] / 3.0,
        );
        // The free surface still belongs to the original VOF state.
        assert!(
            surface
                .tiles
                .values()
                .any(|tile| tile.cells.contains(&Cell3::Interface))
        );
    }
}

#[test]
fn g3_zero_wall_motion_preserves_native_bits_and_rebound_motion_replays() {
    let wall = PlanarWall3::ZMin;
    let mut moving = tank(wall);
    let mut original = moving.clone();
    original.step().unwrap();
    let zero = moving
        .step_with_planar_wall_motion(wall, [0.0; 3], 1e-12, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(moving.grid.state_bits(), original.grid.state_bits());
    assert_eq!(moving.mass_snapshot(), original.mass_snapshot());
    assert_eq!(moving.conversions(), original.conversions());
    assert_eq!(zero.fluid_work, 0.0);
    let mut replay = moving.clone();
    for velocity in [
        [0.005, 0.0, 0.0],
        [-0.005, 0.0, 0.0],
        [0.0; 3],
        [0.005, 0.0, 0.0],
    ] {
        let a = moving
            .step_with_planar_wall_motion(wall, velocity, 1e-12, &CancelGate::new_clock_free())
            .unwrap();
        let b = replay
            .step_with_planar_wall_motion(wall, velocity, 1e-12, &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(moving.grid.state_bits(), replay.grid.state_bits());
        assert_eq!(moving.mass_snapshot(), replay.mass_snapshot());
    }
    assert_ne!(moving.grid.state_bits(), original.grid.state_bits());
}

#[test]
fn g1_slow_wall_receives_mechanical_work_from_faster_liquid() {
    let wall = PlanarWall3::ZMin;
    let mut surface = tank(wall);
    for (&key, tile) in &surface.tiles {
        let slot = surface.grid.slot_of(key).unwrap();
        for (lane, &cell) in tile.cells.iter().enumerate() {
            if cell != Cell3::Gas {
                surface
                    .grid
                    .set_populations(slot, lane, equilibrium3(1.0, [0.02, 0.0, 0.0]));
            }
        }
    }
    let report = surface
        .step_with_planar_wall_motion(wall, [0.01, 0.0, 0.0], 1e-12, &CancelGate::new_clock_free())
        .unwrap();
    close(report.fluid_impulse[0], 64.0 * (0.01 - 0.02) / 3.0);
    assert!(report.fluid_work < 0.0 && report.wall_work > 0.0);
    close(report.wall_work, 64.0 * 0.01 * (0.02 - 0.01) / 3.0);
}

#[test]
fn g0_g4_contact_normal_relative_mach_and_cancel_refuse_without_publication() {
    let wall = PlanarWall3::ZMin;
    let mut surface = tank(wall);
    let before = surface.grid.state_bits();
    let mass = surface.mass_snapshot();
    for velocity in [[0.0, 0.0, 0.001], [f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0]] {
        assert!(
            surface
                .step_with_planar_wall_motion(wall, velocity, 1e-12, &CancelGate::new_clock_free())
                .is_err()
        );
    }
    assert!(
        surface
            .step_with_planar_wall_motion(wall, [0.0; 3], -1.0, &CancelGate::new_clock_free())
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        surface.step_with_planar_wall_motion(wall, [0.01, 0.0, 0.0], 1e-12, &cancelled),
        Err(FreeSurfaceError3::Cancelled)
    );
    assert_eq!(surface.grid.state_bits(), before);
    assert_eq!(surface.mass_snapshot(), mass);

    let mut contact = FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, z| z == 0,
    )
    .unwrap();
    let before = contact.grid.state_bits();
    assert!(
        matches!(contact.step_with_planar_wall_motion(wall,[0.01,0.0,0.0],1e-12,
        &CancelGate::new_clock_free()),
        Err(FreeSurfaceError3::InvalidWallMotion{reason,..}) if reason.contains("partial wetting"))
    );
    assert_eq!(contact.grid.state_bits(), before);

    // Each inertial velocity is low-Mach, but their relative speed is not.
    let fast = 0.75 * MACH_LIMIT * CS2_3.sqrt();
    for (&key, tile) in &surface.tiles {
        let slot = surface.grid.slot_of(key).unwrap();
        for (lane, &cell) in tile.cells.iter().enumerate() {
            if cell != Cell3::Gas {
                surface
                    .grid
                    .set_populations(slot, lane, equilibrium3(1.0, [-fast, 0.0, 0.0]));
            }
        }
    }
    let before = surface.grid.state_bits();
    assert!(
        matches!(surface.step_with_planar_wall_motion(wall,[fast,0.0,0.0],1e-12,
        &CancelGate::new_clock_free()),
        Err(FreeSurfaceError3::InvalidWallMotion{reason,..}) if reason.contains("relative"))
    );
    assert_eq!(surface.grid.state_bits(), before);
}

#[test]
fn g4_late_endpoint_mach_refusal_preserves_original_populations_mass_and_counters() {
    let wall = PlanarWall3::ZMin;
    let mut surface = tank(wall);
    // Zero initial speeds pass; a deliberately sharp pressure pulse creates
    // an excessive velocity only after the ORIGINAL collision/streaming step.
    let (slot, lane) = surface.grid.resolve_source(4, 4, 1).unwrap();
    surface
        .grid
        .set_populations(slot, lane, equilibrium3(10.0, [0.0; 3]));
    let before = surface.grid.state_bits();
    let mass = surface.mass_snapshot();
    let conversions = surface.conversions();
    assert!(
        matches!(surface.step_with_planar_wall_motion(wall,[0.005,0.0,0.0],1e-12,
        &CancelGate::new_clock_free()),
        Err(FreeSurfaceError3::InvalidWallMotion{reason,..}) if reason.contains("Mach"))
    );
    assert_eq!(surface.grid.state_bits(), before);
    assert_eq!(surface.mass_snapshot(), mass);
    assert_eq!(surface.conversions(), conversions);
}
