use super::*;
use fs_exec::{Budget, Cx, ExecMode, StreamKey};

fn controls() -> SubmergedSolidLoadControls3 {
    SubmergedSolidLoadControls3 {
        maximum_solid_tiles: 8,
        maximum_boundary_links: 4_000,
        maximum_relative_mass_defect: 1e-11,
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-12, "{a:.17e} != {b:.17e}");
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// A linear pressure gradient balances a constant force density exactly on
/// the lattice. The free surface is optional; physical initial velocity uses
/// the original full Guo ingress, not a guessed raw population velocity.
fn tank(
    solids: &[(u32, u32, u32)],
    force: [f64; 3],
    velocity: [f64; 3],
    reference_density: f64,
    height: usize,
) -> FreeSurface3 {
    let mut cells = Vec::new();
    for z in 0..height {
        for y in 0..12 {
            for x in 0..20 {
                if solids.contains(&((x / TILE) as u32, (y / TILE) as u32, (z / TILE) as u32)) {
                    continue;
                }
                let p = [x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5];
                let rho = reference_density + (0..3).map(|a| force[a] * p[a] / CS2_3).sum::<f64>();
                cells.push(FreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    density: rho,
                    velocity,
                });
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        FreeSurface3::from_cells(
            [20, 12, 16],
            0.8,
            force,
            0.0,
            ContactModel3::Neutral,
            solids,
            &cells,
            &cx,
        )
        .unwrap()
    })
}

#[test]
fn archimedes_force_and_center_of_buoyancy_follow_actual_hydrostatic_populations() {
    let solid = [(1, 1, 1)];
    let origin = [1.25, -0.5, 2.75];
    let objects = [SubmergedSolid3 {
        solid_tiles: &solid,
        moment_origin: origin,
    }];
    let center = [1.5 * TILE as f64; 3];
    let volume = (TILE * TILE * TILE) as f64;
    for force in [
        [1e-5, 0.0, 0.0],
        [0.0, -2e-5, 0.0],
        [0.0, 0.0, -3e-5],
        [1e-5, -2e-5, -3e-5],
    ] {
        // A different uniform pressure datum must not alter net buoyancy.
        for rho in [0.9, 1.1] {
            let mut flow = tank(&solid, force, [0.0; 3], rho, 16);
            let report = flow
                .step_with_submerged_solid_loads(
                    &objects,
                    controls(),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            let expected = force.map(|g| -volume * g);
            let moment = cross(core::array::from_fn(|a| center[a] - origin[a]), expected);
            for a in 0..3 {
                close(report.bodies[0].wall_impulse[a], expected[a]);
                close(report.bodies[0].wall_angular_impulse[a], moment[a]);
            }
            assert_eq!(report.bodies[0].measured_links, 432);
            assert!(report.maximum_speed < 3e-12);
            assert!(report.relative_mass_defect < 1e-12);
            assert_eq!(report.geometry_epoch, flow.geometry_epoch());
        }
    }
}

#[test]
fn actual_free_surface_drag_keeps_the_unobserved_native_trajectory_bit_exact() {
    let solid = [(1, 1, 1)];
    let velocity = [0.003, -0.004, 0.002];
    let origin = [6.0; 3];
    let objects = [SubmergedSolid3 {
        solid_tiles: &solid,
        moment_origin: origin,
    }];
    let mut flow = tank(&solid, [0.0; 3], velocity, 1.0, 12);
    let mut ordinary = flow.clone();
    assert!(
        flow.tiles
            .values()
            .any(|t| t.cells.contains(&Cell3::Interface))
    );
    for step in 0..3 {
        ordinary.step().unwrap();
        let report = flow
            .step_with_submerged_solid_loads(&objects, controls(), &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(flow.grid.state_bits(), ordinary.grid.state_bits());
        assert_eq!(flow.mass_snapshot(), ordinary.mass_snapshot());
        assert_eq!(flow.conversions(), ordinary.conversions());
        if step == 0 {
            // Cube-side n: opposite face links contribute n² each, diagonals
            // 2n²-n each. The tensor sum gives 2*(5n²-2n)/3 times rho*u.
            let n = TILE as f64;
            let coefficient = 2.0 * (5.0 * n * n - 2.0 * n) / 3.0;
            for a in 0..3 {
                close(report.bodies[0].wall_impulse[a], coefficient * velocity[a]);
                close(report.bodies[0].wall_angular_impulse[a], 0.0);
            }
        }
    }
}

#[test]
fn multiple_objects_share_one_step_and_joined_tiles_form_one_closed_body() {
    let solid = [(1, 1, 1), (3, 1, 1)];
    let objects = [
        SubmergedSolid3 {
            solid_tiles: &solid[..1],
            moment_origin: [6.0; 3],
        },
        SubmergedSolid3 {
            solid_tiles: &solid[1..],
            moment_origin: [14.0, 6.0, 6.0],
        },
    ];
    let input = tank(&solid, [0.0, 0.0, -1e-5], [0.0; 3], 1.0, 16);
    let gate = CancelGate::new_clock_free();
    let (first, report) = input
        .submerged_solid_load_candidate(&objects, controls(), &gate)
        .unwrap();
    let (second, reverse) = input
        .submerged_solid_load_candidate(&[objects[1], objects[0]], controls(), &gate)
        .unwrap();
    assert_eq!(first.grid.state_bits(), second.grid.state_bits());
    assert_eq!(report.bodies[0], reverse.bodies[1]);
    assert_eq!(report.bodies[1], reverse.bodies[0]);
    for body in report.bodies {
        close(body.wall_impulse[2], 64.0e-5);
        for moment in body.wall_angular_impulse {
            close(moment, 0.0);
        }
    }

    let joined = [(1, 1, 1), (2, 1, 1)];
    let mut flow = tank(&joined, [0.0, 0.0, -1e-5], [0.0; 3], 1.0, 16);
    let report = flow
        .step_with_submerged_solid_loads(
            &[SubmergedSolid3 {
                solid_tiles: &joined,
                moment_origin: [8.0, 6.0, 6.0],
            }],
            controls(),
            &gate,
        )
        .unwrap();
    close(report.bodies[0].wall_impulse[2], 128.0e-5);
    for moment in report.bodies[0].wall_angular_impulse {
        close(moment, 0.0);
    }
}

#[test]
fn origin_shift_changes_only_the_moment_and_candidates_do_not_publish() {
    let solid = [(1, 1, 1)];
    let input = tank(&solid, [1e-5, -2e-5, 0.0], [0.0; 3], 1.0, 16);
    let bits = input.grid.state_bits();
    let mass = input.mass_snapshot();
    let origin = [6.0; 3];
    let shift = [2.0, -3.0, 1.0];
    let make = |origin| {
        [SubmergedSolid3 {
            solid_tiles: &solid,
            moment_origin: origin,
        }]
    };
    let gate = CancelGate::new_clock_free();
    let (a, report) = input
        .submerged_solid_load_candidate(&make(origin), controls(), &gate)
        .unwrap();
    let (b, shifted) = input
        .submerged_solid_load_candidate(
            &make(core::array::from_fn(|d| origin[d] + shift[d])),
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(a.grid.state_bits(), b.grid.state_bits());
    let correction = cross(shift, report.bodies[0].wall_impulse);
    for d in 0..3 {
        close(
            shifted.bodies[0].wall_impulse[d],
            report.bodies[0].wall_impulse[d],
        );
        close(
            shifted.bodies[0].wall_angular_impulse[d],
            report.bodies[0].wall_angular_impulse[d] - correction[d],
        );
    }
    drop((a, b));
    assert_eq!(input.grid.state_bits(), bits);
    assert_eq!(input.mass_snapshot(), mass);
}

#[test]
fn wet_surface_budgets_cancellation_and_late_native_failure_are_atomic() {
    let solid = [(1, 1, 1)];
    let objects = [SubmergedSolid3 {
        solid_tiles: &solid,
        moment_origin: [6.0; 3],
    }];
    let gate = CancelGate::new_clock_free();
    let mut flow = tank(&solid, [0.0; 3], [0.0; 3], 1.0, 16);
    let bits = flow.grid.state_bits();
    let mass = flow.mass_snapshot();
    let mut small = controls();
    small.maximum_boundary_links = 431;
    assert!(
        flow.step_with_submerged_solid_loads(&objects, small, &gate)
            .is_err()
    );
    let duplicate = [(1, 1, 1), (1, 1, 1)];
    assert!(
        flow.step_with_submerged_solid_loads(
            &[SubmergedSolid3 {
                solid_tiles: &duplicate,
                moment_origin: [0.0; 3]
            }],
            controls(),
            &gate,
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        flow.step_with_submerged_solid_loads(&objects, controls(), &cancelled,),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert_eq!(flow.grid.state_bits(), bits);
    assert_eq!(flow.mass_snapshot(), mass);

    let mut dry = tank(&solid, [0.0; 3], [0.0; 3], 1.0, 7);
    let dry_bits = dry.grid.state_bits();
    assert!(
        dry.step_with_submerged_solid_loads(&objects, controls(), &gate)
            .is_err()
    );
    assert_eq!(dry.grid.state_bits(), dry_bits);

    // The incoming physical velocity is admissible, but one private native
    // step accelerates interior liquid beyond the existing Mach envelope.
    flow.rebind_body_acceleration([1.5 * MACH_LIMIT * CS2_3.sqrt(), 0.0, 0.0], &gate)
        .unwrap();
    let bits = flow.grid.state_bits();
    let forces = flow.local_force.clone();
    let mass = flow.mass_snapshot();
    assert!(
        flow.step_with_submerged_solid_loads(&objects, controls(), &gate)
            .is_err()
    );
    assert_eq!(flow.grid.state_bits(), bits);
    assert_eq!(flow.local_force, forces);
    assert_eq!(flow.mass_snapshot(), mass);
    flow.rebind_body_acceleration([0.0; 3], &gate).unwrap();
    flow.step_with_submerged_solid_loads(&objects, controls(), &gate)
        .unwrap();
}
