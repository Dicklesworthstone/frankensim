//! G0/G3: finite-cell angular transport is distinct from force and wall torque.
use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, E3, FaceBoundary3, PartialSaturationBody3,
    PartialSaturationBodyExchange3, Q3,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a:e} != {b:e}");
}

fn moments(grid: &BoundaryGrid3, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let [nx, ny, nz] = grid.dimensions();
    let (mut p, mut h) = ([0.0; 3], [0.0; 3]);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let f = grid.populations(x, y, z);
                let mut local = [0.0; 3];
                for q in 0..Q3 {
                    for (axis, e) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[axis] += f64::from(e) * f[q];
                    }
                }
                let r = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                let angular = cross(r, local);
                for i in 0..3 {
                    p[i] += local[i];
                    h[i] += angular[i];
                }
            }
        }
    }
    (p, h)
}

fn packet(case: usize, origin: [f64; 3]) -> PartialSaturationBodyExchange3 {
    let [nx, ny, nz] = [4, 8, 12];
    let mut faces = [FaceBoundary3::Periodic; 6];
    if case == 1 {
        faces[2] = FaceBoundary3::stationary_wall();
        faces[3] = FaceBoundary3::stationary_wall();
    }
    // A very long relaxation time leaves the manufactured packet unchanged to
    // floating-point precision. No collision-kernel helper supplies the oracle.
    let mut grid = BoundaryGrid3::new(nx, ny, nz, 1e30, [0.0; 3], BoundarySpec3::new(faces));
    if case == 2 {
        grid.voxelize_sdf(|p| if p == [0.5, 0.5, 5.5] { -1.0 } else { 1.0 });
    }
    let q = E3.iter().position(|&e| e == (1, 1, 0)).unwrap();
    let delta = 1.0 / 1024.0;
    let mut populations = grid.populations(nx - 1, ny - 1, 5);
    populations[q] += delta;
    populations[0] -= delta;
    grid.set_populations(nx - 1, ny - 1, 5, &populations);
    let before = moments(&grid, origin);
    let mut body = PartialSaturationBody3 {
        center: [2.0, 4.0, 6.0],
        mass: 1.0,
        inertia_world: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    };
    let initial_body = body;
    let report = grid
        .step_with_partial_saturation_body(
            &vec![0.0; nx * ny * nz],
            &mut body,
            [0.0; 6],
            origin,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(body, initial_body);
    let after = moments(&grid, origin);
    for i in 0..3 {
        close(after.0[i] - before.0[i] + report.walls.wall_impulse[i], 0.0);
        close(
            after.1[i] - before.1[i] + report.walls.wall_angular_impulse[i],
            report.periodic.angular_momentum_change[i],
        );
        // The rectangular domain prevents the two jumps from cancelling:
        // (-4,-8,0) cross delta*(1,1,0) = (0,0,4*delta).
        let expected = if case == 0 && i == 2 {
            4.0 * delta
        } else {
            0.0
        };
        close(report.periodic.angular_momentum_change[i], expected);
    }
    assert!(report.periodic.face_crossings > 0);
    if case == 0 {
        assert_eq!(
            report.periodic.face_crossings,
            10 * (nx * ny + nx * nz + ny * nz)
        );
        assert_eq!(report.walls.measured_links, 0);
    } else {
        assert!(report.walls.measured_links > 0);
    }
    report
}

#[test]
fn diagonal_double_wrap_uses_both_domain_lengths_and_is_origin_independent() {
    let a = packet(0, [0.0; 3]);
    let b = packet(0, [-17.0, 9.0, 23.0]);
    assert_eq!(a.periodic, b.periodic);
}

#[test]
fn wall_corner_and_wrapped_solid_links_are_wall_owned_without_periodic_double_counting() {
    packet(1, [1.7, -3.2, 0.9]);
    packet(2, [1.7, -3.2, 0.9]);
}
