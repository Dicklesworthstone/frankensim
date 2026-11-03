//! Prescribed spinning film ports on the original sampled body/fluid clock.
//! These checks do not supply a resolved spinning wall boundary.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, LubricationPlane3,
    MovingPartialSaturationBody3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
const RADIUS: f64 = 0.5;
const MU: f64 = 0.11;
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn grid() -> BoundaryGrid3 {
    BoundaryGrid3::new(8, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic())
}
fn body() -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [4.0, 3.0, 2.505],
        orientation: So3::identity(),
        mass: 2.0,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                0.06, 0.008, -0.004, 0.008, 0.08, 0.003, -0.004, 0.003, 0.095,
            ],
        })
        .unwrap(),
        velocity: [0.004, 0.002, 0.0001],
        angular_momentum_body: Vec3::new(0.000008, -0.000006, 0.000004),
    }
}
fn plane(spin: f64) -> LubricationPlane3 {
    LubricationPlane3::new(
        [3.0, 3.0, 2.0],
        [0.0, 0.0, 1.0],
        [0.003, -0.004, 0.00005],
        0.04,
        0.0001,
    )
    .unwrap()
    .with_tangential(0.01, 0.0001)
    .unwrap()
    .with_axial_spin(spin)
    .unwrap()
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(MU, 0.04, 0.0001, 1).unwrap()
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut values = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                values.push(grid.populations(x, y, z));
            }
        }
    }
    values
}

#[test]
fn sampled_spinning_surface_uses_translated_midpoint_and_native_body_clock() {
    let mut grid = grid();
    let initial = body();
    let mut bodies = [initial];
    let mut samplers = [BodyFrameSampler3::new(&grid, 4).unwrap()];
    let plane = plane(0.01);
    let result = grid
        .step_with_lubricated_spheres_and_planes(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [[0.0; 6]],
            model(),
            &[plane],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let motion = result.motion;
    let coupling = motion.coupling;
    let film = coupling.lubrication;
    assert!(coupling.bodies[0].hydrodynamic.coupled_cells > 0);
    assert_eq!(film.active_resistance_modes, 5);
    let midpoint: [f64; 3] =
        core::array::from_fn(|i| initial.center[i] + 0.5 * initial.velocity[i]);
    let axis: [f64; 3] = core::array::from_fn(|i| plane.point()[i] + 0.5 * plane.velocity()[i]);
    let local = [
        plane.velocity()[0] - 0.01 * (midpoint[1] - axis[1]),
        plane.velocity()[1] + 0.01 * (midpoint[0] - axis[0]),
        plane.velocity()[2],
    ];
    let gap = midpoint[2] - axis[2] - RADIUS;
    let k = 4.0 * core::f64::consts::PI * MU * RADIUS / 5.0 * (0.01 / gap).ln();
    let zeta = 6.0 * core::f64::consts::PI * MU * RADIUS * RADIUS * (1.0 / gap - 1.0 / 0.04);
    let v = coupling.bodies[0].velocity_midpoint;
    let w = coupling.bodies[0].angular_velocity_midpoint;
    let u: [f64; 3] = core::array::from_fn(|i| v[i] - local[i]);
    // The independent physical tangent/parallel-spin block, not its production
    // positive-mode factorization. Both tangent directions and squeeze share v.
    let force = [
        -4.0 * k * u[0] + k * RADIUS * w[1],
        -4.0 * k * u[1] - k * RADIUS * w[0],
        -zeta * u[2],
    ];
    let torque = [
        -k * RADIUS * u[1] - 4.0 * k * RADIUS * RADIUS * w[0],
        k * RADIUS * u[0] - 4.0 * k * RADIUS * RADIUS * w[1],
        0.0,
    ];
    let support = result.planes[0];
    for i in 0..3 {
        close(film.body_impulses[0][i], force[i], 2e-14);
        close(film.body_spin_impulses[0][i], torque[i], 2e-14);
        close(film.body_impulses[0][i] + support.impulse[i], 0.0, 2e-16);
        close(
            film.body_angular_impulses[0][i] + support.angular_impulse[i],
            0.0,
            2e-15,
        );
        close(motion.bodies[0].midpoint_center[i], midpoint[i], 1e-15);
        close(
            bodies[0].center[i],
            midpoint[i] + 0.5 * bodies[0].velocity[i],
            1e-15,
        );
        close(motion.bodies[0].momentum_residual[i], 0.0, 3e-14);
        close(motion.bodies[0].momentum_residual[i + 3], 0.0, 3e-14);
    }
    close(
        support.work,
        dot(support.impulse, local) + 0.01 * support.couple_impulse[2],
        2e-17,
    );
    let reference_moment = cross(
        core::array::from_fn(|i| axis[i] - ORIGIN[i]),
        support.impulse,
    );
    close(
        support.work,
        dot(support.impulse, plane.velocity())
            + 0.01 * (support.angular_impulse[2] - reference_moment[2]),
        2e-17,
    );
    close(film.work_residual, 0.0, 3e-17);
    close(motion.bodies[0].work_residual, 0.0, 3e-16);
}

#[test]
fn spinning_surface_local_speed_refusal_is_atomic_and_retry_matches_fresh() {
    let mut grid = grid();
    let original_fluid = snapshot(&grid);
    let initial = body();
    let mut bodies = [initial];
    let mut samplers = [BodyFrameSampler3::new(&grid, 4).unwrap()];
    assert!(
        grid.step_with_lubricated_spheres_and_planes(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [[0.0; 6]],
            model(),
            &[plane(1.0)],
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, [initial]);
    assert_eq!(snapshot(&grid), original_fluid);
    let accepted = grid
        .step_with_lubricated_spheres_and_planes(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [[0.0; 6]],
            model(),
            &[plane(0.01)],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let mut fresh = BoundaryGrid3::new(8, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic());
    let mut expected = [initial];
    let mut fresh_sampler = [BodyFrameSampler3::new(&fresh, 4).unwrap()];
    let repeated = fresh
        .step_with_lubricated_spheres_and_planes(
            &mut fresh_sampler,
            [RADIUS],
            &mut expected,
            [[0.0; 6]],
            model(),
            &[plane(0.01)],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(accepted, repeated);
    assert_eq!(bodies, expected);
    assert_eq!(snapshot(&grid), snapshot(&fresh));
}
