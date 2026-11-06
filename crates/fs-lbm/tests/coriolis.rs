//! Local f-plane sources through the original D2Q9/D3Q19 flow operators.
use fs_exec::CancelGate;
use fs_lbm::core2::frame::{Coriolis2, FrameForcingError2};
use fs_lbm::d3q19::{
    BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, CollisionModel3, Coriolis3,
    FaceBoundary3, Q3, RigidFrame3,
};
use fs_lbm::{Grid, equilibrium};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
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

#[test]
fn periodic_plane_inertial_oscillations_have_no_position_dependent_force() {
    for omega in [-0.07_f64, 0.0, 0.07] {
        let source = Coriolis2 { angular_velocity: omega };
        let initial = [0.025, -0.017];
        let mut grid = Grid::uniform(5, 3, 0.9);
        grid.f.fill(equilibrium(1.3, initial[0], initial[1]));
        let mut scratch = Vec::new();
        for _ in 0..200 {
            source.prepare_step(&mut grid, [0.0; 2]).unwrap();
            grid.step(&mut scratch);
        }
        let angle = -400.0 * omega.atan();
        let expected = [
            angle.cos() * initial[0] - angle.sin() * initial[1],
            angle.sin() * initial[0] + angle.cos() * initial[1],
        ];
        for cell in 0..15 {
            let m = grid.moments(cell);
            near(m.rho, 1.3, 5e-13);
            near(m.u[0], expected[0], 4e-14);
            near(m.u[1], expected[1], 4e-14);
            near(m.u[0] * m.u[0] + m.u[1] * m.u[1],
                initial[0] * initial[0] + initial[1] * initial[1], 5e-15);
        }
    }
}

#[test]
fn arbitrary_axis_periodic_box_preserves_inertial_energy_and_axial_velocity() {
    let gate = CancelGate::new();
    let initial = [0.023, -0.017, 0.011];
    let omega = [0.04, -0.03, 0.02];
    let source = Coriolis3 { angular_velocity: omega };
    let magnitude = dot(omega, omega).sqrt();
    let axis = omega.map(|v| v / magnitude);
    let parallel = dot(axis, initial);
    let turn = cross(axis, initial);
    for model in [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt { tau: 0.9, magic: 3.0 / 16.0 },
        CollisionModel3::CentralMoment { second_order_rate: 1.1, higher_order_rate: 1.4 },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.1, third_order_rate: 1.3, fourth_order_rate: 1.5,
        },
    ] {
        let mut grid = BoundaryGrid3::with_collision_model(
            4, 4, 4, model, [0.0; 3], BoundarySpec3::periodic(),
        );
        grid.initialize_uniform_equilibrium(1.3, initial, &gate).unwrap();
        for _ in 0..100 {
            let before = grid.velocity(1, 2, 3);
            source.prepare_step(&mut grid, [0.0; 3], &gate).unwrap();
            let force = grid.force_density(1, 2, 3);
            for d in 0..3 {
                near(grid.velocity(1, 2, 3)[d], before[d], 3e-15);
                near(force[d], grid.force_density(3, 0, 1)[d], 3e-15);
            }
            grid.step();
            let after = grid.velocity(1, 2, 3);
            let midpoint = core::array::from_fn(|d| 0.5 * (before[d] + after[d]));
            near(dot(force, midpoint), 0.0, 3e-15);
        }
        // Rodrigues is independent of the implemented analytic midpoint inverse.
        let angle = 200.0 * magnitude.atan();
        let after = grid.velocity(1, 2, 3);
        for d in 0..3 {
            let expected = axis[d] * parallel
                + angle.cos() * (initial[d] - axis[d] * parallel) - angle.sin() * turn[d];
            near(after[d], expected, 5e-13);
        }
        near(dot(after, after), dot(initial, initial), 5e-14);
        near(dot(axis, after), parallel, 5e-13);
    }
}

#[test]
fn homogeneous_geostrophic_balance_and_zero_spin_use_effective_acceleration() {
    let gate = CancelGate::new();
    let velocity = [0.024, -0.013, 0.009];
    let omega = [0.03, 0.02, -0.01];
    let acceleration = cross(omega, velocity).map(|a| 2.0 * a);
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    grid.initialize_uniform_equilibrium(1.7, velocity, &gate).unwrap();
    for _ in 0..40 {
        Coriolis3 { angular_velocity: omega }
            .prepare_step(&mut grid, acceleration, &gate).unwrap();
        grid.step();
    }
    for d in 0..3 {
        near(grid.velocity(2, 1, 3)[d], velocity[d], 5e-14);
    }
    let before = grid.velocity(2, 1, 3);
    let body = [1e-4, -2e-4, 3e-4];
    Coriolis3::default().prepare_step(&mut grid, body, &gate).unwrap();
    grid.step();
    for d in 0..3 {
        near(grid.velocity(2, 1, 3)[d] - before[d], body[d], 3e-15);
    }

    let mut plane = Grid::uniform(4, 3, 0.9);
    plane.f.fill(equilibrium(1.4, velocity[0], velocity[1]));
    let spin = 0.06;
    let acceleration = [-2.0 * spin * velocity[1], 2.0 * spin * velocity[0]];
    for _ in 0..40 {
        Coriolis2 { angular_velocity: spin }.prepare_step(&mut plane, acceleration).unwrap();
        plane.step(&mut Vec::new());
    }
    near(plane.moments(7).u[0], velocity[0], 5e-14);
    near(plane.moments(7).u[1], velocity[1], 5e-14);
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<([f64; Q3], [f64; 3])> {
    let [nx, ny, nz] = grid.dimensions();
    let mut state = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                state.push((grid.populations(x, y, z), grid.force_density(x, y, z)));
            }
        }
    }
    state
}

#[test]
fn local_source_keeps_full_frame_boundary_and_atomic_refusal_contracts() {
    let gate = CancelGate::new();
    let source = Coriolis3 { angular_velocity: [0.03, -0.02, 0.01] };
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    let before = snapshot(&grid);
    assert!(matches!(
        RigidFrame3 { angular_velocity: source.angular_velocity, ..Default::default() }
            .prepare_step(&mut grid, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Invalid(_))
    ));
    assert_eq!(snapshot(&grid), before);
    source.prepare_step(&mut grid, [0.0; 3], &gate).unwrap();
    // Rest has no position-dependent centrifugal load, up to moment roundoff.
    for force in [grid.force_density(0, 0, 0), grid.force_density(3, 2, 1)] {
        for component in force {
            near(component, 0.0, 1e-17);
        }
    }
    let before = snapshot(&grid);

    for omega in [[f64::NAN, 0.0, 0.0], [f64::MAX, 0.0, 0.0]] {
        assert!(Coriolis3 { angular_velocity: omega }
            .prepare_step(&mut grid, [0.0; 3], &gate).is_err());
        assert_eq!(snapshot(&grid), before);
    }
    grid.set_populations(3, 3, 3, &[-1.0; Q3]);
    let invalid = snapshot(&grid);
    assert_eq!(source.prepare_step(&mut grid, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Unphysical { cell: [3, 3, 3] }));
    assert_eq!(snapshot(&grid), invalid);
    gate.request();
    assert_eq!(source.prepare_step(&mut grid, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Cancelled));
    assert_eq!(snapshot(&grid), invalid);

    let mut open = BoundaryGrid3::new(
        4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::pressure_duct_z(1.0, 1.0),
    );
    assert!(source.prepare_step(&mut open, [0.0; 3], &CancelGate::new()).is_err());
    let mut closed = BoundaryGrid3::new(
        4, 4, 4, 0.9, [0.0; 3],
        BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]),
    );
    source.prepare_step(&mut closed, [0.0; 3], &CancelGate::new()).unwrap();
}

#[test]
fn plane_source_refuses_a_late_bad_cell_without_publishing() {
    let mut grid = Grid::uniform(3, 2, 0.9);
    grid.f[5] = [-1.0; 9];
    let before = grid.clone();
    assert_eq!(
        Coriolis2 { angular_velocity: 0.03 }.prepare_step(&mut grid, [0.0; 2]),
        Err(FrameForcingError2::InvalidCell { cell: 5 })
    );
    assert_eq!(grid.f, before.f);
    assert_eq!(grid.fext, before.fext);
    assert_eq!(grid.g, before.g);
}
