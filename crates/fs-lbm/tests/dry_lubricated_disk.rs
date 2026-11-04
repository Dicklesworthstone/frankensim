//! Physical cylinder/plane films with independent momentum and resistance checks.
use fs_exec::CancelGate;
use fs_lbm::core2::moving_disk::{
    DryDisk2, DryDiskCouplingConfig2, DryDiskError2, LubricationPlane2, disk_cell_moments,
};
use fs_lbm::core2::{Cell, Grid, PartialSaturationBody2};
const DIRECTIONS: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; delta={:.3e}, tolerance={tolerance:e}",
        a - b
    );
}
fn vector(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for k in 0..3 {
        close(a[k], b[k], tolerance);
    }
}
fn fixture(density: f64) -> (Grid, DryDisk2) {
    let mut grid = Grid::uniform(26, 24, 0.82);
    grid.periodic_x = true;
    grid.periodic_y = true;
    let radius = 1.6;
    let mass = std::f64::consts::PI * radius * radius * density;
    let body = PartialSaturationBody2 {
        mass,
        inertia: 0.5 * mass * radius * radius,
        center: [12.2, 11.7],
        velocity: [-0.008, -0.004],
        angular_velocity: 0.001,
    };
    let disk = DryDisk2::new(&grid, radius, body, 0.3).unwrap();
    (grid, disk)
}
fn plane(disk: &DryDisk2, normal: [f64; 2], gap: f64, velocity: [f64; 2]) -> LubricationPlane2 {
    let center = disk.body().center;
    LubricationPlane2::new(
        core::array::from_fn(|k| center[k] - (disk.radius() + gap) * normal[k]),
        normal,
        velocity,
        0.2,
        0.005,
    )
    .unwrap()
}
fn physical(grid: &Grid, i: usize) -> (f64, [f64; 2]) {
    let rho: f64 = grid.f[i].iter().sum();
    let mut momentum = [0.; 2];
    for (q, f) in grid.f[i].iter().enumerate() {
        momentum[0] += f * f64::from(DIRECTIONS[q].0);
        momentum[1] += f * f64::from(DIRECTIONS[q].1);
    }
    for k in 0..2 {
        momentum[k] += 0.5 * (rho * grid.g[k] + grid.fext[i][k]);
    }
    (rho, momentum)
}
fn joint(grid: &Grid, disk: &DryDisk2, origin: [f64; 2]) -> [f64; 3] {
    let body = disk.body();
    let mut out = [
        body.mass * body.velocity[0],
        body.mass * body.velocity[1],
        body.inertia * body.angular_velocity,
    ];
    out[2] += cross(
        [body.center[0] - origin[0], body.center[1] - origin[1]],
        [out[0], out[1]],
    );
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            if grid.flags[i] != Cell::Fluid {
                continue;
            }
            let (_, p) = physical(grid, i);
            let m = disk_cell_moments(
                disk.radius(),
                [x as f64 - body.center[0], y as f64 - body.center[1]],
            );
            out[0] += (1. - m.area) * p[0];
            out[1] += (1. - m.area) * p[1];
            out[2] += cross(
                [
                    (1. - m.area) * (x as f64 - origin[0]) - m.first_moment[0],
                    (1. - m.area) * (y as f64 - origin[1]) - m.first_moment[1],
                ],
                p,
            );
        }
    }
    out
}

#[test]
fn empty_and_inactive_planes_preserve_the_exact_dry_disk_step() {
    let (mut grid, mut disk) = fixture(0.4);
    let mut empty_grid = grid.clone();
    let mut empty_disk = disk.clone();
    let mut far_grid = grid.clone();
    let mut far_disk = disk.clone();
    let far = plane(&disk, [1., 0.], 0.8, [0.002, 0.004]);
    let expected = disk
        .step(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7., 26.],
            Default::default(),
            &gate(),
        )
        .unwrap();
    for (g, d, planes) in [
        (&mut empty_grid, &mut empty_disk, vec![]),
        (&mut far_grid, &mut far_disk, vec![far]),
    ] {
        let actual = d
            .step_with_lubrication_planes(
                g,
                [0.0005, -0.0002],
                0.0003,
                [-7., 26.],
                0.001,
                &planes,
                Default::default(),
                &gate(),
            )
            .unwrap();
        assert_eq!(actual.motion, expected);
        assert_eq!(g.f, grid.f);
        assert_eq!(d.body(), disk.body());
        assert_eq!(d.angle(), disk.angle());
        assert_eq!(actual.body_impulse, [0.; 3]);
        assert_eq!(actual.dissipated_energy, 0.);
    }
}

#[test]
fn oblique_planes_share_one_midpoint_and_close_translated_origin_momentum_and_work() {
    let origin = [-7., 26.];
    let viscosity = 0.001;
    let (mut grid, mut disk) = fixture(0.4);
    let initial = disk.body();
    let before = joint(&grid, &disk, origin);
    let planes = [
        plane(&disk, [1., 0.], 0.07, [0.002, 0.004]),
        plane(&disk, [0.6, 0.8], 0.07, [-0.003, 0.001]),
    ];
    let loads = [0.0005, -0.0002];
    let couple = 0.0003;
    let receipt = disk
        .step_with_lubrication_planes(
            &mut grid,
            loads,
            couple,
            origin,
            viscosity,
            &planes,
            Default::default(),
            &gate(),
        )
        .unwrap();
    let midpoint = receipt.motion.native.velocity_midpoint;
    let center = receipt.motion.collision_center;
    let mut total = [0.; 3];
    let mut loss = 0.;
    for (plane, exchange) in planes.into_iter().zip(&receipt.planes) {
        let point = plane.translated(0.5).unwrap().point();
        let normal = plane.normal();
        let gap = dot([center[0] - point[0], center[1] - point[1]], normal) - disk.radius();
        let effective = gap.max(plane.minimum_gap());
        let k = 3.
            * std::f64::consts::PI
            * std::f64::consts::SQRT_2
            * viscosity
            * disk.radius().powf(1.5)
            * (effective.powf(-1.5) - plane.cutoff_gap().powf(-1.5));
        let rate = dot(
            [
                midpoint[0] - plane.velocity()[0],
                midpoint[1] - plane.velocity()[1],
            ],
            normal,
        );
        let j = normal.map(|n| -k * rate * n);
        let torque = cross([center[0] - origin[0], center[1] - origin[1]], j);
        close(exchange.native.resistance_increment, k, 2e-13);
        vector(
            [
                exchange.native.first_solid_impulse[0],
                exchange.native.first_solid_impulse[1],
                exchange.native.first_solid_angular_impulse,
            ],
            [j[0], j[1], torque],
            2e-13,
        );
        for a in 0..2 {
            assert_eq!(
                exchange.native.first_solid_impulse[a],
                -exchange.native.second_solid_impulse[a]
            );
        }
        assert_eq!(
            exchange.native.first_solid_angular_impulse,
            -exchange.native.second_solid_angular_impulse
        );
        close(exchange.body_work, dot(j, midpoint), 2e-16);
        close(
            exchange.plane_work,
            dot(j.map(|v| -v), plane.velocity()),
            2e-16,
        );
        close(
            exchange.body_work + exchange.plane_work,
            -k * rate * rate,
            2e-16,
        );
        total[0] += j[0];
        total[1] += j[1];
        total[2] += torque;
        loss += k * rate * rate;
    }
    vector(receipt.body_impulse, total, 3e-13);
    close(receipt.dissipated_energy, loss, 2e-16);
    close(
        receipt.body_work + receipt.plane_work + receipt.dissipated_energy,
        0.,
        2e-16,
    );
    let after = joint(&grid, &disk, origin);
    vector(
        core::array::from_fn(|a| after[a] - before[a]),
        [
            loads[0] + total[0],
            loads[1] + total[1],
            couple
                + cross([center[0] - origin[0], center[1] - origin[1]], loads)
                + total[2]
                + receipt.motion.periodic_angular_impulse,
        ],
        2e-10,
    );
    close(
        receipt.motion.native.kinetic_energy_change,
        receipt.motion.native.hydrodynamic.solid_work
            + receipt.motion.physical_external_work
            + receipt.motion.interior_correction_work
            + receipt.body_work,
        2e-15,
    );
    assert_eq!(disk.body().mass, initial.mass);
    assert_eq!(disk.body().inertia, initial.inertia);
    assert!(receipt.motion.candidate_evaluations > 4);
}

#[test]
fn advancing_plane_supplies_body_work_with_equal_opposite_reservoir_work() {
    let (grid0, original) = fixture(0.2);
    let mut grid = grid0;
    // Exact transverse reflection symmetry: normal-only films add no spin
    // couple, while an off-center discretized accelerating flow need not be
    // reflection symmetric even with an exact geometric interior moment.
    let initial = PartialSaturationBody2 {
        center: [12.2, 12.],
        velocity: [0.; 2],
        angular_velocity: 0.,
        ..original.body()
    };
    let mut disk = DryDisk2::new(&grid, original.radius(), initial, 0.).unwrap();
    let mut wall = plane(&disk, [1., 0.], 0.07, [0.012, 0.006]);
    for _ in 0..5 {
        let receipt = disk
            .step_with_lubrication_planes(
                &mut grid,
                [0.; 2],
                0.,
                [3., -5.],
                0.003,
                &[wall],
                Default::default(),
                &gate(),
            )
            .unwrap();
        assert!(receipt.body_impulse[0] > 0.);
        assert!(receipt.body_work > 0.);
        assert!(receipt.plane_work < 0.);
        assert!(receipt.dissipated_energy >= 0.);
        close(
            receipt.body_work + receipt.plane_work + receipt.dissipated_energy,
            0.,
            2e-15,
        );
        close(disk.body().angular_velocity, 0., 2e-13);
        wall = wall.translated(1.).unwrap();
        assert!(disk.body().center[0] - wall.point()[0] > disk.radius());
        vector(receipt.motion.body_momentum_residual, [0.; 3], 3e-10);
        vector(receipt.motion.combined_momentum_residual, [0.; 3], 3e-10);
    }
}

#[test]
fn films_leave_density_dependent_guo_source_and_cut_moment_subtraction_separate() {
    let origin = [-5., 29.];
    let (mut grid, mut disk) = fixture(0.7);
    grid.g = [0.00002, -0.00001];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            let rho = 1. + 0.0005 * x as f64 + 0.0002 * y as f64;
            grid.f[i] = fs_lbm::equilibrium(rho, 0.002, -0.001);
            grid.fext[i] = [1e-6 * y as f64, -2e-6 * x as f64];
        }
    }
    let old_grid = grid.clone();
    let before = joint(&grid, &disk, origin);
    let plane = plane(&disk, [0.6, 0.8], 0.08, [0.001, 0.]);
    let loads = [-0.0002, 0.0004];
    let receipt = disk
        .step_with_lubrication_planes(
            &mut grid,
            loads,
            -0.0001,
            origin,
            0.001,
            &[plane],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let center = receipt.motion.collision_center;
    let mut exterior_source = [0.; 3];
    let mut interior_source = [0.; 3];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            let old_rho = physical(&old_grid, i).0;
            let rho = physical(&grid, i).0;
            let f = core::array::from_fn::<_, 2, _>(|a| {
                0.5 * (old_rho + rho) * grid.g[a] + grid.fext[i][a]
            });
            let m = disk_cell_moments(disk.radius(), [x as f64 - center[0], y as f64 - center[1]]);
            let full = [
                f[0],
                f[1],
                cross([x as f64 - origin[0], y as f64 - origin[1]], f),
            ];
            let inside = [
                m.area * f[0],
                m.area * f[1],
                cross(
                    [
                        m.area * (x as f64 - origin[0]) + m.first_moment[0],
                        m.area * (y as f64 - origin[1]) + m.first_moment[1],
                    ],
                    f,
                ),
            ];
            for k in 0..3 {
                interior_source[k] += inside[k];
                exterior_source[k] += full[k] - inside[k];
            }
        }
    }
    vector(
        receipt.motion.interior_forcing_impulse,
        interior_source,
        2e-13,
    );
    vector(
        receipt.motion.represented_fluid_forcing_impulse,
        exterior_source,
        2e-13,
    );
    let after = joint(&grid, &disk, origin);
    let external = [
        loads[0],
        loads[1],
        cross([center[0] - origin[0], center[1] - origin[1]], loads) - 0.0001,
    ];
    vector(
        core::array::from_fn(|k| after[k] - before[k]),
        core::array::from_fn(|k| {
            external[k]
                + exterior_source[k]
                + receipt.body_impulse[k]
                + if k == 2 {
                    receipt.motion.periodic_angular_impulse
                } else {
                    0.
                }
        }),
        3e-10,
    );
}

#[test]
fn refused_film_candidate_preserves_every_accepted_state_and_retry() {
    let (mut grid, mut disk) = fixture(0.4);
    let planes = [
        plane(&disk, [1., 0.], 0.07, [0.002, 0.004]),
        plane(&disk, [0.6, 0.8], 0.07, [-0.003, 0.001]),
    ];
    let original_planes = planes;
    let populations = grid.f.clone();
    let body = disk.body();
    let angle = disk.angle();
    let budget = DryDiskCouplingConfig2 {
        maximum_evaluations: 2,
        ..Default::default()
    };
    assert_eq!(
        disk.step_with_lubrication_planes(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7., 26.],
            0.001,
            &planes,
            budget,
            &gate()
        ),
        Err(DryDiskError2::EvaluationBudget)
    );
    assert_eq!(grid.f, populations);
    assert_eq!(disk.body(), body);
    assert_eq!(disk.angle(), angle);
    assert_eq!(planes, original_planes);
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        disk.step_with_lubrication_planes(
            &mut grid,
            [0.; 2],
            0.,
            [0.; 2],
            0.001,
            &planes,
            Default::default(),
            &cancelled
        ),
        Err(DryDiskError2::Cancelled)
    );
    assert_eq!(grid.f, populations);
    assert_eq!(disk.body(), body);
    let mut clean_grid = grid.clone();
    let mut clean_disk = DryDisk2::new(&grid, disk.radius(), body, angle).unwrap();
    let clean = clean_disk
        .step_with_lubrication_planes(
            &mut clean_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7., 26.],
            0.001,
            &planes,
            Default::default(),
            &gate(),
        )
        .unwrap();
    let retried = disk
        .step_with_lubrication_planes(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7., 26.],
            0.001,
            &planes,
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(retried, clean);
    assert_eq!(grid.f, clean_grid.f);
    assert_eq!(disk.body(), clean_disk.body());
    let (mut grid, mut disk) = fixture(0.4);
    let initial = disk.body();
    let populations = grid.f.clone();
    let crossing = plane(&disk, [1., 0.], 0.005, [0.; 2]);
    assert!(matches!(
        disk.step_with_lubrication_planes(
            &mut grid,
            [0.; 2],
            0.,
            [0.; 2],
            1e-20,
            &[crossing],
            Default::default(),
            &gate()
        ),
        Err(DryDiskError2::Invalid(_))
    ));
    assert_eq!(grid.f, populations);
    assert_eq!(disk.body(), initial);
}

#[test]
fn implicit_stiff_film_prevents_the_crossing_of_the_uncorrected_candidate() {
    let (mut grid, mut disk) = fixture(0.4);
    let wall = plane(&disk, [1., 0.], 0.005, [0.; 2]);
    let original = disk.body();
    let mut free_grid = grid.clone();
    let mut free_disk = disk.clone();
    let _ = free_disk
        .step(
            &mut free_grid,
            [0.; 2],
            0.,
            [0.; 2],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert!(free_disk.body().center[0] - wall.point()[0] < disk.radius());
    let receipt = disk
        .step_with_lubrication_planes(
            &mut grid,
            [0.; 2],
            0.,
            [-4., 20.],
            0.001,
            &[wall],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert!(disk.body().center[0] - wall.point()[0] > disk.radius());
    assert!(receipt.body_impulse[0] > 0.);
    assert!(receipt.planes[0].native.gap_floor_applied);
    assert!(receipt.body_work <= 0.);
    assert_eq!(receipt.plane_work, 0.);
    assert_eq!(disk.body().mass, original.mass);
    assert_eq!(disk.body().inertia, original.inertia);
    vector(receipt.motion.combined_momentum_residual, [0.; 3], 3e-10);
}

fn sliding_plane(disk: &DryDisk2) -> LubricationPlane2 {
    let center = disk.body().center;
    // Normal correction is inactive here; only the independently enabled
    // tangent/spin block can enter the native coupled solve.
    LubricationPlane2::new(
        [center[0], center[1] - disk.radius() - 0.02],
        [0.0, 1.0],
        [0.004, 0.0],
        0.005,
        0.001,
    )
    .unwrap()
    .with_tangential(0.03, 0.002)
    .unwrap()
}

#[test]
fn simultaneous_sliding_and_spin_films_close_native_body_wall_fluid_momentum_and_work() {
    let origin = [-7.0, 26.0];
    let viscosity = 0.04;
    let (mut grid, mut disk) = fixture(0.4);
    let initial = disk.body();
    let before = joint(&grid, &disk, origin);
    let planes = [
        sliding_plane(&disk),
        plane(&disk, [0.6, 0.8], 0.025, [-0.003, 0.001])
            .with_tangential(0.03, 0.002)
            .unwrap(),
    ];
    let external = [0.0005, -0.0002];
    let external_couple = 0.0003;
    let receipt = disk
        .step_with_lubrication_planes(
            &mut grid,
            external,
            external_couple,
            origin,
            viscosity,
            &planes,
            Default::default(),
            &gate(),
        )
        .unwrap();
    let center = receipt.motion.collision_center;
    let velocity = receipt.motion.native.velocity_midpoint;
    let omega = receipt.motion.native.angular_velocity_midpoint;
    assert!(
        (omega - initial.angular_velocity).abs() > 1e-7,
        "oracle must distinguish the solved midpoint spin from the old spin"
    );
    assert_eq!(receipt.planes[0].native.resistance_increment, 0.0);
    assert!(receipt.planes[0].tangential_resistance_increment > 0.0);
    let mut wall_total = [0.0; 3];
    let mut loss = 0.0;
    for (wall, exchange) in planes.into_iter().zip(&receipt.planes) {
        let point = wall.translated(0.5).unwrap().point();
        let normal = wall.normal();
        let tangent = [-normal[1], normal[0]];
        let gap = dot([center[0] - point[0], center[1] - point[1]], normal) - disk.radius();
        let (cutoff, floor) = wall.tangential_scales().unwrap();
        let k = std::f64::consts::PI
            * viscosity
            * (8.0 * disk.radius()).sqrt()
            * (gap.max(floor).powf(-0.5) - cutoff.powf(-0.5));
        let u = dot(
            [
                velocity[0] - wall.velocity()[0],
                velocity[1] - wall.velocity()[1],
            ],
            tangent,
        );
        let tangent_force = tangent.map(|t| -k * u * t);
        let spin = -k * disk.radius().powi(2) * omega;
        let tangent_work = dot(tangent_force, velocity) + spin * omega;
        let expected = [
            exchange.native.first_solid_impulse[0] + tangent_force[0],
            exchange.native.first_solid_impulse[1] + tangent_force[1],
            exchange.native.first_solid_angular_impulse
                + cross(
                    [center[0] - origin[0], center[1] - origin[1]],
                    tangent_force,
                )
                + spin,
        ];
        vector(exchange.body_impulse, expected, 3e-12);
        close(exchange.body_center_couple, spin, 3e-14);
        close(exchange.tangential_resistance_increment, k, 3e-12);
        close(
            exchange.body_work,
            dot(exchange.native.first_solid_impulse, velocity) + tangent_work,
            3e-14,
        );
        close(
            exchange.plane_work,
            dot(
                [exchange.plane_impulse[0], exchange.plane_impulse[1]],
                wall.velocity(),
            ),
            3e-14,
        );
        close(
            exchange.plane_impulse[2],
            cross(
                [
                    exchange.plane_footpoint[0] - origin[0],
                    exchange.plane_footpoint[1] - origin[1],
                ],
                [exchange.plane_impulse[0], exchange.plane_impulse[1]],
            ) + exchange.plane_couple,
            3e-12,
        );
        for k in 0..3 {
            assert_eq!(exchange.plane_impulse[k], -exchange.body_impulse[k]);
            wall_total[k] += exchange.plane_impulse[k];
        }
        let expected_loss =
            -exchange.native.solid_pair_work + k * (u * u + (disk.radius() * omega).powi(2));
        close(exchange.dissipated_energy, expected_loss, 3e-14);
        loss += expected_loss;
    }
    vector(receipt.body_impulse, wall_total.map(|j| -j), 3e-12);
    close(receipt.dissipated_energy, loss, 3e-14);
    close(receipt.body_work + receipt.plane_work + loss, 0.0, 3e-14);
    let after = joint(&grid, &disk, origin);
    vector(
        core::array::from_fn(|k| after[k] - before[k] + wall_total[k]),
        [
            external[0],
            external[1],
            external_couple
                + cross([center[0] - origin[0], center[1] - origin[1]], external)
                + receipt.motion.periodic_angular_impulse,
        ],
        3e-10,
    );
    close(
        receipt.motion.native.kinetic_energy_change,
        receipt.motion.native.hydrodynamic.solid_work
            + receipt.motion.physical_external_work
            + receipt.motion.interior_correction_work
            + receipt.body_work,
        3e-14,
    );
    vector(receipt.motion.body_momentum_residual, [0.0; 3], 3e-10);
    vector(receipt.motion.combined_momentum_residual, [0.0; 3], 3e-10);
    assert_eq!(disk.body().mass, initial.mass);
    assert_eq!(disk.body().inertia, initial.inertia);
}

#[test]
fn inactive_tangent_option_preserves_the_exact_normal_only_tick() {
    let (mut grid, mut disk) = fixture(0.4);
    let mut other_grid = grid.clone();
    let mut other_disk = disk.clone();
    let wall = plane(&disk, [1.0, 0.0], 0.07, [0.002, 0.004]);
    let with_inactive_tangent = wall.with_tangential(0.03, 0.002).unwrap();
    let expected = disk
        .step_with_lubrication_planes(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &[wall],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let actual = other_disk
        .step_with_lubrication_planes(
            &mut other_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &[with_inactive_tangent],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(other_grid.f, grid.f);
    assert_eq!(other_disk.body(), disk.body());
    assert_eq!(other_disk.angle(), disk.angle());
}

#[test]
fn active_spin_film_refusal_preserves_state_and_exact_retry() {
    let (mut grid, mut disk) = fixture(0.4);
    let wall = sliding_plane(&disk);
    let populations = grid.f.clone();
    let body = disk.body();
    let angle = disk.angle();
    for (viscosity, candidate) in [
        (0.04, wall.with_tangential(0.04, 0.002).unwrap()),
        (f64::MAX, wall),
    ] {
        assert!(matches!(
            disk.step_with_lubrication_planes(
                &mut grid,
                [0.0; 2],
                0.0,
                [0.0; 2],
                viscosity,
                &[candidate],
                Default::default(),
                &gate(),
            ),
            Err(DryDiskError2::Invalid(_))
        ));
        assert_eq!(grid.f, populations);
        assert_eq!(disk.body(), body);
        assert_eq!(disk.angle(), angle);
    }
    assert_eq!(
        disk.step_with_lubrication_planes(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.04,
            &[wall],
            DryDiskCouplingConfig2 {
                maximum_evaluations: 2,
                ..Default::default()
            },
            &gate(),
        ),
        Err(DryDiskError2::EvaluationBudget)
    );
    assert_eq!(grid.f, populations);
    assert_eq!(disk.body(), body);
    assert_eq!(disk.angle(), angle);
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        disk.step_with_lubrication_planes(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.04,
            &[wall],
            Default::default(),
            &cancelled,
        ),
        Err(DryDiskError2::Cancelled)
    );
    let mut clean_grid = grid.clone();
    let mut clean_disk = DryDisk2::new(&grid, disk.radius(), body, angle).unwrap();
    let step = |g: &mut Grid, d: &mut DryDisk2| {
        d.step_with_lubrication_planes(
            g,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.04,
            &[wall],
            Default::default(),
            &gate(),
        )
        .unwrap()
    };
    let expected = step(&mut clean_grid, &mut clean_disk);
    let actual = step(&mut grid, &mut disk);
    assert_eq!(actual, expected);
    assert_eq!(grid.f, clean_grid.f);
    assert_eq!(disk.body(), clean_disk.body());
    assert_eq!(disk.angle(), clean_disk.angle());
}

// Finite guided surfaces exercise the same actual dry-body/fluid transaction.
use fs_lbm::core2::moving_disk::FiniteLubricationSlider2;

fn slider_momentum(slider: FiniteLubricationSlider2, origin: [f64; 2]) -> [f64; 3] {
    let p = slider.velocity().map(|v| slider.mass() * v);
    let center = slider.center();
    [
        p[0],
        p[1],
        cross([center[0] - origin[0], center[1] - origin[1]], p),
    ]
}

#[test]
fn finite_sliders_share_native_motion_and_close_actual_momentum_storage_and_guide_work() {
    let origin = [-7.0, 26.0];
    let viscosity = 0.001;
    let (mut grid, mut disk) = fixture(0.4);
    grid.rebind_gravity([1e-6, -2e-6], &mut Vec::new()).unwrap();
    let normals = [[0.0, 1.0], [0.6, 0.8]];
    let mut sliders = core::array::from_fn::<_, 2, _>(|i| {
        let surface = plane(&disk, normals[i], 0.018, [0.0; 2])
            .with_tangential(0.03, 0.001)
            .unwrap();
        FiniteLubricationSlider2::new(
            surface,
            [0.02, 0.7][i],
            [8.1 + i as f64, 6.3],
            [0.003, -0.002][i],
        )
        .unwrap()
    });
    let before_sliders = sliders;
    let mut before = joint(&grid, &disk, origin);
    for slider in sliders {
        let p = slider_momentum(slider, origin);
        for a in 0..3 {
            before[a] += p[a];
        }
    }
    let loads = [0.0002, -0.0003];
    let radius = disk.radius();
    let receipt = disk
        .step_with_lubrication_sliders(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            origin,
            viscosity,
            &mut sliders,
            &loads,
            Default::default(),
            &gate(),
        )
        .unwrap();
    let native = &receipt.motion.motion;
    let mut actual = joint(&grid, &disk, origin);
    let mut expected = core::array::from_fn::<_, 3, _>(|a| {
        native.physical_external_impulse[a]
            + native.represented_fluid_forcing_impulse[a]
            + native.boundary_fluid_impulse[a]
            + if a == 2 {
                native.periodic_angular_impulse
            } else {
                0.0
            }
    });
    for i in 0..2 {
        let before_slider = before_sliders[i];
        let slider = sliders[i];
        let response = receipt.sliders[i];
        let film = receipt.motion.planes[i];
        let normal = normals[i];
        let tangent = [-normal[1], normal[0]];
        let disk_speed = dot(native.native.velocity_midpoint, tangent);
        let initial_speed = before_slider.tangential_velocity();
        let midpoint = 0.5 * (initial_speed + slider.tangential_velocity());
        let gap = dot(
            [
                native.collision_center[0] - slider.plane().point()[0],
                native.collision_center[1] - slider.plane().point()[1],
            ],
            normal,
        ) - radius;
        let k = std::f64::consts::PI
            * viscosity
            * (8.0 * radius).sqrt()
            * (gap.max(0.001).powf(-0.5) - 0.03_f64.powf(-0.5));
        let effective = k / (1.0 + k / (2.0 * slider.mass()));
        let disk_tangent_impulse =
            -effective * (disk_speed - initial_speed - loads[i] / (2.0 * slider.mass()));
        close(
            dot([film.body_impulse[0], film.body_impulse[1]], tangent),
            disk_tangent_impulse,
            3e-16,
        );
        close(response.velocity_midpoint, midpoint, 1e-16);
        close(
            slider.mass() * (slider.tangential_velocity() - initial_speed),
            loads[i] - disk_tangent_impulse,
            3e-16,
        );
        assert_eq!(slider.plane(), before_slider.plane());
        for a in 0..2 {
            close(
                slider.center()[a] - before_slider.center()[a],
                midpoint * tangent[a],
                2e-15,
            );
        }
        let measured_kinetic_change =
            0.5 * slider.mass() * (slider.tangential_velocity().powi(2) - initial_speed.powi(2));
        close(
            measured_kinetic_change,
            loads[i] * midpoint + film.plane_work,
            2e-18,
        );
        close(
            response.kinetic_energy_change,
            measured_kinetic_change,
            2e-18,
        );
        close(response.work_residual, 0.0, 2e-18);
        close(
            dot(
                [response.guide_impulse[0], response.guide_impulse[1]],
                tangent,
            ),
            0.0,
            1e-16,
        );
        let arm = [
            response.collision_center[0] - origin[0],
            response.collision_center[1] - origin[1],
        ];
        let film_tangent = tangent.map(|t| -disk_tangent_impulse * t);
        close(
            response.guide_impulse[2],
            -film.plane_impulse[2] + cross(arm, film_tangent),
            1e-15,
        );
        close(
            film.body_center_couple,
            -k * radius * radius * native.native.angular_velocity_midpoint,
            3e-16,
        );
        vector(response.momentum_residual, [0.0; 3], 2e-15);
        let p = slider_momentum(slider, origin);
        for a in 0..3 {
            actual[a] += p[a];
            expected[a] += response.external_impulse[a] + response.guide_impulse[a];
        }
    }
    for a in 0..3 {
        close(actual[a] - before[a], expected[a], 2e-11);
    }
    vector(receipt.combined_momentum_residual, [0.0; 3], 2e-11);
    close(
        receipt.motion.body_work + receipt.motion.plane_work + receipt.motion.dissipated_energy,
        0.0,
        2e-17,
    );
    assert!(receipt.motion.dissipated_energy >= 0.0);
}

#[test]
fn finite_stiff_slider_recoils_inside_the_actual_fluid_candidate() {
    let (mut grid, disk) = fixture(0.4);
    let mut body = disk.body();
    body.velocity = [0.004, 0.0];
    body.angular_velocity = 0.0;
    let mut disk = DryDisk2::new(&grid, disk.radius(), body, disk.angle()).unwrap();
    let surface = plane(&disk, [0.0, 1.0], 0.008, [0.0; 2])
        .with_tangential(0.03, 0.001)
        .unwrap();
    let mass = 1e-5;
    let initial_speed = 0.003;
    let mut sliders =
        [FiniteLubricationSlider2::new(surface, mass, [8.0, 5.0], initial_speed).unwrap()];
    let receipt = disk
        .step_with_lubrication_sliders(
            &mut grid,
            [0.0; 2],
            0.0,
            [-7.0, 26.0],
            0.002,
            &mut sliders,
            &[0.0],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let film = receipt.motion.planes[0];
    let k = film.tangential_resistance_increment;
    assert!(k / (2.0 * mass) > 1000.0);
    let disk_speed = -receipt.motion.motion.native.velocity_midpoint[0];
    let expected = (mass * initial_speed + 0.5 * k * disk_speed) / (mass + 0.5 * k);
    close(receipt.sliders[0].velocity_midpoint, expected, 2e-17);
    close(
        sliders[0].tangential_velocity(),
        2.0 * expected - initial_speed,
        4e-17,
    );
    assert!(film.plane_impulse[0].abs() <= 2.0 * mass * (disk_speed - initial_speed).abs());
    close(
        receipt.sliders[0].kinetic_energy_change,
        film.plane_work,
        2e-20,
    );
    vector(receipt.combined_momentum_residual, [0.0; 3], 2e-11);
}

#[test]
fn inactive_and_empty_sliders_retain_the_ordinary_native_result_and_free_acceleration() {
    let (mut grid, mut disk) = fixture(0.4);
    let mut expected_grid = grid.clone();
    let mut expected_disk = disk.clone();
    let expected = expected_disk
        .step(
            &mut expected_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let mut empty_grid = grid.clone();
    let mut empty_disk = disk.clone();
    let empty = empty_disk
        .step_with_lubrication_sliders(
            &mut empty_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut [],
            &[],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(empty.motion.motion, expected);
    assert_eq!(empty_grid.f, expected_grid.f);
    let surface = plane(&disk, [0.0, 1.0], 0.8, [0.0; 2])
        .with_tangential(0.03, 0.001)
        .unwrap();
    let mut sliders = [FiniteLubricationSlider2::new(surface, 2.0, [8.0, 5.0], 0.003).unwrap()];
    let receipt = disk
        .step_with_lubrication_sliders(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut sliders,
            &[0.004],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(receipt.motion.motion, expected);
    assert_eq!(grid.f, expected_grid.f);
    assert_eq!(disk.body(), expected_disk.body());
    close(sliders[0].tangential_velocity(), 0.005, 1e-18);
    assert_eq!(receipt.motion.body_impulse, [0.0; 3]);
    vector(receipt.sliders[0].guide_impulse, [0.0; 3], 0.0);
    close(
        receipt.sliders[0].kinetic_energy_change,
        receipt.sliders[0].external_work,
        2e-20,
    );
}

#[test]
fn finite_slider_refusal_and_cancellation_preserve_all_state_and_exact_retry() {
    let (mut grid, mut disk) = fixture(0.4);
    let surface = plane(&disk, [0.0, 1.0], 0.8, [0.0; 2])
        .with_tangential(0.03, 0.001)
        .unwrap();
    let mut sliders = [FiniteLubricationSlider2::new(surface, 1.0, [8.0, 5.0], 0.003).unwrap()];
    let initial_grid = grid.clone();
    let initial_disk = disk.clone();
    let initial_sliders = sliders;
    // Inactive films let the original fluid/disk candidate succeed. Only the
    // late finite-slider endpoint exceeds the speed limit.
    assert!(
        disk.step_with_lubrication_sliders(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut sliders,
            &[1.0],
            Default::default(),
            &gate(),
        )
        .is_err()
    );
    assert_eq!(grid.f, initial_grid.f);
    assert_eq!(disk.body(), initial_disk.body());
    assert_eq!(disk.angle(), initial_disk.angle());
    assert_eq!(sliders, initial_sliders);
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        disk.step_with_lubrication_sliders(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut sliders,
            &[0.0],
            Default::default(),
            &cancelled,
        )
        .unwrap_err(),
        DryDiskError2::Cancelled
    );
    assert_eq!(grid.f, initial_grid.f);
    assert_eq!(sliders, initial_sliders);
    let mut fresh_grid = initial_grid;
    let mut fresh_disk = initial_disk;
    let mut fresh_sliders = initial_sliders;
    let actual = disk
        .step_with_lubrication_sliders(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut sliders,
            &[0.0],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let expected = fresh_disk
        .step_with_lubrication_sliders(
            &mut fresh_grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut fresh_sliders,
            &[0.0],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(grid.f, fresh_grid.f);
    assert_eq!(disk.body(), fresh_disk.body());
    assert_eq!(sliders, fresh_sliders);
}

use fs_lbm::core2::moving_disk::FiniteLubricationPlane2;

fn finite_plane_momentum(plane: FiniteLubricationPlane2, origin: [f64; 2]) -> [f64; 3] {
    let momentum = plane.velocity().map(|v| plane.mass() * v);
    let center = plane.center();
    [
        momentum[0],
        momentum[1],
        cross([center[0] - origin[0], center[1] - origin[1]], momentum),
    ]
}

#[test]
fn translating_planes_close_actual_motion_momentum_and_work_over_changing_gaps() {
    let origin = [-7.0, 26.0];
    let viscosity = 0.001;
    let (mut grid, mut disk) = fixture(0.4);
    grid.rebind_gravity([1e-6, -2e-6], &mut Vec::new()).unwrap();
    let normals = [[0.0, 1.0], [0.6, 0.8]];
    let mut planes = core::array::from_fn::<_, 2, _>(|i| {
        let surface = plane(
            &disk,
            normals[i],
            0.025,
            [[0.003, 0.002], [-0.002, 0.001]][i],
        )
        .with_tangential(0.03, 0.001)
        .unwrap();
        FiniteLubricationPlane2::new(surface, [0.04, 0.7][i], [8.1 + i as f64, 6.3]).unwrap()
    });
    let original = planes;
    let authored = [[0.0002, -0.0001], [-0.0003, 0.0002]];
    let radius = disk.radius();
    for _ in 0..3 {
        let initial = planes;
        let mut before = joint(&grid, &disk, origin);
        for plane in planes {
            let p = finite_plane_momentum(plane, origin);
            for a in 0..3 {
                before[a] += p[a];
            }
        }
        let report = disk
            .step_with_translating_lubrication_planes(
                &mut grid,
                [0.0005, -0.0002],
                0.0003,
                origin,
                viscosity,
                &mut planes,
                &authored,
                Default::default(),
                &gate(),
            )
            .unwrap();
        let native = &report.motion.motion;
        let mut actual = joint(&grid, &disk, origin);
        let mut expected = core::array::from_fn::<_, 3, _>(|a| {
            native.physical_external_impulse[a]
                + native.represented_fluid_forcing_impulse[a]
                + native.boundary_fluid_impulse[a]
                + if a == 2 {
                    native.periodic_angular_impulse
                } else {
                    0.0
                }
        });
        for i in 0..2 {
            let previous = initial[i];
            let after = planes[i];
            let response = report.planes[i];
            let film = report.motion.planes[i];
            let n = normals[i];
            let t = [-n[1], n[0]];
            let midpoint = core::array::from_fn::<_, 2, _>(|a| {
                0.5 * (previous.velocity()[a] + after.velocity()[a])
            });
            let plane_midpoint = core::array::from_fn::<_, 2, _>(|a| {
                previous.plane().point()[a] + 0.5 * previous.velocity()[a]
            });
            let gap = dot(
                [
                    native.collision_center[0] - plane_midpoint[0],
                    native.collision_center[1] - plane_midpoint[1],
                ],
                n,
            ) - radius;
            let normal_k = 3.0
                * std::f64::consts::PI
                * std::f64::consts::SQRT_2
                * viscosity
                * radius.powf(1.5)
                * (gap.max(previous.plane().minimum_gap()).powf(-1.5)
                    - previous.plane().cutoff_gap().powf(-1.5));
            let tangent_k = if gap < 0.03 {
                std::f64::consts::PI
                    * viscosity
                    * (8.0 * radius).sqrt()
                    * (gap.max(0.001).powf(-0.5) - 0.03_f64.powf(-0.5))
            } else {
                0.0
            };
            close(film.native.resistance_increment, normal_k, 2e-12);
            close(film.tangential_resistance_increment, tangent_k, 3e-15);
            close(
                film.native.effective_gap,
                gap.max(previous.plane().minimum_gap()),
                3e-15,
            );
            for (direction, k) in [(n, normal_k), (t, tangent_k)] {
                let free = dot(previous.velocity(), direction)
                    + dot(authored[i], direction) / (2.0 * after.mass());
                let effective = k / (1.0 + k / (2.0 * after.mass()));
                let expected_impulse =
                    -effective * (dot(native.native.velocity_midpoint, direction) - free);
                close(
                    dot([film.body_impulse[0], film.body_impulse[1]], direction),
                    expected_impulse,
                    2e-15,
                );
            }
            close(
                film.native.normal_gap_velocity,
                dot(
                    [
                        native.native.velocity_midpoint[0] - midpoint[0],
                        native.native.velocity_midpoint[1] - midpoint[1],
                    ],
                    n,
                ),
                2e-15,
            );
            for a in 0..2 {
                close(response.velocity_midpoint[a], midpoint[a], 1e-16);
                close(
                    after.mass() * (after.velocity()[a] - previous.velocity()[a]),
                    authored[i][a] + film.plane_impulse[a],
                    2e-15,
                );
                close(after.center()[a] - previous.center()[a], midpoint[a], 2e-15);
                close(
                    after.plane().point()[a] - previous.plane().point()[a],
                    midpoint[a],
                    2e-15,
                );
                close(
                    after.plane().point()[a] - after.center()[a],
                    previous.plane().point()[a] - previous.center()[a],
                    3e-15,
                );
            }
            assert_eq!(after.plane().normal(), previous.plane().normal());
            let actual_energy = 0.5
                * after.mass()
                * (dot(after.velocity(), after.velocity())
                    - dot(previous.velocity(), previous.velocity()));
            close(
                actual_energy,
                dot(authored[i], midpoint) + film.plane_work,
                3e-18,
            );
            close(response.kinetic_energy_change, actual_energy, 3e-18);
            close(response.work_residual, 0.0, 3e-18);
            let arm = [
                response.collision_center[0] - origin[0],
                response.collision_center[1] - origin[1],
            ];
            close(
                response.guide_couple,
                -film.plane_impulse[2] + cross(arm, [film.plane_impulse[0], film.plane_impulse[1]]),
                2e-15,
            );
            close(
                film.body_center_couple,
                -tangent_k * radius * radius * native.native.angular_velocity_midpoint,
                2e-15,
            );
            vector(response.momentum_residual, [0.0; 3], 3e-15);
            let p = finite_plane_momentum(after, origin);
            for a in 0..3 {
                actual[a] += p[a];
                expected[a] +=
                    response.external_impulse[a] + if a == 2 { response.guide_couple } else { 0.0 };
            }
        }
        for a in 0..3 {
            close(actual[a] - before[a], expected[a], 3e-11);
        }
        vector(report.combined_momentum_residual, [0.0; 3], 3e-11);
        close(
            report.motion.body_work + report.motion.plane_work + report.motion.dissipated_energy,
            0.0,
            3e-17,
        );
        assert!(report.motion.dissipated_energy >= 0.0);
    }
    assert!((planes[0].velocity()[1] - original[0].velocity()[1]).abs() > 1e-4);
    assert!((planes[0].plane().point()[1] - original[0].plane().point()[1]).abs() > 1e-4);
}

#[test]
fn finite_normal_recoil_uses_its_solved_gap_when_held_wall_motion_would_cross() {
    let (mut grid, disk) = fixture(0.4);
    let mut body = disk.body();
    body.velocity = [0.0; 2];
    body.angular_velocity = 0.0;
    let mut disk = DryDisk2::new(&grid, disk.radius(), body, disk.angle()).unwrap();
    let surface = plane(&disk, [0.0, 1.0], 0.012, [0.0, 0.02]);
    let mut planes = [FiniteLubricationPlane2::new(surface, 0.02, [8.0, 5.0]).unwrap()];
    let report = disk
        .step_with_translating_lubrication_planes(
            &mut grid,
            [0.0; 2],
            0.0,
            [-7.0, 26.0],
            0.001,
            &mut planes,
            &[[0.0; 2]],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let actual_gap = disk.body().center[1] - planes[0].plane().point()[1] - disk.radius();
    let held_gap =
        disk.body().center[1] - surface.translated(1.0).unwrap().point()[1] - disk.radius();
    assert!(
        held_gap < 0.0,
        "a held-wall endpoint would falsely reject this native response"
    );
    assert!(actual_gap > 0.0);
    assert!(
        planes[0].velocity()[1] < 0.0,
        "finite midpoint recoil reverses this stiff approach"
    );
    let response = report.planes[0];
    let impulse = report.motion.planes[0].plane_impulse;
    close(0.02 * (planes[0].velocity()[1] - 0.02), impulse[1], 2e-16);
    close(
        response.kinetic_energy_change,
        report.motion.planes[0].plane_work,
        2e-18,
    );
    vector(report.combined_momentum_residual, [0.0; 3], 3e-11);
}

#[test]
fn translating_plane_empty_and_inactive_films_keep_native_parity_and_free_vector_acceleration() {
    let (mut grid, mut disk) = fixture(0.4);
    let mut expected_grid = grid.clone();
    let mut expected_disk = disk.clone();
    let expected = expected_disk
        .step(
            &mut expected_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let mut empty_grid = grid.clone();
    let mut empty_disk = disk.clone();
    let empty = empty_disk
        .step_with_translating_lubrication_planes(
            &mut empty_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut [],
            &[],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(empty.motion.motion, expected);
    assert_eq!(empty_grid.f, expected_grid.f);
    let surface = plane(&disk, [0.6, 0.8], 0.8, [0.003, -0.002])
        .with_tangential(0.03, 0.001)
        .unwrap();
    let initial = FiniteLubricationPlane2::new(surface, 2.0, [8.0, 5.0]).unwrap();
    let mut planes = [initial];
    let authored = [0.004, -0.006];
    let report = disk
        .step_with_translating_lubrication_planes(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut planes,
            &[authored],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(report.motion.motion, expected);
    assert_eq!(grid.f, expected_grid.f);
    assert_eq!(disk.body(), expected_disk.body());
    for a in 0..2 {
        close(
            planes[0].velocity()[a],
            initial.velocity()[a] + authored[a] / initial.mass(),
            1e-18,
        );
        let drift = initial.velocity()[a] + 0.5 * authored[a] / initial.mass();
        close(
            planes[0].plane().point()[a] - initial.plane().point()[a],
            drift,
            2e-15,
        );
    }
    assert_eq!(report.motion.body_impulse, [0.0; 3]);
    assert_eq!(report.planes[0].guide_couple, 0.0);
    close(
        report.planes[0].kinetic_energy_change,
        report.planes[0].external_work,
        2e-20,
    );
}

#[test]
fn translating_plane_crossing_speed_and_cancellation_refuse_atomically_then_retry_exactly() {
    let (mut grid, disk) = fixture(0.4);
    let mut body = disk.body();
    body.velocity = [0.0; 2];
    body.angular_velocity = 0.0;
    let mut disk = DryDisk2::new(&grid, disk.radius(), body, disk.angle()).unwrap();
    let point = [body.center[0], body.center[1] - disk.radius() - 0.012];
    // Both film modes are inactive; the late candidate's normal acceleration
    // crosses a gap that the old held-wall endpoint would incorrectly admit.
    let surface = LubricationPlane2::new(point, [0.0, 1.0], [0.0; 2], 0.005, 0.001).unwrap();
    let mut planes = [FiniteLubricationPlane2::new(surface, 1.0, [8.0, 5.0]).unwrap()];
    let initial_grid = grid.clone();
    let initial_disk = disk.clone();
    let initial_planes = planes;
    for authored in [[0.0, 0.04], [0.13, -0.13]] {
        // The second case has individually small components but violates the
        // full vector low-Mach envelope; it moves away and cannot hit the gap guard.
        assert!(
            disk.step_with_translating_lubrication_planes(
                &mut grid,
                [0.0; 2],
                0.0,
                [0.0; 2],
                0.001,
                &mut planes,
                &[authored],
                Default::default(),
                &gate(),
            )
            .is_err()
        );
        assert_eq!(grid.f, initial_grid.f);
        assert_eq!(disk.body(), initial_disk.body());
        assert_eq!(disk.angle(), initial_disk.angle());
        assert_eq!(planes, initial_planes);
    }
    // First-half crossing must refuse even if a later authored impulse would
    // reverse the carrier: the impulse has not occurred during that drift.
    let crossing = LubricationPlane2::new(point, [0.0, 1.0], [0.0, 0.03], 0.005, 0.001).unwrap();
    let mut crossing_planes = [FiniteLubricationPlane2::new(crossing, 1.0, [8.0, 5.0]).unwrap()];
    let crossing_before = crossing_planes;
    assert!(
        disk.step_with_translating_lubrication_planes(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut crossing_planes,
            &[[0.0, -0.1]],
            Default::default(),
            &gate(),
        )
        .is_err()
    );
    assert_eq!(crossing_planes, crossing_before);
    assert_eq!(grid.f, initial_grid.f);
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        disk.step_with_translating_lubrication_planes(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut planes,
            &[[0.0; 2]],
            Default::default(),
            &cancelled,
        )
        .unwrap_err(),
        DryDiskError2::Cancelled
    );
    assert_eq!(grid.f, initial_grid.f);
    assert_eq!(planes, initial_planes);
    let mut fresh_grid = initial_grid;
    let mut fresh_disk = initial_disk;
    let mut fresh_planes = initial_planes;
    let actual = disk
        .step_with_translating_lubrication_planes(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut planes,
            &[[0.0, -0.004]],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let expected = fresh_disk
        .step_with_translating_lubrication_planes(
            &mut fresh_grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut fresh_planes,
            &[[0.0, -0.004]],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(grid.f, fresh_grid.f);
    assert_eq!(disk.body(), fresh_disk.body());
    assert_eq!(planes, fresh_planes);
}

use fs_lbm::core2::moving_disk::{CarrierLubricationPatch2, TranslatingFilmCarrier2};

fn carrier_momentum(carrier: TranslatingFilmCarrier2, origin: [f64; 2]) -> [f64; 3] {
    let p = carrier.velocity().map(|v| carrier.mass() * v);
    [
        p[0],
        p[1],
        cross(
            [
                carrier.center()[0] - origin[0],
                carrier.center()[1] - origin[1],
            ],
            p,
        ),
    ]
}

fn attached_patch(
    disk: &DryDisk2,
    carrier: TranslatingFilmCarrier2,
    owner: usize,
    normal: [f64; 2],
    gap: f64,
    normal_cutoff: f64,
    tangent_cutoff: f64,
) -> CarrierLubricationPatch2 {
    let relative = core::array::from_fn(|a| {
        disk.body().center[a] - (disk.radius() + gap) * normal[a] - carrier.center()[a]
    });
    CarrierLubricationPatch2::new(owner, relative, normal, normal_cutoff, 0.003)
        .unwrap()
        .with_tangential(tangent_cutoff, 0.001)
        .unwrap()
}

fn independent_patch_resistance(
    surface: LubricationPlane2,
    gap: f64,
    radius: f64,
    mu: f64,
) -> (f64, f64) {
    let normal = if gap < surface.cutoff_gap() {
        3.0 * std::f64::consts::PI
            * std::f64::consts::SQRT_2
            * mu
            * radius.powf(1.5)
            * (gap.max(surface.minimum_gap()).powf(-1.5) - surface.cutoff_gap().powf(-1.5))
    } else {
        0.0
    };
    let tangent = match surface.tangential_scales() {
        Some((cutoff, floor)) if gap < cutoff => {
            std::f64::consts::PI
                * mu
                * (8.0 * radius).sqrt()
                * (gap.max(floor).powf(-0.5) - cutoff.powf(-0.5))
        }
        _ => 0.0,
    };
    (normal, tangent)
}

#[test]
fn shared_oblique_patches_use_one_coupled_mass_and_one_kinetic_store() {
    let origin = [-7.0, 26.0];
    let mu = 0.001;
    let (mut grid, mut disk) = fixture(0.4);
    grid.rebind_gravity([1e-6, -2e-6], &mut Vec::new()).unwrap();
    let mut carriers = [
        TranslatingFilmCarrier2::new(0.06, [8.1, 6.3], [-0.002, 0.001]).unwrap(),
        TranslatingFilmCarrier2::new(0.5, [15.1, 6.3], [0.001, -0.001]).unwrap(),
    ];
    let patches = [
        attached_patch(&disk, carriers[0], 0, [0.0, 1.0], 0.021, 0.18, 0.029),
        attached_patch(&disk, carriers[0], 0, [0.6, 0.8], 0.027, 0.24, 0.031),
        attached_patch(&disk, carriers[1], 1, [-0.8, 0.6], 0.025, 0.14, 0.028),
    ];
    let authored = [[0.0004, -0.0003], [-0.0002, 0.0001]];
    let radius = disk.radius();
    let mut differs_from_independent = 0.0_f64;
    for _ in 0..3 {
        let initial = carriers;
        let mut before = joint(&grid, &disk, origin);
        for carrier in carriers {
            let p = carrier_momentum(carrier, origin);
            for a in 0..3 {
                before[a] += p[a];
            }
        }
        let report = disk
            .step_with_lubrication_carriers(
                &mut grid,
                [0.0005, -0.0002],
                0.0003,
                origin,
                mu,
                &mut carriers,
                &authored,
                &patches,
                Default::default(),
                &gate(),
            )
            .unwrap();
        assert_eq!(report.carriers.len(), 2);
        assert_eq!(report.motion.planes.len(), 3);
        let native = &report.motion.motion;
        let vd = native.native.velocity_midpoint;
        let mut resistance = [[0.0; 4]; 2];
        let mut coefficients = [(0.0, 0.0); 3];
        for (j, patch) in patches.iter().enumerate() {
            let owner = patch.carrier();
            let surface = patch.geometry(&initial).unwrap();
            let n = surface.normal();
            let t = [-n[1], n[0]];
            let center = core::array::from_fn::<_, 2, _>(|a| {
                initial[owner].center()[a] + 0.5 * initial[owner].velocity()[a]
            });
            let point = core::array::from_fn::<_, 2, _>(|a| center[a] + patch.point_relative()[a]);
            let gap = dot(
                [
                    native.collision_center[0] - point[0],
                    native.collision_center[1] - point[1],
                ],
                n,
            ) - radius;
            coefficients[j] = independent_patch_resistance(surface, gap, radius, mu);
            let (kn, kt) = coefficients[j];
            close(
                report.motion.planes[j].native.resistance_increment,
                kn,
                3e-12,
            );
            close(
                report.motion.planes[j].tangential_resistance_increment,
                kt,
                4e-15,
            );
            for a in 0..2 {
                for b in 0..2 {
                    resistance[owner][2 * a + b] += kn * n[a] * n[b] + kt * t[a] * t[b];
                }
            }
        }
        assert!(
            resistance[0][1].abs() > 1.0,
            "oracle must exercise off-diagonal coupling"
        );
        assert!((coefficients[0].0 - coefficients[1].0).abs() > 1.0);
        let mut expected_midpoint = [[0.0; 2]; 2];
        for owner in 0..2 {
            let r = resistance[owner];
            let mass = initial[owner].mass();
            let a = [2.0 * mass + r[0], r[1], r[2], 2.0 * mass + r[3]];
            let rhs = [
                2.0 * mass * initial[owner].velocity()[0]
                    + authored[owner][0]
                    + r[0] * vd[0]
                    + r[1] * vd[1],
                2.0 * mass * initial[owner].velocity()[1]
                    + authored[owner][1]
                    + r[2] * vd[0]
                    + r[3] * vd[1],
            ];
            // Independent two-by-two inverse, without production scaling,
            // Cholesky, response matrices or per-face elimination.
            let determinant = a[0] * a[3] - a[1] * a[2];
            expected_midpoint[owner] = [
                (a[3] * rhs[0] - a[1] * rhs[1]) / determinant,
                (a[0] * rhs[1] - a[2] * rhs[0]) / determinant,
            ];
            for axis in 0..2 {
                close(
                    report.carriers[owner].velocity_midpoint[axis],
                    expected_midpoint[owner][axis],
                    4e-15,
                );
                close(
                    carriers[owner].velocity()[axis],
                    2.0 * expected_midpoint[owner][axis] - initial[owner].velocity()[axis],
                    8e-15,
                );
            }
        }
        let mut sums = [[0.0; 5]; 2];
        for (j, patch) in patches.iter().enumerate() {
            let owner = patch.carrier();
            let surface = patch.geometry(&initial).unwrap();
            let after_surface = patch.geometry(&carriers).unwrap();
            let n = surface.normal();
            let t = [-n[1], n[0]];
            let u = expected_midpoint[owner];
            let relative = [vd[0] - u[0], vd[1] - u[1]];
            let (kn, kt) = coefficients[j];
            let normal_force = -kn * dot(relative, n);
            let tangent_force = -kt * dot(relative, t);
            let expected_force =
                core::array::from_fn::<_, 2, _>(|a| normal_force * n[a] + tangent_force * t[a]);
            let film = report.motion.planes[j];
            for a in 0..2 {
                close(film.body_impulse[a], expected_force[a], 3e-13);
                close(after_surface.point()[a] - surface.point()[a], u[a], 4e-15);
            }
            close(film.plane_work, -dot(expected_force, u), 3e-16);
            close(
                film.body_center_couple,
                -kt * radius * radius * native.native.angular_velocity_midpoint,
                3e-15,
            );
            let mass = initial[owner].mass();
            let free = core::array::from_fn::<_, 2, _>(|a| {
                initial[owner].velocity()[a] + authored[owner][a] / (2.0 * mass)
            });
            let independent_n =
                -kn / (1.0 + kn / (2.0 * mass)) * dot([vd[0] - free[0], vd[1] - free[1]], n);
            let independent_t =
                -kt / (1.0 + kt / (2.0 * mass)) * dot([vd[0] - free[0], vd[1] - free[1]], t);
            for a in 0..2 {
                differs_from_independent = differs_from_independent.max(
                    (film.body_impulse[a] - independent_n * n[a] - independent_t * t[a]).abs(),
                );
            }
            let center = report.carriers[owner].collision_center;
            let arm = [center[0] - origin[0], center[1] - origin[1]];
            for a in 0..3 {
                sums[owner][a] += film.plane_impulse[a];
            }
            sums[owner][3] += film.plane_work;
            sums[owner][4] +=
                -film.plane_impulse[2] + cross(arm, [film.plane_impulse[0], film.plane_impulse[1]]);
        }
        let mut actual = joint(&grid, &disk, origin);
        let mut expected = core::array::from_fn::<_, 3, _>(|a| {
            native.physical_external_impulse[a]
                + native.represented_fluid_forcing_impulse[a]
                + native.boundary_fluid_impulse[a]
                + if a == 2 {
                    native.periodic_angular_impulse
                } else {
                    0.0
                }
        });
        for owner in 0..2 {
            let receipt = report.carriers[owner];
            vector(
                receipt.film_impulse,
                [sums[owner][0], sums[owner][1], sums[owner][2]],
                4e-15,
            );
            close(receipt.film_work, sums[owner][3], 3e-18);
            close(receipt.guide_couple, sums[owner][4], 4e-15);
            let before_v = initial[owner].velocity();
            let after_v = carriers[owner].velocity();
            let kinetic =
                0.5 * initial[owner].mass() * (dot(after_v, after_v) - dot(before_v, before_v));
            close(
                kinetic,
                dot(authored[owner], receipt.velocity_midpoint) + sums[owner][3],
                4e-18,
            );
            close(receipt.kinetic_energy_change, kinetic, 3e-18);
            vector(receipt.momentum_residual, [0.0; 3], 4e-15);
            let p = carrier_momentum(carriers[owner], origin);
            for a in 0..3 {
                actual[a] += p[a];
                expected[a] +=
                    receipt.external_impulse[a] + if a == 2 { receipt.guide_couple } else { 0.0 };
            }
        }
        for a in 0..3 {
            close(actual[a] - before[a], expected[a], 3e-11);
        }
        vector(report.combined_momentum_residual, [0.0; 3], 3e-11);
        close(
            report.motion.body_work + report.motion.plane_work + report.motion.dissipated_energy,
            0.0,
            3e-17,
        );
    }
    assert!(
        differs_from_independent > 1e-5,
        "independent per-patch masses must fail this oracle"
    );
}

#[test]
fn one_shared_patch_matches_the_original_finite_translating_plane_owner() {
    let (mut shared_grid, mut shared_disk) = fixture(0.4);
    let mut original_grid = shared_grid.clone();
    let mut original_disk = shared_disk.clone();
    let mut carriers = [TranslatingFilmCarrier2::new(0.4, [8.1, 6.3], [-0.002, 0.001]).unwrap()];
    let patches = [attached_patch(
        &shared_disk,
        carriers[0],
        0,
        [0.6, 0.8],
        0.025,
        0.2,
        0.03,
    )];
    let mut planes = [FiniteLubricationPlane2::new(
        patches[0].geometry(&carriers).unwrap(),
        carriers[0].mass(),
        carriers[0].center(),
    )
    .unwrap()];
    let authored = [[0.0004, -0.0003]];
    let shared = shared_disk
        .step_with_lubrication_carriers(
            &mut shared_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut carriers,
            &authored,
            &patches,
            Default::default(),
            &gate(),
        )
        .unwrap();
    let original = original_disk
        .step_with_translating_lubrication_planes(
            &mut original_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut planes,
            &authored,
            Default::default(),
            &gate(),
        )
        .unwrap();
    for a in 0..2 {
        close(carriers[0].velocity()[a], planes[0].velocity()[a], 2e-13);
        close(carriers[0].center()[a], planes[0].center()[a], 2e-13);
        close(
            shared_disk.body().velocity[a],
            original_disk.body().velocity[a],
            2e-13,
        );
    }
    close(
        shared_disk.body().angular_velocity,
        original_disk.body().angular_velocity,
        2e-13,
    );
    for (left, right) in shared_grid
        .f
        .iter()
        .flatten()
        .zip(original_grid.f.iter().flatten())
    {
        close(*left, *right, 2e-13);
    }
    vector(
        shared.motion.planes[0].body_impulse,
        original.motion.planes[0].body_impulse,
        2e-13,
    );
    close(
        shared.carriers[0].kinetic_energy_change,
        original.planes[0].kinetic_energy_change,
        2e-16,
    );
}

#[test]
fn inactive_shared_faces_apply_authored_impulse_once_and_preserve_native_bits() {
    let (mut grid, mut disk) = fixture(0.4);
    let mut expected_grid = grid.clone();
    let mut expected_disk = disk.clone();
    let expected = expected_disk
        .step(
            &mut expected_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            Default::default(),
            &gate(),
        )
        .unwrap();
    let mut empty_grid = grid.clone();
    let mut empty_disk = disk.clone();
    let empty = empty_disk
        .step_with_lubrication_carriers(
            &mut empty_grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut [],
            &[],
            &[],
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(empty.motion.motion, expected);
    assert_eq!(empty_grid.f, expected_grid.f);
    let mut carriers = [TranslatingFilmCarrier2::new(2.0, [8.0, 5.0], [0.003, -0.002]).unwrap()];
    let initial = carriers;
    let patches = [
        attached_patch(&disk, carriers[0], 0, [0.0, 1.0], 0.8, 0.2, 0.03),
        attached_patch(&disk, carriers[0], 0, [0.6, 0.8], 0.8, 0.2, 0.03),
    ];
    let authored = [[0.004, -0.006]];
    let report = disk
        .step_with_lubrication_carriers(
            &mut grid,
            [0.0005, -0.0002],
            0.0003,
            [-7.0, 26.0],
            0.001,
            &mut carriers,
            &authored,
            &patches,
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(report.motion.motion, expected);
    assert_eq!(grid.f, expected_grid.f);
    assert_eq!(disk.body(), expected_disk.body());
    assert_eq!(report.carriers.len(), 1);
    for a in 0..2 {
        close(
            carriers[0].velocity()[a],
            initial[0].velocity()[a] + authored[0][a] / 2.0,
            1e-18,
        );
    }
    assert_eq!(report.carriers[0].film_impulse, [0.0; 3]);
    assert_eq!(report.carriers[0].guide_couple, 0.0);
    close(
        report.carriers[0].kinetic_energy_change,
        report.carriers[0].external_work,
        2e-20,
    );
}

#[test]
fn shared_patch_refusals_leave_all_state_and_exact_retry_intact() {
    let (mut grid, disk) = fixture(0.4);
    let mut body = disk.body();
    body.velocity = [0.0; 2];
    body.angular_velocity = 0.0;
    let mut disk = DryDisk2::new(&grid, disk.radius(), body, disk.angle()).unwrap();
    let mut carriers = [TranslatingFilmCarrier2::new(1.0, [8.0, 5.0], [0.0; 2]).unwrap()];
    let patch = |owner, n: [f64; 2], gap: f64, cutoff| {
        CarrierLubricationPatch2::new(
            owner,
            core::array::from_fn(|a| {
                body.center[a] - (disk.radius() + gap) * n[a] - carriers[0].center()[a]
            }),
            n,
            cutoff,
            0.001,
        )
        .unwrap()
    };
    let patches = [
        patch(0, [0.0, 1.0], 0.012, 0.005),
        patch(0, [1.0, 0.0], 0.3, 0.005),
    ];
    let bad_index = [patch(2, [0.0, 1.0], 0.012, 0.005)];
    let active = [
        patch(0, [0.0, 1.0], 0.025, 0.2),
        patch(0, [0.6, 0.8], 0.025, 0.2),
    ];
    let initial_grid = grid.clone();
    let initial_disk = disk.clone();
    let initial_carriers = carriers;
    for authored in [[0.0, 0.04], [0.13, -0.13]] {
        assert!(
            disk.step_with_lubrication_carriers(
                &mut grid,
                [0.0; 2],
                0.0,
                [0.0; 2],
                0.001,
                &mut carriers,
                &[authored],
                &patches,
                Default::default(),
                &gate(),
            )
            .is_err()
        );
        assert_eq!(grid.f, initial_grid.f);
        assert_eq!(disk.body(), initial_disk.body());
        assert_eq!(carriers, initial_carriers);
    }
    assert!(
        disk.step_with_lubrication_carriers(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut carriers,
            &[[0.0; 2]],
            &bad_index,
            Default::default(),
            &gate(),
        )
        .is_err()
    );
    let mut unused = [carriers[0], carriers[0]];
    assert!(
        disk.step_with_lubrication_carriers(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut unused,
            &[[0.0; 2]; 2],
            &patches,
            Default::default(),
            &gate(),
        )
        .is_err()
    );
    let mut controls = DryDiskCouplingConfig2::default();
    controls.maximum_condition_number = 1.01;
    assert_eq!(
        disk.step_with_lubrication_carriers(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut carriers,
            &[[0.0; 2]],
            &active,
            controls,
            &gate(),
        )
        .unwrap_err(),
        DryDiskError2::Invalid("shared carrier response is ill-conditioned")
    );
    // A body-relative face near the disk cannot silently lose its motion
    // when a huge COM and local offset cancel. Refuse before any publication.
    let mut remote = [TranslatingFilmCarrier2::new(1.0, [8.0, 1e16], [0.0, 0.004]).unwrap()];
    let remote_before = remote;
    let remote_patch =
        [CarrierLubricationPatch2::new(0, [4.2, -1e16 + 10.0], [0.0, 1.0], 0.2, 0.001).unwrap()];
    assert_eq!(
        disk.step_with_lubrication_carriers(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut remote,
            &[[0.0; 2]],
            &remote_patch,
            Default::default(),
            &gate(),
        )
        .unwrap_err(),
        DryDiskError2::Invalid("shared patch translation is unrepresentable")
    );
    assert_eq!(remote, remote_before);
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        disk.step_with_lubrication_carriers(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut carriers,
            &[[0.0; 2]],
            &patches,
            Default::default(),
            &cancelled,
        )
        .unwrap_err(),
        DryDiskError2::Cancelled
    );
    assert_eq!(grid.f, initial_grid.f);
    assert_eq!(disk.body(), initial_disk.body());
    assert_eq!(carriers, initial_carriers);
    let mut fresh_grid = initial_grid;
    let mut fresh_disk = initial_disk;
    let mut fresh_carriers = initial_carriers;
    let actual = disk
        .step_with_lubrication_carriers(
            &mut grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut carriers,
            &[[0.0, -0.004]],
            &patches,
            Default::default(),
            &gate(),
        )
        .unwrap();
    let expected = fresh_disk
        .step_with_lubrication_carriers(
            &mut fresh_grid,
            [0.0; 2],
            0.0,
            [0.0; 2],
            0.001,
            &mut fresh_carriers,
            &[[0.0, -0.004]],
            &patches,
            Default::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(grid.f, fresh_grid.f);
    assert_eq!(disk.body(), fresh_disk.body());
    assert_eq!(carriers, fresh_carriers);
}
