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
