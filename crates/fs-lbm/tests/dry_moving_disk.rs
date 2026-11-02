//! Physical dry-disk scenarios against independently integrated populations.
use fs_exec::CancelGate;
use fs_lbm::core2::moving_disk::{
    DryDisk2, DryDiskCouplingConfig2, DryDiskError2, disk_cell_moments,
};
use fs_lbm::core2::{Cell, Grid, PartialSaturationBody2};

const C: [[f64; 2]; 9] = [
    [0., 0.],
    [1., 0.],
    [0., 1.],
    [-1., 0.],
    [0., -1.],
    [1., 1.],
    [-1., 1.],
    [-1., -1.],
    [1., -1.],
];
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, delta {:.3e}, tolerance {tolerance:.3e}",
        actual - expected
    );
}
fn vector_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for k in 0..3 {
        close(actual[k], expected[k], tolerance);
    }
}
fn physical(grid: &Grid, i: usize) -> (f64, [f64; 2]) {
    let rho = grid.f[i].iter().sum::<f64>();
    let mut p = [0.; 2];
    for (value, c) in grid.f[i].iter().zip(C) {
        p[0] += value * c[0];
        p[1] += value * c[1];
    }
    for k in 0..2 {
        p[k] += 0.5 * (rho * grid.g[k] + grid.fext[i][k]);
    }
    (rho, p)
}
fn joint(grid: &Grid, disk: &DryDisk2, origin: [f64; 2]) -> [f64; 3] {
    let body = disk.body();
    let mut result = [
        body.mass * body.velocity[0],
        body.mass * body.velocity[1],
        body.inertia * body.angular_velocity,
    ];
    result[2] += cross(
        [body.center[0] - origin[0], body.center[1] - origin[1]],
        [result[0], result[1]],
    );
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            if grid.flags[i] != Cell::Fluid {
                continue;
            }
            let (_, p) = physical(grid, i);
            let moments = disk_cell_moments(
                disk.radius(),
                [x as f64 - body.center[0], y as f64 - body.center[1]],
            );
            let weight = 1. - moments.area;
            result[0] += weight * p[0];
            result[1] += weight * p[1];
            let exterior_arm = [
                weight * (x as f64 - origin[0]) - moments.first_moment[0],
                weight * (y as f64 - origin[1]) - moments.first_moment[1],
            ];
            result[2] += cross(exterior_arm, p);
        }
    }
    result
}
fn body(radius: f64, density: f64, center: [f64; 2]) -> PartialSaturationBody2 {
    let mass = std::f64::consts::PI * radius * radius * density;
    PartialSaturationBody2 {
        mass,
        inertia: 0.5 * mass * radius * radius,
        center,
        velocity: [0.; 2],
        angular_velocity: 0.,
    }
}
fn periodic_grid(nx: usize, ny: usize) -> Grid {
    let mut grid = Grid::uniform(nx, ny, 0.82);
    grid.periodic_x = true;
    grid.periodic_y = true;
    grid
}
fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn periodic_angular_transport(grid: &Grid) -> f64 {
    let mut result = 0.;
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            if grid.flags[i] != Cell::Fluid {
                continue;
            }
            for (q, c) in C.into_iter().enumerate() {
                if let Some(s) = grid.source(x, y, q) {
                    if grid.flags[s] != Cell::Fluid {
                        continue;
                    }
                    let shift = [
                        x as f64 - (s % grid.nx) as f64 - c[0],
                        y as f64 - (s / grid.nx) as f64 - c[1],
                    ];
                    result += grid.f[i][q] * cross(shift, c);
                }
            }
        }
    }
    result
}

#[test]
fn light_neutral_and_heavy_dry_disks_share_uniform_translation_with_original_inertia() {
    // Reflection symmetry fixes the transverse first moment while analytic
    // fractions change in every axial subcell drift.
    let velocity = [0.018, 0.];
    for density in [0.15, 1., 3.] {
        let mut grid = periodic_grid(24, 22);
        grid.f
            .fill(fs_lbm::equilibrium(1., velocity[0], velocity[1]));
        let mut initial = body(1.7, density, [10.2, 11.]);
        initial.velocity = velocity;
        let mut disk = DryDisk2::new(&grid, 1.7, initial, 0.3).unwrap();
        for _ in 0..16 {
            let exchange = disk
                .step(
                    &mut grid,
                    [0.; 2],
                    0.,
                    [-7., 26.],
                    DryDiskCouplingConfig2::default(),
                    &gate(),
                )
                .unwrap();
            vector_close(exchange.combined_momentum_residual, [0.; 3], 2e-11);
            close(disk.body().velocity[0], velocity[0], 3e-14);
            close(disk.body().velocity[1], velocity[1], 3e-14);
            close(disk.body().angular_velocity, 0., 3e-14);
            assert_eq!(disk.body().mass.to_bits(), initial.mass.to_bits());
            assert_eq!(disk.body().inertia.to_bits(), initial.inertia.to_bits());
        }
        close(
            disk.body().center[0],
            initial.center[0] + 16. * velocity[0],
            2e-13,
        );
        close(
            disk.body().center[1],
            initial.center[1] + 16. * velocity[1],
            2e-13,
        );
    }
}

#[test]
fn reduced_pressure_gravity_rises_or_sediments_with_physical_dry_mass() {
    let radius = 2.;
    let displaced_mass = std::f64::consts::PI * radius * radius;
    for density in [0.5, 2.] {
        let mut grid = periodic_grid(28, 28);
        let initial = body(radius, density, [13.2, 13.7]);
        let mut disk = DryDisk2::new(&grid, radius, initial, 0.).unwrap();
        // Hydrostatic pressure is removed: authored load is weight minus
        // buoyancy at reference fluid density one, with no bulk Guo source.
        let load = [0., -(initial.mass - displaced_mass) * 1e-4];
        for _ in 0..8 {
            let before = joint(&grid, &disk, [0.; 2]);
            let receipt = disk
                .step(
                    &mut grid,
                    load,
                    0.,
                    [0.; 2],
                    DryDiskCouplingConfig2::default(),
                    &gate(),
                )
                .unwrap();
            let after = joint(&grid, &disk, [0.; 2]);
            vector_close(
                core::array::from_fn(|k| after[k] - before[k]),
                [
                    load[0],
                    load[1],
                    cross(receipt.collision_center, load) + receipt.periodic_angular_impulse,
                ],
                2e-10,
            );
            close(
                receipt.physical_external_work,
                load[1] * receipt.native.velocity_midpoint[1],
                1e-17,
            );
            close(
                receipt.native.kinetic_energy_change,
                receipt.physical_external_work
                    + receipt.interior_correction_work
                    + receipt.native.hydrodynamic.solid_work,
                2e-12,
            );
        }
        assert!(disk.body().velocity[1] * (1. - density) > 0.);
        assert_eq!(disk.body().mass, initial.mass);
        assert_eq!(disk.body().inertia, initial.inertia);
    }
}

#[test]
fn counterrotating_disk_exchanges_spin_and_orbital_momentum_about_translated_origin() {
    let mut grid = periodic_grid(26, 25);
    let center = [12.3, 11.7];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let r = [x as f64 - center[0], y as f64 - center[1]];
            let taper = (1. - (r[0] * r[0] + r[1] * r[1]) / 30.).max(0.);
            let velocity = [-0.005 * taper * r[1], 0.005 * taper * r[0]];
            let i = grid.idx(x, y);
            grid.f[i] = fs_lbm::equilibrium(1., velocity[0], velocity[1]);
        }
    }
    let mut initial = body(2., 0.7, center);
    initial.velocity = [0.009, -0.004];
    initial.angular_velocity = -0.004;
    let mut disk = DryDisk2::new(&grid, 2., initial, 0.2).unwrap();
    let origin = [41., -17.];
    for _ in 0..3 {
        let before = joint(&grid, &disk, origin);
        let old = disk.body();
        let old_angle = disk.angle();
        let receipt = disk
            .step(
                &mut grid,
                [0.; 2],
                0.,
                origin,
                DryDiskCouplingConfig2::default(),
                &gate(),
            )
            .unwrap();
        let after = joint(&grid, &disk, origin);
        vector_close(
            core::array::from_fn(|k| after[k] - before[k]),
            [0., 0., receipt.periodic_angular_impulse],
            2e-10,
        );
        assert!(disk.body().angular_velocity > old.angular_velocity);
        assert!(receipt.corrected_hydrodynamic_impulse[2] != 0.);
        close(
            disk.angle() - old_angle,
            0.5 * (old.angular_velocity + disk.body().angular_velocity),
            1e-14,
        );
        vector_close(receipt.body_momentum_residual, [0.; 3], 1e-10);
    }
}

fn forced_fixture() -> (Grid, DryDisk2) {
    let mut grid = periodic_grid(25, 23);
    grid.g = [0.0007, -0.0005];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            let rho = 1. + 0.035 * (0.6 * x as f64).sin() + 0.025 * (0.5 * y as f64).cos();
            grid.fext[i] = [3e-5 * (0.3 * y as f64).sin(), 2e-5 * (0.4 * x as f64).cos()];
            let velocity = [
                0.012 + 0.002 * (0.3 * y as f64).sin(),
                -0.007 + 0.003 * (0.4 * x as f64).cos(),
            ];
            let force = [
                rho * grid.g[0] + grid.fext[i][0],
                rho * grid.g[1] + grid.fext[i][1],
            ];
            grid.f[i] = fs_lbm::equilibrium(
                rho,
                velocity[0] - 0.5 * force[0] / rho,
                velocity[1] - 0.5 * force[1] / rho,
            );
        }
    }
    let mut initial = body(2.2, 0.6, [11.35, 10.7]);
    initial.velocity = [0.015, -0.01];
    initial.angular_velocity = 0.001;
    let disk = DryDisk2::new(&grid, 2.2, initial, 0.1).unwrap();
    (grid, disk)
}

#[test]
fn changing_density_uses_endpoint_half_source_in_the_actual_dry_correction() {
    let (mut grid, mut disk) = forced_fixture();
    let origin = [-5.3, 31.2];
    let old = grid.clone();
    let before = joint(&grid, &disk, origin);
    let radius = disk.radius();
    let initial = disk.body();
    let center = core::array::from_fn::<_, 2, _>(|k| initial.center[k] + 0.5 * initial.velocity[k]);
    let load = [0.001, -0.0006];
    let couple = -0.0007;
    let receipt = disk
        .step(
            &mut grid,
            load,
            couple,
            origin,
            DryDiskCouplingConfig2::default(),
            &gate(),
        )
        .unwrap();
    let after = joint(&grid, &disk, origin);
    let mut exterior_source = [0.; 3];
    let mut interior_source = [0.; 3];
    let mut held_interior = [0.; 3];
    let mut changed_density = 0.0f64;
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            let rho0 = old.f[i].iter().sum::<f64>();
            let rho1 = grid.f[i].iter().sum::<f64>();
            changed_density = changed_density.max((rho1 - rho0).abs());
            let held = core::array::from_fn::<_, 2, _>(|k| rho0 * grid.g[k] + grid.fext[i][k]);
            let average = core::array::from_fn::<_, 2, _>(|k| {
                0.5 * (rho0 + rho1) * grid.g[k] + grid.fext[i][k]
            });
            let moments = disk_cell_moments(radius, [x as f64 - center[0], y as f64 - center[1]]);
            let fraction = moments.area;
            let arm = [x as f64 - origin[0], y as f64 - origin[1]];
            let inside_arm = [
                fraction * arm[0] + moments.first_moment[0],
                fraction * arm[1] + moments.first_moment[1],
            ];
            let full = [average[0], average[1], cross(arm, average)];
            let inside = [
                fraction * average[0],
                fraction * average[1],
                cross(inside_arm, average),
            ];
            let held = [
                fraction * held[0],
                fraction * held[1],
                cross(inside_arm, held),
            ];
            for k in 0..3 {
                exterior_source[k] += full[k] - inside[k];
                interior_source[k] += inside[k];
                held_interior[k] += held[k];
            }
        }
    }
    assert!(
        changed_density > 1e-3,
        "the streamed density must actually change"
    );
    assert!(
        (interior_source[1] - held_interior[1]).abs() > 1e-7,
        "fixture must distinguish the old held-force correction"
    );
    vector_close(receipt.interior_forcing_impulse, interior_source, 2e-13);
    vector_close(
        receipt.represented_fluid_forcing_impulse,
        exterior_source,
        2e-12,
    );
    let authored = [
        load[0],
        load[1],
        couple + cross([center[0] - origin[0], center[1] - origin[1]], load),
    ];
    let periodic = periodic_angular_transport(&grid);
    close(receipt.periodic_angular_impulse, periodic, 1e-12);
    assert_eq!(receipt.boundary_fluid_impulse, [0.; 3]);
    let independent_residual = core::array::from_fn(|k| {
        after[k] - before[k] - authored[k] - exterior_source[k] - if k == 2 { periodic } else { 0. }
    });
    vector_close(independent_residual, [0.; 3], 2e-9);
    // Native held-source work stays separately owned; no physical half-source
    // or fictitious fluid correction is relabelled as authored body work.
    close(
        receipt.physical_external_work,
        load[0] * receipt.native.velocity_midpoint[0]
            + load[1] * receipt.native.velocity_midpoint[1]
            + couple * receipt.native.angular_velocity_midpoint,
        2e-18,
    );
    close(
        receipt.native.external_work,
        receipt.physical_external_work + receipt.interior_correction_work,
        2e-17,
    );
    assert!(receipt.candidate_evaluations >= 5);
}

#[test]
fn forcing_only_fictitious_interior_does_not_become_a_physical_external_load() {
    let mut grid = periodic_grid(24, 24);
    let initial = body(2.4, 0.55, [11.2, 11.6]);
    let mut disk = DryDisk2::new(&grid, 2.4, initial, 0.).unwrap();
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            if disk.solid_fraction(x, y) == 1. {
                let i = grid.idx(x, y);
                grid.fext[i] = [0.0002 * (x as f64 - 11.), -0.0003];
            }
        }
    }
    let origin = [-3., 7.];
    let before = joint(&grid, &disk, origin);
    let receipt = disk
        .step(
            &mut grid,
            [0.; 2],
            0.,
            origin,
            DryDiskCouplingConfig2::default(),
            &gate(),
        )
        .unwrap();
    let after = joint(&grid, &disk, origin);
    vector_close(receipt.represented_fluid_forcing_impulse, [0.; 3], 1e-15);
    assert!(receipt.interior_forcing_impulse[1].abs() > 1e-3);
    assert!(
        receipt
            .interior_correction_impulse
            .iter()
            .any(|x| x.abs() > 1e-6)
    );
    vector_close(
        core::array::from_fn(|k| after[k] - before[k]),
        [0., 0., receipt.periodic_angular_impulse],
        2e-10,
    );
    assert_eq!(receipt.physical_external_work, 0.);
}

#[test]
fn real_probe_budget_cancellation_and_late_drift_refuse_atomically_and_retry_exactly() {
    let (mut grid, mut disk) = forced_fixture();
    let initial_grid = grid.clone();
    let initial_disk = disk.clone();
    let limited = DryDiskCouplingConfig2 {
        maximum_evaluations: 3,
        ..DryDiskCouplingConfig2::default()
    };
    assert_eq!(
        disk.step(
            &mut grid,
            [0.001, -0.0006],
            -0.0007,
            [-5., 31.],
            limited,
            &gate()
        ),
        Err(DryDiskError2::EvaluationBudget)
    );
    assert_eq!(grid.f, initial_grid.f);
    assert_eq!(disk.body(), initial_disk.body());
    assert_eq!(disk.angle(), initial_disk.angle());
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        disk.step(
            &mut grid,
            [0.; 2],
            0.,
            [0.; 2],
            DryDiskCouplingConfig2::default(),
            &cancelled
        ),
        Err(DryDiskError2::Cancelled)
    );
    assert_eq!(grid.f, initial_grid.f);
    assert_eq!(disk.body(), initial_disk.body());
    let mut fresh_grid = initial_grid;
    let mut fresh_disk = initial_disk;
    let actual = disk
        .step(
            &mut grid,
            [0.001, -0.0006],
            -0.0007,
            [-5., 31.],
            DryDiskCouplingConfig2::default(),
            &gate(),
        )
        .unwrap();
    let expected = fresh_disk
        .step(
            &mut fresh_grid,
            [0.001, -0.0006],
            -0.0007,
            [-5., 31.],
            DryDiskCouplingConfig2::default(),
            &gate(),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(grid.f, fresh_grid.f);
    assert_eq!(disk.body(), fresh_disk.body());

    let mut grid = periodic_grid(20, 20);
    let mut initial = body(1.5, 1., [2.001, 10.]);
    initial.velocity = [-0.0005, 0.];
    let mut disk = DryDisk2::new(&grid, 1.5, initial, 0.).unwrap();
    let old = grid.f.clone();
    let error = disk
        .step(
            &mut grid,
            [-0.1, 0.],
            0.,
            [0.; 2],
            DryDiskCouplingConfig2::default(),
            &gate(),
        )
        .unwrap_err();
    assert_eq!(
        error,
        DryDiskError2::Invalid("physical dry disk left the interior lattice domain")
    );
    assert_eq!(grid.f, old);
    assert_eq!(disk.body(), initial);
    assert_eq!(disk.angle(), 0.);
}

#[test]
fn stationary_walls_and_exterior_bounces_keep_their_separate_physical_reaction() {
    let mut grid = Grid::uniform(24, 22, 0.82);
    grid.periodic_x = false;
    grid.periodic_y = false;
    let velocity = [0.012, 0.];
    grid.f
        .fill(fs_lbm::equilibrium(1., velocity[0], velocity[1]));
    for x in 0..grid.nx {
        let bottom = grid.idx(x, 0);
        let top = grid.idx(x, grid.ny - 1);
        grid.flags[bottom] = Cell::Wall;
        grid.flags[top] = Cell::Wall;
    }
    let mut initial = body(1.8, 0.8, [11.3, 10.]);
    initial.velocity = velocity;
    let mut disk = DryDisk2::new(&grid, 1.8, initial, 0.).unwrap();
    let origin = [-9., 23.];
    let before = joint(&grid, &disk, origin);
    // The far boundary sees unperturbed uniform equilibrium. Native BGK
    // supplies its actual outgoing populations independently of disk receipts.
    let mut post = Vec::new();
    grid.collide_into(&mut post);
    let opposite = [0, 3, 4, 1, 2, 7, 8, 5, 6];
    let mut expected = [0.; 3];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            if grid.flags[i] != Cell::Fluid {
                continue;
            }
            for (q, c) in C.into_iter().enumerate() {
                if grid
                    .source(x, y, q)
                    .is_none_or(|s| grid.flags[s] == Cell::Wall)
                {
                    let impulse = [
                        2. * post[i][opposite[q]] * c[0],
                        2. * post[i][opposite[q]] * c[1],
                    ];
                    expected[0] += impulse[0];
                    expected[1] += impulse[1];
                    expected[2] += cross([x as f64 - origin[0], y as f64 - origin[1]], impulse);
                }
            }
        }
    }
    let receipt = disk
        .step(
            &mut grid,
            [0.; 2],
            0.,
            origin,
            DryDiskCouplingConfig2::default(),
            &gate(),
        )
        .unwrap();
    let after = joint(&grid, &disk, origin);
    assert!(expected[0].abs() > 0.1);
    vector_close(receipt.boundary_fluid_impulse, expected, 2e-12);
    vector_close(
        core::array::from_fn(|k| after[k] - before[k]),
        expected,
        2e-10,
    );
    assert_eq!(receipt.periodic_angular_impulse, 0.);
    assert_eq!(receipt.physical_external_impulse, [0.; 3]);
}

#[test]
fn oblique_common_translation_does_not_create_fictitious_spin_at_any_dry_density() {
    let velocity = [0.018, -0.007];
    let origin = [-7., 26.];
    for density in [0.15, 1., 3.] {
        let mut grid = periodic_grid(24, 22);
        grid.f
            .fill(fs_lbm::equilibrium(1., velocity[0], velocity[1]));
        let mut initial = body(1.7, density, [10.2, 11.3]);
        initial.velocity = velocity;
        let mut disk = DryDisk2::new(&grid, 1.7, initial, 0.3).unwrap();
        let conserved = joint(&grid, &disk, origin);
        for step in 1..=16 {
            let exchange = disk
                .step(
                    &mut grid,
                    [0.; 2],
                    0.,
                    origin,
                    DryDiskCouplingConfig2::default(),
                    &gate(),
                )
                .unwrap();
            // Previously the first light-body tick generated O(1e-5) angular
            // correction and changed both velocity components despite zero slip.
            close(disk.body().velocity[0], velocity[0], 3e-14);
            close(disk.body().velocity[1], velocity[1], 3e-14);
            close(disk.body().angular_velocity, 0., 3e-14);
            close(disk.angle(), 0.3, 3e-14);
            vector_close(exchange.interior_correction_impulse, [0.; 3], 1e-12);
            vector_close(joint(&grid, &disk, origin), conserved, 3e-10);
            assert_eq!(exchange.iterations, 0);
            assert_eq!(exchange.candidate_evaluations, 1);
            assert_eq!(disk.body().mass.to_bits(), initial.mass.to_bits());
            assert_eq!(disk.body().inertia.to_bits(), initial.inertia.to_bits());
            for k in 0..2 {
                close(
                    disk.body().center[k],
                    initial.center[k] + step as f64 * velocity[k],
                    2e-13,
                );
            }
        }
    }
}
