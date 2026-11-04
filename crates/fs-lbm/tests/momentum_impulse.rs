//! G0/G3: actual momentum jumps, viscous stress, kinetic work and SI invariance.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{CollisionModel3, E3, Q3};
use fs_lbm::impulse::{ImpulseControls, ImpulseError};
use fs_lbm::material::plan_scaling_from_material;
use fs_lbm::{BoundaryGrid3, BoundarySpec3, CS2, Cell, Grid, Q, equilibrium, equilibrium3};
use fs_material::fluid_transport::FluidTransportState;

const E2: [[f64; 2]; Q] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];

fn controls() -> ImpulseControls {
    ImpulseControls {
        max_cells: 1024,
        maximum_speed: 0.15,
        maximum_momentum_defect: 1e-11,
        maximum_work_defect: 1e-12,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual {a:.17e}, expected {b:.17e}"
    );
}

// Independent quadrature of the force-corrected moments; no production kick,
// moment observer, or equilibrium difference is used as the oracle.
fn moments<const D: usize, const N: usize>(
    f: [f64; N],
    e: [[f64; D]; N],
    force: [f64; D],
) -> (f64, [f64; D], [[f64; D]; D], f64) {
    let rho = f.iter().sum::<f64>();
    let p = core::array::from_fn(|a| {
        f.iter().zip(e).map(|(f, e)| f * e[a]).sum::<f64>() + force[a] / 2.0
    });
    let u: [f64; D] = p.map(|p| p / rho);
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter().zip(e).map(|(f, e)| f * e[a] * e[b]).sum::<f64>()
                - rho * u[a] * u[b]
                - if a == b { rho * CS2 } else { 0.0 }
                + (u[a] * force[b] + u[b] * force[a]) / 2.0
        })
    });
    (rho, p, stress, p.iter().map(|p| p * p / (2.0 * rho)).sum())
}
fn moments2(grid: &Grid, i: usize) -> (f64, [f64; 2], [[f64; 2]; 2], f64) {
    let rho = grid.f[i].iter().sum::<f64>();
    moments(
        grid.f[i],
        E2,
        core::array::from_fn(|a| rho * grid.g[a] + grid.fext[i][a]),
    )
}
fn coords3(i: usize) -> [usize; 3] {
    [i % 8, (i / 8) % 4, i / 32]
}
fn moments3(grid: &BoundaryGrid3, i: usize) -> (f64, [f64; 3], [[f64; 3]; 3], f64) {
    let [x, y, z] = coords3(i);
    moments(
        grid.populations(x, y, z),
        E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]),
        grid.force_density(x, y, z),
    )
}
fn bits3(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    (0..128)
        .map(|i| {
            let [x, y, z] = coords3(i);
            grid.populations(x, y, z).map(f64::to_bits)
        })
        .collect()
}

#[test]
fn dense_impulse_preserves_viscous_stress_and_reports_signed_work() {
    let mut grid = Grid::uniform(4, 3, 0.9);
    let mut scratch = Vec::new();
    for i in 0..12 {
        grid.f[i] = equilibrium(0.7 + 0.11 * i as f64, 0.035, -0.025);
        grid.f[i][1] += 0.001;
        grid.f[i][3] += 0.001;
        grid.f[i][0] -= 0.002;
        grid.f[i][5] += 0.0003;
        grid.f[i][7] += 0.0003;
        grid.f[i][6] -= 0.0003;
        grid.f[i][8] -= 0.0003;
    }
    grid.flags[0] = Cell::Wall;
    grid.flags[1] = Cell::Gas;
    grid.flags[2] = Cell::Interface;
    grid.rebind_forcing([0.002, -0.001], &vec![[0.0007, -0.0003]; 12], &mut scratch)
        .unwrap();
    let before = grid.clone();
    let impulses: Vec<_> = (0..12)
        .map(|i| if i < 3 { [0.0; 2] } else { [-0.01, 0.008] })
        .collect();
    let report = grid
        .apply_impulse_density(
            &impulses,
            controls(),
            &mut scratch,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(report.affected_cells, 9);
    assert!(
        report.work < 0.0,
        "a braking impulse removes kinetic energy"
    );
    let mut delta_k = 0.0;
    let mut work = 0.0;
    for i in 0..12 {
        if i < 3 {
            assert_eq!(grid.f[i], before.f[i]);
            continue;
        }
        let old = moments2(&before, i);
        let new = moments2(&grid, i);
        close(new.0, old.0, 8e-16);
        for a in 0..2 {
            close(new.1[a] - old.1[a], impulses[i][a], 2e-16);
            for b in 0..2 {
                close(new.2[a][b], old.2[a][b], 8e-16);
            }
            work += impulses[i][a] * (old.1[a] + 0.5 * impulses[i][a]) / old.0;
        }
        delta_k += new.3 - old.3;
    }
    close(report.kinetic_energy_change, delta_k, 2e-17);
    close(report.work, work, 2e-17);
    close(report.work, report.kinetic_energy_change, 2e-17);
    assert_eq!(grid.g, before.g);
    assert_eq!(grid.fext, before.fext);
    assert_eq!(grid.flags, before.flags);
    assert_eq!(grid.tau, before.tau);
}

#[test]
fn tiled_nonuniform_impulses_use_original_cells_and_preserve_full_stress() {
    let gate = CancelGate::new_clock_free();
    let mut grid = BoundaryGrid3::new(8, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    grid.voxelize_sdf(|p| (p[0] - 3.5).abs() - 0.25);
    for i in 0..128 {
        let [x, y, z] = coords3(i);
        if grid.is_solid(x, y, z) {
            continue;
        }
        let mut f = equilibrium3(1.2 + 0.02 * y as f64, [0.015, -0.02, 0.01]);
        f[1] += 0.002;
        f[2] += 0.002;
        f[0] -= 0.004;
        f[7] += 0.001;
        f[8] += 0.001;
        f[9] -= 0.001;
        f[10] -= 0.001;
        grid.set_populations(x, y, z, &f);
    }
    grid.rebind_uniform_force([0.0002, -0.0003, 0.0001], &gate)
        .unwrap();
    grid.rebind_force_field_with_gate(
        |[x, y, z]| [x as f64 * 1e-5, y as f64 * -2e-5, z as f64 * 3e-5],
        &gate,
    )
    .unwrap();
    let original_bits = bits3(&grid);
    let before: Vec<_> = (0..128)
        .map(|i| {
            let [x, y, z] = coords3(i);
            (!grid.is_solid(x, y, z)).then(|| moments3(&grid, i))
        })
        .collect();
    let impulses: Vec<_> = (0..128)
        .map(|i| {
            let [x, y, z] = coords3(i);
            if grid.is_solid(x, y, z) {
                [0.0; 3]
            } else {
                [
                    0.001 * (x + 1) as f64,
                    -0.002 * (y + 1) as f64,
                    0.003 * (z + 1) as f64,
                ]
            }
        })
        .collect();
    let report = grid
        .apply_impulse_density(&impulses, controls(), &gate)
        .unwrap();
    let mut expected = [0.0; 3];
    let mut delta_k = 0.0;
    for i in 0..128 {
        let [x, y, z] = coords3(i);
        if grid.is_solid(x, y, z) {
            assert_eq!(
                grid.populations(x, y, z).map(f64::to_bits),
                original_bits[i]
            );
            continue;
        }
        let old = before[i].unwrap();
        let new = moments3(&grid, i);
        close(new.0, old.0, 2e-15);
        for a in 0..3 {
            expected[a] += impulses[i][a];
            close(new.1[a] - old.1[a], impulses[i][a], 4e-16);
            for b in 0..3 {
                close(new.2[a][b], old.2[a][b], 2e-15);
            }
        }
        delta_k += new.3 - old.3;
    }
    for a in 0..3 {
        close(report.applied_impulse[a], expected[a], 3e-15);
        close(report.momentum_change[a], expected[a], 8e-15);
    }
    close(report.kinetic_energy_change, delta_k, 1e-16);
    close(report.work, delta_k, 3e-16);
}

#[test]
fn impulse_is_not_reapplied_by_subsequent_native_collision() {
    let gate = CancelGate::new_clock_free();
    for model in [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 0.9,
        },
    ] {
        let mut grid = BoundaryGrid3::with_collision_model(
            8,
            4,
            4,
            model,
            [0.0; 3],
            BoundarySpec3::periodic(),
        );
        grid.initialize_uniform_equilibrium(1.4, [0.011, -0.009, 0.006], &gate)
            .unwrap();
        let force = [0.00002, -0.00001, 0.00003];
        grid.rebind_uniform_force(force, &gate).unwrap();
        let p0 = moments3(&grid, 0).1;
        let impulse = [0.014, -0.007, 0.021];
        grid.apply_impulse_density(&vec![impulse; 128], controls(), &gate)
            .unwrap();
        for step in 0..4 {
            let p = moments3(&grid, 0).1;
            for a in 0..3 {
                close(p[a], p0[a] + impulse[a] + step as f64 * force[a], 2e-15);
            }
            grid.step();
        }
    }
}

#[test]
fn refusals_and_zero_fields_leave_original_populations_untouched() {
    let gate = CancelGate::new_clock_free();
    let mut grid = Grid::uniform(4, 3, 0.9);
    let before = grid.f.clone();
    let mut scratch = Vec::new();
    let mut impulses = vec![[0.001, 0.0]; 12];
    for bad in [f64::NAN, f64::INFINITY, 2.0] {
        impulses[11][0] = bad;
        assert!(
            grid.apply_impulse_density(&impulses, controls(), &mut scratch, &gate)
                .is_err()
        );
        assert_eq!(grid.f, before);
    }
    impulses[11][0] = 0.001;
    grid.flags[11] = Cell::Interface;
    assert!(
        grid.apply_impulse_density(&impulses, controls(), &mut scratch, &gate)
            .is_err()
    );
    assert_eq!(grid.f, before);
    grid.flags[11] = Cell::Fluid;
    let strict = ImpulseControls {
        maximum_momentum_defect: 0.0,
        maximum_work_defect: 0.0,
        ..controls()
    };
    assert_eq!(
        grid.apply_impulse_density(&impulses, strict, &mut scratch, &gate),
        Err(ImpulseError::Balance)
    );
    assert_eq!(grid.f, before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        grid.apply_impulse_density(&impulses, controls(), &mut scratch, &cancelled),
        Err(ImpulseError::Cancelled)
    );
    assert_eq!(grid.f, before);
    let report = grid
        .apply_impulse_density(&vec![[0.0; 2]; 12], strict, &mut scratch, &gate)
        .unwrap();
    assert_eq!(report.affected_cells, 0);
    assert_eq!(report.work, 0.0);
    assert_eq!(grid.f, before);

    let mut grid = BoundaryGrid3::new(8, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    let before = bits3(&grid);
    let mut impulses = vec![[0.001, 0.0, 0.0]; 128];
    impulses[127][2] = f64::NAN;
    assert!(
        grid.apply_impulse_density(&impulses, controls(), &gate)
            .is_err()
    );
    assert_eq!(bits3(&grid), before);
    assert_eq!(
        grid.apply_impulse_density(&vec![[0.001; 3]; 128], controls(), &cancelled),
        Err(ImpulseError::Cancelled)
    );
    assert_eq!(bits3(&grid), before);
    assert!(
        grid.apply_impulse_density(&vec![[0.001; 3]; 128], strict, &gate)
            .is_err()
    );
    assert_eq!(bits3(&grid), before);
    grid.apply_impulse_density(&vec![[0.0; 3]; 128], strict, &gate)
        .unwrap();
    assert_eq!(bits3(&grid), before);
}

#[test]
fn admitted_impulse_locks_original_solid_topology() {
    let mut grid = BoundaryGrid3::new(8, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    grid.apply_impulse_density(
        &vec![[0.001, 0.0, 0.0]; 128],
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let before = bits3(&grid);
    // Removing an already accelerated cell would silently remove its mass
    // and impulse. Geometry changes require the original topology owner.
    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            grid.voxelize_sdf(|point| point[0] - 1.0);
        }))
        .is_err()
    );
    assert_eq!(bits3(&grid), before);
}

#[test]
fn material_impulse_has_same_si_velocity_jump_for_different_dx_dt_and_density() {
    let gate = CancelGate::new_clock_free();
    for density in [1.2, 1000.0] {
        let fluid = FluidTransportState::try_new(density, 1000.0, 0.026, 1.8e-5).unwrap();
        for (dx, dt) in [(0.001, 0.001), (0.0005, 0.0001), (0.002, 0.0002)] {
            let scaling = plan_scaling_from_material(&fluid, dx, dt, 0.1).unwrap();
            let mut grid = scaling.flow_grid(2, 2, 4).unwrap();
            // Actual density is twice the material reference. The same physical
            // impulse must therefore give half the reference-density velocity.
            grid.f.fill(equilibrium(2.0, 0.0, 0.0));
            let j_si = density * 0.08;
            let j = scaling.impulse_density_lattice(j_si).unwrap();
            close(scaling.impulse_density_si(j).unwrap(), j_si, 1e-13);
            let report = grid
                .apply_impulse_density(&[[j, 0.0]; 4], controls(), &mut Vec::new(), &gate)
                .unwrap();
            close(
                scaling.velocity_si(grid.moments(0).u[0]).unwrap(),
                0.04,
                2e-15,
            );
            close(
                scaling
                    .impulse_density_si(report.applied_impulse[0])
                    .unwrap(),
                4.0 * j_si,
                4e-13,
            );
            let velocity_scale = dx / dt;
            let work_per_depth = report.work * density * velocity_scale * velocity_scale * dx * dx;
            close(work_per_depth, 4.0 * dx * dx * j_si * 0.02, 2e-17);
        }
    }
}
