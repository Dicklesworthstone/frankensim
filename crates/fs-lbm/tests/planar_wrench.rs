//! Independent physical quadratures for total planar force-and-couple loading.
use fs_exec::CancelGate;
use fs_lbm::{
    CS2, Cell, ContactModel, FreeSurface, Grid, Q, equilibrium,
    impulse::{
        ImpulseControls, ImpulseError, WrenchImpulse2, WrenchImpulseControls2, WrenchImpulseReport2,
    },
};

const DIRECTIONS: [[f64; 2]; Q] = [
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

fn controls() -> WrenchImpulseControls2 {
    WrenchImpulseControls2 {
        impulse: ImpulseControls {
            max_cells: 6,
            maximum_speed: 0.15,
            maximum_momentum_defect: 1e-12,
            maximum_work_defect: 1e-13,
        },
        maximum_angular_momentum_defect: 1e-12,
    }
}

fn fixture(boost: [f64; 2]) -> Grid {
    let mut grid = Grid::uniform(3, 2, 0.9);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [0.0002, -0.0001];
    for i in 0..6 {
        grid.f[i] = equilibrium(0.8 + 0.17 * i as f64, 0.025 + boost[0], -0.018 + boost[1]);
        // Distinct normal/shear viscous stresses, without changing raw momentum.
        grid.f[i][0] -= 0.0006;
        grid.f[i][1] += 0.0003;
        grid.f[i][3] += 0.0003;
        grid.f[i][5] += 0.0002;
        grid.f[i][7] += 0.0002;
        grid.f[i][6] -= 0.0002;
        grid.f[i][8] -= 0.0002;
        grid.fext[i] = [0.00001 * i as f64, -0.00002 * (i + 1) as f64];
    }
    grid
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn position(grid: &Grid, i: usize) -> [f64; 2] {
    [(i % grid.nx) as f64 + 0.5, (i / grid.nx) as f64 + 0.5]
}

// Independent quadrature, including the Guo half force and physical viscous
// stress. Does not call the production impulse or momentum observer.
fn observe(grid: &Grid, i: usize) -> (f64, [f64; 2], [[f64; 2]; 2]) {
    let rho: f64 = grid.f[i].iter().sum();
    let force: [f64; 2] = core::array::from_fn(|a| rho * grid.g[a] + grid.fext[i][a]);
    let velocity: [f64; 2] = core::array::from_fn(|a| {
        (grid.f[i]
            .iter()
            .zip(DIRECTIONS)
            .map(|(f, e)| f * e[a])
            .sum::<f64>()
            + 0.5 * force[a])
            / rho
    });
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            grid.f[i]
                .iter()
                .zip(DIRECTIONS)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                - rho * velocity[a] * velocity[b]
                - if a == b { rho * CS2 } else { 0.0 }
                + 0.5 * (velocity[a] * force[b] + velocity[b] * force[a])
        })
    });
    (rho, velocity, stress)
}

fn check_physics(
    before: &Grid,
    after: &Grid,
    cells: &[(usize, f64)],
    wrench: WrenchImpulse2,
    report: &WrenchImpulseReport2,
) {
    let stencil_mass = cells.iter().map(|(_, mass)| mass).sum::<f64>();
    let center: [f64; 2] = core::array::from_fn(|a| {
        cells
            .iter()
            .map(|&(i, mass)| mass * position(before, i)[a])
            .sum::<f64>()
            / stencil_mass
    });
    let inertia = cells
        .iter()
        .map(|&(i, mass)| {
            let p = position(before, i);
            mass * ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2))
        })
        .sum::<f64>();
    let center_arm = core::array::from_fn(|a| center[a] - wrench.origin[a]);
    let rotation = (wrench.angular - cross(center_arm, wrench.linear)) / inertia;
    let mut impulse = [0.0; 2];
    let mut angular = 0.0;
    let mut kinetic = 0.0;
    for &(i, mass) in cells {
        let old = observe(before, i);
        let new = observe(after, i);
        let r = position(before, i);
        let expected = [
            wrench.linear[0] / stencil_mass - rotation * (r[1] - center[1]),
            wrench.linear[1] / stencil_mass + rotation * (r[0] - center[0]),
        ];
        let delta = core::array::from_fn(|a| mass * (new.1[a] - old.1[a]));
        close(new.0, old.0, 2e-15);
        for a in 0..2 {
            close(new.1[a] - old.1[a], expected[a], 2e-15);
            impulse[a] += delta[a];
            kinetic += 0.5 * mass * (new.1[a] * new.1[a] - old.1[a] * old.1[a]);
            for b in 0..2 {
                close(new.2[a][b], old.2[a][b], 2e-15);
            }
        }
        angular += cross(core::array::from_fn(|a| r[a] - wrench.origin[a]), delta);
    }
    for i in 0..before.f.len() {
        if cells.iter().all(|&(cell, _)| cell != i) {
            assert_eq!(before.f[i].map(f64::to_bits), after.f[i].map(f64::to_bits));
        }
    }
    close(report.stencil_mass, stencil_mass, 2e-15);
    for a in 0..2 {
        close(report.stencil_center_of_mass[a], center[a], 2e-15);
        close(impulse[a], wrench.linear[a], 2e-15);
        close(report.impulse.momentum_change[a], impulse[a], 2e-15);
    }
    close(angular, wrench.angular, 3e-15);
    close(report.angular_momentum_change, angular, 3e-15);
    close(report.impulse.work, kinetic, 2e-16);
    close(report.impulse.kinetic_energy_change, kinetic, 2e-16);
    assert_eq!(before.flags, after.flags);
    assert_eq!(before.g, after.g);
    assert_eq!(before.fext, after.fext);
    assert_eq!(before.tau, after.tau);
}

#[test]
fn dense_wrench_realizes_translation_and_couple_in_retained_guo_frame() {
    let mut grid = fixture([0.0; 2]);
    let before = grid.clone();
    let indices = [0, 2, 5];
    let cells: Vec<_> = indices
        .iter()
        .map(|&i| (i, observe(&before, i).0))
        .collect();
    let wrench = WrenchImpulse2 {
        origin: [-0.75, 0.2],
        linear: [0.012, -0.006],
        angular: 0.02,
    };
    let report = grid
        .apply_wrench_impulse(
            &indices,
            wrench,
            controls(),
            &mut Vec::new(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    check_physics(&before, &grid, &cells, wrench, &report);
}

#[test]
fn free_surface_wrench_uses_partial_liquid_inertia_and_retains_inventory() {
    let mut liquid = FreeSurface::new(fixture([0.0; 2]), 0.0, ContactModel::Neutral);
    liquid.grid.flags[1] = Cell::Interface;
    liquid.grid.flags[5] = Cell::Interface;
    liquid.mass[1] = 0.25 * observe(&liquid.grid, 1).0;
    liquid.mass[5] = 0.6 * observe(&liquid.grid, 5).0;
    let before = liquid.clone();
    let cells = [
        (0, observe(&before.grid, 0).0),
        (1, before.mass[1]),
        (5, before.mass[5]),
    ];
    let wrench = WrenchImpulse2 {
        origin: [0.2, -0.4],
        linear: [-0.01, 0.006],
        angular: 0.015,
    };
    let report = liquid
        .apply_wrench_impulse(
            &[0, 1, 5],
            wrench,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    check_physics(&before.grid, &liquid.grid, &cells, wrench, &report);
    assert_eq!(liquid.mass, before.mass);
    assert_eq!(liquid.carry, before.carry);
    for i in [0, 1, 5] {
        assert_eq!(
            observe(&liquid.grid, i).0.to_bits(),
            observe(&before.grid, i).0.to_bits()
        );
    }
    assert_eq!(report.impulse.mass_change, 0.0);
}

#[test]
fn wrench_is_origin_covariant_and_work_has_the_galilean_shift() {
    let gate = CancelGate::new_clock_free();
    let indices = [0, 2, 5];
    let wrench = WrenchImpulse2 {
        origin: [-0.75, 0.2],
        linear: [0.012, -0.006],
        angular: 0.02,
    };
    let mut first = fixture([0.0; 2]);
    let mut shifted = first.clone();
    let mut boosted = fixture([0.03, -0.01]);
    let original = first.clone();
    let boosted_original = boosted.clone();
    let origin = [3.25, -1.8];
    let shifted_wrench = WrenchImpulse2 {
        origin,
        angular: wrench.angular
            + cross(
                core::array::from_fn(|a| wrench.origin[a] - origin[a]),
                wrench.linear,
            ),
        ..wrench
    };
    let report = first
        .apply_wrench_impulse(&indices, wrench, controls(), &mut Vec::new(), &gate)
        .unwrap();
    let shifted_report = shifted
        .apply_wrench_impulse(&indices, shifted_wrench, controls(), &mut Vec::new(), &gate)
        .unwrap();
    close(
        shifted_report.angular_momentum_change,
        shifted_wrench.angular,
        3e-15,
    );
    let boost_report = boosted
        .apply_wrench_impulse(&indices, wrench, controls(), &mut Vec::new(), &gate)
        .unwrap();
    for i in 0..6 {
        for q in 0..Q {
            close(first.f[i][q], shifted.f[i][q], 2e-16);
        }
        for a in 0..2 {
            close(
                observe(&first, i).1[a] - observe(&original, i).1[a],
                observe(&boosted, i).1[a] - observe(&boosted_original, i).1[a],
                2e-15,
            );
        }
    }
    close(
        boost_report.impulse.work - report.impulse.work,
        0.03 * wrench.linear[0] - 0.01 * wrench.linear[1],
        2e-16,
    );
}

#[test]
fn a_single_cell_carries_orbital_momentum_but_refuses_a_free_couple() {
    let gate = CancelGate::new_clock_free();
    let mut grid = fixture([0.0; 2]);
    let before = grid.clone();
    let linear = [-0.01, 0.006];
    let origin = [-0.5, -0.5];
    let arm = core::array::from_fn(|a| position(&grid, 2)[a] - origin[a]);
    let wrench = WrenchImpulse2 {
        origin,
        linear,
        angular: cross(arm, linear),
    };
    let report = grid
        .apply_wrench_impulse(&[2], wrench, controls(), &mut Vec::new(), &gate)
        .unwrap();
    assert_eq!(report.impulse.affected_cells, 1);
    assert!(report.impulse.work < 0.0, "braking removes kinetic energy");
    close(report.angular_momentum_change, wrench.angular, 2e-15);
    close(
        report.impulse.work,
        report.impulse.kinetic_energy_change,
        2e-16,
    );
    let accepted = grid.f.clone();
    assert!(
        grid.apply_wrench_impulse(
            &[2],
            WrenchImpulse2 {
                angular: wrench.angular + 0.001,
                ..wrench
            },
            controls(),
            &mut Vec::new(),
            &gate,
        )
        .is_err()
    );
    assert_eq!(grid.f, accepted);
    assert_ne!(grid.f, before.f);
}

#[test]
fn planar_wrench_refusal_and_cancellation_never_publish_a_prefix() {
    let gate = CancelGate::new_clock_free();
    let wrench = WrenchImpulse2 {
        origin: [0.0; 2],
        linear: [0.012, -0.006],
        angular: 0.02,
    };
    let mut grid = fixture([0.0; 2]);
    let before = grid.f.clone();
    for indices in [&[][..], &[0, 0], &[2, 0], &[0, 6]] {
        assert!(
            grid.apply_wrench_impulse(indices, wrench, controls(), &mut Vec::new(), &gate,)
                .is_err()
        );
        assert_eq!(grid.f, before);
    }
    grid.periodic_y = true;
    assert!(
        grid.apply_wrench_impulse(&[0, 2, 5], wrench, controls(), &mut Vec::new(), &gate,)
            .is_err()
    );
    grid.periodic_y = false;
    assert!(
        grid.apply_wrench_impulse(
            &[0, 2, 5],
            WrenchImpulse2 {
                angular: 10.0,
                ..wrench
            },
            controls(),
            &mut Vec::new(),
            &gate,
        )
        .is_err()
    );
    assert_eq!(grid.f, before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        grid.apply_wrench_impulse(&[0, 2, 5], wrench, controls(), &mut Vec::new(), &cancelled,),
        Err(ImpulseError::Cancelled)
    );
    assert_eq!(grid.f, before);

    let mut liquid = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    liquid.grid.flags[5] = Cell::Interface;
    for invalid_mass in [0.0, f64::NAN, 2.0 * observe(&liquid.grid, 5).0] {
        liquid.mass[5] = invalid_mass;
        let before = liquid.grid.f.clone();
        assert!(
            liquid
                .apply_wrench_impulse(&[0, 2, 5], wrench, controls(), &gate)
                .is_err()
        );
        assert_eq!(liquid.grid.f, before);
    }
}

#[test]
fn zero_planar_wrench_preserves_every_population_bit() {
    let gate = CancelGate::new_clock_free();
    let zero = WrenchImpulse2 {
        origin: [0.0; 2],
        linear: [0.0; 2],
        angular: 0.0,
    };
    let mut grid = fixture([0.0; 2]);
    let before = grid.f.clone();
    let report = grid
        .apply_wrench_impulse(&[0, 2, 5], zero, controls(), &mut Vec::new(), &gate)
        .unwrap();
    assert_eq!(report.impulse.affected_cells, 0);
    assert_eq!(report.impulse.work, 0.0);
    assert_eq!(grid.f, before);
    let mut liquid = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let report = liquid
        .apply_wrench_impulse(&[0, 2, 5], zero, controls(), &gate)
        .unwrap();
    assert_eq!(report.impulse.affected_cells, 0);
    assert_eq!(liquid.grid.f, before);
}
