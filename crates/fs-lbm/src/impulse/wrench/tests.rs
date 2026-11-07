use super::*;
use crate::d3q19::{E3, Q3, W3};
use crate::impulse::ImpulseSum;

fn controls() -> WrenchImpulseControls3 {
    WrenchImpulseControls3 {
        impulse: ImpulseControls {
            max_cells: 128,
            maximum_speed: 0.1,
            maximum_momentum_defect: 2e-12,
            maximum_work_defect: 2e-12,
        },
        maximum_angular_momentum_defect: 3e-12,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn cloud() -> Vec<WrenchCell3> {
    [[0.0; 3], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
        .into_iter()
        .enumerate()
        .map(|(index, position)| WrenchCell3 {
            index,
            position,
            mass: (index + 1) as f64,
        })
        .collect()
}

fn totals(cells: &[WrenchCell3], impulses: &[[f64; 3]], origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (cell, impulse) in cells.iter().zip(impulses) {
        let r: [f64; 3] = core::array::from_fn(|a| cell.position[a] - origin[a]);
        // Independent component equations for r cross J.
        let moment = [
            r[1] * impulse[2] - r[2] * impulse[1],
            r[2] * impulse[0] - r[0] * impulse[2],
            r[0] * impulse[1] - r[1] * impulse[0],
        ];
        for a in 0..3 {
            linear[a] += impulse[a];
            angular[a] += moment[a];
        }
    }
    (linear, angular)
}

#[test]
fn wrench_collinear_cells_carry_transverse_couple_and_refuse_axial_spin() {
    let cells = [
        WrenchCell3 {
            index: 0,
            position: [0.5, 0.5, 0.5],
            mass: 1.0,
        },
        WrenchCell3 {
            index: 1,
            position: [1.5, 0.5, 0.5],
            mass: 3.0,
        },
    ];
    let gate = CancelGate::new_clock_free();
    let requested = WrenchImpulse3 {
        origin: [0.0; 3],
        linear: [0.0; 3],
        angular: [0.0, 0.0, 0.003],
    };
    let plan = WrenchPlan3::prepare(&cells, requested, controls(), &gate).unwrap();
    let kicks = cells.map(|cell| plan.impulse(cell));
    close(kicks[0][1], -0.003, 1e-17);
    close(kicks[1][1], 0.003, 1e-17);
    let (linear, angular) = totals(&cells, &kicks, requested.origin);
    for a in 0..3 {
        close(linear[a], requested.linear[a], 1e-17);
        close(angular[a], requested.angular[a], 1e-17);
    }
    let kinetic: f64 = cells
        .iter()
        .zip(kicks)
        .map(|(c, j)| 0.5 * j.into_iter().map(|x| x * x).sum::<f64>() / c.mass)
        .sum();
    close(kinetic, 6e-6, 1e-19);
    assert!(matches!(
        WrenchPlan3::prepare(
            &cells,
            WrenchImpulse3 {
                angular: [0.001, 0.0, 0.003],
                ..requested
            },
            controls(),
            &gate,
        ),
        Err(ImpulseError::Invalid(_))
    ));
}

#[test]
fn wrench_off_center_load_is_origin_covariant_and_minimizes_velocity_change() {
    let cells = cloud();
    let gate = CancelGate::new_clock_free();
    let requested = WrenchImpulse3 {
        origin: [0.25, -0.5, 0.75],
        linear: [0.003, -0.002, 0.001],
        angular: [0.001, 0.004, -0.002],
    };
    let plan = WrenchPlan3::prepare(&cells, requested, controls(), &gate).unwrap();
    let kicks: Vec<_> = cells.iter().map(|&c| plan.impulse(c)).collect();
    let (linear, angular) = totals(&cells, &kicks, requested.origin);
    for a in 0..3 {
        close(linear[a], requested.linear[a], 2e-17);
        close(angular[a], requested.angular[a], 2e-17);
    }

    let shift = [5.0, -2.0, 3.0];
    let orbital = cross(shift, requested.linear);
    let translated = WrenchPlan3::prepare(
        &cells,
        WrenchImpulse3 {
            origin: core::array::from_fn(|a| requested.origin[a] + shift[a]),
            angular: core::array::from_fn(|a| requested.angular[a] - orbital[a]),
            ..requested
        },
        controls(),
        &gate,
    )
    .unwrap();
    for (cell, kick) in cells.iter().zip(&kicks) {
        for a in 0..3 {
            close(translated.impulse(*cell)[a], kick[a], 2e-17);
        }
    }

    // A radial expansion about this independently computed mass center has
    // zero linear and angular impulse. It must increase the quadratic change.
    let mass: f64 = cells.iter().map(|c| c.mass).sum();
    let center: [f64; 3] =
        core::array::from_fn(|a| cells.iter().map(|c| c.mass * c.position[a]).sum::<f64>() / mass);
    let alternative: Vec<[f64; 3]> = cells
        .iter()
        .zip(&kicks)
        .map(|(c, j)| core::array::from_fn(|a| j[a] + 0.001 * c.mass * (c.position[a] - center[a])))
        .collect();
    let (alternative_linear, alternative_angular) = totals(&cells, &alternative, requested.origin);
    for a in 0..3 {
        close(alternative_linear[a], linear[a], 3e-17);
        close(alternative_angular[a], angular[a], 3e-17);
    }
    let objective = |field: &[[f64; 3]]| {
        cells
            .iter()
            .zip(field)
            .map(|(c, j)| j.iter().map(|x| x * x).sum::<f64>() / c.mass)
            .sum::<f64>()
    };
    assert!(objective(&alternative) > objective(&kicks) + 1e-6);
}

#[test]
fn wrench_single_cell_accepts_point_force_but_cannot_invent_a_free_couple() {
    let cells = [WrenchCell3 {
        index: 0,
        position: [0.5; 3],
        mass: 2.0,
    }];
    let gate = CancelGate::new_clock_free();
    let origin = [-2.0, 1.0, 0.25];
    let linear = [0.001, -0.002, 0.0005];
    let arm = core::array::from_fn(|a| cells[0].position[a] - origin[a]);
    let angular = cross(arm, linear);
    let request = WrenchImpulse3 {
        origin,
        linear,
        angular,
    };
    let plan = WrenchPlan3::prepare(&cells, request, controls(), &gate).unwrap();
    assert_eq!(plan.impulse(cells[0]), linear);
    let mut coupled = request;
    coupled.angular[0] += 0.001;
    assert!(WrenchPlan3::prepare(&cells, coupled, controls(), &gate).is_err());
}

fn populations(rho: f64, velocity: [f64; 3], force: [f64; 3]) -> [f64; Q3] {
    let uu: f64 = velocity.into_iter().map(|u| u * u).sum();
    core::array::from_fn(|q| {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let eu: f64 = (0..3).map(|a| e[a] * velocity[a]).sum();
        let ef: f64 = (0..3).map(|a| e[a] * force[a]).sum();
        let uf: f64 = (0..3).map(|a| velocity[a] * force[a]).sum();
        let equilibrium = W3[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uu);
        let half_source = 0.5 * W3[q] * (3.0 * (ef - uf) + 9.0 * eu * ef);
        // Nonzero traceless nonequilibrium stress: the kick must not erase it.
        equilibrium - half_source + 1e-5 * W3[q] * (e[0] * e[0] - e[1] * e[1])
    })
}

fn physical_momentum(f: &[f64; Q3], force: [f64; 3]) -> [f64; 3] {
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    core::array::from_fn(|a| {
        f.iter().zip(directions).map(|(v, e)| v * e[a]).sum::<f64>() + 0.5 * force[a]
    })
}

fn stress(f: &[f64; Q3], force: [f64; 3]) -> [f64; 9] {
    let rho: f64 = f.iter().sum();
    let u = physical_momentum(f, force).map(|p| p / rho);
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    core::array::from_fn(|k| {
        let (a, b) = (k / 3, k % 3);
        f.iter()
            .zip(directions)
            .map(|(v, e)| v * e[a] * e[b])
            .sum::<f64>()
            + 0.5 * (u[a] * force[b] + u[b] * force[a])
            - rho * u[a] * u[b]
            - if a == b { rho / 3.0 } else { 0.0 }
    })
}

#[test]
fn wrench_population_impulse_preserves_guo_stress_and_galilean_work() {
    let gate = CancelGate::new_clock_free();
    let requested = WrenchImpulse3 {
        origin: [-0.25, 0.5, -0.75],
        linear: [0.003, -0.002, 0.001],
        angular: [0.001, 0.004, -0.002],
    };
    let boost = [0.02, -0.01, 0.015];
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let mut work = [0.0; 2];
    for (run, velocity) in [[0.0; 3], boost].into_iter().enumerate() {
        let mut cells = cloud();
        let forces: Vec<_> = (0..cells.len())
            .map(|i| [1e-5 * (i as f64 + 1.0), -2e-5, 3e-5])
            .collect();
        let before: Vec<_> = cells
            .iter()
            .zip(&forces)
            .map(|(c, &f)| populations(c.mass, velocity, f))
            .collect();
        for (cell, pop) in cells.iter_mut().zip(&before) {
            cell.mass = pop.iter().sum();
        }
        let plan = WrenchPlan3::prepare(&cells, requested, controls(), &gate).unwrap();
        let mut after = before.clone();
        let mut sum = ImpulseSum::new();
        for (i, &cell) in cells.iter().enumerate() {
            sum.kick(
                &mut after[i],
                &directions,
                &W3,
                [0.0; 3],
                forces[i],
                plan.impulse(cell),
                controls().impulse,
                i,
            )
            .unwrap();
        }
        let report = plan
            .finish(
                &cells,
                sum.finish(controls().impulse).unwrap(),
                controls(),
                &gate,
                |i| {
                    let old = physical_momentum(&before[i], forces[i]);
                    let new = physical_momentum(&after[i], forces[i]);
                    core::array::from_fn(|a| new[a] - old[a])
                },
            )
            .unwrap();
        for a in 0..3 {
            close(
                report.impulse.momentum_change[a],
                requested.linear[a],
                2e-14,
            );
            close(
                report.angular_momentum_change[a],
                requested.angular[a],
                5e-14,
            );
        }
        close(
            report.impulse.kinetic_energy_change,
            report.impulse.work,
            2e-14,
        );
        close(report.impulse.mass_change, 0.0, 2e-14);
        for i in 0..cells.len() {
            let old = stress(&before[i], forces[i]);
            let new = stress(&after[i], forces[i]);
            for k in 0..9 {
                close(new[k], old[k], 3e-15);
            }
        }
        work[run] = report.impulse.work;
    }
    let observer_work: f64 = (0..3).map(|a| boost[a] * requested.linear[a]).sum();
    close(work[1] - work[0], observer_work, 3e-15);
}

#[test]
fn wrench_rejects_bad_inertia_duplicate_cells_and_cancellation() {
    let gate = CancelGate::new_clock_free();
    let request = WrenchImpulse3 {
        origin: [0.0; 3],
        linear: [0.001; 3],
        angular: [0.0; 3],
    };
    let mut invalid = cloud();
    invalid[1].index = invalid[0].index;
    assert!(WrenchPlan3::prepare(&invalid, request, controls(), &gate).is_err());
    invalid = cloud();
    invalid[2].mass = 0.0;
    assert!(WrenchPlan3::prepare(&invalid, request, controls(), &gate).is_err());
    invalid = cloud();
    invalid[2].position[1] = f64::NAN;
    assert!(WrenchPlan3::prepare(&invalid, request, controls(), &gate).is_err());
    gate.request();
    assert!(matches!(
        WrenchPlan3::prepare(&cloud(), request, controls(), &gate),
        Err(ImpulseError::Cancelled)
    ));
}
