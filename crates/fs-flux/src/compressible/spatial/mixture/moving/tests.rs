use super::*;
use crate::compressible::AxisBoundaries;
use crate::compressible::mixture::tests::{close, gas, source, standard};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

mod forced;

const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const WALLS: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};
const OPENS: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: OPEN,
    z: OPEN,
};

fn field(
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    initial: impl Fn([f64; 3]) -> ([f64; 3], [f64; 3], f64),
) -> MixtureEuler3d {
    let mut states = vec![0.0; grid.len() * 7];
    for (k, row) in states.chunks_exact_mut(7).enumerate() {
        let (species, velocity, pressure) = initial(grid.cell_center(k).unwrap());
        state_into(&gas, &species, velocity, pressure, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}

fn balance(report: StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    assert_eq!(report.source_change, ForceChange3d::default());
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 5e-12);
    }
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 5e-12);
    }
    close(report.conservation_residual.energy, 0.0, 6e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn g0_zero_motion_retains_original_species_step_and_poll_sequence() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let sample = |x: [f64; 3]| {
        (
            [0.6 + 0.1 * x[0], 0.4 - 0.05 * x[0], 0.0],
            [0.2, -0.1, 0.3],
            1.0 + 0.2 * x[1],
        )
    };
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut ordinary = field(standard(0.4), grid, PERIODIC, sample);
        let mut moving = field(standard(0.4), grid, PERIODIC, sample);
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        let mut a_polls = 0;
        let a = ordinary
            .step(Time::new(0.01), controls, || {
                a_polls += 1;
                false
            })
            .unwrap();
        let mut b_polls = 0;
        let b = moving
            .step_moving(Time::new(0.01), controls, MeshMotion::default(), || {
                b_polls += 1;
                false
            })
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(a_polls, b_polls);
        assert_eq!(ordinary.states(), moving.states());
        assert_eq!(ordinary.grid(), moving.grid());
        assert_eq!(ordinary.time(), moving.time());
        assert_eq!(ordinary.boundary_impulses, moving.boundary_impulses);
    }
}

#[test]
fn g1_species_free_stream_obeys_gcl_for_every_piston_axis() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let densities = [0.6, 0.35, 0.0];
    let velocity = [0.3, -0.2, 0.1];
    for axis in 0..3 {
        for dilation in [-0.4, 0.4] {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let mut flow = field(standard(0.4), grid, OPENS, |_| (densities, velocity, 1.7));
                let expected = flow.states()[..7].to_vec();
                let mut motion = MeshMotion {
                    lower_velocity: [-0.1, 0.15, 0.05],
                    upper_velocity: [-0.1, 0.15, 0.05],
                };
                motion.upper_velocity[axis] += dilation;
                for _ in 0..3 {
                    let previous = *flow.grid();
                    let report = flow
                        .step_moving(
                            Time::new(0.005),
                            StepControl {
                                reconstruction,
                                ..StepControl::default()
                            },
                            motion,
                            || false,
                        )
                        .unwrap();
                    balance(report);
                    let dt = report.elapsed.value();
                    let [a, b] = tangents(axis);
                    let swept = dt * dilation * previous.spacing[a] * previous.spacing[b]
                        / previous.shape[axis] as f64;
                    close(
                        flow.grid().cell_volume() - previous.cell_volume(),
                        swept,
                        3e-15,
                    );
                    for d in 0..3 {
                        close(
                            flow.grid().origin()[d],
                            previous.origin()[d] + dt * motion.lower_velocity[d],
                            3e-15,
                        );
                    }
                    for row in flow.states().chunks_exact(7) {
                        for (actual, expected) in row.iter().zip(&expected) {
                            close(*actual, *expected, 3e-12);
                        }
                        assert_eq!(row[2], 0.0);
                    }
                    let change = grid.len() as f64 * swept;
                    for (d, expected) in expected.iter().enumerate() {
                        close(
                            report.boundary_change.component(d),
                            change * expected,
                            4e-12,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g1_comoving_unequal_caloric_contacts_and_both_shears_remain_stationary() {
    let grid = CartesianGrid::new(8, 2, 2, 0.125, 0.3, 0.4).unwrap();
    let speed = [0.3, -0.2, 0.1];
    let motion = MeshMotion {
        lower_velocity: speed,
        upper_velocity: speed,
    };
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut flow = field(standard(0.4), grid, PERIODIC, |x| {
            let a = (std::f64::consts::TAU * x[0]).sin();
            (
                [0.6 + 0.2 * a, 0.4 - 0.1 * a, 0.0],
                [speed[0], speed[1] + 0.1 * a, speed[2] - 0.15 * a],
                1.0,
            )
        });
        let original = flow.states().to_vec();
        let totals = flow.totals();
        for _ in 0..4 {
            balance(
                flow.step_moving(
                    Time::new(0.01),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    motion,
                    || false,
                )
                .unwrap(),
            );
        }
        for (actual, expected) in flow.states().iter().zip(&original) {
            close(*actual, *expected, 5e-12);
        }
        for k in 0..flow.cell_count() {
            let q = flow.cell(k).unwrap();
            assert_eq!(q.partial_densities[2], 0.0);
            close(
                q.primitive(flow.gas(), || false).unwrap().pressure,
                1.0,
                5e-12,
            );
        }
        for d in 0..7 {
            close(flow.totals().component(d), totals.component(d), 5e-12);
        }
    }
}

#[test]
fn g1_piston_preserves_each_species_and_exchanges_the_same_impulse_and_work() {
    let grid = CartesianGrid::new(6, 2, 2, 0.1, 0.3, 0.4).unwrap();
    for axis in 0..3 {
        let mut flow = field(standard(0.4), grid, WALLS, |_| {
            ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
        });
        let initial = flow.totals();
        let mut motion = MeshMotion::default();
        motion.upper_velocity[axis] = -0.1;
        for _ in 0..4 {
            let report = flow
                .step_moving(Time::new(0.002), StepControl::default(), motion, || false)
                .unwrap();
            balance(report);
            for species in 0..3 {
                assert_eq!(report.boundary_change.species_masses[species], 0.0);
                close(
                    report.after.species_masses[species],
                    initial.species_masses[species],
                    3e-13,
                );
            }
            for side in BoundarySide::ALL {
                let speed = if side.sign() < 0.0 {
                    motion.lower_velocity[side.axis()]
                } else {
                    motion.upper_velocity[side.axis()]
                };
                for impulse in flow.boundary_impulses(side) {
                    assert!(impulse.species_masses.iter().all(|x| *x == 0.0));
                    close(impulse.energy, speed * impulse.momentum[side.axis()], 3e-14);
                    for d in 0..3 {
                        if d != side.axis() {
                            assert_eq!(impulse.momentum[d], 0.0);
                        }
                    }
                }
            }
            assert!(report.boundary_change.energy > 0.0);
            assert!(report.after.momentum[axis] < report.before.momentum[axis]);
            close(
                report.after.absolute_energy(flow.gas()).unwrap()
                    - report.before.absolute_energy(flow.gas()).unwrap(),
                report.boundary_change.energy,
                3e-12,
            );
        }
        assert!(flow.totals().energy > initial.energy);
        assert!(flow.grid().spacing()[axis] < grid.spacing()[axis]);
    }
}

#[test]
fn g1_original_mixture_adiabatic_expansion_has_second_order_pressure() {
    // Both active species have Cp/R=3.5, so gamma=1.4 despite different
    // molar masses and formation references. The exact affine solution is
    // rho_i=rho_i0/(1+a*t), u=a*x0, p=p0/(1+a*t)^gamma.
    let grid = CartesianGrid::new(4, 1, 1, 0.25, 1.0, 1.0).unwrap();
    let strain = 0.4;
    let end = 0.2;
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let mut flow = field(standard(0.0), grid, WALLS, |x| {
            ([0.6, 0.4, 0.0], [strain * x[0], 0.0, 0.0], 1.0)
        });
        let motion = MeshMotion {
            lower_velocity: [0.0; 3],
            upper_velocity: [strain, 0.0, 0.0],
        };
        let dt = end / steps as f64;
        for _ in 0..steps {
            let report = flow
                .step_moving(Time::new(dt), StepControl::default(), motion, || false)
                .unwrap();
            close(report.elapsed.value(), dt, 1e-15);
            balance(report);
        }
        let stretch = 1.0 + strain * end;
        let mut error: f64 = 0.0;
        for k in 0..flow.cell_count() {
            let q = flow.cell(k).unwrap();
            for (rho, original) in q.partial_densities.iter().zip([0.6, 0.4, 0.0]) {
                close(*rho, original / stretch, 5e-13);
            }
            let p = q.primitive(flow.gas(), || false).unwrap();
            close(
                p.velocity[0],
                strain * grid.cell_center(k).unwrap()[0],
                5e-12,
            );
            error = error.max((p.pressure - stretch.powf(-1.4)).abs());
        }
        errors.push(error);
    }
    assert!(errors[0] > 3.5 * errors[1], "{errors:?}");
    assert!(errors[1] > 3.5 * errors[2], "{errors:?}");
    assert!(errors[2] < 2e-7, "{errors:?}");
}

#[test]
fn g1_species_advection_refines_on_a_translating_three_dimensional_grid() {
    let velocity = [0.4, 0.25, -0.15];
    let motion = MeshMotion {
        lower_velocity: [0.11, -0.12, 0.08],
        upper_velocity: [0.11, -0.12, 0.08],
    };
    let end = 0.08;
    let mut errors = Vec::new();
    for n in [8, 16] {
        let dx = 1.0 / n as f64;
        let z = std::f64::consts::PI * dx;
        let average = (z.sin() / z).powi(3);
        let exact = |x: [f64; 3], time: f64| {
            0.55 + 0.2 * average * (std::f64::consts::TAU * (x[0] + x[1] + x[2] - 0.5 * time)).sin()
        };
        let grid = CartesianGrid::new(n, n, n, dx, dx, dx).unwrap();
        let mut flow = field(standard(0.0), grid, PERIODIC, |x| {
            let a = exact(x, 0.0);
            ([a, 1.0 - a, 0.0], velocity, 1.0)
        });
        let initial = flow.totals();
        while flow.time().value() < end {
            balance(
                flow.step_moving(
                    Time::new(end - flow.time().value()),
                    StepControl::default(),
                    motion,
                    || false,
                )
                .unwrap(),
            );
        }
        let mut error = 0.0;
        for k in 0..flow.cell_count() {
            let q = flow.cell(k).unwrap();
            let p = q.primitive(flow.gas(), || false).unwrap();
            close(p.pressure, 1.0, 5e-11);
            close(
                p.temperature,
                1.0 / (q.partial_densities[0] + 0.5 * q.partial_densities[1]),
                5e-11,
            );
            assert_eq!(q.partial_densities[2], 0.0);
            error +=
                (q.partial_densities[0] - exact(flow.grid().cell_center(k).unwrap(), end)).abs();
        }
        for d in 0..7 {
            close(flow.totals().component(d), initial.component(d), 4e-12);
        }
        errors.push(error / flow.cell_count() as f64);
    }
    let order = (errors[0] / errors[1]).ln() / 2.0_f64.ln();
    eprintln!("moving mixture contact L1 on 8^3/16^3: {errors:?}, order={order}");
    assert!(order > 1.3, "{errors:?}");
    assert!(errors[1] < 3e-3, "{errors:?}");
}

#[test]
fn g3_moving_species_momentum_and_boundary_energy_transform_under_a_boost() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let boost = [0.7, -0.4, 0.3];
    let initial = |x: [f64; 3]| {
        (
            [0.6 + 0.1 * x[0], 0.4 - 0.03 * x[1], 0.0],
            [0.2 + 0.1 * x[1], -0.1, 0.3],
            1.0 + 0.2 * x[2],
        )
    };
    let mut flow = field(standard(0.4), grid, OPENS, initial);
    let mut boosted = field(standard(0.4), grid, OPENS, |x| {
        let (species, velocity, pressure) = initial(x);
        (
            species,
            std::array::from_fn(|a| velocity[a] + boost[a]),
            pressure,
        )
    });
    let motion = MeshMotion {
        lower_velocity: [-0.1, 0.05, 0.0],
        upper_velocity: [0.15, 0.05, 0.0],
    };
    let shifted = MeshMotion {
        lower_velocity: std::array::from_fn(|a| motion.lower_velocity[a] + boost[a]),
        upper_velocity: std::array::from_fn(|a| motion.upper_velocity[a] + boost[a]),
    };
    let a = flow
        .step_moving(Time::new(0.004), StepControl::default(), motion, || false)
        .unwrap();
    let b = boosted
        .step_moving(Time::new(0.004), StepControl::default(), shifted, || false)
        .unwrap();
    balance(a);
    balance(b);
    close(a.elapsed.value(), b.elapsed.value(), 1e-15);
    for k in 0..flow.cell_count() {
        let q = flow.cell(k).unwrap();
        let shifted = boosted.cell(k).unwrap();
        for (a, b) in q.partial_densities.iter().zip(shifted.partial_densities) {
            close(*a, *b, 5e-12);
        }
        let mass = sum(q.partial_densities.iter().copied());
        for d in 0..3 {
            close(shifted.momentum[d], q.momentum[d] + boost[d] * mass, 5e-12);
        }
        let work: f64 = (0..3)
            .map(|d| boost[d] * q.momentum[d] + 0.5 * boost[d] * boost[d] * mass)
            .sum();
        close(shifted.energy, q.energy + work, 6e-12);
    }
    for side in BoundarySide::ALL {
        for (a, b) in flow
            .boundary_impulses(side)
            .zip(boosted.boundary_impulses(side))
        {
            for (x, y) in a.species_masses.iter().zip(b.species_masses) {
                close(*x, *y, 4e-13);
            }
            let mass = sum(a.species_masses.iter().copied());
            for d in 0..3 {
                close(b.momentum[d], a.momentum[d] + boost[d] * mass, 4e-13);
            }
            let work: f64 = (0..3)
                .map(|d| boost[d] * a.momentum[d] + 0.5 * boost[d] * boost[d] * mass)
                .sum();
            close(b.energy, a.energy + work, 5e-13);
        }
    }
}

#[test]
fn g3_original_formation_references_do_not_drive_mixture_motion() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let models = |shift: f64| {
        gas(vec![
            source("A", R, 3.5, 0.0, 30.0 + shift),
            source("B", 2.0 * R, 3.5, 0.4, -60.0 - 2.0 * shift),
            source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
        ])
    };
    let sample = |x: [f64; 3]| {
        (
            [0.6 + 0.1 * x[0], 0.4 - 0.1 * x[0], 0.0],
            [0.2, -0.1, 0.3],
            1.0,
        )
    };
    let mut original = field(models(0.0), grid, OPENS, sample);
    let mut shifted = field(models(1024.0), grid, OPENS, sample);
    let motion = MeshMotion {
        lower_velocity: [-0.1, 0.05, 0.0],
        upper_velocity: [0.15, 0.05, 0.0],
    };
    for flow in [&mut original, &mut shifted] {
        let report = flow
            .step_moving(Time::new(0.004), StepControl::default(), motion, || false)
            .unwrap();
        close(report.absolute_energy_residual_j, 0.0, 3e-10);
    }
    for k in 0..original.cell_count() {
        let p = original
            .cell(k)
            .unwrap()
            .primitive(original.gas(), || false)
            .unwrap();
        let q = shifted
            .cell(k)
            .unwrap()
            .primitive(shifted.gas(), || false)
            .unwrap();
        close(p.density, q.density, 3e-11);
        close(p.pressure, q.pressure, 3e-11);
        close(p.temperature, q.temperature, 3e-11);
        for d in 0..3 {
            close(p.velocity[d], q.velocity[d], 3e-11);
        }
    }
    // Physical energy retains the different input reference chart even though
    // the dynamics agree. Species-referenced energy is not a second heat source.
    let expected_shift = sum((0..3).map(|d| {
        original.totals().species_masses[d]
            * (shifted.gas().energy_references()[d] - original.gas().energy_references()[d])
    }));
    close(
        shifted.totals().absolute_energy(shifted.gas()).unwrap()
            - original.totals().absolute_energy(original.gas()).unwrap(),
        expected_shift,
        3e-11,
    );
}

#[test]
fn g4_failed_motion_and_last_poll_preserve_the_accepted_transaction() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let make = || {
        field(standard(0.4), grid, WALLS, |_| {
            ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
        })
    };
    let motion = MeshMotion {
        lower_velocity: [0.0; 3],
        upper_velocity: [-0.1, 0.0, 0.0],
    };
    let mut control = make();
    let mut refused = make();
    for flow in [&mut control, &mut refused] {
        flow.step_moving(Time::new(0.001), StepControl::default(), motion, || false)
            .unwrap();
    }
    let states = refused.states().to_vec();
    let geometry = *refused.grid();
    let time = refused.time();
    let report = *refused.last_step().unwrap();
    let loads = refused.boundary_impulses.clone();
    for invalid in [
        MeshMotion {
            lower_velocity: [f64::NAN, 0.0, 0.0],
            ..motion
        },
        MeshMotion {
            upper_velocity: [-0.1, 0.2, 0.0],
            ..motion
        },
    ] {
        assert!(
            refused
                .step_moving(Time::new(0.01), StepControl::default(), invalid, || false)
                .is_err()
        );
        assert_eq!(refused.states(), states);
        assert_eq!(*refused.grid(), geometry);
        assert_eq!(refused.time(), time);
        assert_eq!(*refused.last_step().unwrap(), report);
        assert_eq!(refused.boundary_impulses, loads);
    }
    let mut polls = 0;
    let expected = control
        .step_moving(Time::new(0.002), StepControl::default(), motion, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert_eq!(
        refused.step_moving(Time::new(0.002), StepControl::default(), motion, || {
            count += 1;
            count == polls
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(refused.states(), states);
    assert_eq!(*refused.grid(), geometry);
    assert_eq!(refused.time(), time);
    assert_eq!(*refused.last_step().unwrap(), report);
    assert_eq!(refused.boundary_impulses, loads);
    let replay = refused
        .step_moving(Time::new(0.002), StepControl::default(), motion, || false)
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(refused.states(), control.states());
    assert_eq!(refused.grid(), control.grid());
}

#[test]
fn g4_unsupported_boundary_data_and_lost_translation_refuse_atomically() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let initial = |_| ([0.6, 0.4, 0.0], [0.0; 3], 1.0);
    let motion = MeshMotion {
        lower_velocity: [0.1, 0.0, 0.0],
        upper_velocity: [0.1, 0.0, 0.0],
    };
    let mut periodic = field(standard(0.4), grid, PERIODIC, initial);
    assert!(
        periodic
            .step_moving(
                Time::new(0.01),
                StepControl::default(),
                MeshMotion {
                    upper_velocity: [0.2, 0.0, 0.0],
                    ..motion
                },
                || false
            )
            .is_err()
    );
    let mut profiled = field(standard(0.4), grid, OPENS, initial);
    let sides = [BoundarySide::XMin];
    let bytes = profiled
        .required_reservoir_profile_workspace_bytes(&sides)
        .unwrap();
    profiled
        .prepare_reservoir_profiles(
            &sides,
            bytes,
            |_, _, _, out| {
                out.copy_from_slice(&[0.6, 0.4, 0.0]);
                Ok(ReservoirFaceState {
                    velocity: [0.0; 3],
                    pressure_pa: 1.0,
                })
            },
            || false,
        )
        .unwrap();
    let before = profiled.states().to_vec();
    assert!(
        profiled
            .step_moving(Time::new(0.01), StepControl::default(), motion, || false)
            .is_err()
    );
    assert_eq!(profiled.states(), before);
    assert_eq!(*profiled.grid(), grid);
    let source = standard(0.4);
    let exterior = ReservoirBoundary {
        side: BoundarySide::XMin,
        partial_densities: &[0.6, 0.4, 0.0],
        velocity: [0.0; 3],
        pressure_pa: 1.0,
    };
    let bytes =
        MixtureEuler3d::required_workspace_bytes_with_reservoirs(&source, &grid, 1).unwrap();
    let mut reservoir = MixtureEuler3d::new_with_reservoirs(
        source,
        grid,
        before,
        OPENS,
        &[exterior],
        bytes,
        || false,
    )
    .unwrap();
    assert!(
        reservoir
            .step_moving(Time::new(0.01), StepControl::default(), motion, || false)
            .is_err()
    );
    assert_eq!(*reservoir.grid(), grid);
    assert_eq!(reservoir.time().value(), 0.0);
    let mut lost = field(standard(0.4), grid, WALLS, initial);
    lost.grid.origin[0] = 1e12;
    let old = *lost.grid();
    let states = lost.states().to_vec();
    let tiny = MeshMotion {
        lower_velocity: [1e-8, 0.0, 0.0],
        upper_velocity: [1e-8, 0.0, 0.0],
    };
    assert!(
        lost.step_moving(
            Time::new(0.001),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            tiny,
            || false
        )
        .is_err()
    );
    assert_eq!(lost.states(), states);
    assert_eq!(*lost.grid(), old);
    assert_eq!(lost.time().value(), 0.0);
    assert!(lost.last_step().is_none());
    assert!(lost.boundary_impulses.iter().all(|x| *x == 0.0));
}
