use super::*;
use crate::compressible::tests::{close, gas};

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
    source: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    sample: impl Fn([f64; 3]) -> (f64, [f64; 3], f64),
) -> Euler3d {
    let states = (0..grid.len())
        .map(|cell| {
            let (density, velocity, pressure) = sample(grid.cell_center(cell).unwrap());
            State::from_primitive(&source, density, velocity, pressure).unwrap()
        })
        .collect();
    Euler3d::new(source, grid, states, boundaries, || false).unwrap()
}

fn balance(report: StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    assert_eq!(report.source_change, ForceChange3d::default());
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 8e-12);
    }
}

#[test]
fn g0_zero_mesh_motion_is_the_original_stationary_transaction() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let sample = |p: [f64; 3]| (1.0 + 0.1 * p[0], [0.2, -0.1, 0.3], 1.0 + 0.2 * p[1]);
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut stationary = field(gas(0.4), grid, PERIODIC, sample);
        let mut moving = field(gas(0.4), grid, PERIODIC, sample);
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        let a = stationary
            .step(Time::new(0.01), controls, || false)
            .unwrap();
        let b = moving
            .step_moving(Time::new(0.01), controls, MeshMotion::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(stationary.states(), moving.states());
        assert_eq!(stationary.grid(), moving.grid());
        assert_eq!(stationary.time(), moving.time());
        for side in BoundarySide::ALL {
            assert_eq!(
                stationary.boundary_impulses(side),
                moving.boundary_impulses(side)
            );
        }
    }
}

#[test]
fn g1_deforming_three_dimensional_free_stream_obeys_discrete_gcl() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let velocity = [0.3, -0.2, 0.1];
    for axis in 0..3 {
        for dilation in [-0.4, 0.4] {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let mut flow = field(gas(0.4), grid, OPENS, |_| (1.1, velocity, 1.7));
                let state = flow.states()[0];
                let lower_velocity = [-0.1, 0.15, 0.05];
                let mut upper_velocity = lower_velocity;
                upper_velocity[axis] += dilation;
                let motion = MeshMotion {
                    lower_velocity,
                    upper_velocity,
                };
                for _ in 0..3 {
                    let old = *flow.grid();
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
                    let swept =
                        dt * dilation * old.spacing[a] * old.spacing[b] / old.shape[axis] as f64;
                    close(flow.grid().cell_volume() - old.cell_volume(), swept, 3e-15);
                    for d in 0..3 {
                        close(
                            flow.grid().origin()[d],
                            old.origin()[d] + dt * lower_velocity[d],
                            2e-15,
                        );
                    }
                    for cell in flow.states() {
                        for (actual, expected) in cell.array().into_iter().zip(state.array()) {
                            close(actual, expected, 2e-12);
                        }
                    }
                    let change =
                        grid.len() as f64 * (flow.grid().cell_volume() - old.cell_volume());
                    for (actual, expected) in report
                        .boundary_change
                        .array()
                        .into_iter()
                        .zip(state.array())
                    {
                        close(actual, change * expected, 3e-12);
                    }
                }
            }
        }
    }
}

#[test]
fn g1_translating_periodic_grid_retains_comoving_shear_and_density_contacts() {
    let grid = CartesianGrid::new(8, 2, 2, 0.125, 0.3, 0.4).unwrap();
    let speed = [0.3, -0.2, 0.1];
    let motion = MeshMotion {
        lower_velocity: speed,
        upper_velocity: speed,
    };
    let mut flow = field(gas(0.4), grid, PERIODIC, |p| {
        let wave = (std::f64::consts::TAU * p[0]).sin();
        (
            1.0 + 0.2 * wave,
            [speed[0], speed[1] + 0.1 * wave, speed[2] - 0.15 * wave],
            1.0,
        )
    });
    let original = flow.states().to_vec();
    let totals = flow.totals();
    for _ in 0..4 {
        balance(
            flow.step_moving(Time::new(0.01), StepControl::default(), motion, || false)
                .unwrap(),
        );
    }
    for (actual, expected) in flow.states().iter().zip(original) {
        for (a, e) in actual.array().into_iter().zip(expected.array()) {
            close(a, e, 3e-12);
        }
    }
    for (actual, expected) in flow.totals().array().into_iter().zip(totals.array()) {
        close(actual, expected, 5e-12);
    }
    for d in 0..3 {
        close(
            flow.grid().origin()[d],
            speed[d] * flow.time().value(),
            3e-15,
        );
    }
}

#[test]
fn g1_closed_piston_accounts_for_actual_wall_impulse_and_pressure_work() {
    let grid = CartesianGrid::new(6, 2, 2, 0.1, 0.3, 0.4).unwrap();
    for axis in 0..3 {
        let mut flow = field(gas(0.4), grid, WALLS, |_| (1.0, [0.0; 3], 1.0));
        let mut motion = MeshMotion::default();
        motion.upper_velocity[axis] = -0.1;
        let mass = flow.totals().mass;
        let energy = flow.totals().energy;
        for _ in 0..4 {
            let report = flow
                .step_moving(Time::new(0.002), StepControl::default(), motion, || false)
                .unwrap();
            balance(report);
            close(report.boundary_change.mass, 0.0, 2e-15);
            close(flow.totals().mass, mass, 3e-13);
            for side in BoundarySide::ALL {
                let speed = if side.sign() < 0.0 {
                    motion.lower_velocity[side.axis()]
                } else {
                    motion.upper_velocity[side.axis()]
                };
                for impulse in flow.boundary_impulses(side) {
                    assert_eq!(impulse.mass, 0.0);
                    close(impulse.energy, speed * impulse.momentum[side.axis()], 2e-14);
                    for d in 0..3 {
                        if d != side.axis() {
                            assert_eq!(impulse.momentum[d], 0.0);
                        }
                    }
                }
            }
            assert!(report.boundary_change.energy > 0.0);
            assert!(report.after.momentum[axis] < report.before.momentum[axis]);
        }
        assert!(flow.totals().energy > energy);
        assert!(flow.grid().spacing()[axis] < grid.spacing()[axis]);
    }
}

#[test]
fn g1_piston_rk2_refines_toward_the_exact_homogeneous_adiabatic_expansion() {
    // Euler's exact affine expansion has x=x0*(1+a*t), u=a*x0,
    // rho=rho0/(1+a*t), p=p0/(1+a*t)^gamma. The fixture's original
    // NASA9 model has R=1 and constant Cp/R=3.5, hence gamma=1.4.
    // Linear velocity is exactly reconstructed, including moving wall ghosts.
    let grid = CartesianGrid::new(4, 1, 1, 0.25, 1.0, 1.0).unwrap();
    let strain = 0.4;
    let end = 0.2;
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let mut flow = field(gas(0.0), grid, WALLS, |p| {
            (1.0, [strain * p[0], 0.0, 0.0], 1.0)
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
        let pressure = stretch.powf(-1.4);
        let mut error: f64 = 0.0;
        for (cell, q) in flow.states().iter().enumerate() {
            let p = q.primitive(flow.gas(), || false).unwrap();
            close(p.density, 1.0 / stretch, 5e-13);
            close(
                p.velocity[0],
                strain * grid.cell_center(cell).unwrap()[0],
                5e-12,
            );
            error = error.max((p.pressure - pressure).abs());
        }
        errors.push(error);
    }
    assert!(errors[0] > 3.5 * errors[1], "errors={errors:?}");
    assert!(errors[1] > 3.5 * errors[2], "errors={errors:?}");
    assert!(errors[2] < 2e-7, "errors={errors:?}");
}

#[test]
fn g3_moving_volume_step_is_galilean_covariant_in_all_momentum_and_work() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let boost = [0.7, -0.4, 0.3];
    let initial = |p: [f64; 3]| {
        (
            1.0 + 0.1 * p[0],
            [0.2 + 0.1 * p[1], -0.1, 0.3],
            1.0 + 0.2 * p[2],
        )
    };
    let mut flow = field(gas(0.4), grid, OPENS, initial);
    let mut boosted = field(gas(0.4), grid, OPENS, |p| {
        let (rho, u, pressure) = initial(p);
        (rho, std::array::from_fn(|a| u[a] + boost[a]), pressure)
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
    for (q, shifted) in flow.states().iter().zip(boosted.states()) {
        close(shifted.density, q.density, 5e-12);
        for axis in 0..3 {
            close(
                shifted.momentum[axis],
                q.momentum[axis] + boost[axis] * q.density,
                5e-12,
            );
        }
        let work: f64 = (0..3)
            .map(|axis| {
                boost[axis] * q.momentum[axis] + 0.5 * boost[axis] * boost[axis] * q.density
            })
            .sum();
        close(shifted.energy, q.energy + work, 6e-12);
    }
    for side in BoundarySide::ALL {
        for (f, shifted) in flow
            .boundary_impulses(side)
            .iter()
            .zip(boosted.boundary_impulses(side))
        {
            close(shifted.mass, f.mass, 4e-13);
            for axis in 0..3 {
                close(
                    shifted.momentum[axis],
                    f.momentum[axis] + boost[axis] * f.mass,
                    4e-13,
                );
            }
            let work: f64 = (0..3)
                .map(|axis| {
                    boost[axis] * f.momentum[axis] + 0.5 * boost[axis] * boost[axis] * f.mass
                })
                .sum();
            close(shifted.energy, f.energy + work, 5e-13);
        }
    }
}

#[test]
fn g4_moving_geometry_refusal_and_final_cancellation_are_atomic() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let make = || field(gas(0.4), grid, WALLS, |_| (1.0, [0.0; 3], 1.0));
    let motion = MeshMotion {
        lower_velocity: [0.0; 3],
        upper_velocity: [-0.1, 0.0, 0.0],
    };
    assert!(motion.grid_after(&grid, 7.0).is_err());
    let mut control = make();
    let mut refused = make();
    for flow in [&mut control, &mut refused] {
        flow.step_moving(Time::new(0.001), StepControl::default(), motion, || false)
            .unwrap();
    }
    let initial = refused.states().to_vec();
    let geometry = *refused.grid();
    let time = refused.time();
    let last = *refused.last_step().unwrap();
    let loads: Vec<_> = BoundarySide::ALL
        .into_iter()
        .map(|side| refused.boundary_impulses(side).to_vec())
        .collect();
    for invalid in [
        MeshMotion {
            lower_velocity: [f64::NAN, 0.0, 0.0],
            ..motion
        },
        MeshMotion {
            upper_velocity: [-0.1, 0.2, 0.0],
            ..motion
        },
        MeshMotion {
            upper_velocity: [-1e300, 0.0, 0.0],
            ..motion
        },
    ] {
        assert!(
            refused
                .step_moving(
                    Time::new(0.01),
                    StepControl {
                        maximum_halvings: 0,
                        ..StepControl::default()
                    },
                    invalid,
                    || false
                )
                .is_err()
        );
        assert_eq!(refused.states(), initial);
        assert_eq!(*refused.grid(), geometry);
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
    assert_eq!(refused.states(), initial);
    assert_eq!(*refused.grid(), geometry);
    assert_eq!(refused.time(), time);
    assert_eq!(*refused.last_step().unwrap(), last);
    for (side, previous) in BoundarySide::ALL.into_iter().zip(loads) {
        assert_eq!(refused.boundary_impulses(side), previous);
    }
    let replay = refused
        .step_moving(Time::new(0.002), StepControl::default(), motion, || false)
        .unwrap();
    assert_eq!(expected, replay);
    assert_eq!(control.states(), refused.states());
    assert_eq!(control.grid(), refused.grid());
    let mut periodic = field(gas(0.4), grid, PERIODIC, |_| (1.0, [0.0; 3], 1.0));
    assert!(
        periodic
            .step_moving(Time::new(0.01), StepControl::default(), motion, || false)
            .is_err()
    );
}

#[test]
fn g4_unrepresentable_translation_preserves_field_and_geometry() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let mut flow = field(gas(0.4), grid, WALLS, |_| (1.0, [0.0; 3], 1.0));
    // Valid cell widths at this origin, but the requested sub-ulp translation
    // would leave stored geometry fixed while the ALE flux uses moving faces.
    flow.grid.origin[0] = 1e12;
    let before = flow.states().to_vec();
    let geometry = *flow.grid();
    let motion = MeshMotion {
        lower_velocity: [1e-8, 0.0, 0.0],
        upper_velocity: [1e-8, 0.0, 0.0],
    };
    let result = flow.step_moving(
        Time::new(0.001),
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        motion,
        || false,
    );
    assert!(result.is_err());
    assert_eq!(flow.states(), before);
    assert_eq!(*flow.grid(), geometry);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    for side in BoundarySide::ALL {
        assert!(
            flow.boundary_impulses(side)
                .iter()
                .all(|x| *x == Totals::default())
        );
    }
}
