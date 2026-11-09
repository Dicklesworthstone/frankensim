use super::*;
use crate::compressible::tests::{close, gas};

const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const WALLS: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};

fn controls(velocity: [f64; 2]) -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: velocity.map(|v| [v - 4.0, v + 4.0]),
        momentum_tolerance_n_s: 3e-11,
        energy_tolerance_j: 3e-11,
        maximum_iterations: 24,
    }
}

fn field(grid: CartesianGrid, sample: impl Fn([f64; 3]) -> (f64, [f64; 3], f64)) -> Euler3d {
    let source = gas(0.4);
    let states = (0..grid.len())
        .map(|cell| {
            let (rho, velocity, pressure) = sample(grid.cell_center(cell).unwrap());
            State::from_primitive(&source, rho, velocity, pressure).unwrap()
        })
        .collect();
    Euler3d::new(source, grid, states, WALLS, || false).unwrap()
}

fn small_grid() -> CartesianGrid {
    CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap()
}

fn kinetic(wall: Piston) -> f64 {
    0.5 * wall.mass_kg() * wall.velocity_m_s().powi(2)
}

fn check(report: SpatialInertialWallsStep, tolerance: InertialWallsControl) {
    let r = report.mechanics;
    assert!(r.gas.maximum_cfl <= 0.40000000000001);
    assert_eq!(r.gas.source_change, ForceChange3d::default());
    for residual in r.gas.conservation_residual.array() {
        close(residual, 0.0, 3e-12);
    }
    for residual in r
        .momentum_residual_n_s
        .into_iter()
        .chain([r.combined_momentum_residual_n_s])
    {
        assert!(residual.abs() <= tolerance.momentum_tolerance_n_s);
    }
    for residual in r
        .mechanical_energy_residual_j
        .into_iter()
        .chain([r.total_energy_residual_j])
    {
        assert!(residual.abs() <= tolerance.energy_tolerance_j);
    }
    for end in 0..2 {
        for axis in 0..3 {
            assert_eq!(
                report.guide_angular_impulse_n_m_s[end][axis],
                -report.gas_angular_impulse_n_m_s[end][axis]
            );
        }
    }
}

#[test]
fn g0_comoving_pressure_equilibrium_retains_native_cells_and_both_finite_walls() {
    for axis in 0..3 {
        for velocity in [-0.25, 0.0, 0.35] {
            let grid = small_grid();
            let mut motion = [0.0; 3];
            motion[axis] = velocity;
            let mut flow = field(grid, |_| (1.0, motion, 1.0));
            let initial = flow.states().to_vec();
            let mut walls = [
                Piston::new(0.2, velocity).unwrap(),
                Piston::new(3.0, velocity).unwrap(),
            ];
            let original = walls;
            let control = controls([velocity; 2]);
            let report = flow
                .step_with_inertial_walls(
                    axis,
                    &mut walls,
                    Time::new(0.005),
                    StepControl::default(),
                    [InertialWallLoad {
                        exterior_pressure_pa: 1.0,
                        external_force_n: 0.0,
                    }; 2],
                    control,
                    || false,
                )
                .unwrap();
            check(report, control);
            assert_eq!(walls, original);
            assert_eq!(report.mechanics.iterations, 0);
            assert_eq!(report.mechanics.gas_trials, 1);
            close(
                flow.grid().origin()[axis],
                velocity * report.mechanics.gas.elapsed.value(),
                3e-15,
            );
            assert_eq!(flow.grid().spacing(), grid.spacing());
            for (q, before) in flow.states().iter().zip(&initial) {
                for (actual, expected) in q.array().into_iter().zip(before.array()) {
                    close(actual, expected, 2e-14);
                }
            }
        }
    }
}

#[test]
fn g1_initial_plane_acceleration_matches_pressure_area_over_each_actual_mass() {
    for axis in 0..3 {
        let grid = small_grid();
        let mut flow = field(grid, |_| (1.0, [0.0; 3], 2.0));
        let mut walls = [
            Piston::new(2.0, 0.0).unwrap(),
            Piston::new(5.0, 0.0).unwrap(),
        ];
        let control = controls([0.0; 2]);
        let report = flow
            .step_with_inertial_walls(
                axis,
                &mut walls,
                Time::new(1e-6),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                control,
                || false,
            )
            .unwrap();
        check(report, control);
        let side = BoundarySide::ALL[2 * axis];
        let force = 2.0 * grid.face_area(side) * grid.face_count(side) as f64;
        let dt = report.mechanics.gas.elapsed.value();
        close(
            walls[0].velocity_m_s() / dt,
            -force / walls[0].mass_kg(),
            4e-6,
        );
        close(
            walls[1].velocity_m_s() / dt,
            force / walls[1].mass_kg(),
            4e-6,
        );
        assert!(report.mechanics.midpoint_velocity_m_s[0] < 0.0);
        assert!(report.mechanics.midpoint_velocity_m_s[1] > 0.0);
    }
}

#[test]
fn g1_resolved_unequal_mass_expansion_closes_gas_plus_wall_momentum_and_energy() {
    for axis in 0..3 {
        let mut flow = field(small_grid(), |_| (1.0, [0.0; 3], 2.0));
        let mut walls = [
            Piston::new(0.2, 0.0).unwrap(),
            Piston::new(1.5, 0.0).unwrap(),
        ];
        let initial = flow.totals();
        let mut previous = initial;
        for _ in 0..12 {
            let before = walls;
            let control = controls(walls.map(|wall| wall.velocity_m_s()));
            let report = flow
                .step_with_inertial_walls(
                    axis,
                    &mut walls,
                    Time::new(0.01),
                    StepControl::default(),
                    [InertialWallLoad::default(); 2],
                    control,
                    || false,
                )
                .unwrap();
            check(report, control);
            assert_eq!(flow.last_step(), Some(&report.mechanics.gas));
            let after = flow.totals();
            close(after.mass, initial.mass, 2e-13);
            close(
                after.momentum[axis]
                    + walls
                        .iter()
                        .map(|wall| wall.mass_kg() * wall.velocity_m_s())
                        .sum::<f64>(),
                0.0,
                4e-10,
            );
            close(
                after.energy + walls.into_iter().map(kinetic).sum::<f64>(),
                initial.energy,
                4e-10,
            );
            assert!(after.energy < previous.energy);
            for end in 0..2 {
                close(
                    kinetic(walls[end]) - kinetic(before[end]),
                    report.mechanics.kinetic_energy_change_j[end],
                    1e-14,
                );
                let side = BoundarySide::ALL[2 * axis + end];
                let impulse: f64 = flow
                    .boundary_impulses(side)
                    .iter()
                    .map(|q| q.momentum[axis])
                    .sum();
                let work: f64 = flow.boundary_impulses(side).iter().map(|q| q.energy).sum();
                close(impulse, report.mechanics.gas_impulse_n_s[end], 2e-14);
                close(
                    work,
                    report.mechanics.midpoint_velocity_m_s[end] * impulse,
                    3e-14,
                );
                close(kinetic(walls[end]) - kinetic(before[end]), work, 4e-11);
            }
            previous = after;
        }
        assert!(flow.grid().spacing()[axis] > small_grid().spacing()[axis]);
        assert!(walls[0].velocity_m_s() < 0.0 && walls[1].velocity_m_s() > 0.0);
    }
}

#[test]
fn g3_native_motion_replay_and_nonuniform_pressure_torque_keep_original_face_work() {
    for axis in 0..3 {
        let grid = small_grid();
        let [a, b] = tangents(axis);
        let sample = |p: [f64; 3]| (1.0 + 0.05 * p[b], [0.0; 3], 1.0 + 0.8 * p[a]);
        let mut flow = field(grid, sample);
        let mut replay = field(grid, sample);
        let mut walls = [
            Piston::new(0.3, 0.0).unwrap(),
            Piston::new(2.0, 0.0).unwrap(),
        ];
        let control = controls([0.0; 2]);
        let report = flow
            .step_with_inertial_walls(
                axis,
                &mut walls,
                Time::new(0.01),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                control,
                || false,
            )
            .unwrap();
        check(report, control);
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = report.mechanics.midpoint_velocity_m_s[0];
        motion.upper_velocity[axis] = report.mechanics.midpoint_velocity_m_s[1];
        let prescribed = replay
            .step_moving(
                report.mechanics.gas.elapsed,
                StepControl::default(),
                motion,
                || false,
            )
            .unwrap();
        assert_eq!(prescribed, report.mechanics.gas);
        assert_eq!(flow.states(), replay.states());
        assert_eq!(flow.grid(), replay.grid());
        assert_eq!(flow.time(), replay.time());
        for side in BoundarySide::ALL {
            assert_eq!(flow.boundary_impulses(side), replay.boundary_impulses(side));
        }
        for end in 0..2 {
            let side = BoundarySide::ALL[2 * axis + end];
            let mut torque = [0.0; 3];
            for (face, impulse) in flow.boundary_impulses(side).iter().enumerate() {
                let position = grid.face_center(side, face).unwrap();
                let mut arm = [0.0; 3];
                arm[a] = position[a] - 0.5 * grid.shape()[a] as f64 * grid.spacing()[a];
                arm[b] = position[b] - 0.5 * grid.shape()[b] as f64 * grid.spacing()[b];
                let j = impulse.momentum;
                let cross = [
                    arm[1] * j[2] - arm[2] * j[1],
                    arm[2] * j[0] - arm[0] * j[2],
                    arm[0] * j[1] - arm[1] * j[0],
                ];
                for d in 0..3 {
                    torque[d] += cross[d];
                }
            }
            assert!(torque.into_iter().any(|x| x.abs() > 1e-5));
            for d in 0..3 {
                close(report.gas_angular_impulse_n_m_s[end][d], torque[d], 3e-15);
            }
        }
    }
}

#[test]
fn g3_axial_boost_preserves_reciprocal_motion_and_transforms_total_work() {
    let boost = 1.75;
    for axis in 0..3 {
        let grid = small_grid();
        let length = grid.shape()[axis] as f64 * grid.spacing()[axis];
        let sample = |x: [f64; 3], shift: f64| {
            let s = x[axis] / length;
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.1 * (2.0 * core::f64::consts::PI * s).sin() + shift;
            (1.0 + 0.1 * s, velocity, 1.0 + 0.2 * s)
        };
        let mut base = field(grid, |x| sample(x, 0.0));
        let mut shifted = field(grid, |x| sample(x, boost));
        let mut a = [
            Piston::new(0.5, -0.1).unwrap(),
            Piston::new(2.0, 0.2).unwrap(),
        ];
        let mut b = a.map(|wall| Piston::new(wall.mass_kg(), wall.velocity_m_s() + boost).unwrap());
        let loads = [
            InertialWallLoad {
                exterior_pressure_pa: 0.8,
                external_force_n: 0.1,
            },
            InertialWallLoad {
                exterior_pressure_pa: 1.2,
                external_force_n: -0.2,
            },
        ];
        let ca = controls(a.map(|wall| wall.velocity_m_s()));
        let cb = controls(b.map(|wall| wall.velocity_m_s()));
        let ra = base
            .step_with_inertial_walls(
                axis,
                &mut a,
                Time::new(0.005),
                StepControl::default(),
                loads,
                ca,
                || false,
            )
            .unwrap();
        let rb = shifted
            .step_with_inertial_walls(
                axis,
                &mut b,
                Time::new(0.005),
                StepControl::default(),
                loads,
                cb,
                || false,
            )
            .unwrap();
        check(ra, ca);
        check(rb, cb);
        assert_eq!(ra.mechanics.gas.elapsed, rb.mechanics.gas.elapsed);
        let dt = ra.mechanics.gas.elapsed.value();
        for end in 0..2 {
            close(b[end].velocity_m_s(), a[end].velocity_m_s() + boost, 3e-10);
            close(
                rb.mechanics.gas_impulse_n_s[end],
                ra.mechanics.gas_impulse_n_s[end],
                3e-10,
            );
        }
        for (p, q) in base.states().iter().zip(shifted.states()) {
            close(q.density, p.density, 2e-10);
            for d in 0..3 {
                close(
                    q.momentum[d],
                    p.momentum[d] + if d == axis { p.density * boost } else { 0.0 },
                    2e-10,
                );
            }
            close(
                q.energy,
                p.energy + boost * p.momentum[axis] + 0.5 * boost * boost * p.density,
                2e-10,
            );
        }
        for d in 0..3 {
            close(shifted.grid().spacing()[d], base.grid().spacing()[d], 2e-12);
            close(
                shifted.grid().origin()[d],
                base.grid().origin()[d] + if d == axis { boost * dt } else { 0.0 },
                2e-12,
            );
        }
        let ga = ra.mechanics.gas;
        let gb = rb.mechanics.gas;
        close(
            gb.boundary_change.energy,
            ga.boundary_change.energy + boost * ga.boundary_change.momentum[axis],
            3e-10,
        );
    }
}

#[test]
fn g4_complete_coupled_retry_uses_one_duration_and_exhaustion_preserves_all_owners() {
    let mut flow = field(small_grid(), |_| (1.0, [0.0; 3], 1.0));
    let mut walls = [Piston::new(1.0, 0.0).unwrap(); 2];
    let original = walls;
    let states = flow.states().to_vec();
    let grid = *flow.grid();
    let loads = [
        InertialWallLoad {
            exterior_pressure_pa: 0.0,
            external_force_n: -1000.0,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.0,
            external_force_n: 1000.0,
        },
    ];
    let control = InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[-0.1, 0.1]; 2],
        ..controls([0.0; 2])
    };
    let failed = flow.step_with_inertial_walls(
        0,
        &mut walls,
        Time::new(0.005),
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        loads,
        control,
        || false,
    );
    assert!(matches!(
        failed,
        Err(InertialWallsError::StepRetries { attempts: 1, .. })
    ));
    assert_eq!(walls, original);
    assert_eq!(flow.states(), states);
    assert_eq!(*flow.grid(), grid);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(flow.last_step(), None);
    let report = flow
        .step_with_inertial_walls(
            0,
            &mut walls,
            Time::new(0.005),
            StepControl::default(),
            loads,
            control,
            || false,
        )
        .unwrap();
    check(report, control);
    assert!(report.mechanics.gas.halvings >= 4);
    assert!(report.mechanics.gas_trials > report.mechanics.iterations + 1);
    let dt = report.mechanics.gas.elapsed.value();
    assert_eq!(flow.time().value(), dt);
    for end in 0..2 {
        close(
            report.mechanics.displacement_m[end],
            dt * report.mechanics.midpoint_velocity_m_s[end],
            3e-15,
        );
        close(
            kinetic(walls[end]),
            report.mechanics.external_force_work_j[end]
                + report.mechanics.gas_impulse_n_s[end]
                    * report.mechanics.midpoint_velocity_m_s[end],
            4e-11,
        );
    }
}

#[test]
fn g4_final_checkpoint_cancellation_preserves_previous_field_loads_time_and_walls() {
    let setup = || {
        let mut flow = field(small_grid(), |_| (1.0, [0.0; 3], 1.0));
        let mut walls = [
            Piston::new(0.3, 0.0).unwrap(),
            Piston::new(2.0, 0.0).unwrap(),
        ];
        flow.step_with_inertial_walls(
            0,
            &mut walls,
            Time::new(0.005),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls([0.0; 2]),
            || false,
        )
        .unwrap();
        (flow, walls)
    };
    let (mut reference, mut reference_walls) = setup();
    let control = controls(reference_walls.map(|wall| wall.velocity_m_s()));
    let mut polls = 0;
    let expected = reference
        .step_with_inertial_walls(
            0,
            &mut reference_walls,
            Time::new(0.005),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            control,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, mut walls) = setup();
    let saved_states = flow.states().to_vec();
    let saved_grid = *flow.grid();
    let saved_time = flow.time();
    let saved_report = *flow.last_step().unwrap();
    let saved_walls = walls;
    let saved_loads = BoundarySide::ALL.map(|side| flow.boundary_impulses(side).to_vec());
    let mut seen = 0;
    let refusal = flow.step_with_inertial_walls(
        0,
        &mut walls,
        Time::new(0.005),
        StepControl::default(),
        [InertialWallLoad::default(); 2],
        control,
        || {
            seen += 1;
            seen == polls
        },
    );
    assert_eq!(refusal, Err(InertialWallsError::Cancelled));
    assert_eq!(seen, polls);
    assert_eq!(flow.states(), saved_states);
    assert_eq!(*flow.grid(), saved_grid);
    assert_eq!(flow.time(), saved_time);
    assert_eq!(flow.last_step(), Some(&saved_report));
    assert_eq!(walls, saved_walls);
    for (side, loads) in BoundarySide::ALL.into_iter().zip(saved_loads) {
        assert_eq!(flow.boundary_impulses(side), loads);
    }
    let actual = flow
        .step_with_inertial_walls(
            0,
            &mut walls,
            Time::new(0.005),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            control,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(walls, reference_walls);
}

#[test]
fn g4_unsupported_open_geometry_and_unrepresentable_plane_travel_are_atomic() {
    let mut flow = field(small_grid(), |_| (1.0, [0.0; 3], 1.0));
    flow.boundaries.y = AxisBoundaries::Ends {
        lower: BoundaryCondition::Transmissive,
        upper: BoundaryCondition::Reflecting,
    };
    let mut walls = [Piston::new(1.0, 0.0).unwrap(); 2];
    let cells = flow.states().to_vec();
    let refusal = flow.step_with_inertial_walls(
        0,
        &mut walls,
        Time::new(0.005),
        StepControl::default(),
        [InertialWallLoad::default(); 2],
        controls([0.0; 2]),
        || false,
    );
    assert!(refusal.is_err());
    assert_eq!(flow.states(), cells);
    assert_eq!(flow.time().value(), 0.0);
    // A stored upper endpoint near 1e16 cannot express this microscopic move,
    // even though the minimum origin itself can represent its translation.
    let grid = CartesianGrid::new(1, 1, 1, 1e16, 1.0, 1.0).unwrap();
    let mut enormous = field(grid, |_| (1.0, [0.0; 3], 1.0));
    let mut translating = [Piston::new(1.0, 1.0).unwrap(); 2];
    let initial_walls = translating;
    let cells = enormous.states().to_vec();
    let refusal = enormous.step_with_inertial_walls(
        0,
        &mut translating,
        Time::new(1e-6),
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        [InertialWallLoad {
            exterior_pressure_pa: 1.0,
            external_force_n: 0.0,
        }; 2],
        controls([1.0; 2]),
        || false,
    );
    assert!(refusal.is_err());
    assert_eq!(enormous.states(), cells);
    assert_eq!(*enormous.grid(), grid);
    assert_eq!(enormous.time().value(), 0.0);
    assert_eq!(translating, initial_walls);
}
