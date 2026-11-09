use super::*;
use crate::compressible::spatial::molecular::{InertialWallLoad, InertialWallsControl, Piston};

fn channel(axis: usize) -> Boundaries3d {
    let mut result = PERIODIC;
    match axis {
        0 => result.x = WALL,
        1 => result.y = WALL,
        2 => result.z = WALL,
        _ => unreachable!(),
    }
    result
}

fn owner(boundaries: Boundaries3d, velocity: [f64; 3]) -> MolecularEuler3d<ConstantTransport> {
    let flow = field(
        gas(0.0),
        CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap(),
        boundaries,
        |_| [1.0, velocity[0], velocity[1], velocity[2], 2.0],
    );
    let material = law(flow.gas(), 0.02, 0.009, 0.04);
    MolecularEuler3d::new(flow, material, || false).unwrap()
}

fn wall(kind: usize, velocity: [f64; 3]) -> MolecularWall3d {
    match kind {
        0 => MolecularWall3d::SlipAdiabatic,
        1 => MolecularWall3d::NoSlip {
            velocity,
            thermal: WallThermalCondition::Adiabatic,
        },
        _ => MolecularWall3d::NavierSlip {
            velocity,
            thermal: WallThermalCondition::Adiabatic,
            slip_length_m: 0.12,
            fluid_heat_fraction: 0.35,
        },
    }
}

fn tensors(view: &StressField3d<'_>) -> Vec<[f64; 6]> {
    (0..view.grid().len())
        .map(|cell| view.cell_stress_pa(cell).unwrap())
        .collect()
}

#[test]
fn g1_moving_observer_returns_independent_normal_and_slip_tractions_on_all_sides() {
    let velocity = [0.15, -0.1, 0.07];
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let mut wall_velocity = [0.0; 3];
        wall_velocity[a] = 0.4;
        wall_velocity[b] = -0.3;
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = 0.23;
        motion.upper_velocity[axis] = 0.23;
        for kind in 0..3 {
            let mut flow = owner(channel(axis), velocity)
                .with_wall(side, wall(kind, wall_velocity))
                .unwrap();
            let momentum = flow.flow().states()[0].momentum;
            let distance = 0.5 * flow.flow().grid().spacing[axis];
            let mut face = [1; 3];
            face[axis] = if side.sign() < 0.0 { 0 } else { 3 };
            let moving = {
                let view = flow.stress_field_moving(motion, || false).unwrap();
                assert_eq!(view.momentum(0), Some(momentum));
                assert_eq!(view.time_s(), 0.0);
                view.face_traction_pa(axis, face).unwrap()
            };
            close(
                moving[axis],
                side.sign() * (4.0 * 0.02 / 3.0 + 0.009) * (0.23 - velocity[axis]) / distance,
                3e-14,
            );
            for t in [a, b] {
                let expected = if kind == 0 {
                    0.0
                } else {
                    side.sign() * 0.02 * (wall_velocity[t] - velocity[t])
                        / (distance + if kind == 1 { 0.0 } else { 0.12 })
                };
                close(moving[t], expected, 3e-14);
            }
            // Explicit motion never changes the stationary observer's wall law.
            let stationary = flow
                .stress_field(|| false)
                .unwrap()
                .face_traction_pa(axis, face)
                .unwrap();
            close(
                stationary[axis],
                -side.sign() * (4.0 * 0.02 / 3.0 + 0.009) * velocity[axis] / distance,
                3e-14,
            );
            let zero = flow
                .stress_field_moving(MeshMotion::default(), || false)
                .unwrap()
                .face_traction_pa(axis, face)
                .unwrap();
            assert_eq!(zero, stationary);
        }
    }
}

#[test]
fn g1_moving_observer_sees_zero_stress_in_a_comoving_box_for_all_slip_laws() {
    let velocity = [0.3, -0.2, 0.15];
    let motion = MeshMotion {
        lower_velocity: velocity,
        upper_velocity: velocity,
    };
    for kind in 0..3 {
        let mut flow = owner(WALLS, velocity);
        for side in BoundarySide::ALL {
            let mut tangent_velocity = velocity;
            tangent_velocity[side.axis()] = 0.0;
            flow = flow.with_wall(side, wall(kind, tangent_velocity)).unwrap();
        }
        flow.step_moving(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
        let time = flow.flow().time().value();
        let grid = *flow.flow().grid();
        let view = flow.stress_field_moving(motion, || false).unwrap();
        assert_eq!(view.time_s(), time);
        assert_eq!(*view.grid(), grid);
        for tensor in tensors(&view) {
            for component in tensor {
                close(component, 0.0, 3e-14);
            }
        }
    }
}

#[test]
fn g1_moving_observer_uses_explicit_piston_endpoint_velocity_not_interval_midpoint() {
    for axis in 0..3 {
        let mut flow = owner(channel(axis), [0.0; 3]);
        let mut pistons = [
            Piston::new(0.8, 0.0).unwrap(),
            Piston::new(1.7, 0.0).unwrap(),
        ];
        let report = flow
            .step_with_inertial_walls(
                axis,
                &mut pistons,
                Time::new(0.001),
                MolecularStepControl::default(),
                [InertialWallLoad::default(); 2],
                InertialWallsControl {
                    midpoint_velocity_bracket_m_s: [[-1.0, 1.0]; 2],
                    momentum_tolerance_n_s: 3e-11,
                    energy_tolerance_j: 3e-11,
                    maximum_iterations: 24,
                },
                || false,
            )
            .unwrap();
        let mut endpoint = MeshMotion::default();
        endpoint.lower_velocity[axis] = pistons[0].velocity_m_s();
        endpoint.upper_velocity[axis] = pistons[1].velocity_m_s();
        let mut midpoint = MeshMotion::default();
        midpoint.lower_velocity[axis] = report.mechanics.midpoint_velocity_m_s[0];
        midpoint.upper_velocity[axis] = report.mechanics.midpoint_velocity_m_s[1];
        let grid = *flow.flow().grid();
        let half_width = 0.5 * grid.spacing[axis];
        let observed = {
            let view = flow.stress_field_moving(endpoint, || false).unwrap();
            std::array::from_fn::<_, 2, _>(|end| {
                let mut face = [1; 3];
                face[axis] = 3 * end;
                let mut cell = face;
                cell[axis] = 2 * end;
                let p = view.primitive(grid.index(cell)).unwrap();
                let speed = pistons[end].velocity_m_s();
                let sign = if end == 0 { -1.0 } else { 1.0 };
                let actual = view.face_traction_pa(axis, face).unwrap()[axis];
                close(
                    actual,
                    sign * (4.0 * 0.02 / 3.0 + 0.009) * (speed - p.velocity[axis]) / half_width,
                    4e-14,
                );
                actual
            })
        };
        let held = flow.stress_field_moving(midpoint, || false).unwrap();
        for end in 0..2 {
            let mut face = [1; 3];
            face[axis] = 3 * end;
            let sign = if end == 0 { -1.0 } else { 1.0 };
            let difference = observed[end] - held.face_traction_pa(axis, face).unwrap()[axis];
            close(
                difference,
                sign * (4.0 * 0.02 / 3.0 + 0.009)
                    * (pistons[end].velocity_m_s() - report.mechanics.midpoint_velocity_m_s[end])
                    / half_width,
                4e-14,
            );
            assert!(difference.abs() > 1e-6);
        }
    }
}

#[test]
fn g4_moving_observation_cancellation_and_refusal_preserve_accepted_history() {
    let motion = MeshMotion {
        lower_velocity: [-0.05, 0.0, 0.0],
        upper_velocity: [0.1, 0.0, 0.0],
    };
    let mut flow = owner(channel(0), [0.1, -0.07, 0.04]);
    let mut reference = owner(channel(0), [0.1, -0.07, 0.04]);
    for owner in [&mut flow, &mut reference] {
        owner
            .step_moving(
                Time::new(0.001),
                MolecularStepControl::default(),
                motion,
                || false,
            )
            .unwrap();
    }
    let states = flow.flow().states().to_vec();
    let grid = *flow.flow().grid();
    let time = flow.flow().time();
    let report = flow.last_step().copied();
    let euler_report = flow.flow().last_step().copied();
    let loads = BoundarySide::ALL.map(|side| flow.boundary_impulses(side).to_vec());
    let mut polls = 0;
    let expected = tensors(
        &flow
            .stress_field_moving(motion, || {
                polls += 1;
                false
            })
            .unwrap(),
    );
    for cutoff in [polls / 2, polls] {
        let mut calls = 0;
        assert_eq!(
            flow.stress_field_moving(motion, || {
                calls += 1;
                calls == cutoff
            })
            .err(),
            Some(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), report);
        assert_eq!(flow.flow().last_step().copied(), euler_report);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.boundary_impulses(side), loads[i]);
        }
        assert_eq!(
            tensors(&flow.stress_field_moving(motion, || false).unwrap()),
            expected
        );
    }
    let invalid = MeshMotion {
        lower_velocity: [0.0; 3],
        upper_velocity: [0.1, 0.2, 0.0],
    };
    assert!(flow.stress_field_moving(invalid, || false).is_err());
    assert_eq!(flow.flow().states(), states);
    let a = flow
        .step_moving(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    let b = reference
        .step_moving(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(flow.flow().states(), reference.flow().states());
}
