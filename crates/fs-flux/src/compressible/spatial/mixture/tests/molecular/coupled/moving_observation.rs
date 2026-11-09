use super::*;
use crate::compressible::spatial::mixture::molecular::{
    InertialWallLoad, InertialWallsControl, Piston,
};
use crate::compressible::spatial::{MeshMotion, molecular::StressField3d};

fn channel(axis: usize) -> Boundaries3d {
    let mut result = PERIODIC;
    match axis {
        0 => result.x = WALLS,
        1 => result.y = WALLS,
        2 => result.z = WALLS,
        _ => unreachable!(),
    }
    result
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

fn owner(boundaries: Boundaries3d, velocity: [f64; 3], walls: [MolecularWall3d; 6]) -> Coupled {
    let flow = field(
        equal_calorics(),
        CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap(),
        boundaries,
        |_| ([0.4, 0.6, 0.0], velocity, 2.0),
    );
    let mut molecular = prepare(&flow, law(flow.gas(), 0.01, 0.02, 0.009, 0.04));
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        if boundaries.axis(side.axis()) == WALLS {
            molecular = molecular.with_wall(side, walls[i]).unwrap();
        }
    }
    coupled(flow, molecular).0
}

fn tensors(view: &StressField3d<'_>) -> Vec<[f64; 6]> {
    (0..view.grid().len())
        .map(|cell| view.cell_stress_pa(cell).unwrap())
        .collect()
}

#[test]
fn g1_mixture_moving_observer_returns_original_normal_and_navier_tractions() {
    let velocity = [0.15, -0.1, 0.07];
    for (slot, side) in BoundarySide::ALL.into_iter().enumerate() {
        let axis = side.axis();
        let tangential: Vec<_> = (0..3).filter(|a| *a != axis).collect();
        let mut wall_velocity = [0.0; 3];
        wall_velocity[tangential[0]] = 0.4;
        wall_velocity[tangential[1]] = -0.3;
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = 0.23;
        motion.upper_velocity[axis] = 0.23;
        for kind in 0..3 {
            let mut walls = [MolecularWall3d::SlipAdiabatic; 6];
            walls[slot] = wall(kind, wall_velocity);
            let mut flow = owner(channel(axis), velocity, walls);
            let momentum = flow.flow().cell(0).unwrap().momentum;
            let distance = 0.5 * flow.flow().grid().spacing()[axis];
            let mut face = [1; 3];
            face[axis] = if side.sign() < 0.0 { 0 } else { 3 };
            let moving = {
                let view = flow
                    .stress_field_moving(DissipationControl::default(), motion, || false)
                    .unwrap();
                assert_eq!(view.momentum(0), Some(momentum));
                assert_eq!(view.time_s(), 0.0);
                view.face_traction_pa(axis, face).unwrap()
            };
            close(
                moving[axis],
                side.sign() * (4.0 * 0.02 / 3.0 + 0.009) * (0.23 - velocity[axis]) / distance,
                4e-13,
            );
            for &t in &tangential {
                let expected = if kind == 0 {
                    0.0
                } else {
                    side.sign() * 0.02 * (wall_velocity[t] - velocity[t])
                        / (distance + if kind == 1 { 0.0 } else { 0.12 })
                };
                close(moving[t], expected, 4e-13);
            }
            let stationary = flow
                .stress_field(DissipationControl::default(), || false)
                .unwrap()
                .face_traction_pa(axis, face)
                .unwrap();
            close(
                stationary[axis],
                -side.sign() * (4.0 * 0.02 / 3.0 + 0.009) * velocity[axis] / distance,
                4e-13,
            );
            let zero = flow
                .stress_field_moving(DissipationControl::default(), MeshMotion::default(), || {
                    false
                })
                .unwrap()
                .face_traction_pa(axis, face)
                .unwrap();
            assert_eq!(zero, stationary);
        }
    }
}

#[test]
fn g1_mixture_comoving_stress_observation_respects_each_original_slip_law() {
    let velocity = [0.3, -0.2, 0.15];
    let motion = MeshMotion {
        lower_velocity: velocity,
        upper_velocity: velocity,
    };
    let boundaries = Boundaries3d {
        x: WALLS,
        y: WALLS,
        z: WALLS,
    };
    for kind in 0..3 {
        let walls = BoundarySide::ALL.map(|side| {
            let mut tangent_velocity = velocity;
            tangent_velocity[side.axis()] = 0.0;
            wall(kind, tangent_velocity)
        });
        let mut flow = owner(boundaries, velocity, walls);
        flow.step_moving(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
        let grid = *flow.flow().grid();
        let time = flow.flow().time().value();
        let bytes = flow.workspace_bytes();
        {
            let view = flow
                .stress_field_moving(DissipationControl::default(), motion, || false)
                .unwrap();
            assert_eq!(*view.grid(), grid);
            assert_eq!(view.time_s(), time);
            for tensor in tensors(&view) {
                for component in tensor {
                    close(component, 0.0, 4e-13);
                }
            }
        }
        assert_eq!(flow.workspace_bytes(), bytes);
    }
}

#[test]
fn g1_mixture_piston_stress_observer_accepts_current_endpoint_velocity_explicitly() {
    for axis in 0..3 {
        let mut flow = owner(channel(axis), [0.0; 3], [MolecularWall3d::SlipAdiabatic; 6]);
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
        let half_width = 0.5 * grid.spacing()[axis];
        let observed = {
            let view = flow
                .stress_field_moving(DissipationControl::default(), endpoint, || false)
                .unwrap();
            std::array::from_fn::<_, 2, _>(|end| {
                let mut face = [1; 3];
                face[axis] = 3 * end;
                let mut cell = face;
                cell[axis] = 2 * end;
                let p = view.primitive(grid.index(cell)).unwrap();
                let sign = if end == 0 { -1.0 } else { 1.0 };
                let actual = view.face_traction_pa(axis, face).unwrap()[axis];
                close(
                    actual,
                    sign * (4.0 * 0.02 / 3.0 + 0.009)
                        * (pistons[end].velocity_m_s() - p.velocity[axis])
                        / half_width,
                    5e-13,
                );
                actual
            })
        };
        let held = flow
            .stress_field_moving(DissipationControl::default(), midpoint, || false)
            .unwrap();
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
                5e-13,
            );
            assert!(difference.abs() > 1e-6);
        }
    }
}

#[test]
fn g4_mixture_moving_observation_is_atomic_for_species_geometry_and_all_loads() {
    let motion = MeshMotion {
        lower_velocity: [-0.05, 0.0, 0.0],
        upper_velocity: [0.1, 0.0, 0.0],
    };
    let make = || {
        owner(
            channel(0),
            [0.1, -0.07, 0.04],
            [MolecularWall3d::SlipAdiabatic; 6],
        )
    };
    let mut flow = make();
    let mut reference = make();
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
    let combined = loads(&flow);
    let molecular = molecular_loads(&flow);
    let bytes = flow.workspace_bytes();
    let mut polls = 0;
    let expected = tensors(
        &flow
            .stress_field_moving(DissipationControl::default(), motion, || {
                polls += 1;
                false
            })
            .unwrap(),
    );
    for cutoff in [polls / 2, polls] {
        let mut calls = 0;
        assert_eq!(
            flow.stress_field_moving(DissipationControl::default(), motion, || {
                calls += 1;
                calls == cutoff
            })
            .err(),
            Some(DissipationError::Cancelled)
        );
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), report);
        assert_eq!(flow.flow().last_step().copied(), euler_report);
        assert_eq!(loads(&flow), combined);
        assert_eq!(molecular_loads(&flow), molecular);
        assert_eq!(flow.workspace_bytes(), bytes);
        assert_eq!(
            tensors(
                &flow
                    .stress_field_moving(DissipationControl::default(), motion, || false)
                    .unwrap()
            ),
            expected
        );
    }
    let invalid = MeshMotion {
        lower_velocity: [0.0; 3],
        upper_velocity: [0.1, 0.2, 0.0],
    };
    assert!(
        flow.stress_field_moving(DissipationControl::default(), invalid, || false)
            .is_err()
    );
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
