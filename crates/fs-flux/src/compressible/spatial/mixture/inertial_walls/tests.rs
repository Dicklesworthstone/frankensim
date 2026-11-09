use super::*;
use crate::compressible::mixture::tests::{close, gas, source, standard};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

mod forced;

const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const WALLS: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};

fn field(
    gas: MixtureGas,
    grid: CartesianGrid,
    sample: impl Fn([f64; 3]) -> ([f64; 3], [f64; 3], f64),
) -> MixtureEuler3d {
    let mut states = vec![0.0; grid.len() * 7];
    for (k, row) in states.chunks_exact_mut(7).enumerate() {
        let (species, velocity, pressure) = sample(grid.cell_center(k).unwrap());
        state_into(&gas, &species, velocity, pressure, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, WALLS, bytes, || false).unwrap()
}

fn controls(boost: f64) -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[boost - 1.0, boost + 1.0]; 2],
        momentum_tolerance_n_s: 2e-11,
        energy_tolerance_j: 2e-11,
        maximum_iterations: 32,
    }
}

fn walls(boost: f64) -> [Piston; 2] {
    [
        Piston::new(0.3, boost).unwrap(),
        Piston::new(0.6, boost).unwrap(),
    ]
}

fn kinetic(walls: &[Piston; 2]) -> f64 {
    walls
        .iter()
        .map(|w| 0.5 * w.mass_kg() * w.velocity_m_s().powi(2))
        .sum()
}

fn momentum(walls: &[Piston; 2]) -> f64 {
    walls.iter().map(|w| w.mass_kg() * w.velocity_m_s()).sum()
}

fn balance(step: SpatialInertialWallsStep<StepReport>, controls: InertialWallsControl) {
    let m = step.mechanics;
    assert!(m.combined_momentum_residual_n_s.abs() <= controls.momentum_tolerance_n_s);
    assert!(m.total_energy_residual_j.abs() <= controls.energy_tolerance_j);
    for end in 0..2 {
        assert!(m.momentum_residual_n_s[end].abs() <= controls.momentum_tolerance_n_s);
        assert!(m.mechanical_energy_residual_j[end].abs() <= controls.energy_tolerance_j);
        for d in 0..3 {
            assert_eq!(
                step.gas_angular_impulse_n_m_s[end][d],
                -step.guide_angular_impulse_n_m_s[end][d]
            );
        }
    }
    for x in m.gas.conservation_residual.species_masses {
        close(x, 0.0, 5e-12);
    }
    for x in m.gas.conservation_residual.momentum {
        close(x, 0.0, 5e-12);
    }
    close(m.gas.conservation_residual.energy, 0.0, 6e-12);
    close(m.gas.absolute_energy_residual_j, 0.0, 3e-10);
    let referenced = sum(
        std::iter::once(m.gas.after.energy - m.gas.before.energy).chain((0..2).map(|i| {
            m.kinetic_energy_change_j[i] + m.exterior_pressure_work_j[i]
                - m.external_force_work_j[i]
        })),
    );
    assert!(referenced.abs() <= controls.energy_tolerance_j);
}

#[test]
fn g1_unequal_caloric_species_release_both_finite_planes_with_closed_balances() {
    let grid = CartesianGrid::new(4, 3, 2, 0.2, 0.3, 0.4).unwrap();
    for axis in 0..3 {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut flow = field(standard(0.4), grid, |x| {
                (
                    [0.6 + 0.1 * x[0], 0.4 - 0.05 * x[1], 0.0],
                    [0.0; 3],
                    1.0 + 0.2 * x[1],
                )
            });
            let mut planes = walls(0.0);
            let original = flow.totals();
            let physical_energy = original.absolute_energy(flow.gas()).unwrap() + kinetic(&planes);
            let axial_momentum = original.momentum[axis] + momentum(&planes);
            let mut found_torque = false;
            for _ in 0..6 {
                let report = flow
                    .step_with_inertial_walls(
                        axis,
                        &mut planes,
                        Time::new(0.002),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        [InertialWallLoad::default(); 2],
                        controls(0.0),
                        || false,
                    )
                    .unwrap();
                balance(report, controls(0.0));
                assert_eq!(report.axis, axis);
                assert_eq!(flow.last_step(), Some(&report.mechanics.gas));
                assert!(planes[0].velocity_m_s() < 0.0);
                assert!(planes[1].velocity_m_s() > 0.0);
                for d in 0..3 {
                    close(
                        flow.totals().species_masses[d],
                        original.species_masses[d],
                        5e-13,
                    );
                }
                for k in 0..flow.cell_count() {
                    assert_eq!(flow.cell(k).unwrap().partial_densities[2], 0.0);
                }
                close(
                    flow.totals().momentum[axis] + momentum(&planes),
                    axial_momentum,
                    2e-10,
                );
                close(
                    flow.totals().absolute_energy(flow.gas()).unwrap() + kinetic(&planes),
                    physical_energy,
                    2e-10,
                );
                for side in BoundarySide::ALL {
                    for transfer in flow.boundary_impulses(side) {
                        assert!(transfer.species_masses.iter().all(|x| *x == 0.0));
                    }
                }
                found_torque |= report
                    .gas_angular_impulse_n_m_s
                    .iter()
                    .flatten()
                    .any(|x| x.abs() > 1e-8);
            }
            assert!(flow.grid().spacing()[axis] > grid.spacing()[axis]);
            if axis != 1 {
                assert!(found_torque);
            }
        }
    }
}

#[test]
fn g0_equal_pressure_gives_a_translating_equilibrium_for_each_axis() {
    let grid = CartesianGrid::new(2, 2, 2, 0.5, 0.5, 0.5).unwrap();
    for axis in 0..3 {
        let boost = 0.25;
        let mut velocity = [0.0; 3];
        velocity[axis] = boost;
        let mut flow = field(standard(0.4), grid, |_| ([0.6, 0.4, 0.0], velocity, 1.0));
        let initial = flow.states().to_vec();
        let mut planes = walls(boost);
        let report = flow
            .step_with_inertial_walls(
                axis,
                &mut planes,
                Time::new(0.01),
                StepControl::default(),
                [InertialWallLoad {
                    exterior_pressure_pa: 1.0,
                    external_force_n: 0.0,
                }; 2],
                controls(boost),
                || false,
            )
            .unwrap();
        balance(report, controls(boost));
        assert_eq!(report.mechanics.iterations, 0);
        for plane in planes {
            close(plane.velocity_m_s(), boost, 5e-13);
        }
        for (a, b) in flow.states().iter().zip(initial) {
            close(*a, b, 5e-13);
        }
        close(
            flow.grid().origin()[axis],
            boost * report.mechanics.gas.elapsed.value(),
            3e-15,
        );
        for d in 0..3 {
            close(flow.grid().spacing()[d], grid.spacing()[d], 3e-15);
        }
    }
}

#[test]
fn g3_finite_mixture_planes_are_galilean_covariant_with_external_work() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let boost = 0.4;
    let sample = |x: [f64; 3], speed: f64| {
        (
            [0.6 + 0.1 * x[0], 0.4 - 0.05 * x[1], 0.0],
            [speed + 0.12 + 0.07 * x[0], 0.02 * x[0], -0.01 * x[1]],
            1.0 + 0.2 * x[2],
        )
    };
    let mut original = field(standard(0.4), grid, |x| sample(x, 0.0));
    let mut shifted = field(standard(0.4), grid, |x| sample(x, boost));
    let mut a_walls = walls(0.0);
    let mut b_walls = walls(boost);
    let loads = [
        InertialWallLoad {
            exterior_pressure_pa: 0.4,
            external_force_n: 0.03,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.7,
            external_force_n: -0.02,
        },
    ];
    let a = original
        .step_with_inertial_walls(
            0,
            &mut a_walls,
            Time::new(0.001),
            StepControl::default(),
            loads,
            controls(0.0),
            || false,
        )
        .unwrap();
    let b = shifted
        .step_with_inertial_walls(
            0,
            &mut b_walls,
            Time::new(0.001),
            StepControl::default(),
            loads,
            controls(boost),
            || false,
        )
        .unwrap();
    balance(a, controls(0.0));
    balance(b, controls(boost));
    close(
        a.mechanics.gas.elapsed.value(),
        b.mechanics.gas.elapsed.value(),
        3e-15,
    );
    let dt = a.mechanics.gas.elapsed.value();
    for end in 0..2 {
        close(
            b_walls[end].velocity_m_s(),
            a_walls[end].velocity_m_s() + boost,
            3e-10,
        );
        close(
            b.mechanics.midpoint_velocity_m_s[end],
            a.mechanics.midpoint_velocity_m_s[end] + boost,
            3e-10,
        );
        close(
            b.mechanics.gas_impulse_n_s[end],
            a.mechanics.gas_impulse_n_s[end],
            5e-12,
        );
        close(
            b.mechanics.external_force_work_j[end] - a.mechanics.external_force_work_j[end],
            loads[end].external_force_n * boost * dt,
            3e-12,
        );
    }
    for k in 0..original.cell_count() {
        let q = original.cell(k).unwrap();
        let shifted_q = shifted.cell(k).unwrap();
        for (a, b) in q.partial_densities.iter().zip(shifted_q.partial_densities) {
            close(*a, *b, 5e-10);
        }
        let mass = sum(q.partial_densities.iter().copied());
        close(shifted_q.momentum[0], q.momentum[0] + boost * mass, 5e-10);
        for d in 1..3 {
            close(shifted_q.momentum[d], q.momentum[d], 5e-10);
        }
        close(
            shifted_q.energy,
            q.energy + boost * q.momentum[0] + 0.5 * boost * boost * mass,
            6e-10,
        );
    }
    close(
        shifted.grid().origin()[0],
        original.grid().origin()[0] + boost * dt,
        5e-12,
    );
}

#[test]
fn g3_reference_gauges_do_not_change_wall_recoil_or_species_transport() {
    let grid = CartesianGrid::new(4, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let models = |shift: f64| {
        gas(vec![
            source("A", R, 3.5, 0.0, 30.0 + shift),
            source("B", 2.0 * R, 3.5, 0.4, -60.0 - 2.0 * shift),
            source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
        ])
    };
    let sample = |x: [f64; 3]| ([0.6 + 0.1 * x[0], 0.4 - 0.1 * x[0], 0.0], [0.0; 3], 1.0);
    let mut original = field(models(0.0), grid, sample);
    let mut shifted = field(models(1024.0), grid, sample);
    let mut a_walls = walls(0.0);
    let mut b_walls = walls(0.0);
    for _ in 0..3 {
        let a = original
            .step_with_inertial_walls(
                0,
                &mut a_walls,
                Time::new(0.002),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                controls(0.0),
                || false,
            )
            .unwrap();
        let b = shifted
            .step_with_inertial_walls(
                0,
                &mut b_walls,
                Time::new(0.002),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                controls(0.0),
                || false,
            )
            .unwrap();
        balance(a, controls(0.0));
        balance(b, controls(0.0));
        close(
            a.mechanics.gas.elapsed.value(),
            b.mechanics.gas.elapsed.value(),
            3e-15,
        );
        for end in 0..2 {
            close(
                a_walls[end].velocity_m_s(),
                b_walls[end].velocity_m_s(),
                3e-10,
            );
        }
    }
    for k in 0..original.cell_count() {
        let a = original.cell(k).unwrap();
        let b = shifted.cell(k).unwrap();
        for (x, y) in a.partial_densities.iter().zip(b.partial_densities) {
            close(*x, *y, 3e-10);
        }
        let p = a.primitive(original.gas(), || false).unwrap();
        let q = b.primitive(shifted.gas(), || false).unwrap();
        close(p.pressure, q.pressure, 3e-10);
        close(p.temperature, q.temperature, 3e-10);
    }
}

#[test]
fn g4_cancelled_coupled_proposals_and_final_poll_leave_both_owners_unchanged() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let make = || {
        field(standard(0.4), grid, |x| {
            ([0.6 + 0.05 * x[0], 0.4, 0.0], [0.0; 3], 1.0)
        })
    };
    let mut control = make();
    let mut refused = make();
    let mut a_walls = walls(0.0);
    let mut b_walls = walls(0.0);
    for (flow, planes) in [(&mut control, &mut a_walls), (&mut refused, &mut b_walls)] {
        flow.step_with_inertial_walls(
            0,
            planes,
            Time::new(0.001),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || false,
        )
        .unwrap();
    }
    let state = refused.states().to_vec();
    let grid_before = *refused.grid();
    let time = refused.time();
    let last = *refused.last_step().unwrap();
    let transfers = refused.boundary_impulses.clone();
    let planes_before = b_walls;
    let mut polls = 0;
    let expected = control
        .step_with_inertial_walls(
            0,
            &mut a_walls,
            Time::new(0.002),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            refused.step_with_inertial_walls(
                0,
                &mut b_walls,
                Time::new(0.002),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                controls(0.0),
                || {
                    count += 1;
                    count == stop
                }
            ),
            Err(PistonError::Cancelled)
        );
        assert_eq!(b_walls, planes_before);
        assert_eq!(refused.states(), state);
        assert_eq!(*refused.grid(), grid_before);
        assert_eq!(refused.time(), time);
        assert_eq!(*refused.last_step().unwrap(), last);
        assert_eq!(refused.boundary_impulses, transfers);
    }
    let replay = refused
        .step_with_inertial_walls(
            0,
            &mut b_walls,
            Time::new(0.002),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(refused.states(), control.states());
    assert_eq!(refused.grid(), control.grid());
    assert_eq!(a_walls, b_walls);
}

#[test]
fn g4_invalid_axes_open_sides_and_unrepresentable_wall_travel_refuse() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let make = || field(standard(0.4), grid, |_| ([0.6, 0.4, 0.0], [0.0; 3], 1.0));
    let mut flow = make();
    let mut planes = walls(0.0);
    let initial = flow.states().to_vec();
    let old_planes = planes;
    assert!(
        flow.step_with_inertial_walls(
            3,
            &mut planes,
            Time::new(0.001),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || false
        )
        .is_err()
    );
    flow.boundaries.y = AxisBoundaries::Ends {
        lower: BoundaryCondition::Transmissive,
        upper: BoundaryCondition::Reflecting,
    };
    assert!(
        flow.step_with_inertial_walls(
            0,
            &mut planes,
            Time::new(0.001),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), initial);
    assert_eq!(*flow.grid(), grid);
    assert_eq!(planes, old_planes);
    let mut lost = make();
    lost.grid.origin[0] = 1e12;
    let old_grid = *lost.grid();
    let old_state = lost.states().to_vec();
    let mut moving_planes = walls(1e-8);
    let old_walls = moving_planes;
    assert!(
        lost.step_with_inertial_walls(
            0,
            &mut moving_planes,
            Time::new(0.001),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            [InertialWallLoad::default(); 2],
            controls(1e-8),
            || false
        )
        .is_err()
    );
    assert_eq!(*lost.grid(), old_grid);
    assert_eq!(lost.states(), old_state);
    assert_eq!(moving_planes, old_walls);
    assert_eq!(lost.time().value(), 0.0);
    assert!(lost.last_step().is_none());
}
