use super::*;

mod forced;
use crate::compressible::spatial::mixture::molecular::{
    InertialWallLoad, InertialWallsControl, MolecularInertialWallsError,
    MolecularInertialWallsStep, Piston,
};

fn controls() -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[-1.0, 1.0]; 2],
        momentum_tolerance_n_s: 1e-10,
        energy_tolerance_j: 1e-10,
        maximum_iterations: 32,
    }
}
fn pistons() -> [Piston; 2] {
    [
        Piston::new(0.3, 0.0).unwrap(),
        Piston::new(0.6, 0.0).unwrap(),
    ]
}
fn external_loads() -> [InertialWallLoad; 2] {
    [
        InertialWallLoad {
            exterior_pressure_pa: 0.8,
            external_force_n: 0.02,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.7,
            external_force_n: -0.01,
        },
    ]
}
fn case(axis: usize, reference_shift: f64) -> Coupled {
    let model = gas(vec![
        source("A", R, 3.5, 0.0, 30.0 + reference_shift),
        source("B", 2.0 * R, 3.5, 0.4, -60.0 - 2.0 * reference_shift),
        source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
    ]);
    let flow = field(
        model,
        CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap(),
        Boundaries3d {
            x: WALLS,
            y: WALLS,
            z: WALLS,
        },
        |x| {
            let a = 0.5 + 0.1 * x[0] - 0.04 * x[1];
            (
                [a, 1.0 - a, 0.0],
                [0.07 + 0.03 * x[1], -0.05 + 0.02 * x[2], 0.04 - 0.02 * x[0]],
                1.0 + 0.03 * x[2],
            )
        },
    );
    let mut molecular = prepare(&flow, law(flow.gas(), 0.015, 0.07, 0.02, 0.04));
    for side in BoundarySide::ALL {
        let thermal = if side == BoundarySide::ALL[2 * axis + 1] {
            WallThermalCondition::Isothermal(Temperature::new(1.8))
        } else if side == BoundarySide::ALL[2 * ((axis + 1) % 3)] {
            WallThermalCondition::Isothermal(Temperature::new(0.9))
        } else {
            WallThermalCondition::Adiabatic
        };
        let wall = if side == BoundarySide::ALL[2 * axis + 1] {
            MolecularWall3d::NavierSlip {
                velocity: [0.0; 3],
                thermal,
                slip_length_m: 0.08,
                fluid_heat_fraction: 0.6,
            }
        } else {
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal,
            }
        };
        molecular = molecular.with_wall(side, wall).unwrap();
    }
    coupled(flow, molecular).0
}

fn check_pistons(report: MolecularInertialWallsStep, gas: &MixtureGas) {
    let m = report.mechanics;
    check(&m.gas.transport);
    let c = controls();
    assert!(m.combined_momentum_residual_n_s.abs() <= c.momentum_tolerance_n_s);
    assert!(m.total_energy_residual_j.abs() <= c.energy_tolerance_j);
    for end in 0..2 {
        assert!(m.momentum_residual_n_s[end].abs() <= c.momentum_tolerance_n_s);
        assert!(m.mechanical_energy_residual_j[end].abs() <= c.energy_tolerance_j);
        assert_eq!(report.guide_impulse_n_s[end][report.axis], 0.0);
        for a in 0..3 {
            assert_eq!(
                report.guide_angular_impulse_n_m_s[end][a],
                -report.gas_angular_impulse_n_m_s[end][a]
            );
        }
    }
    let heat = m.gas.wall_heat_j.iter().sum::<f64>();
    close(
        m.gas.transport.first.boundary_change.energy
            + m.gas.transport.second.boundary_change.energy,
        heat + m.gas.molecular_wall_work_j.iter().sum::<f64>(),
        4e-11,
    );
    let before = m.gas.transport.before;
    let after = m.gas.transport.after;
    let wall_change = (0..2)
        .map(|i| {
            m.kinetic_energy_change_j[i] + m.exterior_pressure_work_j[i]
                - m.external_force_work_j[i]
        })
        .sum::<f64>();
    let referenced = after.energy - before.energy + wall_change - heat;
    assert!(referenced.abs() <= c.energy_tolerance_j);
    let chemical = (0..gas.species_count())
        .map(|k| (after.species_masses[k] - before.species_masses[k]) * gas.energy_references()[k])
        .sum::<f64>();
    assert!((referenced + chemical).abs() <= c.energy_tolerance_j);
    for axis in 0..3 {
        if axis != report.axis {
            close(
                after.momentum[axis]
                    - before.momentum[axis]
                    - report.side_boundary_impulse_n_s[axis]
                    - report.guide_impulse_n_s[0][axis]
                    - report.guide_impulse_n_s[1][axis],
                0.0,
                2e-11,
            );
        }
    }
}

#[test]
fn g1_finite_mixture_pistons_close_all_momenta_reference_energy_heat_and_guides() {
    for axis in 0..3 {
        let mut flow = case(axis, 0.0);
        let bytes = flow.workspace_bytes();
        let initial_grid = *flow.flow().grid();
        let initial = flow.flow().states().to_vec();
        let mut walls = pistons();
        let report: MolecularInertialWallsStep = flow
            .step_with_inertial_walls(
                axis,
                &mut walls,
                Time::new(0.002),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || false,
            )
            .unwrap();
        check_pistons(report, flow.flow().gas());
        assert_eq!(flow.last_step(), Some(&report.mechanics.gas.transport));
        assert_eq!(flow.workspace_bytes(), bytes);
        assert_ne!(flow.flow().states(), initial);
        assert!(report.side_boundary_impulse_n_s[axis].abs() > 1e-7);
        assert!(
            report
                .guide_impulse_n_s
                .iter()
                .flatten()
                .any(|x| x.abs() > 1e-7)
        );
        assert!(
            report
                .gas_angular_impulse_n_m_s
                .iter()
                .flatten()
                .any(|x| x.abs() > 1e-8)
        );
        assert!(
            report
                .mechanics
                .gas
                .molecular_wall_work_j
                .iter()
                .any(|x| x.abs() > 1e-12)
        );
        assert!(report.mechanics.gas.wall_heat_j[2 * axis + 1] > 0.0);
        assert!(report.mechanics.gas.wall_heat_j[2 * ((axis + 1) % 3)] < 0.0);
        for side in BoundarySide::ALL {
            for q in flow.boundary_impulses(side) {
                assert!(q.species_masses.iter().all(|x| *x == 0.0));
            }
        }
        let length = |grid: &CartesianGrid| grid.shape()[axis] as f64 * grid.spacing()[axis];
        close(
            length(flow.flow().grid()) - length(&initial_grid),
            report.mechanics.displacement_m[1] - report.mechanics.displacement_m[0],
            3e-14,
        );
        for end in 0..2 {
            close(
                walls[end].velocity_m_s(),
                report.mechanics.velocity_after_m_s[end],
                1e-14,
            );
        }
    }
}

#[test]
fn g1_original_formation_references_do_not_change_viscous_piston_response() {
    let mut a = case(0, 0.0);
    let mut b = case(0, 1024.0);
    let mut wa = pistons();
    let mut wb = pistons();
    for _ in 0..2 {
        let ra = a
            .step_with_inertial_walls(
                0,
                &mut wa,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || false,
            )
            .unwrap();
        let rb = b
            .step_with_inertial_walls(
                0,
                &mut wb,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || false,
            )
            .unwrap();
        check_pistons(ra, a.flow().gas());
        check_pistons(rb, b.flow().gas());
        for end in 0..2 {
            close(wa[end].velocity_m_s(), wb[end].velocity_m_s(), 3e-10);
            close(
                ra.mechanics.gas_impulse_n_s[end],
                rb.mechanics.gas_impulse_n_s[end],
                3e-11,
            );
        }
    }
    for index in 0..a.flow().cell_count() {
        let qa = a.flow().cell(index).unwrap();
        let qb = b.flow().cell(index).unwrap();
        for species in 0..3 {
            close(
                qa.partial_densities[species],
                qb.partial_densities[species],
                3e-10,
            );
        }
        let pa = qa.primitive(a.flow().gas(), || false).unwrap();
        let pb = qb.primitive(b.flow().gas(), || false).unwrap();
        close(pa.pressure, pb.pressure, 3e-10);
        close(pa.temperature, pb.temperature, 3e-10);
        for axis in 0..3 {
            close(pa.velocity[axis], pb.velocity[axis], 3e-10);
        }
    }
}

#[test]
fn g0_zero_molecular_coefficients_recover_original_finite_piston_dynamics() {
    for axis in 0..3 {
        let make = || {
            field(
                standard(0.4),
                CartesianGrid::new(2, 2, 2, 0.4, 0.4, 0.4).unwrap(),
                Boundaries3d {
                    x: WALLS,
                    y: WALLS,
                    z: WALLS,
                },
                |_| ([0.6, 0.4, 0.0], [0.0; 3], 1.0),
            )
        };
        let mut inviscid = make();
        let original = make();
        let molecular = prepare(&original, law(original.gas(), 1e-20, 0.0, 0.0, 0.0));
        let mut flow = coupled(original, molecular).0;
        let mut wa = pistons();
        let mut wb = pistons();
        let a = inviscid
            .step_with_inertial_walls(
                axis,
                &mut wa,
                Time::new(0.002),
                StepControl::default(),
                external_loads(),
                controls(),
                || false,
            )
            .unwrap();
        let b = flow
            .step_with_inertial_walls(
                axis,
                &mut wb,
                Time::new(0.002),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || false,
            )
            .unwrap();
        check_pistons(b, flow.flow().gas());
        for end in 0..2 {
            close(wa[end].velocity_m_s(), wb[end].velocity_m_s(), 3e-11);
        }
        for (x, y) in inviscid.states().iter().zip(flow.flow().states()) {
            close(*x, *y, 3e-11);
        }
        close(
            a.mechanics.gas.elapsed.value(),
            b.mechanics.gas.transport.elapsed.value(),
            2e-15,
        );
        assert_eq!(b.mechanics.gas.molecular_wall_work_j, [0.0; 6]);
        assert_eq!(b.mechanics.gas.wall_heat_j, [0.0; 6]);
    }
}

#[test]
fn g3_failed_or_cancelled_mixture_piston_trials_publish_no_partial_owner() {
    let mut expected = case(0, 0.0);
    let mut expected_walls = pistons();
    let mut polls = 0;
    let accepted = expected
        .step_with_inertial_walls(
            0,
            &mut expected_walls,
            Time::new(0.002),
            MolecularStepControl::default(),
            external_loads(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut flow = case(0, 0.0);
        let mut walls = pistons();
        let original_walls = walls;
        let states = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let boundary = loads(&flow);
        let molecular_boundary = molecular_loads(&flow);
        let bytes = flow.workspace_bytes();
        let mut seen = 0;
        let result: Result<MolecularInertialWallsStep, MolecularInertialWallsError> = flow
            .step_with_inertial_walls(
                0,
                &mut walls,
                Time::new(0.002),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || {
                    seen += 1;
                    seen == stop
                },
            );
        assert!(matches!(
            result,
            Err(MolecularInertialWallsError::Cancelled)
        ));
        assert_eq!(walls, original_walls);
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time().value(), 0.0);
        assert!(flow.last_step().is_none());
        assert!(flow.molecular().last_step().is_none());
        assert_eq!(loads(&flow), boundary);
        assert_eq!(molecular_loads(&flow), molecular_boundary);
        assert_eq!(flow.workspace_bytes(), bytes);
        let recovered = flow
            .step_with_inertial_walls(
                0,
                &mut walls,
                Time::new(0.002),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(recovered, accepted);
        assert_eq!(walls, expected_walls);
        assert_eq!(flow.flow().states(), expected.flow().states());
    }
    let (mut driven, _) = wall_case();
    let original = driven.flow().states().to_vec();
    let mut walls = pistons();
    assert!(
        driven
            .step_with_inertial_walls(
                0,
                &mut walls,
                Time::new(0.002),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(driven.flow().states(), original);
    assert_eq!(walls, pistons());
}
