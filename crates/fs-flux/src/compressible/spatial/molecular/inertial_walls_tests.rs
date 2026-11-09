use super::*;
use crate::compressible::spatial::molecular::{
    InertialWallLoad, InertialWallsControl, MolecularInertialWallsError,
    MolecularInertialWallsStep, Piston,
};

fn controls(velocity: [f64; 2]) -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: velocity.map(|v| [v - 2.0, v + 2.0]),
        momentum_tolerance_n_s: 3e-11,
        energy_tolerance_j: 3e-11,
        maximum_iterations: 24,
    }
}

fn channel(axis: usize) -> Boundaries3d {
    let mut boundaries = PERIODIC;
    match axis {
        0 => boundaries.x = WALL,
        1 => boundaries.y = WALL,
        2 => boundaries.z = WALL,
        _ => unreachable!(),
    }
    boundaries
}

fn grid() -> CartesianGrid {
    CartesianGrid::new(4, 3, 2, 0.25, 0.3, 0.4).unwrap()
}

fn owner(
    boundaries: Boundaries3d,
    mu: f64,
    bulk: f64,
    conductivity: f64,
    sample: impl Fn([f64; 3]) -> [f64; 5],
) -> MolecularEuler3d<ConstantTransport> {
    let flow = field(gas(0.0), grid(), boundaries, sample);
    let transport = law(flow.gas(), mu, bulk, conductivity);
    MolecularEuler3d::new(flow, transport, || false).unwrap()
}

fn check(report: MolecularInertialWallsStep, control: InertialWallsControl) {
    let mechanics = report.mechanics;
    let gas = mechanics.gas.transport;
    assert_eq!(gas.first.elapsed.value() * 2.0, gas.elapsed.value());
    assert_eq!(gas.second.elapsed, gas.first.elapsed);
    assert_eq!(gas.advection.elapsed, gas.elapsed);
    assert!(gas.first.maximum_parabolic_number <= 0.40000000000001);
    assert!(gas.second.maximum_parabolic_number <= 0.40000000000001);
    assert!(gas.advection.maximum_cfl <= 0.40000000000001);
    for residual in gas.conservation_residual.array() {
        close(residual, 0.0, 2e-12);
    }
    for residual in mechanics
        .momentum_residual_n_s
        .into_iter()
        .chain([mechanics.combined_momentum_residual_n_s])
    {
        assert!(residual.abs() <= control.momentum_tolerance_n_s);
    }
    for residual in mechanics
        .mechanical_energy_residual_j
        .into_iter()
        .chain([mechanics.total_energy_residual_j])
    {
        assert!(residual.abs() <= control.energy_tolerance_j);
    }
}

fn kinetic(wall: Piston) -> f64 {
    0.5 * wall.mass_kg() * wall.velocity_m_s().powi(2)
}

#[test]
fn g0_zero_transport_finite_planes_reproduce_the_native_inviscid_owner_on_all_axes() {
    for axis in 0..3 {
        let sample = |x: [f64; 3]| {
            let mut u = [0.0; 3];
            u[axis] = 0.08 + 0.06 * x[axis];
            [1.0 + 0.02 * x[axis], u[0], u[1], u[2], 2.0 + 0.1 * x[axis]]
        };
        let mut viscous = owner(channel(axis), 0.0, 0.0, 0.0, sample);
        let mut inviscid = field(gas(0.0), grid(), channel(axis), sample);
        let mut a = [
            Piston::new(0.8, -0.1).unwrap(),
            Piston::new(1.7, 0.15).unwrap(),
        ];
        let mut b = a;
        let control = controls(a.map(|wall| wall.velocity_m_s()));
        let loads = [InertialWallLoad {
            exterior_pressure_pa: 1.0,
            external_force_n: 0.0,
        }; 2];
        let report = viscous
            .step_with_inertial_walls(
                axis,
                &mut a,
                Time::new(0.001),
                MolecularStepControl::default(),
                loads,
                control,
                || false,
            )
            .unwrap();
        let baseline = inviscid
            .step_with_inertial_walls(
                axis,
                &mut b,
                Time::new(0.001),
                StepControl::default(),
                loads,
                control,
                || false,
            )
            .unwrap();
        check(report, control);
        assert_eq!(a, b);
        assert_eq!(viscous.flow().states(), inviscid.states());
        assert_eq!(viscous.flow().grid(), inviscid.grid());
        assert_eq!(
            report.mechanics.gas.transport.advection,
            baseline.mechanics.gas
        );
        assert_eq!(
            report.mechanics.gas_impulse_n_s,
            baseline.mechanics.gas_impulse_n_s
        );
        assert_eq!(report.mechanics.gas.wall_heat_j, [0.0; 6]);
        assert_eq!(report.mechanics.gas.molecular_wall_work_j, [0.0; 6]);
    }
}

#[test]
fn g1_finite_piston_recoil_includes_original_longitudinal_viscous_traction() {
    let (mu, bulk, speed) = (0.04, 0.02, 0.3);
    for axis in 0..3 {
        let sample = |_| {
            let mut u = [0.0; 3];
            u[axis] = speed;
            [1.0, u[0], u[1], u[2], 2.0]
        };
        let mut native = owner(channel(axis), mu, bulk, 0.0, sample);
        let mut pressure_only = owner(channel(axis), 0.0, 0.0, 0.0, sample);
        let mut walls = [
            Piston::new(2.0, 0.0).unwrap(),
            Piston::new(5.0, 0.0).unwrap(),
        ];
        let mut baseline_walls = walls;
        let control = controls([0.0; 2]);
        let report = native
            .step_with_inertial_walls(
                axis,
                &mut walls,
                Time::new(1e-5),
                MolecularStepControl::default(),
                [InertialWallLoad::default(); 2],
                control,
                || false,
            )
            .unwrap();
        let pressure = pressure_only
            .step_with_inertial_walls(
                axis,
                &mut baseline_walls,
                Time::new(1e-5),
                MolecularStepControl::default(),
                [InertialWallLoad::default(); 2],
                control,
                || false,
            )
            .unwrap();
        check(report, control);
        let side = BoundarySide::ALL[2 * axis];
        let area = grid().face_area(side) * grid().face_count(side) as f64;
        // At t=0, both planes receive +eta*u/(h/2)*A from the resolved
        // normal viscous stress. The pressure-only run removes that term.
        let expected_force = (4.0 * mu / 3.0 + bulk) * speed / (0.5 * grid().spacing[axis]) * area;
        let dt = report.mechanics.gas.transport.elapsed.value();
        for end in 0..2 {
            let force = (report.mechanics.gas_impulse_n_s[end]
                - pressure.mechanics.gas_impulse_n_s[end])
                / dt;
            assert!(force > 0.0);
            assert!(
                (force - expected_force).abs() < 0.005 * expected_force,
                "axis={axis} end={end} force={force} expected={expected_force}"
            );
            assert!(walls[end].velocity_m_s() > baseline_walls[end].velocity_m_s());
        }
    }
}

fn heated_channel() -> MolecularEuler3d<ConstantTransport> {
    let mut result = owner(channel(0), 0.035, 0.01, 0.07, |_| {
        [1.0, 0.05, 0.08, -0.03, 2.0]
    });
    for side in [BoundarySide::XMin, BoundarySide::XMax] {
        result = result
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                    slip_length_m: 0.08,
                    fluid_heat_fraction: 0.4,
                },
            )
            .unwrap();
    }
    result
}

#[test]
fn g1_isothermal_heat_remains_distinct_from_reciprocal_finite_piston_work() {
    let mut native = heated_channel();
    let mut walls = [
        Piston::new(0.8, -0.1).unwrap(),
        Piston::new(1.7, 0.15).unwrap(),
    ];
    let original = walls;
    let before = native.flow().totals();
    let control = controls(walls.map(|wall| wall.velocity_m_s()));
    let loads = [InertialWallLoad {
        exterior_pressure_pa: 1.0,
        external_force_n: 0.0,
    }; 2];
    let report = native
        .step_with_inertial_walls(
            0,
            &mut walls,
            Time::new(0.002),
            MolecularStepControl::default(),
            loads,
            control,
            || false,
        )
        .unwrap();
    check(report, control);
    let heat: f64 = report.mechanics.gas.wall_heat_j.into_iter().sum();
    assert!(heat > 1e-5);
    let physical_exchange = native.flow().totals().energy - before.energy
        + (0..2)
            .map(|end| {
                kinetic(walls[end]) - kinetic(original[end])
                    + report.mechanics.exterior_pressure_work_j[end]
                    - report.mechanics.external_force_work_j[end]
            })
            .sum::<f64>();
    close(physical_exchange, heat, 8e-11);
    assert!(physical_exchange.abs() > 1000.0 * control.energy_tolerance_j);
    let gas = report.mechanics.gas;
    close(
        gas.molecular_wall_work_j.into_iter().sum::<f64>() + heat,
        gas.transport.first.boundary_change.energy + gas.transport.second.boundary_change.energy,
        2e-13,
    );
    for end in 0..2 {
        let side = BoundarySide::ALL[end];
        let pressure_work: f64 = native
            .flow()
            .boundary_impulses(side)
            .iter()
            .map(|q| -q.energy)
            .sum();
        let wall_work = pressure_work + gas.molecular_wall_work_j[end];
        close(
            kinetic(walls[end]) - kinetic(original[end])
                + wall_work
                + report.mechanics.exterior_pressure_work_j[end],
            0.0,
            4e-11,
        );
    }
}

#[test]
fn g1_fixed_side_shear_and_piston_guides_close_vector_and_angular_momentum() {
    let mut native = owner(WALLS, 0.04, 0.015, 0.0, |x| {
        [
            1.0,
            0.25 + 0.03 * x[1],
            0.12 + 0.04 * x[2],
            -0.05 + 0.03 * x[1],
            2.0 + 0.4 * x[1] + 0.2 * x[2],
        ]
    });
    for side in BoundarySide::ALL {
        native = native
            .with_wall(
                side,
                MolecularWall3d::NoSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                },
            )
            .unwrap();
    }
    let mut walls = [
        Piston::new(2.0, 0.0).unwrap(),
        Piston::new(5.0, 0.0).unwrap(),
    ];
    let original = walls;
    let before = native.flow().totals();
    let initial_grid = *native.flow().grid();
    let control = controls([0.0; 2]);
    let report = native
        .step_with_inertial_walls(
            0,
            &mut walls,
            Time::new(0.001),
            MolecularStepControl::default(),
            [InertialWallLoad::default(); 2],
            control,
            || false,
        )
        .unwrap();
    check(report, control);
    assert!(report.side_boundary_impulse_n_s[0] < -1e-6);
    assert!(
        report
            .guide_impulse_n_s
            .iter()
            .flatten()
            .any(|q| q.abs() > 1e-6)
    );
    let mut side_impulse = [0.0; 3];
    for side in [
        BoundarySide::YMin,
        BoundarySide::YMax,
        BoundarySide::ZMin,
        BoundarySide::ZMax,
    ] {
        for impulse in native.boundary_impulses(side) {
            for d in 0..3 {
                side_impulse[d] -= impulse.momentum[d];
            }
        }
    }
    for d in 0..3 {
        close(report.side_boundary_impulse_n_s[d], side_impulse[d], 2e-14);
    }
    let after = native.flow().totals();
    let piston_impulse: f64 = (0..2)
        .map(|end| {
            walls[end].mass_kg() * (walls[end].velocity_m_s() - original[end].velocity_m_s())
        })
        .sum();
    close(
        after.momentum[0] - before.momentum[0] + piston_impulse - side_impulse[0],
        0.0,
        4e-11,
    );
    for d in [1, 2] {
        close(
            after.momentum[d]
                - before.momentum[d]
                - side_impulse[d]
                - report.guide_impulse_n_s[0][d]
                - report.guide_impulse_n_s[1][d],
            0.0,
            4e-12,
        );
    }
    for end in 0..2 {
        let side = BoundarySide::ALL[end];
        let mut torque = [0.0; 3];
        let mut inward = [0.0; 3];
        for (face, impulse) in native.boundary_impulses(side).iter().enumerate() {
            let position = initial_grid.face_center(side, face).unwrap();
            let arm = [
                0.0,
                position[1]
                    - initial_grid.origin[1]
                    - 0.5 * initial_grid.shape[1] as f64 * initial_grid.spacing[1],
                position[2]
                    - initial_grid.origin[2]
                    - 0.5 * initial_grid.shape[2] as f64 * initial_grid.spacing[2],
            ];
            let j = impulse.momentum;
            let cross = [arm[1] * j[2] - arm[2] * j[1], arm[2] * j[0], -arm[1] * j[0]];
            for d in 0..3 {
                torque[d] += cross[d];
                inward[d] -= j[d];
            }
        }
        assert!(torque.into_iter().any(|q| q.abs() > 1e-6));
        for d in 0..3 {
            close(report.gas_angular_impulse_n_m_s[end][d], torque[d], 2e-14);
            close(
                report.guide_angular_impulse_n_m_s[end][d],
                -torque[d],
                2e-14,
            );
            close(
                report.guide_impulse_n_s[end][d],
                if d == 0 { 0.0 } else { inward[d] },
                2e-14,
            );
        }
    }
}

#[test]
fn g2_g4_finite_viscous_pistons_retry_whole_intervals_and_cancel_before_publication() {
    let transport = MolecularStepControl {
        molecular: MolecularControl {
            maximum_step_s: 0.001,
            maximum_steps: 1,
            ..MolecularControl::default()
        },
        ..MolecularStepControl::default()
    };
    let initial = [
        Piston::new(0.8, -0.1).unwrap(),
        Piston::new(1.7, 0.15).unwrap(),
    ];
    let control = controls(initial.map(|wall| wall.velocity_m_s()));
    let loads = [InertialWallLoad {
        exterior_pressure_pa: 1.0,
        external_force_n: 0.0,
    }; 2];
    let mut reference = heated_channel();
    let mut reference_walls = initial;
    let mut polls = 0;
    let report = reference
        .step_with_inertial_walls(
            0,
            &mut reference_walls,
            Time::new(0.008),
            transport,
            loads,
            control,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    check(report, control);
    assert!(report.mechanics.gas.transport.coupling_retries >= 2);
    let mut direct = heated_channel();
    let mut direct_walls = initial;
    let direct_report = direct
        .step_with_inertial_walls(
            0,
            &mut direct_walls,
            report.mechanics.gas.transport.elapsed,
            transport,
            loads,
            control,
            || false,
        )
        .unwrap();
    assert_eq!(reference.flow().states(), direct.flow().states());
    assert_eq!(reference.flow().grid(), direct.flow().grid());
    assert_eq!(reference_walls, direct_walls);
    assert_eq!(
        report.mechanics.gas.transport.first,
        direct_report.mechanics.gas.transport.first
    );
    assert_eq!(
        report.mechanics.gas.transport.second,
        direct_report.mechanics.gas.transport.second
    );
    for cutoff in [polls / 2, polls] {
        let mut cancelled = heated_channel();
        let mut walls = initial;
        let state = cancelled.flow().states().to_vec();
        let grid = *cancelled.flow().grid();
        let previous = cancelled.last_step().copied();
        let loads_before = BoundarySide::ALL.map(|side| cancelled.boundary_impulses(side).to_vec());
        let mut calls = 0;
        assert_eq!(
            cancelled.step_with_inertial_walls(
                0,
                &mut walls,
                Time::new(0.008),
                transport,
                loads,
                control,
                || {
                    calls += 1;
                    calls == cutoff
                },
            ),
            Err(MolecularInertialWallsError::Cancelled)
        );
        assert_eq!(walls, initial);
        assert_eq!(cancelled.flow().states(), state);
        assert_eq!(*cancelled.flow().grid(), grid);
        assert_eq!(cancelled.flow().time(), Time::new(0.0));
        assert_eq!(cancelled.last_step().copied(), previous);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(cancelled.boundary_impulses(side), loads_before[i]);
        }
        let replay = cancelled
            .step_with_inertial_walls(
                0,
                &mut walls,
                Time::new(0.008),
                transport,
                loads,
                control,
                || false,
            )
            .unwrap();
        assert_eq!(replay, report);
        assert_eq!(walls, reference_walls);
        assert_eq!(cancelled.flow().states(), reference.flow().states());
    }
}

#[test]
fn g3_axial_boost_preserves_viscous_recoil_heat_and_reciprocal_work() {
    let boost = 0.8;
    let make = |shift: f64| {
        let mut result = owner(channel(0), 0.03, 0.01, 0.05, |x| {
            [
                1.0 + 0.03 * x[0],
                0.1 + 0.04 * x[0] + shift,
                0.07,
                -0.03,
                2.0 + 0.1 * x[0],
            ]
        });
        for side in [BoundarySide::XMin, BoundarySide::XMax] {
            result = result
                .with_wall(
                    side,
                    MolecularWall3d::NavierSlip {
                        velocity: [0.0; 3],
                        thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                        slip_length_m: 0.08,
                        fluid_heat_fraction: 0.4,
                    },
                )
                .unwrap();
        }
        result
    };
    let mut a = make(0.0);
    let mut b = make(boost);
    let mut wa = [
        Piston::new(0.8, -0.1).unwrap(),
        Piston::new(1.7, 0.15).unwrap(),
    ];
    let mut wb = wa.map(|wall| Piston::new(wall.mass_kg(), wall.velocity_m_s() + boost).unwrap());
    let ca = controls(wa.map(|wall| wall.velocity_m_s()));
    let cb = controls(wb.map(|wall| wall.velocity_m_s()));
    let loads = [InertialWallLoad {
        exterior_pressure_pa: 1.0,
        external_force_n: 0.0,
    }; 2];
    let ra = a
        .step_with_inertial_walls(
            0,
            &mut wa,
            Time::new(0.001),
            MolecularStepControl::default(),
            loads,
            ca,
            || false,
        )
        .unwrap();
    let rb = b
        .step_with_inertial_walls(
            0,
            &mut wb,
            Time::new(0.001),
            MolecularStepControl::default(),
            loads,
            cb,
            || false,
        )
        .unwrap();
    check(ra, ca);
    check(rb, cb);
    for end in 0..2 {
        close(
            wb[end].velocity_m_s() - wa[end].velocity_m_s(),
            boost,
            3e-10,
        );
        close(
            rb.mechanics.gas_impulse_n_s[end],
            ra.mechanics.gas_impulse_n_s[end],
            3e-10,
        );
    }
    for (qa, qb) in a.flow().states().iter().zip(b.flow().states()) {
        close(qb.density, qa.density, 4e-11);
        close(qb.momentum[0], qa.momentum[0] + boost * qa.density, 4e-10);
        close(
            qb.energy,
            qa.energy + boost * qa.momentum[0] + 0.5 * boost.powi(2) * qa.density,
            4e-10,
        );
    }
    for side in 0..6 {
        close(
            rb.mechanics.gas.wall_heat_j[side],
            ra.mechanics.gas.wall_heat_j[side],
            3e-11,
        );
    }
    // The same complete boundary impulse transforms laboratory work.
    let ja: f64 = ra.mechanics.gas_impulse_n_s.into_iter().sum();
    let da = ra.mechanics.gas.transport.after.energy - ra.mechanics.gas.transport.before.energy;
    let db = rb.mechanics.gas.transport.after.energy - rb.mechanics.gas.transport.before.energy;
    close(db, da - boost * ja, 4e-10);
}
