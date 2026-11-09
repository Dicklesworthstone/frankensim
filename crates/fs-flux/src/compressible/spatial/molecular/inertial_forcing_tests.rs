use super::*;
use crate::compressible::spatial::molecular::{
    InertialWallLoad, InertialWallsControl, MolecularInertialWallsError, Piston,
};
use crate::compressible::spatial::{BodyForce3d, ForceControl3d};

fn mechanics(velocity: [f64; 2]) -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: velocity.map(|u| [u - 2.0, u + 2.0]),
        momentum_tolerance_n_s: 5e-11,
        energy_tolerance_j: 5e-11,
        maximum_iterations: 24,
    }
}
fn cap() -> ForceControl3d {
    ForceControl3d { maximum_mach: 10.0 }
}
fn loads() -> [InertialWallLoad; 2] {
    [0.1, -0.04].map(|external_force_n| InertialWallLoad {
        exterior_pressure_pa: 1.0,
        external_force_n,
    })
}
fn pistons() -> [Piston; 2] {
    [
        Piston::new(0.7, -0.08).unwrap(),
        Piston::new(1.4, 0.12).unwrap(),
    ]
}
fn grid() -> CartesianGrid {
    CartesianGrid::new(3, 2, 2, 0.3, 0.4, 0.5).unwrap()
}
fn channel(axis: usize) -> Boundaries3d {
    let mut b = PERIODIC;
    match axis {
        0 => b.x = WALL,
        1 => b.y = WALL,
        _ => b.z = WALL,
    }
    b
}
fn initial(boundaries: Boundaries3d) -> Euler3d {
    field(gas(0.0), grid(), boundaries, |x| {
        [
            1.0 + 0.02 * x[0],
            0.07 + 0.03 * x[1],
            -0.05 + 0.02 * x[2],
            0.09 - 0.03 * x[0],
            2.0 + 0.1 * x[1] + 0.08 * x[2],
        ]
    })
}
fn viscous() -> MolecularEuler3d<ConstantTransport> {
    let flow = initial(WALLS);
    let source = law(flow.gas(), 0.03, 0.01, 0.05);
    let mut owner = MolecularEuler3d::new(flow, source, || false).unwrap();
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        owner =
            owner
                .with_wall(
                    side,
                    MolecularWall3d::NavierSlip {
                        velocity: [0.0; 3],
                        slip_length_m: 0.05,
                        fluid_heat_fraction: 0.4,
                        thermal: WallThermalCondition::Isothermal(Temperature::new(
                            if i % 2 == 0 { 3.0 } else { 1.0 },
                        )),
                    },
                )
                .unwrap();
    }
    owner
}
fn kinetic(wall: Piston) -> f64 {
    0.5 * wall.mass_kg() * wall.velocity_m_s().powi(2)
}

#[test]
fn g1_forced_finite_pistons_keep_body_impulse_external_and_zero_transport_matches_on_all_axes() {
    let a = [0.7, -0.4, 0.2];
    let force = BodyForce3d::UniformAcceleration(a);
    for axis in 0..3 {
        let flow = initial(channel(axis));
        let law = law(flow.gas(), 0.0, 0.0, 0.0);
        let mut molecular = MolecularEuler3d::new(flow, law, || false).unwrap();
        let mut inviscid = initial(channel(axis));
        let mut mw = pistons();
        let original = mw;
        let mut iw = mw;
        let control = mechanics(mw.map(|w| w.velocity_m_s()));
        let mut mwork = molecular.flow().forcing_workspace(usize::MAX).unwrap();
        let mut iwork = inviscid.forcing_workspace(usize::MAX).unwrap();
        let report = molecular
            .step_with_inertial_walls_and_forcing(
                axis,
                &mut mw,
                Time::new(0.001),
                MolecularStepControl::default(),
                loads(),
                control,
                force,
                cap(),
                &mut mwork,
                || false,
            )
            .unwrap();
        let baseline = inviscid
            .step_with_inertial_walls_and_forcing(
                axis,
                &mut iw,
                Time::new(0.001),
                StepControl::default(),
                loads(),
                control,
                force,
                cap(),
                &mut iwork,
                || false,
            )
            .unwrap();
        let gas = report.mechanics.gas.transport;
        assert_eq!(gas.elapsed, Time::new(0.001));
        assert_eq!(mw, iw);
        assert_eq!(molecular.flow().states(), inviscid.states());
        assert_eq!(molecular.flow().grid(), inviscid.grid());
        assert_eq!(gas.advection, baseline.mechanics.gas);
        assert_eq!(
            report.mechanics.gas_impulse_n_s,
            baseline.mechanics.gas_impulse_n_s
        );
        let dt = gas.elapsed.value();
        for d in 0..3 {
            close(
                gas.source_change.impulse[d],
                gas.before.mass * a[d] * dt,
                2e-12,
            );
        }
        let body_dp = (0..2)
            .map(|i| original[i].mass_kg() * (mw[i].velocity_m_s() - original[i].velocity_m_s()))
            .sum::<f64>();
        close(
            gas.after.momentum[axis] - gas.before.momentum[axis] + body_dp
                - gas.source_change.impulse[axis]
                - dt * (loads()[0].external_force_n + loads()[1].external_force_n),
            0.0,
            control.momentum_tolerance_n_s,
        );
        assert!(gas.source_change.impulse[axis].abs() > 1e-5);
        assert!(
            report.mechanics.combined_momentum_residual_n_s.abs() < control.momentum_tolerance_n_s
        );
        assert!(report.mechanics.total_energy_residual_j.abs() < control.energy_tolerance_j);
        for r in gas.conservation_residual.array() {
            close(r, 0.0, 5e-12);
        }
    }
}

#[test]
fn g1_forced_viscous_pistons_retain_heat_side_support_and_full_guide_momentum_accounts() {
    let f = [0.7, -0.4, 0.2];
    for axis in 0..3 {
        let mut flow = viscous();
        let old_grid = *flow.flow().grid();
        let mut walls = pistons();
        let original = walls;
        let control = mechanics(walls.map(|w| w.velocity_m_s()));
        let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
        let report = flow
            .step_with_inertial_walls_and_forcing(
                axis,
                &mut walls,
                Time::new(0.001),
                MolecularStepControl::default(),
                loads(),
                control,
                BodyForce3d::UniformForceDensity(f),
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        let gas = report.mechanics.gas.transport;
        assert_eq!(gas.elapsed, Time::new(0.001));
        let dt = gas.elapsed.value();
        let v0 = old_grid.cell_volume() * flow.flow().states().len() as f64;
        let v1 = flow.flow().grid().cell_volume() * flow.flow().states().len() as f64;
        for d in 0..3 {
            close(
                gas.source_change.impulse[d],
                0.5 * dt * f[d] * (v0 + v1),
                2e-12,
            );
            let body_dp = if d == axis {
                (0..2)
                    .map(|i| {
                        original[i].mass_kg()
                            * (walls[i].velocity_m_s() - original[i].velocity_m_s())
                    })
                    .sum()
            } else {
                0.0
            };
            let external = if d == axis {
                dt * (loads()[0].external_force_n + loads()[1].external_force_n)
            } else {
                0.0
            };
            let guide = report.guide_impulse_n_s.iter().map(|q| q[d]).sum::<f64>();
            close(
                gas.after.momentum[d] - gas.before.momentum[d] + body_dp
                    - gas.source_change.impulse[d]
                    - report.side_boundary_impulse_n_s[d]
                    - guide
                    - external,
                0.0,
                control.momentum_tolerance_n_s,
            );
            let side: f64 = BoundarySide::ALL
                .into_iter()
                .filter(|s| s.axis() != axis)
                .flat_map(|s| flow.boundary_impulses(s))
                .map(|q| -q.momentum[d])
                .sum();
            close(report.side_boundary_impulse_n_s[d], side, 3e-13);
        }
        let heat = report.mechanics.gas.wall_heat_j.iter().sum::<f64>();
        let wall_energy = (0..2)
            .map(|i| {
                kinetic(walls[i]) - kinetic(original[i])
                    + report.mechanics.exterior_pressure_work_j[i]
                    - report.mechanics.external_force_work_j[i]
            })
            .sum::<f64>();
        close(
            gas.after.energy - gas.before.energy
                + (gas.after.mass - gas.before.mass) * flow.flow().gas().energy_reference()
                + wall_energy
                - heat
                - gas.source_change.work,
            0.0,
            control.energy_tolerance_j,
        );
        assert!(
            report
                .mechanics
                .gas
                .wall_heat_j
                .iter()
                .any(|q| q.abs() > 1e-6)
        );
        assert!(
            report
                .mechanics
                .gas
                .molecular_wall_work_j
                .iter()
                .any(|q| q.abs() > 1e-8)
        );
        assert!(gas.source_change.work.abs() > 1e-7);
        assert!(report.side_boundary_impulse_n_s[axis].abs() > 1e-6);
        assert!(
            report
                .guide_impulse_n_s
                .iter()
                .flatten()
                .any(|q| q.abs() > 1e-6)
        );
        assert!(
            report
                .gas_angular_impulse_n_m_s
                .iter()
                .flatten()
                .any(|q| q.abs() > 1e-8)
        );
        for r in gas.conservation_residual.array() {
            close(r, 0.0, 5e-12);
        }
    }
}

#[test]
fn g0_forced_finite_piston_zero_source_is_the_original_native_transaction() {
    let zero = BodyForce3d::UniformAcceleration([0.0; 3]);
    let mut a = viscous();
    let mut b = viscous();
    let mut aw = pistons();
    let mut bw = aw;
    let control = mechanics(aw.map(|w| w.velocity_m_s()));
    let mut work = b.flow().forcing_workspace(usize::MAX).unwrap();
    let expected = a
        .step_with_inertial_walls(
            0,
            &mut aw,
            Time::new(0.001),
            MolecularStepControl::default(),
            loads(),
            control,
            || false,
        )
        .unwrap();
    let actual = b
        .step_with_inertial_walls_and_forcing(
            0,
            &mut bw,
            Time::new(0.001),
            MolecularStepControl::default(),
            loads(),
            control,
            zero,
            ForceControl3d { maximum_mach: 1e-6 },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(expected, actual);
    assert_eq!(aw, bw);
    assert_eq!(a.flow().states(), b.flow().states());
    let mut a = initial(channel(0));
    let mut b = initial(channel(0));
    let mut aw = pistons();
    let mut bw = aw;
    let mut work = b.forcing_workspace(usize::MAX).unwrap();
    let expected = a
        .step_with_inertial_walls(
            0,
            &mut aw,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            control,
            || false,
        )
        .unwrap();
    let actual = b
        .step_with_inertial_walls_and_forcing(
            0,
            &mut bw,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            control,
            zero,
            ForceControl3d { maximum_mach: 1e-6 },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(expected, actual);
    assert_eq!(aw, bw);
    assert_eq!(a.states(), b.states());
}

#[test]
fn g4_forced_molecular_piston_final_cancel_and_root_exhaustion_preserve_all_owners() {
    let force = BodyForce3d::UniformForceDensity([0.7, -0.4, 0.2]);
    let make = || {
        let mut flow = viscous();
        let mut walls = pistons();
        let control = mechanics(walls.map(|w| w.velocity_m_s()));
        let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
        flow.step_with_inertial_walls_and_forcing(
            0,
            &mut walls,
            Time::new(0.001),
            MolecularStepControl::default(),
            loads(),
            control,
            force,
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
        (flow, walls, work)
    };
    let (mut probe, mut pw, mut work) = make();
    let control = mechanics(pw.map(|w| w.velocity_m_s()));
    let mut polls = 0;
    let expected = probe
        .step_with_inertial_walls_and_forcing(
            0,
            &mut pw,
            Time::new(0.001),
            MolecularStepControl::default(),
            loads(),
            control,
            force,
            cap(),
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for failure in 0..3 {
        let (mut flow, mut walls, mut work) = make();
        let state = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let time = flow.flow().time();
        let previous = flow.last_step().copied();
        let old = walls;
        let boundary = BoundarySide::ALL.map(|s| flow.boundary_impulses(s).to_vec());
        let mut calls = 0;
        let (transport, mechanics) = if failure == 2 {
            (
                MolecularStepControl {
                    maximum_coupling_retries: 0,
                    ..MolecularStepControl::default()
                },
                InertialWallsControl {
                    maximum_iterations: 1,
                    ..control
                },
            )
        } else {
            (MolecularStepControl::default(), control)
        };
        let result = flow.step_with_inertial_walls_and_forcing(
            0,
            &mut walls,
            Time::new(0.001),
            transport,
            loads(),
            mechanics,
            force,
            cap(),
            &mut work,
            || {
                calls += 1;
                failure < 2 && calls == if failure == 0 { polls / 2 } else { polls }
            },
        );
        assert!(result.is_err());
        if failure < 2 {
            assert_eq!(result, Err(MolecularInertialWallsError::Cancelled));
        }
        assert_eq!(flow.flow().states(), state);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), previous);
        assert_eq!(walls, old);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.boundary_impulses(side), boundary[i]);
        }
        let accepted = flow
            .step_with_inertial_walls_and_forcing(
                0,
                &mut walls,
                Time::new(0.001),
                MolecularStepControl::default(),
                loads(),
                control,
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(accepted, expected);
        assert_eq!(walls, pw);
        assert_eq!(flow.flow().states(), probe.flow().states());
    }
}

#[test]
fn g4_forced_inviscid_piston_final_poll_and_mach_refusal_are_atomic() {
    let force = BodyForce3d::UniformAcceleration([0.7, -0.4, 0.2]);
    let make = || initial(channel(1));
    let mut probe = make();
    let mut pw = pistons();
    let control = mechanics(pw.map(|w| w.velocity_m_s()));
    let mut work = probe.forcing_workspace(usize::MAX).unwrap();
    let mut polls = 0;
    let expected = probe
        .step_with_inertial_walls_and_forcing(
            1,
            &mut pw,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            control,
            force,
            cap(),
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cancel in [false, true] {
        let mut flow = make();
        let mut walls = pistons();
        let states = flow.states().to_vec();
        let old_grid = *flow.grid();
        let original = walls;
        let mut work = flow.forcing_workspace(usize::MAX).unwrap();
        let mut calls = 0;
        let result = flow.step_with_inertial_walls_and_forcing(
            1,
            &mut walls,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            control,
            force,
            if cancel {
                cap()
            } else {
                ForceControl3d { maximum_mach: 1e-6 }
            },
            &mut work,
            || {
                calls += 1;
                cancel && calls == polls
            },
        );
        assert!(result.is_err());
        assert_eq!(flow.states(), states);
        assert_eq!(*flow.grid(), old_grid);
        assert_eq!(flow.time(), Time::new(0.0));
        assert!(flow.last_step().is_none());
        assert_eq!(walls, original);
        let accepted = flow
            .step_with_inertial_walls_and_forcing(
                1,
                &mut walls,
                Time::new(0.001),
                StepControl::default(),
                loads(),
                control,
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(accepted, expected);
        assert_eq!(walls, pw);
        assert_eq!(flow.states(), probe.states());
    }
}
