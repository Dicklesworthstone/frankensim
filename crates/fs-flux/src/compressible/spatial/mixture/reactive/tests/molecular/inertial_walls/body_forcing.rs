use super::*;
use crate::compressible::spatial::mixture::{BodyForce3d, ForceControl3d, ForcingWorkspace3d};

fn work(flow: &ReactingEuler3d) -> ForcingWorkspace3d {
    flow.flow()
        .forcing_workspace(flow.flow().required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

#[test]
fn g4_reacting_piston_admits_initial_mach_before_any_molecular_or_mechanical_proposal() {
    let (mut flow, mut op) = viscous_case();
    let mut scratch = work(&flow);
    let bytes = budget(&flow, &op) + scratch.workspace_bytes();
    let mut pistons = walls([-0.02, 0.03]);
    let original_walls = pistons;
    let state = flow.flow().states().to_vec();
    let grid = *flow.flow().grid();
    let time = flow.flow().time();
    let boundary = loads(&flow);
    let molecular_boundary = molecular_loads(&op);
    let original_laws = BoundarySide::ALL.map(|side| op.wall(side));
    let p = flow
        .flow()
        .cell(0)
        .unwrap()
        .primitive(flow.flow().gas(), || false)
        .unwrap();
    let initial_mach = p.velocity.iter().map(|u| u * u).sum::<f64>().sqrt() / p.sound_speed;
    assert!(initial_mach > 0.001);
    let error = flow.step_with_dissipation_forcing_and_inertial_walls(
        &mut op,
        0,
        &mut pistons,
        Time::new(0.001),
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        DissipationControl::default(),
        bytes,
        [InertialWallLoad::default(); 2],
        controls(),
        BodyForce3d::UniformAcceleration([0.4, -0.3, 0.2]),
        ForceControl3d {
            maximum_mach: 0.001,
        },
        &mut scratch,
        || false,
    );
    // Initial physical admission must fail directly, before the first D half
    // can change this state and before a root proposal consumes retry budget.
    assert!(matches!(
        error,
        Err(ReactiveInertialWallsError::Transport(_))
    ));
    assert_eq!(flow.flow().states(), state);
    assert_eq!(*flow.flow().grid(), grid);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(pistons, original_walls);
    assert_eq!(flow.last_step(), None);
    assert_eq!(op.last_step(), None);
    assert_eq!(loads(&flow), boundary);
    assert_eq!(molecular_loads(&op), molecular_boundary);
    assert_eq!(BoundarySide::ALL.map(|side| op.wall(side)), original_laws);
}

fn force_control() -> ForceControl3d {
    ForceControl3d { maximum_mach: 4.0 }
}

// Compute the assembly balance independently from the accepted gas states,
// actual piston velocities, prescribed loads and source-owned gas impulse/work.
fn assembly_balance(
    report: &ReactiveInertialWallsStep,
    gas: &MixtureGas,
    before: [Piston; 2],
    after: [Piston; 2],
    loads: [InertialWallLoad; 2],
    face_area: f64,
) {
    let m = &report.mechanics;
    let r = &m.gas.reactive;
    balance(r);
    let dt = r.elapsed.value();
    let momentum: f64 = (0..2)
        .map(|i| after[i].mass_kg() * (after[i].velocity_m_s() - before[i].velocity_m_s()))
        .sum();
    let applied_impulse = dt
        * (loads[0].external_force_n
            + loads[1].external_force_n
            + face_area * (loads[0].exterior_pressure_pa - loads[1].exterior_pressure_pa));
    close(
        r.after.momentum[report.axis] - r.before.momentum[report.axis] + momentum,
        r.source_change.impulse[report.axis]
            + report.side_boundary_impulse_n_s[report.axis]
            + applied_impulse,
        5e-10,
    );
    for axis in 0..3 {
        if axis != report.axis {
            close(
                r.after.momentum[axis] - r.before.momentum[axis],
                r.source_change.impulse[axis]
                    + report.side_boundary_impulse_n_s[axis]
                    + report.guide_impulse_n_s[0][axis]
                    + report.guide_impulse_n_s[1][axis],
                3e-11,
            );
        }
    }
    let mut body_energy = 0.0;
    for i in 0..2 {
        let v0 = before[i].velocity_m_s();
        let v1 = after[i].velocity_m_s();
        let kinetic = 0.5 * after[i].mass_kg() * (v1 * v1 - v0 * v0);
        let displacement = 0.5 * (v0 + v1) * dt;
        let exterior_work = if i == 0 { -1.0 } else { 1.0 }
            * loads[i].exterior_pressure_pa
            * face_area
            * displacement;
        let applied_work = loads[i].external_force_n * displacement;
        close(m.kinetic_energy_change_j[i], kinetic, 3e-13);
        close(m.exterior_pressure_work_j[i], exterior_work, 3e-13);
        close(m.external_force_work_j[i], applied_work, 3e-13);
        body_energy += kinetic + exterior_work - applied_work;
    }
    let chemical: f64 = (0..gas.species_count())
        .map(|k| {
            (r.after.species_masses[k] - r.before.species_masses[k]) * gas.energy_references()[k]
        })
        .sum();
    let heat: f64 = m.gas.wall_heat_j.iter().sum();
    let stored_change = r.after.energy - r.before.energy;
    assert!(
        (stored_change + chemical + body_energy - heat - r.source_change.work).abs()
            <= controls().energy_tolerance_j
    );
    assert!(
        (stored_change - r.chemistry.referenced_energy_change_j + body_energy
            - heat
            - r.source_change.work)
            .abs()
            <= controls().energy_tolerance_j
    );
    for residual in m
        .momentum_residual_n_s
        .into_iter()
        .chain([m.combined_momentum_residual_n_s])
    {
        assert!(residual.abs() <= controls().momentum_tolerance_n_s);
    }
    for residual in m
        .mechanical_energy_residual_j
        .into_iter()
        .chain([m.total_energy_residual_j])
    {
        assert!(residual.abs() <= controls().energy_tolerance_j);
    }
}

#[test]
fn g1_reacting_pistons_accept_gas_body_impulse_and_work_on_every_axis() {
    let acceleration = [0.7, -0.4, 0.25];
    let area = [0.8, 0.6, 0.48];
    for axis in 0..3 {
        let mut flow = case(axis, 16.0, 12.0);
        let mut scratch = work(&flow);
        let mut pistons = walls([-0.02, 0.03]);
        let before = pistons;
        let external = [InertialWallLoad::default(); 2];
        let report = flow
            .step_with_inertial_walls_and_forcing(
                axis,
                &mut pistons,
                Time::new(0.001),
                ReactiveStepControl::default(),
                external,
                controls(),
                BodyForce3d::UniformAcceleration(acceleration),
                force_control(),
                &mut scratch,
                || false,
            )
            .unwrap();
        assembly_balance(
            &report,
            flow.flow().gas(),
            before,
            pistons,
            external,
            area[axis],
        );
        let r = &report.mechanics.gas.reactive;
        assert!(r.chemistry.species_change_kg[1] > 0.0);
        assert!(r.source_change.work.abs() > 1e-9);
        for a in 0..3 {
            close(
                r.source_change.impulse[a],
                r.before.mass() * acceleration[a] * r.elapsed.value(),
                3e-13,
            );
        }
        assert_eq!(flow.last_step(), Some(r));
    }
}

#[test]
fn g2_reacting_viscous_pistons_keep_force_density_heat_guides_and_external_loads_distinct() {
    let (mut flow, mut op) = viscous_case();
    let mut scratch = work(&flow);
    let bytes = budget(&flow, &op) + scratch.workspace_bytes();
    let mut pistons = walls([-0.02, 0.03]);
    let before = pistons;
    let volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
    let density = [0.6, -0.2, 0.3];
    let external = [
        InertialWallLoad {
            exterior_pressure_pa: 0.7,
            external_force_n: 0.2,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.7,
            external_force_n: -0.1,
        },
    ];
    let report = flow
        .step_with_dissipation_forcing_and_inertial_walls(
            &mut op,
            0,
            &mut pistons,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            external,
            controls(),
            BodyForce3d::UniformForceDensity(density),
            force_control(),
            &mut scratch,
            || false,
        )
        .unwrap();
    assembly_balance(&report, flow.flow().gas(), before, pistons, external, 0.8);
    let gas = &report.mechanics.gas;
    let r = &gas.reactive;
    let next_volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
    for a in 0..3 {
        close(
            r.source_change.impulse[a],
            density[a] * r.elapsed.value() * 0.5 * (volume + next_volume),
            3e-13,
        );
    }
    assert!(r.source_change.work.abs() > 1e-8);
    assert!(r.chemistry.species_change_kg[1] > 0.0);
    assert!(gas.wall_heat_j.iter().sum::<f64>() > 1e-5);
    assert!(report.side_boundary_impulse_n_s[0].abs() > 1e-6);
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
            .any(|q| q.abs() > 1e-6)
    );
    let split = r.dissipation.unwrap();
    close(
        gas.wall_heat_j.iter().sum::<f64>() + gas.molecular_wall_work_j.iter().sum::<f64>(),
        split.first.boundary_change.energy + split.second.boundary_change.energy,
        3e-12,
    );
    assert_eq!(op.last_step(), Some(&split.second));
    assert_eq!(budget(&flow, &op) + scratch.workspace_bytes(), bytes);
}

#[test]
fn g0_zero_body_source_preserves_original_reacting_piston_transactions() {
    for molecular in [false, true] {
        let (mut expected, mut expected_op) = viscous_case();
        let (mut flow, mut op) = viscous_case();
        let mut reference_walls = walls([-0.02, 0.03]);
        let mut pistons = reference_walls;
        let mut scratch = work(&flow);
        let bytes = budget(&flow, &op);
        let reference = if molecular {
            expected
                .step_with_dissipation_and_inertial_walls(
                    &mut expected_op,
                    0,
                    &mut reference_walls,
                    Time::new(0.001),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    [InertialWallLoad::default(); 2],
                    controls(),
                    || false,
                )
                .unwrap()
        } else {
            expected
                .step_with_inertial_walls(
                    0,
                    &mut reference_walls,
                    Time::new(0.001),
                    ReactiveStepControl::default(),
                    [InertialWallLoad::default(); 2],
                    controls(),
                    || false,
                )
                .unwrap()
        };
        let actual = if molecular {
            flow.step_with_dissipation_forcing_and_inertial_walls(
                &mut op,
                0,
                &mut pistons,
                Time::new(0.001),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes + scratch.workspace_bytes(),
                [InertialWallLoad::default(); 2],
                controls(),
                BodyForce3d::UniformAcceleration([0.0; 3]),
                force_control(),
                &mut scratch,
                || false,
            )
            .unwrap()
        } else {
            flow.step_with_inertial_walls_and_forcing(
                0,
                &mut pistons,
                Time::new(0.001),
                ReactiveStepControl::default(),
                [InertialWallLoad::default(); 2],
                controls(),
                BodyForce3d::UniformAcceleration([0.0; 3]),
                force_control(),
                &mut scratch,
                || false,
            )
            .unwrap()
        };
        assert_eq!(actual, reference);
        assert_eq!(pistons, reference_walls);
        assert_eq!(flow.flow().states(), expected.flow().states());
        assert_eq!(flow.flow().grid(), expected.flow().grid());
        assert_eq!(loads(&flow), loads(&expected));
        assert_eq!(molecular_loads(&op), molecular_loads(&expected_op));
        assert_eq!(op.last_step(), expected_op.last_step());
    }
}

fn forced_viscous(
    flow: &mut ReactingEuler3d,
    op: &mut MolecularTransport3d<ConstantTransport>,
    pistons: &mut [Piston; 2],
    scratch: &mut ForcingWorkspace3d,
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveInertialWallsStep, ReactiveInertialWallsError> {
    let bytes = budget(flow, op) + scratch.workspace_bytes();
    flow.step_with_dissipation_forcing_and_inertial_walls(
        op,
        0,
        pistons,
        Time::new(0.001),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        bytes,
        [InertialWallLoad::default(); 2],
        controls(),
        BodyForce3d::UniformAcceleration([0.4, -0.3, 0.2]),
        force_control(),
        scratch,
        cancelled,
    )
}

#[test]
fn g4_forced_reacting_piston_final_checkpoint_replays_and_failed_root_publishes_nothing() {
    let make = || {
        let (mut flow, mut op) = viscous_case();
        let mut pistons = walls([-0.02, 0.03]);
        let mut scratch = work(&flow);
        forced_viscous(&mut flow, &mut op, &mut pistons, &mut scratch, || false).unwrap();
        (flow, op, pistons, scratch)
    };
    let (mut reference, mut reference_op, mut reference_walls, mut reference_work) = make();
    let mut polls = 0;
    let expected = forced_viscous(
        &mut reference,
        &mut reference_op,
        &mut reference_walls,
        &mut reference_work,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let (mut flow, mut op, mut pistons, mut scratch) = make();
    let state = flow.flow().states().to_vec();
    let grid = *flow.flow().grid();
    let time = flow.flow().time();
    let previous = flow.last_step().copied();
    let operator_report = op.last_step().copied();
    let original_walls = pistons;
    let wall_data = BoundarySide::ALL.map(|side| op.wall(side));
    let original_loads = loads(&flow);
    let original_molecular = molecular_loads(&op);
    let mut calls = 0;
    assert_eq!(
        forced_viscous(&mut flow, &mut op, &mut pistons, &mut scratch, || {
            calls += 1;
            calls == polls
        },),
        Err(ReactiveInertialWallsError::Cancelled)
    );
    assert_eq!(calls, polls);
    assert_eq!(flow.flow().states(), state);
    assert_eq!(*flow.flow().grid(), grid);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step().copied(), previous);
    assert_eq!(op.last_step().copied(), operator_report);
    assert_eq!(pistons, original_walls);
    assert_eq!(loads(&flow), original_loads);
    assert_eq!(molecular_loads(&op), original_molecular);
    assert_eq!(BoundarySide::ALL.map(|side| op.wall(side)), wall_data);
    let recovered =
        forced_viscous(&mut flow, &mut op, &mut pistons, &mut scratch, || false).unwrap();
    assert_eq!(recovered, expected);
    assert_eq!(pistons, reference_walls);
    assert_eq!(flow.flow().states(), reference.flow().states());

    let (mut refused, mut refused_op) = viscous_case();
    let mut refused_work = work(&refused);
    let bytes = budget(&refused, &refused_op) + refused_work.workspace_bytes();
    let mut refused_walls = walls([0.0; 2]);
    let original = refused_walls;
    let state = refused.flow().states().to_vec();
    let grid = *refused.flow().grid();
    let boundary = loads(&refused);
    let molecular_boundary = molecular_loads(&refused_op);
    let error = refused.step_with_dissipation_forcing_and_inertial_walls(
        &mut refused_op,
        0,
        &mut refused_walls,
        Time::new(0.001),
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        DissipationControl::default(),
        bytes,
        [InertialWallLoad::default(); 2],
        InertialWallsControl {
            midpoint_velocity_bracket_m_s: [[-1e-9, 1e-9]; 2],
            maximum_iterations: 1,
            ..controls()
        },
        BodyForce3d::UniformAcceleration([0.4, -0.3, 0.2]),
        force_control(),
        &mut refused_work,
        || false,
    );
    assert!(matches!(
        error,
        Err(ReactiveInertialWallsError::StepRetries { attempts: 1, .. })
    ));
    assert_eq!(refused.flow().states(), state);
    assert_eq!(*refused.flow().grid(), grid);
    assert_eq!(refused_walls, original);
    assert_eq!(refused.last_step(), None);
    assert_eq!(refused_op.last_step(), None);
    assert_eq!(loads(&refused), boundary);
    assert_eq!(molecular_loads(&refused_op), molecular_boundary);
}
