//! Independent volumetric heat/shear and wall-conservation fixtures.
use super::*;

#[test]
fn rkl_volumetric_heat_resolves_every_axis_and_retains_source_cv() {
    let n = 8;
    let lengths: [f64; 3] = [1.0, 1.5, 0.75];
    let h = lengths.map(|l| l / n as f64);
    let k = lengths.map(|l| 2.0 * std::f64::consts::PI / l);
    let mode = |x: [f64; 3]| (k[0] * x[0]).sin() * (k[1] * x[1]).cos() * (k[2] * x[2]).sin();
    let end = 0.2;
    let lambda = 4.0 * (0.05 / 2.5) * (std::f64::consts::PI / n as f64).sin().powi(2)
        * h.iter().map(|v| v.recip().powi(2)).sum::<f64>();
    let expected = 0.2 * (-lambda * end).exp();
    let mut errors = Vec::new();
    for step in [0.2, 0.1] {
        let mut flow = field(
            gas(0.0), CartesianGrid::new(n, n, n, h[0], h[1], h[2]).unwrap(),
            PERIODIC, |x| [1.0, 0.0, 0.0, 0.0, 2.0 + 0.2 * mode(x)],
        );
        let mut op = MolecularTransport3d::new(
            &flow, law(flow.gas(), 0.0, 0.0, 0.05), || false,
        ).unwrap();
        let report = op.advance_stabilized(
            &mut flow, Time::new(end),
            MolecularControl { maximum_step_s: step, ..MolecularControl::default() },
            MolecularStabilization::default(), || false,
        ).unwrap();
        let mut projection = 0.0;
        let mut norm = 0.0;
        for (i, state) in flow.states().iter().enumerate() {
            let phi = mode(flow.grid().cell_center(i).unwrap());
            let p = state.primitive(flow.gas(), || false).unwrap();
            projection += (p.temperature - 2.0) * phi;
            norm += phi * phi;
            assert_eq!(state.density, 1.0);
            assert_eq!(state.momentum, [0.0; 3]);
        }
        errors.push((projection / norm - expected).abs());
        assert!(report.maximum_stages_used > 2);
        assert!(report.transport.maximum_parabolic_number > 0.4);
        assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
        for residual in report.transport.conservation_residual.array() {
            close(residual, 0.0, 4e-12);
        }
        assert_eq!(flow.time().value(), 0.0);
    }
    assert!(errors[0] < 0.002, "{errors:?}");
    assert!(errors[1] < 0.35 * errors[0], "{errors:?}");
}

#[test]
fn rkl_shear_decay_turns_lost_kinetic_energy_into_original_internal_energy() {
    let n = 10;
    let end = 0.1;
    let wave = 2.0 * std::f64::consts::PI;
    let mut flow = cube(n, |x| [1.0, 0.0, 0.1 * (wave * x[0]).sin(), 0.0, 2.0]);
    let (kinetic0, internal0) = partition(&flow);
    let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), 0.05, 0.0, 0.0), || false).unwrap();
    let report = op.advance_stabilized(
        &mut flow, Time::new(end), MolecularControl::default(),
        MolecularStabilization::default(), || false,
    ).unwrap();
    let lambda = 4.0 * 0.05 * (n as f64).powi(2)
        * (std::f64::consts::PI / n as f64).sin().powi(2);
    let amplitude = 0.1 * (-lambda * end).exp();
    for (i, q) in flow.states().iter().enumerate() {
        let expected = amplitude * (wave * flow.grid().cell_center(i).unwrap()[0]).sin();
        close(q.momentum[1], expected, 8e-5);
        close(q.momentum[0], 0.0, 1e-13);
        close(q.momentum[2], 0.0, 1e-13);
    }
    let (kinetic1, internal1) = partition(&flow);
    assert!(kinetic1 < kinetic0);
    assert!(internal1 > internal0);
    close(kinetic0 - kinetic1, internal1 - internal0, 4e-12);
    assert!(report.maximum_stages_used > 2);
    for residual in report.transport.conservation_residual.array() {
        close(residual, 0.0, 4e-12);
    }
}

#[test]
fn rkl_heated_walls_use_the_same_signed_recurrence_and_keep_last_accepted_loads() {
    let mut flow = field(
        gas(0.0), CartesianGrid::new(4, 3, 2, 0.25, 0.3, 0.4).unwrap(),
        WALLS, |_| [1.0, 0.0, 0.0, 0.0, 2.0],
    );
    let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), 0.0, 0.0, 0.5), || false).unwrap();
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        op = op.with_wall(side, MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(fs_qty::Temperature::new(3.0 + 0.1 * i as f64)),
        }).unwrap();
    }
    let selected = MolecularStabilization::default();
    let report = op.advance_stabilized(
        &mut flow, Time::new(0.2), MolecularControl::default(), selected, || false,
    ).unwrap();
    assert!(report.maximum_stages_used > 2);
    assert!(report.transport.after.energy > report.transport.before.energy);
    let outward: f64 = BoundarySide::ALL.into_iter()
        .flat_map(|s| op.boundary_impulses(s)).map(|q| q.energy).sum();
    close(outward + report.transport.boundary_change.energy, 0.0, 2e-12);
    for residual in report.transport.conservation_residual.array() {
        close(residual, 0.0, 4e-12);
    }
    let accepted = flow.states().to_vec();
    let loads = BoundarySide::ALL.map(|s| op.boundary_impulses(s).to_vec());
    assert!(matches!(
        op.advance_stabilized(
            &mut flow, Time::new(0.2), MolecularControl::default(),
            MolecularStabilization { maximum_operator_evaluations: 2, ..selected }, || false,
        ),
        Err(MolecularError::Law(
            fs_gas::multicomponent::dissipation::DissipationError::OperatorWorkLimit { .. }
        ))
    ));
    assert_eq!(flow.states(), accepted);
    assert_eq!(flow.time().value(), 0.0);
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        assert_eq!(op.boundary_impulses(side), loads[i]);
    }
}

#[test]
fn rkl_forced_spatial_retries_keep_spent_work_and_replay_one_physical_state() {
    use crate::compressible::spatial::{BodyForce3d, ForceControl3d};
    let force = BodyForce3d::UniformAcceleration([20.0, 0.0, 0.0]);
    let cap = ForceControl3d { maximum_mach: 10.0 };
    let controls = MolecularStepControl::default();
    let selected = MolecularStabilization::default();
    let make = || {
        let initial = cube(3, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
        let source = law(initial.gas(), 0.0, 0.0, 0.0);
        let mut flow = MolecularEuler3d::new(initial, source, || false).unwrap();
        let mut work = flow.flow().forcing_workspace(40 * 27).unwrap();
        flow.step_stabilized_with_forcing(
            Time::new(0.001), controls, selected, force, cap, &mut work, || false,
        ).unwrap();
        (flow, work)
    };
    let (mut probe, mut work) = make();
    let before = probe.flow().totals();
    let mut polls = 0;
    let report = probe.step_stabilized_with_forcing(
        Time::new(1.0), controls, selected, force, cap, &mut work,
        || { polls += 1; false },
    ).unwrap();
    let physical = report.transport;
    assert!(physical.coupling_retries > 0);
    assert_eq!(physical.first.elapsed.value() * 2.0, physical.elapsed.value());
    assert_eq!(physical.second.elapsed, physical.first.elapsed);
    assert_eq!(physical.advection.elapsed, physical.elapsed);
    let speed0 = before.momentum[0] / before.mass;
    let speed1 = speed0 + 20.0 * physical.elapsed.value();
    close(physical.source_change.impulse[0], before.mass * (speed1 - speed0), 3e-13);
    close(physical.source_change.work, 0.5 * before.mass * (speed1 * speed1 - speed0 * speed0), 3e-13);
    for residual in physical.conservation_residual.array() {
        close(residual, 0.0, 4e-12);
    }
    let (mut direct, mut direct_work) = make();
    let replay = direct.step_stabilized_with_forcing(
        physical.elapsed, controls, selected, force, cap, &mut direct_work, || false,
    ).unwrap();
    assert_eq!(replay.transport.coupling_retries, 0);
    assert_eq!(replay.transport.first, physical.first);
    assert_eq!(replay.transport.advection, physical.advection);
    assert_eq!(replay.transport.second, physical.second);
    assert_eq!(direct.flow().states(), probe.flow().states());
    assert!(report.operator_evaluations > replay.operator_evaluations);

    // A cap sufficient for the final duration cannot pay for discarded probes.
    let (mut limited, mut limited_work) = make();
    let states = limited.flow().states().to_vec();
    let time = limited.flow().time();
    let previous = limited.last_step().copied();
    let loads = BoundarySide::ALL.map(|s| limited.boundary_impulses(s).to_vec());
    assert!(matches!(
        limited.step_stabilized_with_forcing(
            Time::new(1.0), controls,
            MolecularStabilization { maximum_operator_evaluations: replay.operator_evaluations, ..selected },
            force, cap, &mut limited_work, || false,
        ),
        Err(MolecularError::Law(
            fs_gas::multicomponent::dissipation::DissipationError::OperatorWorkLimit { .. }
        ))
    ));
    assert_eq!(limited.flow().states(), states);
    assert_eq!(limited.flow().time(), time);
    assert_eq!(limited.last_step().copied(), previous);
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        assert_eq!(limited.boundary_impulses(side), loads[i]);
    }
    let (mut interrupted, mut interrupted_work) = make();
    let mut count = 0;
    assert_eq!(
        interrupted.step_stabilized_with_forcing(
            Time::new(1.0), controls, selected, force, cap, &mut interrupted_work,
            || { count += 1; count == polls },
        ),
        Err(MolecularError::Cancelled)
    );
    assert_eq!(interrupted.flow().states(), states);
    assert_eq!(interrupted.flow().time(), time);
    assert_eq!(interrupted.last_step().copied(), previous);
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        assert_eq!(interrupted.boundary_impulses(side), loads[i]);
    }
    let completed = interrupted.step_stabilized_with_forcing(
        Time::new(1.0), controls, selected, force, cap, &mut interrupted_work, || false,
    ).unwrap();
    assert_eq!(completed, report);
    assert_eq!(interrupted.flow().states(), probe.flow().states());
}

#[test]
fn rkl_spatial_coupling_closes_heated_wall_energy_across_both_halves() {
    let initial = field(
        gas(0.0), CartesianGrid::new(4, 3, 2, 0.25, 0.3, 0.4).unwrap(),
        WALLS, |_| [1.0, 0.0, 0.0, 0.0, 2.0],
    );
    let source = law(initial.gas(), 0.01, 0.003, 5.0);
    let mut flow = MolecularEuler3d::new(initial, source, || false).unwrap();
    for side in BoundarySide::ALL {
        flow = flow.with_wall(side, MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(fs_qty::Temperature::new(3.0)),
        }).unwrap();
    }
    let report = flow.step_stabilized(
        Time::new(0.005), MolecularStepControl::default(),
        MolecularStabilization::default(), || false,
    ).unwrap();
    let physical = report.transport;
    assert!(report.maximum_stages_used > 2);
    assert!(physical.after.energy > physical.before.energy);
    assert!(physical.first.boundary_change.energy > 0.0);
    assert!(physical.second.boundary_change.energy > 0.0);
    let outward: f64 = BoundarySide::ALL.into_iter()
        .flat_map(|s| flow.boundary_impulses(s)).map(|q| q.energy).sum();
    close(outward + physical.boundary_change.energy, 0.0, 3e-12);
    close(physical.after.energy - physical.before.energy, physical.boundary_change.energy, 4e-12);
    assert_eq!(flow.last_step().copied(), Some(physical));
}
