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
