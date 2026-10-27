//! Actual regional source state refreshes a continuing native gas/wall solve.
use super::*;
use fs_duct::dg::{DgSidePort, DgWallResolution, PortId};

fn acoustic(field: &AxialGasField) -> DgDuct {
    DgDuct::new_with_gas_field(
        &fs_duct::Duct { segments: vec![fs_duct::Segment::Cylinder {
            radius: 0.001, length: 0.3,
        }] },
        field, DgMesh { elements: 2, order: 1 },
        DgBoundary::Closed, DgBoundary::Closed,
    ).unwrap().with_side_ports(&[DgSidePort { id: PortId(7), downstream_element: 1 }]).unwrap()
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 16, quadrature_points: 2,
            band_rad_s: [20.0, 4000.0], maximum_relative_error: 0.001,
        }).unwrap()
}

#[test]
fn regional_source_refresh_continues_native_diffusion_wall_motion_and_atomic_publication() {
    let initial_field = gas_field(0.8, 350.0).unwrap();
    let replacement = gas_field(0.1, 450.0).unwrap();
    // The inlet is unchanged: this must consume the actual varying field.
    assert_eq!(initial_field.state_at(0.0).unwrap().state(), replacement.state_at(0.0).unwrap().state());
    assert_ne!(initial_field.state_at(0.15).unwrap().state(), replacement.state_at(0.15).unwrap().state());
    assert_eq!(initial_field.model().source_packs(), replacement.model().source_packs());
    let source = source_wall(1800.0, 4e6);
    let right = LocalMaterialWall::new(
        source.material().clone(), 2e-5, source.axial_constraint(),
        source.loss().clone(), source.band_hz(),
    ).unwrap();
    let walls = vec![kelvin_wall(1000.0, 1e6, 2.0, 2e-5), right];
    let mut bound = RegionalMaterialWallDg::new(walls.clone(), acoustic(&initial_field)).unwrap();
    let mut step = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|x| 10.0 + x, |_| 1e-6);
    for _ in 0..12 { step.step(&mut state, &[]).unwrap(); }
    let before = state.clone();
    let old_motion = bound.duct().compliant_wall_state(&state).unwrap();
    assert!(old_motion.relaxation_stored_energy_j > 0.0);
    let old_step = step.clone();
    let mut control = state.clone();
    let report = bound.rebind_spatial_viscothermal_gas_field_controlled(
        &replacement, &mut state, &mut step, 0.001, 1.0, &mut || false,
    ).unwrap();
    assert_eq!(bound.walls(), walls.as_slice());
    assert_eq!(bound.band_hz(), BAND);
    assert_eq!(bound.duct().gas_field(), Some(&replacement));
    assert_eq!(bound.duct().gas_field().unwrap().model().source_packs().len(), 2);
    assert_eq!(bound.duct().side_port_source(PortId(7)).unwrap().gas, *replacement.state_at(0.15).unwrap().state());
    assert_eq!(bound.duct().compliant_wall_state(&state).unwrap(), old_motion);
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(state.cumulative_dissipated_j(), before.cumulative_dissipated_j());
    assert_eq!(step.dt(), old_step.dt());
    assert_eq!(step.tolerance(), old_step.tolerance());
    let mut difference = 0.0;
    let mut signal = 0.0;
    for _ in 0..16 {
        let frame = step.step(&mut state, &[]).unwrap();
        old_step.step(&mut control, &[]).unwrap();
        step.sync_nodal(&mut state).unwrap();
        old_step.sync_nodal(&mut control).unwrap();
        for (&a, &b) in state.pressure().iter().zip(control.pressure()) {
            difference += (a-b).powi(2);
            signal += b*b;
        }
        assert_eq!(frame.ledger.cumulative_parameter_work_j, report.parameter_work_j);
        close(state.stored_energy_j() + state.cumulative_dissipated_j(),
              before.initial_energy_j() + report.parameter_work_j);
    }
    assert!(difference > 1e-7 * signal);

    let accepted = bound.clone();
    let accepted_state = state.clone();
    let accepted_step = step.clone();
    for (error, budget) in [(1e-16, 1.0), (0.001, 0.0)] {
        assert!(bound.rebind_spatial_viscothermal_gas_field_controlled(
            &initial_field, &mut state, &mut step, error, budget, &mut || false,
        ).is_err());
        assert_eq!(state, accepted_state);
        assert_eq!(bound.walls(), accepted.walls());
        assert_eq!(bound.duct().gas_field(), accepted.duct().gas_field());
        assert_eq!(bound.duct().operator_identity(), accepted.duct().operator_identity());
    }
    let mut polls = 0;
    accepted.clone().rebind_spatial_viscothermal_gas_field_controlled(
        &initial_field, &mut accepted_state.clone(), &mut accepted_step.clone(),
        0.001, 1.0, &mut || { polls += 1; false },
    ).unwrap();
    let mut count = 0;
    assert!(matches!(bound.rebind_spatial_viscothermal_gas_field_controlled(
        &initial_field, &mut state, &mut step, 0.001, 1.0,
        &mut || { count += 1; count == polls },
    ), Err(AcousticRealizeError::Cancelled)));
    assert_eq!(state, accepted_state);
    assert_eq!(bound.walls(), accepted.walls());
    assert_eq!(bound.duct().gas_field(), accepted.duct().gas_field());
    let mut foreign = acoustic(&replacement).propagator(1e-6, 1e-13).unwrap();
    assert!(bound.rebind_spatial_viscothermal_gas_field_controlled(
        &initial_field, &mut state, &mut foreign, 0.001, 1.0, &mut || false,
    ).is_err());
    assert_eq!(state, accepted_state);
    let mut expected = accepted_state.clone();
    assert_eq!(step.step(&mut state, &[]).unwrap(), accepted_step.step(&mut expected, &[]).unwrap());

    // The homogeneous wall source wrapper uses the same prepared transaction.
    let mut homogeneous = walls[1].clone().bind_dg(acoustic(&initial_field)).unwrap();
    let mut stepping = homogeneous.duct().propagator(1e-6, 1e-13).unwrap();
    let mut vibration = homogeneous.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..3 { stepping.step(&mut vibration, &[]).unwrap(); }
    homogeneous.rebind_spatial_viscothermal_gas_field_controlled(
        &replacement, &mut vibration, &mut stepping, 0.001, 1.0, &mut || false,
    ).unwrap();
    assert_eq!(homogeneous.wall(), &walls[1]);
    assert_eq!(homogeneous.duct().gas_field(), Some(&replacement));
    assert_eq!(vibration.step_count(), 3);
    stepping.step(&mut vibration, &[]).unwrap();
}
