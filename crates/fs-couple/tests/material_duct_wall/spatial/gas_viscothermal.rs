//! Synthetic complete gas sources refresh the continuing native diffusion owner.
use super::*;
use fs_duct::dg::{DgSidePort, PortId};

fn uniform_field(temperature: f64, fraction: f64, pressure: f64) -> AxialGasField {
    let catalog = gas_field(1.0, 300.0).unwrap();
    let primitive = PrimitiveGasState {
        temperature_k: temperature,
        composition: composition(fraction),
    };
    AxialGasField::new(
        vec![0.0, 0.3], pressure, vec![primitive; 2], catalog.model().clone(),
        QueryPoint::new(), ThermochemicalPhase::FrozenGas,
    ).unwrap()
}
fn acoustic(field: &AxialGasField) -> DgDuct {
    DgDuct::new_with_gas_field(
        &fs_duct::Duct { segments: vec![fs_duct::Segment::Cylinder {
            radius: 0.001, length: 0.3,
        }] },
        field, DgMesh { elements: 2, order: 1 },
        DgBoundary::Closed, DgBoundary::Closed,
    ).unwrap().with_side_ports(&[DgSidePort { id: PortId(7), downstream_element: 1 }]).unwrap()
        .with_viscothermal(16, [50.0, 2000.0], 0.001).unwrap()
}

#[test]
fn source_refresh_keeps_live_diffusion_wall_receipts_and_certified_clock_atomic() {
    let initial_field = uniform_field(300.0, 1.0, 101_325.0);
    let replacement = uniform_field(420.0, 0.2, 150_000.0);
    let old_gas = *initial_field.uniform_state().unwrap().state();
    let new_gas = *replacement.uniform_state().unwrap().state();
    assert_eq!(initial_field.model().source_packs().len(), 2);
    assert_eq!(initial_field.uniform_state().unwrap().source_packs().len(), 1);
    assert_eq!(replacement.uniform_state().unwrap().source_packs().len(), 2);
    assert_eq!(initial_field.model().source_packs(), replacement.model().source_packs());
    let source = source_wall(1800.0, 4e6);
    // Preserve the original source law while explicitly authoring a thin
    // membrane for this narrow radius: thickness/radius=0.02 <= 0.05.
    let right = LocalMaterialWall::new(
        source.material().clone(), 2e-5, source.axial_constraint(),
        source.loss().clone(), source.band_hz(),
    ).unwrap();
    let walls = vec![kelvin_wall(1000.0, 1e6, 2.0, 2e-5), right];
    let mut bound = RegionalMaterialWallDg::new(walls.clone(), acoustic(&initial_field)).unwrap();
    let mut propagation = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|_| 10.0, |_| 1e-6);
    for _ in 0..12 { propagation.step(&mut state, &[]).unwrap(); }
    propagation.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let old_motion = bound.duct().compliant_wall_state(&state).unwrap();
    assert!(old_motion.relaxation_stored_energy_j > 0.0);
    let old_model = bound.duct().viscothermal_model().unwrap().clone();
    let old_propagation = propagation.clone();
    let mut control = state.clone();
    let report = bound.rebind_viscothermal_gas_field_controlled(
        &replacement, &mut state, &mut propagation, 0.001, 1.0, &mut || false,
    ).unwrap();
    assert_eq!(bound.walls(), walls.as_slice());
    assert_eq!(bound.duct().gas_field(), Some(&replacement));
    assert_eq!(bound.duct().side_port_source(PortId(7)).unwrap().gas, new_gas);
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(state.cumulative_dissipated_j(), before.cumulative_dissipated_j());
    assert_eq!(propagation.dt(), old_propagation.dt());
    assert_eq!(propagation.tolerance(), old_propagation.tolerance());
    let new_model = bound.duct().viscothermal_model().unwrap();
    assert_eq!(new_model.band_rad_s(), old_model.band_rad_s());
    assert_eq!(new_model.series_terms().len(), old_model.series_terms().len());
    assert_ne!(new_model.series_terms()[0].rate_per_s, old_model.series_terms()[0].rate_per_s);
    assert_ne!(new_model.shunt_terms()[0].rate_per_s, old_model.shunt_terms()[0].rate_per_s);
    let motion = bound.duct().compliant_wall_state(&state).unwrap();
    close(motion.stored_energy_j, old_motion.stored_energy_j);
    close(motion.relaxation_stored_energy_j, old_motion.relaxation_stored_energy_j);
    let old_bulk = old_gas.density * old_gas.sound_speed.powi(2);
    let new_bulk = new_gas.density * new_gas.sound_speed.powi(2);
    for i in 0..bound.duct().field_dofs() {
        close(state.pressure()[i] / new_bulk, before.pressure()[i] / old_bulk);
        close(new_gas.density * state.volume_flow()[i], old_gas.density * before.volume_flow()[i]);
        close(motion.displacement_m[i], old_motion.displacement_m[i]);
        close(motion.velocity_m_s[i], old_motion.velocity_m_s[i]);
        for (a, b) in motion.relaxation_pressure_pa.iter().zip(&old_motion.relaxation_pressure_pa) {
            close(a[i], b[i]);
        }
    }
    let mut difference = 0.0;
    let mut signal = 0.0;
    for _ in 0..15 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        old_propagation.step(&mut control, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        old_propagation.sync_nodal(&mut control).unwrap();
        for (&a, &b) in state.pressure().iter().zip(control.pressure()) {
            difference += (a - b).powi(2);
            signal += b * b;
        }
        assert_eq!(frame.ledger.cumulative_parameter_work_j, report.parameter_work_j);
        close(state.stored_energy_j() + state.cumulative_dissipated_j(),
              before.initial_energy_j() + report.parameter_work_j);
    }
    assert!(difference > 1e-6 * signal);

    // Failed new-gas admission, nonuniform source and a cancellation AFTER the
    // replacement propagator is prepared must preserve the accepted trajectory.
    let accepted = bound.clone();
    let accepted_state = state.clone();
    let accepted_propagation = propagation.clone();
    for (field, error) in [
        (initial_field.clone(), 1e-16),
        (gas_field(0.2, 350.0).unwrap(), 0.001),
    ] {
        assert!(bound.rebind_viscothermal_gas_field_controlled(
            &field, &mut state, &mut propagation, error, 1.0, &mut || false,
        ).is_err());
        assert_eq!(state, accepted_state);
        assert_eq!(bound.walls(), accepted.walls());
        assert_eq!(bound.duct().gas_field(), accepted.duct().gas_field());
        assert_eq!(bound.duct().operator_identity(), accepted.duct().operator_identity());
    }
    let mut polls = 0;
    accepted.clone().rebind_viscothermal_gas_field_controlled(
        &initial_field, &mut accepted_state.clone(), &mut accepted_propagation.clone(),
        0.001, 1.0, &mut || { polls += 1; false },
    ).unwrap();
    let mut count = 0;
    assert!(matches!(bound.rebind_viscothermal_gas_field_controlled(
        &initial_field, &mut state, &mut propagation, 0.001, 1.0,
        &mut || { count += 1; count == polls },
    ), Err(AcousticRealizeError::Cancelled)));
    assert_eq!(state, accepted_state);
    assert_eq!(bound.walls(), accepted.walls());
    assert_eq!(bound.duct().gas_field(), accepted.duct().gas_field());
    let mut foreign = acoustic(&replacement).propagator(1e-6, 1e-13).unwrap();
    assert!(bound.rebind_viscothermal_gas_field_controlled(
        &initial_field, &mut state, &mut foreign, 0.001, 1.0, &mut || false,
    ).is_err());
    assert_eq!(state, accepted_state);
    let mut expected = accepted_state.clone();
    assert_eq!(propagation.step(&mut state, &[]).unwrap(),
               accepted_propagation.step(&mut expected, &[]).unwrap());
}
