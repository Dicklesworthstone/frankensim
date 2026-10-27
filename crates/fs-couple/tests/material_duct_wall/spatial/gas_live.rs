//! A new shared gas source changes the continuing native duct, not a new voice.
use super::*;
use fs_duct::dg::{DgSidePort, PortId};

fn uniform_field(temperature: f64, fraction: f64, pressure: f64, length: f64) -> AxialGasField {
    let catalog = gas_field(1.0, 300.0).unwrap();
    let primitive = PrimitiveGasState {
        temperature_k: temperature,
        composition: composition(fraction),
    };
    AxialGasField::new(
        vec![0.0, length],
        pressure,
        vec![primitive; 2],
        catalog.model().clone(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}
fn bulk(gas: &GasState) -> f64 {
    gas.density * gas.sound_speed.powi(2)
}

#[test]
fn source_gas_refresh_reaches_live_pressure_flow_ports_and_regional_wall_memory() {
    let initial_field = uniform_field(300.0, 1.0, 101_325.0, 0.3);
    let replacement = uniform_field(420.0, 0.2, 150_000.0, 0.3);
    let old_gas = *initial_field.uniform_state().unwrap().state();
    let new_gas = *replacement.uniform_state().unwrap().state();
    assert_eq!(initial_field.model().source_packs().len(), 2);
    assert_eq!(
        initial_field.uniform_state().unwrap().source_packs().len(),
        1
    );
    assert_eq!(replacement.uniform_state().unwrap().source_packs().len(), 2);
    assert_eq!(
        initial_field.model().source_packs(),
        replacement.model().source_packs()
    );
    let left = kelvin_wall(1000.0, 1e6, 2.0, 0.0002);
    let right = source_wall(1800.0, 4e6);
    let sources = vec![left, right.clone()];
    let side = DgSidePort {
        id: PortId(7),
        downstream_element: 1,
    };
    let mut bound = RegionalMaterialWallDg::new(
        sources.clone(),
        duct(&initial_field).with_side_ports(&[side]).unwrap(),
    )
    .unwrap();
    let mut propagation = bound.duct().propagator(2e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..30 {
        propagation.step(&mut state, &[]).unwrap();
    }
    propagation.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let motion = bound.duct().compliant_wall_state(&state).unwrap();
    assert!(motion.relaxation_stored_energy_j > 0.0);
    let mut control = state.clone();
    let control_propagation = propagation.clone();
    let account = bound
        .rebind_gas_field_controlled(&replacement, &mut state, &mut propagation, 1.0, &mut || {
            false
        })
        .unwrap();
    assert_eq!(bound.walls(), sources.as_slice());
    assert_eq!(bound.duct().gas_field(), Some(&replacement));
    assert_eq!(bound.duct().side_port_source(side.id).unwrap().gas, new_gas);
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    assert_eq!(
        state.cumulative_parameter_work_j(),
        account.parameter_work_j
    );
    assert_ne!(old_gas.dynamic_viscosity, new_gas.dynamic_viscosity);
    assert_ne!(old_gas.thermal_conductivity, new_gas.thermal_conductivity);
    let after = bound.duct().compliant_wall_state(&state).unwrap();
    assert_eq!(
        after.relaxation_stored_energy_j,
        motion.relaxation_stored_energy_j
    );
    let momentum_scale = old_gas.density
        * before
            .volume_flow()
            .iter()
            .map(|q| q.abs())
            .fold(1e-10, f64::max);
    for i in 0..bound.duct().field_dofs() {
        close(
            state.pressure()[i] / bulk(&new_gas),
            before.pressure()[i] / bulk(&old_gas),
        );
        assert!(
            (new_gas.density * state.volume_flow()[i] - old_gas.density * before.volume_flow()[i])
                .abs()
                <= 3e-10 * momentum_scale
        );
        close(after.displacement_m[i], motion.displacement_m[i]);
        close(after.velocity_m_s[i], motion.velocity_m_s[i]);
        for (a, b) in after
            .relaxation_pressure_pa
            .iter()
            .zip(&motion.relaxation_pressure_pa)
        {
            close(a[i], b[i]);
        }
    }
    let mut difference = 0.0;
    let mut signal = 0.0;
    for _ in 0..25 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        control_propagation.step(&mut control, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        control_propagation.sync_nodal(&mut control).unwrap();
        for (&p, &reference) in state.pressure().iter().zip(control.pressure()) {
            difference += (p - reference).powi(2);
            signal += reference * reference;
        }
        assert_eq!(
            frame.ledger.cumulative_parameter_work_j,
            account.parameter_work_j
        );
        close(
            state.stored_energy_j() + state.cumulative_dissipated_j(),
            before.initial_energy_j() + account.parameter_work_j,
        );
    }
    assert!(difference > 1e-6 * signal);
    assert_eq!(state.step_count(), 55);
    assert_eq!(propagation.dt(), control_propagation.dt());
    assert_eq!(propagation.tolerance(), control_propagation.tolerance());

    // The existing live wall-material owner must still work on the new gas
    // coordinates, retaining the other region's original Maxwell histories.
    let stiffer = kelvin_wall(1000.0, 2e6, 2.0, 0.0002);
    let wall_work = bound
        .rebind_controlled(
            vec![stiffer, right],
            &mut state,
            &mut propagation,
            1.0,
            &mut || false,
        )
        .unwrap();
    let total_work = account.parameter_work_j + wall_work.parameter_work_j;
    assert_eq!(bound.duct().gas_field(), Some(&replacement));
    assert_eq!(state.cumulative_parameter_work_j(), total_work);
    let frame = propagation.step(&mut state, &[]).unwrap();
    assert_eq!(frame.ledger.cumulative_parameter_work_j, total_work);
    close(
        state.stored_energy_j() + state.cumulative_dissipated_j(),
        before.initial_energy_j() + total_work,
    );
}

#[test]
fn nonuniform_gas_shape_budget_stale_propagator_and_late_cancel_never_publish_partial_gas() {
    let field = uniform_field(300.0, 1.0, 101_325.0, 0.3);
    let replacement = uniform_field(420.0, 0.2, 150_000.0, 0.3);
    let original = source_wall(1000.0, 1e6).bind_dg(duct(&field)).unwrap();
    let stepping = original.duct().propagator(2e-6, 1e-13).unwrap();
    let mut initial = original.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..20 {
        stepping.step(&mut initial, &[]).unwrap();
    }
    stepping.sync_nodal(&mut initial).unwrap();
    assert!(
        original
            .duct()
            .compliant_wall_state(&initial)
            .unwrap()
            .relaxation_stored_energy_j
            > 0.0
    );
    for (field, budget) in [
        (replacement.clone(), 0.0),
        (gas_field(0.2, 350.0).unwrap(), 1.0),
        (uniform_field(420.0, 0.2, 150_000.0, 0.4), 1.0),
    ] {
        let mut bound = original.clone();
        let mut state = initial.clone();
        let mut propagation = stepping.clone();
        assert!(
            bound
                .rebind_gas_field_controlled(
                    &field,
                    &mut state,
                    &mut propagation,
                    budget,
                    &mut || false
                )
                .is_err()
        );
        assert_eq!(state, initial);
        assert_eq!(bound.wall(), original.wall());
        assert_eq!(bound.duct().gas_field(), original.duct().gas_field());
        assert_eq!(
            bound.duct().operator_identity(),
            original.duct().operator_identity()
        );
        let mut control = initial.clone();
        assert_eq!(
            propagation.step(&mut state, &[]).unwrap(),
            stepping.step(&mut control, &[]).unwrap()
        );
    }
    let mut polls = 0;
    original
        .clone()
        .rebind_gas_field_controlled(
            &replacement,
            &mut initial.clone(),
            &mut stepping.clone(),
            1.0,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    let mut bound = original.clone();
    let mut state = initial.clone();
    let mut propagation = stepping.clone();
    assert!(matches!(
        bound.rebind_gas_field_controlled(
            &replacement,
            &mut state,
            &mut propagation,
            1.0,
            &mut || {
                count += 1;
                count == polls
            }
        ),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(state, initial);
    assert_eq!(bound.duct().gas_field(), original.duct().gas_field());
    assert_eq!(
        bound.duct().operator_identity(),
        original.duct().operator_identity()
    );
    let mut foreign = duct(&field).propagator(2e-6, 1e-13).unwrap();
    assert!(
        bound
            .rebind_gas_field_controlled(&replacement, &mut state, &mut foreign, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(state, initial);
    let mut control = initial.clone();
    assert_eq!(
        propagation.step(&mut state, &[]).unwrap(),
        stepping.step(&mut control, &[]).unwrap()
    );
}
