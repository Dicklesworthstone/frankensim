//! The original moving source's solved flow heats its transported gas once.
use super::*;

fn moving_heat_network() -> (Parcels, MixedNodes, MixedInputs, Vec<SnapshotNozzleLink>) {
    let (parcels, nodes, mut inputs, links) = mixed();
    let FiniteGasAcousticInput::Aperture(aperture) = inputs.remove(&20).unwrap() else {
        unreachable!("original moving-source fixture")
    };
    inputs.insert(
        20,
        FiniteGasAcousticInput::ViscousMolecularAperture {
            aperture,
            heat_budgets: absent_bore_sources(),
        },
    );
    (parcels, nodes, inputs, links)
}

#[test]
fn g1_moving_source_shear_heats_transported_gas_without_repeating_mechanics() {
    let (mut parcels, mut nodes, inputs, links) = moving_heat_network();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let endpoints = inputs
        .iter()
        .map(|(&id, row)| (id, thermal(row).gas))
        .collect();
    let transported = material_transport::step_candidate_with_paths(
        &parcels,
        paths,
        DT,
        &endpoints,
        solve().transport,
        || false,
    )
    .unwrap();
    let FiniteGasAcousticNode::Aperture(original) = &old_nodes[&20] else {
        unreachable!()
    };
    let FiniteGasAcousticInput::ViscousMolecularAperture { aperture, .. } = &inputs[&20]
    else {
        unreachable!()
    };
    let source = aperture.source.unwrap();
    let mut native = original.clone();
    let frozen_frame = native
        .step_with_forces_and_pressures(source.drive, source.forces_n, source.pressures_pa)
        .unwrap();
    // This integral uses the actual velocity profile from the ORIGINAL medium
    // and the flow solved by the moving source, before any composition refresh.
    let expected_heat = passage_heat(native.load());
    assert!(expected_heat > 0.0);
    assert!(native.load().last_terminal_channel_frames()[0].step.unwrap().numerical_dissipation_j > 0.0);
    assert!(frozen_frame.aperture.mechanical_work_j.abs() > 0.0);

    let report = advance_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport, *transported.step());
    let FiniteGasAcousticStep::Aperture(step) = &report.acoustics[&20] else {
        unreachable!()
    };
    let FiniteGasAcousticNode::Aperture(next) = &nodes[&20] else {
        unreachable!()
    };
    assert_eq!(step.coupled.aperture, frozen_frame.aperture);
    assert_eq!(step.gas.molecular_heat_j, 0.0);
    assert_eq!(step.gas.viscous_heat_j, 0.0);
    assert!((step.gas.terminal_viscous_heat_j / expected_heat - 1.0).abs() < 2e-13);
    assert_eq!(report.molecular_heat_j, 0.0);
    assert_eq!(report.viscous_heat_j, step.gas.terminal_viscous_heat_j);
    assert_eq!(next.duct_state().step_count(), original.duct_state().step_count() + 1);
    assert_eq!(next.aperture().accepted_steps(), original.aperture().accepted_steps() + 1);
    assert_eq!(next.load().last_terminal_frames(), native.load().last_terminal_frames());
    assert_eq!(
        next.load().last_terminal_channel_frames(),
        native.load().last_terminal_channel_frames()
    );
    assert_eq!(
        next.load().terminal().unwrap().pole_states(),
        native.load().terminal().unwrap().pole_states()
    );
    let (LoadedAperture::Paired(actual), LoadedAperture::Paired(expected)) =
        (next.aperture(), native.aperture())
    else {
        unreachable!("original paired-source fixture")
    };
    assert_eq!(actual.blade_states(), expected.blade_states());
    assert_eq!(
        next.load().duct().gas_field().unwrap().uniform_state(),
        Some(parcels[&20].state())
    );
    assert_ne!(
        next.load().terminal_channels()[0].gas().dynamic_viscosity,
        original.load().terminal_channels()[0].gas().dynamic_viscosity
    );

    let external_heat: f64 = inputs.values().map(|row| thermal(row).supplied_heat_j).sum();
    assert!((report.supplied_heat_j - external_heat).abs() < 1e-18);
    let gas_change: f64 = parcels.iter()
        .map(|(id, parcel)| gas_energy(parcel) - gas_energy(&old_gas[id]))
        .sum();
    let storage_change: f64 = nodes.iter()
        .map(|(id, node)| node.stored_energy_j() - old_nodes[id].stored_energy_j())
        .sum();
    let residual = gas_change + storage_change + report.acoustic_dissipation_j
        - expected_heat - external_heat - report.port_work_j;
    assert!(residual.abs() <= solve().energy_tolerance_j, "{residual:e}");
    assert!(step.combined_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(
        (step.gas.exported_acoustic_dissipation_j + expected_heat
            - step.gas.acoustic_dissipation_j).abs() < 1e-20
    );
    for species in 0..2 {
        let before: f64 = old_gas.values().map(|gas| species_masses(gas)[species]).sum();
        let after: f64 = parcels.values().map(|gas| species_masses(gas)[species]).sum();
        assert!((after - before).abs() <= solve().transport.exchange.mass_tolerance_kg);
    }
}

#[test]
fn g4_moving_source_heat_rejects_wrong_owners_and_cancels_before_joint_publication() {
    let (parcels, nodes, inputs, links) = moving_heat_network();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    for id in [10, 30, 90] {
        let mut wrong = inputs.clone();
        wrong.insert(id, inputs[&20].clone());
        let (mut attempted_gas, mut attempted_nodes) = (parcels.clone(), nodes.clone());
        assert!(advance_finite_gas_acoustic_transport_with_paths(
            &mut attempted_gas, &mut attempted_nodes, &wrong, paths, DT, solve(),
            &mut || false,
        ).is_err());
        assert_eq!(attempted_gas, parcels);
        equal_nodes(&attempted_nodes, &nodes);
    }

    let (mut expected_gas, mut expected_nodes) = (parcels.clone(), nodes.clone());
    let mut polls = 0;
    let expected = advance_finite_gas_acoustic_transport_with_paths(
        &mut expected_gas, &mut expected_nodes, &inputs, paths, DT, solve(),
        &mut || { polls += 1; false },
    ).unwrap();
    assert!(polls > 1);
    let (mut attempted_gas, mut attempted_nodes) = (parcels.clone(), nodes.clone());
    let mut count = 0;
    let result = advance_finite_gas_acoustic_transport_with_paths(
        &mut attempted_gas, &mut attempted_nodes, &inputs, paths, DT, solve(),
        &mut || { count += 1; count == polls },
    );
    assert!(matches!(result, Err(FiniteGasDgTransportError::Cancelled)));
    assert_eq!(count, polls);
    assert_eq!(attempted_gas, parcels);
    equal_nodes(&attempted_nodes, &nodes);
    let retry = advance_finite_gas_acoustic_transport_with_paths(
        &mut attempted_gas, &mut attempted_nodes, &inputs, paths, DT, solve(),
        &mut || false,
    ).unwrap();
    assert_eq!(attempted_gas, expected_gas);
    equal_nodes(&attempted_nodes, &expected_nodes);
    assert_eq!(retry.transport, expected.transport);
    assert_eq!(retry.viscous_heat_j, expected.viscous_heat_j);
    assert_eq!(retry.exported_acoustic_dissipation_j, expected.exported_acoustic_dissipation_j);
}
