//! Resolved main-passage shear heats the same inventory during species transport.
mod aperture;

use super::*;
use fs_couple::bernoulli_aperture::wall::material::FiniteGasDgHeatBudget;
use fs_phs::linear_exact::DenseExponentialBudget;

fn absent_bore_sources() -> FiniteGasDgHeatBudget {
    let none = DenseExponentialBudget {
        max_dimension: 0,
        max_dense_bytes: 0,
        max_multiply_terms: 0,
    };
    FiniteGasDgHeatBudget {
        molecular: none,
        viscous: none,
    }
}

fn heated_network() -> (Parcels, MixedNodes, MixedInputs, Vec<SnapshotNozzleLink>) {
    let (parcels, nodes, mut inputs, links) = mixed();
    let FiniteGasAcousticInput::Passive(thermal) = inputs.remove(&10).unwrap() else {
        unreachable!("original complete-load fixture")
    };
    inputs.insert(
        10,
        FiniteGasAcousticInput::ViscousMolecularPort {
            thermal,
            heat_budgets: absent_bore_sources(),
        },
    );
    (parcels, nodes, inputs, links)
}

// Independent integral of mu*(du/dy)^2 on the original cell-centred no-slip
// field, including the half-cell wall distances. This is neither pressure
// work nor the backward-Euler numerical energy decrement.
fn passage_heat(load: &DgPortLoad) -> f64 {
    let channel = &load.terminal_channels()[0];
    let DgTerminalChannelKind::Slit {
        width_m,
        height_m,
        length_m,
        transverse_cells,
        ..
    } = channel.input().kind
    else {
        unreachable!("original slit fixture")
    };
    let velocity = channel.channel().unwrap().velocities_m_s();
    assert_eq!(velocity.len(), transverse_cells);
    let squares = 2.0 * velocity[0].powi(2)
        + velocity
            .windows(2)
            .map(|pair| (pair[1] - pair[0]).powi(2))
            .sum::<f64>()
        + 2.0 * velocity[transverse_cells - 1].powi(2);
    DT * channel.gas().dynamic_viscosity * width_m * length_m / (height_m / transverse_cells as f64)
        * squares
}

#[test]
fn g1_original_channel_shear_heats_the_transported_inventory_and_closes_energy() {
    let (mut parcels, mut nodes, inputs, links) = heated_network();
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

    // The independent reference retains the original medium while advancing
    // the SAME pressure/flow problem. Later composition may not reprice its heat.
    let FiniteGasAcousticNode::Port(original) = &old_nodes[&10] else {
        unreachable!()
    };
    let mut native = original.clone();
    let mut native_gas = old_gas[&10].clone();
    let input = thermal(&inputs[&10]);
    native
        .step_with_finite_gas_heat_controlled(
            &mut native_gas,
            input.ports[0].effort,
            0.0,
            input.gas,
            input.acoustic.unwrap(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(native_gas, old_gas[&10]);
    let expected_heat = passage_heat(&native);
    assert!(expected_heat > 0.0);
    let native_channel = native.last_terminal_channel_frames()[0].step.unwrap();
    assert!(native_channel.numerical_dissipation_j > 0.0);

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
    let FiniteGasAcousticStep::Port(step) = &report.acoustics[&10] else {
        unreachable!()
    };
    let FiniteGasAcousticNode::Port(next) = &nodes[&10] else {
        unreachable!()
    };
    assert_eq!(step.molecular_heat_j, 0.0);
    assert_eq!(step.viscous_heat_j, 0.0);
    assert!((step.terminal_viscous_heat_j / expected_heat - 1.0).abs() < 2e-13);
    assert_eq!(report.molecular_heat_j, 0.0);
    assert_eq!(report.viscous_heat_j, step.terminal_viscous_heat_j);
    assert_eq!(
        next.last_terminal_channel_frames(),
        native.last_terminal_channel_frames()
    );
    assert_eq!(next.last_terminal_frames(), native.last_terminal_frames());
    assert_eq!(
        next.terminal().unwrap().pole_states(),
        native.terminal().unwrap().pole_states()
    );
    assert_eq!(next.state().step_count(), original.state().step_count() + 1);
    assert_ne!(
        next.terminal_channels()[0].gas().dynamic_viscosity,
        original.terminal_channels()[0].gas().dynamic_viscosity
    );
    assert_eq!(
        next.duct().gas_field().unwrap().uniform_state(),
        Some(parcels[&10].state())
    );

    let external_heat: f64 = inputs
        .values()
        .map(|row| thermal(row).supplied_heat_j)
        .sum();
    assert!((report.supplied_heat_j - external_heat).abs() < 1e-18);
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, parcel)| gas_energy(parcel) - gas_energy(&old_gas[id]))
        .sum();
    let storage_change: f64 = nodes
        .iter()
        .map(|(id, node)| node.stored_energy_j() - old_nodes[id].stored_energy_j())
        .sum();
    let residual = gas_change + storage_change + report.acoustic_dissipation_j
        - expected_heat
        - external_heat
        - report.port_work_j;
    assert!(residual.abs() <= solve().energy_tolerance_j, "{residual:e}");
    assert!(
        (step.exported_acoustic_dissipation_j + expected_heat - step.acoustic_dissipation_j).abs()
            < 1e-20
    );
    for species in 0..2 {
        let before: f64 = old_gas
            .values()
            .map(|gas| species_masses(gas)[species])
            .sum();
        let after: f64 = parcels
            .values()
            .map(|gas| species_masses(gas)[species])
            .sum();
        assert!((after - before).abs() <= solve().transport.exchange.mass_tolerance_kg);
    }
}

#[test]
fn g4_complete_load_heat_refuses_wrong_owners_and_late_cancel_without_partial_transport() {
    let (parcels, nodes, inputs, links) = heated_network();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    for id in [30, 90] {
        let mut wrong = inputs.clone();
        let thermal = thermal(&wrong[&id]).clone();
        wrong.insert(
            id,
            FiniteGasAcousticInput::ViscousMolecularPort {
                thermal,
                heat_budgets: absent_bore_sources(),
            },
        );
        let (mut attempted_gas, mut attempted_nodes) = (parcels.clone(), nodes.clone());
        assert!(
            advance_finite_gas_acoustic_transport_with_paths(
                &mut attempted_gas,
                &mut attempted_nodes,
                &wrong,
                paths,
                DT,
                solve(),
                &mut || false,
            )
            .is_err()
        );
        assert_eq!(attempted_gas, parcels);
        equal_nodes(&attempted_nodes, &nodes);
    }

    let (mut expected_gas, mut expected_nodes) = (parcels.clone(), nodes.clone());
    let mut polls = 0;
    let expected = advance_finite_gas_acoustic_transport_with_paths(
        &mut expected_gas,
        &mut expected_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 1);
    let (mut attempted_gas, mut attempted_nodes) = (parcels.clone(), nodes.clone());
    let mut count = 0;
    let result = advance_finite_gas_acoustic_transport_with_paths(
        &mut attempted_gas,
        &mut attempted_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || {
            count += 1;
            count == polls
        },
    );
    assert!(matches!(result, Err(FiniteGasDgTransportError::Cancelled)));
    assert_eq!(count, polls);
    assert_eq!(attempted_gas, parcels);
    equal_nodes(&attempted_nodes, &nodes);
    let retry = advance_finite_gas_acoustic_transport_with_paths(
        &mut attempted_gas,
        &mut attempted_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(attempted_gas, expected_gas);
    equal_nodes(&attempted_nodes, &expected_nodes);
    assert_eq!(retry.transport, expected.transport);
    assert_eq!(retry.viscous_heat_j, expected.viscous_heat_j);
    assert_eq!(
        retry.exported_acoustic_dissipation_j,
        expected.exported_acoustic_dissipation_j
    );
}
