//! Both physical source heats compose with one original material-transport stage.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{
    FiniteGasDgBinding, FiniteGasDgNode, FiniteGasDgTransportError, FiniteGasDgTransportInput,
    FiniteGasDgTransportSolve,
};
use fs_couple::dg_port::finite_gas::transport::assembly::{
    FiniteGasAcousticInput, FiniteGasAcousticNode, FiniteGasAcousticStep,
    advance_finite_gas_acoustic_transport_with_paths,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::network::{
        SpeciesNetworkSolve,
        transport::{GasTransportPaths, GasTransportSolve},
    },
};
use std::collections::BTreeMap;

type Parcels = BTreeMap<u64, ClosedGasParcel>;
type Nodes = BTreeMap<u64, FiniteGasAcousticNode>;
type Inputs = BTreeMap<u64, FiniteGasAcousticInput>;

fn transport_solve() -> FiniteGasDgTransportSolve {
    FiniteGasDgTransportSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.05,
                mass_tolerance_kg: 1e-16,
                energy_tolerance_j: 1e-11,
            },
            entropy_tolerance_j_k: 1e-12,
        },
        energy_tolerance_j: 5e-11,
    }
}

fn fixture() -> (Parcels, Nodes, Inputs, SnapshotNozzleLink) {
    let (binding, state, propagator, gas) = regional_both();
    let parcels = BTreeMap::from([(0, gas), (1, parcel(1.0, 300.0, 140000.0))]);
    let nodes = BTreeMap::from([(
        0,
        FiniteGasAcousticNode::Wall(FiniteGasDgNode {
            binding: FiniteGasDgBinding::Regional(binding),
            state,
            propagator,
        }),
    )]);
    let inputs = BTreeMap::from([
        (
            0,
            FiniteGasAcousticInput::ViscousMolecularWall {
                thermal: FiniteGasDgTransportInput {
                    gas: energy_solve(),
                    supplied_heat_j: 0.0,
                    acoustic: Some(solve()),
                    ports: vec![PortInput {
                        id: PortId(7),
                        effort: 1e-7,
                    }],
                },
                heat_budgets: budgets(),
            },
        ),
        (
            1,
            FiniteGasAcousticInput::Passive(FiniteGasDgTransportInput {
                gas: energy_solve(),
                supplied_heat_j: 0.0,
                acoustic: None,
                ports: vec![],
            }),
        ),
    ]);
    let nozzle = SnapshotNozzleLink {
        id: 9,
        parcels: [1, 0],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: energy_solve(),
    };
    (parcels, nodes, inputs, nozzle)
}

fn wall(nodes: &Nodes) -> &FiniteGasDgNode {
    let FiniteGasAcousticNode::Wall(node) = &nodes[&0] else {
        panic!("source wall required")
    };
    node
}

fn unchanged(parcels: &Parcels, nodes: &Nodes, original: &(Parcels, Nodes)) {
    assert_eq!(parcels, &original.0);
    let (node, old) = (wall(nodes), wall(&original.1));
    assert_eq!(node.state, old.state);
    assert_eq!(
        node.binding.duct().operator_identity(),
        old.binding.duct().operator_identity()
    );
    assert_eq!(
        node.binding.duct().gas_field(),
        old.binding.duct().gas_field()
    );
    let (mut left, mut right) = (node.state.clone(), old.state.clone());
    assert_eq!(
        node.propagator.step(&mut left, &[]).unwrap(),
        old.propagator.step(&mut right, &[]).unwrap()
    );
    assert_eq!(left, right);
}

#[test]
fn joint_gas_heat_and_nozzle_transport_share_original_sources_and_cancel_once() {
    let (mut parcels, mut nodes, inputs, nozzle) = fixture();
    let (original_parcels, original_nodes) = (parcels.clone(), nodes.clone());
    let node = wall(&nodes);
    let ports = [PortInput {
        id: PortId(7),
        effort: 1e-7,
    }];
    let molecular = node
        .binding
        .duct()
        .prepare_molecular_heat(&node.propagator, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&node.state, &ports, &mut || false)
        .unwrap();
    let viscous = node
        .binding
        .duct()
        .prepare_viscous_heat(&node.propagator, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&node.state, &ports, &mut || false)
        .unwrap();
    assert!(molecular > 0.0 && viscous > 0.0);
    let mut marched = node.state.clone();
    let native = node.propagator.step(&mut marched, &ports).unwrap();
    let motion = node.binding.duct().compliant_wall_state(&marched).unwrap();
    let flow = parcels[&1]
        .nozzle_flow_to(
            parcels[&0].state().state().pressure,
            nozzle.nozzle,
            energy_solve(),
            || false,
        )
        .unwrap();
    let report = advance_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        GasTransportPaths {
            nozzles: std::slice::from_ref(&nozzle),
            ..Default::default()
        },
        1e-6,
        transport_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport.nozzles[0].flow, flow);
    assert_eq!(
        report.transport.nozzles[0].transferred_mass_kg,
        flow.mass_flow_kg_s * 1e-6
    );
    assert!(parcels[&0].mass_kg() > original_parcels[&0].mass_kg());
    assert_eq!(report.molecular_heat_j, molecular);
    assert_eq!(report.viscous_heat_j, viscous);
    assert_eq!(report.supplied_heat_j, 0.0);
    let FiniteGasAcousticStep::Wall(local) = &report.acoustics[&0] else {
        panic!("wall report required")
    };
    assert_eq!(local.molecular_heat_j, molecular);
    assert_eq!(local.viscous_heat_j, viscous);
    assert_eq!(
        report.exported_acoustic_dissipation_j,
        report.acoustic_dissipation_j - (molecular + viscous)
    );
    assert_eq!(
        local.acoustic.ledger.port_work_step_j,
        native.ledger.port_work_step_j
    );
    assert_eq!(
        local.acoustic.port_outputs[0].delta_x,
        native.port_outputs[0].delta_x
    );
    let current = wall(&nodes);
    assert_eq!(
        current.state.step_count(),
        wall(&original_nodes).state.step_count() + 1
    );
    same_motion(
        &current
            .binding
            .duct()
            .compliant_wall_state(&current.state)
            .unwrap(),
        &motion,
    );
    assert_eq!(
        current.binding.duct().gas_field().unwrap().uniform_state(),
        Some(parcels[&0].state())
    );
    assert_eq!(
        current.binding.duct().molecular_relaxation_law(),
        wall(&original_nodes)
            .binding
            .duct()
            .molecular_relaxation_law()
    );
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy_change(&original_parcels[id], gas))
        .sum();
    let wave_change =
        current.state.stored_energy_j() - wall(&original_nodes).state.stored_energy_j();
    let residual =
        gas_change + wave_change + report.exported_acoustic_dissipation_j - report.port_work_j;
    assert!(residual.abs() <= transport_solve().energy_tolerance_j);
    assert!(report.combined_energy_residual_j.abs() <= transport_solve().energy_tolerance_j);
    assert!(
        report.total_mass_residual_kg.abs()
            <= transport_solve().transport.exchange.mass_tolerance_kg
    );
}

#[test]
fn joint_transport_refuses_missing_source_allowance_or_late_cancellation_atomically() {
    let (mut parcels, mut nodes, inputs, nozzle) = fixture();
    let original = (parcels.clone(), nodes.clone());
    let mut missing = inputs.clone();
    let FiniteGasAcousticInput::ViscousMolecularWall { heat_budgets, .. } =
        missing.get_mut(&0).unwrap()
    else {
        panic!("joint source input required")
    };
    heat_budgets.viscous = no_budget();
    assert!(matches!(
        advance_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &missing,
            GasTransportPaths {
                nozzles: std::slice::from_ref(&nozzle),
                ..Default::default()
            },
            1e-6,
            transport_solve(),
            &mut || false,
        ),
        Err(FiniteGasDgTransportError::Acoustic { parcel: 0, .. })
    ));
    unchanged(&parcels, &nodes, &original);
    let mut wrong_owner = inputs.clone();
    let FiniteGasAcousticInput::Passive(thermal) = inputs[&1].clone() else {
        unreachable!()
    };
    wrong_owner.insert(
        1,
        FiniteGasAcousticInput::ViscousMolecularWall {
            thermal,
            heat_budgets: budgets(),
        },
    );
    assert!(matches!(
        advance_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &wrong_owner,
            GasTransportPaths::default(),
            1e-6,
            transport_solve(),
            &mut || false,
        ),
        Err(FiniteGasDgTransportError::Invalid(_))
    ));
    unchanged(&parcels, &nodes, &original);

    let (mut expected_parcels, mut expected_nodes) = original.clone();
    let mut successful_polls = 0;
    let expected = advance_finite_gas_acoustic_transport_with_paths(
        &mut expected_parcels,
        &mut expected_nodes,
        &inputs,
        GasTransportPaths {
            nozzles: std::slice::from_ref(&nozzle),
            ..Default::default()
        },
        1e-6,
        transport_solve(),
        &mut || {
            successful_polls += 1;
            false
        },
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        advance_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &inputs,
            GasTransportPaths {
                nozzles: std::slice::from_ref(&nozzle),
                ..Default::default()
            },
            1e-6,
            transport_solve(),
            &mut || {
                polls += 1;
                polls == successful_polls
            },
        ),
        Err(FiniteGasDgTransportError::Cancelled)
    ));
    assert_eq!(polls, successful_polls);
    unchanged(&parcels, &nodes, &original);
    let retry = advance_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        GasTransportPaths {
            nozzles: std::slice::from_ref(&nozzle),
            ..Default::default()
        },
        1e-6,
        transport_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, expected_parcels);
    assert_eq!(wall(&nodes).state, wall(&expected_nodes).state);
    assert_eq!(retry.molecular_heat_j, expected.molecular_heat_j);
    assert_eq!(retry.viscous_heat_j, expected.viscous_heat_j);
}
