//! Different original acoustic owners share one conservative material snapshot.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, FiniteGasDgBackground, FiniteGasDgBinding, FiniteGasDgNode, LocalMaterialWall,
    WallMaterialLoss,
};
use fs_couple::dg_port::finite_gas::transport::assembly::{
    FiniteGasAcousticInput, FiniteGasAcousticNode, FiniteGasAcousticStep,
    advance_finite_gas_acoustic_transport_with_paths,
    finite_gas_acoustic_transport_candidate_with_paths,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_qty::{Density, Dims, Pressure};

type MixedNodes = BTreeMap<u64, FiniteGasAcousticNode>;
type MixedInputs = BTreeMap<u64, FiniteGasAcousticInput>;

fn wall_node() -> (FiniteGasDgNode, ClosedGasParcel, FiniteGasDgTransportInput) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 1e6),
        ("poisson_ratio", Dims::NONE, 0.3),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "authored elastic wall for transport conformance, not measured material"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-wall".into(),
            phase: "solid".into(),
            process: "manufactured integration fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let field = field(300.0, 1.0, 101_325.0, LENGTH);
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &field,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap();
    let wall = LocalMaterialWall::new(
        material,
        2e-5,
        AxialConstraint::StressFree,
        WallMaterialLoss::Elastic,
        [20.0, 100_000.0],
    )
    .unwrap()
    .bind_dg(duct)
    .unwrap();
    let propagator = wall.duct().propagator(DT, 1e-11).unwrap();
    let mut state = wall
        .duct()
        .project_initial_state(|x| 10.0 + 200.0 * x, |_| 1e-9);
    let ports = vec![fs_duct::dg::PortInput {
        id: PortId(7),
        effort: 1e-9,
    }];
    for _ in 0..24 {
        propagator.step(&mut state, &ports).unwrap();
    }
    propagator.sync_nodal(&mut state).unwrap();
    let gas = ClosedGasParcel::new(
        field.model().with_composition(composition(1.0)).unwrap(),
        point(300.0, 101_325.0),
        core::f64::consts::PI * RADIUS * RADIUS * LENGTH,
    )
    .unwrap();
    (
        FiniteGasDgNode {
            binding: FiniteGasDgBinding::Uniform(wall),
            state,
            propagator,
        },
        gas,
        FiniteGasDgTransportInput {
            gas: energy_solve(),
            supplied_heat_j: 1e-5,
            acoustic: Some(controls()),
            ports,
        },
    )
}

fn mixed() -> (Parcels, MixedNodes, MixedInputs, Vec<SnapshotNozzleLink>) {
    let (mut parcels, mut sources, mut input) = network();
    let passive = sources.remove(&10).unwrap().load().clone();
    let mut passive_input = input.remove(&10).unwrap().thermal;
    passive_input.ports = vec![fs_duct::dg::PortInput {
        id: passive.inlet(),
        effort: 1e-8,
    }];
    let paired = sources.remove(&20).unwrap();
    let (wall, wall_gas, wall_input) = wall_node();
    parcels.insert(30, wall_gas);
    let nodes = BTreeMap::from([
        (10, FiniteGasAcousticNode::Port(passive)),
        (20, FiniteGasAcousticNode::Aperture(paired)),
        (30, FiniteGasAcousticNode::Wall(wall)),
    ]);
    let inputs = BTreeMap::from([
        (10, FiniteGasAcousticInput::Passive(passive_input)),
        (
            20,
            FiniteGasAcousticInput::Aperture(input.remove(&20).unwrap()),
        ),
        (30, FiniteGasAcousticInput::Passive(wall_input)),
        (
            90,
            FiniteGasAcousticInput::Passive(input.remove(&90).unwrap().thermal),
        ),
    ]);
    let mut links = nozzles().to_vec();
    links.push(SnapshotNozzleLink {
        id: 42,
        parcels: [90, 30],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: energy_solve(),
    });
    (parcels, nodes, inputs, links)
}

fn thermal(input: &FiniteGasAcousticInput) -> &FiniteGasDgTransportInput {
    match input {
        FiniteGasAcousticInput::Passive(input) => input,
        FiniteGasAcousticInput::MolecularWall { thermal, .. } => thermal,
        FiniteGasAcousticInput::Aperture(input) => &input.thermal,
    }
}

fn equal_nodes(actual: &MixedNodes, expected: &MixedNodes) {
    assert!(actual.keys().eq(expected.keys()));
    for (&id, node) in actual {
        assert_eq!(node.stored_energy_j(), expected[&id].stored_energy_j());
        assert_eq!(node.accepted_steps(), expected[&id].accepted_steps());
        assert_eq!(node.time_step_s(), expected[&id].time_step_s());
        match (node, &expected[&id]) {
            (FiniteGasAcousticNode::Wall(a), FiniteGasAcousticNode::Wall(b)) => {
                assert_eq!(a.state, b.state);
                assert_eq!(
                    a.binding.duct().operator_identity(),
                    b.binding.duct().operator_identity()
                );
                assert_eq!(a.binding.duct().gas_field(), b.binding.duct().gas_field());
            }
            (FiniteGasAcousticNode::Port(a), FiniteGasAcousticNode::Port(b)) => {
                assert_eq!(a.state(), b.state());
                assert_eq!(a.duct().gas_field(), b.duct().gas_field());
                assert_eq!(
                    a.terminal().unwrap().pole_states(),
                    b.terminal().unwrap().pole_states()
                );
                assert_eq!(a.last_terminal_frames(), b.last_terminal_frames());
                assert_eq!(
                    a.last_terminal_channel_frames(),
                    b.last_terminal_channel_frames()
                );
                assert_eq!(
                    a.terminal_channels()[0].channel().unwrap().velocities_m_s(),
                    b.terminal_channels()[0].channel().unwrap().velocities_m_s()
                );
            }
            (FiniteGasAcousticNode::Aperture(a), FiniteGasAcousticNode::Aperture(b)) => {
                unchanged(
                    &BTreeMap::from([(id, a.clone())]),
                    &BTreeMap::from([(id, b.clone())]),
                );
            }
            _ => panic!("original acoustic owner changed"),
        }
    }
}

#[test]
fn g1_all_three_original_owners_share_species_transport_and_one_energy_account() {
    let (mut parcels, mut nodes, inputs, links) = mixed();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let endpoints = inputs
        .iter()
        .map(|(&id, input)| (id, thermal(input).gas))
        .collect();
    let transport = material_transport::step_candidate_with_paths(
        &parcels,
        paths,
        DT,
        &endpoints,
        solve().transport,
        || false,
    )
    .unwrap();
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
    assert_eq!(report.transport, *transport.step());
    assert_eq!(report.acoustics.len(), 3);
    assert_eq!(report.reservoirs.len(), 1);
    let mut native = old_nodes.clone();
    let (mut net_work, mut losses) = (0.0, 0.0);
    for (&id, node) in &mut native {
        let input = thermal(&inputs[&id]);
        let transported = &transport.parcels()[&id];
        let enthalpy = transport
            .step()
            .exchange
            .parcels
            .iter()
            .find(|row| row.id == id)
            .unwrap()
            .net_enthalpy_j;
        let mut gas = old_gas[&id].clone();
        let evolve = |_: &ClosedGasParcel, _: f64, work: f64, cancel: &mut fn() -> bool| {
            let mut next = transported.clone();
            let source = next.replace_model_with_heat_work_volume(
                transported.model().clone(),
                transported.volume_m3(),
                input.supplied_heat_j,
                work,
                input.gas,
                cancel,
            )?;
            Ok::<_, GasEnergyError>(FiniteGasDgBackground {
                gas: next,
                supplied_heat_j: input.supplied_heat_j,
                transported_enthalpy_j: enthalpy,
                source,
            })
        };
        let mut cancel: fn() -> bool = || false;
        match (node, &report.acoustics[&id]) {
            (FiniteGasAcousticNode::Wall(node), FiniteGasAcousticStep::Wall(actual)) => {
                let FiniteGasDgBinding::Uniform(wall) = &mut node.binding else {
                    unreachable!()
                };
                let expected = wall
                    .step_with_finite_gas_controlled(
                        &mut node.state,
                        &mut node.propagator,
                        &mut gas,
                        &input.ports,
                        input.acoustic.unwrap(),
                        evolve,
                        &mut cancel,
                    )
                    .unwrap();
                assert_eq!(actual.acoustic, expected.acoustic);
                assert_eq!(actual.refresh, expected.refresh);
                net_work += input
                    .ports
                    .iter()
                    .map(|port| {
                        port.effort
                            * actual
                                .acoustic
                                .port_outputs
                                .iter()
                                .find(|p| p.id == port.id)
                                .unwrap()
                                .delta_x
                    })
                    .sum::<f64>();
                losses += actual.acoustic_dissipation_j;
                assert_eq!(
                    node.binding.duct().gas_field().unwrap().uniform_state(),
                    Some(gas.state())
                );
            }
            (FiniteGasAcousticNode::Port(node), FiniteGasAcousticStep::Port(actual)) => {
                let expected = node
                    .step_with_finite_gas_controlled(
                        &mut gas,
                        input.ports[0].effort,
                        input.acoustic.unwrap(),
                        evolve,
                        &mut cancel,
                    )
                    .unwrap();
                assert_eq!(actual.acoustic, expected.acoustic);
                assert_eq!(actual.refresh, expected.refresh);
                let inlet = actual
                    .acoustic
                    .port_outputs
                    .iter()
                    .find(|p| p.id == node.inlet())
                    .unwrap();
                net_work += input.ports[0].effort * inlet.delta_x
                    - node
                        .last_terminal_frames()
                        .iter()
                        .map(|t| t.load_work_j)
                        .sum::<f64>();
                losses += actual.acoustic_dissipation_j;
                assert_eq!(
                    node.duct().gas_field().unwrap().uniform_state(),
                    Some(gas.state())
                );
            }
            (FiniteGasAcousticNode::Aperture(node), FiniteGasAcousticStep::Aperture(actual)) => {
                let FiniteGasAcousticInput::Aperture(source) = &inputs[&id] else {
                    unreachable!()
                };
                let source = source.source.as_ref().unwrap();
                let expected = node
                    .step_with_finite_gas_controlled(
                        &mut gas,
                        source.drive,
                        source.forces_n,
                        source.pressures_pa,
                        input.acoustic.unwrap(),
                        evolve,
                        &mut cancel,
                    )
                    .unwrap();
                assert_eq!(actual.coupled, expected.coupled);
                assert_eq!(actual.gas.refresh, expected.gas.refresh);
                let frame = &actual.coupled;
                net_work +=
                    frame.upstream_work_j + frame.body_work_j + frame.aperture.mechanical_work_j
                        - frame.terminals.iter().map(|t| t.load_work_j).sum::<f64>();
                losses += frame.inferred_total_dissipation_j;
                assert_eq!(
                    node.load().duct().gas_field().unwrap().uniform_state(),
                    Some(gas.state())
                );
            }
            _ => panic!("mixed report was relabeled to another physical owner"),
        }
        assert_eq!(gas, parcels[&id]);
        assert!(species_masses(&gas)[1] > 0.0);
    }
    equal_nodes(&nodes, &native);
    for species in 0..2 {
        let before: f64 = old_gas.values().map(|p| species_masses(p)[species]).sum();
        let after: f64 = parcels.values().map(|p| species_masses(p)[species]).sum();
        assert!((after - before).abs() <= solve().transport.exchange.mass_tolerance_kg);
    }
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, p)| gas_energy(p) - gas_energy(&old_gas[id]))
        .sum();
    let storage_change: f64 = nodes
        .iter()
        .map(|(id, node)| node.stored_energy_j() - old_nodes[id].stored_energy_j())
        .sum();
    let heat: f64 = inputs
        .values()
        .map(|input| thermal(input).supplied_heat_j)
        .sum();
    let residual = gas_change + storage_change + losses - heat - net_work;
    assert!(residual.abs() <= solve().energy_tolerance_j, "{residual:e}");
    close(report.port_work_j, net_work, net_work.abs());
    close(
        report.acoustic_storage_change_j,
        storage_change,
        storage_change.abs(),
    );
    assert!(report.combined_energy_residual_j.abs() <= solve().energy_tolerance_j);
    assert!(nodes.values().all(|node| node.accepted_steps() == 25));
}

#[test]
fn g4_mistyped_or_cancelled_assembly_preserves_every_original_owner_and_retries() {
    let (mut parcels, mut nodes, inputs, links) = mixed();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    for case in 0..3 {
        let mut bad = inputs.clone();
        match case {
            0 => {
                bad.insert(10, inputs[&20].clone());
            }
            1 => {
                bad.insert(
                    20,
                    FiniteGasAcousticInput::Passive(thermal(&inputs[&20]).clone()),
                );
            }
            _ => {
                bad.insert(
                    90,
                    FiniteGasAcousticInput::Aperture(FiniteGasApertureTransportInput {
                        thermal: thermal(&inputs[&90]).clone(),
                        source: None,
                    }),
                );
            }
        }
        assert!(
            advance_finite_gas_acoustic_transport_with_paths(
                &mut parcels,
                &mut nodes,
                &bad,
                paths,
                DT,
                solve(),
                &mut || false
            )
            .is_err()
        );
        assert_eq!(parcels, old_gas);
        equal_nodes(&nodes, &old_nodes);
    }
    let preview = finite_gas_acoustic_transport_candidate_with_paths(
        &parcels,
        &nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, old_gas);
    equal_nodes(&nodes, &old_nodes);
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
    assert_eq!(preview.parcels(), &expected_gas);
    equal_nodes(preview.nodes(), &expected_nodes);
    let mut seen = 0;
    assert!(matches!(
        advance_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &inputs,
            paths,
            DT,
            solve(),
            &mut || {
                seen += 1;
                seen == polls
            }
        ),
        Err(FiniteGasDgTransportError::Cancelled)
    ));
    assert_eq!(seen, polls);
    assert_eq!(parcels, old_gas);
    equal_nodes(&nodes, &old_nodes);
    let retry = advance_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, expected_gas);
    equal_nodes(&nodes, &expected_nodes);
    assert_eq!(retry.transport, expected.transport);
    assert_eq!(
        retry.combined_energy_residual_j,
        expected.combined_energy_residual_j
    );
    // Continue another full mixed physical step, retaining the native clocks.
    advance_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    advance_finite_gas_acoustic_transport_with_paths(
        &mut expected_gas,
        &mut expected_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, expected_gas);
    equal_nodes(&nodes, &expected_nodes);
}
