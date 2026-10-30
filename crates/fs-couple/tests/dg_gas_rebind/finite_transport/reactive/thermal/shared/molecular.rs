//! Source-owned molecular heat enters the original transported gas and one body.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, FiniteGasDgBinding, FiniteGasDgLoss, FiniteGasDgNode, LocalMaterialWall,
    WallMaterialLoss,
};
use fs_duct::dg::PortInput;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_qty::Dims;

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        // The source owner admits a conservative 2048*(2n)^3 preparation bound.
        max_multiply_terms: 40_000_000_000,
    }
}

fn molecular_wall(law: &GasReaction, active: bool) -> (FiniteGasDgNode, ClosedGasParcel) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 4000.0),
        ("young_modulus", Pressure::DIMS, 1e9),
        ("poisson_ratio", Dims::NONE, 0.3),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored elastic wall for molecular-transfer conformance".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                name,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic molecular-wall fixture".into(),
            phase: "solid".into(),
            process: "elastic".into(),
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
    // The native molecular history owner currently admits fixed species fractions.
    // Pure I permits real same-composition inflow and the mechanism's inert branch;
    // the other native volumes continue their nonzero chemical reactions.
    let model = law
        .catalog()
        .with_composition(mixture([0.0, 0.0, 1.0]))
        .unwrap();
    let field = AxialGasField::new(
        vec![0.0, LENGTH],
        101_325.0,
        vec![
            PrimitiveGasState {
                temperature_k: 300.0,
                composition: mixture([0.0, 0.0, 1.0]),
            };
            2
        ],
        model.clone(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let relaxation = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: SpeciesId::new("I").unwrap(),
        processes: if active {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 900.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 3.0,
                    t_ref_k: 300.0,
                    exponent: 1.2,
                    t_min_k: 200.0,
                    t_max_k: 900.0,
                },
            }]
        } else {
            Vec::new()
        },
    }])
    .unwrap();
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
    .unwrap()
    .with_viscothermal(8, [50.0, 2000.0], 0.02)
    .unwrap()
    .with_molecular_relaxation_field(relaxation)
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
        .project_initial_state(|x| 1000.0 + 3000.0 * (1.0 - 2.0 * x / LENGTH), |_| 1e-9);
    for _ in 0..24 {
        propagator.step(&mut state, &wall_ports()).unwrap();
    }
    propagator.sync_nodal(&mut state).unwrap();
    if active {
        assert!(wall.duct().molecular_dissipated_power_w(&state).unwrap() > 0.0);
    }
    let gas = ClosedGasParcel::new(
        model,
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
    )
}

fn wall_ports() -> Vec<PortInput> {
    vec![PortInput {
        id: PortId(7),
        effort: 1e-9,
    }]
}

fn molecular_fixture(law: &GasReaction, active: bool) -> (Parcels, MixedNodes, MixedInputs) {
    let (mut parcels, mut nodes, mut inputs) = fixture(law);
    let (wall, gas) = molecular_wall(law, active);
    parcels.insert(30, gas);
    nodes.insert(30, FiniteGasAcousticNode::Wall(wall));
    inputs.insert(
        30,
        FiniteGasAcousticInput::MolecularWall {
            thermal: FiniteGasDgTransportInput {
                gas: energy_solve(),
                supplied_heat_j: 0.0,
                acoustic: Some(
                    fs_couple::bernoulli_aperture::wall::material::FiniteGasDgSolve {
                        loss: FiniteGasDgLoss::UniformMolecularDiffusion {
                            maximum_relative_error: 0.02,
                        },
                        ..controls()
                    },
                ),
                ports: wall_ports(),
            },
            heat_budget: heat_budget(),
        },
    );
    (parcels, nodes, inputs)
}

fn reactions(law: &GasReaction) -> BTreeMap<u64, FiniteGasAcousticReaction<'_>> {
    let mut result = bindings(law);
    result.insert(
        30,
        FiniteGasAcousticReaction {
            reaction: law,
            controls: chemistry(),
        },
    );
    result
}

fn links() -> Vec<SnapshotNozzleLink> {
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
    links
}

fn export_only(inputs: &MixedInputs) -> MixedInputs {
    let mut result = inputs.clone();
    result.insert(
        30,
        FiniteGasAcousticInput::Passive(thermal(&inputs[&30]).clone()),
    );
    result
}

fn original_heat(nodes: &MixedNodes) -> f64 {
    let FiniteGasAcousticNode::Wall(wall) = &nodes[&30] else {
        unreachable!()
    };
    wall.binding
        .duct()
        .prepare_molecular_heat(&wall.propagator, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&wall.state, &wall_ports(), &mut || false)
        .unwrap()
}

// Independent time quadrature of the actual molecular source, not the prepared
// quadratic integral or an allocation of the total acoustic energy decrement.
fn quadrature_heat(nodes: &MixedNodes) -> f64 {
    let FiniteGasAcousticNode::Wall(wall) = &nodes[&30] else {
        unreachable!()
    };
    let duct = wall.binding.duct();
    let count = 128;
    let step = duct.propagator(DT / count as f64, 1e-12).unwrap();
    let mut state = wall.state.clone();
    let mut sum = duct.molecular_dissipated_power_w(&state).unwrap();
    for i in 1..=count {
        step.step(&mut state, &wall_ports()).unwrap();
        let weight = if i == count {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * duct.molecular_dissipated_power_w(&state).unwrap();
    }
    DT * sum / (3.0 * count as f64)
}

fn all_equal(actual: &MixedNodes, expected: &MixedNodes) {
    let FiniteGasAcousticNode::Wall(a) = &actual[&30] else {
        unreachable!()
    };
    let FiniteGasAcousticNode::Wall(b) = &expected[&30] else {
        unreachable!()
    };
    assert_eq!(a.state, b.state);
    assert_eq!(a.propagator.dt(), b.propagator.dt());
    assert_eq!(
        a.binding.duct().operator_identity(),
        b.binding.duct().operator_identity()
    );
    assert_eq!(a.binding.duct().gas_field(), b.binding.duct().gas_field());
    assert_eq!(
        a.binding.duct().molecular_relaxation_law(),
        b.binding.duct().molecular_relaxation_law()
    );
    let without_wall = |map: &MixedNodes| {
        map.iter()
            .filter(|(id, _)| **id != 30)
            .map(|(&id, node)| (id, node.clone()))
            .collect()
    };
    equal_nodes(&without_wall(actual), &without_wall(expected));
}

fn complete_storage(node: &FiniteGasAcousticNode) -> f64 {
    match node {
        FiniteGasAcousticNode::Wall(wall) => wall.state.stored_energy_j(),
        other => node_energy(other),
    }
}

fn work_and_total_loss<S>(
    reports: &BTreeMap<u64, FiniteGasAcousticStep<S>>,
    nodes: &MixedNodes,
    inputs: &MixedInputs,
) -> (f64, f64) {
    let (mut work, mut loss) = (0.0, 0.0);
    for (&id, report) in reports {
        match (report, &nodes[&id]) {
            (FiniteGasAcousticStep::Wall(report), FiniteGasAcousticNode::Wall(_)) => {
                work += report
                    .acoustic
                    .port_outputs
                    .iter()
                    .map(|port| {
                        thermal(&inputs[&id])
                            .ports
                            .iter()
                            .find(|p| p.id == port.id)
                            .unwrap()
                            .effort
                            * port.delta_x
                    })
                    .sum::<f64>();
                loss += report.acoustic_dissipation_j;
            }
            (FiniteGasAcousticStep::Port(report), FiniteGasAcousticNode::Port(load)) => {
                let port = report
                    .acoustic
                    .port_outputs
                    .iter()
                    .find(|p| p.id == load.inlet())
                    .unwrap();
                work += thermal(&inputs[&id]).ports[0].effort * port.delta_x
                    - load
                        .last_terminal_frames()
                        .iter()
                        .map(|t| t.load_work_j)
                        .sum::<f64>();
                loss += report.acoustic_dissipation_j;
            }
            (FiniteGasAcousticStep::Aperture(report), FiniteGasAcousticNode::Aperture(node)) => {
                let frame = &report.coupled;
                work +=
                    frame.upstream_work_j + frame.body_work_j + frame.aperture.mechanical_work_j
                        - node
                            .load()
                            .last_terminal_frames()
                            .iter()
                            .map(|t| t.load_work_j)
                            .sum::<f64>();
                loss += frame.inferred_total_dissipation_j;
            }
            _ => unreachable!(),
        }
    }
    (work, loss)
}

#[test]
fn g1_original_molecular_heat_enters_transported_chemical_boundary_once() {
    let law = reaction(50_000.0);
    let (mut parcels, mut nodes, inputs) = molecular_fixture(&law, true);
    let reactions = reactions(&law);
    let (before_gas, before_nodes) = (parcels.clone(), nodes.clone());
    let heat = original_heat(&nodes);
    let integrated = quadrature_heat(&nodes);
    assert!(heat > 1e-12);
    assert!((heat - integrated).abs() < 1e-6 * integrated + 1e-16);
    let links = links();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let caloric = inputs
        .iter()
        .map(|(&id, input)| (id, thermal(input).gas))
        .collect();
    let transported = material_transport::step_candidate_with_paths(
        &parcels,
        paths,
        DT,
        &caloric,
        solve().transport,
        || false,
    )
    .unwrap();
    let report = advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport.transport, *transported.step());
    assert!(parcels[&30].mass_kg() > before_gas[&30].mass_kg());
    assert_eq!(report.transport.molecular_heat_j, heat);
    let FiniteGasAcousticStep::Wall(wall) = &report.transport.acoustics[&30] else {
        unreachable!()
    };
    assert_eq!(wall.molecular_heat_j, heat);
    assert_eq!(wall.supplied_heat_j, 0.0);
    assert_eq!(
        wall.exported_acoustic_dissipation_j,
        wall.acoustic_dissipation_j - heat
    );
    let chemical = report.reaction(30).unwrap();
    assert_eq!(chemical.boundary.heat_j, heat);
    assert_eq!(chemical.boundary.work_by_gas_j, wall.work_by_gas_j);
    assert_eq!(chemical.reaction.extent_mol, 0.0);
    assert!(report.reaction(20).unwrap().reaction.extent_mol > 1e-9);
    let direct = law
        .boundary_implicit_candidate(
            &transported.parcels()[&30],
            DT,
            chemical.boundary,
            chemistry(),
            || false,
        )
        .unwrap();
    assert_eq!(direct.gas(), &parcels[&30]);
    assert_eq!(direct.step(), chemical);
    let external: f64 = inputs
        .values()
        .map(|input| thermal(input).supplied_heat_j)
        .sum();
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy(gas) - energy(&before_gas[id]))
        .sum();
    let native_change: f64 = nodes
        .iter()
        .map(|(id, node)| complete_storage(node) - complete_storage(&before_nodes[id]))
        .sum();
    let (work, loss) = work_and_total_loss(&report.transport.acoustics, &nodes, &inputs);
    assert!(
        (gas_change + native_change + loss - heat - external - work).abs()
            < solve().energy_tolerance_j
    );
    assert!((report.transport.supplied_heat_j - external).abs() < 1e-18);
    assert!((report.transport.exported_acoustic_dissipation_j - (loss - heat)).abs() < 1e-18);
    assert!(report.maximum_species_residual_mol <= reactive_solve().species_tolerance_mol);
    assert_eq!(nodes[&30].accepted_steps(), 25);
}

#[test]
fn g3_molecular_heat_changes_shared_body_and_zero_source_retains_exact_old_path() {
    let law = reaction(50_000.0);
    let reactions = reactions(&law);
    let (mut parcels, mut nodes, inputs) = molecular_fixture(&law, true);
    let (mut old_gas, mut old_nodes) = (parcels.clone(), nodes.clone());
    let no_heat = export_only(&inputs);
    let (card, key, caloric) = body_material(true);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let (mut body_state, mut old_body) = (initial, initial);
    let mut body = shared_body(caloric.curve(), initial, 0.0);
    body.mass_kg /= 10.0;
    body.volume_m3 /= 10.0;
    let mut ports = shared_ports(&card, &key, &context, [20, 30]);
    for port in ports.values_mut() {
        port.exchange.source_tolerance_j = 1e-14;
    }
    let mut solve = shared_solve();
    solve.body_exchange.source_tolerance_j = 1e-14;
    let links = links();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let before_gas = parcels.clone();
    let before_nodes = nodes.clone();
    let heat = original_heat(&nodes);
    let report = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut body_state,
        body,
        &inputs,
        &reactions,
        &ports,
        paths,
        DT,
        solve,
        &mut || false,
    )
    .unwrap();
    let old = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut old_gas,
        &mut old_nodes,
        &mut old_body,
        body,
        &no_heat,
        &reactions,
        &ports,
        paths,
        DT,
        solve,
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport.molecular_heat_j, heat);
    assert_eq!(old.transport.molecular_heat_j, 0.0);
    assert_eq!(
        old.transport.exported_acoustic_dissipation_j,
        old.transport.acoustic_dissipation_j
    );
    assert!(parcels[&30].state().state().temperature > old_gas[&30].state().state().temperature);
    // Linear constant-capacity feedback gives a useful independent signal scale:
    // beta*M/(C_body + alpha_30 + alpha_20). The small conductivity and reaction
    // changes here do not justify exact linear equality, but the effect must be
    // resolved well above the declared root allowance, not merely differ by an ulp.
    let gas_capacity = 1.5 * R * parcels[&30].mass_kg() / 0.004;
    let body_capacity = body.mass_kg * BODY_CP;
    let a = DT * (AREA / 2.0) / (GAS_DISTANCE / 0.026 + BODY_DISTANCE / 0.2);
    let expected_scale = (a / (gas_capacity + a)) * heat / (body_capacity + 2.0 * a);
    assert!(expected_scale > 16.0 * solve.body_exchange.source_tolerance_j / body_capacity);
    assert!(body_state.temperature_k() - old_body.temperature_k() > 0.25 * expected_scale);
    let contact = report.port(30).unwrap().solid_to_gas_heat_j;
    assert!((report.reaction(30).unwrap().boundary.heat_j - (heat + contact)).abs() < 1e-18);
    assert!(report.reaction(20).unwrap().reaction.extent_mol > 1e-9);
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy(gas) - energy(&before_gas[id]))
        .sum();
    let native_change: f64 = nodes
        .iter()
        .map(|(id, node)| complete_storage(node) - complete_storage(&before_nodes[id]))
        .sum();
    let body_change =
        body.mass_kg * BODY_CP * (body_state.temperature_k() - initial.temperature_k());
    let external: f64 = inputs
        .values()
        .map(|input| thermal(input).supplied_heat_j)
        .sum();
    let (work, loss) = work_and_total_loss(&report.transport.acoustics, &nodes, &inputs);
    assert!(
        (gas_change + native_change + body_change + loss - heat - external - work).abs()
            < solve.reactive.transport.energy_tolerance_j
    );
    let next_heat = original_heat(&nodes);
    body.initial = body_state;
    let next = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut body_state,
        body,
        &inputs,
        &reactions,
        &ports,
        paths,
        DT,
        solve,
        &mut || false,
    )
    .unwrap();
    assert_eq!(next.transport.molecular_heat_j, next_heat);
    for id in [10, 20, 30] {
        assert_eq!(nodes[&id].accepted_steps(), 26);
    }

    // The explicitly empty constitutive source is exactly the former export-only
    // path: no changed gas, body, acoustic history or source work from opting in.
    let (mut zero_gas, mut zero_nodes, zero_input) = molecular_fixture(&law, false);
    let (mut plain_gas, mut plain_nodes) = (zero_gas.clone(), zero_nodes.clone());
    let (mut zero_body, mut plain_body) = (initial, initial);
    body.initial = initial;
    let zero = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut zero_gas,
        &mut zero_nodes,
        &mut zero_body,
        body,
        &zero_input,
        &reactions,
        &ports,
        paths,
        DT,
        solve,
        &mut || false,
    )
    .unwrap();
    let plain = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut plain_gas,
        &mut plain_nodes,
        &mut plain_body,
        body,
        &export_only(&zero_input),
        &reactions,
        &ports,
        paths,
        DT,
        solve,
        &mut || false,
    )
    .unwrap();
    assert_eq!(zero.transport.molecular_heat_j, 0.0);
    assert_eq!(zero_gas, plain_gas);
    assert_eq!(zero_body, plain_body);
    all_equal(&zero_nodes, &plain_nodes);
    assert_eq!(zero.transport.transport, plain.transport.transport);
    assert_eq!(
        zero.combined_energy_residual_j,
        plain.combined_energy_residual_j
    );
    assert_eq!(zero.endpoint_evaluations, plain.endpoint_evaluations);
}

#[test]
fn g4_molecular_source_refusals_and_final_cancellation_preserve_all_owners() {
    let law = reaction(50_000.0);
    let reactions = reactions(&law);
    let (mut parcels, mut nodes, inputs) = molecular_fixture(&law, true);
    let (card, key, caloric) = body_material(true);
    let context = body_point(300.0);
    let mut state = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let body = shared_body(caloric.curve(), state, 0.0);
    let ports = shared_ports(&card, &key, &context, [20, 30]);
    let links = links();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let (before_gas, before_nodes, before_body) = (parcels.clone(), nodes.clone(), state);
    let (mut accepted_gas, mut accepted_nodes, mut accepted_body) =
        (parcels.clone(), nodes.clone(), state);
    let mut polls = 0;
    let accepted = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut accepted_gas,
        &mut accepted_nodes,
        &mut accepted_body,
        body,
        &inputs,
        &reactions,
        &ports,
        paths,
        DT,
        shared_solve(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut visited = 0;
    let cancelled = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut state,
        body,
        &inputs,
        &reactions,
        &ports,
        paths,
        DT,
        shared_solve(),
        &mut || {
            visited += 1;
            visited == polls
        },
    );
    assert!(matches!(
        cancelled,
        Err(ThermalReactiveFiniteGasAcousticTransportError::Cancelled)
    ));
    assert_eq!(parcels, before_gas);
    all_equal(&nodes, &before_nodes);
    assert_eq!(state, before_body);
    for bad in 0..3 {
        let mut wrong = inputs.clone();
        match bad {
            0 => {
                wrong.insert(
                    10,
                    FiniteGasAcousticInput::MolecularWall {
                        thermal: thermal(&inputs[&10]).clone(),
                        heat_budget: heat_budget(),
                    },
                );
            }
            1 => {
                let FiniteGasAcousticInput::MolecularWall { heat_budget, .. } =
                    wrong.get_mut(&30).unwrap()
                else {
                    unreachable!()
                };
                heat_budget.max_dimension = 1;
            }
            _ => {
                thermal_mut(wrong.get_mut(&30).unwrap())
                    .acoustic
                    .as_mut()
                    .unwrap()
                    .loss = FiniteGasDgLoss::Lossless;
            }
        }
        assert!(
            advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
                &mut parcels,
                &mut nodes,
                &mut state,
                body,
                &wrong,
                &reactions,
                &ports,
                paths,
                DT,
                shared_solve(),
                &mut || false,
            )
            .is_err()
        );
        assert_eq!(parcels, before_gas);
        all_equal(&nodes, &before_nodes);
        assert_eq!(state, before_body);
    }
    // The source integral is not authority to invent acoustic species-history
    // transfer: different incoming fractions still fail at the original owner.
    let mut changed_mix = parcels.clone();
    changed_mix.insert(
        90,
        ClosedGasParcel::new(
            law.catalog().clone(),
            point(300.0, 180_000.0),
            parcels[&90].volume_m3(),
        )
        .unwrap(),
    );
    let changed_before = changed_mix.clone();
    assert!(
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut changed_mix,
            &mut nodes,
            &mut state,
            body,
            &inputs,
            &reactions,
            &ports,
            paths,
            DT,
            shared_solve(),
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(changed_mix, changed_before);
    all_equal(&nodes, &before_nodes);
    assert_eq!(state, before_body);
    let retry = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut state,
        body,
        &inputs,
        &reactions,
        &ports,
        paths,
        DT,
        shared_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, accepted_gas);
    all_equal(&nodes, &accepted_nodes);
    assert_eq!(state, accepted_body);
    assert_eq!(
        retry.transport.molecular_heat_j,
        accepted.transport.molecular_heat_j
    );
    assert_eq!(
        retry.combined_energy_residual_j,
        accepted.combined_energy_residual_j
    );
    assert_eq!(retry.endpoint_evaluations, accepted.endpoint_evaluations);
}
