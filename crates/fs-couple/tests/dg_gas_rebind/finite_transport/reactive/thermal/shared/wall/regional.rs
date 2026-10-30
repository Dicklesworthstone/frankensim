//! One thermal inventory preserves distinct original DG wall regions.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{RegionalMaterialWallDg, WallRelaxation};
use fs_couple::dg_port::finite_gas::transport::assembly::reactive::thermal::shared::{
    SharedThermalReactiveFiniteGasAcousticTransportCandidate,
    SharedThermalReactiveFiniteGasAcousticTransportStep,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_qty::DynViscosity;

const ARM: f64 = 4e7;
const ARM_SLOPE: f64 = 0.002;
const BAND: [f64; 2] = [20.0, 100_000.0];

#[derive(Clone, Copy)]
enum Loss {
    Elastic,
    Kelvin,
    Maxwell,
}

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn source_card(viscosity_upper: f64) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let (base, key, _) = material(0.0, 800.0, false);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.validity = claim.validity.with("lot", 7.0, 8.0);
        claims.insert_claim(claim).unwrap();
    }
    for (name, dims, values, upper) in [
        (
            "arm_modulus",
            Pressure::DIMS,
            [0.8 * ARM, ARM, 2.0 * ARM],
            800.0,
        ),
        ("arm_time", Dims([0, 0, 1, 0, 0, 0]), [4e-5; 3], 800.0),
        (
            "young_viscosity",
            DynViscosity::DIMS,
            [70.0, 100.0, 250.0],
            viscosity_upper,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [200.0, 300.0, 800.0].into_iter().zip(values).collect(),
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        200.0,
                        upper,
                    )
                    .with("lot", 7.0, 8.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic regional wall heat fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "one synthetic regional wall inventory".into(),
            phase: "solid".into(),
            process: "fixed density, differing regional geometry and loss".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &body_point(200.0),
            upper: &body_point(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap();
    (card, key, caloric)
}

fn local_wall(
    card: &MaterialCard,
    temperature: f64,
    lot: f64,
    thickness: f64,
    loss: Loss,
) -> LocalMaterialWall {
    let mut properties = vec![
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("poisson_ratio", Dims::NONE),
    ];
    match loss {
        Loss::Elastic => {}
        Loss::Kelvin => properties.push(("young_viscosity", DynViscosity::DIMS)),
        Loss::Maxwell => properties.extend([
            ("arm_modulus", Pressure::DIMS),
            ("arm_time", Dims([0, 0, 1, 0, 0, 0])),
        ]),
    }
    let requirements = properties
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_with_quantity(
                name,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let material = resolve_material_state_point(
        card,
        &body_point(temperature).with("lot", lot).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let loss = match loss {
        Loss::Elastic => WallMaterialLoss::Elastic,
        Loss::Kelvin => WallMaterialLoss::KelvinVoigt,
        Loss::Maxwell => WallMaterialLoss::Maxwell(
            WallRelaxation::from_state_point(
                &material,
                &[("arm_modulus".into(), "arm_time".into())],
                BAND,
                1.0,
                0.8,
            )
            .unwrap(),
        ),
    };
    LocalMaterialWall::new(material, thickness, AxialConstraint::StressFree, loss, BAND).unwrap()
}

fn regional_fixture(
    law: &GasReaction,
    walls: Vec<LocalMaterialWall>,
    uniform: bool,
) -> (Parcels, MixedNodes, MixedInputs) {
    let field = field_for(law);
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &field,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap();
    let binding = if uniform {
        assert_eq!(walls[0], walls[1]);
        FiniteGasDgBinding::Uniform(walls[0].clone().bind_dg(duct).unwrap())
    } else {
        FiniteGasDgBinding::Regional(RegionalMaterialWallDg::new(walls, duct).unwrap())
    };
    let propagator = binding.duct().propagator(DT, 1e-11).unwrap();
    let mut state = binding
        .duct()
        .project_initial_state(|x| 1000.0 + 3000.0 * (1.0 - 2.0 * x / LENGTH), |_| 1e-9);
    for _ in 0..24 {
        propagator.step(&mut state, &held_ports()).unwrap();
    }
    propagator.sync_nodal(&mut state).unwrap();
    let gas = ClosedGasParcel::new(
        field.model().clone(),
        point(300.0, 101_325.0),
        core::f64::consts::PI * RADIUS * RADIUS * LENGTH,
    )
    .unwrap();
    (
        [(30, gas)].into(),
        [(
            30,
            FiniteGasAcousticNode::Wall(FiniteGasDgNode {
                binding,
                state,
                propagator,
            }),
        )]
        .into(),
        [(
            30,
            FiniteGasAcousticInput::Passive(FiniteGasDgTransportInput {
                gas: energy_solve(),
                supplied_heat_j: 0.0,
                acoustic: Some(controls()),
                ports: held_ports(),
            }),
        )]
        .into(),
    )
}

fn regions(nodes: &MixedNodes) -> &[LocalMaterialWall] {
    let FiniteGasDgBinding::Regional(wall) = &wall_node(nodes).binding else {
        unreachable!()
    };
    wall.walls()
}

fn equal_regional(actual: &MixedNodes, expected: &MixedNodes) {
    let (a, b) = (wall_node(actual), wall_node(expected));
    assert_eq!(a.state, b.state);
    assert_eq!(a.propagator.dt(), b.propagator.dt());
    assert_eq!(a.propagator.tolerance(), b.propagator.tolerance());
    assert_eq!(
        a.binding.duct().operator_identity(),
        b.binding.duct().operator_identity()
    );
    assert_eq!(a.binding.duct().gas_field(), b.binding.duct().gas_field());
    assert_eq!(regions(actual), regions(expected));
}

fn thermal_body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    initial: EquilibriumPhaseState,
    mean_thickness: f64,
) -> UniformEnthalpyStepInput<'_> {
    let mut body = wall_body(curve, initial, 0.0);
    body.volume_m3 *= mean_thickness / THICKNESS;
    body.mass_kg *= mean_thickness / THICKNESS;
    body
}

fn preview(
    gas: &Parcels,
    nodes: &MixedNodes,
    body: UniformEnthalpyStepInput<'_>,
    inputs: &MixedInputs,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
) -> SharedThermalReactiveFiniteGasAcousticTransportCandidate {
    finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths(
        gas,
        nodes,
        &body.initial,
        body,
        inputs,
        &BTreeMap::new(),
        ports,
        GasTransportPaths::default(),
        DT,
        shared_solve(),
        &mut || false,
    )
    .unwrap()
}

fn advance(
    gas: &mut Parcels,
    nodes: &mut MixedNodes,
    state: &mut EquilibriumPhaseState,
    body: UniformEnthalpyStepInput<'_>,
    inputs: &MixedInputs,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    SharedThermalReactiveFiniteGasAcousticTransportStep,
    ThermalReactiveFiniteGasAcousticTransportError,
> {
    advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        gas,
        nodes,
        state,
        body,
        inputs,
        &BTreeMap::new(),
        ports,
        GasTransportPaths::default(),
        DT,
        shared_solve(),
        cancelled,
    )
}

#[test]
fn g1_one_body_preserves_distinct_regional_material_heat_and_maxwell_history() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = source_card(800.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let original = vec![
        local_wall(&card, 300.0, 7.0, THICKNESS, Loss::Maxwell),
        local_wall(&card, 300.0, 8.0, 2.0 * THICKNESS, Loss::Kelvin),
    ];
    let (gas, nodes, inputs) = regional_fixture(&law, original.clone(), false);
    let context = body_point(300.0);
    let mut ports = wall_ports(&card, &key, &context, 0.0);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(heat_budget());
    let mut body = thermal_body(caloric.curve(), initial, 1.5 * THICKNESS);
    let capacity = body.mass_kg * BODY_CP;
    body.internal_heat_j = capacity * 20.0;
    let owner = wall_node(&nodes);
    let heat = owner
        .binding
        .duct()
        .prepare_compliant_wall_heat(&owner.propagator, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&owner.state, &held_ports(), &mut || false)
        .unwrap();
    assert!(heat > 0.0);
    let frozen = finite_gas_reacting_acoustic_transport_candidate_with_paths(
        &gas,
        &nodes,
        &inputs,
        &BTreeMap::new(),
        GasTransportPaths::default(),
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    let before = wall_node(frozen.nodes())
        .binding
        .duct()
        .compliant_wall_state(&wall_node(frozen.nodes()).state)
        .unwrap();
    assert_eq!(before.displacement_m.len(), 6);
    assert_eq!(before.relaxation_pressure_pa.len(), 1);
    let arm_k = ARM * THICKNESS / RADIUS.powi(2);
    let extension: Vec<_> = before.relaxation_pressure_pa[0][..3]
        .iter()
        .map(|p| p / arm_k)
        .collect();
    let arm_energy = 0.5 * spring_energy(&extension, ARM);
    assert!((arm_energy - before.relaxation_stored_energy_j).abs() < 1e-13);
    let equilibrium = 0.5 * spring_energy(&before.displacement_m[..3], YOUNG_300)
        + spring_energy(&before.displacement_m[3..], YOUNG_300);
    let work_per_k = equilibrium * YOUNG_SLOPE + arm_energy * ARM_SLOPE;
    let expected_delta = (body.internal_heat_j + heat) / (capacity + work_per_k);
    let candidate = preview(&gas, &nodes, body, &inputs, &ports);
    let report = candidate.step();
    let t = candidate.body_state().temperature_k();
    let update = &report.wall_updates[&30];
    assert!(update.uniform_material().is_none());
    let receipts = update.regional_materials().unwrap();
    assert_eq!(receipts.len(), 2);
    assert!((t - 300.0 - expected_delta).abs() < 1e-9);
    assert_eq!(update.dissipated_heat_j, heat);
    assert_eq!(report.wall_heat_j, heat);
    assert!((report.wall_parameter_work_j - work_per_k * (t - 300.0)).abs() < 1e-13);
    for (i, receipt) in receipts.iter().enumerate() {
        assert_eq!(receipt, regions(candidate.nodes())[i].material());
        assert_eq!(receipt.query().axes()["lot"], 7.0 + i as f64);
        assert_eq!(receipt.query().axes()["T"], t);
        for name in ["density", "young_modulus", "poisson_ratio"] {
            assert_eq!(
                receipt.property(name).unwrap().answer().receipt.selected,
                original[i]
                    .material()
                    .property(name)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    let accepted = wall_node(candidate.nodes());
    let after = accepted
        .binding
        .duct()
        .compliant_wall_state(&accepted.state)
        .unwrap();
    for i in 0..6 {
        assert!((after.displacement_m[i] - before.displacement_m[i]).abs() < 1e-15);
        assert!((after.velocity_m_s[i] - before.velocity_m_s[i]).abs() < 1e-12);
        if i < 3 {
            assert!(
                (after.relaxation_pressure_pa[0][i] / (arm_k * (1.0 + ARM_SLOPE * (t - 300.0)))
                    - extension[i])
                    .abs()
                    < 1e-15
            );
        } else {
            assert_eq!(after.relaxation_pressure_pa[0][i], 0.0);
        }
    }
    let body_storage = capacity * (t - 300.0);
    let acoustic_storage = accepted.state.stored_energy_j() - owner.state.stored_energy_j();
    let balance = energy(&candidate.parcels()[&30]) - energy(&gas[&30])
        + acoustic_storage
        + body_storage
        + report.final_exported_acoustic_dissipation_j
        - body.internal_heat_j
        - report.transport.port_work_j;
    assert!(balance.abs() < shared_solve().reactive.transport.energy_tolerance_j);
    assert!(
        (body_storage + report.wall_parameter_work_j - body.internal_heat_j - heat).abs() < 1e-11
    );
    let next = preview(
        candidate.parcels(),
        candidate.nodes(),
        UniformEnthalpyStepInput {
            initial: *candidate.body_state(),
            internal_heat_j: 0.0,
            ..body
        },
        &inputs,
        &ports,
    );
    assert_eq!(
        wall_node(next.nodes()).state.step_count(),
        accepted.state.step_count() + 1
    );
    assert_eq!(
        next.step().wall_updates[&30]
            .regional_materials()
            .unwrap()
            .len(),
        2
    );
    assert!(next.step().wall_heat_j > 0.0);
}

#[test]
fn g3_homogeneous_regional_binding_retains_the_uniform_trajectory_and_receipt_layout() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = source_card(800.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let source = local_wall(&card, 300.0, 7.0, THICKNESS, Loss::Elastic);
    let walls = vec![source.clone(), source];
    let (gas_u, nodes_u, inputs_u) = regional_fixture(&law, walls.clone(), true);
    let (gas_r, nodes_r, inputs_r) = regional_fixture(&law, walls, false);
    assert_eq!(wall_node(&nodes_u).state, wall_node(&nodes_r).state);
    let context = body_point(300.0);
    let mut ports = wall_ports(&card, &key, &context, 0.0);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(heat_budget());
    let mut body = thermal_body(caloric.curve(), initial, THICKNESS);
    body.internal_heat_j = body.mass_kg * BODY_CP * 20.0;
    let u = preview(&gas_u, &nodes_u, body, &inputs_u, &ports);
    let r = preview(&gas_r, &nodes_r, body, &inputs_r, &ports);
    assert_eq!(u.parcels(), r.parcels());
    assert_eq!(u.body_state(), r.body_state());
    assert_eq!(wall_node(u.nodes()).state, wall_node(r.nodes()).state);
    assert_eq!(
        wall_node(u.nodes()).binding.duct().operator_identity(),
        wall_node(r.nodes()).binding.duct().operator_identity()
    );
    assert_eq!(
        u.step().wall_parameter_work_j,
        r.step().wall_parameter_work_j
    );
    assert_eq!(u.step().wall_heat_j, 0.0);
    assert_eq!(r.step().wall_heat_j, 0.0);
    let uniform = u.step().wall_updates[&30].uniform_material().unwrap();
    assert!(u.step().wall_updates[&30].regional_materials().is_none());
    assert_eq!(
        r.step().wall_updates[&30].regional_materials().unwrap(),
        &[uniform.clone(), uniform.clone()]
    );
}

#[test]
fn g4_later_regional_source_inventory_and_cancellation_fail_without_partial_publication() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = source_card(800.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let walls = vec![
        local_wall(&card, 300.0, 7.0, THICKNESS, Loss::Maxwell),
        local_wall(&card, 300.0, 8.0, 2.0 * THICKNESS, Loss::Kelvin),
    ];
    let (mut gas, mut nodes, inputs) = regional_fixture(&law, walls, false);
    let (old_gas, old_nodes) = (gas.clone(), nodes.clone());
    let context = body_point(300.0);
    let mut ports = wall_ports(&card, &key, &context, 0.0);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(heat_budget());
    let mut body = thermal_body(caloric.curve(), initial, 1.5 * THICKNESS);
    body.internal_heat_j = body.mass_kg * BODY_CP * 20.0;
    let mut state = initial;
    let (mut accepted_gas, mut accepted_nodes, mut accepted_state) =
        (gas.clone(), nodes.clone(), state);
    let mut polls = 0;
    let accepted = advance(
        &mut accepted_gas,
        &mut accepted_nodes,
        &mut accepted_state,
        body,
        &inputs,
        &ports,
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut visited = 0;
    assert!(matches!(
        advance(
            &mut gas,
            &mut nodes,
            &mut state,
            body,
            &inputs,
            &ports,
            &mut || {
                visited += 1;
                visited == polls
            }
        ),
        Err(ThermalReactiveFiniteGasAcousticTransportError::Cancelled)
    ));
    assert_eq!(visited, polls);
    assert_eq!(gas, old_gas);
    equal_regional(&nodes, &old_nodes);
    assert_eq!(state, initial);
    // First region fits; the second's actual extra thickness exceeds inventory.
    let bad_body = UniformEnthalpyStepInput {
        volume_m3: body.volume_m3 * 0.5,
        mass_kg: body.mass_kg * 0.5,
        ..body
    };
    assert!(
        advance(
            &mut gas,
            &mut nodes,
            &mut state,
            bad_body,
            &inputs,
            &ports,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(gas, old_gas);
    equal_regional(&nodes, &old_nodes);
    assert_eq!(state, initial);
    // The first Maxwell source remains valid at 320 K. Only the second region's
    // own Kelvin viscosity source refuses the private heated endpoint.
    let (limited, limited_key, limited_caloric) = source_card(315.0);
    let limited_initial = limited_caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let (mut limited_gas, mut limited_nodes, limited_inputs) = regional_fixture(
        &law,
        vec![
            local_wall(&limited, 300.0, 7.0, THICKNESS, Loss::Maxwell),
            local_wall(&limited, 300.0, 8.0, 2.0 * THICKNESS, Loss::Kelvin),
        ],
        false,
    );
    let (before_gas, before_nodes) = (limited_gas.clone(), limited_nodes.clone());
    let mut limited_state = limited_initial;
    let mut limited_body = thermal_body(limited_caloric.curve(), limited_initial, 1.5 * THICKNESS);
    limited_body.internal_heat_j = limited_body.mass_kg * BODY_CP * 20.0;
    let limited_ports = wall_ports(&limited, &limited_key, &context, 0.0);
    assert!(
        advance(
            &mut limited_gas,
            &mut limited_nodes,
            &mut limited_state,
            limited_body,
            &limited_inputs,
            &limited_ports,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(limited_gas, before_gas);
    equal_regional(&limited_nodes, &before_nodes);
    assert_eq!(limited_state, limited_initial);
    let retry = advance(
        &mut gas,
        &mut nodes,
        &mut state,
        body,
        &inputs,
        &ports,
        &mut || false,
    )
    .unwrap();
    assert_eq!(gas, accepted_gas);
    equal_regional(&nodes, &accepted_nodes);
    assert_eq!(state, accepted_state);
    assert_eq!(retry.wall_updates, accepted.wall_updates);
}
