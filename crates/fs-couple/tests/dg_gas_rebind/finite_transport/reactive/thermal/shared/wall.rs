//! Finite body enthalpy pays the original vibrating wall's elastic material work.
//! These authored laws are regression fixtures, not measured material data.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, FiniteGasDgBinding, FiniteGasDgNode, LocalMaterialWall, WallMaterialLoss,
};
use fs_couple::dg_port::finite_gas::transport::assembly::reactive::thermal::shared::SharedFiniteGasAcousticWallSolve;
use fs_duct::dg::PortInput;
use fs_material::phase::EnthalpyPhaseKnot;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::Dims;

const WALL_DENSITY: f64 = 4000.0;
const THICKNESS: f64 = 2e-5;
const YOUNG_300: f64 = 1e8;
const YOUNG_SLOPE: f64 = -0.001;

fn young(t: f64) -> f64 {
    YOUNG_300 * (1.0 + YOUNG_SLOPE * (t - 300.0))
}

fn material(
    density_slope: f64,
    elastic_max_temperature: f64,
    competing_claim: bool,
) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let elastic = PropertyKey::new("young_modulus", Pressure::DIMS);
    let mut claims = ClaimSet::new();
    for (key, values, upper) in [
        (cp.clone(), [BODY_CP; 3], 800.0),
        (
            density.clone(),
            [
                WALL_DENSITY - 100.0 * density_slope,
                WALL_DENSITY,
                WALL_DENSITY + 500.0 * density_slope,
            ],
            800.0,
        ),
        (conductivity.clone(), [0.2; 3], 800.0),
        (
            elastic.clone(),
            [young(200.0), young(300.0), young(800.0)],
            elastic_max_temperature,
        ),
        (
            PropertyKey::new("poisson_ratio", Dims::NONE),
            [0.3; 3],
            800.0,
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
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
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite elastic thermal inventory".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    if competing_claim {
        claims
            .insert_claim(PropertyClaim {
                key: elastic,
                value: PropertyValue::Scalar {
                    value: 3.0 * YOUNG_300,
                    dims: Pressure::DIMS,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        310.0,
                        800.0,
                    )
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "deliberately competing hot-wall claim; must not replace selection"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic elastic thermal body".into(),
            phase: "solid".into(),
            process: "constant heat capacity, temperature-dependent stiffness".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &density,
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
    (card, conductivity, caloric)
}

fn held_ports() -> Vec<PortInput> {
    vec![PortInput {
        id: PortId(7),
        effort: 1e-9,
    }]
}

fn wall_fixture(law: &GasReaction, card: &MaterialCard) -> (Parcels, MixedNodes, MixedInputs) {
    let requirements = [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("poisson_ratio", Dims::NONE),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_with_quantity(
            name,
            QuantitySpec::dimensional(dims),
            ScalarAdmissibility::Finite,
        )
        .unwrap()
    });
    let resolved = resolve_material_state_point(
        card,
        &body_point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
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
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap();
    let wall = LocalMaterialWall::new(
        resolved,
        THICKNESS,
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
                binding: FiniteGasDgBinding::Uniform(wall),
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

fn wall_body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    initial: EquilibriumPhaseState,
    heat_j: f64,
) -> UniformEnthalpyStepInput<'_> {
    let volume = core::f64::consts::TAU * RADIUS * LENGTH * THICKNESS;
    UniformEnthalpyStepInput {
        curve,
        initial,
        mass_kg: WALL_DENSITY * volume,
        volume_m3: volume,
        surface_area_m2: AREA,
        internal_heat_j: heat_j,
        duration_s: DT,
    }
}

fn wall_ports<'a>(
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    context: &'a QueryPoint,
    area_m2: f64,
) -> BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'a>> {
    let mut ports = shared_ports(card, key, context, [20, 30]);
    ports.remove(&20);
    let port = ports.get_mut(&30).unwrap();
    port.layers.area_m2 = area_m2;
    port.wall = Some(SharedFiniteGasAcousticWallSolve {
        maximum_absolute_parameter_work_j: 1e-4,
    });
    ports
}

fn wall_node(nodes: &MixedNodes) -> &FiniteGasDgNode {
    let FiniteGasAcousticNode::Wall(node) = &nodes[&30] else {
        unreachable!()
    };
    node
}

fn wall_source(nodes: &MixedNodes) -> &LocalMaterialWall {
    let FiniteGasDgBinding::Uniform(wall) = &wall_node(nodes).binding else {
        unreachable!()
    };
    wall.wall()
}

fn same_owners(actual: &MixedNodes, expected: &MixedNodes) {
    let (a, b) = (wall_node(actual), wall_node(expected));
    assert_eq!(a.state, b.state);
    assert_eq!(a.propagator.dt(), b.propagator.dt());
    assert_eq!(a.propagator.tolerance(), b.propagator.tolerance());
    assert_eq!(
        a.binding.duct().operator_identity(),
        b.binding.duct().operator_identity()
    );
    assert_eq!(a.binding.duct().gas_field(), b.binding.duct().gas_field());
    assert_eq!(wall_source(actual), wall_source(expected));
    let others = |nodes: &MixedNodes| {
        nodes
            .iter()
            .filter(|(id, _)| **id != 30)
            .map(|(&id, node)| (id, node.clone()))
            .collect()
    };
    equal_nodes(&others(actual), &others(expected));
}

// Independent exact Gauss integral of the degree-two physical displacement,
// with stress-free membrane stiffness E*t/r^2 and surface measure 2*pi*r*dx.
fn spring_energy(displacement: &[f64], modulus: f64) -> f64 {
    assert_eq!(displacement.len(), 3);
    let a = (3.0_f64 / 5.0).sqrt();
    let square_integral = [(-a, 5.0 / 9.0), (0.0, 8.0 / 9.0), (a, 5.0 / 9.0)]
        .into_iter()
        .map(|(x, weight)| {
            let u = 0.5 * x * (x - 1.0) * displacement[0]
                + (1.0 - x * x) * displacement[1]
                + 0.5 * x * (x + 1.0) * displacement[2];
            weight * u * u
        })
        .sum::<f64>()
        * LENGTH
        / 2.0;
    0.5 * modulus * THICKNESS / (RADIUS * RADIUS)
        * core::f64::consts::TAU
        * RADIUS
        * square_integral
}

#[test]
fn g1_finite_body_pays_signed_elastic_work_without_inventing_contact_heat() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = material(0.0, 800.0, true);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let (parcels, nodes, inputs) = wall_fixture(&law, &card);
    let ports = wall_ports(&card, &key, &context, 0.0);
    // Independently march the accepted frozen native owner once. The thermal
    // intervention has zero duration and must retain this physical wall motion.
    let mut native = wall_node(&nodes).clone();
    native
        .propagator
        .step(&mut native.state, &held_ports())
        .unwrap();
    native.propagator.sync_nodal(&mut native.state).unwrap();
    let motion = native
        .binding
        .duct()
        .compliant_wall_state(&native.state)
        .unwrap();
    let spring = spring_energy(&motion.displacement_m, YOUNG_300);
    assert!(spring > 1e-12);
    for delta_t in [-20.0, 20.0] {
        let mut body = wall_body(caloric.curve(), initial, 0.0);
        let capacity = body.mass_kg * BODY_CP;
        body.internal_heat_j = capacity * delta_t;
        let candidate = finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths(
            &parcels,
            &nodes,
            &initial,
            body,
            &inputs,
            &BTreeMap::new(),
            &ports,
            GasTransportPaths::default(),
            DT,
            shared_solve(),
            &mut || false,
        )
        .unwrap();
        let report = candidate.step();
        let expected_delta = body.internal_heat_j / (capacity + YOUNG_SLOPE * spring);
        let expected_work = YOUNG_SLOPE * expected_delta * spring;
        let updated = wall_node(candidate.nodes());
        let final_motion = updated
            .binding
            .duct()
            .compliant_wall_state(&updated.state)
            .unwrap();
        let update = &report.wall_updates[&30];
        let temperature = candidate.body_state().temperature_k();
        let body_storage = capacity * (temperature - 300.0);
        let acoustic_storage =
            updated.state.stored_energy_j() - wall_node(&nodes).state.stored_energy_j();
        assert!((temperature - 300.0 - expected_delta).abs() < 1e-9);
        assert!(report.wall_parameter_work_j * delta_t < 0.0);
        assert!((report.wall_parameter_work_j - expected_work).abs() < 1e-13);
        assert!(
            (update.native.wall_energy_after_j
                - update.native.wall_energy_before_j
                - expected_work)
                .abs()
                < 1e-13
        );
        assert!(
            (final_motion.stored_energy_j - motion.stored_energy_j - expected_work).abs() < 1e-12
        );
        assert!((body_storage + expected_work - body.internal_heat_j).abs() < 2e-11);
        assert_eq!(report.port(30).unwrap().solid_to_gas_heat_j, 0.0);
        assert_eq!(report.solid_to_gas_heat_j, 0.0);
        assert!(report.body_boundary_heat_j.abs() < 1e-13);
        assert_eq!(report.final_acoustic_storage_change_j, acoustic_storage);
        assert!(
            (report.transport.acoustic_storage_change_j + expected_work - acoustic_storage).abs()
                < 1e-12
        );
        let gas_storage = energy(&candidate.parcels()[&30]) - energy(&parcels[&30]);
        let global = gas_storage
            + acoustic_storage
            + body_storage
            + report.transport.exported_acoustic_dissipation_j
            - body.internal_heat_j
            - report.transport.port_work_j;
        assert!(global.abs() < shared_solve().reactive.transport.energy_tolerance_j);
        for (before, after) in motion
            .displacement_m
            .iter()
            .zip(&final_motion.displacement_m)
        {
            assert!((before - after).abs() < 1e-15);
        }
        for (before, after) in motion.velocity_m_s.iter().zip(&final_motion.velocity_m_s) {
            assert!((before - after).abs() < 1e-12);
        }
        assert_eq!(updated.state.step_count(), native.state.step_count());
        assert_eq!(updated.propagator.dt(), native.propagator.dt());
        assert_eq!(update.material, *wall_source(candidate.nodes()).material());
        assert_eq!(
            update.material.query_point(),
            body_point(temperature)
                .axes()
                .iter()
                .map(|(a, &v)| (a.clone(), v))
                .collect::<Vec<_>>()
        );
        for name in ["density", "young_modulus", "poisson_ratio"] {
            assert_eq!(
                update
                    .material
                    .property(name)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                wall_source(&nodes)
                    .material()
                    .property(name)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        assert!(
            (update
                .material
                .property("young_modulus")
                .unwrap()
                .value_si()
                - young(temperature))
            .abs()
                < 1e-6
        );
        assert_eq!(update.material.card_identity(), card.content_hash());
        if delta_t > 0.0 {
            assert_eq!(
                update
                    .material
                    .property("young_modulus")
                    .unwrap()
                    .answer()
                    .receipt
                    .in_domain
                    .len(),
                2
            );
        }
    }
}

#[test]
fn g3_refreshed_wall_continues_with_shared_reacting_gas_and_none_retains_native_path() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = material(0.0, 800.0, false);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let (wall_gas, wall_nodes, wall_inputs) = wall_fixture(&law, &card);
    let mut disconnected = wall_ports(&card, &key, &context, 0.0);
    disconnected.get_mut(&30).unwrap().wall = None;
    let no_wall = finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths(
        &wall_gas,
        &wall_nodes,
        &initial,
        wall_body(caloric.curve(), initial, 0.0),
        &wall_inputs,
        &BTreeMap::new(),
        &disconnected,
        GasTransportPaths::default(),
        DT,
        shared_solve(),
        &mut || false,
    )
    .unwrap();
    let old = finite_gas_reacting_acoustic_transport_candidate_with_paths(
        &wall_gas,
        &wall_nodes,
        &wall_inputs,
        &BTreeMap::new(),
        GasTransportPaths::default(),
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(no_wall.parcels(), old.parcels());
    same_owners(no_wall.nodes(), old.nodes());
    assert_eq!(*no_wall.body_state(), initial);
    assert!(no_wall.step().wall_updates.is_empty());
    assert_eq!(no_wall.step().wall_parameter_work_j, 0.0);
    assert_eq!(
        no_wall.step().final_acoustic_storage_change_j,
        no_wall.step().transport.acoustic_storage_change_j
    );

    let (mut parcels, mut nodes, mut inputs) = fixture(&law);
    parcels.extend(wall_gas);
    nodes.extend(wall_nodes);
    inputs.extend(wall_inputs);
    let (mut frozen_gas, mut frozen_nodes, mut frozen_body) =
        (parcels.clone(), nodes.clone(), initial);
    let mut state = initial;
    let mut ports = shared_ports(&card, &key, &context, [20, 30]);
    let frozen_ports = shared_ports(&card, &key, &context, [20, 30]);
    ports.get_mut(&30).unwrap().wall = Some(SharedFiniteGasAcousticWallSolve {
        maximum_absolute_parameter_work_j: 1e-4,
    });
    let mut reactions = bindings(&law);
    reactions.insert(
        30,
        FiniteGasAcousticReaction {
            reaction: &law,
            controls: chemistry(),
        },
    );
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    for sample in 0..2 {
        let mut body = wall_body(caloric.curve(), state, 0.0);
        body.internal_heat_j = if sample == 0 {
            body.mass_kg * BODY_CP * 20.0
        } else {
            0.0
        };
        let frozen_input = UniformEnthalpyStepInput {
            initial: frozen_body,
            ..body
        };
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut frozen_gas,
            &mut frozen_nodes,
            &mut frozen_body,
            frozen_input,
            &inputs,
            &reactions,
            &frozen_ports,
            paths,
            DT,
            shared_solve(),
            &mut || false,
        )
        .unwrap();
        let report = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
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
        assert_eq!(wall_node(&nodes).state.step_count(), 25 + sample);
        assert!(report.port(30).unwrap().solid_to_gas_heat_j > 0.0);
        assert!(report.reaction(30).is_some());
        assert!(
            (wall_source(&nodes)
                .material()
                .property("young_modulus")
                .unwrap()
                .value_si()
                - young(state.temperature_k()))
            .abs()
                < 1e-6
        );
        assert!(
            report.body_energy_residual_j.abs() <= shared_solve().body_exchange.energy_tolerance_j
        );
        assert!(
            report.combined_energy_residual_j.abs()
                <= shared_solve().reactive.transport.energy_tolerance_j
        );
    }
    let pressure_difference = wall_node(&nodes)
        .state
        .pressure()
        .iter()
        .zip(wall_node(&frozen_nodes).state.pressure())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        pressure_difference > 1e-8,
        "the continuing wall must use its new stiffness"
    );
    assert_eq!(
        wall_source(&frozen_nodes).material().query().axes()["T"],
        300.0
    );
    assert_eq!(
        wall_source(&nodes).material().query().axes()["T"],
        state.temperature_k()
    );
}

#[test]
fn g4_wall_work_source_inventory_and_late_cancellation_refuse_atomically() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = material(0.0, 800.0, false);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let (mut parcels, mut nodes, inputs) = wall_fixture(&law, &card);
    let (before_gas, before_nodes) = (parcels.clone(), nodes.clone());
    let mut state = initial;
    let mut body = wall_body(caloric.curve(), initial, 0.0);
    body.internal_heat_j = body.mass_kg * BODY_CP * 20.0;
    let ports = wall_ports(&card, &key, &context, 0.0);
    let (mut accepted_gas, mut accepted_nodes, mut accepted_body) =
        (parcels.clone(), nodes.clone(), state);
    let mut polls = 0;
    let accepted = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut accepted_gas,
        &mut accepted_nodes,
        &mut accepted_body,
        body,
        &inputs,
        &BTreeMap::new(),
        &ports,
        GasTransportPaths::default(),
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
        &BTreeMap::new(),
        &ports,
        GasTransportPaths::default(),
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
    assert_eq!(visited, polls);
    assert_eq!(parcels, before_gas);
    same_owners(&nodes, &before_nodes);
    assert_eq!(state, initial);
    let (foreign_card, foreign_key, _) = material(0.0, 790.0, false);
    for refusal in 0..4 {
        let mut bad_body = body;
        let mut bad_ports = wall_ports(&card, &key, &context, 0.0);
        match refusal {
            0 => {
                bad_ports
                    .get_mut(&30)
                    .unwrap()
                    .wall
                    .as_mut()
                    .unwrap()
                    .maximum_absolute_parameter_work_j = 0.0
            }
            1 => {
                bad_body.mass_kg *= 0.5;
                bad_body.volume_m3 *= 0.5;
            }
            2 => bad_body.volume_m3 *= 2.0,
            _ => bad_ports = wall_ports(&foreign_card, &foreign_key, &context, 0.0),
        }
        assert!(
            advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
                &mut parcels,
                &mut nodes,
                &mut state,
                bad_body,
                &inputs,
                &BTreeMap::new(),
                &bad_ports,
                GasTransportPaths::default(),
                DT,
                shared_solve(),
                &mut || false,
            )
            .is_err(),
            "wall refusal {refusal}"
        );
        assert_eq!(parcels, before_gas);
        same_owners(&nodes, &before_nodes);
        assert_eq!(state, initial);
    }
    for wrong_owner in [10, 20] {
        let (mut gas, mut owners, inputs) = fixture(&law);
        let (old_gas, old_owners) = (gas.clone(), owners.clone());
        let mut body_state = initial;
        let mut ports = shared_ports(&card, &key, &context, [10, 20]);
        ports.get_mut(&wrong_owner).unwrap().wall = Some(SharedFiniteGasAcousticWallSolve {
            maximum_absolute_parameter_work_j: 1e-4,
        });
        assert!(
            advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
                &mut gas,
                &mut owners,
                &mut body_state,
                body,
                &inputs,
                &BTreeMap::new(),
                &ports,
                GasTransportPaths::default(),
                DT,
                shared_solve(),
                &mut || false,
            )
            .is_err()
        );
        assert_eq!(gas, old_gas);
        equal_nodes(&owners, &old_owners);
        assert_eq!(body_state, initial);
    }
    // Initial points are admissible. Only the private heated endpoint leaves
    // the native wall's source range or changes its fixed-density inventory.
    for (density_slope, elastic_upper) in [(0.0, 315.0), (-0.1, 800.0)] {
        let (card, key, caloric) = material(density_slope, elastic_upper, false);
        let initial = caloric
            .curve()
            .state_at_specific_enthalpy(BODY_CP * 100.0)
            .unwrap();
        let (mut gas, mut owners, inputs) = wall_fixture(&law, &card);
        let (old_gas, old_owners) = (gas.clone(), owners.clone());
        let mut body_state = initial;
        let body = wall_body(caloric.curve(), initial, body.internal_heat_j);
        let ports = wall_ports(&card, &key, &context, 0.0);
        assert!(
            advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
                &mut gas,
                &mut owners,
                &mut body_state,
                body,
                &inputs,
                &BTreeMap::new(),
                &ports,
                GasTransportPaths::default(),
                DT,
                shared_solve(),
                &mut || false,
            )
            .is_err()
        );
        assert_eq!(gas, old_gas);
        same_owners(&owners, &old_owners);
        assert_eq!(body_state, initial);
    }
    // Equal density does not license liquid material to retain an elastic
    // solid wall. Check both an accepted mixed phase and a private heated
    // endpoint crossing the latent plateau on the same original card.
    let phase_curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 200.0, 0.0),
            (50_000.0, 300.0, 0.0),
            (55_000.0, 310.0, 0.0),
            (65_000.0, 310.0, 1.0),
            (310_000.0, 800.0, 1.0),
        ]
        .map(
            |(specific_enthalpy_j_kg, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: WALL_DENSITY,
            },
        )
        .to_vec(),
    )
    .unwrap();
    for enthalpy in [50_000.0, 60_000.0] {
        let phase_initial = phase_curve.state_at_specific_enthalpy(enthalpy).unwrap();
        let mut phase_state = phase_initial;
        let phase_body = wall_body(&phase_curve, phase_initial, body.internal_heat_j);
        let refusal = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut phase_state,
            phase_body,
            &inputs,
            &BTreeMap::new(),
            &ports,
            GasTransportPaths::default(),
            DT,
            shared_solve(),
            &mut || false,
        );
        assert!(matches!(
            refusal,
            Err(ThermalReactiveFiniteGasAcousticTransportError::Invalid(what))
                if what.contains("all-solid")
        ));
        assert_eq!(parcels, before_gas);
        same_owners(&nodes, &before_nodes);
        assert_eq!(phase_state, phase_initial);
    }
    let retry = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut state,
        body,
        &inputs,
        &BTreeMap::new(),
        &ports,
        GasTransportPaths::default(),
        DT,
        shared_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, accepted_gas);
    same_owners(&nodes, &accepted_nodes);
    assert_eq!(state, accepted_body);
    assert_eq!(retry.wall_updates, accepted.wall_updates);
    assert_eq!(retry.wall_parameter_work_j, accepted.wall_parameter_work_j);
}
