//! Original solid loss heats its finite body while native wall memory continues.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::WallRelaxation;
use fs_couple::dg_port::finite_gas::transport::assembly::reactive::thermal::shared::{
    SharedThermalReactiveFiniteGasAcousticTransportCandidate,
    SharedThermalReactiveFiniteGasAcousticTransportStep,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_qty::DynViscosity;

const ARM: f64 = 4e7;
const ARM_SLOPE: f64 = 0.002;
const TAU: f64 = 4e-5;
const TAU_SLOPE: f64 = 0.005;
const BAND: [f64; 2] = [20.0, 100_000.0];

#[derive(Clone, Copy)]
enum Loss {
    Elastic,
    Kelvin,
    Maxwell,
    Authored,
}

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn source(
    vary_moduli: bool,
    extinguish_arm: bool,
) -> (
    MaterialCard,
    PropertyKey,
    ResolvedHeatCapacityCurve,
    ResolvedMaterialStatePoint,
) {
    let temperatures = [200.0, 300.0, 320.0, 800.0];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, values) in [
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [BODY_CP; 4],
        ),
        ("density", Density::DIMS, [WALL_DENSITY; 4]),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, [0.2; 4]),
        (
            "young_modulus",
            Pressure::DIMS,
            temperatures.map(|t| if vary_moduli { young(t) } else { YOUNG_300 }),
        ),
        ("poisson_ratio", Dims::NONE, [0.3; 4]),
        (
            "arm_modulus",
            Pressure::DIMS,
            temperatures.map(|t| {
                if extinguish_arm && t >= 320.0 {
                    0.0
                } else if vary_moduli {
                    ARM * (1.0 + ARM_SLOPE * (t - 300.0))
                } else {
                    ARM
                }
            }),
        ),
        (
            "arm_time",
            Dims([0, 0, 1, 0, 0, 0]),
            temperatures.map(|t| TAU * (1.0 + TAU_SLOPE * (t - 300.0))),
        ),
        (
            "young_viscosity",
            DynViscosity::DIMS,
            temperatures.map(|t| 100.0 * (1.0 + 0.003 * (t - 300.0))),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: temperatures.into_iter().zip(values).collect(),
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
                        800.0,
                    )
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic finite-body wall memory; not measured validation".into(),
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
            chemistry: "synthetic thermal viscoelastic inventory".into(),
            phase: "solid".into(),
            process: "fixed density and complete original source identities".into(),
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
    let point = resolve_material_state_point(
        &card,
        &body_point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        card,
        PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        caloric,
        point,
    )
}

fn memory_fixture(
    law: &GasReaction,
    source: &ResolvedMaterialStatePoint,
    loss: Loss,
) -> (Parcels, MixedNodes, MixedInputs) {
    let loss = match loss {
        Loss::Elastic => WallMaterialLoss::Elastic,
        Loss::Kelvin => WallMaterialLoss::KelvinVoigt,
        Loss::Maxwell | Loss::Authored => {
            let pinned = WallRelaxation::from_state_point(
                source,
                &[("arm_modulus".into(), "arm_time".into())],
                BAND,
                1.0,
                0.8,
            )
            .unwrap();
            WallMaterialLoss::Maxwell(if matches!(loss, Loss::Authored) {
                WallRelaxation::new(
                    pinned.material().clone(),
                    "same coefficients without original source-arm identity".into(),
                    BAND,
                    1.0,
                    0.8,
                )
                .unwrap()
            } else {
                pinned
            })
        }
    };
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
        source.clone(),
        THICKNESS,
        AxialConstraint::StressFree,
        loss,
        BAND,
    )
    .unwrap()
    .bind_dg(duct)
    .unwrap();
    let propagator = wall.duct().propagator(DT, 1e-11).unwrap();
    let mut state = wall
        .duct()
        .project_initial_state(|x| 1000.0 + 3000.0 * (1.0 - 2.0 * x / LENGTH), |_| 1e-9);
    for _ in 0..4 {
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
fn g1_original_maxwell_heat_warms_one_body_and_preserves_physical_memory() {
    let law = reaction(50_000.0);
    let (card, key, caloric, source) = source(true, false);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let body = wall_body(caloric.curve(), initial, 0.0);
    let capacity = body.mass_kg * BODY_CP;
    let context = body_point(300.0);
    let (gas, nodes, inputs) = memory_fixture(&law, &source, Loss::Maxwell);
    let mut ports = wall_ports(&card, &key, &context, 0.0);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(heat_budget());
    let original = wall_node(&nodes);
    let heat = original
        .binding
        .duct()
        .prepare_compliant_wall_heat(&original.propagator, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&original.state, &held_ports(), &mut || false)
        .unwrap();
    assert!(heat > 1e-11);
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
    assert_eq!(before.relaxation_pressure_pa.len(), 1);
    let arm_stiffness = ARM * THICKNESS / RADIUS.powi(2);
    let extension: Vec<_> = before.relaxation_pressure_pa[0]
        .iter()
        .map(|p| p / arm_stiffness)
        .collect();
    let arm_energy = spring_energy(&extension, ARM);
    let work_per_kelvin =
        spring_energy(&before.displacement_m, YOUNG_300) * YOUNG_SLOPE + arm_energy * ARM_SLOPE;
    assert!((arm_energy - before.relaxation_stored_energy_j).abs() < 1e-13);
    let expected_delta = heat / (capacity + work_per_kelvin);
    let candidate = preview(&gas, &nodes, body, &inputs, &ports);
    let step = candidate.step();
    let new_temperature = candidate.body_state().temperature_k();
    let after_node = wall_node(candidate.nodes());
    let after = after_node
        .binding
        .duct()
        .compliant_wall_state(&after_node.state)
        .unwrap();
    let expected_work = work_per_kelvin * (new_temperature - 300.0);
    assert!(new_temperature > 300.0);
    assert!((new_temperature - 300.0 - expected_delta).abs() < 2e-10);
    assert_eq!(step.wall_heat_j, heat);
    assert_eq!(step.wall_updates[&30].dissipated_heat_j, heat);
    assert_eq!(step.external_body_heat_j, 0.0);
    assert_eq!(step.external_gas_heat_j, 0.0);
    assert_eq!(step.solid_to_gas_heat_j, 0.0);
    assert!((step.wall_parameter_work_j - expected_work).abs() < 1e-13);
    let body_storage = capacity * (new_temperature - 300.0);
    assert!((body_storage + expected_work - heat).abs() < 3e-12);
    assert_eq!(
        step.final_exported_acoustic_dissipation_j,
        step.transport.exported_acoustic_dissipation_j - heat
    );
    let acoustic_storage = after_node.state.stored_energy_j() - original.state.stored_energy_j();
    let combined = energy(&candidate.parcels()[&30]) - energy(&gas[&30])
        + acoustic_storage
        + body_storage
        + step.final_exported_acoustic_dissipation_j
        - step.transport.port_work_j;
    assert!(combined.abs() < shared_solve().reactive.transport.energy_tolerance_j);
    let new_arm = ARM * (1.0 + ARM_SLOPE * (new_temperature - 300.0));
    for i in 0..before.displacement_m.len() {
        assert!((after.displacement_m[i] - before.displacement_m[i]).abs() < 1e-15);
        assert!((after.velocity_m_s[i] - before.velocity_m_s[i]).abs() < 1e-12);
        assert!(
            (after.relaxation_pressure_pa[0][i] / (new_arm * THICKNESS / RADIUS.powi(2))
                - extension[i])
                .abs()
                < 1e-15
        );
    }
    assert_eq!(
        after_node.state.step_count(),
        original.state.step_count() + 1
    );
    assert_eq!(
        after_node.state.cumulative_dissipated_j(),
        wall_node(frozen.nodes()).state.cumulative_dissipated_j()
    );
    let material = &step.wall_updates[&30].material;
    assert_eq!(material.query(), body_point(new_temperature));
    for name in ["young_modulus", "arm_modulus", "arm_time"] {
        assert_eq!(
            material.property(name).unwrap().answer().receipt.selected,
            source.property(name).unwrap().answer().receipt.selected
        );
    }
    assert!(
        (material.property("arm_time").unwrap().value_si()
            - TAU * (1.0 + TAU_SLOPE * (new_temperature - 300.0)))
            .abs()
            < 1e-15
    );
    let continued = preview(
        candidate.parcels(),
        candidate.nodes(),
        UniformEnthalpyStepInput {
            initial: *candidate.body_state(),
            ..body
        },
        &inputs,
        &ports,
    );
    assert_eq!(
        wall_node(continued.nodes()).state.step_count(),
        after_node.state.step_count() + 1
    );
    assert!(continued.step().wall_heat_j > 0.0);
    assert!(continued.body_state().temperature_k() > new_temperature);
    assert_eq!(
        wall_source(continued.nodes()).material().query().axes()["T"],
        continued.body_state().temperature_k()
    );
}

#[test]
fn g3_kelvin_heat_elastic_zero_and_tau_only_refresh_keep_separate_accounts() {
    let law = reaction(50_000.0);
    let (card, key, caloric, source) = source(false, false);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let context = body_point(300.0);
    for loss in [Loss::Kelvin, Loss::Elastic] {
        let (gas, nodes, inputs) = memory_fixture(&law, &source, loss);
        let mut ports = wall_ports(&card, &key, &context, 0.0);
        let body = wall_body(caloric.curve(), initial, 0.0);
        let exported = preview(&gas, &nodes, body, &inputs, &ports);
        ports
            .get_mut(&30)
            .unwrap()
            .wall
            .as_mut()
            .unwrap()
            .dissipation_heat = Some(heat_budget());
        let deposited = preview(&gas, &nodes, body, &inputs, &ports);
        assert_eq!(exported.step().wall_heat_j, 0.0);
        assert_eq!(
            exported.step().final_exported_acoustic_dissipation_j,
            exported.step().transport.exported_acoustic_dissipation_j
        );
        assert_eq!(*exported.body_state(), initial);
        assert_eq!(deposited.step().wall_parameter_work_j, 0.0);
        if matches!(loss, Loss::Kelvin) {
            assert!(deposited.step().wall_heat_j > 0.0);
            assert!(deposited.body_state().temperature_k() > 300.0);
            assert!(
                wall_source(deposited.nodes())
                    .material()
                    .property("young_viscosity")
                    .unwrap()
                    .value_si()
                    > 100.0
            );
            assert!(
                (deposited.step().body_storage_change_j - deposited.step().wall_heat_j).abs()
                    < 1e-11
            );
        } else {
            assert_eq!(deposited.step().wall_heat_j, 0.0);
            assert_eq!(*deposited.body_state(), initial);
            same_owners(deposited.nodes(), exported.nodes());
            assert_eq!(deposited.parcels(), exported.parcels());
        }
    }
    // Changing tau changes future damping, without instantaneous spring work.
    let (gas, nodes, inputs) = memory_fixture(&law, &source, Loss::Maxwell);
    let ports = wall_ports(&card, &key, &context, 0.0);
    let mut body = wall_body(caloric.curve(), initial, 0.0);
    body.internal_heat_j = body.mass_kg * BODY_CP * 20.0;
    let candidate = preview(&gas, &nodes, body, &inputs, &ports);
    assert!((candidate.body_state().temperature_k() - 320.0).abs() < 1e-11);
    assert_eq!(candidate.step().wall_parameter_work_j, 0.0);
    assert_eq!(candidate.step().wall_heat_j, 0.0);
    assert!(
        (wall_source(candidate.nodes())
            .material()
            .property("arm_time")
            .unwrap()
            .value_si()
            - 1.1 * TAU)
            .abs()
            < 1e-15
    );
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
    let power = |nodes: &MixedNodes| {
        let node = wall_node(nodes);
        node.binding
            .duct()
            .compliant_wall_state(&node.state)
            .unwrap()
            .dissipated_power_w
    };
    assert!(power(frozen.nodes()) > 0.0);
    assert!((power(candidate.nodes()) / power(frozen.nodes()) - 1.0 / 1.1).abs() < 1e-11);
}

#[test]
fn g4_original_memory_source_heat_budget_and_late_cancel_are_atomic() {
    let law = reaction(50_000.0);
    let (card, key, caloric, source) = source(true, false);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let context = body_point(300.0);
    let body = wall_body(caloric.curve(), initial, 0.0);
    let (mut gas, mut nodes, inputs) = memory_fixture(&law, &source, Loss::Maxwell);
    let before_gas = gas.clone();
    let before_nodes = nodes.clone();
    let mut state = initial;
    let mut ports = wall_ports(&card, &key, &context, 0.0);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(heat_budget());
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
    let refused = advance(
        &mut gas,
        &mut nodes,
        &mut state,
        body,
        &inputs,
        &ports,
        &mut || {
            visited += 1;
            visited == polls
        },
    );
    assert!(matches!(
        refused,
        Err(ThermalReactiveFiniteGasAcousticTransportError::Cancelled)
    ));
    assert_eq!(visited, polls);
    assert_eq!(gas, before_gas);
    same_owners(&nodes, &before_nodes);
    assert_eq!(state, initial);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(DenseExponentialBudget {
        max_dimension: 1,
        ..heat_budget()
    });
    assert!(
        advance(
            &mut gas,
            &mut nodes,
            &mut state,
            body,
            &inputs,
            &ports,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(gas, before_gas);
    same_owners(&nodes, &before_nodes);
    assert_eq!(state, initial);
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(heat_budget());
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
    same_owners(&nodes, &accepted_nodes);
    assert_eq!(state, accepted_state);
    assert_eq!(retry.wall_heat_j, accepted.wall_heat_j);
    assert_eq!(retry.wall_updates, accepted.wall_updates);

    for (authored, extinguish) in [(true, false), (false, true)] {
        let (card, key, caloric, source) = self::source(true, extinguish);
        let initial = caloric
            .curve()
            .state_at_specific_enthalpy(BODY_CP * 100.0)
            .unwrap();
        let (mut gas, mut nodes, inputs) = memory_fixture(
            &law,
            &source,
            if authored {
                Loss::Authored
            } else {
                Loss::Maxwell
            },
        );
        let (old_gas, old_nodes) = (gas.clone(), nodes.clone());
        let mut state = initial;
        let mut body = wall_body(caloric.curve(), initial, 0.0);
        body.internal_heat_j = body.mass_kg * BODY_CP * 25.0;
        let mut ports = wall_ports(&card, &key, &context, 0.0);
        ports
            .get_mut(&30)
            .unwrap()
            .wall
            .as_mut()
            .unwrap()
            .dissipation_heat = Some(heat_budget());
        assert!(
            advance(
                &mut gas,
                &mut nodes,
                &mut state,
                body,
                &inputs,
                &ports,
                &mut || false
            )
            .is_err()
        );
        assert_eq!(gas, old_gas);
        same_owners(&nodes, &old_nodes);
        assert_eq!(state, initial);
    }
}
