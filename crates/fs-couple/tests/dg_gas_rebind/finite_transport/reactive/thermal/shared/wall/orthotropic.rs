//! Principal material axes and proportional viscosity on the finite-body path.
//! These complete source curves are manufactured regression data.
use super::*;
use fs_couple::dg_port::finite_gas::transport::assembly::reactive::thermal::shared::{
    SharedThermalReactiveFiniteGasAcousticTransportCandidate,
    SharedThermalReactiveFiniteGasAcousticTransportStep,
};
use fs_phs::linear_exact::DenseExponentialBudget;

const MODULI: [f64; 2] = [1e8, 4e7];
const NU: f64 = 0.25;
const TIME: &str = "normal_plane_viscosity_time";
const VISCOSITY_TIME: f64 = 2e-6;
const VISCOSITY_SLOPE: f64 = 0.003;
const BAND: [f64; 2] = [20.0, 100_000.0];

#[derive(Clone, Copy, Debug)]
enum SourceCase {
    Supported,
    ShortDomain,
    IndefiniteCompliance,
    ChangingDensity,
}

fn source(
    case: SourceCase,
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
        (
            "density",
            Density::DIMS,
            temperatures.map(|t| {
                WALL_DENSITY
                    + if matches!(case, SourceCase::ChangingDensity) {
                        0.1 * (t - 300.0)
                    } else {
                        0.0
                    }
            }),
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, [0.2; 4]),
        (
            "young_modulus_1",
            Pressure::DIMS,
            temperatures.map(|t| MODULI[0] * (1.0 + YOUNG_SLOPE * (t - 300.0))),
        ),
        (
            "young_modulus_2",
            Pressure::DIMS,
            temperatures.map(|t| MODULI[1] * (1.0 + YOUNG_SLOPE * (t - 300.0))),
        ),
        (
            "poisson_ratio_12",
            Dims::NONE,
            temperatures.map(|t| {
                if matches!(case, SourceCase::IndefiniteCompliance) && t >= 320.0 {
                    2.0
                } else {
                    NU
                }
            }),
        ),
        (
            TIME,
            Dims([0, 0, 1, 0, 0, 0]),
            temperatures.map(|t| VISCOSITY_TIME * (1.0 + VISCOSITY_SLOPE * (t - 300.0))),
        ),
    ] {
        let upper = if matches!(case, SourceCase::ShortDomain) && name == "young_modulus_2" {
            315.0
        } else {
            800.0
        };
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
                        upper,
                    )
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured directional finite-body membrane".into(),
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
            chemistry: "synthetic directional wall".into(),
            phase: "solid".into(),
            process: "principal plane and proportional viscosity".into(),
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
    let material = resolve_material_state_point(
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
        material,
    )
}

fn fixture(
    law: &GasReaction,
    material: ResolvedMaterialStatePoint,
    axes: [usize; 2],
    axial: AxialConstraint,
    viscous: bool,
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
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap();
    let wall = LocalMaterialWall::new_orthotropic(
        material,
        THICKNESS,
        axial,
        axes,
        viscous.then_some(TIME),
        BAND,
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

fn frozen(gas: &Parcels, nodes: &MixedNodes, inputs: &MixedInputs) -> MixedNodes {
    finite_gas_reacting_acoustic_transport_candidate_with_paths(
        gas,
        nodes,
        inputs,
        &BTreeMap::new(),
        GasTransportPaths::default(),
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap()
    .nodes()
    .clone()
}

// Independent dimensional normal-compliance reduction, with nu12/E1 retained
// when the hoop and axial directions are reversed.
fn effective_modulus(axes: [usize; 2], axial: AxialConstraint) -> f64 {
    let shh = 1.0 / MODULI[axes[0]];
    let saa = 1.0 / MODULI[axes[1]];
    let sha = -NU / MODULI[0];
    match axial {
        AxialConstraint::StressFree => shh.recip(),
        AxialConstraint::StrainFixed => (shh - sha * sha / saa).recip(),
    }
}

#[test]
fn g1_directional_temperature_work_keeps_axes_restraint_and_source_claims() {
    let law = reaction(50_000.0);
    let (card, key, caloric, source) = source(SourceCase::Supported);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let mut body = wall_body(caloric.curve(), initial, 0.0);
    let capacity = body.mass_kg * BODY_CP;
    body.internal_heat_j = capacity * 20.0;
    let ports = wall_ports(&card, &key, &context, 0.0);
    for axes in [[0, 1], [1, 0]] {
        for axial in [AxialConstraint::StressFree, AxialConstraint::StrainFixed] {
            let (gas, nodes, inputs) = fixture(&law, source.clone(), axes, axial, false);
            let unchanged = frozen(&gas, &nodes, &inputs);
            let old_node = wall_node(&unchanged);
            let old_motion = old_node
                .binding
                .duct()
                .compliant_wall_state(&old_node.state)
                .unwrap();
            let spring = spring_energy(&old_motion.displacement_m, effective_modulus(axes, axial));
            assert!(spring > 1e-12);
            let expected_delta = body.internal_heat_j / (capacity + YOUNG_SLOPE * spring);
            let candidate = preview(&gas, &nodes, body, &inputs, &ports);
            let report = candidate.step();
            let temperature = candidate.body_state().temperature_k();
            let expected_work = YOUNG_SLOPE * spring * (temperature - 300.0);
            let node = wall_node(candidate.nodes());
            let motion = node
                .binding
                .duct()
                .compliant_wall_state(&node.state)
                .unwrap();
            assert!((temperature - 300.0 - expected_delta).abs() < 1e-9);
            assert!(report.wall_parameter_work_j < 0.0);
            assert!((report.wall_parameter_work_j - expected_work).abs() < 1e-13);
            assert!(
                (motion.stored_energy_j - old_motion.stored_energy_j - expected_work).abs() < 1e-12
            );
            assert!(
                (capacity * (temperature - 300.0) + expected_work - body.internal_heat_j).abs()
                    < 2e-11
            );
            assert_eq!(report.wall_heat_j, 0.0);
            assert_eq!(report.solid_to_gas_heat_j, 0.0);
            assert_eq!(
                report.final_acoustic_storage_change_j,
                node.state.stored_energy_j() - wall_node(&nodes).state.stored_energy_j()
            );
            assert_eq!(
                wall_source(candidate.nodes()).orthotropic_axes(),
                Some(axes)
            );
            assert_eq!(wall_source(candidate.nodes()).axial_constraint(), axial);
            assert_eq!(
                wall_source(candidate.nodes()).proportional_viscosity_time_property(),
                None
            );
            let material = report.wall_updates[&30].uniform_material().unwrap();
            assert_eq!(material.query(), body_point(temperature));
            for name in [
                "density",
                "young_modulus_1",
                "young_modulus_2",
                "poisson_ratio_12",
            ] {
                assert_eq!(
                    material.property(name).unwrap().answer().receipt.selected,
                    source.property(name).unwrap().answer().receipt.selected
                );
            }
            assert!(material.property("young_modulus").is_none());
            for i in 0..motion.displacement_m.len() {
                assert!((motion.displacement_m[i] - old_motion.displacement_m[i]).abs() < 1e-15);
                assert!((motion.velocity_m_s[i] - old_motion.velocity_m_s[i]).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn g3_proportional_directional_wall_heat_changes_the_same_continuing_operator() {
    let law = reaction(50_000.0);
    let (card, key, caloric, source) = source(SourceCase::Supported);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let body = wall_body(caloric.curve(), initial, 0.0);
    let (gas, nodes, inputs) = fixture(&law, source, [1, 0], AxialConstraint::StrainFixed, true);
    let mut ports = wall_ports(&card, &key, &context, 0.0);
    let budget = DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    };
    ports
        .get_mut(&30)
        .unwrap()
        .wall
        .as_mut()
        .unwrap()
        .dissipation_heat = Some(budget);
    let original = wall_node(&nodes);
    let heat = original
        .binding
        .duct()
        .prepare_compliant_wall_heat(&original.propagator, budget, &mut || false)
        .unwrap()
        .integrated_heat_j(&original.state, &held_ports(), &mut || false)
        .unwrap();
    assert!(heat > 1e-11);
    let unchanged = frozen(&gas, &nodes, &inputs);
    let before = wall_node(&unchanged)
        .binding
        .duct()
        .compliant_wall_state(&wall_node(&unchanged).state)
        .unwrap();
    let candidate = preview(&gas, &nodes, body, &inputs, &ports);
    let report = candidate.step();
    let temperature = candidate.body_state().temperature_k();
    let after = wall_node(candidate.nodes())
        .binding
        .duct()
        .compliant_wall_state(&wall_node(candidate.nodes()).state)
        .unwrap();
    let expected_work = spring_energy(
        &before.displacement_m,
        effective_modulus([1, 0], AxialConstraint::StrainFixed),
    ) * YOUNG_SLOPE
        * (temperature - 300.0);
    assert!(temperature > 300.0);
    assert_eq!(report.wall_heat_j, heat);
    assert_eq!(report.wall_updates[&30].dissipated_heat_j, heat);
    assert!((report.wall_parameter_work_j - expected_work).abs() < 1e-13);
    assert!((report.body_storage_change_j + expected_work - heat).abs() < 2e-11);
    assert_eq!(
        report.final_exported_acoustic_dissipation_j,
        report.transport.exported_acoustic_dissipation_j - heat
    );
    assert_eq!(
        wall_source(candidate.nodes()).orthotropic_axes(),
        Some([1, 0])
    );
    assert_eq!(
        wall_source(candidate.nodes()).proportional_viscosity_time_property(),
        Some(TIME)
    );
    let material = report.wall_updates[&30].uniform_material().unwrap();
    assert!(
        (material.property(TIME).unwrap().value_si()
            - VISCOSITY_TIME * (1.0 + VISCOSITY_SLOPE * (temperature - 300.0)))
            .abs()
            < 1e-16
    );
    let expected_power_ratio = (1.0 + YOUNG_SLOPE * (temperature - 300.0))
        * (1.0 + VISCOSITY_SLOPE * (temperature - 300.0));
    assert!(before.dissipated_power_w > 0.0);
    assert!(
        (after.dissipated_power_w / before.dissipated_power_w - expected_power_ratio).abs() < 2e-11
    );
    let mut continuing = wall_node(candidate.nodes()).clone();
    let mut frozen_next = wall_node(&unchanged).clone();
    for node in [&mut continuing, &mut frozen_next] {
        node.propagator
            .step(&mut node.state, &held_ports())
            .unwrap();
        node.propagator.sync_nodal(&mut node.state).unwrap();
    }
    assert_eq!(
        continuing.state.step_count(),
        original.state.step_count() + 2
    );
    let change = continuing
        .state
        .pressure()
        .iter()
        .zip(frozen_next.state.pressure())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        change > 1e-12,
        "the retained native operator must use the refreshed directional law"
    );
}

#[test]
fn g4_directional_source_compliance_density_and_late_cancel_preserve_owners() {
    let law = reaction(50_000.0);
    for case in [
        SourceCase::ShortDomain,
        SourceCase::IndefiniteCompliance,
        SourceCase::ChangingDensity,
    ] {
        let (card, key, caloric, material) = source(case);
        let context = body_point(300.0);
        let initial = caloric
            .curve()
            .state_at_specific_enthalpy(BODY_CP * 100.0)
            .unwrap();
        let mut state = initial;
        let mut body = wall_body(caloric.curve(), initial, 0.0);
        body.internal_heat_j = body.mass_kg * BODY_CP * 20.0;
        let (mut gas, mut nodes, inputs) =
            fixture(&law, material, [0, 1], AxialConstraint::StrainFixed, false);
        let (before_gas, before_nodes) = (gas.clone(), nodes.clone());
        let ports = wall_ports(&card, &key, &context, 0.0);
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
            .is_err(),
            "{case:?}"
        );
        assert_eq!(gas, before_gas);
        same_owners(&nodes, &before_nodes);
        assert_eq!(state, initial);
    }
    let (card, key, caloric, material) = source(SourceCase::Supported);
    let context = body_point(300.0);
    let initial = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let mut state = initial;
    let mut body = wall_body(caloric.curve(), initial, 0.0);
    body.internal_heat_j = body.mass_kg * BODY_CP * 20.0;
    let ports = wall_ports(&card, &key, &context, 0.0);
    let (mut gas, mut nodes, inputs) =
        fixture(&law, material, [1, 0], AxialConstraint::StrainFixed, true);
    let (before_gas, before_nodes) = (gas.clone(), nodes.clone());
    let (mut accepted_gas, mut accepted_nodes, mut accepted_body) =
        (gas.clone(), nodes.clone(), initial);
    let mut polls = 0;
    let accepted = advance(
        &mut accepted_gas,
        &mut accepted_nodes,
        &mut accepted_body,
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
    let failure = advance(
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
        failure,
        Err(ThermalReactiveFiniteGasAcousticTransportError::Cancelled)
    ));
    assert_eq!(visited, polls);
    assert_eq!(gas, before_gas);
    same_owners(&nodes, &before_nodes);
    assert_eq!(state, initial);
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
    assert_eq!(state, accepted_body);
    assert_eq!(retry.wall_updates, accepted.wall_updates);
    assert_eq!(retry.wall_parameter_work_j, accepted.wall_parameter_work_j);
}
