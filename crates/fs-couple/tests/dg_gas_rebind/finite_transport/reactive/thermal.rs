//! One material store, one reciprocal heat transfer, and continuing acoustic state.
use super::*;
use fs_couple::cavity_phs::piston::thermal::{HeatExchangeSolve, layers::MaterialConductionLayers};
use fs_couple::dg_port::finite_gas::transport::assembly::reactive::thermal::{
    FiniteGasAcousticThermalBinding, ThermalReactiveFiniteGasAcousticTransportError,
    advance_thermal_reacting_finite_gas_acoustic_transport_with_paths,
    finite_gas_thermal_reacting_acoustic_transport_candidate_with_paths,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_qty::{Density, Temperature};

const BODY_MASS: f64 = 1e-5;
const BODY_CP: f64 = 500.0;
const AREA: f64 = 2e-4;
const GAS_DISTANCE: f64 = 0.001;
const BODY_DISTANCE: f64 = 0.0001;
type Bodies = BTreeMap<u64, EquilibriumPhaseState>;

fn body_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}

fn body_material(
    variable_conductivity: bool,
) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut claims = ClaimSet::new();
    for (key, low, high) in [
        (cp.clone(), BODY_CP, BODY_CP),
        (density.clone(), 4000.0, 4000.0),
        (
            conductivity.clone(),
            if variable_conductivity {
                0.2 * 200.0 / 300.0
            } else {
                0.2
            },
            if variable_conductivity {
                0.2 * 800.0 / 300.0
            } else {
                0.2
            },
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, low), (800.0, high)],
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
                    source: "authored finite thermal body; not measured material".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal store".into(),
            phase: "solid".into(),
            process: "constant-capacity regression".into(),
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

fn body_states(curve: &EquilibriumEnthalpyPhaseCurve) -> Bodies {
    [(10, 500.0), (20, 250.0)]
        .map(|(id, t)| {
            (
                id,
                curve
                    .state_at_specific_enthalpy(BODY_CP * (t - 200.0))
                    .unwrap(),
            )
        })
        .into()
}

fn heat_exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-2e-4, 2e-4],
        source_tolerance_j: 1e-13,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 80,
    }
}

fn body_bindings<'a>(
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    states: &Bodies,
) -> BTreeMap<u64, FiniteGasAcousticThermalBinding<'a>> {
    states
        .iter()
        .map(|(&id, &initial)| {
            (
                id,
                FiniteGasAcousticThermalBinding {
                    body: UniformEnthalpyStepInput {
                        curve,
                        initial,
                        mass_kg: BODY_MASS,
                        volume_m3: BODY_MASS / 4000.0,
                        surface_area_m2: AREA,
                        internal_heat_j: if id == 20 { 1e-6 } else { 3e-6 },
                        duration_s: DT,
                    },
                    layers: MaterialConductionLayers {
                        material: card,
                        context,
                        temperature_axis: "T",
                        liquid_fraction_axis: None,
                        conductivity: key,
                        selection: MaterialPropertySelection::SingleClaimOnly,
                        area_m2: AREA,
                        gas_distance_m: GAS_DISTANCE,
                        body_distance_m: BODY_DISTANCE,
                    },
                    exchange: heat_exchange(),
                },
            )
        })
        .collect()
}

// External input work is reconstructed from original physical ports and source
// motion, rather than using the assembly's reported sum as its own oracle.
fn work_and_losses<S>(
    steps: &BTreeMap<u64, FiniteGasAcousticStep<S>>,
    nodes: &MixedNodes,
    inputs: &MixedInputs,
) -> (f64, f64) {
    let (mut work, mut losses) = (0.0, 0.0);
    for (&id, step) in steps {
        match (step, &nodes[&id]) {
            (FiniteGasAcousticStep::Port(step), FiniteGasAcousticNode::Port(load)) => {
                let inlet = step
                    .acoustic
                    .port_outputs
                    .iter()
                    .find(|p| p.id == load.inlet())
                    .unwrap();
                work += thermal(&inputs[&id]).ports[0].effort * inlet.delta_x
                    - load
                        .last_terminal_frames()
                        .iter()
                        .map(|t| t.load_work_j)
                        .sum::<f64>();
                losses += step.acoustic_dissipation_j;
            }
            (FiniteGasAcousticStep::Aperture(step), FiniteGasAcousticNode::Aperture(node)) => {
                let frame = &step.coupled;
                work +=
                    frame.upstream_work_j + frame.body_work_j + frame.aperture.mechanical_work_j
                        - node
                            .load()
                            .last_terminal_frames()
                            .iter()
                            .map(|t| t.load_work_j)
                            .sum::<f64>();
                losses += frame.inferred_total_dissipation_j;
            }
            _ => panic!("fixture changed its physical owner"),
        }
    }
    (work, losses)
}

#[test]
fn g1_finite_material_heat_closes_transport_chemistry_and_acoustic_work() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (card, key, caloric) = body_material(true);
    let context = body_point(300.0);
    let mut bodies = body_states(caloric.curve());
    let original_bodies = bodies.clone();
    let mut links = body_bindings(&card, &key, caloric.curve(), &context, &bodies);
    let (original_gas, original_nodes) = (parcels.clone(), nodes.clone());
    let nozzles = nozzles();
    let paths = GasTransportPaths {
        nozzles: &nozzles,
        ..Default::default()
    };
    let endpoints = inputs
        .iter()
        .map(|(&id, input)| (id, thermal(input).gas))
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
    let (mut no_wall_gas, mut no_wall_nodes) = (parcels.clone(), nodes.clone());
    advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut no_wall_gas,
        &mut no_wall_nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    let report = advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut bodies,
        &inputs,
        &reactions,
        &links,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport.transport, *transported.step());
    let mut transferred = 0.0;
    for id in [10, 20] {
        let wall = report.thermal(id).unwrap();
        let q = wall.solid_to_gas_heat_j;
        let temperature = bodies[&id].temperature_k();
        let gas_temperature = parcels[&id].state().state().temperature;
        let body_k = 0.2 * temperature / 300.0;
        let gas_k = parcels[&id].state().state().thermal_conductivity;
        let conductance = AREA / (GAS_DISTANCE / gas_k + BODY_DISTANCE / body_k);
        assert!(q * (temperature - gas_temperature) > 0.0);
        assert!(
            (q - DT * conductance * (temperature - gas_temperature)).abs()
                <= heat_exchange().source_tolerance_j
        );
        close(wall.conductance_w_k, conductance, conductance);
        close(
            wall.thermal
                .report
                .property_by_key(&key)
                .unwrap()
                .value_si(),
            body_k,
            body_k,
        );
        assert_eq!(wall.thermal.state, bodies[&id]);
        assert_eq!(wall.thermal.external_heat_j, -q);
        assert_eq!(
            wall.thermal.state.material_card_identity(),
            card.content_hash()
        );
        let body_change =
            BODY_MASS * BODY_CP * (temperature - original_bodies[&id].temperature_k());
        assert!((body_change + q - links[&id].body.internal_heat_j).abs() < 1e-11);
        assert!((wall.body_storage_change_j - body_change).abs() < 1e-12);
        let chemical = report.reaction(id).unwrap();
        let work = accepted_work(&report.transport.acoustics[&id]);
        assert!(chemical.reaction.extent_mol > 1e-9);
        assert!(work.abs() > 1e-18);
        assert_eq!(
            chemical.boundary.heat_j,
            thermal(&inputs[&id]).supplied_heat_j + q
        );
        assert_eq!(chemical.boundary.work_by_gas_j, work);
        assert_eq!(
            chemical.endpoint_heat_entropy_j_k,
            chemical.boundary.heat_j / gas_temperature
        );
        let native = law
            .boundary_implicit_candidate(
                &transported.parcels()[&id],
                DT,
                chemical.boundary,
                chemistry(),
                || false,
            )
            .unwrap();
        assert_eq!(native.gas(), &parcels[&id]);
        assert_eq!(native.step(), chemical);
        assert_eq!(nodes[&id].accepted_steps(), 25);
        assert_ne!(parcels[&id], no_wall_gas[&id]);
        transferred += q;
    }
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, g)| energy(g) - energy(&original_gas[id]))
        .sum();
    let acoustic_change: f64 = nodes
        .iter()
        .map(|(id, n)| node_energy(n) - node_energy(&original_nodes[id]))
        .sum();
    let body_change: f64 = bodies
        .iter()
        .map(|(id, state)| {
            BODY_MASS * BODY_CP * (state.temperature_k() - original_bodies[id].temperature_k())
        })
        .sum();
    let external_gas: f64 = inputs
        .values()
        .map(|input| thermal(input).supplied_heat_j)
        .sum();
    let external_body: f64 = links.values().map(|link| link.body.internal_heat_j).sum();
    let (work, losses) = work_and_losses(&report.transport.acoustics, &nodes, &inputs);
    let residual =
        gas_change + body_change + acoustic_change + losses - external_gas - external_body - work;
    assert!(residual.abs() < solve().energy_tolerance_j);
    assert!((report.combined_energy_residual_j - residual).abs() < 2e-12);
    assert!((report.solid_to_gas_heat_j - transferred).abs() < 1e-15);
    assert!((report.transport.supplied_heat_j - external_gas - transferred).abs() < 1e-15);
    assert_eq!(report.external_gas_heat_j, external_gas);
    assert_eq!(report.external_body_heat_j, external_body);
    assert!(report.maximum_species_residual_mol <= reactive_solve().species_tolerance_mol);

    // Continue accepted body and original source histories; do not recreate a
    // cold acoustic owner or repeat heat from the first accepted material state.
    for (&id, binding) in &mut links {
        binding.body.initial = bodies[&id];
    }
    advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut bodies,
        &inputs,
        &reactions,
        &links,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut no_wall_gas,
        &mut no_wall_nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert!(bodies[&10].temperature_k() < links[&10].body.initial.temperature_k());
    assert!(bodies[&20].temperature_k() > links[&20].body.initial.temperature_k());
    assert_eq!(nodes[&20].accepted_steps(), 26);
    let (FiniteGasAcousticNode::Aperture(a), FiniteGasAcousticNode::Aperture(b)) =
        (&nodes[&20], &no_wall_nodes[&20])
    else {
        unreachable!()
    };
    let (LoadedAperture::Paired(a), LoadedAperture::Paired(b)) = (a.aperture(), b.aperture())
    else {
        unreachable!()
    };
    assert_ne!(a.blade_states(), b.blade_states());
}

#[test]
fn g1_finite_reservoir_and_body_match_two_capacity_backward_euler() {
    let law = reaction(50_000.0);
    let (card, key, caloric) = body_material(false);
    let context = body_point(300.0);
    for temperature in [250.0, 500.0] {
        let (mut parcels, mut nodes, inputs) = fixture(&law);
        let mut bodies = BTreeMap::from([(
            90,
            caloric
                .curve()
                .state_at_specific_enthalpy(BODY_CP * (temperature - 200.0))
                .unwrap(),
        )]);
        let links = body_bindings(&card, &key, caloric.curve(), &context, &bodies);
        let gas_capacity = 1.5 * R * parcels[&90].mass_kg() / 0.004;
        let body_capacity = BODY_CP * BODY_MASS;
        // Pure inert source has constant k=0.026 and Cp/R=2.5. With no mass
        // paths or reaction this is an independent closed-form BE heat solve.
        let conductance = AREA / (GAS_DISTANCE / 0.026 + BODY_DISTANCE / 0.2);
        let body_heat = links[&90].body.internal_heat_j;
        let gas_heat = thermal(&inputs[&90]).supplied_heat_j;
        let expected_q = DT
            * conductance
            * (temperature + body_heat / body_capacity - 300.0 - gas_heat / gas_capacity)
            / (1.0 + DT * conductance * (body_capacity.recip() + gas_capacity.recip()));
        let report = advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut bodies,
            &inputs,
            &BTreeMap::new(),
            &links,
            GasTransportPaths::default(),
            DT,
            reactive_solve(),
            &mut || false,
        )
        .unwrap();
        let actual = report.thermal(90).unwrap();
        assert!((actual.solid_to_gas_heat_j - expected_q).abs() < 2e-13);
        assert!(
            (bodies[&90].temperature_k() - temperature - (body_heat - expected_q) / body_capacity)
                .abs()
                < 1e-9
        );
        assert!(
            (parcels[&90].state().state().temperature
                - 300.0
                - (gas_heat + expected_q) / gas_capacity)
                .abs()
                < 1e-9
        );
        assert!(actual.solid_to_gas_heat_j * (temperature - 300.0) > 0.0);
        assert!(report.reaction(90).is_none());
    }
}

#[test]
fn g3_empty_thermal_bindings_preserve_original_reacting_owners() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (mut original_gas, mut original_nodes) = (parcels.clone(), nodes.clone());
    let nozzles = nozzles();
    let paths = GasTransportPaths {
        nozzles: &nozzles,
        ..Default::default()
    };
    let original = advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut original_gas,
        &mut original_nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    let report = advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut BTreeMap::new(),
        &inputs,
        &reactions,
        &BTreeMap::new(),
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, original_gas);
    equal_nodes(&nodes, &original_nodes);
    assert_eq!(report.transport.transport, original.transport.transport);
    assert_eq!(
        report.maximum_species_residual_mol,
        original.maximum_species_residual_mol
    );
    assert_eq!(report.body_storage_change_j, 0.0);
    assert_eq!(report.external_body_heat_j, 0.0);
    assert_eq!(report.solid_to_gas_heat_j, 0.0);
    assert_eq!(
        report.transport.combined_energy_residual_j,
        original.transport.combined_energy_residual_j
    );
    // The enclosing report uses a separate compensated summation. Its roundoff
    // need not equal the original report's grouping when adding zero bodies.
    assert!(report.combined_energy_residual_j.abs() < solve().energy_tolerance_j);
    for id in [10, 20, 90] {
        assert!(report.thermal(id).is_none());
        assert_eq!(report.reaction(id), original.reaction(id));
    }
    for (&id, step) in &report.transport.acoustics {
        match (step, &original.transport.acoustics[&id]) {
            (FiniteGasAcousticStep::Port(a), FiniteGasAcousticStep::Port(b)) => {
                assert_eq!(a.acoustic, b.acoustic);
                assert_eq!(a.refresh, b.refresh);
                assert_eq!(a.background.native, b.background);
            }
            (FiniteGasAcousticStep::Aperture(a), FiniteGasAcousticStep::Aperture(b)) => {
                assert_eq!(a.coupled, b.coupled);
                assert_eq!(a.gas.acoustic, b.gas.acoustic);
                assert_eq!(a.gas.refresh, b.gas.refresh);
                assert_eq!(a.gas.background.native, b.gas.background);
            }
            _ => unreachable!(),
        }
    }
    assert_eq!(
        report.transport.reservoirs[&90].native,
        original.transport.reservoirs[&90]
    );
}

#[test]
fn g4_finite_thermal_refusal_and_final_cancellation_preserve_all_three_owners() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (card, key, caloric) = body_material(true);
    let context = body_point(300.0);
    let mut bodies = body_states(caloric.curve());
    let links = body_bindings(&card, &key, caloric.curve(), &context, &bodies);
    let (before_gas, before_nodes, before_bodies) =
        (parcels.clone(), nodes.clone(), bodies.clone());
    let nozzles = nozzles();
    let paths = GasTransportPaths {
        nozzles: &nozzles,
        ..Default::default()
    };
    let preview = finite_gas_thermal_reacting_acoustic_transport_candidate_with_paths(
        &parcels,
        &nodes,
        &bodies,
        &inputs,
        &reactions,
        &links,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, before_gas);
    assert_eq!(bodies, before_bodies);
    equal_nodes(&nodes, &before_nodes);
    let (mut accepted_gas, mut accepted_nodes, mut accepted_bodies) =
        (parcels.clone(), nodes.clone(), bodies.clone());
    let mut polls = 0;
    let accepted = advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut accepted_gas,
        &mut accepted_nodes,
        &mut accepted_bodies,
        &inputs,
        &reactions,
        &links,
        paths,
        DT,
        reactive_solve(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(preview.parcels(), &accepted_gas);
    assert_eq!(preview.body_states(), &accepted_bodies);
    equal_nodes(preview.nodes(), &accepted_nodes);
    let mut visited = 0;
    let cancelled = advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut bodies,
        &inputs,
        &reactions,
        &links,
        paths,
        DT,
        reactive_solve(),
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
    assert_eq!(bodies, before_bodies);
    equal_nodes(&nodes, &before_nodes);
    // Initial sources admit; the second owner's heated endpoint leaves the
    // original chart only after the first physical owner has been proposed.
    let mut excessive = body_bindings(&card, &key, caloric.curve(), &context, &bodies);
    excessive.get_mut(&20).unwrap().body.internal_heat_j = 100.0;
    assert!(
        advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut bodies,
            &inputs,
            &reactions,
            &excessive,
            paths,
            DT,
            reactive_solve(),
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(parcels, before_gas);
    assert_eq!(bodies, before_bodies);
    equal_nodes(&nodes, &before_nodes);
    let retry = advance_thermal_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &mut bodies,
        &inputs,
        &reactions,
        &links,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, accepted_gas);
    assert_eq!(bodies, accepted_bodies);
    equal_nodes(&nodes, &accepted_nodes);
    assert_eq!(retry.transport.transport, accepted.transport.transport);
    assert_eq!(
        retry.combined_energy_residual_j,
        accepted.combined_energy_residual_j
    );
    for id in [10, 20] {
        assert_eq!(
            retry.thermal(id).unwrap().thermal,
            accepted.thermal(id).unwrap().thermal
        );
        assert_eq!(retry.reaction(id), accepted.reaction(id));
    }
}
