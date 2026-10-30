//! One shared material inventory couples several original finite gas owners.
mod molecular;
mod wall;

use super::*;
use fs_couple::dg_port::finite_gas::transport::assembly::reactive::thermal::shared::{
    SharedFiniteGasAcousticThermalPort, SharedFiniteGasAcousticThermalSolve,
    advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths,
    finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths,
};

fn shared_solve() -> SharedFiniteGasAcousticThermalSolve {
    SharedFiniteGasAcousticThermalSolve {
        reactive: reactive_solve(),
        body_exchange: heat_exchange(),
        maximum_endpoint_evaluations: 4096,
    }
}

fn shared_body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    initial: EquilibriumPhaseState,
    heat_j: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial,
        mass_kg: BODY_MASS,
        volume_m3: BODY_MASS / 4000.0,
        surface_area_m2: AREA,
        internal_heat_j: heat_j,
        duration_s: DT,
    }
}

fn shared_ports<'a>(
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    context: &'a QueryPoint,
    ids: [u64; 2],
) -> BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'a>> {
    ids.map(|id| {
        (
            id,
            SharedFiniteGasAcousticThermalPort {
                layers: MaterialConductionLayers {
                    material: card,
                    context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                    conductivity: key,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    area_m2: AREA / 2.0,
                    gas_distance_m: GAS_DISTANCE,
                    body_distance_m: BODY_DISTANCE,
                },
                exchange: heat_exchange(),
                wall: None,
            },
        )
    })
    .into()
}

#[test]
fn g1_shared_body_resolves_opposing_heat_even_when_net_exchange_is_zero() {
    let law = reaction(50_000.0);
    let (mut parcels, mut nodes, mut inputs) = fixture(&law);
    let (card, key, caloric) = body_material(false);
    let context = body_point(300.0);
    let mut state = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 100.0)
        .unwrap();
    let body = shared_body(caloric.curve(), state, 0.0);
    let ports = shared_ports(&card, &key, &context, [90, 91]);
    // The same number of pure-I molecules occupies two fixed volumes. Their
    // temperatures straddle the single body, so each exchange is nonzero while
    // their net transfer vanishes. The existing warm acoustic owners continue.
    let original = parcels[&90].clone();
    for (id, temperature) in [(90, 250.0), (91, 350.0)] {
        parcels.insert(
            id,
            ClosedGasParcel::new(
                original.model().clone(),
                point(temperature, 180_000.0),
                original.volume_m3() * temperature / 300.0,
            )
            .unwrap(),
        );
    }
    inputs.insert(
        91,
        FiniteGasAcousticInput::Passive(FiniteGasDgTransportInput {
            gas: energy_solve(),
            supplied_heat_j: 0.0,
            acoustic: None,
            ports: Vec::new(),
        }),
    );
    for input in inputs.values_mut() {
        thermal_mut(input).supplied_heat_j = 0.0;
    }
    let before = parcels.clone();
    let native_before = nodes.clone();
    let gas_capacity = 1.5 * R * original.mass_kg() / 0.004;
    let body_capacity = BODY_MASS * BODY_CP;
    let conductance = (AREA / 2.0) / (GAS_DISTANCE / 0.026 + BODY_DISTANCE / 0.2);
    let alpha = DT * conductance / (1.0 + DT * conductance / gas_capacity);
    let expected_body =
        (body_capacity * 300.0 + alpha * (250.0 + 350.0)) / (body_capacity + 2.0 * alpha);
    let expected_q = [
        alpha * (expected_body - 250.0),
        alpha * (expected_body - 350.0),
    ];
    let report = advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
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
    let mut exchanged = 0.0;
    let mut absolute = 0.0;
    for (id, expected) in [(90, expected_q[0]), (91, expected_q[1])] {
        let port = report.port(id).unwrap();
        let initial_temperature = before[&id].state().state().temperature;
        let final_temperature = parcels[&id].state().state().temperature;
        assert!((port.solid_to_gas_heat_j - expected).abs() < 2e-13);
        assert!((final_temperature - initial_temperature - expected / gas_capacity).abs() < 1e-8);
        assert!(port.solid_to_gas_heat_j * (300.0 - initial_temperature) > 0.0);
        assert!((final_temperature - initial_temperature).abs() > 1e-4);
        assert_eq!(port.body_temperature_k, state.temperature_k());
        close(port.conductance_w_k, conductance, conductance);
        assert_eq!(parcels[&id].mass_kg(), before[&id].mass_kg());
        exchanged += port.solid_to_gas_heat_j;
        absolute += port.solid_to_gas_heat_j.abs();
    }
    assert!(absolute > 1e-6);
    assert!(exchanged.abs() < 1e-6 * absolute);
    assert!((state.temperature_k() - expected_body).abs() < 1e-9);
    assert!(report.solid_to_gas_heat_j.abs() < 1e-12);
    assert!(report.body_storage_change_j.abs() < 1e-11);
    assert!(report.body_energy_residual_j.abs() < shared_solve().body_exchange.energy_tolerance_j);
    assert!(report.endpoint_evaluations >= 4);
    for id in [10, 20] {
        assert_eq!(
            nodes[&id].accepted_steps(),
            native_before[&id].accepted_steps() + 1
        );
        assert!(report.port(id).is_none());
    }
    // Enclosing the shared body once cancels its two signed exchanges.
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy(gas) - energy(&before[id]))
        .sum();
    let native_change: f64 = nodes
        .iter()
        .map(|(id, node)| node_energy(node) - node_energy(&native_before[id]))
        .sum();
    let body_change = body_capacity * (state.temperature_k() - body.initial.temperature_k());
    let (work, losses) = work_and_losses(&report.transport.acoustics, &nodes, &inputs);
    assert!(
        (gas_change + native_change + body_change + losses - work).abs()
            < solve().energy_tolerance_j
    );
}

#[test]
fn g1_shared_material_storage_closes_reacting_transport_and_continuing_sources() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (card, key, caloric) = body_material(true);
    let context = body_point(300.0);
    let mut state = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 300.0)
        .unwrap();
    let mut body = shared_body(caloric.curve(), state, 5e-6);
    let ports = shared_ports(&card, &key, &context, [10, 20]);
    let (before_gas, before_nodes, before_body) = (parcels.clone(), nodes.clone(), state);
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
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
    let (mut isolated_gas, mut isolated_nodes) = (parcels.clone(), nodes.clone());
    advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut isolated_gas,
        &mut isolated_nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
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
    assert_eq!(report.transport.transport, *transported.step());
    assert_eq!(report.body_state, state);
    let mut exchanged = 0.0;
    for id in [10, 20] {
        let port = report.port(id).unwrap();
        let chemical = report.reaction(id).unwrap();
        let q = port.solid_to_gas_heat_j;
        let gas_temperature = parcels[&id].state().state().temperature;
        let body_k = 0.2 * state.temperature_k() / 300.0;
        let conductance = (AREA / 2.0)
            / (GAS_DISTANCE / parcels[&id].state().state().thermal_conductivity
                + BODY_DISTANCE / body_k);
        assert!(q > 0.0);
        assert_eq!(port.body_temperature_k, state.temperature_k());
        close(
            port.material.property_by_key(&key).unwrap().value_si(),
            body_k,
            body_k,
        );
        assert_eq!(port.material.card_identity(), card.content_hash());
        assert!(
            (q - DT * conductance * (state.temperature_k() - gas_temperature)).abs()
                <= heat_exchange().source_tolerance_j
        );
        let work = accepted_work(&report.transport.acoustics[&id]);
        assert!(work.abs() > 1e-18);
        assert!(chemical.reaction.extent_mol > 1e-9);
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
        exchanged += q;
    }
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy(gas) - energy(&before_gas[id]))
        .sum();
    let native_change: f64 = nodes
        .iter()
        .map(|(id, node)| node_energy(node) - node_energy(&before_nodes[id]))
        .sum();
    let body_change = BODY_MASS * BODY_CP * (state.temperature_k() - before_body.temperature_k());
    let external_gas: f64 = inputs
        .values()
        .map(|input| thermal(input).supplied_heat_j)
        .sum();
    let (work, losses) = work_and_losses(&report.transport.acoustics, &nodes, &inputs);
    let residual = gas_change + native_change + body_change + losses
        - external_gas
        - body.internal_heat_j
        - work;
    assert!(residual.abs() < solve().energy_tolerance_j);
    assert!((report.combined_energy_residual_j - residual).abs() < 2e-12);
    assert!((report.body_storage_change_j - body_change).abs() < 1e-12);
    assert!(
        (body_change + exchanged - body.internal_heat_j).abs() < heat_exchange().energy_tolerance_j
    );
    assert_eq!(report.external_body_heat_j, body.internal_heat_j);
    assert!((report.solid_to_gas_heat_j - exchanged).abs() < 1e-15);
    assert!((report.transport.supplied_heat_j - external_gas - exchanged).abs() < 1e-15);
    assert!(report.source_residual_j.abs() <= shared_solve().body_exchange.source_tolerance_j);
    assert!(report.maximum_species_residual_mol <= reactive_solve().species_tolerance_mol);
    assert!(report.endpoint_evaluations <= shared_solve().maximum_endpoint_evaluations);

    body.initial = state;
    advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
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
    advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut isolated_gas,
        &mut isolated_nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert!(state.temperature_k() < body.initial.temperature_k());
    assert_eq!(nodes[&20].accepted_steps(), 26);
    let (FiniteGasAcousticNode::Aperture(a), FiniteGasAcousticNode::Aperture(b)) =
        (&nodes[&20], &isolated_nodes[&20])
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
fn g4_shared_body_refusals_budget_and_cancellation_preserve_one_accepted_state() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (card, key, caloric) = body_material(true);
    let context = body_point(300.0);
    let mut state = caloric
        .curve()
        .state_at_specific_enthalpy(BODY_CP * 300.0)
        .unwrap();
    let body = shared_body(caloric.curve(), state, 5e-6);
    let ports = shared_ports(&card, &key, &context, [10, 20]);
    let (before_gas, before_nodes, before_body) = (parcels.clone(), nodes.clone(), state);
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let preview = finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths(
        &parcels,
        &nodes,
        &state,
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
    assert_eq!(parcels, before_gas);
    equal_nodes(&nodes, &before_nodes);
    assert_eq!(state, before_body);
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
    assert_eq!(preview.parcels(), &accepted_gas);
    equal_nodes(preview.nodes(), &accepted_nodes);
    assert_eq!(preview.body_state(), &accepted_body);
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
    equal_nodes(&nodes, &before_nodes);
    assert_eq!(state, before_body);

    let mut coverage = shared_ports(&card, &key, &context, [10, 20]);
    for port in coverage.values_mut() {
        // Each patch separately fits; their sum does not fit the one body.
        port.layers.area_m2 = 0.75 * AREA;
    }
    assert!(
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut state,
            body,
            &inputs,
            &reactions,
            &coverage,
            paths,
            DT,
            shared_solve(),
            &mut || false,
        )
        .is_err()
    );
    let (foreign_card, foreign_key, _) = body_material(false);
    let mut foreign = shared_ports(&card, &key, &context, [10, 20]);
    foreign.get_mut(&20).unwrap().layers.material = &foreign_card;
    foreign.get_mut(&20).unwrap().layers.conductivity = &foreign_key;
    assert!(
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut state,
            body,
            &inputs,
            &reactions,
            &foreign,
            paths,
            DT,
            shared_solve(),
            &mut || false,
        )
        .is_err()
    );
    let mismatch = UniformEnthalpyStepInput {
        initial: caloric
            .curve()
            .state_at_specific_enthalpy(BODY_CP * 200.0)
            .unwrap(),
        ..body
    };
    assert!(
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut state,
            mismatch,
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
    let exhausted = SharedFiniteGasAcousticThermalSolve {
        maximum_endpoint_evaluations: 1,
        ..shared_solve()
    };
    assert!(
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut state,
            body,
            &inputs,
            &reactions,
            &ports,
            paths,
            DT,
            exhausted,
            &mut || false,
        )
        .is_err()
    );
    // A valid original chart/source cannot support this proposed heated state.
    let overheated = UniformEnthalpyStepInput {
        internal_heat_j: 100.0,
        ..body
    };
    assert!(
        advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &mut state,
            overheated,
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
    assert_eq!(parcels, before_gas);
    equal_nodes(&nodes, &before_nodes);
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
    equal_nodes(&nodes, &accepted_nodes);
    assert_eq!(state, accepted_body);
    assert_eq!(retry.transport.transport, accepted.transport.transport);
    assert_eq!(retry.body_storage_change_j, accepted.body_storage_change_j);
    assert_eq!(
        retry.combined_energy_residual_j,
        accepted.combined_energy_residual_j
    );
    assert_eq!(retry.endpoint_evaluations, accepted.endpoint_evaluations);
    for id in [10, 20] {
        assert_eq!(retry.reaction(id), accepted.reaction(id));
        assert_eq!(
            retry.port(id).unwrap().material,
            accepted.port(id).unwrap().material
        );
        assert_eq!(
            retry.port(id).unwrap().solid_to_gas_heat_j,
            accepted.port(id).unwrap().solid_to_gas_heat_j
        );
    }
}
