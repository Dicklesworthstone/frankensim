//! Native transport, finite channel momentum and causal acoustics share one endpoint.
use super::*;
use fs_couple::dg_port::finite_gas::transport::{
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    advance_finite_gas_port_transport, advance_finite_gas_port_transport_with_paths,
    finite_gas_port_transport_candidate_with_paths,
};
use fs_duct::dg::PortInput;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::{SnapshotNozzleError, SnapshotNozzleLink},
    species::{
        diffusion::{
            BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
            GasDiffusionGeometry,
            network::{DiffusionFaceSolve, DiffusionNetworkLink},
        },
        network::{
            SpeciesNetworkSolve,
            transport::{
                self as material_transport, GasTransportError, GasTransportPaths, GasTransportSolve,
            },
        },
    },
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
use std::collections::BTreeMap;

type Parcels = BTreeMap<u64, ClosedGasParcel>;
type Loads = BTreeMap<u64, DgPortLoad>;

fn network() -> (Parcels, Loads, BTreeMap<u64, FiniteGasDgTransportInput>) {
    let (slit, first) = primed(false);
    let (round, second) = primed(true);
    let reservoir = ClosedGasParcel::new(
        first.model().with_composition(composition(0.1)).unwrap(),
        point(300.0, 180_000.0),
        6.0 * volume(false),
    )
    .unwrap();
    let parcels = BTreeMap::from([(10, first), (20, second), (90, reservoir)]);
    let nodes = BTreeMap::from([(10, slit), (20, round)]);
    let inputs = parcels
        .keys()
        .map(|&id| {
            (
                id,
                FiniteGasDgTransportInput {
                    gas: gas_solve(),
                    supplied_heat_j: match id {
                        10 => 2e-5,
                        20 => -1e-5,
                        _ => 3e-6,
                    },
                    acoustic: nodes.contains_key(&id).then_some(controls()),
                    ports: nodes.get(&id).map_or_else(Vec::new, |load| {
                        vec![PortInput {
                            id: load.inlet(),
                            effort: if id == 10 { FLOW } else { -0.2 * FLOW },
                        }]
                    }),
                },
            )
        })
        .collect();
    (parcels, nodes, inputs)
}

fn solve() -> FiniteGasDgTransportSolve {
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

fn nozzles() -> [SnapshotNozzleLink; 2] {
    [(9, 10), (31, 20)].map(|(id, receiver)| SnapshotNozzleLink {
        id,
        parcels: [90, receiver],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: gas_solve(),
    })
}

fn species_mass(parcels: &Parcels, name: &str) -> f64 {
    let species = SpeciesId::new(name).unwrap();
    parcels
        .values()
        .map(|parcel| {
            let law = parcel.model().thermodynamics();
            let masses: Vec<_> = law.components().iter().map(|c| c.molar_mass()).collect();
            let fractions = law.mole_composition().to_mass_fractions(&masses).unwrap();
            law.components()
                .iter()
                .zip(fractions.fractions())
                .filter(|(component, _)| component.species() == &species)
                .map(|(_, fraction)| parcel.mass_kg() * fraction)
                .sum::<f64>()
        })
        .sum()
}

fn unchanged_loads(nodes: &Loads, original: &Loads) {
    assert!(nodes.keys().eq(original.keys()));
    for (&id, load) in nodes {
        unchanged(load, &original[&id]);
    }
}

#[test]
fn g1_one_native_snapshot_drives_two_complete_loads_and_a_finite_reservoir() {
    let (mut parcels, mut nodes, inputs) = network();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
    let endpoints = inputs.iter().map(|(&id, input)| (id, input.gas)).collect();
    let reference = material_transport::step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        },
        DT,
        &endpoints,
        solve().transport,
        || false,
    )
    .unwrap();
    let report = advance_finite_gas_port_transport(
        &mut parcels,
        &mut nodes,
        &inputs,
        &links,
        &[],
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport, *reference.step());
    assert_eq!(
        report.reservoirs[&90].supplied_heat_j,
        inputs[&90].supplied_heat_j
    );
    for link in &links {
        let row = report
            .transport
            .nozzles
            .iter()
            .find(|row| row.id == link.id)
            .unwrap();
        let native = old_gas[&90]
            .nozzle_flow_to(
                old_gas[&link.parcels[1]].state().state().pressure,
                link.nozzle,
                gas_solve(),
                || false,
            )
            .unwrap();
        assert_eq!(row.flow, native);
        assert_eq!(row.transferred_mass_kg, native.mass_flow_kg_s * DT);
    }
    for name in ["A", "B"] {
        assert!((species_mass(&parcels, name) - species_mass(&old_gas, name)).abs() < 1e-16);
    }
    let mut inlet_work = 0.0;
    let mut exterior_work = 0.0;
    let mut losses = 0.0;
    for (&id, load) in &nodes {
        let (mut marched, mut frozen_gas) = (old_nodes[&id].clone(), old_gas[&id].clone());
        let flow = inputs[&id].ports[0].effort;
        let frozen = frozen_step(&mut marched, &mut frozen_gas, flow);
        let step = &report.acoustics[&id];
        assert_eq!(
            load.terminal().unwrap().pole_states(),
            marched.terminal().unwrap().pole_states()
        );
        assert_eq!(load.last_terminal_frames(), marched.last_terminal_frames());
        assert_eq!(
            load.last_terminal_channel_frames(),
            marched.last_terminal_channel_frames()
        );
        assert_eq!(
            step.acoustic.ledger.port_work_step_j,
            frozen.acoustic.ledger.port_work_step_j
        );
        assert_eq!(
            load.state().step_count(),
            old_nodes[&id].state().step_count() + 1
        );
        assert_eq!(
            load.duct().gas_field().unwrap().uniform_state(),
            Some(parcels[&id].state())
        );
        assert_eq!(
            load.terminal_channels()[0].gas(),
            parcels[&id].state().state()
        );
        assert_ne!(
            parcels[&id].model().thermodynamics().mole_composition(),
            old_gas[&id].model().thermodynamics().mole_composition()
        );
        // Momentum in each actual channel cell is invariant through the source
        // refresh. Its independently integrated kinetic work belongs here too.
        let old = marched.terminal_channels()[0].channel().unwrap();
        let next = load.terminal_channels()[0].channel().unwrap();
        let (rho0, rho1) = (old.gas().density, next.gas().density);
        let (mut kinetic_before, mut kinetic_after) = (0.0, 0.0);
        for (i, (&a, &b)) in old
            .velocities_m_s()
            .iter()
            .zip(next.velocities_m_s())
            .enumerate()
        {
            let cell = if id == 20 {
                let inner = CHANNEL_RADIUS * i as f64 / CELLS as f64;
                let outer = CHANNEL_RADIUS * (i + 1) as f64 / CELLS as f64;
                core::f64::consts::PI * CHANNEL_LENGTH * (outer * outer - inner * inner)
            } else {
                channel_volume(false) / CELLS as f64
            };
            near(rho1 * cell * b, rho0 * cell * a, (rho0 * cell * a).abs());
            kinetic_before += 0.5 * rho0 * cell * a * a;
            kinetic_after += 0.5 * rho1 * cell * b * b;
        }
        near(kinetic_before, old.stored_energy_j(), kinetic_before);
        near(kinetic_after, kinetic_before * rho0 / rho1, kinetic_before);
        let channel_work = step.refresh.terminal_channels[0]
            .transfer
            .unwrap()
            .parameter_work_j;
        near(channel_work, kinetic_after - kinetic_before, kinetic_before);
        assert!(channel_work != 0.0);
        inlet_work += flow
            * frozen
                .acoustic
                .port_outputs
                .iter()
                .find(|port| port.id == load.inlet())
                .unwrap()
                .delta_x;
        exterior_work += load
            .last_terminal_frames()
            .iter()
            .map(|frame| frame.load_work_j)
            .sum::<f64>();
        losses += load.cumulative_dissipated_j() - old_nodes[&id].cumulative_dissipated_j();
    }
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy_change(&old_gas[id], gas))
        .sum();
    let acoustic_change: f64 = nodes
        .iter()
        .map(|(id, load)| load.stored_energy_j() - old_nodes[id].stored_energy_j())
        .sum();
    let heat: f64 = inputs.values().map(|input| input.supplied_heat_j).sum();
    let residual = gas_change + acoustic_change + losses + exterior_work - inlet_work - heat;
    assert!(residual.abs() <= solve().energy_tolerance_j, "{residual:e}");
    assert!(report.combined_energy_residual_j.abs() <= solve().energy_tolerance_j);
    near(
        report.acoustic_storage_change_j,
        acoustic_change,
        acoustic_change.abs(),
    );
    assert!(report.total_mass_residual_kg.abs() <= solve().transport.exchange.mass_tolerance_kg);
}

fn diffusivity(value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                200.0,
                900.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic binary diffusion reference".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "A/B pair".into(),
            phase: "gas".into(),
            process: "manufactured transport fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn g3_source_diffusion_matches_binary_limit_and_changes_the_continuing_sound() {
    let mut traces = Vec::new();
    for coefficient in [2e-5, 4e-5] {
        let (mut parcels, mut nodes, mut inputs) = network();
        // Retain both real loads; only the slit chamber exchanges species here.
        let reservoir = ClosedGasParcel::new(
            parcels[&90].model().clone(),
            point(300.0, 101325.0),
            parcels[&90].volume_m3(),
        )
        .unwrap();
        parcels.insert(90, reservoir);
        for input in inputs.values_mut() {
            input.supplied_heat_j = 0.0;
        }
        let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
        let card = diffusivity(coefficient);
        let context = QueryPoint::new();
        let pairs = [BinaryDiffusionMaterial {
            species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            material: &card,
            context: &context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }];
        let faces = [DiffusionNetworkLink {
            id: 41,
            parcels: [10, 90],
            pairs: &pairs,
            geometry: GasDiffusionGeometry {
                area_m2: 2e-6,
                distance_m: 1e-5,
            },
            solve: DiffusionFaceSolve {
                maximum_temperature_difference_k: 1.0,
                maximum_pressure_difference_pa: 200.0,
                gradient_tolerance_per_m: 1e-7,
                molar_flux_tolerance_mol_m2_s: 1e-10,
            },
        }];
        let paths = GasTransportPaths {
            diffusion: &faces,
            ..Default::default()
        };
        let preview = finite_gas_port_transport_candidate_with_paths(
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
        unchanged_loads(&nodes, &old_nodes);
        let flux = &preview.step().transport.diffusion[0].flux;
        let expected =
            101325.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 300.0) * coefficient * (0.8 - 0.1)
                / faces[0].geometry.distance_m;
        near(flux.species[0].molar_flux_mol_m2_s, expected, expected);
        near(flux.species[1].molar_flux_mol_m2_s, -expected, expected);
        assert_eq!(flux.pairs[0].diffusivity_m2_s, coefficient);
        assert_eq!(flux.pairs[0].material.card_identity(), card.content_hash());
        let (expected_gas, expected_nodes, expected_step) = preview.into_parts();
        let report = advance_finite_gas_port_transport_with_paths(
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
        unchanged_loads(&nodes, &expected_nodes);
        assert_eq!(report.transport, expected_step.transport);
        for name in ["A", "B"] {
            assert!((species_mass(&parcels, name) - species_mass(&old_gas, name)).abs() < 1e-16);
        }
        let gas = parcels.get_mut(&10).unwrap();
        let load = nodes.get_mut(&10).unwrap();
        assert_ne!(
            gas.state().state().sound_speed,
            old_gas[&10].state().state().sound_speed
        );
        assert_eq!(load.terminal_channels()[0].gas(), gas.state().state());
        let mut trace = Vec::new();
        for _ in 0..8 {
            let step = frozen_step(load, gas, 0.0);
            trace.extend(
                step.acoustic
                    .port_outputs
                    .iter()
                    .map(|port| port.end_conjugate),
            );
        }
        traces.push(trace);
    }
    let difference: f64 = traces[0]
        .iter()
        .zip(&traces[1])
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    let signal: f64 = traces[0].iter().map(|p| p * p).sum();
    assert!(difference > 1e-20 * signal);
}

#[test]
fn g4_late_transport_refusal_and_final_cancel_preserve_all_owners_and_exact_retry() {
    let (mut parcels, mut nodes, inputs) = network();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
    for ports in [
        vec![],
        vec![PortInput {
            id: PortId(999),
            effort: FLOW,
        }],
    ] {
        let mut bad = inputs.clone();
        bad.get_mut(&20).unwrap().ports = ports;
        assert!(matches!(
            advance_finite_gas_port_transport(
                &mut parcels,
                &mut nodes,
                &bad,
                &links,
                &[],
                DT,
                solve(),
                &mut || false,
            ),
            Err(FiniteGasDgTransportError::Acoustic {
                parcel: 20,
                error: FiniteGasDgError::Invalid(_),
            })
        ));
        assert_eq!(parcels, old_gas);
        unchanged_loads(&nodes, &old_nodes);
    }
    let mut crossing_gas = parcels.clone();
    crossing_gas.insert(
        90,
        ClosedGasParcel::new(
            parcels[&90].model().clone(),
            point(300.0, 100000.0),
            parcels[&90].volume_m3(),
        )
        .unwrap(),
    );
    let original_crossing = crossing_gas.clone();
    let mut hot = inputs.clone();
    hot.get_mut(&90).unwrap().supplied_heat_j = 0.05;
    let endpoints = hot.iter().map(|(&id, input)| (id, input.gas)).collect();
    material_transport::step_candidate_with_paths(
        &crossing_gas,
        GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        },
        DT,
        &endpoints,
        solve().transport,
        || false,
    )
    .unwrap();
    let error = advance_finite_gas_port_transport(
        &mut crossing_gas,
        &mut nodes,
        &hot,
        &links,
        &[],
        DT,
        solve(),
        &mut || false,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            FiniteGasDgTransportError::Transport(GasTransportError::Nozzle(
                SnapshotNozzleError::EquilibriumCrossed { .. }
            ))
        ),
        "{error:?}"
    );
    assert_eq!(crossing_gas, original_crossing);
    unchanged_loads(&nodes, &old_nodes);
    let (mut expected_gas, mut expected_nodes) = (parcels.clone(), nodes.clone());
    let mut polls = 0;
    let expected = advance_finite_gas_port_transport(
        &mut expected_gas,
        &mut expected_nodes,
        &inputs,
        &links,
        &[],
        DT,
        solve(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        advance_finite_gas_port_transport(
            &mut parcels,
            &mut nodes,
            &inputs,
            &links,
            &[],
            DT,
            solve(),
            &mut || {
                count += 1;
                count == polls
            },
        ),
        Err(FiniteGasDgTransportError::Cancelled)
    ));
    assert_eq!(count, polls);
    assert_eq!(parcels, old_gas);
    unchanged_loads(&nodes, &old_nodes);
    let retry = advance_finite_gas_port_transport(
        &mut parcels,
        &mut nodes,
        &inputs,
        &links,
        &[],
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, expected_gas);
    unchanged_loads(&nodes, &expected_nodes);
    assert_eq!(retry.transport, expected.transport);
    for (&id, step) in &retry.acoustics {
        assert_eq!(step.acoustic, expected.acoustics[&id].acoustic);
        assert_eq!(step.refresh, expected.acoustics[&id].refresh);
        assert_eq!(step.work_by_gas_j, expected.acoustics[&id].work_by_gas_j);
    }
}
