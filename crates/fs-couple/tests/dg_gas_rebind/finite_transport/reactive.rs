//! Original chemistry changes the same transported gas and continuing acoustic owners.
//! Authored isomer/transport laws test implementation, not measured combustion.
#[path = "reactive/thermal.rs"]
mod thermal;

use super::*;
use fs_couple::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, boundary::GasReactionBoundary,
    implicit::coupled::ImplicitMechanismSolve,
};
use fs_couple::dg_port::finite_gas::transport::assembly::{
    FiniteGasAcousticInput, FiniteGasAcousticNode, FiniteGasAcousticStep,
    advance_finite_gas_acoustic_transport_with_paths,
    reactive::{
        FiniteGasAcousticReaction, FiniteGasAcousticReactionReservoir,
        FiniteGasAcousticReactionSource, ReactiveFiniteGasAcousticTransportError,
        ReactiveFiniteGasAcousticTransportSolve,
        advance_reacting_finite_gas_acoustic_transport_with_paths,
        finite_gas_reacting_acoustic_transport_candidate_with_paths,
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::{MolarMass, Pressure};
use fs_thermochem::equilibrium::{
    DeclaredStoichiometricMassActionV1, IdealGasReactionEquilibriumV1,
    MassActionKineticOrderConventionV1, ReactionProgressRateQuantityV1,
    ReactionProgressRateScaleV1, ThermodynamicReverseRateClosureV1,
};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReactionId, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, StoichiometricMatrix, required_nasa9_parameter_dims, verify_conservation,
};

type MixedNodes = BTreeMap<u64, FiniteGasAcousticNode>;
type MixedInputs = BTreeMap<u64, FiniteGasAcousticInput>;
const RELEASE: f64 = 2000.0;

fn mixture(fractions: [f64; 3]) -> Composition {
    Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap()
}

fn reaction(scale: f64) -> GasReaction {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.028, 3.5, 200.0);
    let inert = gas_support::packs("I", 0.004, 2.5, 200.0);
    let transport = from_packs(
        &[a.binding(), b.binding(), inert.binding()],
        mixture([0.7, 0.05, 0.25]),
        None,
    )
    .unwrap();
    // Constant Cp and a nonzero B formation energy make the total-U oracle
    // independent of the production NASA evaluator and catch duplicate heat.
    let species = [
        ("A", 0.028, 3.5, 0.0),
        ("B", 0.028, 3.5, -RELEASE),
        ("I", 0.004, 2.5, 0.0),
    ]
    .map(|(name, mass, cp_r, formation)| {
        let mut coefficients = [0.0; 9];
        coefficients[2] = cp_r;
        coefficients[7] = formation / R;
        let mut parameters: BTreeMap<_, _> = required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .zip(coefficients)
            .map(|((name, dims), value)| (name.to_owned(), LawParameter { value, dims }))
            .collect();
        parameters.insert(
            "reference_pressure".into(),
            LawParameter {
                value: 100_000.0,
                dims: Pressure::DIMS,
            },
        );
        let card = ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
            sources: vec![],
            provenance: Provenance {
                source: "authored balanced isomers; not empirical kinetics".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        Nasa9StandardStateModelV1::new(
            SpeciesId::new(name).unwrap(),
            MolarMass::new(mass),
            StandardStateConventionV1::new(
                StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas,
                Pressure::new(100_000.0),
                ElementalReferenceIdV1::new("synthetic-isomer-reference").unwrap(),
            )
            .unwrap(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap()
    });
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(species.to_vec(), mixture([0.7, 0.05, 0.25])).unwrap(),
        transport.transport_materials().to_vec(),
        None,
    )
    .unwrap();
    let ids = vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let id = ReactionId::new("balanced-acoustic-isomer").unwrap();
    let atoms = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        ids.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let stoichiometry =
        StoichiometricMatrix::new(ids.clone(), vec![id.clone()], vec![vec![-1], vec![1]]).unwrap();
    let charge = ChargeVector::new(ids, vec![0, 0]).unwrap();
    let balance = verify_conservation(&atoms, &stoichiometry, &charge).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(
        stoichiometry,
        balance,
        id.clone(),
        catalog.thermodynamics().components()[..2].to_vec(),
    )
    .unwrap();
    let kinetics = DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(scale)).unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    );
    GasReaction::new(catalog, kinetics)
        .unwrap()
        .with_implicit_mechanism(ImplicitMechanismSolve {
            extent_brackets_mol: BTreeMap::from([(id, [-1e-6, 1e-6])]),
            extent_tolerance_mol: 1e-15,
            interior_fraction: 0.05,
            jacobian_relative_step: 1e-5,
            maximum_iterations: 32,
            maximum_evaluations: 1024,
            newton: fs_solver::NewtonKrylovConfig {
                absolute_tolerance: 1e-13,
                relative_tolerance: 1e-12,
                linear_restart: 8,
                max_linear_cycles: 8,
                forcing_maximum: 1e-6,
                ..Default::default()
            },
        })
        .unwrap()
}

fn chemistry() -> GasReactionSolve {
    GasReactionSolve {
        gas: energy_solve(),
        maximum_consumed_fraction: 0.15,
        species_tolerance_mol: 1e-14,
        entropy_tolerance_j_k: 1e-12,
    }
}
fn reactive_solve() -> ReactiveFiniteGasAcousticTransportSolve {
    ReactiveFiniteGasAcousticTransportSolve {
        transport: solve(),
        species_tolerance_mol: 1e-13,
    }
}
fn bindings(law: &GasReaction) -> BTreeMap<u64, FiniteGasAcousticReaction<'_>> {
    // The pure-I reservoir deliberately exercises an inert zero-extent endpoint.
    [10, 20, 90]
        .map(|id| {
            (
                id,
                FiniteGasAcousticReaction {
                    reaction: law,
                    controls: chemistry(),
                },
            )
        })
        .into()
}
fn field_for(law: &GasReaction) -> AxialGasField {
    AxialGasField::new(
        vec![0.0, LENGTH],
        101_325.0,
        vec![
            PrimitiveGasState {
                temperature_k: 300.0,
                composition: mixture([0.7, 0.05, 0.25])
            };
            2
        ],
        law.catalog().clone(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}
fn fixture(law: &GasReaction) -> (Parcels, MixedNodes, MixedInputs) {
    let field = field_for(law);
    let load = channel_load(&field);
    let mut first =
        LoadedApertureDg::new_loaded(aperture(&load, load.impedance_pa_s_m3()), load).unwrap();
    let load = channel_load(&field);
    let mut second =
        LoadedApertureDg::new_loaded(paired_source(&load, load.impedance_pa_s_m3()), load).unwrap();
    for sample in 0..24 {
        first.step(drive(sample)).unwrap();
        second.step(drive(sample)).unwrap();
    }
    let volume = core::f64::consts::PI * RADIUS * RADIUS * LENGTH + 0.002 * 4e-5 * 0.001;
    let gas = ClosedGasParcel::new(law.catalog().clone(), point(300.0, 101_325.0), volume).unwrap();
    let reservoir = ClosedGasParcel::new(
        law.catalog()
            .with_composition(mixture([0.0, 0.0, 1.0]))
            .unwrap(),
        point(300.0, 180_000.0),
        6.0 * volume,
    )
    .unwrap();
    let thermal = |acoustic: bool| FiniteGasDgTransportInput {
        gas: energy_solve(),
        supplied_heat_j: 2e-5,
        acoustic: acoustic.then_some(controls()),
        ports: Vec::new(),
    };
    let mut passive = thermal(true);
    passive.ports.push(fs_duct::dg::PortInput {
        id: first.load().inlet(),
        effort: 1e-8,
    });
    (
        BTreeMap::from([(10, gas.clone()), (20, gas), (90, reservoir)]),
        BTreeMap::from([
            (10, FiniteGasAcousticNode::Port(first.load().clone())),
            (20, FiniteGasAcousticNode::Aperture(second)),
        ]),
        BTreeMap::from([
            (10, FiniteGasAcousticInput::Passive(passive)),
            (
                20,
                FiniteGasAcousticInput::Aperture(FiniteGasApertureTransportInput {
                    thermal: thermal(true),
                    source: Some(FiniteGasApertureDrive {
                        drive: TubeDrive {
                            upstream_pressure_pa: 20.0,
                            body_flow_m3_s: -1e-9,
                        },
                        forces_n: [2e-5, -3e-5],
                        pressures_pa: [0.0; 2],
                    }),
                }),
            ),
            (90, FiniteGasAcousticInput::Passive(thermal(false))),
        ]),
    )
}
fn thermal(input: &FiniteGasAcousticInput) -> &FiniteGasDgTransportInput {
    match input {
        FiniteGasAcousticInput::Passive(row) => row,
        FiniteGasAcousticInput::MolecularWall { thermal, .. } => thermal,
        FiniteGasAcousticInput::Aperture(row) => &row.thermal,
    }
}
fn thermal_mut(input: &mut FiniteGasAcousticInput) -> &mut FiniteGasDgTransportInput {
    match input {
        FiniteGasAcousticInput::Passive(row) => row,
        FiniteGasAcousticInput::MolecularWall { thermal, .. } => thermal,
        FiniteGasAcousticInput::Aperture(row) => &mut row.thermal,
    }
}
fn equal_nodes(actual: &MixedNodes, expected: &MixedNodes) {
    assert!(actual.keys().eq(expected.keys()));
    for (&id, node) in actual {
        assert_eq!(node.accepted_steps(), expected[&id].accepted_steps());
        assert_eq!(node.stored_energy_j(), expected[&id].stored_energy_j());
        match (node, &expected[&id]) {
            (FiniteGasAcousticNode::Port(a), FiniteGasAcousticNode::Port(b)) => {
                assert_eq!(a.state(), b.state());
                assert_eq!(a.duct().gas_field(), b.duct().gas_field());
                assert_eq!(
                    a.terminal().map(|x| x.pole_states()),
                    b.terminal().map(|x| x.pole_states())
                );
                assert_eq!(a.last_terminal_frames(), b.last_terminal_frames());
                assert_eq!(
                    a.last_terminal_channel_frames(),
                    b.last_terminal_channel_frames()
                );
                assert_eq!(a.terminal_channels().len(), b.terminal_channels().len());
                for (x, y) in a.terminal_channels().iter().zip(b.terminal_channels()) {
                    assert_eq!(
                        x.channel().map(|c| c.velocities_m_s()),
                        y.channel().map(|c| c.velocities_m_s())
                    );
                }
            }
            (FiniteGasAcousticNode::Aperture(a), FiniteGasAcousticNode::Aperture(b)) => unchanged(
                &BTreeMap::from([(id, a.clone())]),
                &BTreeMap::from([(id, b.clone())]),
            ),
            _ => panic!("fixture changed native owner type"),
        }
    }
}
// Independent molar inventory and formation-inclusive constant-Cp NASA integral.
fn moles(gas: &ClosedGasParcel) -> [f64; 3] {
    let model = gas.model().thermodynamics();
    let mut fractions = [0.0; 3];
    for (index, name) in ["A", "B", "I"].into_iter().enumerate() {
        if let Some(position) = model
            .components()
            .iter()
            .position(|s| s.species() == &SpeciesId::new(name).unwrap())
        {
            fractions[index] = model.mole_composition().fractions()[position];
        }
    }
    let mean_mass = 0.028 * (fractions[0] + fractions[1]) + 0.004 * fractions[2];
    fractions.map(|x| gas.mass_kg() * x / mean_mass)
}
fn energy(gas: &ClosedGasParcel) -> f64 {
    let n = moles(gas);
    R * (2.5 * (n[0] + n[1]) + 1.5 * n[2]) * gas.state().state().temperature - RELEASE * n[1]
}
fn node_energy(node: &FiniteGasAcousticNode) -> f64 {
    match node {
        FiniteGasAcousticNode::Aperture(node) => complete_energy(node),
        FiniteGasAcousticNode::Port(load) => {
            let gas = load.duct().inlet_gas();
            let kinetic: f64 = load.terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s()
                .iter()
                .map(|v| 0.5 * gas.density * (0.002 * 4e-5 * 0.001) / 8.0 * v * v)
                .sum();
            wave_energy(load.state(), gas, LENGTH) + kinetic
        }
        _ => unreachable!(),
    }
}
fn accepted_work<S>(step: &FiniteGasAcousticStep<S>) -> f64 {
    match step {
        FiniteGasAcousticStep::Port(step) => step.work_by_gas_j,
        FiniteGasAcousticStep::Aperture(step) => step.gas.work_by_gas_j,
        FiniteGasAcousticStep::Wall(step) => step.work_by_gas_j,
    }
}

#[test]
fn g1_reacting_transport_drives_original_sources_and_closes_formation_energy() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (original_gas, original_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
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
    let (mut frozen_gas, mut frozen_nodes) = (parcels.clone(), nodes.clone());
    advance_finite_gas_acoustic_transport_with_paths(
        &mut frozen_gas,
        &mut frozen_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
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
    let mut extent = 0.0;
    let mut net_work = 0.0;
    let mut losses = 0.0;
    for id in [10, 20, 90] {
        let chemical = report.reaction(id).unwrap();
        let work = report
            .transport
            .acoustics
            .get(&id)
            .map_or(0.0, accepted_work);
        assert_eq!(
            chemical.boundary.heat_j,
            thermal(&inputs[&id]).supplied_heat_j
        );
        assert_eq!(chemical.boundary.work_by_gas_j, work);
        let original = law
            .boundary_implicit_candidate(
                &transported.parcels()[&id],
                DT,
                GasReactionBoundary {
                    final_volume_m3: transported.parcels()[&id].volume_m3(),
                    heat_j: thermal(&inputs[&id]).supplied_heat_j,
                    work_by_gas_j: work,
                },
                chemistry(),
                || false,
            )
            .unwrap();
        assert_eq!(original.gas(), &parcels[&id]);
        assert_eq!(original.step(), chemical);
        let t = parcels[&id].state().state().temperature;
        assert_eq!(
            chemical.endpoint_heat_entropy_j_k,
            chemical.boundary.heat_j / t
        );
        assert!(chemical.entropy_residual_j_k >= -chemistry().entropy_tolerance_j_k);
        if id == 90 {
            assert_eq!(chemical.reaction.extent_mol, 0.0);
            assert_eq!(moles(&parcels[&id])[0], 0.0);
            continue;
        }
        assert!(chemical.reaction.extent_mol > 1e-9);
        assert!(work.abs() > 1e-18);
        extent += chemical.reaction.extent_mol;
        assert_ne!(parcels[&id], frozen_gas[&id]);
        assert_eq!(nodes[&id].accepted_steps(), 25);
        match (&nodes[&id], &report.transport.acoustics[&id]) {
            (FiniteGasAcousticNode::Port(load), FiniteGasAcousticStep::Port(step)) => {
                let inlet = step
                    .acoustic
                    .port_outputs
                    .iter()
                    .find(|p| p.id == load.inlet())
                    .unwrap();
                net_work += thermal(&inputs[&id]).ports[0].effort * inlet.delta_x
                    - load
                        .last_terminal_frames()
                        .iter()
                        .map(|t| t.load_work_j)
                        .sum::<f64>();
                losses += step.acoustic_dissipation_j;
                assert_eq!(
                    load.duct().gas_field().unwrap().uniform_state(),
                    Some(parcels[&id].state())
                );
            }
            (FiniteGasAcousticNode::Aperture(node), FiniteGasAcousticStep::Aperture(step)) => {
                let frame = &step.coupled;
                net_work +=
                    frame.upstream_work_j + frame.body_work_j + frame.aperture.mechanical_work_j
                        - node
                            .load()
                            .last_terminal_frames()
                            .iter()
                            .map(|t| t.load_work_j)
                            .sum::<f64>();
                losses += frame.inferred_total_dissipation_j;
                assert_eq!(
                    node.load().duct().gas_field().unwrap().uniform_state(),
                    Some(parcels[&id].state())
                );
            }
            _ => panic!("native physical report changed owner"),
        }
    }
    for (species, chemical_change) in [-extent, extent, 0.0].into_iter().enumerate() {
        let before: f64 = original_gas.values().map(|g| moles(g)[species]).sum();
        let after: f64 = parcels.values().map(|g| moles(g)[species]).sum();
        assert!((after - before - chemical_change).abs() < 1e-14);
    }
    assert!(report.maximum_species_residual_mol < reactive_solve().species_tolerance_mol);
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, g)| energy(g) - energy(&original_gas[id]))
        .sum();
    let acoustic_change: f64 = nodes
        .iter()
        .map(|(id, n)| node_energy(n) - node_energy(&original_nodes[id]))
        .sum();
    let heat: f64 = inputs.values().map(|i| thermal(i).supplied_heat_j).sum();
    let independent = gas_change + acoustic_change + losses - heat - net_work;
    assert!(independent.abs() < solve().energy_tolerance_j);
    close(report.transport.port_work_j, net_work, net_work.abs());
    close(
        report.transport.acoustic_storage_change_j,
        acoustic_change,
        acoustic_change.abs(),
    );
    assert!((report.transport.combined_energy_residual_j - independent).abs() < 2e-12);
    // A continuing second sample responds to the already changed native source
    // gas; a one-sample frozen-medium result alone cannot establish that feedback.
    advance_reacting_finite_gas_acoustic_transport_with_paths(
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
    advance_finite_gas_acoustic_transport_with_paths(
        &mut frozen_gas,
        &mut frozen_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    let (FiniteGasAcousticNode::Aperture(a), FiniteGasAcousticNode::Aperture(b)) =
        (&nodes[&20], &frozen_nodes[&20])
    else {
        panic!("paired source required")
    };
    let (LoadedAperture::Paired(a), LoadedAperture::Paired(b)) = (a.aperture(), b.aperture())
    else {
        panic!("original paired owner required")
    };
    assert_ne!(a.blade_states(), b.blade_states());
}

#[test]
fn g3_empty_reaction_bindings_reduce_exactly_to_original_mixed_transport() {
    let law = reaction(50_000.0);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (mut old_parcels, mut old_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let original = advance_finite_gas_acoustic_transport_with_paths(
        &mut old_parcels,
        &mut old_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    let reactive = advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        &BTreeMap::new(),
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, old_parcels);
    equal_nodes(&nodes, &old_nodes);
    let report = &reactive.transport;
    assert_eq!(report.transport, original.transport);
    assert_eq!(
        [
            report.gas_storage_change_j,
            report.acoustic_storage_change_j,
            report.supplied_heat_j,
            report.port_work_j,
            report.acoustic_dissipation_j,
            report.parameter_work_j,
            report.total_mass_residual_kg,
            report.combined_energy_residual_j
        ],
        [
            original.gas_storage_change_j,
            original.acoustic_storage_change_j,
            original.supplied_heat_j,
            original.port_work_j,
            original.acoustic_dissipation_j,
            original.parameter_work_j,
            original.total_mass_residual_kg,
            original.combined_energy_residual_j
        ],
    );
    for (&id, step) in &report.acoustics {
        assert!(reactive.reaction(id).is_none());
        match (step, &original.acoustics[&id]) {
            (FiniteGasAcousticStep::Port(a), FiniteGasAcousticStep::Port(b)) => {
                assert_eq!(a.acoustic, b.acoustic);
                assert_eq!(a.refresh, b.refresh);
                assert_eq!(
                    a.background,
                    FiniteGasAcousticReactionSource::Caloric(b.background.clone())
                );
            }
            (FiniteGasAcousticStep::Aperture(a), FiniteGasAcousticStep::Aperture(b)) => {
                assert_eq!(a.coupled, b.coupled);
                assert_eq!(a.gas.acoustic, b.gas.acoustic);
                assert_eq!(a.gas.refresh, b.gas.refresh);
                assert_eq!(
                    a.gas.background,
                    FiniteGasAcousticReactionSource::Caloric(b.gas.background.clone())
                );
            }
            _ => panic!("native reduction changed owner"),
        }
    }
    assert_eq!(
        report.reservoirs[&90],
        FiniteGasAcousticReactionReservoir::Caloric(original.reservoirs[&90].clone())
    );
}

#[test]
fn g4_reacting_network_refusals_cancel_and_retry_all_original_histories() {
    let law = reaction(50_000.0);
    let reactions = bindings(&law);
    let (mut parcels, mut nodes, inputs) = fixture(&law);
    let (before_gas, before_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let preview = finite_gas_reacting_acoustic_transport_candidate_with_paths(
        &parcels,
        &nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, before_gas);
    equal_nodes(&nodes, &before_nodes);
    let (mut accepted_gas, mut accepted_nodes) = (parcels.clone(), nodes.clone());
    let mut polls = 0;
    let accepted = advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut accepted_gas,
        &mut accepted_nodes,
        &inputs,
        &reactions,
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
    equal_nodes(preview.nodes(), &accepted_nodes);
    let mut visited = 0;
    let cancelled = advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        &reactions,
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
        Err(ReactiveFiniteGasAcousticTransportError::Cancelled)
    ));
    assert_eq!(parcels, before_gas);
    equal_nodes(&nodes, &before_nodes);
    // ID20 refuses after the earlier complete load has already been proposed.
    let mut overheated = inputs.clone();
    thermal_mut(overheated.get_mut(&20).unwrap()).supplied_heat_j = 1000.0;
    assert!(
        advance_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &overheated,
            &reactions,
            paths,
            DT,
            reactive_solve(),
            &mut || false,
        )
        .is_err()
    );
    let mut foreign = reactions.clone();
    foreign.insert(
        999,
        FiniteGasAcousticReaction {
            reaction: &law,
            controls: chemistry(),
        },
    );
    assert!(matches!(
        advance_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &inputs,
            &foreign,
            paths,
            DT,
            reactive_solve(),
            &mut || false,
        ),
        Err(ReactiveFiniteGasAcousticTransportError::Invalid(_))
    ));
    let mut incompatible = reactions.clone();
    incompatible
        .get_mut(&20)
        .unwrap()
        .controls
        .gas
        .minimum_temperature_k = 250.0;
    assert!(
        advance_reacting_finite_gas_acoustic_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &inputs,
            &incompatible,
            paths,
            DT,
            reactive_solve(),
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(parcels, before_gas);
    equal_nodes(&nodes, &before_nodes);
    let retry = advance_reacting_finite_gas_acoustic_transport_with_paths(
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
    assert_eq!(parcels, accepted_gas);
    equal_nodes(&nodes, &accepted_nodes);
    for id in [10, 20, 90] {
        assert_eq!(retry.reaction(id), accepted.reaction(id));
    }
    assert_final_entropy_refusal();
}

fn assert_final_entropy_refusal() {
    let law = reaction(1e-12);
    let field = field_for(&law);
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
        DgBoundary::Port(PortId(0)),
        DgBoundary::Closed,
    )
    .unwrap();
    let mut load = DgPortLoad::new(duct, DT, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap();
    let volume = core::f64::consts::PI * RADIUS * RADIUS * LENGTH;
    let mut gas =
        ClosedGasParcel::new(law.catalog().clone(), point(300.0, 101_325.0), volume).unwrap();
    // Closed-end constant inflow accumulates compression; oscillatory kinetic
    // storage stays bounded. This makes injection's material work positive.
    for _ in 0..128 {
        load.step_with_finite_gas_heat_controlled(
            &mut gas,
            1e-7,
            0.0,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    }
    let reservoir = ClosedGasParcel::new(
        law.catalog()
            .with_composition(mixture([0.0, 0.0, 1.0]))
            .unwrap(),
        point(300.0, 180_000.0),
        6.0 * volume,
    )
    .unwrap();
    let mut parcels = BTreeMap::from([(10, gas), (90, reservoir)]);
    let mut nodes = BTreeMap::from([(10, FiniteGasAcousticNode::Port(load))]);
    let inputs = [10, 90]
        .map(|id| {
            (
                id,
                FiniteGasAcousticInput::Passive(FiniteGasDgTransportInput {
                    gas: energy_solve(),
                    supplied_heat_j: 0.0,
                    acoustic: (id == 10).then_some(controls()),
                    ports: if id == 10 {
                        vec![fs_duct::dg::PortInput {
                            id: PortId(0),
                            effort: 0.0,
                        }]
                    } else {
                        vec![]
                    },
                }),
            )
        })
        .into();
    let links = [nozzles()[0]];
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let (before_gas, before_nodes) = (parcels.clone(), nodes.clone());
    let (mut frozen_gas, mut frozen_nodes) = (parcels.clone(), nodes.clone());
    let frozen = advance_finite_gas_acoustic_transport_with_paths(
        &mut frozen_gas,
        &mut frozen_nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert!(accepted_work(&frozen.acoustics[&10]) > 1e-15);
    let reactions = BTreeMap::from([(
        10,
        FiniteGasAcousticReaction {
            reaction: &law,
            controls: GasReactionSolve {
                entropy_tolerance_j_k: 1e-18,
                ..chemistry()
            },
        },
    )]);
    let failure = advance_reacting_finite_gas_acoustic_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        &reactions,
        paths,
        DT,
        reactive_solve(),
        &mut || false,
    )
    .err()
    .expect("uncompensated fixed-volume work must fail final entropy admission");
    assert!(matches!(
        failure,
        ReactiveFiniteGasAcousticTransportError::Native(FiniteGasDgTransportError::Gas {
            parcel: 10,
            error: GasReactionError::EntropyDecrease { .. },
        })
    ));
    assert_eq!(parcels, before_gas);
    equal_nodes(&nodes, &before_nodes);
    // Refusal did not spend transport or advance the native clock.
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
    assert_eq!(parcels, frozen_gas);
    equal_nodes(&nodes, &frozen_nodes);
}
