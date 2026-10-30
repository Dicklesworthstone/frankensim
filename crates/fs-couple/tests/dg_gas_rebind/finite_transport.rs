//! Species transport, moving mechanics and original acoustic histories share one commit.
//! Synthetic source laws independently check implementation, not measured instruments.
use super::finite_gas::{channel_load, controls, energy_solve, gas};
use super::*;
use fs_couple::bernoulli_aperture::dg::finite_gas::transport::{
    FiniteGasApertureDrive, FiniteGasApertureTransportInput,
    advance_finite_gas_aperture_transport_with_paths,
    finite_gas_aperture_transport_candidate_with_paths,
};
use fs_couple::dg_port::finite_gas::transport::{
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
    species::network::{
        SpeciesNetworkSolve,
        transport::{self as material_transport, GasTransportPaths, GasTransportSolve},
    },
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::collections::BTreeMap;

type Parcels = BTreeMap<u64, ClosedGasParcel>;
type Nodes = BTreeMap<u64, LoadedApertureDg>;
type Inputs = BTreeMap<u64, FiniteGasApertureTransportInput>;

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure,
        )
        .unwrap()
}

fn network() -> (Parcels, Nodes, Inputs) {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let first = channel_load(&cold);
    let single =
        LoadedApertureDg::new_loaded(aperture(&first, first.impedance_pa_s_m3()), first).unwrap();
    let second = channel_load(&cold);
    let paired =
        LoadedApertureDg::new_loaded(paired_source(&second, second.impedance_pa_s_m3()), second)
            .unwrap();
    let mut nodes = BTreeMap::from([(10, single), (20, paired)]);
    // Clone/preview tests begin with actual retained motion, channel momentum
    // and exterior poles, not an all-zero cold state.
    for node in nodes.values_mut() {
        for sample in 0..24 {
            node.step(drive(sample)).unwrap();
        }
    }
    let first = gas(&cold);
    let reservoir = ClosedGasParcel::new(
        first.model().with_composition(composition(0.1)).unwrap(),
        point(300.0, 180_000.0),
        6.0 * first.volume_m3(),
    )
    .unwrap();
    let parcels = BTreeMap::from([(10, first.clone()), (20, first), (90, reservoir)]);
    let inputs = parcels
        .keys()
        .map(|&id| {
            (
                id,
                FiniteGasApertureTransportInput {
                    thermal: FiniteGasDgTransportInput {
                        gas: energy_solve(),
                        supplied_heat_j: if id == 20 { -1e-5 } else { 2e-5 },
                        acoustic: nodes.contains_key(&id).then_some(controls()),
                        ports: Vec::new(),
                    },
                    source: nodes.contains_key(&id).then_some(FiniteGasApertureDrive {
                        drive: TubeDrive {
                            upstream_pressure_pa: if id == 10 { 30.0 } else { 20.0 },
                            body_flow_m3_s: if id == 10 { 2e-9 } else { -1e-9 },
                        },
                        forces_n: if id == 10 { [2e-5, 0.0] } else { [2e-5, -3e-5] },
                        pressures_pa: [0.0; 2],
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
        nozzle_solve: energy_solve(),
    })
}

// Independent ideal-mixture masses and NASA integrals for the declared fixture:
// species A/B have molar masses 0.028/0.044 and Cp/R = 3.5/4.5 + 0.001*T.
fn species_masses(gas: &ClosedGasParcel) -> [f64; 2] {
    let fractions = gas.model().thermodynamics().mole_composition().fractions();
    let mean_mass = fractions[0] * 0.028 + fractions[1] * 0.044;
    [
        gas.mass_kg() * fractions[0] * 0.028 / mean_mass,
        gas.mass_kg() * fractions[1] * 0.044 / mean_mass,
    ]
}
fn gas_energy(gas: &ClosedGasParcel) -> f64 {
    let mass = species_masses(gas);
    let t = gas.state().state().temperature;
    R * (mass[0] / 0.028 * (2.5 * t + 0.0005 * t * t)
        + mass[1] / 0.044 * (3.5 * t + 0.0005 * t * t))
}
fn mechanical_energy(node: &LoadedApertureDg) -> f64 {
    match node.aperture() {
        LoadedAperture::Single(source) => {
            let state = source.state();
            let spec = source.spec();
            0.5 * spec.mass_kg * state.opening_velocity_m_s.powi(2)
                + 0.5
                    * spec.stiffness_n_m
                    * (state.opening_m - spec.aperture.rest_opening_m).powi(2)
        }
        LoadedAperture::Paired(source) => source
            .spec()
            .blades
            .iter()
            .zip(source.blade_states())
            .map(|(blade, state)| {
                0.5 * blade.mass_kg * state.velocity_m_s.powi(2)
                    + 0.5 * blade.stiffness_n_m * state.displacement_m.powi(2)
            })
            .sum(),
        _ => unreachable!("fixture has exactly the two original mechanical layouts"),
    }
}
fn complete_energy(node: &LoadedApertureDg) -> f64 {
    let gas = node.load().duct().inlet_gas();
    let channel = node.load().terminal_channels()[0].channel().unwrap();
    let channel_volume = 0.002 * 4e-5 * 0.001;
    let kinetic = channel
        .velocities_m_s()
        .iter()
        .map(|v| 0.5 * gas.density * channel_volume / 8.0 * v * v)
        .sum::<f64>();
    mechanical_energy(node) + wave_energy(node.duct_state(), gas, LENGTH) + kinetic
}

fn unchanged(nodes: &Nodes, before: &Nodes) {
    assert!(nodes.keys().eq(before.keys()));
    for (&id, node) in nodes {
        let old = &before[&id];
        assert_eq!(node.duct_state(), old.duct_state());
        assert_eq!(
            node.load().duct().gas_field(),
            old.load().duct().gas_field()
        );
        assert_eq!(
            node.aperture().accepted_steps(),
            old.aperture().accepted_steps()
        );
        assert_eq!(node.stored_energy_j(), old.stored_energy_j());
        assert_eq!(
            node.load().terminal().unwrap().pole_states(),
            old.load().terminal().unwrap().pole_states()
        );
        assert_eq!(
            node.load().last_terminal_frames(),
            old.load().last_terminal_frames()
        );
        assert_eq!(
            node.load().last_terminal_channel_frames(),
            old.load().last_terminal_channel_frames()
        );
        assert_eq!(
            node.load().terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s(),
            old.load().terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s()
        );
        match (node.aperture(), old.aperture()) {
            (LoadedAperture::Single(a), LoadedAperture::Single(b)) => {
                assert_eq!(a.state(), b.state());
                assert_eq!(a.spec(), b.spec());
                assert_eq!(
                    a.viscous_slit().unwrap().gas(),
                    b.viscous_slit().unwrap().gas()
                );
            }
            (LoadedAperture::Paired(a), LoadedAperture::Paired(b)) => {
                assert_eq!(a.blade_states(), b.blade_states());
                assert_eq!(a.spec(), b.spec());
                assert_eq!(a.gas(), b.gas());
            }
            _ => panic!("mechanical ownership changed"),
        }
    }
}

#[test]
fn g1_transport_changes_the_same_moving_sources_and_closes_independent_energy() {
    let (mut parcels, mut nodes, inputs) = network();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let endpoints = inputs
        .iter()
        .map(|(&id, input)| (id, input.thermal.gas))
        .collect();
    let original_transport = material_transport::step_candidate_with_paths(
        &parcels,
        paths,
        DT,
        &endpoints,
        solve().transport,
        || false,
    )
    .unwrap();
    let report = advance_finite_gas_aperture_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        paths,
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport, *original_transport.step());
    for species in 0..2 {
        let before: f64 = old_gas.values().map(|p| species_masses(p)[species]).sum();
        let after: f64 = parcels.values().map(|p| species_masses(p)[species]).sum();
        assert!((after - before).abs() < 1e-16);
    }
    let (mut work, mut loss, mut storage_change) = (0.0, 0.0, 0.0);
    let mut frozen = old_nodes.clone();
    for (&id, node) in &nodes {
        let source = inputs[&id].source.as_ref().unwrap();
        let frame = frozen
            .get_mut(&id)
            .unwrap()
            .step_with_forces_and_pressures(source.drive, source.forces_n, source.pressures_pa)
            .unwrap();
        let step = &report.acoustics[&id];
        assert_eq!(step.coupled.aperture, frame.aperture);
        assert_eq!(
            node.load().terminal().unwrap().pole_states(),
            frozen[&id].load().terminal().unwrap().pole_states()
        );
        assert_eq!(
            node.load().last_terminal_frames(),
            frozen[&id].load().last_terminal_frames()
        );
        assert_eq!(node.aperture().accepted_steps(), 25);
        assert_eq!(node.duct_state().step_count(), 25);
        assert_eq!(
            node.load().duct().gas_field().unwrap().uniform_state(),
            Some(parcels[&id].state())
        );
        assert_eq!(
            node.load().terminal_channels()[0].gas(),
            parcels[&id].state().state()
        );
        assert_ne!(species_masses(&parcels[&id])[1], 0.0);
        assert!(frame.body_work_j != 0.0);
        assert!(frame.aperture.mechanical_work_j != 0.0);
        assert!(frame.terminals[0].channel_work_j != 0.0);
        assert!(frame.terminals[0].load_work_j != 0.0);
        // Only the actual exterior exchange leaves the enclosure. Channel
        // pressure work already contains represented finite-channel storage.
        let exterior_load: f64 = frame
            .terminals
            .iter()
            .map(|terminal| terminal.load_work_j)
            .sum();
        work += frame.upstream_work_j + frame.body_work_j + frame.aperture.mechanical_work_j
            - exterior_load;
        loss += frame.inferred_total_dissipation_j;
        close(
            complete_energy(node),
            node.stored_energy_j(),
            node.stored_energy_j(),
        );
        storage_change += complete_energy(node) - complete_energy(&old_nodes[&id]);
        let density = parcels[&id].state().state().density;
        match node.aperture() {
            LoadedAperture::Single(a) => {
                assert_eq!(a.spec().density_kg_m3, density);
                assert_eq!(a.spec().impedance_pa_s_m3, node.load().impedance_pa_s_m3());
                assert_eq!(
                    a.viscous_slit().unwrap().gas(),
                    Some(parcels[&id].state().state())
                );
            }
            LoadedAperture::Paired(a) => {
                assert_eq!(a.gas(), Some(parcels[&id].state().state()));
                assert_ne!(a.blade_states()[0], a.blade_states()[1]);
            }
            _ => unreachable!(),
        }
    }
    let caloric_change: f64 = parcels
        .iter()
        .map(|(id, p)| gas_energy(p) - gas_energy(&old_gas[id]))
        .sum();
    let heat: f64 = inputs.values().map(|i| i.thermal.supplied_heat_j).sum();
    let residual = caloric_change + storage_change + loss - heat - work;
    assert!(residual.abs() <= solve().energy_tolerance_j, "{residual:e}");
    close(report.port_work_j, work, work.abs());
    close(
        report.acoustic_storage_change_j,
        storage_change,
        storage_change.abs(),
    );
    assert!(report.combined_energy_residual_j.abs() <= solve().energy_tolerance_j);
    assert!(report.total_mass_residual_kg.abs() <= solve().transport.exchange.mass_tolerance_kg);
    // The next step actually feels the accepted gas; transport is not a side report.
    let mut changed = 0.0;
    for (&id, node) in &mut nodes {
        let source = inputs[&id].source.as_ref().unwrap();
        let next = node
            .step_with_forces_and_pressures(source.drive, source.forces_n, source.pressures_pa)
            .unwrap();
        let reference = frozen
            .get_mut(&id)
            .unwrap()
            .step_with_forces_and_pressures(source.drive, source.forces_n, source.pressures_pa)
            .unwrap();
        changed += (next.aperture.bore_pressure_pa - reference.aperture.bore_pressure_pa).powi(2);
    }
    assert!(changed > 0.0);
}

#[test]
fn g3_no_transport_reduces_to_the_original_finite_gas_mechanics_and_heat() {
    let (mut parcels, mut nodes, inputs) = network();
    let (mut expected_gas, mut expected_nodes) = (parcels.clone(), nodes.clone());
    let mut expected = BTreeMap::new();
    for (&id, node) in &mut expected_nodes {
        let input = &inputs[&id];
        let source = input.source.as_ref().unwrap();
        let step = node
            .step_with_finite_gas_heat_controlled(
                expected_gas.get_mut(&id).unwrap(),
                source.drive,
                source.forces_n,
                source.pressures_pa,
                input.thermal.supplied_heat_j,
                input.thermal.gas,
                input.thermal.acoustic.unwrap(),
                &mut || false,
            )
            .unwrap();
        expected.insert(id, step);
    }
    expected_gas
        .get_mut(&90)
        .unwrap()
        .add_heat(inputs[&90].thermal.supplied_heat_j, energy_solve(), || {
            false
        })
        .unwrap();
    let report = advance_finite_gas_aperture_transport_with_paths(
        &mut parcels,
        &mut nodes,
        &inputs,
        GasTransportPaths::default(),
        DT,
        solve(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, expected_gas);
    unchanged(&nodes, &expected_nodes);
    for (&id, step) in &report.acoustics {
        assert_eq!(step.coupled, expected[&id].coupled);
        assert_eq!(step.gas.refresh, expected[&id].gas.refresh);
        assert_eq!(step.gas.work_by_gas_j, expected[&id].gas.work_by_gas_j);
    }
}

#[test]
fn g4_final_cancellation_and_late_source_refusal_preserve_all_owners_and_retry() {
    let (mut parcels, mut nodes, inputs) = network();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let preview = finite_gas_aperture_transport_candidate_with_paths(
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
    unchanged(&nodes, &old_nodes);
    let mut bad = inputs.clone();
    bad.get_mut(&20).unwrap().thermal.supplied_heat_j = 1000.0;
    assert!(matches!(
        advance_finite_gas_aperture_transport_with_paths(
            &mut parcels,
            &mut nodes,
            &bad,
            paths,
            DT,
            solve(),
            &mut || false
        ),
        Err(FiniteGasDgTransportError::Acoustic { parcel: 20, .. })
    ));
    assert_eq!(parcels, old_gas);
    unchanged(&nodes, &old_nodes);
    let (mut expected_gas, mut expected_nodes) = (parcels.clone(), nodes.clone());
    let mut polls = 0;
    let expected = advance_finite_gas_aperture_transport_with_paths(
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
    unchanged(preview.nodes(), &expected_nodes);
    let mut seen = 0;
    assert!(matches!(
        advance_finite_gas_aperture_transport_with_paths(
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
    unchanged(&nodes, &old_nodes);
    let retry = advance_finite_gas_aperture_transport_with_paths(
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
    unchanged(&nodes, &expected_nodes);
    assert_eq!(retry.transport, expected.transport);
    for (&id, step) in &retry.acoustics {
        assert_eq!(step.coupled, expected.acoustics[&id].coupled);
        assert_eq!(step.gas.refresh, expected.acoustics[&id].gas.refresh);
    }
    for sample in 25..32 {
        for (&id, node) in &mut nodes {
            assert_eq!(
                node.step(drive(sample)).unwrap(),
                expected_nodes
                    .get_mut(&id)
                    .unwrap()
                    .step(drive(sample))
                    .unwrap()
            );
        }
    }
}

#[test]
fn g4_transport_refuses_ambiguous_source_ports_reservoir_drives_and_clock() {
    let (mut parcels, mut nodes, inputs) = network();
    let (old_gas, old_nodes) = (parcels.clone(), nodes.clone());
    for case in 0..5 {
        let mut bad = inputs.clone();
        match case {
            0 => bad.get_mut(&10).unwrap().source = None,
            1 => bad.get_mut(&90).unwrap().source = inputs[&10].source.clone(),
            2 => {
                bad.get_mut(&10).unwrap().thermal.ports = vec![fs_duct::dg::PortInput {
                    id: PortId(0),
                    effort: 1e-8,
                }]
            }
            3 => bad.get_mut(&20).unwrap().source.as_mut().unwrap().forces_n[0] = f64::NAN,
            _ => {}
        }
        let duration = if case == 4 { DT * 2.0 } else { DT };
        assert!(
            advance_finite_gas_aperture_transport_with_paths(
                &mut parcels,
                &mut nodes,
                &bad,
                GasTransportPaths::default(),
                duration,
                solve(),
                &mut || false
            )
            .is_err()
        );
        assert_eq!(parcels, old_gas);
        unchanged(&nodes, &old_nodes);
    }
}
