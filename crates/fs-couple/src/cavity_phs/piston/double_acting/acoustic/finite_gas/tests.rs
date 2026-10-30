//! Shared native momentum, conservative species transfer and paid gas/wave work.
//! Synthetic source laws check implementation and conservation, not measured gases.

use super::super::AcousticMotionFace;
use super::super::tests as acoustic_fixture;
use super::*;
use crate::cavity_phs::piston::double_acting::DoublePistonInput;
use crate::cavity_phs::piston::tests::{controls as single_controls, gas};
use crate::dg_port::finite_gas::transport::{FiniteGasDgTransportInput, FiniteGasDgTransportSolve};
use crate::dg_port::finite_gas::{FiniteGasDgLoss, FiniteGasDgSolve};
use crate::dg_port::{DgPortLoad, DgTerminalChannelInput, DgTerminalChannelKind};
use fs_duct::dg::{PortId, PortInput};
use fs_matdb::QueryPoint;
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
use std::collections::BTreeMap;

pub(super) type Parcels = BTreeMap<u64, ClosedGasParcel>;
pub(super) type Loads = BTreeMap<u64, DgPortLoad>;
pub(super) type Inputs = BTreeMap<u64, FiniteGasDgTransportInput>;

const BORE_VOLUME: f64 = core::f64::consts::PI * 0.001 * 0.001 * 0.05;
const AREAS: [f64; 2] = [0.001, -0.0007];

pub(super) fn finite_controls() -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        work_bracket_j: [-1e-5, 1e-5],
        work_tolerance_j: 1e-14,
        energy_tolerance_j: 2e-11,
        maximum_iterations: 80,
        maximum_temperature_increment_k: 80.0,
        loss: FiniteGasDgLoss::Lossless,
    }
}

pub(super) fn transport_controls() -> FiniteGasDgTransportSolve {
    FiniteGasDgTransportSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.05,
                mass_tolerance_kg: 1e-16,
                energy_tolerance_j: 1e-11,
            },
            entropy_tolerance_j_k: 1e-12,
        },
        // This gate encloses the original piston as well as regional gas work.
        energy_tolerance_j: acoustic_fixture::motion().energy_tolerance_j,
    }
}

pub(super) fn fixture(dt: f64) -> (Parcels, Loads, Inputs) {
    let mut parcels = BTreeMap::new();
    let mut loads = BTreeMap::new();
    for (id, fraction) in [(10, 0.15), (20, 0.85)] {
        let source = gas(fraction, 0.001);
        // The piston fixture's gas has the right sources, but its chamber volume
        // is not the represented acoustic inventory. Bind the actual bore here.
        let parcel =
            ClosedGasParcel::new(source.model().clone(), source.point().clone(), BORE_VOLUME)
                .unwrap();
        let load = acoustic_fixture::load(fraction, dt);
        assert_eq!(
            load.duct().gas_field().unwrap().uniform_state(),
            Some(parcel.state())
        );
        parcels.insert(id, parcel);
        loads.insert(id, load);
    }
    let reservoir_source = gas(0.1, 0.001);
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            400.0,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            180_000.0,
        )
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    parcels.insert(
        90,
        ClosedGasParcel::new(reservoir_source.model().clone(), point, 6.0 * BORE_VOLUME).unwrap(),
    );
    let inputs = parcels
        .keys()
        .map(|&id| {
            (
                id,
                FiniteGasDgTransportInput {
                    gas: single_controls().gas,
                    supplied_heat_j: match id {
                        10 => 2e-5,
                        20 => -1e-5,
                        _ => 3e-6,
                    },
                    acoustic: loads.contains_key(&id).then_some(finite_controls()),
                    // The shared mechanical owner supplies the only inlet flow.
                    ports: Vec::new(),
                },
            )
        })
        .collect();
    (parcels, loads, inputs)
}

pub(super) fn faces(loads: &Loads) -> [FiniteGasAcousticFace<'_>; 2] {
    [(10, AREAS[0]), (20, AREAS[1])].map(|(parcel, signed_area_m2)| FiniteGasAcousticFace {
        parcel,
        load: &loads[&parcel],
        signed_area_m2,
    })
}

pub(super) fn original_faces(loads: &Loads) -> [AcousticMotionFace<'_>; 2] {
    [(10, AREAS[0]), (20, AREAS[1])].map(|(parcel, signed_area_m2)| AcousticMotionFace {
        load: &loads[&parcel],
        signed_area_m2,
    })
}

pub(super) fn nozzles() -> [SnapshotNozzleLink; 2] {
    [(9, 10), (31, 20)].map(|(id, receiver)| SnapshotNozzleLink {
        id,
        parcels: [90, receiver],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: single_controls().gas,
    })
}

pub(super) fn assert_loads_same(actual: &Loads, expected: &Loads) {
    assert!(actual.keys().eq(expected.keys()));
    for (&id, load) in actual {
        acoustic_fixture::assert_load_same(load, &expected[&id]);
        assert_eq!(
            load.duct().operator_identity(),
            expected[&id].duct().operator_identity()
        );
        assert_eq!(
            load.cumulative_parameter_work_j(),
            expected[&id].cumulative_parameter_work_j()
        );
    }
}

fn species_mass(parcels: &Parcels, index: usize) -> f64 {
    // Independent conversion for the fixture's original A/B species masses.
    let molar = [0.028, 0.004];
    parcels
        .values()
        .map(|parcel| {
            let x = parcel
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions();
            parcel.mass_kg() * x[index] * molar[index] / (x[0] * molar[0] + x[1] * molar[1])
        })
        .sum()
}

#[test]
fn one_transport_snapshot_and_one_shared_motion_close_the_actual_enclosure() {
    let input = acoustic_fixture::forcing();
    let live = acoustic_fixture::pair(0.03);
    let (parcels, loads, inputs) = fixture(input.duration_s);
    let (before, old_parcels, old_loads) = (live.clone(), parcels.clone(), loads.clone());
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let endpoints = inputs.iter().map(|(&id, row)| (id, row.gas)).collect();
    let material = material_transport::step_candidate_with_paths(
        &parcels,
        paths,
        input.duration_s,
        &endpoints,
        transport_controls().transport,
        || false,
    )
    .unwrap();
    let frozen = live
        .acoustic_candidate(
            input,
            acoustic_fixture::mechanics(),
            &original_faces(&loads),
            acoustic_fixture::motion(),
            || false,
        )
        .unwrap();
    let candidate = live
        .transported_acoustic_candidate(
            &faces(&loads),
            &parcels,
            &inputs,
            paths,
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || false,
        )
        .unwrap();
    let report = candidate.step();
    acoustic_fixture::assert_piston_same(&candidate.native().0, frozen.piston());
    assert_eq!(&candidate.native().1, frozen.piston_step());
    assert_eq!(report.transport.transport, *material.step());
    assert_eq!(report.face_parcels, vec![10, 20]);
    assert_ne!(
        loads[&10].duct().inlet_gas().density,
        loads[&20].duct().inlet_gas().density
    );
    for row in &report.transport.transport.nozzles {
        let link = links.iter().find(|link| link.id == row.id).unwrap();
        let original = parcels[&90]
            .nozzle_flow_to(
                parcels[&row.receiver].state().state().pressure,
                link.nozzle,
                link.nozzle_solve,
                || false,
            )
            .unwrap();
        assert_eq!(row.flow, original);
        assert_eq!(
            row.transferred_mass_kg,
            original.mass_flow_kg_s * input.duration_s
        );
        assert!(row.transferred_mass_kg > 0.0);
    }
    for index in 0..2 {
        acoustic_fixture::near(
            species_mass(candidate.parcels(), index),
            species_mass(&parcels, index),
            transport_controls().transport.exchange.mass_tolerance_kg,
        );
    }
    let (mut storage, mut loss, mut exterior, mut parameter) = (0.0, 0.0, 0.0, 0.0);
    for (index, id) in [10, 20].into_iter().enumerate() {
        let next = &candidate.loads()[&id];
        let original = &loads[&id];
        let native = &frozen.loads()[index];
        let frame = &report.motion.faces[index];
        let initial_frame = &frozen.motion().faces[index];
        assert_eq!(next.state().step_count(), original.state().step_count() + 1);
        assert_eq!(
            next.terminal().unwrap().pole_states(),
            native.terminal().unwrap().pole_states()
        );
        assert_eq!(next.last_terminal_frames(), native.last_terminal_frames());
        assert_eq!(frame.terminals, initial_frame.terminals);
        assert_eq!(frame.inlet_flow_m3_s, initial_frame.inlet_flow_m3_s);
        assert_eq!(frame.mean_pressure_pa, initial_frame.mean_pressure_pa);
        assert_eq!(frame.inlet_work_j, initial_frame.inlet_work_j);
        assert_eq!(frame.duct, report.transport.acoustics[&id].acoustic);
        assert_eq!(
            next.duct().gas_field().unwrap().uniform_state(),
            Some(candidate.parcels()[&id].state())
        );
        assert_eq!(
            candidate.parcels()[&id].volume_m3(),
            parcels[&id].volume_m3()
        );
        assert_eq!(
            candidate.parcels()[&id].mass_kg(),
            material.parcels()[&id].mass_kg()
        );
        assert_ne!(
            candidate.parcels()[&id]
                .model()
                .thermodynamics()
                .mole_composition(),
            parcels[&id].model().thermodynamics().mole_composition()
        );
        storage += next.stored_energy_j() - original.stored_energy_j();
        loss += next.cumulative_dissipated_j() - original.cumulative_dissipated_j();
        exterior += next
            .last_terminal_frames()
            .iter()
            .map(|row| row.load_work_j)
            .sum::<f64>();
        parameter += report.transport.acoustics[&id].refresh.parameter_work_j;
    }
    acoustic_fixture::near(report.motion.acoustic_storage_change_j, storage, 1e-12);
    acoustic_fixture::near(report.motion.acoustic_dissipated_j, loss, 1e-12);
    acoustic_fixture::near(report.motion.parameter_work_j, parameter, 1e-12);
    assert!(parameter.abs() > 1e-12);

    // Use actual endpoint inventories and the independently integrated fixture
    // caloric laws, rather than recombining the reported local residuals.
    let next = &candidate.native().0;
    let gas_change = acoustic_fixture::caloric_energy(next.primary().gas())
        - acoustic_fixture::caloric_energy(live.primary().gas())
        + acoustic_fixture::caloric_energy(next.opposing())
        - acoustic_fixture::caloric_energy(live.opposing())
        + candidate
            .parcels()
            .iter()
            .map(|(id, parcel)| {
                acoustic_fixture::caloric_energy(parcel)
                    - acoustic_fixture::caloric_energy(&parcels[id])
            })
            .sum::<f64>();
    let kinetic = 0.5
        * live.primary().mass_kg()
        * (next.primary().velocity_m_s().powi(2) - live.primary().velocity_m_s().powi(2));
    let supplied = input.heat_j.iter().sum::<f64>()
        + inputs.values().map(|row| row.supplied_heat_j).sum::<f64>()
        + input.external_force_n * candidate.native().1.displacement_m;
    let residual = gas_change + kinetic + storage + loss + exterior - supplied;
    acoustic_fixture::near(residual, 0.0, acoustic_fixture::motion().energy_tolerance_j);
    acoustic_fixture::near(candidate.total_energy_residual_j(), residual, 5e-10);
    acoustic_fixture::near(
        report.motion.acoustic_energy_residual_j,
        storage + loss + exterior
            - report
                .motion
                .faces
                .iter()
                .map(|row| row.inlet_work_j)
                .sum::<f64>()
            - parameter,
        1e-12,
    );
    acoustic_fixture::assert_piston_same(&live, &before);
    assert_eq!(parcels, old_parcels);
    assert_loads_same(&loads, &old_loads);
}

#[test]
fn heated_regions_change_the_next_pressure_reaction_without_replaying_the_current_motion() {
    let input = acoustic_fixture::forcing();
    let live = acoustic_fixture::pair(0.03);
    let (parcels, loads, mut cold_inputs) = fixture(input.duration_s);
    for row in cold_inputs.values_mut() {
        row.supplied_heat_j = 0.0;
    }
    let mut hot_inputs = cold_inputs.clone();
    hot_inputs.get_mut(&10).unwrap().supplied_heat_j = 0.002;
    hot_inputs.get_mut(&20).unwrap().supplied_heat_j = 0.001;
    let run = |inputs: &Inputs| {
        live.transported_acoustic_candidate(
            &faces(&loads),
            &parcels,
            inputs,
            GasTransportPaths::default(),
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || false,
        )
        .unwrap()
    };
    let cold = run(&cold_inputs);
    let hot = run(&hot_inputs);
    acoustic_fixture::assert_piston_same(&hot.native().0, &cold.native().0);
    assert_eq!(hot.native().1, cold.native().1);
    let frozen = live
        .acoustic_candidate(
            input,
            acoustic_fixture::mechanics(),
            &original_faces(&loads),
            acoustic_fixture::motion(),
            || false,
        )
        .unwrap();
    acoustic_fixture::assert_piston_same(&cold.native().0, frozen.piston());
    for (index, id) in [10, 20].into_iter().enumerate() {
        let a = &hot.step().motion.faces[index];
        let b = &cold.step().motion.faces[index];
        assert_eq!(a.inlet_flow_m3_s, b.inlet_flow_m3_s);
        assert_eq!(a.mean_pressure_pa, b.mean_pressure_pa);
        assert_eq!(a.inlet_work_j, b.inlet_work_j);
        assert_eq!(a.terminals, b.terminals);
        assert_eq!(hot.parcels()[&id].mass_kg(), parcels[&id].mass_kg());
        assert!(
            hot.parcels()[&id].state().state().temperature
                > cold.parcels()[&id].state().state().temperature + 5.0
        );
        assert_eq!(
            cold.loads()[&id].terminal().unwrap().pole_states(),
            frozen.loads()[index].terminal().unwrap().pole_states()
        );
        // An identical-model caloric inversion preserves the state to its
        // admitted physical tolerance; it does not promise identical bits.
        acoustic_fixture::near(
            cold.loads()[&id].stored_energy_j(),
            frozen.loads()[index].stored_energy_j(),
            finite_controls().energy_tolerance_j,
        );
    }
    let next_input = DoublePistonInput {
        heat_j: [0.0; 2],
        external_force_n: 0.0,
        ..input
    };
    let next_cold = cold
        .native()
        .0
        .acoustic_candidate(
            next_input,
            acoustic_fixture::mechanics(),
            &original_faces(cold.loads()),
            acoustic_fixture::motion(),
            || false,
        )
        .unwrap();
    let next_hot = hot
        .native()
        .0
        .acoustic_candidate(
            next_input,
            acoustic_fixture::mechanics(),
            &original_faces(hot.loads()),
            acoustic_fixture::motion(),
            || false,
        )
        .unwrap();
    assert_ne!(next_hot.motion().force_n, next_cold.motion().force_n);
    let velocity_resolution =
        acoustic_fixture::mechanics().piston.momentum_tolerance_n_s / live.primary().mass_kg();
    assert!(
        (next_hot.piston().primary().velocity_m_s() - next_cold.piston().primary().velocity_m_s())
            .abs()
            > velocity_resolution
    );
    assert_eq!(next_hot.loads()[0].state().step_count(), 2);
    assert_eq!(next_cold.loads()[0].state().step_count(), 2);
}

#[test]
fn finite_terminal_gas_volume_and_native_channel_momentum_are_kept_in_the_same_account() {
    let input = acoustic_fixture::forcing();
    let live = acoustic_fixture::pair(0.03);
    let (mut parcels, mut loads, inputs) = fixture(input.duration_s);
    let radius = 0.0007;
    let length = 0.004;
    let cells = 8;
    let channel_volume = core::f64::consts::PI * radius * radius * length;
    let original = &parcels[&10];
    let gas = ClosedGasParcel::new(
        original.model().clone(),
        original.point().clone(),
        BORE_VOLUME + channel_volume,
    )
    .unwrap();
    let load = loads[&10]
        .clone()
        .with_terminal_channels(
            vec![DgTerminalChannelInput {
                port: PortId(1),
                kind: DgTerminalChannelKind::Circular {
                    radius_m: radius,
                    length_m: length,
                    radial_cells: cells,
                    maximum_mach: 0.3,
                },
            }],
            cells,
        )
        .unwrap();
    parcels.insert(10, gas);
    loads.insert(10, load);
    let frozen = live
        .acoustic_candidate(
            input,
            acoustic_fixture::mechanics(),
            &original_faces(&loads),
            acoustic_fixture::motion(),
            || false,
        )
        .unwrap();
    let links = nozzles();
    let candidate = live
        .transported_acoustic_candidate(
            &faces(&loads),
            &parcels,
            &inputs,
            GasTransportPaths {
                nozzles: &links,
                ..Default::default()
            },
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || false,
        )
        .unwrap();
    let old = frozen.loads()[0].terminal_channels()[0].channel().unwrap();
    let next = candidate.loads()[&10].terminal_channels()[0]
        .channel()
        .unwrap();
    let (rho0, rho1) = (old.gas().density, next.gas().density);
    assert_ne!(rho0, rho1);
    let (mut kinetic0, mut kinetic1) = (0.0, 0.0);
    for (i, (&v0, &v1)) in old
        .velocities_m_s()
        .iter()
        .zip(next.velocities_m_s())
        .enumerate()
    {
        let inner = radius * i as f64 / cells as f64;
        let outer = radius * (i + 1) as f64 / cells as f64;
        let volume = core::f64::consts::PI * length * (outer * outer - inner * inner);
        acoustic_fixture::near(rho1 * volume * v1, rho0 * volume * v0, 1e-20);
        kinetic0 += 0.5 * rho0 * volume * v0 * v0;
        kinetic1 += 0.5 * rho1 * volume * v1 * v1;
    }
    assert!(kinetic0 > 0.0);
    acoustic_fixture::near(old.stored_energy_j(), kinetic0, 1e-15);
    acoustic_fixture::near(next.stored_energy_j(), kinetic1, 1e-15);
    let work = candidate.step().transport.acoustics[&10]
        .refresh
        .terminal_channels[0]
        .transfer
        .unwrap()
        .parameter_work_j;
    acoustic_fixture::near(work, kinetic1 - kinetic0, 1e-15);
    assert_eq!(
        candidate.loads()[&10].last_terminal_channel_frames(),
        frozen.loads()[0].last_terminal_channel_frames()
    );
    assert_eq!(candidate.loads()[&10].state().step_count(), 1);
    assert_eq!(
        candidate.loads()[&10].terminal_channels()[0].gas(),
        candidate.parcels()[&10].state().state()
    );
    acoustic_fixture::near(
        candidate.total_energy_residual_j(),
        0.0,
        acoustic_fixture::motion().energy_tolerance_j,
    );
}

#[test]
fn invalid_bindings_and_late_cancellation_leave_every_owner_replayable() {
    let input = acoustic_fixture::forcing();
    let live = acoustic_fixture::pair(0.03);
    let (parcels, loads, inputs) = fixture(input.duration_s);
    let (before, before_parcels, before_loads) = (live.clone(), parcels.clone(), loads.clone());
    let links = nozzles();
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let run = |faces: &[FiniteGasAcousticFace<'_>], gas: &Parcels, rows: &Inputs| {
        live.transported_acoustic_candidate(
            faces,
            gas,
            rows,
            paths,
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || false,
        )
    };
    let mut supplied_ports = inputs.clone();
    supplied_ports.get_mut(&10).unwrap().ports.push(PortInput {
        id: loads[&10].inlet(),
        effort: 0.0,
    });
    assert!(run(&faces(&loads), &parcels, &supplied_ports).is_err());
    let mut incomplete = inputs.clone();
    incomplete.remove(&90);
    assert!(run(&faces(&loads), &parcels, &incomplete).is_err());
    let duplicate = [
        FiniteGasAcousticFace {
            parcel: 10,
            load: &loads[&10],
            signed_area_m2: AREAS[0],
        },
        FiniteGasAcousticFace {
            parcel: 10,
            load: &loads[&20],
            signed_area_m2: AREAS[1],
        },
    ];
    assert!(run(&duplicate, &parcels, &inputs).is_err());
    let mut wrong_volume = parcels.clone();
    let source = &parcels[&10];
    wrong_volume.insert(
        10,
        ClosedGasParcel::new(
            source.model().clone(),
            source.point().clone(),
            source.volume_m3() * 1.01,
        )
        .unwrap(),
    );
    assert!(run(&faces(&loads), &wrong_volume, &inputs).is_err());
    let mut exhausted = inputs.clone();
    exhausted.get_mut(&20).unwrap().acoustic = Some(FiniteGasDgSolve {
        maximum_temperature_increment_k: 1e-12,
        ..finite_controls()
    });
    assert!(run(&faces(&loads), &parcels, &exhausted).is_err());

    let mut polls = 0usize;
    let expected = live
        .transported_acoustic_candidate(
            &faces(&loads),
            &parcels,
            &inputs,
            paths,
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);
    let final_poll = polls;
    polls = 0;
    let refused = live.transported_acoustic_candidate(
        &faces(&loads),
        &parcels,
        &inputs,
        paths,
        input,
        acoustic_fixture::mechanics(),
        acoustic_fixture::motion(),
        transport_controls(),
        &mut || {
            polls += 1;
            polls == final_poll
        },
    );
    assert!(matches!(refused, Err(FiniteGasAcousticError::Cancelled)));
    assert_eq!(polls, final_poll);
    acoustic_fixture::assert_piston_same(&live, &before);
    assert_eq!(parcels, before_parcels);
    assert_loads_same(&loads, &before_loads);

    let retried = run(&faces(&loads), &parcels, &inputs).unwrap();
    acoustic_fixture::assert_piston_same(&retried.native().0, &expected.native().0);
    assert_eq!(retried.native().1, expected.native().1);
    assert_eq!(retried.parcels(), expected.parcels());
    assert_loads_same(retried.loads(), expected.loads());
    assert_eq!(retried.step().motion, expected.step().motion);
    assert_eq!(
        retried.step().transport.transport,
        expected.step().transport.transport
    );
    let (a, gas_a, loads_a, _) = retried.into_parts();
    let (b, gas_b, loads_b, _) = expected.into_parts();
    let next_a =
        a.0.transported_acoustic_candidate(
            &faces(&loads_a),
            &gas_a,
            &inputs,
            paths,
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || false,
        )
        .unwrap();
    let next_b =
        b.0.transported_acoustic_candidate(
            &faces(&loads_b),
            &gas_b,
            &inputs,
            paths,
            input,
            acoustic_fixture::mechanics(),
            acoustic_fixture::motion(),
            transport_controls(),
            &mut || false,
        )
        .unwrap();
    acoustic_fixture::assert_piston_same(&next_a.native().0, &next_b.native().0);
    assert_eq!(next_a.native().1, next_b.native().1);
    assert_eq!(next_a.parcels(), next_b.parcels());
    assert_loads_same(next_a.loads(), next_b.loads());
}
