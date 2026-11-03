//! Native owner replay and a failed scheduled material exchange on a live source.
use super::*;
use fs_couple::bernoulli_aperture::performance::finite_transport::{
    FiniteGasReservoir, FiniteGasTransportProgram, GasNozzleControl,
    FiniteGasEffusionLink, FiniteGasEffusionProgram, GasEffusionControl,
};
use fs_couple::bernoulli_aperture::dg::finite_gas::transport::{
    FiniteGasApertureDrive, FiniteGasApertureTransportInput,
    FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    scalar::advance_finite_gas_scalar_aperture_transport_with_paths,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
    species::effusion::{EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
    species::network::transport::effusion::EffusionNetworkLink,
    species::network::{
        SpeciesNetworkSolve, transport::{GasTransportPaths, GasTransportSolve},
    },
};
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use std::collections::BTreeMap;

#[test]
fn reservoir_playback_replays_original_owner_and_late_refusal_preserves_every_participant() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let source_gas = gas(&cold);
    let point = QueryPoint::new()
        .with_quantity("temperature", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static)), 300.0).unwrap()
        .with_quantity("pressure", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Pressure, ValueForm::Static)), 180_000.0).unwrap();
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let reservoir_model = from_packs(
        &[a.binding(), b.binding()], composition(0.1), None,
    ).unwrap();
    let reservoir = ClosedGasParcel::new(
        reservoir_model,
        point, 6.0 * source_gas.volume_m3(),
    ).unwrap();
    let transport = FiniteGasDgTransportSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.05, mass_tolerance_kg: 1e-16,
                energy_tolerance_j: 1e-11,
            },
            entropy_tolerance_j_k: 1e-12,
        },
        energy_tolerance_j: 5e-11,
    };
    let heat = program();
    let authored = FiniteGasTransportProgram {
        source_id: 10,
        source: heat.clone(),
        reservoirs: BTreeMap::from([(90, FiniteGasReservoir {
            gas: reservoir.clone(), caloric: energy_solve(),
        })]),
        nozzles: vec![SnapshotNozzleLink {
            id: 9, parcels: [90, 10],
            nozzle: GasNozzle { area_m2: 0.0, discharge_coefficient: 0.8 },
            nozzle_solve: energy_solve(),
        }],
        controls: vec![
            GasNozzleControl { sample: 7, id: 9,
                nozzle: GasNozzle { area_m2: 1e-8, discharge_coefficient: 0.8 } },
            // This finite geometry is a valid authored control, but its actual
            // withdrawal exceeds the original inventory limit at sample 19.
            GasNozzleControl { sample: 19, id: 9,
                nozzle: GasNozzle { area_m2: 1.0, discharge_coefficient: 0.8 } },
        ],
        transport,
    };
    let mut player = performer(&cold, AperturePerformanceConfig {
        max_controls: 16, max_compile_work: 256, ..config(1)
    }).with_finite_gas_transport(source_gas.clone(), authored.clone(),
        Arc::new(CancelGate::new())).unwrap();
    // Cold admission previews privately; it cannot consume a donor or a clock.
    let mut parcels = BTreeMap::from([(10, source_gas.clone()), (90, reservoir)]);
    assert_eq!(player.finite_gas_parcels(), Some(&parcels));
    assert_eq!(player.samples_rendered(), 0);
    assert!(player.applied_nozzle_controls().is_empty());
    assert!(player.applied_heat_controls().is_empty());
    let mut nodes = BTreeMap::from([(10, source(&cold))]);
    let mut nozzles = authored.nozzles.clone();
    for sample in 0..19 {
        for event in authored.controls.iter().filter(|event| event.sample == sample) {
            nozzles[0].nozzle = event.nozzle;
        }
        let power = heat.controls.iter().rev()
            .find(|event| event.sample <= sample).unwrap().power_w;
        let inputs = BTreeMap::from([
            (10, FiniteGasApertureTransportInput {
                thermal: FiniteGasDgTransportInput {
                    gas: energy_solve(), supplied_heat_j: power * DT,
                    acoustic: Some(controls()), ports: vec![],
                },
                source: Some(FiniteGasApertureDrive {
                    drive: TubeDrive { upstream_pressure_pa: 30.0, body_flow_m3_s: 0.0 },
                    forces_n: [0.0; 2], pressures_pa: [0.0; 2],
                }),
            }),
            (90, FiniteGasApertureTransportInput {
                thermal: FiniteGasDgTransportInput {
                    gas: energy_solve(), supplied_heat_j: 0.0,
                    acoustic: None, ports: vec![],
                }, source: None,
            }),
        ]);
        let direct = advance_finite_gas_scalar_aperture_transport_with_paths(
            &mut parcels, &mut nodes, &inputs,
            GasTransportPaths { nozzles: &nozzles, ..Default::default() },
            DT, transport, &mut || false,
        ).unwrap();
        let mut pressure = [f64::NAN];
        player.block(&mut pressure).unwrap();
        let actual = player.last_finite_gas_transport_step().unwrap();
        assert_eq!(pressure[0].to_bits(),
            direct.acoustics[&10].coupled.aperture.bore_pressure_pa.to_bits());
        assert_eq!(actual.transport, direct.transport);
        assert_eq!(actual.acoustics[&10].coupled, direct.acoustics[&10].coupled);
        assert_eq!(actual.acoustics[&10].gas.background, direct.acoustics[&10].gas.background);
        assert_eq!(actual.combined_energy_residual_j, direct.combined_energy_residual_j);
        assert!(actual.combined_energy_residual_j.abs() <= transport.energy_tolerance_j);
        assert_eq!(player.finite_gas_parcels(), Some(&parcels));
        assert_eq!(player.finite_gas(), parcels.get(&10));
        same_state(retained(&player), &nodes[&10]);
        assert_eq!(player.applied_heat_controls(),
            &heat.controls[..heat.controls.iter().take_while(|event| event.sample <= sample).count()]);
        assert_eq!(player.applied_nozzle_controls(),
            &authored.controls[..authored.controls.iter().take_while(|event| event.sample <= sample).count()]);
    }
    assert!(parcels[&10].mass_kg() > source_gas.mass_kg());
    assert_ne!(parcels[&10].state().state().sound_speed, source_gas.state().state().sound_speed);
    let previous = player.last_finite_gas_transport_step().unwrap();
    let previous_transport = previous.transport.clone();
    let previous_frame = previous.acoustics[&10].coupled.clone();
    let heat_cursor = player.applied_heat_controls().to_vec();
    let nozzle_cursor = player.applied_nozzle_controls().to_vec();
    let pressure_cursor = player.applied_controls().to_vec();
    let mut refused_output = [123.0];
    assert!(player.block(&mut refused_output).is_err());
    assert_eq!(refused_output, [123.0]);
    assert_eq!(player.samples_rendered(), 19);
    assert_eq!(player.finite_gas_parcels(), Some(&parcels));
    assert_eq!(player.applied_heat_controls(), heat_cursor);
    assert_eq!(player.applied_nozzle_controls(), nozzle_cursor);
    assert_eq!(player.applied_controls(), pressure_cursor);
    assert_eq!(player.held_gas_nozzles(), nozzles);
    same_state(retained(&player), &nodes[&10]);
    let last = player.last_finite_gas_transport_step().unwrap();
    assert_eq!(last.transport, previous_transport);
    assert_eq!(last.acoustics[&10].coupled, previous_frame);
    assert!(matches!(player.block(&mut refused_output), Err(RenderError::Poisoned)));
}

fn effusion_fixture(cold: &AxialGasField)
    -> (ClosedGasParcel, FiniteGasTransportProgram, FiniteGasEffusionProgram)
{
    let source = gas(cold);
    let point = QueryPoint::new()
        .with_quantity("temperature", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static)), 300.0).unwrap()
        .with_quantity("pressure", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Pressure, ValueForm::Static)), 40_000.0).unwrap();
    let reservoir = ClosedGasParcel::new(
        cold.model().with_composition(composition(0.1)).unwrap(),
        point, 6.0 * source.volume_m3(),
    ).unwrap();
    let mut heat = program();
    heat.controls.truncate(1);
    let transport = FiniteGasTransportProgram {
        source_id: 10, source: heat,
        reservoirs: BTreeMap::from([(90, FiniteGasReservoir {
            gas: reservoir, caloric: energy_solve(),
        })]),
        nozzles: vec![], controls: vec![],
        transport: FiniteGasDgTransportSolve {
            transport: GasTransportSolve {
                exchange: SpeciesNetworkSolve {
                    maximum_species_fraction: 0.05, mass_tolerance_kg: 1e-16,
                    energy_tolerance_j: 1e-11,
                },
                entropy_tolerance_j_k: 1e-12,
            },
            energy_tolerance_j: 5e-11,
        },
    };
    // The source starts as pure A; qualification must retain the full original
    // A/B catalog because the other reservoir can introduce B during playback.
    let transmissions = cold.model().thermodynamics().components().iter()
        .map(|component| EffusiveSpeciesTransmission::try_new(
            component.species().clone(), 1.0, "synthetic reciprocal thin aperture",
        ).unwrap()).collect();
    let interface = EffusiveInterface::try_new(cold.model().clone(), transmissions,
        EffusiveSupport { temperature_k: (200.0, 900.0), pressure_pa: (1000.0, 200_000.0) },
        EffusiveRegime::try_new(1e-9, 0.0, 1e-7,
            "synthetic all-mixture collisionless qualification; not measured gas validation").unwrap(),
    ).unwrap();
    assert_eq!(interface.transmissions().len(), 2);
    let effusion = FiniteGasEffusionProgram {
        links: vec![FiniteGasEffusionLink { id: 12, parcels: [10, 90], interface, area_m2: 0.0 }],
        controls: [(3, 1e-8), (11, 0.0), (14, 2e-8), (19, 1.0)].into_iter()
            .map(|(sample, area_m2)| GasEffusionControl { sample, id: 12, area_m2 }).collect(),
    };
    (source, transport, effusion)
}

#[test]
fn effusive_playback_replays_species_cooling_and_rejects_late_area_atomically() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let (source_gas, mut transport, effusion) = effusion_fixture(&cold);
    // Distinct continuum and molecular exposures read one original inventory.
    transport.nozzles.push(SnapshotNozzleLink {
        id: 5, parcels: [10, 90],
        nozzle: GasNozzle { area_m2: 1e-10, discharge_coefficient: 0.8 },
        nozzle_solve: energy_solve(),
    });
    let mut player = performer(&cold, AperturePerformanceConfig {
        max_controls: 16, max_compile_work: 256, ..config(1)
    }).with_finite_gas_transport(source_gas.clone(), transport.clone(),
        Arc::new(CancelGate::new())).unwrap()
        .with_finite_gas_effusion(effusion.clone()).unwrap();
    let mut parcels = BTreeMap::from([
        (10, source_gas.clone()), (90, transport.reservoirs[&90].gas.clone()),
    ]);
    let mut nodes = BTreeMap::from([(10, source(&cold))]);
    assert_eq!(player.finite_gas_parcels(), Some(&parcels));
    assert_eq!(player.finite_gas_effusion_program(), Some(&effusion));
    assert_eq!(player.held_effusion_areas_m2(), &[0.0]);
    assert!(player.applied_effusion_controls().is_empty());
    let inputs = BTreeMap::from([
        (10, FiniteGasApertureTransportInput {
            thermal: FiniteGasDgTransportInput {
                gas: energy_solve(), supplied_heat_j: 0.0,
                acoustic: Some(controls()), ports: vec![],
            },
            source: Some(FiniteGasApertureDrive {
                drive: TubeDrive { upstream_pressure_pa: 30.0, body_flow_m3_s: 0.0 },
                forces_n: [0.0; 2], pressures_pa: [0.0; 2],
            }),
        }),
        (90, FiniteGasApertureTransportInput {
            thermal: FiniteGasDgTransportInput {
                gas: energy_solve(), supplied_heat_j: 0.0,
                acoustic: None, ports: vec![],
            }, source: None,
        }),
    ]);
    let mut area_m2 = 0.0;
    for sample in 0..19 {
        if let Some(event) = effusion.controls.iter().find(|event| event.sample == sample) {
            area_m2 = event.area_m2;
        }
        let path = EffusionNetworkLink {
            id: 12, parcels: [10, 90], interface: &effusion.links[0].interface, area_m2,
        };
        let direct = advance_finite_gas_scalar_aperture_transport_with_paths(
            &mut parcels, &mut nodes, &inputs,
            GasTransportPaths { nozzles: &transport.nozzles, effusion: &[path], ..Default::default() },
            DT, transport.transport, &mut || false,
        ).unwrap();
        let mut output = [f64::NAN];
        player.block(&mut output).unwrap();
        let actual = player.last_finite_gas_transport_step().unwrap();
        assert_eq!(actual.transport, direct.transport);
        assert_eq!(actual.acoustics[&10].coupled, direct.acoustics[&10].coupled);
        assert_eq!(actual.combined_energy_residual_j, direct.combined_energy_residual_j);
        assert!(actual.combined_energy_residual_j.abs() <= transport.transport.energy_tolerance_j);
        assert_eq!(output[0].to_bits(), direct.acoustics[&10].coupled.aperture.bore_pressure_pa.to_bits());
        assert_eq!(player.finite_gas_parcels(), Some(&parcels));
        assert_eq!(player.held_effusion_areas_m2(), &[area_m2]);
        assert_eq!(player.applied_effusion_controls(),
            &effusion.controls[..effusion.controls.iter().take_while(|event| event.sample <= sample).count()]);
        assert_eq!(actual.transport.effusion[0].flux.species.is_empty(), area_m2 == 0.0);
        same_state(retained(&player), &nodes[&10]);
    }
    assert!(parcels[&10].mass_kg() < source_gas.mass_kg());
    assert!(parcels[&10].state().state().temperature < 300.0);
    assert!(parcels[&10].model().thermodynamics().mole_composition().fractions()[1] > 0.0);
    let previous = player.last_finite_gas_transport_step().unwrap();
    let previous_transport = previous.transport.clone();
    let previous_frame = previous.acoustics[&10].coupled.clone();
    let cursor = player.applied_effusion_controls().to_vec();
    let mut output = [123.0];
    assert!(player.block(&mut output).is_err());
    assert_eq!(output, [123.0]);
    assert_eq!(player.samples_rendered(), 19);
    assert_eq!(player.finite_gas_parcels(), Some(&parcels));
    assert_eq!(player.held_effusion_areas_m2(), &[area_m2]);
    assert_eq!(player.applied_effusion_controls(), cursor);
    let last = player.last_finite_gas_transport_step().unwrap();
    assert_eq!(last.transport, previous_transport);
    assert_eq!(last.acoustics[&10].coupled, previous_frame);
    same_state(retained(&player), &nodes[&10]);
    assert!(matches!(player.block(&mut output), Err(RenderError::Poisoned)));
}

#[test]
fn effusion_controls_share_budgets_and_require_original_preplayback_bindings() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let (source, transport, effusion) = effusion_fixture(&cold);
    let attach = |config, program| performer(&cold, config)
        .with_finite_gas_transport(source.clone(), transport.clone(), Arc::new(CancelGate::new()))
        .unwrap().with_finite_gas_effusion(program);
    assert!(attach(config(1), effusion.clone()).is_ok());
    let mut active = effusion.clone();
    active.controls[0].sample = 0;
    let active_player = attach(config(1), active).unwrap();
    assert_eq!(active_player.finite_gas(), Some(&source));
    assert_eq!(active_player.samples_rendered(), 0);
    assert!(active_player.applied_effusion_controls().is_empty());
    assert_eq!(active_player.held_effusion_areas_m2(), &[0.0]);
    for limited in [AperturePerformanceConfig { max_controls: 7, ..config(1) },
        AperturePerformanceConfig { max_compile_work: 61, ..config(1) }] {
        assert!(attach(limited, effusion.clone()).is_err());
    }
    for defect in 0..5 {
        let mut bad = effusion.clone();
        match defect {
            0 => bad.links[0].parcels[1] = 91,
            1 => bad.controls[0].id = 13,
            2 => bad.controls[0].sample = 48,
            3 => bad.controls[0].area_m2 = f64::NAN,
            _ => bad.controls[1] = bad.controls[0],
        }
        assert!(attach(config(1), bad).is_err());
    }
    assert!(performer(&cold, config(1)).with_finite_gas_effusion(effusion.clone()).is_err());
    let player = attach(config(1), effusion.clone()).unwrap();
    assert!(player.with_finite_gas_effusion(effusion.clone()).is_err());
    let mut cancelled = attach(config(1), effusion.clone()).unwrap();
    let before = cancelled.finite_gas_parcels().unwrap().clone();
    let before_source = retained(&cancelled).clone();
    cancelled.finite_gas_transport_cancel_gate().unwrap().request();
    let mut output = [123.0];
    assert!(matches!(cancelled.block(&mut output), Err(RenderError::Voice(
        fs_couple::acoustic_realize::AcousticRealizeError::Cancelled))));
    assert_eq!(output, [123.0]);
    assert_eq!(cancelled.finite_gas_parcels(), Some(&before));
    assert_eq!(cancelled.samples_rendered(), 0);
    assert!(cancelled.applied_effusion_controls().is_empty());
    assert_eq!(cancelled.held_effusion_areas_m2(), &[0.0]);
    same_state(retained(&cancelled), &before_source);
    let mut player = performer(&cold, config(1)).with_finite_gas_transport(
        source, transport, Arc::new(CancelGate::new())).unwrap();
    player.block(&mut [0.0]).unwrap();
    assert!(player.with_finite_gas_effusion(effusion).is_err());
}
