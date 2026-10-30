//! Native owner replay and a failed scheduled material exchange on a live source.
use super::*;
use fs_couple::bernoulli_aperture::performance::finite_transport::{
    FiniteGasReservoir, FiniteGasTransportProgram, GasNozzleControl,
};
use fs_couple::bernoulli_aperture::dg::finite_gas::transport::{
    FiniteGasApertureDrive, FiniteGasApertureTransportInput,
    FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    scalar::advance_finite_gas_scalar_aperture_transport_with_paths,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
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
