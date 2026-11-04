//! Original primitive gas and complete load memory on the shared exact clock.
use super::*;
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_thermochem::{Composition, CompositionBasis};

#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

#[path = "gas_field/molecular.rs"]
mod molecular;

fn field(length: f64, uniform: bool) -> AxialGasField {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let primitive = |t, f| PrimitiveGasState {
        temperature_k: t,
        composition: Composition::new(CompositionBasis::MoleFraction, vec![f, 1.0 - f]).unwrap(),
    };
    let inlet = primitive(300.0, 0.8);
    let outlet = if uniform {
        inlet.clone()
    } else {
        primitive(600.0, 0.2)
    };
    let model = from_packs(&[a.binding(), b.binding()], inlet.composition.clone(), None).unwrap();
    AxialGasField::new(
        vec![0.0, length],
        101325.0,
        vec![inlet, outlet],
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

#[test]
fn uniform_source_field_preserves_complete_time_history_and_face_work() {
    let profile = RadiusProfile::exponential(0.025, 4.0, 0.04).unwrap();
    let original = field(profile.length(), true);
    let gas = original.uniform_state().unwrap().state();
    for end in [ModalEnd::Rigid, ModalEnd::Ports] {
        let mut uniform = ContinuousModalLine::new(
            &profile,
            gas,
            ModeCount(2),
            3,
            end,
            1.0 / 96000.0,
            budget(),
            || false,
        )
        .unwrap();
        let mut sourced = ContinuousModalLine::new_field(
            &profile,
            &original,
            ModeCount(2),
            3,
            end,
            1.0 / 96000.0,
            budget(),
            || false,
        )
        .unwrap();
        assert_eq!(sourced.source().gas_field(), Some(&original));
        assert_eq!(uniform.state().len(), sourced.state().len());
        for step in 0..24 {
            let mut input = vec![0.0; sourced.port_count()];
            input[0] = if step < 8 { 70.0 } else { 0.0 };
            input[1] = 40.0 * det::sin(step as f64 * 0.13);
            if input.len() == 4 {
                input[3] = -25.0;
            }
            let a = uniform.push(&input).unwrap();
            let b = sourced.push(&input).unwrap();
            assert_eq!(a.physical_heat.to_bits(), b.physical_heat.to_bits());
            assert_eq!(
                snapshot(&uniform)
                    .into_iter()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>(),
                snapshot(&sourced)
                    .into_iter()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>()
            );
        }
    }
}

#[test]
fn spatial_source_and_original_mouth_memory_share_physical_work_and_clock_subdivision() {
    let profile = RadiusProfile::exponential(0.025, 0.0, 0.04).unwrap();
    let original = field(profile.length(), false);
    let outlet = *original.state_at(profile.length()).unwrap().state();
    let duration = 1.0 / 2400.0;
    let mut results = Vec::new();
    for subdivision in [1, 2, 4] {
        let load = coupled_load(profile.outlet_radius(), outlet, &[1.0, 0.6]);
        let mut line = ContinuousModalLine::new_field_loaded(
            &profile,
            &original,
            ModeCount(2),
            3,
            load,
            duration / subdivision as f64,
            budget(),
            || false,
        )
        .unwrap();
        assert_eq!(line.source().gas_field(), Some(&original));
        assert_eq!(line.source().outlet_gas(), &outlet);
        assert_eq!(line.load_state().len(), 2);
        assert_eq!(line.state().len(), 18);
        let (mut incident, mut outgoing, mut heat) = (0.0, 0.0, 0.0);
        for _ in 0..subdivision {
            let record = line.push(&[75.0, 0.0]).unwrap();
            incident += record.incident_energy;
            outgoing += record.outgoing_energy;
            heat += record.physical_heat;
            assert_eq!(record.external_work, 0.0);
        }
        assert!(heat > 0.0, "the original resistive load must absorb energy");
        assert!(line.load_state().iter().any(|v| v.abs() > 1e-12));
        close(
            line.stored_energy_j() + heat + outgoing,
            incident,
            3e-8,
            1e-14,
        );
        results.push((
            line.state().to_vec(),
            line.stored_energy_j(),
            heat,
            outgoing,
        ));
    }
    for result in &results[1..] {
        for (actual, expected) in result.0.iter().zip(&results[0].0) {
            close(*actual, *expected, 2e-7, 1e-12);
        }
        close(result.1, results[0].1, 2e-7, 1e-14);
        close(result.2, results[0].2, 2e-7, 1e-14);
        close(result.3, results[0].3, 2e-7, 1e-14);
    }
}
