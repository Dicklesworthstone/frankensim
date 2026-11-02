//! Manufactured composition checks, not measured flute or breath validation.
use super::*;
use crate::nsf_dg::{NsfDgDrive, NsfFiniteReservoirInventory, NsfFiniteReservoirSource};
use fs_lbm::nsf3::NsfGasMedium;
use fs_material::gas_species::input::ThermochemicalPhase;

#[path = "../../../../render/schedule/reed/field_input/tests/support.rs"]
mod support;

const INPUT: &str = include_str!("../../../../../examples/resolved-thermoviscous-breath.performance");
const MIXTURE: &str = "mixture 300 95000 mole 4 N2 0.75 O2 0.20 CO2 0.04 H2O 0.01";
const MOUTH: &str = "reservoir_observation aperture 1 baffled -0.08 0 0\n";

fn input() -> String {
    INPUT.replace("output 1920 1", "output 384 1")
        .replace("boundary finite-reservoir-actuation 0 0.01 4\nramp 384 0 0 -1e-7 0.0005\nhold 768 -1e-7 0.0005\nramp 1152 -1e-7 0.0005 1e-7 0\nhold 1920 1e-7 0",
            "boundary finite-reservoir-actuation 0 0.01 1\nhold 384 -1e-7 0.0005")
}

fn read(text: &str, block: usize) -> ResolvedThermoviscousDgInput {
    ResolvedThermoviscousDgInput::from_bytes(text.as_bytes(), block,
        Arc::new(CancelGate::new_clock_free())).unwrap()
}

fn render(source: ResolvedThermoviscousDgInput, block: usize)
    -> (Vec<f64>, u64, NsfFiniteReservoirInventory)
{
    source.with_physical_renderer(|renderer, info| {
        let mut pressure = vec![0.0; info.samples as usize];
        for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        assert!(pressure.iter().all(|p| p.is_finite()));
        assert!(renderer.validate_sample_count(1).is_err());
        Ok((pressure, renderer.system().zone().state_hash(),
            renderer.system().finite_reservoir_inventory().unwrap()))
    }).unwrap()
}

#[test]
fn original_inventory_and_every_material_microphone_survive_the_native_field_handoff() {
    let text = input();
    let source = read(&text, 37);
    let original = source.source().renderer();
    let inventory = original.system().finite_reservoir_inventory().unwrap();
    let zone_hash = original.system().zone().state_hash();
    let duct_state = original.system().load().state().clone();
    let reference = *original.system().reference_gas();
    let propagation = original.observer_gas().unwrap();
    let gas = propagation.gas;
    let law = propagation.law;
    let latency = propagation.additional_latency_samples;
    let mut bore = original.baffled_receiver().unwrap().clone();
    let mut escape = original.reservoir_receiver_at(1).unwrap().clone();
    let result = source.with_physical_renderer(|renderer, info| {
        assert_eq!(renderer.system().finite_reservoir_inventory(), Some(inventory));
        assert_eq!(renderer.system().zone().state_hash(), zone_hash);
        assert_eq!(renderer.system().load().system().load().state(), &duct_state);
        assert_eq!(renderer.system().load().system().source().gas(), &reference);
        assert_eq!(renderer.observer_gas().unwrap().gas, gas);
        assert_eq!(renderer.observer_gas().unwrap().law, law);
        assert_eq!(renderer.observer_gas().unwrap().additional_latency_samples, latency);
        assert_eq!(renderer.reservoir_observer_gas_at(1).unwrap().additional_latency_samples, latency);
        assert_eq!(renderer.reservoir_observation_count(), 1);
        assert!(renderer.reservoir_receiver_at(0).is_none());
        assert!(renderer.drive_history().iter().all(|drive| matches!(drive,
            NsfDgDrive::FiniteReservoirActuation { aperture: 0, .. })));

        let mut pressure = vec![0.0; info.samples as usize];
        let mut escape_peak = 0.0_f64;
        for sample in &mut pressure {
            renderer.block(core::slice::from_mut(sample)).unwrap();
            // These are the independently retained original physical receiver
            // kernels, driven by the accepted interval-overlap flow projections.
            let bore_pa = bore.step_trace(&renderer.terminal_flows).unwrap();
            let escape_pa = escape.step_trace(&renderer.reservoirs[0].flows).unwrap();
            assert_eq!(sample.to_bits(), (bore_pa + escape_pa).to_bits());
            escape_peak = escape_peak.max(escape_pa.abs());
            let frame = renderer.last_frame().unwrap();
            let stored = renderer.system().finite_reservoir_inventory().unwrap();
            let account = frame.finite_reservoir.as_ref().unwrap();
            assert_eq!(account.source_exchange, frame.reservoir_exchanges[0]);
            assert_eq!(account.final_mass_kg, stored.mass_kg);
            assert_eq!(account.final_volume_m3, stored.volume_m3);
            assert_eq!(account.final_pressure_pa, stored.gas.pressure);
            assert!(account.mass_residual_kg.abs() <= 1e-18);
            assert!(account.energy_residual_j.abs() <= 1e-10);
            assert!(account.actuation.as_ref().unwrap().step.work_by_gas_j < 0.0);
        }
        assert!(escape_peak > 0.0, "actual escape flow must reach the microphone");
        let final_inventory = renderer.system().finite_reservoir_inventory().unwrap();
        assert!(final_inventory.volume_m3 < inventory.volume_m3);
        assert_ne!(final_inventory.gas.pressure, inventory.gas.pressure);
        Ok((pressure, renderer.system().zone().state_hash(), final_inventory))
    }).unwrap();
    assert_eq!(result, render(read(&text, 23), 23));
    let bore_only = render(read(&text.replace(MOUTH, ""), 23), 23);
    assert_eq!(result.1, bore_only.1, "observation must not change source state");
    assert_eq!(result.2, bore_only.2, "observation must not change material storage");
    assert_ne!(result.0, bore_only.0, "escape pressure is independent of bore flow");
}

fn sourced(cp_r: f64) -> ResolvedThermoviscousDgInput {
    let packs = support::packs("arbitrary-source", 0.028, cp_r, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let text = input().replace(MIXTURE, "thermochemical 300 95000 mole 1 arbitrary-source 1")
        .replace("\nbulk_viscosity", "\nenergy_solve 1e-8 96\nbulk_viscosity");
    ResolvedThermoviscousDgInput::from_bytes_with_sources(text.as_bytes(), 29,
        Arc::new(CancelGate::new_clock_free()), Some(&sources),
        |_| Err("this manufactured source has no material wall".into())).unwrap()
} // The source text and borrowed pack catalogue are gone before playback.

#[test]
fn supplied_caloric_source_changes_the_complete_field_and_sound_without_geometry_retuning() {
    let first = sourced(3.5);
    let second = sourced(4.5);
    let a = first.source().renderer().system();
    let b = second.source().renderer().system();
    assert_eq!(a.zone().law().dx_m, b.zone().law().dx_m);
    assert_eq!(a.zone().law().dt_s, b.zone().law().dt_s);
    assert_eq!(a.reference_gas().density, b.reference_gas().density);
    assert_ne!(a.reference_gas().specific_heat_cp, b.reference_gas().specific_heat_cp);
    assert_ne!(a.reference_gas().sound_speed, b.reference_gas().sound_speed);
    for source in [&first, &second] {
        let system = source.source().renderer().system();
        assert!(matches!(system.zone().medium(), NsfGasMedium::Thermochemical(_)));
        assert_eq!(system.zone().thermochemical_state([5, 1, 1]).unwrap().unwrap().source_packs().len(), 1);
        assert_eq!(system.finite_reservoir_inventory().unwrap().source_kind,
            NsfFiniteReservoirSource::Thermochemical);
    }
    let a = render(first, 29);
    let b = render(second, 29);
    assert_ne!(a.0, b.0);
    assert_ne!(a.2.gas.temperature, b.2.gas.temperature);
    assert_eq!(a.2.source_kind, NsfFiniteReservoirSource::Thermochemical);
    assert_eq!(b.2.source_kind, NsfFiniteReservoirSource::Thermochemical);
}

#[test]
fn source_loader_refusals_do_not_fall_back_to_a_different_wall_or_gas() {
    let text = input().replace("\nobservation baffled",
        "\ndg_wall isotropic missing.pack 300 1000 0.00002 stress-free 20 2000 elastic\nobservation baffled");
    let mut requested = Vec::new();
    let refused = ResolvedThermoviscousDgInput::from_bytes_with_sources(text.as_bytes(), 8,
        Arc::new(CancelGate::new_clock_free()), None, |name| {
            requested.push(name.to_owned());
            Err("original-material-loader-refusal".into())
        });
    assert_eq!(requested, ["missing.pack"]);
    match refused {
        Err(error) => assert!(error.to_string().contains("original-material-loader-refusal")),
        Ok(_) => panic!("missing wall data must not select a replacement material"),
    }
    let missing = input().replace(MIXTURE, "thermochemical 300 95000 mole 1 arbitrary-source 1")
        .replace("\nbulk_viscosity", "\nenergy_solve 1e-8 96\nbulk_viscosity");
    assert!(ResolvedThermoviscousDgInput::from_bytes(missing.as_bytes(), 8,
        Arc::new(CancelGate::new_clock_free())).is_err());
}

#[test]
fn cancellation_and_cold_admission_preserve_all_accepted_physical_owners() {
    let text = input();
    let gate = Arc::new(CancelGate::new_clock_free());
    let source = ResolvedThermoviscousDgInput::from_bytes(text.as_bytes(), 8, Arc::clone(&gate)).unwrap();
    source.with_physical_renderer(|renderer, _| {
        renderer.block(&mut [0.0; 8]).unwrap();
        let before = (renderer.system().zone().state_hash(),
            renderer.system().finite_reservoir_inventory(), renderer.samples_rendered(),
            renderer.system().load().system().state().clone(),
            renderer.system().load().system().load().state().clone());
        gate.request();
        let mut output = [12345.0; 1];
        assert!(renderer.block(&mut output).is_err());
        assert_eq!(output, [12345.0]);
        assert_eq!(before, (renderer.system().zone().state_hash(),
            renderer.system().finite_reservoir_inventory(), renderer.samples_rendered(),
            renderer.system().load().system().state().clone(),
            renderer.system().load().system().load().state().clone()));
        Ok(())
    }).unwrap();

    let mut advanced = read(&text, 1).source.into_renderer();
    advanced.block(&mut [0.0]).unwrap();
    assert!(advanced.into_cold_parts().is_err());
    let upstream = text.replace("observation baffled",
        "upstream_dg dg_tube 0.03 0.0008462843753216344 matched 2 2 12 20 2000 0.01 1e-11 512\nobservation baffled");
    let source = read(&upstream, 8);
    let invoked = std::cell::Cell::new(false);
    assert!(source.with_renderer(|_, _| { invoked.set(true); Ok(()) }).is_err());
    assert!(!invoked.get(), "unsupported acoustic topology cannot be discarded");
}
