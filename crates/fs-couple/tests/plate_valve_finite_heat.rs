//! Authored heat follows the original material inventory and pressure clock.
use fs_couple::bernoulli_aperture::performance::{AperturePerformance, CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::{PlateValveInputError, PlateValvePerformance};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

fn performance() -> String {
    include_str!("../examples/plate-valve-finite-gas.performance")
        .replace("audio 96000 9600 20", "audio 96000 64 20")
        .replace("95000", "101325")
        .replace("N2", "A").replace("CO2", "B")
        .replace("heat 0 0.01", "heat 0 1")
        .replace("heat 3200 0", "heat 16 0")
        .replace("heat 6400 -0.005", "heat 32 -0.2")
}

fn load(text: &str, sources: Option<&ThermochemicalSources<'_>>, block: usize)
    -> Result<AperturePerformance, PlateValveInputError>
{
    PlateValvePerformance::from_bytes_with_gas_sources(
        text.as_bytes(), block, &CancelGate::new_clock_free(), sources,
    ).map(PlateValvePerformance::into_renderer)
}
fn render(text: &str, sources: &ThermochemicalSources<'_>, block: usize)
    -> (Vec<f64>, AperturePerformance)
{
    let mut renderer = load(text, Some(sources), block).unwrap();
    let mut pressure = vec![0.0; 64];
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    assert_eq!(renderer.samples_rendered(), 64);
    assert!(pressure.iter().all(|p| p.is_finite()));
    (pressure, renderer)
}

#[test]
fn heat_warms_the_same_inventory_changes_pressure_and_preserves_callback_identity() {
    let a = support::packs("A", 0.028, 3.5, 200.0);
    let b = support::packs("B", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let text = performance();
    let original = load(&text, Some(&sources), 7).unwrap();
    let initial = original.finite_gas().unwrap();
    assert!((initial.volume_m3() / (core::f64::consts::PI * 1e-6 * 0.05) - 1.0).abs() < 1e-14);
    assert!(original.applied_heat_controls().is_empty());
    let (pressure, heated) = render(&text, &sources, 7);
    let gas = heated.finite_gas().unwrap();
    assert_eq!(gas.mass_kg(), initial.mass_kg());
    assert_eq!(gas.volume_m3(), initial.volume_m3());
    assert!(gas.state().state().temperature > 300.1);
    assert!(gas.state().state().pressure > initial.state().state().pressure);
    let CoupledAperture::Dg(system) = heated.system() else { panic!("original DG owner required") };
    assert_eq!(system.load().duct().gas_field().unwrap().uniform_state(), Some(gas.state()));
    assert_eq!(system.aperture().accepted_steps(), 64);
    assert_eq!(heated.applied_heat_controls(), heated.finite_gas_heat_program().unwrap().controls.as_slice());
    let last = heated.last_finite_gas_step().unwrap();
    assert_eq!(last.gas.supplied_heat_j, -0.2 * (1.0 / 96000.0));
    assert!(last.gas.gas_energy_residual_j.abs() <= 1e-11);
    let cold = text.replace("heat 0 1", "heat 0 0").replace("heat 32 -0.2", "heat 32 0");
    let (reference, _) = render(&cold, &sources, 7);
    assert!(pressure.iter().zip(reference).any(|(hot, cold)| (hot - cold).abs() > 1e-10));
    for block in [1, 64] {
        let (other, renderer) = render(&text, &sources, block);
        assert_eq!(pressure, other);
        assert_eq!(heated.finite_gas(), renderer.finite_gas());
        assert_eq!(heated.applied_heat_controls(), renderer.applied_heat_controls());
    }
}

#[test]
fn heat_requires_the_original_uniform_source_frozen_phase_and_bounded_controls() {
    let a = support::packs("A", 0.028, 3.5, 200.0);
    let b = support::packs("B", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let text = performance();
    assert!(load(&text, None, 7).is_err());
    for input in [
        text.replace("finite_gas thermochemical 300", "finite_gas thermochemical 301"),
        text.replace("finite_gas thermochemical 300 101325 mole 1 A 1",
            "finite_gas thermochemical 300 101325 mole 1 B 1"),
        text.replace("gas_station 0.05 300", "gas_station 0.05 301"),
        text.replace("dg_segment cylinder 0.001 0.05", "dg_segment cone 0.001 0.0011 0.05"),
        text.replace("gas_energy 200 900", "gas_energy 301 900"),
        text.replace("gas_feedback -1e-5 1e-5", "gas_feedback 1e-5 1e-5"),
        text.replace("heat 0 1", "heat 1 1"),
        text.replace("heat 16 0", "heat 0 0"),
        text.replace("heat 32 -0.2", "heat 64 -0.2"),
        text.replace("heat 32 -0.2", "heat 32 NaN"),
    ] {
        assert!(load(&input, Some(&sources), 7).is_err(), "{input}");
    }
    let water = fs_qty::chemistry::SpeciesId::new("H2O").unwrap();
    let w = support::packs("H2O", 0.018, 4.0, 200.0);
    let phase_bindings = [a.binding(), w.binding()];
    let phase = ThermochemicalSources { species: &phase_bindings, water_species: Some(&water), phase: ThermochemicalPhase::LiquidWaterLimit };
    let error = load(&text.replace(" B", " H2O"), Some(&phase), 7).err().expect("future phase must be explicit");
    assert!(error.to_string().contains("FrozenGas"), "{error}");
    let exterior = text.replace("dg_bore 1 open", "dg_bore 1 baffled")
        .replace("\nfinite_gas", concat!(
            "\ndg_radiation computed-baffled 20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257\n",
            "dg_exterior thermochemical 300 101325 mole 1 A 1\nfinite_gas"));
    let error = load(&exterior, Some(&sources), 7).err().expect("mean-pressure exchange needs transport");
    assert!(error.to_string().contains("closed-inventory"), "{error}");
}

#[test]
fn a_rejected_heating_sample_keeps_the_original_inventory_and_source_clock() {
    let a = support::packs("A", 0.028, 3.5, 200.0);
    let b = support::packs("B", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let mut renderer = load(&performance(), Some(&sources), 7).unwrap();
    let before = renderer.finite_gas().unwrap().clone();
    renderer.finite_gas_cancel_gate().unwrap().request();
    assert!(renderer.block(&mut [0.0]).is_err());
    assert_eq!(renderer.finite_gas(), Some(&before));
    assert!(renderer.applied_heat_controls().is_empty());
    assert!(renderer.last_finite_gas_step().is_none());
    let CoupledAperture::Dg(system) = renderer.system() else { panic!("original DG owner required") };
    assert_eq!(system.aperture().accepted_steps(), 0);
}
