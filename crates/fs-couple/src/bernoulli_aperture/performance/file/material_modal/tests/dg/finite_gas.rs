//! Heat changes actual retained material modes and the original DG pressure.
use super::{pack, render, valve, CancelGate, LoadedAperture, MaterialModalDgValvePerformance,
    PlateValveInputError, PressureRenderer};
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};

#[path = "../../../../../../render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const INPUT: &str = include_str!("../../../../../../../examples/material_modal_valve/paired-dg-finite-gas.performance");

fn input() -> String {
    let nodes = "pressure_node 1\n".repeat(21);
    let loads = format!(
        "pressure_shapes 1\npressure_shape 0 1000000 finite-gas-pressure\npressure_nodes 21\n{nodes}\
source_forces 2\nsource_force 3 1 -0.0001\nsource_force 31 1 0\n\
source_pressures 2\nsource_pressure 5 0 1\nsource_pressure 32 0 0\n"
    );
    INPUT.replace("audio 96000 96000 4800 20", "audio 96000 96000 64 20")
        .replace("95000", "101325").replace("N2", "A").replace("CO2", "B")
        .replace("heat 0 0.01", "heat 0 1")
        .replace("heat 1600 0", "heat 16 0")
        .replace("heat 3200 -0.005", "heat 32 -0.2")
        .replace("compile_limits", &format!("{loads}compile_limits"))
}

fn build(text: &str, block: usize, sources: Option<&ThermochemicalSources<'_>>)
    -> Result<MaterialModalDgValvePerformance, PlateValveInputError>
{
    MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), block, &CancelGate::new_clock_free(), sources,
        |_| Ok(pack(900.0, false)),
    )
}

#[test]
fn heat_changes_modal_motion_flow_and_pressure_on_the_same_material_clock() {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None,
        phase: ThermochemicalPhase::FrozenGas };
    let text = input();
    let source = build(&text, 7, Some(&sources)).unwrap();
    assert_eq!(source.info().modes, 6);
    assert_eq!(source.material_records().len(), 4);
    let original_modes = valve(&source).spec().modes.clone();
    let original_materials = source.material_records().iter()
        .map(|record| record.pack().content_hash()).collect::<Vec<_>>();
    let initial = source.renderer().finite_gas().unwrap().clone();
    let mut heated = source.into_renderer();
    let mut pressure = Vec::new();
    let mut force_work = 0.0;
    let mut spatial_work = 0.0;
    for _ in 0..64 {
        let mut sample = [0.0];
        heated.block(&mut sample).unwrap();
        pressure.push(sample[0]);
        let frame = heated.last_frame().unwrap();
        force_work += frame.aperture.mechanical_work_j.abs();
        spatial_work += frame.spatial_pressure[0].work_j.abs();
        let report = heated.last_finite_gas_step().unwrap();
        assert!(report.gas.gas_energy_residual_j.abs() <= 1e-11);
        assert_eq!(report.coupled.aperture.bore_pressure_pa, sample[0]);
    }
    let gas = heated.finite_gas().unwrap();
    assert_eq!(gas.mass_kg(), initial.mass_kg());
    assert_eq!(gas.volume_m3(), initial.volume_m3());
    assert!(gas.state().state().temperature > 300.1);
    assert!(gas.state().state().pressure > initial.state().state().pressure);
    assert!(force_work > 0.0 && spatial_work > 0.0);
    assert_eq!(heated.applied_heat_controls(),
        heated.finite_gas_heat_program().unwrap().controls.as_slice());
    assert_eq!(heated.last_finite_gas_step().unwrap().gas.supplied_heat_j, -0.2 * (1.0 / 96000.0));
    assert_eq!(heated.system().load().state().step_count(), 64);
    assert_eq!(heated.system().load().duct().gas_field().unwrap().uniform_state(), Some(gas.state()));
    let LoadedAperture::Modal(material) = heated.system().aperture() else {
        panic!("thermal playback must retain the original material-modal valve");
    };
    assert_eq!(material.spec().modes, original_modes);
    assert_eq!(material.gas(), Some(heated.system().load().duct().inlet_gas()));
    assert_eq!(material.spec().impedance_pa_s_m3.to_bits(),
        heated.system().load().impedance_pa_s_m3().to_bits());
    let reference = text.replace("heat 0 1", "heat 0 0").replace("heat 32 -0.2", "heat 32 0");
    let (cold_pressure, cold) = render(build(&reference, 7, Some(&sources)).unwrap());
    let LoadedAperture::Modal(cold_material) = cold.system().aperture() else {
        panic!("original reference modal source");
    };
    assert_ne!(material.modal_states(), cold_material.modal_states());
    assert_ne!(heated.last_frame().unwrap().aperture.bore_flow_m3_s,
        cold.last_frame().unwrap().aperture.bore_flow_m3_s);
    assert!(pressure.iter().zip(cold_pressure).any(|(hot, cold)| (hot - cold).abs() > 1e-10));
    for block in [1, 7, 64] {
        let source = build(&text, block, Some(&sources)).unwrap();
        assert_eq!(source.material_records().iter()
            .map(|record| record.pack().content_hash()).collect::<Vec<_>>(), original_materials);
        let (other, renderer) = render(source);
        assert_eq!(pressure, other);
        assert_eq!(heated.finite_gas(), renderer.finite_gas());
        assert_eq!(heated.system().load().state(), renderer.system().load().state());
        let LoadedAperture::Modal(other_material) = renderer.system().aperture() else {
            panic!("same original modal source");
        };
        assert_eq!(material.modal_states(), other_material.modal_states());
        assert_eq!(heated.applied_heat_controls(), renderer.applied_heat_controls());
    }
}

#[test]
fn rejected_heat_keeps_the_original_material_inventory_wave_and_controls() {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None,
        phase: ThermochemicalPhase::FrozenGas };
    let source = build(&input(), 7, Some(&sources)).unwrap();
    let modes = valve(&source).modal_states().to_vec();
    let mut renderer = source.into_renderer();
    let gas = renderer.finite_gas().unwrap().clone();
    let wave = renderer.system().load().state().clone();
    renderer.finite_gas_cancel_gate().unwrap().request();
    assert!(renderer.block(&mut [0.0]).is_err());
    assert_eq!(renderer.finite_gas(), Some(&gas));
    assert_eq!(renderer.system().load().state(), &wave);
    assert_eq!(renderer.samples_rendered(), 0);
    assert!(renderer.applied_controls().is_empty());
    assert!(renderer.applied_heat_controls().is_empty());
    assert!(renderer.last_frame().is_none());
    assert!(renderer.last_finite_gas_step().is_none());
    let LoadedAperture::Modal(material) = renderer.system().aperture() else {
        panic!("same original modal source");
    };
    assert_eq!(material.modal_states(), modes);
    assert!(matches!(renderer.block(&mut [0.0]), Err(crate::render::RenderError::Poisoned)));
}

#[test]
fn native_heat_requires_its_original_inventory_and_shared_control_budget() {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None,
        phase: ThermochemicalPhase::FrozenGas };
    let text = input();
    assert!(build(&text, 7, None).is_err());
    for refused in [
        text.replace("finite_gas thermochemical 300", "finite_gas thermochemical 301"),
        text.replace("gas_station 0.05 300", "gas_station 0.05 301"),
        text.replace("heat 0 1", "heat 1 1"),
        text.replace("heat 16 0", "heat 0 0"),
        text.replace("heat 32 -0.2", "heat 64 -0.2"),
        text.replace("compile_limits 100000 8192", "compile_limits 100000 7"),
        text.replace("finite_gas thermochemical", "supply slit 0.00002 0.01 0.001 16\nfinite_gas thermochemical"),
    ] {
        assert!(build(&refused, 7, Some(&sources)).is_err(), "{refused}");
    }
    let exterior = text.replace("dg_bore 1 open", "dg_bore 1 baffled")
        .replace("\nfinite_gas", concat!(
            "\ndg_radiation computed-baffled 20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257\n",
            "dg_exterior thermochemical 300 101325 mole 1 A 1\nfinite_gas"));
    let error = build(&exterior, 7, Some(&sources)).err()
        .expect("thermal inventory cannot replace open-flow mean-pressure exchange");
    assert!(error.to_string().contains("closed-inventory"), "{error}");
}
