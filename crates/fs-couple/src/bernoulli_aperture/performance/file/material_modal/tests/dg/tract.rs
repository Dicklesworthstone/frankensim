//! Full native source and microphone retain original upstream physical owners.
use super::{build, input, pack, render, valve, CancelGate, LoadedAperture,
    MaterialModalDgValvePerformance, NormalizedMaterialCardPack};
use crate::bernoulli_aperture::dg::loaded::performance::LoadedDgPerformance;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};

#[path = "../../../../../../render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const TRACT: &str = "upstream_tract\ndg_tube 0.04 0.001 closed 1 2 8 50 2000 0.01 1e-11 256\n";
const WALL: &str = "dg_wall isotropic blade.matpack 293.15 100 0.00002 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time\n";

fn states(renderer: &LoadedDgPerformance) -> Vec<crate::bernoulli_aperture::paired::BladeState> {
    let LoadedAperture::Modal(source) = renderer.system().aperture() else { unreachable!() };
    source.modal_states().to_vec()
}

#[test]
fn material_tract_changes_reed_feedback_and_bore_microphone_on_the_original_clock() {
    let text = include_str!("../../../../../../../examples/material_modal_valve/paired-dg-chimney.performance")
        .replace("audio 96000 96000 4800 20", "audio 96000 96000 512 20")
        .replace("slit_length", &format!("{TRACT}{WALL}slit_length"));
    let mut loads = 0;
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        },
    ).unwrap();
    assert_eq!(loads, 1, "tract wall and reed must share the immutable material artifact");
    let wall = source.upstream_wall_layout().unwrap().uniform_record().unwrap();
    assert!(std::ptr::eq(wall.pack(), source.material_records()[0].pack()));
    let system = source.renderer().system();
    let tract = system.upstream_load().unwrap();
    assert_eq!(tract.duct().inlet_gas(), system.load().duct().inlet_gas());
    assert_eq!(valve(&source).spec().impedance_pa_s_m3,
        tract.impedance_pa_s_m3() + system.load().impedance_pa_s_m3());
    assert_eq!(tract.propagator().dt(), system.load().propagator().dt());
    assert!(tract.duct().has_viscothermal_memory());
    assert!((tract.duct().compliant_wall().unwrap().surface_density - 0.018).abs() < 1e-15);
    let (reference, rendered) = render(source);
    assert!(reference.iter().any(|p| p.abs() > 1e-15));
    let tract = rendered.system().upstream_load().unwrap();
    let motion = tract.duct().compliant_wall_state(tract.state()).unwrap();
    assert!(motion.stored_energy_j > 0.0 && motion.relaxation_stored_energy_j > 0.0);
    assert_eq!(tract.state().step_count(), 512);
    assert_eq!(rendered.system().load().state().step_count(), 512);
    assert!(rendered.last_frame().unwrap().upstream.is_some());
    for block in [1, 512] {
        let (other, source) = render(build(&text, block));
        assert_eq!(reference, other);
        assert_eq!(tract.state(), source.system().upstream_load().unwrap().state());
        assert_eq!(states(&rendered), states(&source));
    }
    let changed = text.replace("dg_wall isotropic blade.matpack", "dg_wall isotropic altered.matpack");
    let changed = build(&changed, 37);
    assert_eq!(changed.material_records()[0].pack(), &pack(900.0, false),
        "only the tract wall changed; every reed material and control stays fixed");
    let (other, changed) = render(changed);
    assert_ne!(reference, other);
    assert_ne!(states(&rendered), states(&changed));
    let (other, changed) = render(build(&text.replace("dg_tube 0.04", "dg_tube 0.08"), 37));
    assert_ne!(reference, other);
    assert_ne!(states(&rendered), states(&changed));
}

#[test]
fn remote_tract_composition_and_temperature_reach_the_same_material_reed() {
    let a = gas_support::packs("tracta", 0.028, 3.5, 200.0);
    let b = gas_support::packs("tractb", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let build_gas = |text: &str| MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), Some(&sources),
        |_| Ok(pack(900.0, false)),
    );
    let declaration = concat!(
        "upstream_tract\n",
        "dg_bore 1 closed 2 2 12 50 2000 0.01 1e-11 256 4 gas_field 3 mole 2 tracta tractb\n",
        "dg_segment cylinder 0.001 0.04\n",
        "gas_station 0 293.15 0.7 0.3\n",
        "gas_station 0.02 293.15 0.7 0.3\n",
        "gas_station 0.04 293.15 0.7 0.3\n",
    );
    let text = input().replace("ambient 293.15 101325 0",
        "ambient thermochemical 293.15 101325 mole 2 tracta 0.7 tractb 0.3")
        .replace("slit_length", &format!("{declaration}slit_length"));
    let base = build_gas(&text).unwrap();
    assert_eq!(base.renderer().system().upstream_load().unwrap().duct().inlet_gas(),
        base.renderer().system().load().duct().inlet_gas());
    let (reference, source) = render(base);
    for changed in [
        text.replace("gas_station 0.04 293.15 0.7 0.3", "gas_station 0.04 303.15 0.7 0.3"),
        text.replace("gas_station 0.04 293.15 0.7 0.3", "gas_station 0.04 293.15 0.5 0.5"),
    ] {
        let changed = build_gas(&changed).unwrap();
        assert_eq!(valve(&changed).gas(), Some(source.system().load().duct().inlet_gas()));
        let (other, changed) = render(changed);
        assert_ne!(reference, other, "only the remote upstream material-gas field changed");
        assert_ne!(states(&source), states(&changed));
    }
    let mismatched = text.replace("gas_station 0 293.15", "gas_station 0 303.15");
    assert!(build_gas(&mismatched).is_err(),
        "one small-signal through-flow requires matching complete inlet gases");
}

#[test]
fn tract_is_optional_and_retains_original_dg_source_budgets() {
    assert!(build(&input(), 37).renderer().system().upstream_load().is_none());
    for declaration in [
        "upstream_tract unexpected\n".to_owned(),
        TRACT.replace("1e-11 256", "1e-11 1"),
        format!("{TRACT}dg_wall regions 1\ndg_wall_region 0 0.05 rigid\n"),
    ] {
        let text = input().replace("slit_length", &format!("{declaration}slit_length"));
        assert!(MaterialModalDgValvePerformance::from_bytes_with_sources(
            text.as_bytes(), 37, &CancelGate::new_clock_free(), None,
            |_| -> Result<NormalizedMaterialCardPack, String> { Ok(pack(900.0, false)) },
        ).is_err());
    }
}
