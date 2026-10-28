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
    let nonrelaxing = text.replace("slit_length", concat!(
        "dg_molecular declared 2 tracta 0 tractb 0\nslit_length",
    ));
    let (same, _) = render(build_gas(&nonrelaxing).unwrap());
    assert_eq!(reference, same, "explicit zero-reservoir field preserves the original response");
}

#[test]
fn molecular_tract_keeps_remote_temperature_and_composition_in_reed_feedback() {
    let a = gas_support::packs("tracta", 0.028, 3.5, 200.0);
    let b = gas_support::packs("tractb", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let gate = CancelGate::new_clock_free();
    let build_gas = |text: &str, block| MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), block, &gate, Some(&sources), |_| Ok(pack(900.0, false)),
    ).unwrap();
    // Manufactured constitutive data isolate the native shared-law connection.
    let law = "dg_molecular declared 2 tracta 1 vibration 1 1000 1 power 10 293.15 0 273.15 323.15 tractb 0\n";
    let tract = format!(concat!(
        "upstream_tract\n",
        "dg_bore 1 closed 1 2 12 50 2000 0.01 1e-11 256 4 gas_field 2 mole 2 tracta tractb\n",
        "dg_segment cone 0.001 0.0015 0.04\n",
        "gas_station 0 293.15 0.7 0.3\n",
        "gas_station 0.04 303.15 0.5 0.5\n{}"), law);
    let text = input().replace("ambient 293.15 101325 0",
        "ambient thermochemical 293.15 101325 mole 2 tracta 0.7 tractb 0.3")
        .replace("slit_length", &format!("{tract}slit_length"));
    let original = build_gas(&text, 37);
    let load = original.renderer().system().upstream_load().unwrap();
    assert!(load.duct().molecular_relaxation_field().is_some());
    assert!(load.duct().molecular_relaxation().is_none());
    assert_eq!(load.duct().total_dofs(), 150);
    let (reference, accepted) = render(original);
    let (partitioned, other) = render(build_gas(&text, 512));
    assert_eq!(reference, partitioned);
    assert_eq!(states(&accepted), states(&other));
    assert_eq!(accepted.system().upstream_load().unwrap().state(),
        other.system().upstream_load().unwrap().state());
    for changed in [
        text.replace(law, ""),
        text.replace("power 10 ", "power 1 "),
        text.replace("gas_station 0.04 303.15 0.5 0.5", "gas_station 0.04 313.15 0.5 0.5"),
        text.replace("gas_station 0.04 303.15 0.5 0.5", "gas_station 0.04 303.15 0.3 0.7"),
    ] {
        let (pressure, other) = render(build_gas(&changed, 37));
        assert_ne!(reference, pressure);
        assert_ne!(states(&accepted), states(&other));
        assert_eq!(accepted.system().load().duct().inlet_gas(), other.system().load().duct().inlet_gas());
    }
}

#[test]
fn conical_bore_and_tract_molecular_states_reach_both_material_blades() {
    // Synthetic constitutive data exercise the full declared material/source
    // connection; they do not represent measured cane or gas relaxation rates.
    let law = "dg_molecular declared 2 N2 1 vibration 1 1000 1 power 10 293.15 0 273.15 323.15 O2 0\n";
    let bore = format!(concat!(
        "dg_bore 1 open 1 2 12 50 2000 0.01 1e-11 256 4\n",
        "dg_segment cone 0.001 0.002 0.05\n{}"), law);
    let tract = format!(concat!(
        "upstream_tract\n",
        "dg_bore 1 closed 1 2 12 50 2000 0.01 1e-11 256 4\n",
        "dg_segment cone 0.001 0.0015 0.04\n{}"), law);
    let text = input().replace("ambient 293.15 101325 0",
        "ambient mixture 293.15 101325 mole 2 N2 0.8 O2 0.2")
        .replace("dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256\n", &bore)
        .replace("slit_length", &format!("{tract}slit_length"));
    let source = build(&text, 37);
    assert_eq!(source.info().blades, 2);
    let main = source.renderer().system().load().duct();
    let upstream = source.renderer().system().upstream_load().unwrap().duct();
    for duct in [main, upstream] {
        assert!(duct.spatial_viscothermal_admission().is_some());
        assert_eq!(duct.molecular_relaxation().unwrap().gas(), duct.inlet_gas());
        assert_eq!(duct.total_dofs(), 150);
    }
    let (reference, accepted) = render(source);
    let (partitioned, other) = render(build(&text, 512));
    assert_eq!(reference, partitioned);
    assert_eq!(states(&accepted), states(&other));
    assert_eq!(accepted.system().load().state(), other.system().load().state());
    assert_eq!(accepted.system().upstream_load().unwrap().state(),
        other.system().upstream_load().unwrap().state());
    for without in [
        text.replacen(law, "", 1),
        text.replace(&format!("{law}slit_length"), "slit_length"),
    ] {
        let (changed, state) = render(build(&without, 37));
        assert_ne!(reference, changed,
            "each region's molecular histories must change the accepted bore pressure");
        assert_ne!(states(&accepted), states(&state));
    }
    assert_eq!(accepted.system().load().state().step_count(), 512);
    assert_eq!(accepted.system().upstream_load().unwrap().state().step_count(), 512);
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
