//! Data-selected supply momentum stays inside the original material/DG solve.
use super::{
    CancelGate, LoadedAperture, MaterialModalDgValvePerformance, build, input, pack, render, valve,
};

const SLIT: &str = "supply slit 0.004 0.002 0.0002 16\n";
const CIRCULAR: &str = "supply circular 0.004 0.0002 16\n";
const TRACT: &str = "upstream_tract\ndg_tube 0.04 0.001 closed 1 2 8 50 2000 0.01 1e-11 256\n";

#[test]
fn optional_supply_none_preserves_existing_native_samples_and_histories() {
    let original = render(build(&input(), 37));
    let explicit = render(build(
        &input().replace("slit_length", "supply none\nslit_length"),
        37,
    ));
    assert_eq!(original.0, explicit.0);
    assert_eq!(
        original.1.system().load().state(),
        explicit.1.system().load().state()
    );
    assert!(explicit.1.system().supply_channel().is_none());
    assert!(explicit.1.last_frame().unwrap().supply.is_none());
}

#[test]
fn both_native_geometries_reach_original_material_modes_and_are_block_invariant() {
    let plain = render(build(&input(), 37));
    for declaration in [SLIT, CIRCULAR] {
        let text = input().replace("slit_length", &format!("{declaration}slit_length"));
        let source = build(&text, 37);
        let system = source.renderer().system();
        let channel = system.supply_channel().unwrap();
        assert_eq!(channel.gas(), system.load().duct().inlet_gas());
        assert_eq!(channel.gas(), valve(&source).gas().unwrap());
        assert_eq!(channel.sample_period_s(), system.load().propagator().dt());
        assert_eq!(
            valve(&source).spec().impedance_pa_s_m3,
            system.load().impedance_pa_s_m3() + channel.pressure_load().unwrap().impedance_pa_s_m3
        );
        assert_eq!(
            source.info().modes,
            6,
            "the original two independent material plates survive"
        );
        let (pressure, reference) = render(source);
        assert!(pressure.iter().all(|p| p.is_finite()));
        assert_ne!(
            plain.0, pressure,
            "the native fluid momentum must feed the acoustic source"
        );
        assert!(
            reference
                .system()
                .supply_channel()
                .unwrap()
                .stored_energy_j()
                > 0.0
        );
        let LoadedAperture::Modal(plain_valve) = plain.1.system().aperture() else {
            unreachable!()
        };
        let LoadedAperture::Modal(reference_valve) = reference.system().aperture() else {
            unreachable!()
        };
        assert_ne!(plain_valve.modal_states(), reference_valve.modal_states());
        for block in [1, 512] {
            let (partitioned, other) = render(build(&text, block));
            assert_eq!(pressure, partitioned);
            assert_eq!(
                reference.system().load().state(),
                other.system().load().state()
            );
            assert_eq!(
                reference
                    .system()
                    .supply_channel()
                    .unwrap()
                    .velocities_m_s(),
                other.system().supply_channel().unwrap().velocities_m_s()
            );
            let LoadedAperture::Modal(other_valve) = other.system().aperture() else {
                unreachable!()
            };
            assert_eq!(reference_valve.modal_states(), other_valve.modal_states());
        }
    }
}

#[test]
fn geometry_pressure_and_species_data_change_one_coupled_tract_channel_reed_bore() {
    let text = input()
        .replace("slit_length", &format!("{TRACT}{SLIT}slit_length"))
        .replace(
            "ambient 293.15 101325 0",
            "ambient mixture 293.15 101325 mole 2 N2 0.79 O2 0.21",
        );
    let source = build(&text, 37);
    let original_gas = *source.renderer().system().supply_channel().unwrap().gas();
    let (reference, original) = render(source);
    assert_eq!(
        original
            .system()
            .upstream_load()
            .unwrap()
            .state()
            .step_count(),
        512
    );
    assert!(original.last_frame().unwrap().upstream.is_some());
    assert!(original.last_frame().unwrap().supply.is_some());
    for changed in [
        text.replace("supply slit 0.004", "supply slit 0.008"),
        text.replace("101325", "81060"),
        text.replace("N2 0.79 O2 0.21", "He 0.79 O2 0.21"),
    ] {
        let source = build(&changed, 37);
        let system = source.renderer().system();
        let channel = system.supply_channel().unwrap();
        assert_eq!(channel.gas(), system.load().duct().inlet_gas());
        assert_eq!(
            channel.gas(),
            system.upstream_load().unwrap().duct().inlet_gas()
        );
        assert_eq!(channel.gas(), valve(&source).gas().unwrap());
        if changed.contains("He 0.79") {
            assert!(channel.gas().density < original_gas.density);
            assert!(channel.gas().sound_speed > original_gas.sound_speed);
            assert_ne!(
                channel.gas().dynamic_viscosity,
                original_gas.dynamic_viscosity
            );
        }
        let (pressure, other) = render(source);
        assert_ne!(reference, pressure);
        assert_ne!(
            original.system().supply_channel().unwrap().velocities_m_s(),
            other.system().supply_channel().unwrap().velocities_m_s()
        );
        let LoadedAperture::Modal(a) = original.system().aperture() else {
            unreachable!()
        };
        let LoadedAperture::Modal(b) = other.system().aperture() else {
            unreachable!()
        };
        assert_ne!(a.modal_states(), b.modal_states());
    }
}

#[test]
fn malformed_supply_records_refuse_without_a_substitute_passage() {
    for declaration in [
        "supply other 0.004 0.002 0.0002 16\n",
        "supply slit 0.004 0.002 0.0002 1\n",
        "supply circular 0.004 0.0002 4097\n",
        "supply slit 0 0.002 0.0002 16\n",
        "supply circular 0.004 -0.0002 16\n",
        "supply none extra\n",
    ] {
        let text = input().replace("slit_length", &format!("{declaration}slit_length"));
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(pack(900.0, false)),
            )
            .is_err(),
            "accepted {declaration}"
        );
    }
}
