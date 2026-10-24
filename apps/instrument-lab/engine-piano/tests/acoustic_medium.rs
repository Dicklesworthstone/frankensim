//! The Rayleigh receiver takes the caller's admitted gas as given (no dry-air
//! substitute), and a gas whose derived fields were mutated independently is
//! refused before anything is built.
use fs_live_piano::{Config, Description, LivePiano};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};

#[path = "../../../../crates/fs-couple/src/render/schedule/reed/field_input/tests/support.rs"]
#[allow(dead_code)] // the shared fixture also exposes a borrowed-registry helper
mod gas_sources;

fn source_packet(molar_mass: f64) -> (Vec<u8>, Vec<u32>) {
    let source = gas_sources::packs("piano-test-gas", molar_mass, 3.5, 200.0);
    let parts = [source.species.to_bytes(), source.nasa9.to_bytes(), source.transport.to_bytes()];
    let lengths = parts.iter().map(|part| u32::try_from(part.len()).unwrap()).collect();
    (parts.concat(), lengths)
}

const SOURCE_GAS: &str = "thermochemical 300 100000 mole 1 piano-test-gas 1";

/// G3: source bytes, without a named-species branch, reach actual material
/// string drag, bridge motion and observer pressure at both output clocks.
#[test]
fn canonical_gas_sources_change_material_string_motion_and_sound() {
    use fs_live_piano::{PianoAudio, preparation::{self, materials::{StringMaterials, operating_point}}};
    use fs_material::gas_species::input::GasRecord;
    use std::sync::Arc;
    let atmospheres = [0.028, 0.056].map(|mass| {
        let (bytes, lengths) = source_packet(mass);
        preparation::resolve_atmosphere_sources(SOURCE_GAS, &bytes, &lengths, "", "frozen").unwrap()
    }); // the original byte buffers and borrowed packs are gone
    let source_ids: Vec<_> = atmospheres.iter().map(|record| {
        let GasRecord::Thermochemical(gas) = record else { panic!("source ownership was discarded") };
        assert_eq!(gas.source_packs().len(), 1);
        gas.source_packs()[0].species_pack
    }).collect();
    assert_ne!(source_ids[0], source_ids[1]);
    let gases = atmospheres.each_ref().map(|record| *record.state());
    assert!((gases[1].density / gases[0].density - 2.0).abs() < 1e-13);
    assert!(gases[1].sound_speed < gases[0].sound_speed);

    let mut desc = Description::embedded().unwrap().truncated(4).unwrap();
    desc.courses.retain(|course| course.midi == 69);
    let construction = "frankensim-piano-materials-v2\ntension,installed-tension\n\
        strength,minimum_tensile_strength,pack=1\nloss,thermoelastic\n\
        key,layer,outer_radius_m,wire_diameter_m,pitch_m,material_index\n69,core,0.0005,,,0\n";
    let strings = Arc::new(StringMaterials::from_packs(&[
        include_bytes!("fixtures/packs/steel-music-wire-string-heptonstall-ligo-p1400229.fsmatpk"),
        include_bytes!("fixtures/packs/steel-music-wire-min-strength-qq-w-470b.fsmatpk"),
    ], construction, &operating_point("T", 293.15, "frequency", 1000.0).unwrap(),
        "bending_loss_tangent").unwrap());
    preparation::prepare_with_materials(&mut desc, None, None, Some(&strings)).unwrap();
    let mut cfg = Config::live(48_000);
    cfg.max_partials = 8;
    cfg.silent_bank_hz = 0.0;
    cfg.duplex = false;
    for substeps in [1, 2] {
        let mut a = PianoAudio::new_with_materials(&desc, cfg.clone(), substeps, gases[0], strings.clone()).unwrap();
        let mut b = PianoAudio::new_with_materials(&desc, cfg.clone(), substeps, gases[1], strings.clone()).unwrap();
        assert_eq!(a.piano().acoustic_gas(), gases[0]);
        assert_eq!(b.piano().acoustic_gas(), gases[1]);
        assert_ne!(a.piano().string_gas_drag_rates_per_s(69), b.piano().string_gas_drag_rates_per_s(69));
        assert!(b.microphone_delay_frames().0 > a.microphone_delay_frames().0);
        a.note_on(69, 0.25).unwrap(); b.note_on(69, 0.25).unwrap();
        let (mut peak, mut sound_change, mut motion_change) = (0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..4096 {
            let x = a.step_pressure().unwrap();
            let y = b.step_pressure().unwrap();
            assert!(x.is_finite() && y.is_finite());
            peak = peak.max(x.abs());
            sound_change = sound_change.max((x - y).abs());
            motion_change = motion_change.max((a.piano().bridge_velocity(69).unwrap()
                - b.piano().bridge_velocity(69).unwrap()).abs());
        }
        assert!(peak > 1e-9 && sound_change > 1e-8 * peak);
        assert!(motion_change > 1e-14, "source gas must affect mechanics, not only receiver gain");
        assert_ne!(a.close_ledger().string_loss_j, b.close_ledger().string_loss_j);
    }
}

#[test]
fn canonical_gas_ingress_preserves_source_domain_phase_and_packet_refusals() {
    use fs_live_piano::preparation::resolve_atmosphere_sources;
    let (bytes, lengths) = source_packet(0.028);
    let read = |record, data: &[u8], sizes: &[u32], water, phase|
        resolve_atmosphere_sources(record, data, sizes, water, phase);
    assert!(read(SOURCE_GAS, &bytes, &lengths, "", "frozen").is_ok());
    assert!(read(SOURCE_GAS, &bytes, &lengths, "", "equilibrium").is_err());
    assert!(read(SOURCE_GAS, &bytes, &lengths, "", "liquid-water").is_err());
    assert!(read(SOURCE_GAS, &bytes, &lengths[..2], "", "frozen").is_err());
    assert!(read(SOURCE_GAS, &bytes[..bytes.len() - 1], &lengths, "", "frozen").is_err());
    for record in [format!("{SOURCE_GAS} trailing"), SOURCE_GAS.replace("300", "199"),
        SOURCE_GAS.replace("piano-test-gas", "missing-source"), "293.15 100000 0".into()] {
        assert!(read(&record, &bytes, &lengths, "", "frozen").is_err(), "{record}");
    }
}

fn helium_like() -> GasState {
    GasState::try_new(
        &GasSpec {
            molar_mass: 4.002_602e-3,
            gamma: 5.0 / 3.0,
            sutherland_beta: 1.05e-6,
            sutherland_s: 79.4,
            conductivity: ConductivityModel::Eucken,
        },
        293.15,
        101_325.0,
    )
    .unwrap()
}

#[test]
fn arbitrary_admitted_gas_sets_the_receiver_without_dry_air_replacement() {
    let desc = Description::embedded().unwrap();
    let gas = helium_like();
    assert!(gas.sound_speed > 900.0 && gas.density < 0.3);
    let air = LivePiano::new(&desc, Config::live(48_000)).unwrap();
    let helium = LivePiano::new_with_gas(&desc, Config::live(48_000), gas).unwrap();
    // The same board-to-microphone paths (about 1 m), travelled about 2.9x faster.
    let (a, h) = (air.microphone_delay_samples().0, helium.microphone_delay_samples().0);
    assert!(a > 100 && 2 * h < a, "bulk travel delay: helium {h} samples, air {a}");
}

#[test]
fn independently_mutated_derived_gas_fields_refuse() {
    let desc = Description::embedded().unwrap();
    let mut gas = helium_like();
    gas.characteristic_impedance *= 1.1;
    assert!(LivePiano::new_with_gas(&desc, Config::live(48_000), gas).is_err());
    let mut gas = helium_like();
    gas.prandtl = f64::NAN;
    assert!(LivePiano::new_with_gas(&desc, Config::live(48_000), gas).is_err());
    // Repairing impedance alone must not bypass the shared ideal-gas closure.
    let mut gas = helium_like();
    gas.density *= 1.1;
    gas.characteristic_impedance = gas.density * gas.sound_speed;
    assert!(LivePiano::new_with_gas(&desc, Config::live(48_000), gas).is_err());
}

#[test]
fn atmosphere_text_reaches_real_pressure_at_both_audio_clocks() {
    use fs_live_piano::{PianoAudio, preparation::resolve_atmosphere};
    let mut desc = Description::embedded().unwrap().truncated(4).unwrap();
    desc.courses.retain(|c| c.midi == 69);
    let mut cfg = Config::live(48_000);
    cfg.max_partials = 8;
    cfg.silent_bank_hz = 0.0;
    cfg.duplex = false;
    let low = resolve_atmosphere("293.15 50000 0").unwrap();
    let high = resolve_atmosphere("293.15 100000 0").unwrap();
    let carbon = resolve_atmosphere("mixture 293.15 100000 mole 1 CO2 1").unwrap();
    assert_eq!(low.sound_speed, high.sound_speed);
    for substeps in [1, 2] {
        let mut a = PianoAudio::new_with_gas(&desc, cfg.clone(), substeps, low).unwrap();
        let mut b = PianoAudio::new_with_gas(&desc, cfg.clone(), substeps, high).unwrap();
        let c = PianoAudio::new_with_gas(&desc, cfg.clone(), substeps, carbon).unwrap();
        assert_eq!(a.piano().partials_hz(69), b.piano().partials_hz(69));
        assert_eq!(a.piano().partials_hz(69), c.piano().partials_hz(69));
        assert!(c.microphone_delay_frames().0 > a.microphone_delay_frames().0);
        a.note_on(69, 1.0).unwrap(); b.note_on(69, 1.0).unwrap();
        let mut peak = 0.0_f64;
        for _ in 0..1024 {
            let x = a.step_pressure().unwrap();
            let y = b.step_pressure().unwrap();
            peak = peak.max(x.abs());
            assert!((y - 2.0*x).abs() < 1e-12 * y.abs().max(1.0));
        }
        assert!(peak > 1e-7, "the comparison must exercise audible source pressure");
        assert_eq!(format!("{:?}", a.close_ledger()), format!("{:?}", b.close_ledger()));
        assert!(fs_live_piano::describe_audio(&b).contains("\"pressurePa\":100000"));
    }
}

#[test]
fn browser_atmosphere_admission_refuses_bad_and_incomplete_declarations() {
    use fs_live_piano::preparation::resolve_atmosphere;
    for text in ["", "293.15 0 0", "293.15 100000 0 trailing",
        "mixture 293.15 100000 mole 2 N2 0.8 O2 0.8",
        "mixture 293.15 100000 mole 1 unobtanium 1"] {
        assert!(resolve_atmosphere(text).is_err(), "{text}");
    }
}
