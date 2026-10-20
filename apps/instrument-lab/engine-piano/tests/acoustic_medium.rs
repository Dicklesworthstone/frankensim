//! The Rayleigh receiver takes the caller's admitted gas as given (no dry-air
//! substitute), and a gas whose derived fields were mutated independently is
//! refused before anything is built.
use fs_live_piano::{Config, Description, LivePiano};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};

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
