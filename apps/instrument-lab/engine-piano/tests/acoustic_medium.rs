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
}
