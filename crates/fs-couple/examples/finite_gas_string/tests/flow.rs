//! Source bytes and authored nozzle data reach the original layered owners.
use super::*;
use std::collections::BTreeMap;
use std::path::PathBuf;

const NOZZLE: &str = "nozzle reservoir.gas 0.000002 0.00000001 1 0.00000000000001 0.000000001 0.05";
const FEED: &str = "thermochemical 300 95000 mole 1 light 1";

fn scenario() -> input::Spec {
    let text = include_str!("../bonded.string")
        .replace("\nconduction ", &format!("\n{NOZZLE}\nconduction "));
    let mut spec = input::parse(&text).unwrap();
    spec.samples = 32;
    spec.modes = 1;
    spec
}

fn inventories(
    spec: &input::Spec,
    feed: &str,
) -> Result<(IsochoricGasParcel, Option<ClosedGasParcel>)> {
    let mut args = [
        "--gas-record",
        "local.gas",
        "--gas-water",
        "none",
        "--gas-phase",
        "frozen",
    ]
    .map(str::to_owned)
    .to_vec();
    let mut files = BTreeMap::<PathBuf, Vec<u8>>::new();
    for (name, mass, cp_r) in [("heavy", 0.028, 3.5), ("light", 0.004, 2.5)] {
        let packs = gas_support::packs(name, mass, cp_r, 200.0);
        let paths = [
            format!("{name}.species"),
            format!("{name}.nasa9"),
            format!("{name}.transport"),
        ];
        args.push("--gas-component".into());
        args.extend(paths.iter().cloned());
        for (path, bytes) in paths.into_iter().zip([
            packs.species.to_bytes(),
            packs.nasa9.to_bytes(),
            packs.transport.to_bytes(),
        ]) {
            files.insert(PathBuf::from(path), bytes);
        }
    }
    files.insert(
        PathBuf::from("local.gas"),
        b"thermochemical 300 85000 mole 1 heavy 1".to_vec(),
    );
    files.insert(
        PathBuf::from("scenario/reservoir.gas"),
        feed.as_bytes().to_vec(),
    );
    let (remaining, options) = thermochemical::Options::extract(&args).unwrap();
    assert!(remaining.is_empty());
    // Production source loader, canonical pack decoders and original GasRecord
    // parser, with only file I/O replaced by an immutable attachment map.
    load_gas(&options, spec, Path::new("scenario"), |path, limit| {
        let bytes = files
            .get(path)
            .ok_or_else(|| format!("missing fixture {}", path.display()))?;
        assert!(bytes.len() <= limit);
        Ok(bytes.clone())
    })
}

#[test]
fn nozzle_data_is_explicit_bounded_and_never_ignored_by_uniform_modes() {
    let text = include_str!("../bonded.string")
        .replace("\nconduction ", &format!("\n{NOZZLE}\nconduction "));
    let spec = input::parse(&text).unwrap();
    let nozzle = spec.nozzle.as_ref().unwrap();
    assert_eq!(nozzle.geometry.discharge_coefficient, 1.0);
    assert_eq!(
        nozzle.controls(spec.gas_solve).transfer.donor,
        spec.gas_solve
    );
    for bad in [
        text.replace("reservoir.gas", "../reservoir.gas"),
        text.replace("reservoir.gas", "/reservoir.gas"),
        text.replace("reservoir.gas 0.000002", "reservoir.gas 0"),
        text.replace("0.00000001 1", "-0.00000001 1"),
        text.replace("0.00000001 1", "0.00000001 NaN"),
        text.replace("0.00000001 1", "0.00000001 1.01"),
        text.replace("0.00000000000001 0.000000001 0.05", "0 0.000000001 0.05"),
        text.replace("0.000000001 0.05", "0.000000001 1"),
        include_str!("../heated.string")
            .replace("\nconduction ", &format!("\n{NOZZLE}\nconduction ")),
        include_str!("../expanding.string")
            .replace("\nconduction ", &format!("\n{NOZZLE}\nconduction ")),
    ] {
        assert!(input::parse(&bad).is_err());
    }
    assert!(inventories(&spec, "thermochemical 300 95000 mole 1 missing 1").is_err());
    assert!(inventories(&spec, &format!("{FEED} trailing")).is_err());
}

// Independent fixture NASA integral: cp/R = cp0/R + .001 T, zero reference
// integration constant. Composition changes molecular mass AND specific u.
fn energy(mass: f64, temperature: f64, heavy: f64, light: f64) -> f64 {
    let molar_mass = 0.028 * heavy + 0.004 * light;
    mass * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass
        * ((2.5 * heavy + 1.5 * light) * temperature + 0.0005 * temperature.powi(2))
}

#[test]
fn source_loaded_species_exchange_changes_pressure_and_closes_extensive_inventory_energy() {
    scoped(|cx, _| {
        let spec = scenario();
        let cards = [source_law(1.0, Some(2e-5)), source_law(2.0, Some(4e-5))];
        let (gas, reservoir) = inventories(&spec, FEED).unwrap();
        let reservoir = reservoir.unwrap();
        let masses = [gas.mass_kg(), reservoir.mass_kg()];
        let initial_energy = [
            energy(masses[0], 300.0, 1.0, 0.0),
            energy(masses[1], 300.0, 0.0, 1.0),
        ];
        let a = layered::simulate(cx, &spec, &cards, gas.clone(), Some(reservoir.clone())).unwrap();
        let b = layered::simulate(cx, &spec, &cards, gas, Some(reservoir)).unwrap();
        assert_eq!(a.pressure_pa, b.pressure_pa);
        assert_eq!(a.records, b.records);
        let flow = a.gas_flow.as_ref().unwrap();
        assert_eq!(flow.records, b.gas_flow.as_ref().unwrap().records);
        assert_eq!(flow.species, ["heavy", "light"]);
        for (row, flux) in a.records.iter().zip(&flow.records) {
            assert!(flux[7] > 0.0 && flux[8] > 0.0);
            assert!((flux[0] + flux[2] - masses.iter().sum::<f64>()).abs() < 1e-14);
            let molecular = [
                0.028 * flux[12] + 0.004 * flux[14],
                0.028 * flux[13] + 0.004 * flux[15],
            ];
            for (species, molar_mass) in [0.028, 0.004].into_iter().enumerate() {
                let inventory = flux[0] * flux[12 + 2 * species] * molar_mass / molecular[0]
                    + flux[2] * flux[13 + 2 * species] * molar_mass / molecular[1];
                assert!((inventory - masses[species]).abs() < 1e-14);
            }
            let stored = [
                energy(flux[0], row[5], flux[12], flux[14]),
                energy(flux[2], flux[4], flux[13], flux[15]),
            ];
            assert!((flux[1] - stored[0]).abs() < 1e-11);
            assert!((flux[3] - stored[1]).abs() < 1e-11);
            assert!((flux[6] - stored.iter().sum::<f64>()).abs() < 2e-11);
            assert!((flux[10] - (stored[0] - initial_energy[0])).abs() < 2e-11);
            assert!((flux[11] - (stored[1] - initial_energy[1])).abs() < 2e-11);
            assert!((flux[0] - masses[0] - flux[8]).abs() < 1e-14);
            assert!((flux[11] + flux[9]).abs() < 2e-9);
            assert!(row[15].abs() <= row[16]);
        }
        let (gas, reservoir) =
            inventories(&spec, "thermochemical 300 95000 mole 2 heavy 0.5 light 0.5").unwrap();
        let changed = layered::simulate(cx, &spec, &cards, gas, reservoir).unwrap();
        assert_eq!(a.records[0][2], changed.records[0][2]); // same first mechanical trial
        assert_ne!(a.pressure_pa, changed.pressure_pa);
        assert_ne!(
            a.records.last().unwrap()[21],
            changed.records.last().unwrap()[21]
        );
        assert!(flow.records.last().unwrap()[14] > flow.records.first().unwrap()[14]);
    });
}

#[test]
fn closed_path_preserves_heat_only_output_and_pressure_reversal_uses_same_native_route() {
    scoped(|cx, _| {
        let mut closed = scenario();
        closed.nozzle.as_mut().unwrap().geometry.area_m2 = 0.0;
        let cards = [source_law(1.0, Some(2e-5)), source_law(2.0, Some(4e-5))];
        let (gas, reservoir) = inventories(&closed, FEED).unwrap();
        let mut heat_only = closed.clone();
        heat_only.nozzle = None;
        let original = layered::simulate(cx, &heat_only, &cards, gas.clone(), None).unwrap();
        let sealed = layered::simulate(cx, &closed, &cards, gas, reservoir).unwrap();
        assert_eq!(original.pressure_pa, sealed.pressure_pa);
        assert_eq!(original.layer_temperatures, sealed.layer_temperatures);
        for (original, sealed) in original.records.iter().zip(&sealed.records) {
            for index in 0..27 {
                if index != 16 {
                    // composition also admits its explicit transfer tolerance
                    assert_eq!(original[index], sealed[index]);
                }
            }
        }
        assert!(
            sealed
                .gas_flow
                .as_ref()
                .unwrap()
                .records
                .iter()
                .all(|row| row[7] == 0.0 && row[8] == 0.0 && row[11] == 0.0)
        );
        let spec = scenario();
        let (gas, reservoir) =
            inventories(&spec, "thermochemical 300 75000 mole 1 light 1").unwrap();
        assert!(layered::simulate(cx, &spec, &cards, gas.clone(), None).is_err());
        let mut wrong_volume = spec.clone();
        wrong_volume.nozzle.as_mut().unwrap().reservoir_volume_m3 *= 2.0;
        assert!(
            layered::simulate(cx, &wrong_volume, &cards, gas.clone(), reservoir.clone()).is_err()
        );
        let outward = layered::simulate(cx, &spec, &cards, gas, reservoir).unwrap();
        assert!(
            outward
                .gas_flow
                .as_ref()
                .unwrap()
                .records
                .iter()
                .all(|row| row[7] < 0.0)
        );
        assert!(outward.records.iter().all(|row| row[15].abs() <= row[16]));
    });
}
