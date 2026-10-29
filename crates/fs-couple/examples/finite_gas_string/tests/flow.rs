//! Source bytes and nozzle data reach uniform and layered physical owners.
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

fn uniform_scenario(thermoelastic: bool) -> input::Spec {
    let text = if thermoelastic {
        include_str!("../expanding.string")
    } else {
        include_str!("../heated.string")
    };
    let text = text.replace("\nconduction ", &format!("\n{NOZZLE}\nconduction "));
    let mut spec = input::parse(&text).unwrap();
    spec.samples = 32;
    spec.modes = 1;
    if thermoelastic {
        spec.drive = Some(fs_material::conductor::OhmicDrive::Voltage(-0.25));
    }
    spec
}

fn simulate_scenario(
    cx: &Cx<'_>,
    spec: &input::Spec,
    gas: IsochoricGasParcel,
    reservoir: Option<ClosedGasParcel>,
) -> Result<History> {
    if spec.bonded.is_some() {
        let cards = [source_law(1.0, Some(2e-5)), source_law(2.0, Some(4e-5))];
        layered::simulate(cx, spec, &cards, gas, reservoir)
    } else {
        let card = if spec.thermoelastic.is_some() {
            source_law(1.0, Some(2e-5))
        } else {
            source(1.0)
        };
        simulate(cx, spec, &card, gas, reservoir)
    }
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
fn nozzle_data_is_explicit_bounded_and_shared_by_all_supported_strings() {
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
    ] {
        assert!(input::parse(&bad).is_err());
    }
    for uniform in [uniform_scenario(false), uniform_scenario(true)] {
        assert!(uniform.bonded.is_none());
        assert_eq!(uniform.nozzle.unwrap().geometry, nozzle.geometry);
    }
    let supplied = input::parse(include_str!("../heated-flow.string")).unwrap();
    assert!(supplied.bonded.is_none() && supplied.nozzle.is_some());
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
    for spec in [scenario(), uniform_scenario(false), uniform_scenario(true)] {
        scoped(|cx, _| {
            let (gas, reservoir) = inventories(&spec, FEED).unwrap();
            let reservoir = reservoir.unwrap();
            let masses = [gas.mass_kg(), reservoir.mass_kg()];
            let initial_energy = [
                energy(masses[0], 300.0, 1.0, 0.0),
                energy(masses[1], 300.0, 0.0, 1.0),
            ];
            let a = simulate_scenario(cx, &spec, gas.clone(), Some(reservoir.clone())).unwrap();
            let b = simulate_scenario(cx, &spec, gas, Some(reservoir)).unwrap();
            assert_eq!(a.pressure_pa, b.pressure_pa);
            assert_eq!(a.records, b.records);
            let flow = a.gas_flow.as_ref().unwrap();
            assert_eq!(flow.records, b.gas_flow.as_ref().unwrap().records);
            assert_eq!(flow.species, ["heavy", "light"]);
            let last = a.records.last().unwrap();
            assert!(
                last[10] > 0.0,
                "electrical source must remain active during flow"
            );
            assert!(
                last[9].abs() > 1e-12,
                "reciprocal conduction must remain active during flow"
            );
            if matches!(
                spec.drive,
                Some(fs_material::conductor::OhmicDrive::Voltage(_))
            ) {
                assert!(a.records.iter().all(|row| row[17] < 0.0 && row[18] < 0.0));
            }
            let mut previous_joule = 0.0;
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
                if let input::Prestress::FixedTension(tension) = spec.prestress {
                    // The uniform one-mode fixture also admits an independent
                    // complete storage calculation, including BOTH changing masses.
                    let area = std::f64::consts::PI * spec.radius_m.powi(2);
                    let solid_mass = 1000.0 * area * spec.length_m;
                    let k = std::f64::consts::PI / spec.length_m;
                    let moment = area * spec.radius_m.powi(2) / 4.0;
                    let omega2 = (tension * k.powi(2) + 2e9 * moment * k.powi(4)) / (1000.0 * area);
                    let phi2 =
                        2.0 * (std::f64::consts::PI * spec.station).sin().powi(2) / solid_mass;
                    let initial_mechanical = 0.5 * spec.force_n.powi(2) * phi2 / omega2;
                    let mechanical = 0.5 * (row[3].powi(2) + omega2 * row[2].powi(2)) / phi2;
                    let expected = mechanical - initial_mechanical
                        + solid_mass * 0.1 * (row[4] - spec.solid_temperature_k)
                        + stored[0]
                        - initial_energy[0]
                        + stored[1]
                        - initial_energy[1];
                    assert!((row[14] - expected).abs() < 5e-11);
                }
                if spec.bonded.is_none() {
                    // Uniform electrical quadrature uses its original initial law;
                    // the transported gas enthalpy must not become Joule work.
                    let expected = row[17].powi(2) * row[19] / f64::from(spec.sample_rate_hz);
                    assert!((row[10] - previous_joule - expected).abs() < 1e-14);
                    previous_joule = row[10];
                }
            }
            let (gas, reservoir) =
                inventories(&spec, "thermochemical 300 95000 mole 2 heavy 0.5 light 0.5").unwrap();
            let changed = simulate_scenario(cx, &spec, gas, reservoir).unwrap();
            assert_eq!(a.records[0][2], changed.records[0][2]); // same first mechanical trial
            assert_ne!(a.pressure_pa, changed.pressure_pa);
            assert_ne!(
                a.records.last().unwrap()[21],
                changed.records.last().unwrap()[21]
            );
            assert!(flow.records.last().unwrap()[14] > flow.records.first().unwrap()[14]);
            let mut unpowered = spec.clone();
            unpowered.drive = None;
            let (gas, reservoir) = inventories(&unpowered, FEED).unwrap();
            let unpowered = simulate_scenario(cx, &unpowered, gas, reservoir).unwrap();
            assert!(unpowered.records.iter().all(|row| row[10] == 0.0));
            assert!(
                unpowered
                    .gas_flow
                    .as_ref()
                    .unwrap()
                    .records
                    .iter()
                    .all(|row| row[7] > 0.0)
            );
            assert_ne!(unpowered.pressure_pa, a.pressure_pa);
        });
    }
}

#[test]
fn uniform_flow_refuses_invalid_reservoir_sources_and_cancelled_runs() {
    scoped(|cx, gate| {
        let spec = uniform_scenario(false);
        let (gas, reservoir) =
            inventories(&spec, "thermochemical 700 95000 mole 1 light 1").unwrap();
        assert!(simulate_scenario(cx, &spec, gas, reservoir).is_err());
        let (gas, reservoir) = inventories(&spec, FEED).unwrap();
        let mut withdrawal = spec.clone();
        withdrawal.nozzle.as_mut().unwrap().geometry.area_m2 = 1.0;
        assert!(simulate_scenario(cx, &withdrawal, gas.clone(), reservoir.clone()).is_err());
        let mut heating = spec.clone();
        heating.solid_power_w = 1e12;
        assert!(simulate_scenario(cx, &heating, gas.clone(), reservoir.clone()).is_err());
        gate.request();
        assert!(simulate_scenario(cx, &spec, gas, reservoir).is_err());
    });
}

#[test]
fn closed_path_preserves_heat_only_output_and_pressure_reversal_uses_same_native_route() {
    for spec in [scenario(), uniform_scenario(false), uniform_scenario(true)] {
        scoped(|cx, _| {
            let mut closed = spec.clone();
            closed.nozzle.as_mut().unwrap().geometry.area_m2 = 0.0;
            let (gas, reservoir) = inventories(&closed, FEED).unwrap();
            let mut heat_only = closed.clone();
            heat_only.nozzle = None;
            let original = simulate_scenario(cx, &heat_only, gas.clone(), None).unwrap();
            let sealed = simulate_scenario(cx, &closed, gas, reservoir).unwrap();
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
            let (gas, reservoir) =
                inventories(&spec, "thermochemical 300 75000 mole 1 light 1").unwrap();
            assert!(simulate_scenario(cx, &spec, gas.clone(), None).is_err());
            assert!(simulate_scenario(cx, &heat_only, gas.clone(), reservoir.clone()).is_err());
            let mut wrong_volume = spec.clone();
            wrong_volume.nozzle.as_mut().unwrap().reservoir_volume_m3 *= 2.0;
            assert!(simulate_scenario(cx, &wrong_volume, gas.clone(), reservoir.clone()).is_err());
            let outward = simulate_scenario(cx, &spec, gas, reservoir).unwrap();
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
}
