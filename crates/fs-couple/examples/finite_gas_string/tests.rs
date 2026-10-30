//! G1/G3/G4: source packs -> continuous electrical/thermal/gas state -> pressure.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_material::gas_species::input::{
    ThermochemicalPhase, ThermochemicalSources, read_isochoric_parcel,
};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::{Density, DynViscosity, Pressure, Temperature};

#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

fn scoped<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 0x464753,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn source(viscosity_scale: f64) -> MaterialCard {
    source_law(viscosity_scale, None)
}

fn source_law(viscosity_scale: f64, expansion_per_k: Option<f64>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut properties = vec![
        ("density", Density::DIMS, 1000.0, 1000.0, true),
        if expansion_per_k.is_some() {
            ("young_modulus", Pressure::DIMS, 2.02e9, 1.9e9, false)
        } else { ("young_modulus", Pressure::DIMS, 2e9, 2e9, true) },
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            if expansion_per_k.is_some() { 10.0 } else { 0.1 },
            if expansion_per_k.is_some() { 10.0 } else { 0.1 },
            true,
        ),
        (
            "kelvin_voigt_bending_viscosity",
            DynViscosity::DIMS,
            1e6 * viscosity_scale,
            2e7 * viscosity_scale,
            false,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            20.0,
            40.0,
            false,
        ),
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            8e-7,
            2e-6,
            false,
        ),
    ];
    if let Some(alpha) = expansion_per_k {
        properties.push(("linear_thermal_expansion_coefficient",
            fs_qty::Dims([0, 0, 0, -1, 0, 0]), alpha, alpha, true));
    }
    for (name, dims, low, high, constant) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if constant {
                    PropertyValue::Scalar { value: low, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(280.0, low), (400.0, high)],
                        dims,
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 400.0)
                    .with("omega", 0.0, 200000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if constant {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured live-string integration fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let pack = NormalizedPack::new(
        "manufactured-live-string",
        "integration-v1",
        fs_blake3::ContentHash([17; 32]),
        "CC0-1.0",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    // Exercise the same canonical source bytes and material association as the CLI.
    let decoded = NormalizedPack::from_bytes(&pack.to_bytes()).unwrap();
    material(&decoded).unwrap()
}

fn parcel(name: &str, molar_mass: f64, cp_r: f64, volume: f64) -> IsochoricGasParcel {
    let sources = gas_support::packs(name, molar_mass, cp_r, 200.0);
    let bindings = [sources.binding()];
    let registry = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let record = format!("thermochemical 300 85000 mole 1 {name} 1");
    read_isochoric_parcel(&mut record.split_ascii_whitespace(), &registry, volume).unwrap()
}

fn scenario() -> input::Spec {
    let mut spec = input::parse(include_str!("heated.string")).unwrap();
    spec.samples = 128;
    spec.modes = 2;
    spec.validate().unwrap();
    spec
}

#[test]
fn strict_input_retains_signed_drives_and_bounds_physical_work() {
    let text = include_str!("heated.string");
    let signed = input::parse(
        &text
            .replace("drive current 0.1", "drive voltage -0.25")
            .replace("release 0.4 0.001 16", "release 0.4 -0.001 16")
            .replace("heat 0 0", "heat -0.001 0.002"),
    )
    .unwrap();
    assert_eq!(
        signed.drive,
        Some(fs_material::conductor::OhmicDrive::Voltage(-0.25))
    );
    assert_eq!(signed.force_n, -0.001);
    assert_eq!(signed.solid_power_w, -0.001);
    assert_eq!(signed.gas_power_w, 0.002);
    assert!(
        input::parse(&text.replace("drive current 0.1", "drive none"))
            .unwrap()
            .drive
            .is_none()
    );
    for invalid in [
        String::new(),
        format!("{text}\nheat 0 0"),
        text.replace(
            "frankensim-finite-gas-string-v1",
            "frankensim-finite-gas-string-v2",
        ),
        text.replace("audio 16000 3200", "audio 0 3200"),
        text.replace("audio 16000 3200", "audio 16000 262145"),
        text.replace("20 4", "20 33"),
        text.replace("release 0.4 0.001 16", "release 0.4 0.001 3201"),
        text.replace("gas 0.000001", "gas 0"),
        text.replace("conduction 0.001", "conduction 1"),
        text.replace("drive current 0.1", "drive current NaN"),
        text.replace("heat 0 0", "heat inf 0"),
        " ".repeat(input::MAX_INPUT_BYTES + 1),
    ] {
        assert!(
            input::parse(&invalid).is_err(),
            "invalid scenario was accepted"
        );
    }
}

#[test]
fn continuous_pressure_replays_and_matches_independent_three_storage_energy() {
    scoped(|cx, _| {
        let spec = scenario();
        let card = source(1.0);
        let gas = parcel("surroundings", 0.028, 3.5, spec.gas_volume_m3);
        let initial = build(cx, &spec, &card, &gas).unwrap();
        assert!(
            initial
                .mechanical()
                .states()
                .iter()
                .all(|state| state.velocity_m_sqrt_kg_per_s == 0.0)
        );
        assert!(
            initial
                .mechanical()
                .states()
                .iter()
                .any(|state| state.displacement_m_sqrt_kg != 0.0)
        );
        let a = simulate(cx, &spec, &card, gas.clone(), None).unwrap();
        let b = simulate(cx, &spec, &card, gas, None).unwrap();
        assert_eq!(a.pressure_pa, b.pressure_pa);
        assert_eq!(a.records, b.records);
        assert_eq!(a.records.len(), spec.samples);
        assert_eq!(a.pressure_pa.len(), spec.samples);
        for (index, row) in a.records.iter().enumerate() {
            assert!(row.iter().all(|value| value.is_finite()));
            assert_eq!(row[1], a.pressure_pa[index]);
            assert!((row[0] - (index + 1) as f64 / f64::from(spec.sample_rate_hz)).abs() < 1e-15);
            assert!(row[15].abs() <= row[16]);
        }
        assert!(
            a.records[spec.release_sample..]
                .iter()
                .any(|row| row[2].abs() > 1e-10)
        );
        assert!(
            a.pressure_pa[spec.release_sample..]
                .iter()
                .any(|value| value.abs() > 1e-12)
        );

        // With one sine mode, a single physical station reconstructs its entire
        // energy independently of the driver's modal-energy/report functions.
        let mut one = spec.clone();
        one.modes = 1;
        one.force_n = 0.01;
        let history = simulate(
            cx,
            &one,
            &card,
            parcel("surroundings", 0.028, 3.5, one.gas_volume_m3),
            None,
        )
        .unwrap();
        let area = std::f64::consts::PI * one.radius_m.powi(2);
        let second_moment = std::f64::consts::PI * one.radius_m.powi(4) / 4.0;
        let solid_mass = 1000.0 * area * one.length_m;
        let wavenumber = std::f64::consts::PI / one.length_m;
        let input::Prestress::FixedTension(tension) = one.prestress else { panic!("fixed tension fixture") };
        let omega_squared = (tension * wavenumber.powi(2)
            + 2e9 * second_moment * wavenumber.powi(4))
            / (1000.0 * area);
        let phi_squared = 2.0 * (std::f64::consts::PI * one.station).sin().powi(2) / solid_mass;
        let initial_mechanical = 0.5 * one.force_n.powi(2) * phi_squared / omega_squared;
        let universal_r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let gas_mass = 85000.0 * one.gas_volume_m3 * 0.028 / (universal_r * 300.0);
        let dt = 1.0 / f64::from(one.sample_rate_hz);
        let mut previous_joule = 0.0;
        for row in &history.records {
            let mechanical = 0.5 * (row[3].powi(2) + omega_squared * row[2].powi(2)) / phi_squared;
            let solid = solid_mass * 0.1 * (row[4] - 300.0);
            // Fixture NASA law has cp/R = 3.5 + 0.001 T, so cv/R = 2.5 + 0.001 T.
            let gas = gas_mass * universal_r / 0.028
                * (2.5 * (row[5] - 300.0) + 0.0005 * (row[5].powi(2) - 300.0_f64.powi(2)));
            let stored = mechanical - initial_mechanical + solid + gas;
            assert!(
                (row[14] - stored).abs() < 5e-12,
                "independent stored energy {stored}, row {}",
                row[14]
            );
            let supplied = row[10] + row[11] + row[12] + row[13];
            assert!((stored - supplied).abs() <= row[16] + 5e-12);
            assert!((row[10] - previous_joule - row[17].powi(2) * row[19] * dt).abs() < 1e-14);
            previous_joule = row[10];
            let resistivity = 8e-7 + (2e-6 - 8e-7) * (row[4] - 280.0) / 120.0;
            assert!((row[20] / (resistivity * one.length_m / area) - 1.0).abs() < 1e-12);
            assert!((row[6] / (85000.0 * row[5] / 300.0) - 1.0).abs() < 1e-12);
        }
    });
}

#[test]
fn sourced_expansion_changes_live_tension_pitch_and_pressure_through_the_same_input() {
    let text = include_str!("expanding.string");
    assert!(input::parse(&text.replace("fixed_support", "20")).is_err());
    let missing = text.split("thermoelastic").next().unwrap();
    assert!(input::parse(missing).is_err());
    scoped(|cx, _| {
        let mut spec = input::parse(text).unwrap();
        spec.samples = 128;
        spec.modes = 2;
        let thermal = spec.thermoelastic.unwrap();
        let card = source_law(1.0, Some(2e-5));
        let gas = parcel("surroundings", 0.028, 3.5, spec.gas_volume_m3);
        let history = simulate(cx, &spec, &card, gas.clone(), None).unwrap();
        let without_expansion = simulate(cx, &spec, &source_law(1.0, Some(0.0)), gas, None).unwrap();
        assert_ne!(history.pressure_pa, without_expansion.pressure_pa);
        let area = std::f64::consts::PI * spec.radius_m.powi(2);
        let second_moment = area * spec.radius_m.powi(2) / 4.0;
        let k = std::f64::consts::PI / spec.length_m;
        for row in &history.records {
            let young = 2e9 - 1e6 * (row[4] - 300.0);
            let strain = (spec.length_m - thermal.reference_length_m) / thermal.reference_length_m
                - 2e-5 * (row[4] - thermal.reference_temperature_k);
            let tension = young * area * strain;
            let frequency = ((tension * k.powi(2) + young * second_moment * k.powi(4))
                / (1000.0 * area)).sqrt() / (2.0 * std::f64::consts::PI);
            assert!((row[22] - tension).abs() < 1e-12);
            assert!((row[23] - frequency).abs() < 1e-10);
            assert!(row[15].abs() <= row[16]);
        }
        let first = history.records.first().unwrap();
        let last = history.records.last().unwrap();
        assert!(last[4] > first[4]);
        assert!(last[22] < first[22]);
        assert!(last[23] < first[23]);
        assert!(last[24].abs() > 0.0);
        assert!(last[25].abs() > 0.0);
        assert!(last[26].abs() > 0.0);
    });
}

#[test]
fn material_species_inventory_and_drive_changes_reach_the_same_sound() {
    scoped(|cx, _| {
        let spec = scenario();
        let card = source(1.0);
        let baseline = simulate(
            cx,
            &spec,
            &card,
            parcel("surroundings", 0.028, 3.5, spec.gas_volume_m3),
            None,
        )
        .unwrap();
        let changed_card = source(2.0);
        assert_ne!(card.content_hash(), changed_card.content_hash());
        let viscosity = simulate(
            cx,
            &spec,
            &changed_card,
            parcel("surroundings", 0.028, 3.5, spec.gas_volume_m3),
            None,
        )
        .unwrap();
        let species = simulate(
            cx,
            &spec,
            &card,
            parcel("light-source", 0.004, 2.5, spec.gas_volume_m3),
            None,
        )
        .unwrap();
        let mut large = spec.clone();
        large.gas_volume_m3 *= 100.0;
        let inventory = simulate(
            cx,
            &large,
            &card,
            parcel("surroundings", 0.028, 3.5, large.gas_volume_m3),
            None,
        )
        .unwrap();
        let mut no_drive = spec.clone();
        no_drive.drive = None;
        let unpowered = simulate(
            cx,
            &no_drive,
            &card,
            parcel("surroundings", 0.028, 3.5, no_drive.gas_volume_m3),
            None,
        )
        .unwrap();
        for different in [&viscosity, &species, &inventory, &unpowered] {
            assert_ne!(baseline.pressure_pa, different.pressure_pa);
            assert_ne!(
                baseline.records.last().unwrap()[2],
                different.records.last().unwrap()[2]
            );
        }
        let last = baseline.records.last().unwrap();
        assert!(last[10] > 0.0);
        assert!(unpowered.records.iter().all(|row| row[10] == 0.0));
        assert!(last[4] > unpowered.records.last().unwrap()[4]);
        assert!(last[5] > inventory.records.last().unwrap()[5]);
        assert!(viscosity.records.last().unwrap()[21] > last[21]);
    });
}

#[test]
fn source_exit_excess_heat_and_cancellation_return_no_partial_history() {
    scoped(|cx, gate| {
        let spec = scenario();
        let card = source(1.0);
        let gas = parcel("surroundings", 0.028, 3.5, spec.gas_volume_m3);
        let mut unsupported = spec.clone();
        unsupported.solid_maximum_k = 450.0;
        assert!(build(cx, &unsupported, &card, &gas).is_err());
        let mut overheat = spec.clone();
        overheat.solid_power_w = 1e12;
        assert!(simulate(cx, &overheat, &card, gas.clone(), None).is_err());
        let mut excluded_gas = spec.clone();
        excluded_gas.gas_solve.minimum_temperature_k = 310.0;
        assert!(simulate(cx, &excluded_gas, &card, gas.clone(), None).is_err());
        gate.request();
        assert!(build(cx, &spec, &card, &gas).is_err());
        assert!(simulate(cx, &spec, &card, gas, None).is_err());
    });
}

#[test]
fn bonded_input_preserves_layer_order_and_rejects_ignored_geometry() {
    let text = include_str!("bonded.string");
    let spec = input::parse(text).unwrap();
    let bonded = spec.bonded.as_ref().unwrap();
    assert!(spec.thermoelastic.is_none());
    assert_eq!(bonded.sleeves[0].pack, "sleeve.pack");
    assert_eq!(bonded.sleeves[0].outer_radius_m, spec.radius_m);
    for bad in [
        text.replace("fixed_support", "20"),
        text.replace("bonded 0.00015 1", "bonded 0.00015 0"),
        text.replace("bonded 0.00015 1", "bonded 0.00015 16"),
        text.replace("sleeve.pack", "../sleeve.pack"),
        text.replace("sleeve.pack", "/sleeve.pack"),
        text.replace("layer sleeve.pack 0.00025", "layer sleeve.pack 0.00024"),
        text.replace("layered_solve 4096", "layered_solve 0"),
        format!("{text}\nthermoelastic 0.648 300 0.01 0.00001 0.000000000001 64"),
    ] { assert!(input::parse(&bad).is_err()); }
}

#[test]
fn bonded_cards_drive_layer_heat_parallel_current_pitch_and_continuing_pressure() {
    scoped(|cx, gate| {
        let mut spec = input::parse(include_str!("bonded.string")).unwrap();
        spec.samples = 32;
        spec.modes = 1;
        let cards = [source_law(1.0, Some(2e-5)), source_law(2.0, Some(4e-5))];
        let gas = parcel("surroundings", 0.028, 3.5, spec.gas_volume_m3);
        let a = layered::simulate(cx, &spec, &cards, gas.clone(), None).unwrap();
        let b = layered::simulate(cx, &spec, &cards, gas.clone(), None).unwrap();
        assert_eq!(a.records, b.records);
        assert_eq!(a.pressure_pa, b.pressure_pa);
        assert_eq!(a.layer_temperatures, b.layer_temperatures);
        assert_eq!(a.layer_temperatures.len(), spec.samples);
        let bonded = spec.bonded.as_ref().unwrap();
        let radii = [bonded.core_radius_m, spec.radius_m];
        let pi = std::f64::consts::PI;
        let k = pi/spec.length_m;
        let mut previous_joule = 0.0;
        for (row, temperatures) in a.records.iter().zip(&a.layer_temperatures) {
            assert_eq!(temperatures.len(), 2);
            let (mut tension, mut bending, mut conductance) = (0.0, 0.0, 0.0);
            let mut inner = 0.0_f64;
            for (i, (&outer, &t)) in radii.iter().zip(temperatures).enumerate() {
                let area = pi*(outer.powi(2)-inner.powi(2));
                let moment = pi/4.0*(outer.powi(4)-inner.powi(4));
                let young = 2e9-1e6*(t-300.0);
                let strain = (spec.length_m-bonded.reference_length_m)/bonded.reference_length_m
                    -[2e-5, 4e-5][i]*(t-bonded.reference_temperature_k);
                tension += young*area*strain;
                bending += young*moment;
                let resistivity = 8e-7+(2e-6-8e-7)*(t-280.0)/120.0;
                conductance += area/(spec.length_m*resistivity);
                inner = outer;
            }
            let frequency = ((tension*k*k+bending*k.powi(4))/(1000.0*pi*spec.radius_m.powi(2)))
                .sqrt()/(2.0*pi);
            assert!((row[22]-tension).abs() < 1e-11);
            assert!((row[23]-frequency).abs() < 1e-10);
            assert!((row[20]*conductance-1.0).abs() < 1e-12);
            let joule = row[18].powi(2)*row[20]/f64::from(spec.sample_rate_hz);
            assert!((row[10]-previous_joule-joule).abs() < 1e-13);
            previous_joule = row[10];
            assert!(row[15].abs() <= row[16]);
        }
        assert!(a.layer_temperatures.iter().any(|t| t[0] != t[1]));
        assert!(a.pressure_pa[spec.release_sample..].iter().any(|p| p.abs() > 1e-12));
        let swapped = [cards[0].clone(), source_law(4.0, Some(0.0))];
        let other = layered::simulate(cx, &spec, &swapped, gas.clone(), None).unwrap();
        assert_ne!(a.pressure_pa, other.pressure_pa);
        assert_ne!(a.records.last().unwrap()[22], other.records.last().unwrap()[22]);
        let other_gas = layered::simulate(cx, &spec, &cards,
            parcel("light-source", 0.004, 2.5, spec.gas_volume_m3), None).unwrap();
        assert_ne!(a.pressure_pa, other_gas.pressure_pa);
        let mut cold = spec.clone();
        cold.drive = None;
        cold.solid_power_w = 0.0;
        let unpowered = layered::simulate(cx, &cold, &cards, gas.clone(), None).unwrap();
        assert!(unpowered.records.iter().all(|row| row[10] == 0.0));
        assert!(a.records.last().unwrap()[4] > unpowered.records.last().unwrap()[4]);
        assert!(layered::simulate(cx, &spec, &cards[..1], gas.clone(), None).is_err());
        assert!(layered::simulate(cx, &spec, &[source(1.0), cards[1].clone()], gas.clone(), None).is_err());
        let mut excess = spec.clone();
        excess.solid_power_w = 1e12;
        assert!(layered::simulate(cx, &excess, &cards, gas.clone(), None).is_err());
        gate.request();
        assert!(layered::simulate(cx, &spec, &cards, gas, None).is_err());
    });
}

mod flow;
mod electrical_feedback;
mod attached_sources;
