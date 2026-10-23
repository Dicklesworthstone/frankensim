//! G1/G3/G4: actual source ingress -> finite reservoir -> inertial feedback.
use super::*;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};

#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

fn scoped<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 0x544d44,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn parcel(volume: f64, cp_r: f64) -> IsochoricGasParcel {
    let packs = support::packs("surroundings", 0.028, cp_r, 200.0);
    let bindings = [packs.binding()];
    let registry = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    read_isochoric_parcel(
        &mut "thermochemical 300 100000 mole 1 surroundings 1".split_ascii_whitespace(),
        &registry,
        volume,
    )
    .unwrap()
}

fn reservoir(volume: f64) -> Reservoir {
    Reservoir {
        volume_m3: volume,
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 600.0,
        maximum_exchange_j: 1e-5,
        power_w: 0.0,
        solid_temperature_k: 300.0,
    }
}
fn motion() -> Motion {
    Motion {
        velocity: 0.5,
        payload: 0.0,
        duration: 0.01,
        count: 200,
        film: 100.0,
    }
}
fn strings(words: &[&str]) -> Vec<String> {
    words.iter().map(|word| (*word).to_owned()).collect()
}

#[test]
fn finite_options_preserve_the_original_command_and_require_explicit_sources() {
    let old = strings(&["--table", "material.txt", "0.5", "0", "0.04", "800", "0"]);
    let (same, options) = Options::extract(&old).unwrap();
    assert_eq!(same, old);
    assert!(options.is_none());

    let extra = strings(&[
        "--finite-gas",
        "1e-8",
        "200",
        "600",
        "1e-5",
        "--gas-record",
        "mixture.gas",
        "--gas-component",
        "species.pack",
        "nasa9.pack",
        "transport.pack",
        "--gas-water",
        "none",
        "--gas-phase",
        "frozen",
        "--gas-power",
        "-0.001",
        "--solid-temperature",
        "310",
    ]);
    let args = old.iter().chain(&extra).cloned().collect::<Vec<_>>();
    let (same, options) = Options::extract(&args).unwrap();
    assert_eq!(same, old);
    let options = options.unwrap();
    assert_eq!(options.reservoir.volume_m3, 1e-8);
    assert_eq!(options.reservoir.power_w, -0.001);
    assert_eq!(options.reservoir.solid_temperature_k, 310.0);
    // None of these admission tests opens the deliberately nonexistent files.
    for invalid in [
        strings(&[
            "--finite-gas",
            "0",
            "200",
            "600",
            "1e-5",
            "--gas-record",
            "missing",
        ]),
        strings(&[
            "--finite-gas",
            "1e-8",
            "600",
            "200",
            "1e-5",
            "--gas-record",
            "missing",
        ]),
        strings(&["--finite-gas", "1e-8", "200", "600", "1e-5"]),
        strings(&["--gas-record", "missing"]),
        strings(&["--gas-power", "NaN"]),
        strings(&["--solid-temperature", "310"]),
        strings(&["--finite-gas", "1e-8", "200", "600"]),
    ] {
        assert!(Options::extract(&invalid).is_err());
    }
}

#[test]
fn finite_gas_volume_and_caloric_source_change_actual_material_motion_and_close_energy() {
    scoped(|cx, gate| {
        let input = source::Source::table_text(include_str!("material.txt"), gate).unwrap();
        let motion = motion();
        let small_gas = parcel(1e-8, 3.5);
        let mass = small_gas.mass_kg();
        let source_pack = small_gas.state().source_packs()[0].nasa9_pack.to_hex();
        let small = simulate(cx, gate, &input, small_gas.clone(), motion, reservoir(1e-8)).unwrap();
        let large = simulate(cx, gate, &input, parcel(1e-4, 3.5), motion, reservoir(1e-4)).unwrap();
        let changed =
            simulate(cx, gate, &input, parcel(1e-8, 4.5), motion, reservoir(1e-8)).unwrap();
        let replay = simulate(cx, gate, &input, small_gas, motion, reservoir(1e-8)).unwrap();
        assert_eq!(small.records, replay.records);
        assert!(small.description.contains(&source_pack));
        assert_eq!(small.records.len(), motion.count);
        for run in [&small, &large, &changed] {
            for row in &run.records {
                assert!(row.iter().all(|value| value.is_finite()));
                assert!(row[12].abs() <= row[13]);
                assert_eq!(row[9], 0.0);
                assert_eq!(row[10], 0.0);
                assert!((row[5] / 100000.0 - row[4] / 300.0).abs() < 1e-12);
            }
            let last = run.records.last().unwrap();
            assert_eq!(last[0], motion.duration);
            assert!(last[6] > 0.0 && last[7] >= 0.0 && last[8] > 0.0);
            assert!(last[4] > 300.0);
        }
        let end = small.records.last().unwrap();
        assert!(end[4] > large.records.last().unwrap()[4]);
        assert!(
            small
                .records
                .iter()
                .zip(&large.records)
                .any(|(a, b)| (a[2] - b[2]).abs() > 1e-12)
        );
        assert!(
            small
                .records
                .iter()
                .zip(&changed.records)
                .any(|(a, b)| (a[2] - b[2]).abs() > 1e-12)
        );
        // Independent integral of the fixture's NASA heat capacity:
        // cv(T) = R/M * (2.5 + 0.001 T). No solver report supplies this oracle.
        let temperature = end[4];
        let stored_gas = mass * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
            * (2.5 * (temperature - 300.0) + 0.0005 * (temperature * temperature - 300.0 * 300.0));
        let gas_allowance = motion.count as f64 * mass * 1e-7 + 1e-11;
        assert!((stored_gas - end[8]).abs() <= gas_allowance);
    });
}

#[test]
fn gas_heating_reverses_transfer_and_refusals_return_no_partial_history() {
    scoped(|cx, gate| {
        let input = source::Source::table_text(include_str!("material.txt"), gate).unwrap();
        let gas = parcel(1e-8, 3.5);
        let config = Reservoir {
            power_w: 0.001,
            ..reservoir(1e-8)
        };
        let quiet = Motion {
            velocity: 0.0,
            duration: 0.001,
            count: 20,
            ..motion()
        };
        let warmed = simulate(cx, gate, &input, gas.clone(), quiet, config).unwrap();
        let end = warmed.records.last().unwrap();
        assert!(end[8] < 0.0);
        assert!(end[3] > 300.0 && end[4] >= end[3]);
        assert!((end[9] - config.power_w * quiet.duration).abs() <= 64.0 * f64::EPSILON * end[9]);
        assert_eq!(end[6], 0.0);
        assert_eq!(end[7], 0.0);
        assert!(end[12].abs() <= end[13]);
        let invalid = Reservoir {
            solid_temperature_k: 350.0,
            ..config
        };
        assert!(simulate(cx, gate, &input, gas.clone(), quiet, invalid).is_err());
        let outside = Reservoir {
            minimum_temperature_k: 310.0,
            ..config
        };
        assert!(simulate(cx, gate, &input, gas.clone(), quiet, outside).is_err());
        gate.request();
        assert!(simulate(cx, gate, &input, gas, quiet, config).is_err());
    });
}
