//! Authored interior reservoirs bind the original shell/string/heat owner.
//! Collision coefficients below are synthetic integration data, not gas fits.
use super::*;
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, ThermochemicalSources};

#[path = "../../../../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

const AIR: &str = "mixture 293.15 101325 mole 2 N2 0.8 O2 0.2";
const MOLECULAR: &str =
    "cavity_molecular 2 96 declared 2 N2 1 rotation 2 power 2533.125 293.15 0 280 320 O2 0";
const THERMAL_BODY: &str = include_str!("../../closed-shell-thermal.body");

fn selected() -> String {
    // The outer helium has no entry in the interior law. Borrowing its record
    // would therefore refuse rather than accidentally pass this source test.
    PERFORMANCE
        .replace(
            &format!("gas {AIR}"),
            "gas mixture 293.15 101325 mole 1 He 1",
        )
        .replacen(
            "\nbody ",
            &format!("\ncavity_gas {AIR}\n{MOLECULAR}\nbody "),
            1,
        )
}

#[test]
fn authored_cavity_molecular_uses_the_original_interior_species_and_rates() {
    let spec = source(&selected());
    let declared = spec.cavity_molecular.as_ref().unwrap();
    assert_eq!(declared.relaxation.gas(), spec.cavity_gas());
    assert_ne!(declared.relaxation.gas(), &spec.gas);
    let processes = declared.relaxation.processes();
    assert_eq!(processes.len(), 1);
    assert_eq!(processes[0].species.as_str(), "N2");
    let expected_capacity = 0.8 * spec.cavity_gas().specific_gas_constant;
    assert!(
        (processes[0].heat_capacity_j_kg_k - expected_capacity).abs() < 1e-12 * expected_capacity
    );
    assert!((processes[0].relaxation_time_s - 0.025).abs() < 1e-15);
    let faster = source(&selected().replace("2533.125", "253.3125"));
    assert!(
        (faster
            .cavity_molecular
            .as_ref()
            .unwrap()
            .relaxation
            .processes()[0]
            .relaxation_time_s
            - 0.0025)
            .abs()
            < 1e-15
    );
    let inherited = source(&PERFORMANCE.replacen("\nbody ", &format!("\n{MOLECULAR}\nbody "), 1));
    assert!(inherited.cavity_gas_record.is_none());
    assert_eq!(
        inherited.cavity_molecular.as_ref().unwrap().relaxation,
        declared.relaxation
    );
}

#[test]
fn authored_cavity_molecular_keeps_original_thermochemical_source_authority() {
    let original = support::packs("interior-reservoir", 0.028, 3.5, 200.0);
    let bindings = [original.binding()];
    let registry = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let text = selected()
        .replace(
            &format!("cavity_gas {AIR}"),
            "cavity_gas thermochemical 293.15 101325 mole 1 interior-reservoir 1",
        )
        .replace("declared 2 N2 1", "declared 1 interior-reservoir 1")
        .replace(" O2 0", "");
    assert!(input::Spec::parse(text.as_bytes()).is_err());
    let spec = input::Spec::parse_with_gas_sources(text.as_bytes(), Some(&registry)).unwrap();
    let GasRecord::Thermochemical(record) = spec.cavity_gas_record() else {
        panic!("the interior thermochemical source record was lost");
    };
    assert_eq!(
        record.source_packs()[0].species_pack,
        original.species.content_hash()
    );
    assert_eq!(
        record.source_packs()[0].nasa9_pack,
        original.nasa9.content_hash()
    );
    assert_eq!(
        record.source_packs()[0].transport_pack,
        original.transport.content_hash()
    );
    let molecular = &spec.cavity_molecular.as_ref().unwrap().relaxation;
    assert_eq!(molecular.gas(), record.state());
    assert_eq!(
        molecular.processes()[0].species.as_str(),
        "interior-reservoir"
    );
}

#[test]
fn authored_cavity_molecular_binding_changes_the_native_constitutive_pressure() {
    // The existing debug frontend nests large by-value bank builders; its
    // compile frame alone is about 1.5 MiB. Keep the required construction
    // stack local to this full binding test rather than setting a global limit.
    std::thread::Builder::new()
        .name("authored-cavity-molecular-binding".into())
        .stack_size(8 * 1024 * 1024)
        .spawn(|| {
            scoped(|cx, gate| {
                let spec = source(&selected());
                let declared = spec.cavity_molecular.as_ref().unwrap().clone();
                let profile = body::Profile::parse(THERMAL_BODY.as_bytes()).unwrap();
                assert_eq!(profile.material_memory_states(), 43);
                assert_eq!(profile.material_memory_states_for(&spec), 96);
                let first = ready(cx, gate, spec, THERMAL_BODY, 700.0);
                let second = ready(
                    cx,
                    gate,
                    source(&selected().replace("2533.125", "253.3125")),
                    THERMAL_BODY,
                    700.0,
                );
                let first_bank = first.performance.bank().source().bank();
                let second_bank = second.performance.bank().source().bank();
                let interior = first_bank.spatial_cavity().unwrap();
                let a = interior.thermal_loss().unwrap();
                let b = second_bank
                    .spatial_cavity()
                    .unwrap()
                    .thermal_loss()
                    .unwrap();
                assert_eq!(first_bank.gas(), &first.spec.gas);
                assert_eq!(interior.gas(), first.spec.cavity_gas());
                assert_ne!(interior.gas(), first_bank.gas());
                assert_eq!(a.molecular_relaxation(), Some(&declared.relaxation));
                assert!(a.memory_states() > 0 && a.memory_states() <= 96);
                assert_eq!(
                    first.performance.bank().system().material_memory().len(),
                    a.memory_states()
                );
                assert!(a.molecular_complement_memory_states() > 0);
                assert_eq!(
                    a.instantaneous_bulk_modulus_pa(),
                    b.instantaneous_bulk_modulus_pa()
                );
                assert_ne!(a.decay_rates_per_s(), b.decay_rates_per_s());
                assert!(
                    (a.normalized_bulk_response(100.0).unwrap()[0]
                        - b.normalized_bulk_response(100.0).unwrap()[0])
                        .abs()
                        > 1.0
                );
                // Native integration tests already exercise accepted force/sound/heat
                // trajectories. This frontend seam proves authored source selection and
                // complete binding without repeating long dense time integrations.
            })
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn authored_cavity_molecular_refuses_missing_sources_and_unsupported_owners() {
    let text = selected();
    for invalid in [
        text.replace("cavity_molecular 2 96", "cavity_molecular 0 96"),
        text.replace("cavity_molecular 2 96", "cavity_molecular 1 96"),
        text.replace("cavity_molecular 2 96", "cavity_molecular 129 96"),
        text.replace("cavity_molecular 2 96", "cavity_molecular 2 0"),
        text.replace("cavity_molecular 2 96", "cavity_molecular 2 4097"),
        text.replace("declared 2", "declared 1")
            .replace(" O2 0", ""),
        text.replace("0 280 320 O2", "0 300 320 O2"),
        text.replace(&format!("cavity_gas {AIR}"), "cavity_gas 293.15 101325 0"),
        text.replace(MOLECULAR, &format!("{MOLECULAR} ignored")),
        text.replace(MOLECULAR, &format!("{MOLECULAR}\n{MOLECULAR}")),
        PERFORMANCE
            .replacen("\nbody ", &format!("\n{MOLECULAR}\nbody "), 1)
            .replace("cavity none", "cavity neck 0.003 0.01 0.02 100 1"),
    ] {
        assert!(
            input::Spec::parse(invalid.as_bytes()).is_err(),
            "accepted {invalid}"
        );
    }
    scoped(|cx, gate| {
        let spec = source(&text);
        let thermal = body::Profile::parse(THERMAL_BODY.as_bytes()).unwrap();
        let nonthermal = body::Profile::parse(BODY.as_bytes()).unwrap();
        let mut prepared = spec.clone();
        prepared.prepared_linear = Some(PreparedModalBudget {
            max_modes: 128,
            max_dense_bytes: 67_108_864,
            max_multiply_terms: 200_000_000_000,
        });
        let mut stale = spec.clone();
        stale.cavity_gas_record = None;
        let mut excessive = spec.clone();
        excessive
            .cavity_molecular
            .as_mut()
            .unwrap()
            .options
            .maximum_memory_states = 257;
        for (spec, body) in [
            (spec.clone(), None),
            (spec, Some(nonthermal)),
            (prepared, Some(thermal.clone())),
            (stale, Some(thermal.clone())),
            (excessive, Some(thermal)),
        ] {
            let mut reads = 0;
            assert!(
                spec.compile_with_body::<NormalizedPack>(cx, gate, body, |_| {
                    reads += 1;
                    Err("an unconsumed caloric source must refuse before material I/O".into())
                })
                .is_err()
            );
            assert_eq!(reads, 0);
        }
    });
}
