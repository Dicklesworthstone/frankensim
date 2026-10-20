use super::*;
use fs_material::gas_species::input::ThermochemicalPhase;

// Reuse the existing normalized synthetic species/NASA9/transport fixtures.
#[path = "../../../../../../render/schedule/reed/field_input/tests/support.rs"]
mod support;

const GAS: &str = "thermochemical 293.15 101325 mole 1 specimen-gas 1";

fn original_artifact(sources: &ThermochemicalSources<'_>) -> String {
    let record = read_gas_record(&mut GAS.split_ascii_whitespace(), None, Some(sources)).unwrap();
    let GasRecord::Thermochemical(gas) = record else { panic!("original gas receipt"); };
    let mut exterior = format!("exterior {GAS}\ngas_sources {}\n", gas.source_packs().len());
    for source in gas.source_packs() {
        exterior.push_str(&format!("gas_source {} {} {} {}\n",
            source.species, source.species_pack, source.nasa9_pack, source.transport_pack));
    }
    exterior.push_str("gas_state");
    for value in gas_values(gas.state()) { exterior.push_str(&format!(" {value}")); }
    SOURCE.replace("exterior 293.15 101325 0", &exterior)
}

#[test]
fn original_thermochemical_packs_and_complete_numerical_state_are_bound_to_the_artifact() {
    let packs = support::packs("specimen-gas", 0.028, 3.5, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let text = original_artifact(&sources);
    let original = MimoRadiationArtifact::from_bytes_with_sources(text.as_bytes(),
        Some(&sources), &CancelGate::new_clock_free()).unwrap();
    let GasRecord::Thermochemical(gas) = original.gas_record() else { panic!("retained original receipt"); };
    assert_eq!(gas.source_packs()[0].nasa9_pack, packs.nasa9.content_hash());
    assert_eq!(gas.source_packs()[0].species_pack, packs.species.content_hash());
    assert_eq!(gas.source_packs()[0].transport_pack, packs.transport.content_hash());
    assert_eq!(*gas.state(), original.source().exterior);

    // Same canonical species id, changed NASA9 heat-capacity coefficients.
    // Identical radiation bytes must not acquire this newly resolved gas.
    let replacement = support::packs("specimen-gas", 0.028, 4.5, 200.0);
    let replacements = [replacement.binding()];
    let changed = ThermochemicalSources {
        species: &replacements, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let changed_gas = read_gas_record(&mut GAS.split_ascii_whitespace(), None, Some(&changed)).unwrap();
    assert_ne!(changed_gas.state().sound_speed, gas.state().sound_speed);
    assert_ne!(changed_gas.state().specific_heat_cp, gas.state().specific_heat_cp);
    let refused = MimoRadiationArtifact::from_bytes_with_sources(text.as_bytes(),
        Some(&changed), &CancelGate::new_clock_free()).unwrap_err();
    assert!(refused.to_string().contains("source pack/order"));

    // Even updating a receipt cannot silently relabel the original snapshot.
    let changed_receipt = text.replace(&packs.nasa9.content_hash().to_hex(),
        &replacement.nasa9.content_hash().to_hex());
    let refused = MimoRadiationArtifact::from_bytes_with_sources(changed_receipt.as_bytes(),
        Some(&changed), &CancelGate::new_clock_free()).unwrap_err();
    assert!(refused.to_string().contains("gas snapshot"));

    let snapshot = text.lines().find(|line| line.starts_with("gas_state ")).unwrap();
    for (field, _) in gas_values(gas.state()).iter().enumerate() {
        let mut words = snapshot.split_ascii_whitespace().map(String::from).collect::<Vec<_>>();
        let value: f64 = words[field + 1].parse().unwrap();
        words[field + 1] = f64::from_bits(value.to_bits() + 1).to_string();
        let altered = text.replace(snapshot, &words.join(" "));
        assert!(MimoRadiationArtifact::from_bytes_with_sources(altered.as_bytes(),
            Some(&sources), &CancelGate::new_clock_free()).unwrap_err().to_string().contains("gas snapshot"));
    }
    for altered in [
        text.replace("gas_sources 1", "gas_sources 0"),
        text.replace("gas_sources 1", "gas_sources 129"),
        text.replace("gas_source specimen-gas", "gas_source another-species"),
        SOURCE.replace("exterior 293.15 101325 0", &format!("exterior {GAS}")),
    ] {
        assert!(MimoRadiationArtifact::from_bytes_with_sources(altered.as_bytes(),
            Some(&sources), &CancelGate::new_clock_free()).is_err());
    }
}
