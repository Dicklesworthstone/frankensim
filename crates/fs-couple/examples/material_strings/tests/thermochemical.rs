//! Original species/NASA-9/transport data reach the coupled string/body/cavity.
use super::*;
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, ThermochemicalSources};
use fs_material::thermochemical_gas::packs::SpeciesPackRef;

#[path = "../../../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

fn source() -> String {
    String::from_utf8(spec().bytes).unwrap().replace(
        "gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
        "gas thermochemical 293.15 101325 mole 1 bow-atmosphere 1")
}

fn registry<'a>(species: &'a [SpeciesPackRef<'a>]) -> ThermochemicalSources<'a> {
    ThermochemicalSources {
        species, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    }
}

#[test]
fn thermochemical_source_changes_reach_actual_coupled_pressure_and_keep_original_receipts() {
    let original = support::packs("bow-atmosphere", 0.028, 3.5, 200.0);
    let changed = support::packs("bow-atmosphere", 0.028, 4.5, 200.0);
    let original_bindings = [original.binding()]; let replacement_bindings = [changed.binding()];
    let a = input::Spec::parse_with_gas_sources(source().as_bytes(), Some(&registry(&original_bindings))).unwrap();
    let b = input::Spec::parse_with_gas_sources(source().as_bytes(), Some(&registry(&replacement_bindings))).unwrap();
    assert_eq!(a.bytes, b.bytes);
    assert_eq!(a.gas.density, b.gas.density);
    assert_ne!(a.gas.sound_speed, b.gas.sound_speed);
    assert_ne!(a.gas.specific_heat_cp, b.gas.specific_heat_cp);
    let GasRecord::Thermochemical(record) = &a.gas_record else { panic!("source receipt lost"); };
    assert_eq!(record.source_packs()[0].species_pack, original.species.content_hash());
    assert_eq!(record.source_packs()[0].nasa9_pack, original.nasa9.content_hash());
    assert_eq!(record.source_packs()[0].transport_pack, original.transport.content_hash());
    scoped(|cx, gate| {
        let mut first = compile(cx, gate, a.clone(), 7800.0);
        let mut second = compile(cx, gate, b, 7800.0);
        assert_eq!(first.catalog["wire.fsmatpk"].pack.content_hash(),
            second.catalog["wire.fsmatpk"].pack.content_hash());
        let bank = first.performance.bank().source().bank();
        assert_eq!(bank.gas(), a.gas_record.state());
        assert_eq!(bank.material_cavity().unwrap().cavity.gas(), a.gas_record.state());
        let p = pressure(&mut first, 1); let q = pressure(&mut second, 37);
        assert!(p.iter().zip(q).any(|(p, q)| (p - q).abs() > 1e-12));
        let mut replay = compile(cx, gate, a.clone(), 7800.0);
        assert_eq!(p, pressure(&mut replay, 128));
        let mut lower_pressure = input::Spec::parse_with_gas_sources(
            source().replace("293.15 101325", "293.15 60000").as_bytes(),
            Some(&registry(&original_bindings))).unwrap();
        assert!(lower_pressure.gas.density < a.gas.density);
        // A numeric overwrite cannot retain the authority of a different source.
        lower_pressure.gas = a.gas;
        let mut calls = 0;
        assert!(lower_pressure.compile(cx, gate, |name| {
            calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
        }).is_err());
        assert_eq!(calls, 0);
    });
}

#[test]
fn thermochemical_admission_requires_complete_unambiguous_in_domain_sources() {
    assert!(input::Spec::parse(source().as_bytes()).is_err());
    let original = support::packs("bow-atmosphere", 0.028, 3.5, 200.0);
    let late = support::packs("bow-atmosphere", 0.028, 3.5, 350.0);
    let mut missing = original.binding(); missing.nasa9 = None;
    for bindings in [vec![missing], vec![late.binding()],
        vec![original.binding(), original.binding()], vec![]] {
        let registry = ThermochemicalSources {
            species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
        };
        assert!(input::Spec::parse_with_gas_sources(source().as_bytes(), Some(&registry)).is_err());
    }
    let bindings = [original.binding()];
    let registry = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let legacy = spec();
    let parsed = input::Spec::parse_with_gas_sources(&legacy.bytes, Some(&registry)).unwrap();
    assert_eq!(legacy.gas_record, parsed.gas_record);
    assert_eq!(legacy.gas, parsed.gas);
    assert_eq!(legacy.bytes, parsed.bytes);
}
