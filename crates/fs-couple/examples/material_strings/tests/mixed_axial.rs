//! Per-string longitudinal selection: `measured-strings-transverse` leaves measured whole strings transverse
//! (they carry no axial law) while the bonded wire keeps its axial coordinates; without the word the record
//! refuses a measured string exactly as before. A selection that names no string refuses with its own text.
use super::*;

const RECORD: &str = "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2";
const WORD: &str = " measured-strings-transverse";
const MEASURED_STRING: &str = "string cello-g-whole.fsmatpk 0.7 0.000575 tension 116.3 4 0 0.8 3000 0.1\nsection measured 1 temperature_known 0 playing_length diameter_dominant tension_dominant linear_density_dominant bending_stiffness_dominant no-intrinsic-loss\npolarizations two 0.6 0.8 1 0";

/// The bonded wire alone, with the axial record (the existing longitudinal fixture).
fn wire_text() -> String {
    include_str!("../longitudinal-cylinder-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 129 0.01 128")
}
/// The bonded wire (string 0) followed by a measured whole string (string 1), the axial record with `word`.
fn mixed_text(word: &str, gas: bool) -> String {
    let t = wire_text()
        .replace("strings 1\n", "strings 2\n")
        .replace("polarizations two 0.6 0.8 1 0\n", &format!("polarizations two 0.6 0.8 1 0\n{MEASURED_STRING}\n"))
        .replace(RECORD, &format!("{RECORD}{word}"))
        .replace("prepared-linear 64 67108864 200000000000", "prepared-linear 256 268435456 800000000000");
    if gas { t } else { t.replace("cylinder-gas-loading 0.001 30 4\n", "") }
}
fn build(cx: &Cx<'_>, gate: &CancelGate, text: &str) -> Result<Ready> {
    input::Spec::parse(text.as_bytes())?.compile(cx, gate, |name| {
        Ok(if name == "cello-g-whole.fsmatpk" {
            super::measured::measured_pack(1.0)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        })
    })
}

#[test]
fn without_the_word_a_measured_string_still_refuses_the_axial_record_with_todays_text() {
    scoped(|cx, gate| {
        let err = build(cx, gate, &mixed_text("", false)).err().expect("the record without the word must refuse a measured string");
        assert!(
            err.to_string().contains("measured whole strings supply transverse mass and bending stiffness, not torsion, thermal, axial"),
            "{err}"
        );
    });
}

#[test]
fn the_word_leaves_the_measured_string_transverse_and_the_wire_keeps_its_axial_modes() {
    scoped(|cx, gate| {
        let ready = build(cx, gate, &mixed_text(WORD, false)).unwrap();
        let source = ready.performance.bank().source().bank();
        let axial = source.longitudinal_waves().unwrap();
        assert_eq!(axial.strings().len(), 1, "only the wire carries the axial law");
        let law = &axial.strings()[0];
        assert_eq!(law.string_index(), 0);
        assert!(axial.string(0).is_some() && axial.string(1).is_none());
        assert_eq!(law.model().axial_modes(), &[1, 2]);
        let components = ready.performance.bank().system().components();
        assert_eq!(components[law.axial_component()].modes().len(), 2);
        assert!(law.axial_component() >= 2, "the axial component follows every original string");
        // The measured string keeps its own transverse component and no axial one.
        assert_eq!(components[1].modes().len(), 8, "four modes per polarization, untouched");
    });
}

#[test]
fn the_wire_axial_law_is_the_single_string_build_s_own() {
    scoped(|cx, gate| {
        // Both builds carry no cylinder gas loading: the gas row adds the air's mass to the transverse density
        // (1.2 * pi * r^2 = 3.4e-7 kg/m on the 0.3 mm wire), so a build with it is a different transverse card.
        let wire_no_gas = wire_text().replace("cylinder-gas-loading 0.001 30 4\n", "");
        assert!(!wire_no_gas.contains("cylinder-gas-loading"), "the gas row must be gone from the single build");
        let single = compile(cx, gate, input::Spec::parse(wire_no_gas.as_bytes()).unwrap(), 7800.0);
        let mixed = build(cx, gate, &mixed_text(WORD, false)).unwrap();
        let a = single.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0].clone();
        let b = mixed.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0].clone();
        assert_eq!(a.model().axial_modes(), b.model().axial_modes());
        assert_eq!(a.model().axial_omegas(), b.model().axial_omegas());
        assert_eq!(a.solid_linear_density_kg_m(), b.solid_linear_density_kg_m());
        assert_eq!(a.axial_rigidity_n(), b.axial_rigidity_n());
        assert_eq!(a.transverse_linear_density_kg_m(), b.transverse_linear_density_kg_m());
    });
}

#[test]
fn the_mixed_bank_renders_bit_stably_across_block_sizes_and_closes_its_energy_balance() {
    scoped(|cx, gate| {
        let mut a = build(cx, gate, &mixed_text(WORD, false)).unwrap();
        let first = pressure(&mut a, 37);
        assert!(first.iter().any(|p| p.abs() > 1e-15), "the mixed bank must sound");
        for block in [1, 128] {
            let mut b = build(cx, gate, &mixed_text(WORD, false)).unwrap();
            assert_eq!(pressure(&mut b, block), first);
        }
        let frame = a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn axial_gas_memory_follows_only_the_selected_string_and_its_component_addresses_line_up() {
    scoped(|cx, gate| {
        let ready = build(cx, gate, &mixed_text(WORD, true)).unwrap();
        let source = ready.performance.bank().source().bank();
        assert!(source.cylinder_gas_loading().is_some());
        let law = &source.longitudinal_waves().unwrap().strings()[0];
        let fluid = law.axial_fluid_component().expect("the selected wire has axial gas memory");
        assert_eq!(fluid, law.axial_component() + 1, "one selected string: its fluid component follows its axial one");
        let components = ready.performance.bank().system().components();
        let branches = source.cylinder_gas_loading().unwrap().strings()[0].flow().branches().len();
        assert_eq!(components[fluid].modes().len(), 2 * branches);
        assert_eq!(fluid, components.len() - 1, "no axial gas component for the measured string");
    });
}

#[test]
fn the_word_with_only_measured_strings_names_no_string_and_refuses_with_its_own_text() {
    let text = include_str!("../measured-cello-g.material-strings")
        .replace("prepared-linear 64 67108864 200000000000", &format!("prepared-linear 64 67108864 200000000000\n{RECORD}{WORD}"));
    let err = input::Spec::parse(text.as_bytes()).err().expect("a selection naming no string must refuse");
    // Refused while parsing, before any bank or zero-axial network exists.
    let text = err.to_string();
    assert!(text.contains("names no string: every string is a measured whole string"), "{text}");
}
