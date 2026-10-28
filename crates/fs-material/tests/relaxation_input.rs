//! Source adapter regressions: the native data reaches the same material law.
use fs_material::gas_species::input::read_gas_record;
use fs_material::relaxation::input::{medium_from_record, read_relaxation};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

// Synthetic rates exercise data connectivity; they are not CO2 measurements.
const DECLARED: &str =
    "declared 2 CO2 1 vibration 1 960 2 power 1.01325 293.15 0 273.15 323.15 N2 0";

fn medium(text: &str) -> fs_material::relaxation::RelaxingMedium {
    let mut fields = text.split_ascii_whitespace();
    let record = read_gas_record(&mut fields, None, None).unwrap();
    assert!(fields.next().is_none());
    medium_from_record(&record).unwrap()
}

#[test]
fn native_reservoirs_equal_the_original_programmatic_law() {
    let gas = medium("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2");
    let mut fields = DECLARED.split_ascii_whitespace();
    let parsed = read_relaxation(&mut fields).unwrap();
    assert!(fields.next().is_none());
    let direct = RelaxationAbsorption::new(vec![
        SpeciesRelaxation {
            species: SpeciesId::new("N2").unwrap(),
            processes: vec![],
        },
        SpeciesRelaxation {
            species: SpeciesId::new("CO2").unwrap(),
            processes: vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 960.0,
                    degeneracy: 2,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 1.01325,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 273.15,
                    t_max_k: 323.15,
                },
            }],
        },
    ])
    .unwrap();
    assert_eq!(parsed, direct);
    let a = parsed.evaluate(&gas).unwrap();
    let b = direct.evaluate(&gas).unwrap();
    assert_eq!(a, b);
    for hz in [20.0, 440.0, 1000.0, 20000.0] {
        assert_eq!(
            a.causal_wavenumber(core::f64::consts::TAU * hz),
            b.causal_wavenumber(core::f64::consts::TAU * hz)
        );
    }
}

#[test]
fn native_pressure_and_composition_reach_reservoir_rates_and_strengths() {
    let law = read_relaxation(&mut DECLARED.split_ascii_whitespace()).unwrap();
    let a = law
        .evaluate(&medium("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2"))
        .unwrap();
    let b = law
        .evaluate(&medium("mixture 293.15 50662.5 mole 2 N2 0.8 CO2 0.2"))
        .unwrap();
    assert_eq!(
        b.processes()[0].relaxation_time_s,
        2.0 * a.processes()[0].relaxation_time_s
    );
    let c = law
        .evaluate(&medium("mixture 293.15 101325 mole 2 N2 0.6 CO2 0.4"))
        .unwrap();
    assert!(c.processes()[0].heat_capacity_j_kg_k > a.processes()[0].heat_capacity_j_kg_k);
    assert_ne!(
        c.relaxing_wavenumber(10000.0).unwrap(),
        a.relaxing_wavenumber(10000.0).unwrap()
    );
}

#[test]
fn collision_partners_and_temperature_domains_remain_mandatory() {
    let text = "declared 2 CO2 1 vibration 1 960 2 partners 2 N2 0 0 290 300 CO2 0 0 290 300 N2 0";
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let state = law
        .evaluate(&medium("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2"))
        .unwrap();
    assert_eq!(state.processes()[0].relaxation_time_s, 1.0 / 101325.0);
    assert!(
        law.evaluate(&medium("mixture 310 101325 mole 2 N2 0.8 CO2 0.2"))
            .is_err()
    );
    let incomplete = text.replace(
        "partners 2 N2 0 0 290 300 CO2 0 0 290 300",
        "partners 1 CO2 0 0 290 300",
    );
    let law = read_relaxation(&mut incomplete.split_ascii_whitespace()).unwrap();
    assert!(
        law.evaluate(&medium("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2"))
            .is_err()
    );
}

#[test]
fn molecular_records_do_not_infer_species_or_iso_eligibility() {
    let iso = read_relaxation(&mut "iso9613-terrestrial-air".split_ascii_whitespace()).unwrap();
    assert!(
        iso.evaluate(&medium("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2"))
            .is_err()
    );
    let legacy = read_gas_record(
        &mut "293.15 101325 0.4".split_ascii_whitespace(),
        None,
        None,
    )
    .unwrap();
    assert!(medium_from_record(&legacy).is_err());
    let helium = read_relaxation(&mut "declared 1 He 0".split_ascii_whitespace()).unwrap();
    assert!(
        helium
            .evaluate(&medium("mixture 293.15 101325 mole 1 He 1"))
            .unwrap()
            .processes()
            .is_empty()
    );
}

#[test]
fn bounded_native_admission_is_atomic_and_leaves_caller_suffix() {
    for text in [
        "declared 129",
        "declared 1 N2 33",
        "declared 1 N2 1 vibration 17",
        "declared 1 N2 1 rotation 1 power 1 300 0 200 400",
        "declared 1 N2 1 vibration 1 NaN 1 power 1 300 0 200 400",
        "declared 1 N2 1 vibration 1 1000 1 partners 17",
        "declared 2 N2 0 N2 0",
        "iso9613",
        "declared 0",
    ] {
        let mut fields = text.split_ascii_whitespace();
        assert!(read_relaxation(&mut fields).is_err(), "{text}");
        assert_eq!(
            fields.collect::<Vec<_>>(),
            text.split_ascii_whitespace().collect::<Vec<_>>()
        );
    }
    let mut fields = "declared 1 He 0 caller-suffix".split_ascii_whitespace();
    read_relaxation(&mut fields).unwrap();
    assert_eq!(fields.next(), Some("caller-suffix"));
}
