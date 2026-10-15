use super::*;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_qty::semantic::{Composition, CompositionBasis};

const INPUT: &str = include_str!("../../../../examples/reed-duct.performance");
const BASE: &str = "mixture 293.15 100000 mole 2 N2 0.7 O2 0.3";
const CO2: &str = "mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4";

fn source(ambient: &str) -> String {
    // Explicit common geometry, admitted at Nyquist even in the slower mixture.
    INPUT.replace("ambient 293.15 101325 0", &format!("ambient {ambient}"))
        .replace("cylinder 0.0022", "cylinder 0.0015")
}
fn render(text: &str, block: usize) -> Vec<f64> {
    let mut source = ReedPerformance::from_bytes(text.as_bytes(), block).unwrap();
    let mut pressure = vec![0.0; source.info().samples as usize];
    for part in pressure.chunks_mut(block) { source.block(part).unwrap(); }
    assert_eq!(source.remaining_samples(), 0);
    pressure
}
fn changed(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    assert!(a.iter().chain(b).all(|p| p.is_finite()));
    let norm: f64 = a.iter().map(|p| p*p).sum();
    let difference: f64 = a.iter().zip(b).map(|(p,q)| (p-q)*(p-q)).sum();
    assert!(norm > 0.0 && difference > 1e-8*norm, "gas change must reach physical pressure");
}

#[test]
fn complete_reed_input_matches_direct_physics_for_legacy_air_and_all_five_species() {
    let species = [GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::Argon,
        GasSpecies::CarbonDioxide, GasSpecies::WaterVapor];
    let mixture = resolve_gas_mixture(&species,
        &Composition::new(CompositionBasis::MoleFraction, vec![0.70,0.20,0.05,0.04,0.01]).unwrap(),
        300.0, 100000.0).unwrap();
    for (declaration, gas) in [
        ("293.15 101325 0", GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap()),
        ("303.15 90000 0.5", GasState::try_new_moist_air(303.15, 90000.0, 0.5).unwrap()),
        ("mixture 300 100000 mole 5 nitrogen 0.70 oxygen 0.20 argon 0.05 carbon-dioxide 0.04 water-vapor 0.01", mixture),
    ] {
        for mass in [0.0, 1e-5] {
            let stiffness = if mass == 0.0 { 0.0 } else { 100.0 };
            let text = source(declaration).replace("reed 0.0004 0.013 6000 0 0 0.35",
                &format!("reed 0.0004 0.013 6000 {mass} {stiffness} 0.35"));
            let loaded = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
            assert_eq!(*loaded.ambient_gas(), gas);
            assert_eq!(loaded.info().massive_reed, mass > 0.0);
            let duct = Duct { segments: vec![Segment::Cylinder { radius: 0.0015, length: 0.5 }] };
            let mut voice = ReedBoreVoice::new(&duct, &gas, BeatingReed {
                rest_opening_m: 0.0004, width_m: 0.013, closing_pressure_pa: 6000.0,
                mass_kg: mass, stiffness_n_m: stiffness, damping_ratio: 0.35,
                blowing_pressure_pa: 0.0, attack_s: 0.0,
            }, Termination::UnflangedOpen, PlateBank::default(), 1.0, 48000, 4801, None).unwrap();
            let schedule = GestureSchedule::from_canonical_bytes(text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
            let actual = render(&text, 37);
            for (n, pressure) in actual.iter().enumerate() {
                voice.set_blowing_pressure(schedule.sample("blow", n as u64*700/48000).unwrap());
                let mut expected = [0.0];
                voice.step_block(&mut expected).unwrap();
                assert_eq!(pressure.to_bits(), expected[0].to_bits(), "{declaration} mass={mass} sample={n}");
            }
        }
    }
}

#[test]
fn composition_and_pressure_reach_sound_without_geometry_or_gesture_retuning() {
    let base = source(BASE);
    let enriched = source(CO2);
    let low = source(&BASE.replace("100000", "50000"));
    let a = ReedPerformance::from_bytes(base.as_bytes(), 37).unwrap();
    let b = ReedPerformance::from_bytes(enriched.as_bytes(), 37).unwrap();
    let c = ReedPerformance::from_bytes(low.as_bytes(), 37).unwrap();
    assert!(b.ambient_gas().sound_speed < a.ambient_gas().sound_speed);
    assert_eq!(a.ambient_gas().sound_speed, c.ambient_gas().sound_speed);
    assert_eq!(a.ambient_gas().density, 2.0*c.ambient_gas().density);
    assert_eq!(a.ambient_gas().characteristic_impedance, 2.0*c.ambient_gas().characteristic_impedance);
    assert_eq!(a.info().full_scale_pa, b.info().full_scale_pa);
    assert_eq!(a.info().compiled_controls, b.info().compiled_controls);
    let pa = render(&base, 37);
    let pb = render(&enriched, 37);
    let pc = render(&low, 37);
    changed(&pa, &pb);
    changed(&pa, &pc);
    assert_eq!(pb, render(&enriched, 113));
}

#[test]
fn mass_basis_keeps_material_semantics_in_the_reed_input() {
    let text = source("mixture 300 100000 mass 3 N2 0.7 O2 0.2 Ar 0.1");
    let expected = resolve_gas_mixture(&[GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::Argon],
        &Composition::new(CompositionBasis::MassFraction, vec![0.7,0.2,0.1]).unwrap(), 300.0, 100000.0).unwrap();
    assert_eq!(*ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap().ambient_gas(), expected);
}

#[test]
fn source_line_and_material_failures_are_retained_before_rendering() {
    for ambient in ["mixture 300 100000 mole 1 He 1",
        "mixture 300 100000 mole 2 N2 0.5 nitrogen 0.5",
        "mixture 300 100000 mole 2 N2 70 O2 30",
        "mixture 300 100000 mole 1 H2O 1", "mixture 400 100000 mole 1 N2 1"] {
        let error = match ReedPerformance::from_bytes(source(ambient).as_bytes(), 37) {
            Err(error) => error, Ok(_) => panic!("accepted {ambient}"),
        };
        assert!(matches!(error, ReedPerformanceError::Gas { line: 3, .. }));
        assert!(std::error::Error::source(&error).is_some());
    }
    for ambient in ["mixture", "mixture 300 100000 mole 0", "mixture 300 100000 mole 99999999999999999999",
        "mixture 300 100000 volume 1 N2 1", "mixture NaN 100000 mole 1 N2 1",
        "mixture 300 100000 mole 1 N2 1 extra"] {
        assert!(ReedPerformance::from_bytes(source(ambient).as_bytes(), 37).is_err(), "accepted {ambient}");
    }
}

#[test]
fn changed_gas_does_not_bypass_the_owners_radiation_band_refusal() {
    // This geometry is valid for the old dry air, but ka exceeds one in CO2
    // at the same 24 kHz Nyquist. Do not shrink the radius or drop the load.
    assert!(ReedPerformance::from_bytes(INPUT.as_bytes(), 37).is_ok());
    let text = INPUT.replace("ambient 293.15 101325 0", "ambient mixture 293.15 100000 mole 1 CO2 1");
    assert!(matches!(ReedPerformance::from_bytes(text.as_bytes(), 37), Err(ReedPerformanceError::Model(_))));
}
