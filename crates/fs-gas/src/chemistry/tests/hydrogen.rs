use super::*;
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, MOLAR_MASSES_KG_PER_MOL, cantera_h2o2};

// Independent Cantera 3.2.0, ideal-gas ohmech, adiabatic constant volume.
// Pinned source 4a8358eb80cfeb50474386b5f9ec0b3a83519889; exact YAML blob
// 65e6321a60bafb1feb4f3b1be5347e6d1e89c937. The source's upper h/s integration
// constants receive the SAME explicitly selected continuity adjustment as the
// original FrankenSim caloric owner; Cp and every reaction law are unchanged.
// Reproduce with data/reference/thermochem/hydrogen_cantera32.py in an isolated
// development environment containing Cantera 3.2.0. No production dependency.
// Cantera controls: rtol=1e-11, atol=1e-24, max_dt=duration/50, max_steps=200000.
// H2:O2:N2 = 1:1:3.76 (phi=0.5); 101325 Pa, 1 litre; NO radical seed.
struct Reference {
    initial_temperature: f64,
    duration: f64,
    temperature: f64,
    pressure: f64,
    energy: f64,
    amounts: [f64; 10],
}
const REFERENCES: [Reference; 2] = [
    Reference {
        initial_temperature: 1000.0,
        duration: 1e-3,
        temperature: 2458.4913945555377,
        pressure: 228725.14079886756,
        energy: 161.2793063841124,
        amounts: [
            0.0,
            6.840225672803266e-6,
            2.1604586538868265e-5,
            0.002027389303305023,
            5.9750418141946714e-9,
            1.2735377507825906e-7,
            0.007955139300154143,
            2.651343322649553e-5,
            0.0010253989575270015,
            0.0001264897653307908,
        ],
    },
    Reference {
        initial_temperature: 1200.0,
        duration: 1e-4,
        temperature: 2490.1870891852154,
        pressure: 194714.2134198162,
        energy: 184.82023276918432,
        amounts: [
            0.0,
            3.0514265490861635e-5,
            5.43936166029164e-5,
            0.001601192986638355,
            4.0181079890848725e-9,
            1.1665338714300007e-7,
            0.006629282750128453,
            6.862675824435062e-5,
            0.000835875536851778,
            0.00018440206721522228,
        ],
    },
];

fn hydrogen_controls(relative: f64) -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 300.0,
        maximum_temperature_k: 3500.0,
        energy_tolerance_j: 1e-8,
        relative_tolerance: relative,
        absolute_amount_tolerance_mol: 1e-16,
        initial_step_s: 1e-9,
        maximum_step_s: 2e-5,
        maximum_steps: 100000,
        ..IsochoricKineticsSolve::default()
    }
}

fn hydrogen_run(case: &Reference, relative: f64) -> IsochoricKineticsResult {
    let mechanism = cantera_h2o2(HydrogenCalorics::ContinuousAt300K).unwrap();
    let unit = 101325.0 * 0.001 / (R * case.initial_temperature * 5.76);
    advance_adiabatic(
        &mechanism,
        &[
            (id("H2"), Amount::new(unit)),
            (id("O2"), Amount::new(unit)),
            (id("N2"), Amount::new(3.76 * unit)),
        ],
        Volume::new(0.001),
        Temperature::new(case.initial_temperature),
        Time::new(case.duration),
        hydrogen_controls(relative),
        || false,
    )
    .unwrap()
}

fn endpoint_error(result: &IsochoricKineticsResult, case: &Reference) -> f64 {
    let total: f64 = case.amounts.iter().sum();
    let species = result
        .amounts()
        .iter()
        .zip(case.amounts)
        .map(|(a, b)| (a.value() - b).abs())
        .sum::<f64>()
        / total;
    species + (result.temperature().value() - case.temperature).abs() / case.temperature
}

#[test]
fn sourced_hydrogen_induction_and_hot_oxidation_match_cantera_and_refine() {
    for case in &REFERENCES {
        let fine = hydrogen_run(case, 2e-6);
        let t_error = (fine.temperature().value() - case.temperature).abs();
        assert!(t_error < 3.0, "T error={t_error}, {:?}", fine.balance());
        assert!((fine.pressure().value() / case.pressure - 1.0).abs() < 0.002);
        for (i, (a, b)) in fine.amounts().iter().zip(case.amounts).enumerate() {
            assert!(
                (a.value() - b).abs() < 1e-12 + 0.005 * b,
                "species {i}, actual={}, Cantera={b}, {:?}",
                a.value(),
                fine.balance()
            );
        }
        assert_eq!(fine.elapsed_time(), Time::new(case.duration));
        assert!(fine.amounts().iter().all(|n| n.value() >= 0.0));
        assert_eq!(fine.amounts()[0].value(), 0.0); // No argon creation.
        let unit = 101325.0 * 0.001 / (R * case.initial_temperature * 5.76);
        let n = fine.amounts();
        // Independent H/O atom and N2 checks, not only the solver's receipt.
        let hydrogen = n[1].value()
            + 2.0 * (n[2].value() + n[3].value() + n[4].value())
            + n[5].value()
            + n[9].value();
        let oxygen = n[3].value()
            + 2.0 * (n[4].value() + n[5].value() + n[8].value())
            + n[7].value()
            + n[9].value();
        assert!((hydrogen - 2.0 * unit).abs() < 1e-14);
        assert!((oxygen - 2.0 * unit).abs() < 1e-14);
        assert!((n[6].value() - 3.76 * unit).abs() < 1e-16);
        let mass: f64 = n
            .iter()
            .zip(MOLAR_MASSES_KG_PER_MOL)
            .map(|(n, m)| n.value() * m)
            .sum();
        let initial_mass = unit
            * (MOLAR_MASSES_KG_PER_MOL[2]
                + MOLAR_MASSES_KG_PER_MOL[8]
                + 3.76 * MOLAR_MASSES_KG_PER_MOL[6]);
        assert!((mass - initial_mass).abs() < 1e-16);
        assert!((fine.internal_energy().value() - case.energy).abs() < 2e-8);
        assert!(
            fine.balance().energy_residual_j.abs() <= hydrogen_controls(2e-6).energy_tolerance_j
        );
        assert!(fine.balance().entropy_change_j_k > 0.1);
        assert!(fine.balance().maximum_estimated_local_error <= 1.0);
        if case.initial_temperature == 1200.0 {
            let coarse = hydrogen_run(case, 2e-4);
            let fine_error = endpoint_error(&fine, case);
            let coarse_error = endpoint_error(&coarse, case);
            assert!(
                fine_error < 0.3 * coarse_error,
                "coarse={coarse_error:e}, fine={fine_error:e}"
            );
        }
    }
}
