use super::*;
use crate::kinetics::KineticWorkspaceV1;

#[path = "oracle.rs"]
mod oracle;

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn complete_source_identity_calorics_and_directionality() {
    let mechanism = cantera_gri30(Gri30Calorics::Original).unwrap();
    assert_eq!(mechanism.reactions().len(), 325);
    assert_eq!(mechanism.thermodynamics().models().len(), 53);
    let irreversible: Vec<_> = mechanism
        .reactions()
        .iter()
        .enumerate()
        .filter_map(|(i, r)| (!r.is_reversible()).then_some(i + 1))
        .collect();
    assert_eq!(
        irreversible,
        [
            135, 284, 288, 290, 292, 293, 297, 298, 300, 301, 302, 303, 305, 306, 307, 324
        ]
    );
    let mut counts = [0; 4];
    for (i, reaction) in mechanism.reactions().iter().enumerate() {
        assert_eq!(
            reaction.reaction(),
            &ReactionId::new(format!("cantera-gri30-r{:03}", i + 1)).unwrap()
        );
        counts[match reaction.pressure_dependence() {
            PressureDependenceV1::Independent => 0,
            PressureDependenceV1::ThirdBody(_) => 1,
            PressureDependenceV1::Falloff {
                broadening: FalloffBroadeningV1::Lindemann,
                ..
            } => 2,
            PressureDependenceV1::Falloff {
                broadening: FalloffBroadeningV1::Troe(_),
                ..
            } => 3,
        }] += 1;
    }
    assert_eq!(counts, [284, 12, 3, 26]);
    for (i, model) in mechanism.thermodynamics().models().iter().enumerate() {
        assert_eq!(model.species().as_str(), CANONICAL_SPECIES[i]);
        assert_eq!(model.molar_mass().value(), MOLAR_MASSES_KG_PER_MOL[i]);
        assert!(model.caloric_reference().is_none());
        for (case, t) in [700.0, 1800.0].into_iter().enumerate() {
            let actual = model.evaluate_properties(Temperature::new(t)).unwrap();
            for (value, reference) in [actual.cp().value(), actual.h().value(), actual.s().value()]
                .into_iter()
                .zip(oracle::THERMO[case][i])
            {
                close(value, reference, 8e-13, 1e-9);
            }
        }
    }
    // Channel34's explicit O2 collider remains on both sides and is not
    // multiplied by a second generic-M pressure factor.
    let collider = &mechanism.reactions()[33];
    assert_eq!(
        collider
            .reactants()
            .iter()
            .find(|(id, _)| id.as_str() == "O2")
            .unwrap()
            .1,
        2
    );
    assert_eq!(
        collider
            .products()
            .iter()
            .find(|(id, _)| id.as_str() == "O2")
            .unwrap()
            .1,
        1
    );
    assert_eq!(
        collider.pressure_dependence(),
        &PressureDependenceV1::Independent
    );
}

#[test]
fn all_325_rates_agree_with_independent_cantera_and_conserve_every_element() {
    let mechanism = cantera_gri30(Gri30Calorics::Original).unwrap();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let concentrations: Vec<_> = (0..53).map(|i| 0.01 * (i + 1) as f64).collect();
    for (case, t) in [700.0, 1800.0].into_iter().enumerate() {
        mechanism
            .evaluate(Temperature::new(t), &concentrations, &mut work, || false)
            .unwrap();
        for (j, (&actual, &expected)) in work
            .reaction_rates()
            .iter()
            .zip(&oracle::RATES[case])
            .enumerate()
        {
            assert!(
                (actual - expected).abs() <= 3e-11 * expected.abs() + 1e-25,
                "channel{},T={t}: actual={actual:.17e}, expected={expected:.17e}",
                j + 1
            );
        }
        let elements = mechanism.thermodynamics().elemental_matrix();
        for row in 0..elements.row_count() {
            let terms: Vec<_> = work
                .species_rates()
                .iter()
                .enumerate()
                .map(|(species, rate)| elements.get(row, species).unwrap() as f64 * rate)
                .collect();
            assert!(
                terms.iter().sum::<f64>().abs()
                    <= 3e-14 * terms.iter().map(|x| x.abs()).sum::<f64>()
            );
        }
        let terms: Vec<_> = work
            .species_rates()
            .iter()
            .zip(MOLAR_MASSES_KG_PER_MOL)
            .map(|(rate, mass)| rate * mass)
            .collect();
        assert!(
            terms.iter().sum::<f64>().abs() <= 3e-14 * terms.iter().map(|x| x.abs()).sum::<f64>()
        );
    }
}

#[test]
fn continuous_calorics_are_explicit_bounded_and_keep_every_original_cp() {
    let original = cantera_gri30(Gri30Calorics::Original).unwrap();
    let continuous = cantera_gri30(Gri30Calorics::ContinuousAt300K).unwrap();
    assert!(
        original
            .admit_temperature_interval(Temperature::new(300.0), Temperature::new(3000.0), || false)
            .is_err()
    );
    continuous
        .admit_temperature_interval(Temperature::new(300.0), Temperature::new(3000.0), || false)
        .unwrap();
    assert!(
        continuous
            .admit_temperature_interval(Temperature::new(299.0), Temperature::new(3000.0), || false)
            .is_err()
    );
    assert!(
        continuous
            .admit_temperature_interval(Temperature::new(300.0), Temperature::new(3001.0), || false)
            .is_err()
    );
    for (i, (raw, adjusted)) in original
        .thermodynamics()
        .models()
        .iter()
        .zip(continuous.thermodynamics().models())
        .enumerate()
    {
        let offsets = adjusted.caloric_reference().unwrap();
        assert_eq!(offsets.anchor_temperature().value(), 300.0);
        assert_eq!(offsets.regions()[0].enthalpy_offset().value(), 0.0);
        assert_eq!(offsets.regions()[0].entropy_offset().value(), 0.0);
        assert!(offsets.regions()[1].enthalpy_offset().value().abs() <= 0.2);
        assert!(offsets.regions()[1].entropy_offset().value().abs() <= 2e-4);
        for t in [300.0, 700.0, 999.0, 1001.0, 1400.0, 1800.0, 3000.0] {
            let a = raw.evaluate_properties(Temperature::new(t)).unwrap();
            let b = adjusted.evaluate_properties(Temperature::new(t)).unwrap();
            assert_eq!(a.cp().value().to_bits(), b.cp().value().to_bits());
            if t < data::SPECIES[i].temperatures[1] {
                assert_eq!(a, b);
            }
        }
    }
    // The source has three seams other than1000K; none is rounded away.
    for (name, seam) in [("HCNO", 1382.0), ("HOCN", 1368.0), ("HNCO", 1478.0)] {
        let i = CANONICAL_SPECIES.iter().position(|&s| s == name).unwrap();
        assert_eq!(data::SPECIES[i].temperatures[1], seam);
        let model = &continuous.thermodynamics().models()[i];
        let left = model
            .evaluate_properties(Temperature::new(seam - 1e-7))
            .unwrap();
        let right = model
            .evaluate_properties(Temperature::new(seam + 1e-7))
            .unwrap();
        close(left.h().value(), right.h().value(), 0.0, 5e-5);
        close(left.s().value(), right.s().value(), 0.0, 5e-8);
    }
}

#[test]
fn source_jacobians_cover_irreversible_duplicates_colliders_and_falloff() {
    let mechanism = cantera_gri30(Gri30Calorics::Original).unwrap();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let concentrations: Vec<_> = (0..53).map(|i| 0.01 * (i + 1) as f64).collect();
    let t = 1800.0;
    mechanism
        .evaluate(Temperature::new(t), &concentrations, &mut work, || false)
        .unwrap();
    let dc = work.rate_jacobian().to_vec();
    let dt = work.temperature_derivatives().to_vec();
    for name in ["CH4", "H2O", "CH3", "CH2(S)", "O2", "N2"] {
        let species = CANONICAL_SPECIES.iter().position(|&s| s == name).unwrap();
        let h = concentrations[species] * 2e-4;
        let mut plus = concentrations.clone();
        plus[species] += h;
        let mut minus = concentrations.clone();
        minus[species] -= h;
        mechanism
            .evaluate(Temperature::new(t), &plus, &mut work, || false)
            .unwrap();
        let a = work.reaction_rates().to_vec();
        mechanism
            .evaluate(Temperature::new(t), &minus, &mut work, || false)
            .unwrap();
        for (j, (&a, &b)) in a.iter().zip(work.reaction_rates()).enumerate() {
            let numeric = (a - b) / (2.0 * h);
            let analytic = dc[j * 53 + species];
            close(
                numeric,
                analytic,
                2e-5,
                (a.abs() + b.abs()) * 2e-10 / h + 1e-20,
            );
        }
    }
    let h = 0.02;
    mechanism
        .evaluate(Temperature::new(t + h), &concentrations, &mut work, || {
            false
        })
        .unwrap();
    let plus = work.reaction_rates().to_vec();
    mechanism
        .evaluate(Temperature::new(t - h), &concentrations, &mut work, || {
            false
        })
        .unwrap();
    for (j, (&a, &b)) in plus.iter().zip(work.reaction_rates()).enumerate() {
        close(
            (a - b) / (2.0 * h),
            dt[j],
            2e-5,
            (a.abs() + b.abs()) * 2e-10 / h + 1e-20,
        );
    }
}

#[test]
fn pure_inert_has_no_manufactured_radicals_and_cancellation_is_reusable() {
    let mechanism = cantera_gri30(Gri30Calorics::ContinuousAt300K).unwrap();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let mut c = vec![0.0; 53];
    c[0] = 40.0;
    let mut polls = 0;
    assert_eq!(
        mechanism.evaluate(Temperature::new(300.0), &c, &mut work, || {
            polls += 1;
            polls == 80
        }),
        Err(KineticsErrorV1::Cancelled)
    );
    let report = mechanism
        .evaluate(Temperature::new(300.0), &c, &mut work, || false)
        .unwrap();
    assert!(work.reaction_rates().iter().all(|&x| x == 0.0));
    assert!(work.species_rates().iter().all(|&x| x == 0.0));
    assert_eq!(report.entropy_production_w_m3_k, 0.0);
    assert_eq!(report.adiabatic_temperature_rate_k_s, 0.0);
}
