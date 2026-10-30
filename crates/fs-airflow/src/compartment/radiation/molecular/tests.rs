use super::*;
use fs_gas::multicomponent::CaloricDomain;
use fs_qty::{Area, Length};
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, cantera_h2o2};

fn gas() -> MixtureGas {
    let source = cantera_h2o2(HydrogenCalorics::ContinuousAt300K).unwrap();
    MixtureGas::new(
        source.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(3500.0)),
        || false,
    )
    .unwrap()
}

fn room(gas: MixtureGas, t: f64, p: f64, water_mole_fraction: f64) -> TwoZone {
    let mut fractions = vec![0.0; gas.species_count()];
    let mass = water_mole_fraction * 0.018015 + (1.0 - water_mole_fraction) * 0.028014;
    fractions[3] = water_mole_fraction * 0.018015 / mass;
    fractions[6] = (1.0 - water_mole_fraction) * 0.028014 / mass;
    TwoZone::new(
        gas,
        Area::new(1.0),
        Length::new(1.0),
        Length::new(0.0),
        Pressure::new(p),
        [Temperature::new(t); 2],
        [fractions.clone(), fractions],
        || false,
    )
    .unwrap()
}

#[test]
fn published_values_units_piecewise_branch_and_domain_are_preserved() {
    // Independently summed decimal source values at 1000 K.
    for (species, expected) in [5.575471316, 27.3741, 2.422757, 3.311034]
        .into_iter()
        .enumerate()
    {
        assert!((fit(species, 1000.0).unwrap() - expected).abs() < 2e-12);
        for temperature in [300.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 2500.0] {
            assert!(fit(species, temperature).unwrap() >= 0.0);
        }
        for temperature in [299.999, 2500.001, f64::NAN] {
            assert!(matches!(
                fit(species, temperature),
                Err(RadiationError::MolecularTemperature { .. })
            ));
        }
    }
    let low = 4.7869 - 0.06953 * 750.0 + 2.95775e-4 * 750.0_f64.powi(2)
        - 4.25732e-7 * 750.0_f64.powi(3)
        + 2.02894e-10 * 750.0_f64.powi(4);
    assert!((fit(2, 750.0).unwrap() - low).abs() < 2e-13);
    let high_limit = polynomial(
        750.0,
        &[10.09, -0.01183, 4.7753e-6, -5.87209e-10, -2.5334e-14],
    );
    assert!((fit(2, 750.000001).unwrap() - high_limit).abs() < 1e-8);
    // The supplied piecewise fit is not silently joined or refitted.
    assert!((low - high_limit).abs() > 1e-4);
}

#[test]
fn original_moles_drive_partial_pressure_and_opacity_without_mass_fraction_substitution() {
    let gas = gas();
    let model =
        MolecularPlanck::new(&gas, "N2 and radical species treated as nonparticipating").unwrap();
    let base = model
        .evaluate_room(&room(gas.clone(), 1000.0, ATM_PA, 0.1), || false)
        .unwrap();
    assert!(base[0].is_none());
    let base = base[1].unwrap();
    assert!((base.partial_pressures[0].value() - 0.1 * ATM_PA).abs() < 2e-10);
    assert!((base.absorption.value() - 0.5575471316).abs() < 1e-13);
    let doubled = model
        .evaluate_room(&room(gas, 1000.0, 2.0 * ATM_PA, 0.1), || false)
        .unwrap()[1]
        .unwrap();
    assert!((doubled.absorption.value() - 2.0 * base.absorption.value()).abs() < 2e-13);
    let gray = base
        .gray_absorption("one local Planck mean held gray over the authored enclosure")
        .unwrap();
    assert_eq!(gray.coefficient(), base.absorption);
    assert!(gray.source().contains("TNF/RADCAL"));
    assert!(base.gray_absorption("").is_err());
}

#[test]
fn absence_domain_source_identity_and_cancellation_do_not_invent_opacity() {
    let gas = gas();
    assert!(MolecularPlanck::new(&gas, " ").is_err());
    let model =
        MolecularPlanck::new(&gas, "only the supplied four molecular species participate").unwrap();
    let inert = room(gas.clone(), 2800.0, ATM_PA, 0.0);
    assert_eq!(
        model.evaluate_room(&inert, || false).unwrap()[1]
            .unwrap()
            .absorption
            .value(),
        0.0
    );
    let participating = room(gas.clone(), 2800.0, ATM_PA, 0.1);
    assert!(matches!(
        model.evaluate_room(&participating, || false),
        Err(RadiationError::MolecularTemperature {
            temperature_k: 2800.0
        })
    ));
    assert_eq!(
        model.evaluate_room(&inert, || true),
        Err(RadiationError::Cancelled)
    );
    assert_eq!(
        model.evaluate_room(&inert, || false).unwrap()[1]
            .unwrap()
            .absorption
            .value(),
        0.0
    );
    let changed = MixtureGas::new(
        gas.source().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(3400.0)),
        || false,
    )
    .unwrap();
    assert!(
        model
            .evaluate_room(&room(changed, 1000.0, ATM_PA, 0.1), || false)
            .is_err()
    );
}
