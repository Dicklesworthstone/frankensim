//! Conditional saturation reuses the original finite-stock trajectories.

use super::*;
use crate::SpeciesId;

fn model(saturation: f64, particles: bool) -> SolubilityKineticsV1 {
    let law = IdealSolubilityV1::new(
        SpeciesId::new("solvent").unwrap(),
        SpeciesId::new("solute").unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(saturation),
        0.0,
        [Temperature::new(280.0), Temperature::new(350.0)],
    )
    .unwrap();
    let kinetics = SolubilityKineticsV1::new(
        law,
        VolumetricFlowRate::new(0.4),
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap();
    if particles {
        kinetics.with_particle_surface(Amount::new(1.0)).unwrap()
    } else {
        kinetics
    }
}

#[test]
fn conditional_targets_reuse_both_native_trajectories_without_rebinding_sources() {
    for particles in [false, true] {
        let original = model(10.0, particles);
        let before = original.clone();
        for target in [0.3, 2.0, 8.0] {
            let separately_authored = model(target, particles);
            for (dissolved, solid) in [(0.0, 5.0), (1.0, 4.0), (4.0, 1.0), (4.0, 0.0)] {
                let state = SolubilityKineticStateV1::from_phases(
                    Amount::new(dissolved),
                    Amount::new(solid),
                )
                .unwrap();
                let actual = original
                    .advance_at_saturation(
                        Temperature::new(300.0),
                        Volume::new(1.0),
                        &state,
                        Time::new(2.0),
                        AmountConcentration::new(target),
                    )
                    .unwrap();
                let expected = separately_authored
                    .advance(
                        Temperature::new(300.0),
                        Volume::new(1.0),
                        &state,
                        Time::new(2.0),
                    )
                    .unwrap();
                assert_eq!(actual, expected);
                assert_eq!(
                    actual.equilibrium().saturation_concentration.value(),
                    target
                );
                assert_eq!(original, before);
            }
        }
    }
}

#[test]
fn conditional_saturation_keeps_original_temperature_policy_and_range_guards() {
    let original = model(10.0, true);
    let state = SolubilityKineticStateV1::from_phases(Amount::new(1.0), Amount::new(1.0)).unwrap();
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            original
                .advance_at_saturation(
                    Temperature::new(300.0),
                    Volume::new(1.0),
                    &state,
                    Time::new(0.0),
                    AmountConcentration::new(invalid),
                )
                .is_err()
        );
    }
    assert_eq!(
        original.advance_at_saturation(
            Temperature::new(351.0),
            Volume::new(1.0),
            &state,
            Time::new(1.0),
            AmountConcentration::new(2.0),
        ),
        Err(SolubilityErrorV1::OutsideTemperatureDomain)
    );
    let no_crystals =
        SolubilityKineticStateV1::from_phases(Amount::new(4.0), Amount::new(0.0)).unwrap();
    let blocked = original
        .advance_at_saturation(
            Temperature::new(300.0),
            Volume::new(1.0),
            &no_crystals,
            Time::new(100.0),
            AmountConcentration::new(2.0),
        )
        .unwrap();
    assert!(blocked.precipitation_blocked());
    assert_eq!(blocked.state(), &no_crystals);
}
