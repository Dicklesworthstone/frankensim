use super::*;
use crate::cavitation::{BubbleGas, BubbleLimits, KellerMiksisLimits};

fn liquid() -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}
fn bubble(radius: f64, ambient: f64, exponent: f64) -> KellerMiksisBubble {
    let material = liquid();
    KellerMiksisBubble::try_new(
        material,
        BubbleGas::equilibrium(material, radius, ambient, exponent).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.5 * radius,
                maximum_radius_m: 2.0 * radius,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
}
fn class(radius: f64, number: f64) -> BubblePopulationClass {
    BubblePopulationClass::try_new(bubble(radius, 100_000.0, 1.4), number).unwrap()
}
fn volume(radius: f64) -> f64 {
    4.0 * core::f64::consts::PI * radius.powi(3) / 3.0
}

#[test]
fn physical_number_and_volume_moments_keep_each_native_gas_and_interface() {
    let a = BubblePopulationClass::try_new(bubble(20e-6, 100_000.0, 1.0), 2e8).unwrap();
    let b = BubblePopulationClass::try_new(bubble(100e-6, 100_000.0, 1.4), 3e7).unwrap();
    let population = BubblePopulation::try_new(liquid(), vec![a, b], 2).unwrap();
    assert_eq!(population.number_density_per_m3(), 2.3e8);
    let expected = 2e8 * volume(20e-6) + 3e7 * volume(100e-6);
    assert!((population.void_fraction() / expected - 1.0).abs() < 4e-16);
    assert_eq!(population.classes(), &[a, b]);
    assert_ne!(
        a.bubble().gas().reference_pressure_pa(),
        b.bubble().gas().reference_pressure_pa()
    );
    assert_eq!(a.bubble().gas().polytropic_exponent(), 1.0);
    assert_eq!(b.bubble().gas().polytropic_exponent(), 1.4);
    assert_eq!(
        population.equilibrium_absolute_pressure_pa(),
        Some(100_000.0)
    );
    assert_eq!(population.occupied_classes(), 2);
    assert_eq!(population.maximum_classes(), 2);
    assert_eq!(
        a.linearization().unwrap(),
        a.bubble().linearization().unwrap()
    );
}

#[test]
fn splitting_identical_classes_preserves_number_and_occupied_volume() {
    let original = class(100e-6, 1e8);
    let merged = BubblePopulation::try_new(liquid(), vec![original], 1).unwrap();
    let split = BubblePopulation::try_new(liquid(), vec![class(100e-6, 1e8 / 8.0); 8], 8).unwrap();
    assert_eq!(
        merged.number_density_per_m3(),
        split.number_density_per_m3()
    );
    assert_eq!(merged.void_fraction(), split.void_fraction());
    assert_eq!(
        merged.equilibrium_absolute_pressure_pa(),
        split.equilibrium_absolute_pressure_pa()
    );
}

#[test]
fn compensated_number_moment_retains_small_classes_without_normalizing() {
    let population = BubblePopulation::try_new(
        liquid(),
        vec![class(1e-6, 1e16), class(1e-6, 1.0), class(1e-6, 1.0)],
        3,
    )
    .unwrap();
    assert_eq!(population.number_density_per_m3(), 10_000_000_000_000_002.0);
    assert_eq!(population.classes()[1].number_density_per_m3(), 1.0);
    assert_eq!(population.classes()[2].number_density_per_m3(), 1.0);
    assert!(population.void_fraction() > 0.04);
}

#[test]
fn empty_population_and_zero_weight_unstable_class_have_no_bubble_response() {
    let empty = BubblePopulation::try_new(liquid(), vec![], 0).unwrap();
    assert_eq!(empty.number_density_per_m3(), 0.0);
    assert_eq!(empty.void_fraction(), 0.0);
    assert_eq!(empty.equilibrium_absolute_pressure_pa(), None);
    let reference = bubble(100e-6, 100_000.0, 1.4);
    let unstable = KellerMiksisBubble::try_new(
        liquid(),
        BubbleGas::try_new(100e-6, 0.0, 1.0).unwrap(),
        reference.compressibility_limits(),
    )
    .unwrap();
    assert!(unstable.linearization().is_err());
    let absent = BubblePopulationClass::try_new(unstable, 0.0).unwrap();
    let population = BubblePopulation::try_new(liquid(), vec![absent], 1).unwrap();
    assert_eq!(population.occupied_classes(), 0);
    assert_eq!(population.void_fraction(), 0.0);
    assert_eq!(absent.linearization(), None);
}

#[test]
fn class_budget_liquid_identity_and_common_ambient_are_real_admission_limits() {
    let a = class(100e-6, 1e7);
    assert!(matches!(
        BubblePopulation::try_new(liquid(), vec![a, a], 1),
        Err(BubblePopulationError::ClassBudgetExceeded {
            requested: 2,
            maximum: 1
        })
    ));
    let changed = CavitationState::try_new(300.0, 999.0, 0.001, 0.072, 2000.0, 1500.0).unwrap();
    assert!(matches!(
        BubblePopulation::try_new(changed, vec![a], 1),
        Err(BubblePopulationError::MaterialMismatch { class_index: 0 })
    ));
    let incompatible = BubblePopulationClass::try_new(bubble(20e-6, 100_001.0, 1.4), 1e7).unwrap();
    assert!(matches!(
        BubblePopulation::try_new(liquid(), vec![a, incompatible], 2),
        Err(BubblePopulationError::EquilibriumMismatch { class_index: 1, .. })
    ));
}

#[test]
fn malformed_number_density_and_excess_volume_are_not_silently_clamped() {
    let model = bubble(100e-6, 100_000.0, 1.4);
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(BubblePopulationClass::try_new(model, invalid).is_err());
    }
    let overfilled = BubblePopulationClass::try_new(model, 1e12).unwrap();
    assert!(overfilled.void_fraction() > 1.0);
    assert!(BubblePopulation::try_new(liquid(), vec![overfilled], 1).is_err());
}
