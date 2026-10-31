use super::*;
use crate::electron_collision::{AVOGADRO_PER_MOL, ELECTRON_MASS_KG, MaxwellianRateControlsV1, MaxwellianRateWeightingV1, nist};

fn impact_bank() -> IdealGasEquilibriumV1 {
    let me = ELECTRON_MASS_KG * AVOGADRO_PER_MOL;
    let ids = vec![id("H"), id("H+"), id("e")];
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![vec![1, 1, 0]]).unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        vec![source("H", 0.001 + me, 2.5, 0.0, 0.0),
             source("H+", 0.001, 2.5, 13.598_434_599_702 * 1.602_176_634e-19 * AVOGADRO_PER_MOL, 0.0),
             source("e", me, 2.5, 0.0, 0.0)],
    ).unwrap()
}
fn impact() -> ElectronImpactRateV1 {
    ElectronImpactRateV1::new(id("e"), nist::hydrogen(id("H")),
        MaxwellianRateControlsV1 { relative_tolerance: 1e-9, maximum_subdivisions: 4096 }, 0.1).unwrap()
}
fn channel() -> ElementaryReactionV1 {
    ElementaryReactionV1::new(ReactionId::new("H-electron-impact").unwrap(),
        vec![(id("H"), 1), (id("e"), 1)], vec![(id("H+"), 1), (id("e"), 2)], impact()).unwrap()
}

#[test]
fn source_cross_section_reaches_existing_rates_derivatives_and_reverse_balance() {
    let mechanism = IdealGasKineticsV1::new(impact_bank(), "source-H-ionization", 1, vec![channel()], 1e-12).unwrap();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let t = Temperature::new(5000.0);
    let concentrations = [2.0, 0.0, 0.01];
    let evaluated = mechanism.evaluate(t, &concentrations, &mut work, || false).unwrap();
    let rate_source = impact();
    let integral = rate_source.cross_section().evaluate_properties(t, MaxwellianRateWeightingV1::Event, rate_source.controls(), || false).unwrap();
    let expected = AVOGADRO_PER_MOL * integral.coefficient_m3_s().unwrap() * concentrations[0] * concentrations[2];
    close(work.reaction_rates()[0], expected, 2e-12 * expected.abs());
    let expected_derivative = AVOGADRO_PER_MOL * integral.derivative_m3_s_k().unwrap() * concentrations[0] * concentrations[2];
    close(work.temperature_derivatives()[0], expected_derivative, 3e-12 * expected_derivative.abs());
    assert!(evaluated.adiabatic_temperature_rate_k_s < 0.0);
    assert_eq!(work.species_rates()[1], work.species_rates()[2]);
    assert_eq!(work.species_rates()[0], -work.species_rates()[1]);
    let equilibrium = fs_math::det::exp(0.5 * work.log_equilibrium_constants()[0]);
    mechanism.evaluate(t, &[1.0, equilibrium, equilibrium], &mut work, || false).unwrap();
    assert!(work.reaction_rates()[0].abs() <= expected * 1e-10);
    mechanism.evaluate(t, &[1.0, 0.0, 0.0], &mut work, || false).unwrap();
    assert_eq!(work.species_rates(), &[0.0, 0.0, 0.0]);
}

#[test]
fn electron_impact_refuses_wrong_electron_collider_and_nonbinary_orders() {
    let invalid_collider = channel().with_pressure_dependence(PressureDependenceV1::ThirdBody(ThirdBodyEfficienciesV1::default()));
    assert!(matches!(IdealGasKineticsV1::new(impact_bank(), "bad-collider", 1, vec![invalid_collider], 1e-12), Err(KineticsErrorV1::Invalid(_))));
    let wrong_electron = ElectronImpactRateV1::new(id("H+"), nist::hydrogen(id("H")), impact().controls(), 0.1).unwrap();
    let invalid = ElementaryReactionV1::new(ReactionId::new("bad-electron").unwrap(), vec![(id("H"),1), (id("e"),1)], vec![(id("H+"),1), (id("e"),2)], wrong_electron).unwrap();
    assert!(matches!(IdealGasKineticsV1::new(impact_bank(), "bad-electron", 1, vec![invalid], 1e-12), Err(KineticsErrorV1::Invalid(_))));
    let invalid = ElementaryReactionV1::new(ReactionId::new("bad-order").unwrap(), vec![(id("H"),2), (id("e"),1)], vec![(id("H+"),2), (id("e"),3)], impact()).unwrap();
    assert!(matches!(IdealGasKineticsV1::new(impact_bank(), "bad-order", 1, vec![invalid], 1e-12), Err(KineticsErrorV1::Invalid(_))));
}
