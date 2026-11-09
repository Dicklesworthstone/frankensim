//! Native Fourier adapters consume the original amount and caloric owners.
use super::*;

#[test]
fn amount_bound_caloric_law_retains_original_conversion_and_full_bank_identity() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let reference = GasEnthalpyLaw::from_mole(
        &bank,
        mole(vec![0.75, 0.25]),
        point(300.0, 100000.0),
        controls(),
    )
    .unwrap();
    for scale in [1e-280, 1.0, 1e280] {
        let amounts = [fs_qty::Amount::new(3.0 * scale), fs_qty::Amount::new(scale)];
        let law =
            GasEnthalpyLaw::from_mole_amounts(&bank, &amounts, point(300.0, 100000.0), controls())
                .unwrap();
        assert_eq!(law, reference);
        assert_eq!(law.model(), &bank.with_mole_amounts(&amounts).unwrap());
    }
    let pure = GasEnthalpyLaw::from_mole_amounts(
        &bank,
        &[fs_qty::Amount::new(0.0), fs_qty::Amount::new(2.0)],
        point(300.0, 100000.0),
        controls(),
    )
    .unwrap();
    assert_eq!(pure.model().thermodynamics().components().len(), 1);
    assert_eq!(pure.mole_composition().fractions(), &[0.0, 1.0]);
    let descriptor = pure.advected_composition().unwrap();
    assert_eq!(descriptor.species.len(), 2);
    assert_eq!(descriptor.mass_fractions, &[0.0, 1.0]);
    assert_eq!(
        descriptor.enthalpy_reference,
        reference.advected_composition().unwrap().enthalpy_reference
    );
    for values in [[0.0, 0.0], [-1.0, 2.0], [f64::from_bits(1), f64::MAX]] {
        assert!(
            GasEnthalpyLaw::from_mole_amounts(
                &bank,
                &values.map(fs_qty::Amount::new),
                point(300.0, 100000.0),
                controls(),
            )
            .is_err()
        );
    }
}

#[test]
fn fused_caloric_query_keeps_original_curved_cp_bounds_and_cancellation() {
    let bank = model(vec![0.5, 0.5], 0.001);
    // Source transport is supported only to 1e7 Pa. This caloric-only context
    // must never turn a successful inverse into a complete material state.
    let law = GasEnthalpyLaw::from_mole_amounts(
        &bank,
        &[fs_qty::Amount::new(3.0), fs_qty::Amount::new(1.0)],
        point(300.0, 1e9),
        controls(),
    )
    .unwrap();
    let molar_mass = 0.75 * 0.028 + 0.25 * 0.004;
    let exact_h =
        |t: f64| UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (3.25 * t + 0.000375 * t * t) / molar_mass;
    let bounds = law.enthalpy_bounds_j_kg().unwrap();
    for (actual, t) in bounds.into_iter().zip([280.0, 360.0]) {
        assert!((actual - exact_h(t)).abs() < 1e-8);
    }
    for t in [280.0, 300.0, 337.0, 360.0] {
        let state = law.caloric_at_enthalpy(exact_h(t), || false).unwrap();
        let cp = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (3.25 + 0.00075 * t) / molar_mass;
        assert!((state.temperature_k - t).abs() < 1e-9);
        assert!((state.heat_capacity_j_kg_k - cp).abs() < 1e-8);
        assert!(state.residual_j_kg.abs() <= controls().energy_tolerance_j_kg);
        assert!(state.iterations <= controls().maximum_iterations);
        assert!(law.resolve_enthalpy(exact_h(t), || false).is_err());
    }
    assert!(matches!(
        law.caloric_at_enthalpy(exact_h(320.0), || true),
        Err(GasEnergyError::Cancelled)
    ));
    assert!(law.caloric_at_enthalpy(bounds[0] - 1.0, || false).is_err());
    assert!(law.caloric_at_enthalpy(bounds[1] + 1.0, || false).is_err());
    assert!(law.caloric_at_enthalpy(f64::NAN, || false).is_err());
}
