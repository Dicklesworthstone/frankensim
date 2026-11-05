//! G0/G3 detailed balance, near-equilibrium accuracy and ideal-mixture limits.
use super::*;

fn symmetric_isomer(scale: f64) -> DeclaredStoichiometricMassActionV1 {
    let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(1);
    let equilibrium = IdealGasReactionEquilibriumV1::new(
        stoichiometric,
        certificate,
        reaction_id,
        vec![
            model("A", 0.01, 3.0, 0.0, 100_000.0),
            model("B", 0.01, 3.0, 0.0, 100_000.0),
        ],
    )
    .unwrap();
    DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(equilibrium, progress_rate_scale(scale)),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    )
}

fn relative_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_near_equilibrium_keeps_tiny_net_rates_and_nonnegative_entropy() {
    // Equal standard-state potentials give K=1 and the independent linear
    // rate law r=k*(a_A-a_B), including a tiny difference of huge gross rates.
    let delta = 2.0_f64.powi(-40);
    for scale in [1e-200, 1.0, 1e200] {
        let law = symmetric_isomer(scale);
        for sign in [-1.0, 1.0] {
            let (a, b) = if sign > 0.0 {
                (1.0 + delta, 1.0)
            } else {
                (1.0, 1.0 + delta)
            };
            let evaluation = law
                .evaluate(
                    Temperature::new(500.0),
                    vec![activity("A", a), activity("B", b)],
                )
                .unwrap();
            let net = evaluation.net_progress_rate().unwrap().value();
            relative_close(net, sign * scale * delta);
            let affinity = evaluation.dimensionless_affinity().unwrap().value();
            relative_close(affinity, sign * delta.ln_1p());
            let expected = (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * delta.ln_1p()) * (scale * delta);
            let entropy = evaluation.entropy_production_w_m3_k().unwrap().value();
            relative_close(entropy, expected);
            assert!(entropy > 0.0);
        }
    }
}

#[test]
fn g0_exact_detailed_balance_has_zero_progress_and_entropy() {
    let law = symmetric_isomer(1e200);
    for a in [1e-200, 1.0, 1e100] {
        let result = law
            .evaluate(
                Temperature::new(500.0),
                vec![activity("A", a), activity("B", a)],
            )
            .unwrap();
        assert_eq!(result.dimensionless_affinity().unwrap().value(), 0.0);
        assert_eq!(result.net_progress_rate().unwrap().value(), 0.0);
        assert_eq!(result.entropy_production_w_m3_k().unwrap().value(), 0.0);
    }
}

#[test]
fn g3_kinetic_multiplier_scales_progress_and_entropy_without_changing_affinity() {
    let law = symmetric_isomer(7.0);
    let activities = vec![activity("A", 1.0 + 2.0_f64.powi(-40)), activity("B", 1.0)];
    let base = law
        .evaluate(Temperature::new(500.0), activities.clone())
        .unwrap();
    let scaled = law
        .evaluate_with_kinetic_multiplier(
            Temperature::new(500.0),
            activities.clone(),
            Dimensionless::new(3.0),
        )
        .unwrap();
    assert_eq!(
        scaled.dimensionless_affinity(),
        base.dimensionless_affinity()
    );
    relative_close(
        scaled.net_progress_rate().unwrap().value(),
        3.0 * base.net_progress_rate().unwrap().value(),
    );
    relative_close(
        scaled.entropy_production_w_m3_k().unwrap().value(),
        3.0 * base.entropy_production_w_m3_k().unwrap().value(),
    );
    let stopped = law
        .evaluate_with_kinetic_multiplier(
            Temperature::new(500.0),
            activities,
            Dimensionless::new(0.0),
        )
        .unwrap();
    assert_eq!(
        stopped.dimensionless_affinity(),
        base.dimensionless_affinity()
    );
    assert_eq!(stopped.entropy_production_w_m3_k().unwrap().value(), 0.0);
}

#[test]
fn g3_exact_absence_retains_singular_affinity_without_inventing_trace_gas() {
    let law = symmetric_isomer(1.0);
    for (a, b, sign) in [(0.0, 1.0, -1.0), (1.0, 0.0, 1.0)] {
        let activities = vec![activity("A", a), activity("B", b)];
        let result = law
            .evaluate(Temperature::new(500.0), activities.clone())
            .unwrap();
        assert_eq!(result.net_progress_rate().unwrap().value(), sign);
        assert_eq!(
            result.dimensionless_affinity().unwrap().value(),
            sign * f64::INFINITY
        );
        assert_eq!(
            result.entropy_production_w_m3_k().unwrap().value(),
            f64::INFINITY
        );
        let stopped = law
            .evaluate_with_kinetic_multiplier(
                Temperature::new(500.0),
                activities,
                Dimensionless::new(0.0),
            )
            .unwrap();
        assert_eq!(stopped.entropy_production_w_m3_k().unwrap().value(), 0.0);
    }
    let both = law
        .evaluate(
            Temperature::new(500.0),
            vec![activity("A", 0.0), activity("B", 0.0)],
        )
        .unwrap();
    assert!(both.dimensionless_affinity().is_none());
    assert_eq!(both.entropy_production_w_m3_k().unwrap().value(), 0.0);
}

#[test]
fn g3_numerical_range_loss_is_not_physical_equilibrium() {
    let law = symmetric_isomer(1e-310);
    let result = law
        .evaluate(
            Temperature::new(500.0),
            vec![activity("A", 1.0 + 2.0_f64.powi(-48)), activity("B", 1.0)],
        )
        .unwrap();
    assert!(result.forward_progress_rate().is_some());
    assert!(result.reverse_progress_rate().is_some());
    assert!(result.dimensionless_affinity().unwrap().value() > 0.0);
    assert!(result.net_progress_rate().is_none());
    assert!(result.entropy_production_w_m3_k().is_none());
}

#[test]
fn g0_affinity_preserves_small_remainders_of_large_stoichiometric_logs() {
    assert_eq!(mass_action_affinity(1e16, 1.0, 1e16), Some(-1.0));
    assert_eq!(mass_action_affinity(-1e16, -1.0, -1e16), Some(1.0));
}

fn reference_shifted_water(
    reversed: bool,
    enthalpy_reference_r: [f64; 2],
    entropy_reference_r: [f64; 2],
    binding_r: f64,
) -> IdealGasReactionEquilibriumV1 {
    let models = [
        ("H2", 0.002_015_88, 3.0, 2.0, 0.0, 0.0, 1.0),
        ("O2", 0.031_998_8, 3.5, 0.0, 2.0, 0.0, 2.0),
        ("H2O", 0.018_015_28, 4.0, 2.0, 1.0, -binding_r, 4.0),
    ]
    .into_iter()
    .map(|(name, mass, cp_r, hydrogen, oxygen, formation_r, entropy_r)| {
        let mut coefficients = [0.0; 9];
        coefficients[2] = cp_r;
        coefficients[7] = hydrogen * enthalpy_reference_r[0]
            + oxygen * enthalpy_reference_r[1]
            + formation_r;
        coefficients[8] = hydrogen * entropy_reference_r[0]
            + oxygen * entropy_reference_r[1]
            + entropy_r;
        Nasa9StandardStateModelV1::new(
            species(name),
            MolarMass::new(mass),
            convention(100_000.0),
            vec![crate::Nasa9RegionV1::from_card(card(coefficients, 100_000.0)).unwrap()],
        )
        .unwrap()
    })
    .collect();
    let (stoichiometric, certificate, reaction_id) = water_stoichiometry(reversed);
    IdealGasReactionEquilibriumV1::new(stoichiometric, certificate, reaction_id, models).unwrap()
}

#[test]
fn g3_elemental_caloric_reference_shifts_preserve_equilibrium_and_reverse_rates() {
    // These dyadic shifts preserve the authored binding/entropy differences
    // exactly in the source coefficients. They must cancel under A*N=0 before
    // absolute species properties are rounded; lost input bits are not repaired.
    let h_shift = [2.0_f64.powi(50), -2.0_f64.powi(48)];
    let s_shift = [2.0_f64.powi(40), -2.0_f64.powi(39)];
    for t in [250.0, 800.0, 5000.0] {
        let temperature = Temperature::new(t);
        let baseline = ThermodynamicReverseRateClosureV1::new(
            reference_shifted_water(false, [0.0; 2], [0.0; 2], 1024.0),
            progress_rate_scale(7.0),
        )
        .evaluate(temperature)
        .unwrap();
        // Independent constant-Cp solution for 2 H2 + O2 -> 2 H2O.
        let h = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (-1.5 * t - 2048.0);
        let s = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (-1.5 * fs_math::det::ln(t) + 4.0);
        relative_close(baseline.equilibrium().standard_reaction_gibbs().value(), h - t * s);
        for (h_reference, s_reference) in [
            (h_shift, [0.0; 2]),
            ([0.0; 2], s_shift),
            (h_shift.map(|v| -v), s_shift.map(|v| -v)),
        ] {
            let closure = ThermodynamicReverseRateClosureV1::new(
                reference_shifted_water(false, h_reference, s_reference, 1024.0),
                progress_rate_scale(7.0),
            );
            let shifted = closure.evaluate(temperature).unwrap();
            assert_eq!(
                shifted.equilibrium().standard_reaction_gibbs(),
                baseline.equilibrium().standard_reaction_gibbs()
            );
            assert_eq!(
                shifted.equilibrium().log_equilibrium_constant(),
                baseline.equilibrium().log_equilibrium_constant()
            );
            assert_eq!(
                shifted.reverse_progress_rate_scale(),
                baseline.reverse_progress_rate_scale()
            );
            assert_eq!(shifted.receipt().equilibrium_receipt().evaluator_version(), 2);
            assert_eq!(
                shifted.receipt().equilibrium_receipt(),
                shifted.equilibrium().receipt()
            );
            for (term, receipt) in closure
                .equilibrium()
                .terms()
                .iter()
                .zip(shifted.equilibrium().receipt().terms())
            {
                let source = term.model().evaluate(temperature).unwrap();
                assert_eq!(receipt.standard_state_receipt(), source.receipt());
                assert_eq!(receipt.standard_gibbs_bits(), source.properties().g().value().to_bits());
                assert_eq!(
                    receipt.contribution_bits(),
                    ((term.coefficient() as f64) * source.properties().g().value()).to_bits()
                );
            }
            let reversed = reference_shifted_water(true, h_reference, s_reference, 1024.0)
                .evaluate(temperature)
                .unwrap();
            assert_eq!(
                reversed.standard_reaction_gibbs().value().to_bits(),
                (-shifted.equilibrium().standard_reaction_gibbs().value()).to_bits()
            );
            assert_eq!(
                reversed.log_equilibrium_constant().value().to_bits(),
                (-shifted.equilibrium().log_equilibrium_constant().value()).to_bits()
            );
        }
    }
}

#[test]
fn g0_elemental_reference_cancellation_retains_physical_binding_energy() {
    let h_reference = [2.0_f64.powi(50), -2.0_f64.powi(48)];
    let s_reference = [2.0_f64.powi(40), -2.0_f64.powi(39)];
    let temperature = Temperature::new(800.0);
    let original = reference_shifted_water(false, h_reference, s_reference, 1024.0)
        .evaluate(temperature)
        .unwrap();
    let stronger_binding = reference_shifted_water(false, h_reference, s_reference, 2048.0)
        .evaluate(temperature)
        .unwrap();
    relative_close(
        stronger_binding.standard_reaction_gibbs().value()
            - original.standard_reaction_gibbs().value(),
        -2048.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    );
    relative_close(
        stronger_binding.log_equilibrium_constant().value()
            - original.log_equilibrium_constant().value(),
        2048.0 / temperature.value(),
    );
}
