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
