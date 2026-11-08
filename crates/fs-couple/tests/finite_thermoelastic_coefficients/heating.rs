//! G1/G3 finite heat: actual original coefficients, expansion, Cp and waves.

#[path = "heating/reservoir.rs"]
mod reservoir;

#[path = "heating/plane_stress.rs"]
mod plane_stress;

use super::*;
use fs_couple::finite_solid::thermal::thermodynamics::heating::{
    FiniteSolidCaloricState, FiniteThermoelasticHeatError, FiniteThermoelasticHeatSolve,
};
use fs_material::hyper::thermoelastic::FiniteThermoelasticState;
use fs_thermal::enthalpy::ThermalStorageBasis;
use fs_thermochem::solution::caloric::temperature::TemperatureControl;

fn controls() -> FiniteThermoelasticHeatSolve {
    FiniteThermoelasticHeatSolve {
        thermal: TemperatureControl {
            bracket: [291.0, 309.0],
            temperature_tolerance_k: 1e-9,
            energy_tolerance_j: 1e-8,
            maximum_iterations: 128,
        },
        maximum_temperature_step_k: 2.0,
        maximum_state_evaluations: 256,
        entropy_tolerance_j_per_k: 1e-9,
    }
}

fn source_curve(card: &MaterialCard, datum: f64) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new_exact(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(310.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        datum,
    )
    .unwrap()
}

fn prepared(card: &MaterialCard, caloric: &ResolvedHeatCapacityCurve) -> FiniteSolidCaloricState {
    let material = material(card);
    let coefficients = coefficients(card, &material);
    material
        .fixed_deformation_caloric_state(
            F,
            ReferenceDeformation::identity(),
            2.0,
            &coefficients,
            caloric,
            CaloricReference::StressFreeZeroPressure,
            || false,
        )
        .unwrap()
}

// Independent scalar invariant potential for this fixture. No native energy,
// temperature derivative, chart integration or solver is used by the oracle.
fn mechanical(t: f64) -> (f64, f64) {
    let young = 3000.0 - 12.0 * (t - 300.0);
    let nu = 0.25 + 0.004 * (t - 300.0);
    let s = 1.0 + 0.001 * (t - 300.0);
    let mu = young / (2.0 * (1.0 + nu));
    let mu_t = -12.0 / (2.0 * (1.0 + nu)) - mu * 0.004 / (1.0 + nu);
    let ratio = nu / (1.0 - 2.0 * nu);
    let ratio_t = 0.004 / (1.0 - 2.0 * nu).powi(2);
    let lambda = 2.0 * mu * ratio;
    let lambda_t = 2.0 * (mu_t * ratio + mu * ratio_t);
    let det_f = F[0] * (F[4] * F[8] - F[5] * F[7]) - F[1] * (F[3] * F[8] - F[5] * F[6])
        + F[2] * (F[3] * F[7] - F[4] * F[6]);
    let i1 = F.iter().map(|x| x * x).sum::<f64>() / (s * s);
    let log_j = (det_f / s.powi(3)).ln();
    let i1_t = -2.0 * 0.001 * i1 / s;
    let log_j_t = -3.0 * 0.001 / s;
    let w = 0.5 * mu * (i1 - 3.0) - mu * log_j + 0.5 * lambda * log_j * log_j;
    let w_t = 0.5 * mu_t * (i1 - 3.0) + 0.5 * mu * i1_t - mu_t * log_j - mu * log_j_t
        + 0.5 * lambda_t * log_j * log_j
        + lambda * log_j * log_j_t;
    let psi = s.powi(3) * w;
    let entropy = -(3.0 * s * s * 0.001 * w + s.powi(3) * w_t);
    (psi + t * entropy, entropy)
}
fn base(t: f64) -> (f64, f64) {
    let slope = if t < 300.0 { 0.2 } else { 0.4 };
    let dt = t - 300.0;
    (
        12.0 * dt + 0.5 * slope * dt * dt,
        slope * dt + (12.0 - 300.0 * slope) * (t / 300.0).ln(),
    )
}
fn expected_change(a: f64, b: f64) -> (f64, f64) {
    let (ua, sa) = mechanical(a);
    let (ub, sb) = mechanical(b);
    let (ha, ca) = base(a);
    let (hb, cb) = base(b);
    (
        2.0 * (hb - ha) + (ub - ua) / 3.0,
        2.0 * (cb - ca) + (sb - sa) / 3.0,
    )
}

#[test]
fn finite_heat_crosses_original_cp_knot_updates_stress_and_cools_conservatively() {
    let card = card([-12.0, 0.004], 1.0);
    let curve = source_curve(&card, 0.0);
    let anchor = prepared(&card, &curve);
    let initial = anchor.at_temperature(296.0, || false).unwrap();
    assert_eq!(initial.storage_basis(), ThermalStorageBasis::InternalEnergy);
    assert_eq!(initial.reference_mass_kg(), 2.0);
    near(initial.reference_volume_m3(), 1.0 / 3.0);
    let (heat, entropy) = expected_change(296.0, 304.0);
    let heated = initial.add_heat(heat, controls(), || false).unwrap();
    let final_state = heated.final_state();
    near(final_state.temperature_k(), 304.0);
    near(heated.internal_energy_change_j(), heat);
    near(heated.entropy_change_j_per_k(), entropy);
    assert!(heated.energy_residual_j().abs() <= 1e-8);
    assert!(heated.state_evaluations() > 9 && heated.state_evaluations() <= 256);
    near(final_state.base_specific_capacity_j_kg_k(), 13.6);
    let epsilon = 0.002;
    let derivative = expected_change(304.0 - epsilon, 304.0 + epsilon).0 / (2.0 * epsilon);
    assert!((final_state.heat_capacity_j_per_k() - derivative).abs() < 2e-7);
    assert_ne!(
        final_state.thermodynamics().first_piola_pa(),
        anchor.thermodynamics().first_piola_pa()
    );
    let waves = final_state.adiabatic_acoustics().unwrap();
    near(
        waves.fixed_deformation_heat_capacity_j_m3_k(),
        3.0 * final_state.heat_capacity_j_per_k(),
    );
    assert!(waves.acoustic().directional_waves([1.0, 0.0, 0.0]).is_ok());
    let cooled = final_state.add_heat(-heat, controls(), || false).unwrap();
    near(cooled.final_state().temperature_k(), 296.0);
    assert!((heated.entropy_change_j_per_k() + cooled.entropy_change_j_per_k()).abs() < 1e-9);
    assert_eq!(initial.temperature_k(), 296.0);
    assert_eq!(anchor.temperature_k(), 300.0);
}

#[test]
fn finite_heat_below_one_absolute_caloric_ulp_retains_original_temperature_and_entropy() {
    let card = card([-12.0, 0.004], 1.0);
    let plain = prepared(&card, &source_curve(&card, 0.0));
    let shifted = prepared(&card, &source_curve(&card, 1e17));
    let target = plain.at_temperature(300.01, || false).unwrap();
    let heat = plain.internal_energy_change_to(&target).unwrap();
    assert!(heat > 0.0 && heat < 2.0 * (1e17_f64.next_up() - 1e17));
    let mut settings = controls();
    settings.thermal.bracket = [300.0, 300.1];
    settings.maximum_temperature_step_k = 0.02;
    let a = plain.add_heat(heat, settings, || false).unwrap();
    let b = shifted.add_heat(heat, settings, || false).unwrap();
    assert_eq!(
        a.final_state().temperature_k(),
        b.final_state().temperature_k()
    );
    near(b.final_state().temperature_k(), 300.01);
    assert_eq!(a.entropy_change_j_per_k(), b.entropy_change_j_per_k());
    assert_eq!(a.energy_residual_j(), b.energy_residual_j());
    assert_eq!(
        a.final_state().thermodynamics(),
        b.final_state().thermodynamics()
    );
    assert_ne!(
        plain.caloric_source().curve().identity(),
        shifted.caloric_source().curve().identity()
    );
}

#[test]
fn constant_coefficients_retain_native_response_and_zero_expansion_reduces_to_original_cp() {
    for expansion in [[-0.01, 0.0, 0.01], [0.0; 3]] {
        let card = card_with_expansion([0.0; 2], 1.0, expansion);
        let curve = source_curve(&card, 0.0);
        let initial = prepared(&card, &curve);
        let endpoint = initial.at_temperature(304.0, || false).unwrap();
        let native = FiniteThermoelasticState::try_new(
            initial.material().law().law(),
            &F,
            304.0,
            *endpoint.thermodynamics().reference_jet(),
        )
        .unwrap();
        assert_eq!(endpoint.thermodynamics(), &native);
        if expansion == [0.0; 3] {
            assert_eq!(
                endpoint.heat_capacity_j_per_k(),
                2.0 * curve.secant_heat_capacity(304.0, 304.0).unwrap()
            );
            assert_eq!(
                initial.internal_energy_change_to(&endpoint).unwrap(),
                2.0 * curve.specific_enthalpy_change(300.0, 304.0).unwrap()
            );
            assert_eq!(
                initial.entropy_change_to(&endpoint).unwrap(),
                2.0 * curve.specific_entropy_change(300.0, 304.0).unwrap()
            );
        }
    }

    // Exercise the shared original finite-pair owner against a large mechanical
    // baseline. This deformation change is work, not a fixed-F heat request.
    // Its energy survives even when both absolute U values round identically.
    let law = fs_material::hyper::Hyperelastic::new(
        fs_material::hyper::HyperelasticModel::NeoHookean {
            mu: 1e20,
            lambda: 1e20,
        },
        4.0,
    )
    .unwrap();
    let before_f = [2.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0];
    let mut after_f = before_f;
    after_f[1] = 1e-8;
    let jet =
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [0.0; 9], [0.0; 9]).unwrap();
    let before = FiniteThermoelasticState::try_new(&law, &before_f, 300.0, jet).unwrap();
    let after = FiniteThermoelasticState::try_new(&law, &after_f, 310.0, jet).unwrap();
    assert_eq!(before.internal_energy_j_m3(), after.internal_energy_j_m3());
    near(
        before.internal_energy_increment_to(&after).unwrap(),
        0.5 * 1e20 * 1e-8 * 1e-8,
    );
    near(
        after.internal_energy_increment_to(&before).unwrap(),
        -0.5 * 1e20 * 1e-8 * 1e-8,
    );
}

#[test]
fn finite_heat_refuses_unresolved_sources_budgets_and_cancellation_without_mutation() {
    let card = card([-12.0, 0.004], 1.0);
    let curve = source_curve(&card, 0.0);
    let initial = prepared(&card, &curve);
    let original = initial.thermodynamics().clone();
    let (heat, _) = expected_change(300.0, 304.0);
    let mut short = controls();
    short.maximum_state_evaluations = 2;
    assert!(matches!(
        initial.add_heat(heat, short, || false),
        Err(FiniteThermoelasticHeatError::StateBudget { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        initial.add_heat(heat, controls(), || {
            polls += 1;
            polls > 8
        }),
        Err(FiniteThermoelasticHeatError::Cancelled)
    ));
    assert!(initial.add_heat(1e9, controls(), || false).is_err());
    assert!(initial.at_temperature(311.0, || false).is_err());
    let zero = initial.add_heat(0.0, controls(), || false).unwrap();
    assert_eq!(zero.state_evaluations(), 0);
    assert_eq!(zero.final_state().thermodynamics(), &original);
    assert!(
        initial
            .internal_energy_change_to(&prepared(&card, &curve))
            .is_err()
    );
    for card in [
        card_with_expansion([-12.0, 0.004], 1.0, [-0.01, 0.0, 0.02]),
        card_with_expansion([-12.0, 0.004], 1e-6, [-0.01, 0.0, 0.01]),
    ] {
        let material = material(&card);
        let coefficients = coefficients(&card, &material);
        assert!(
            material
                .fixed_deformation_caloric_state(
                    F,
                    ReferenceDeformation::identity(),
                    2.0,
                    &coefficients,
                    &source_curve(&card, 0.0),
                    CaloricReference::StressFreeZeroPressure,
                    || false
                )
                .is_err()
        );
    }
    assert_eq!(initial.thermodynamics(), &original);
}
