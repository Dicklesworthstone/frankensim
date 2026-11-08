//! G1/G3: original heat sources drive actual exact-volume thermal deformation.

use super::*;
use fs_couple::finite_solid::thermal::thermodynamics::heating::plane_stress::IncompressiblePlaneCaloricState;
use fs_material::hyper::plane_stress::PlaneStressSettings;

const A: [f64; 4] = [1.22, 0.16, -0.05, 0.91];

fn plane_state(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
) -> IncompressiblePlaneCaloricState {
    let material = material(card);
    let coefficients = coefficients(card, &material);
    material
        .incompressible_plane_stress_caloric_state(
            A,
            ReferenceDeformation::identity(),
            2.0,
            &coefficients,
            curve,
            CaloricReference::StressFreeZeroPressure,
            PlaneStressSettings::default(),
            || false,
        )
        .unwrap()
}

// Independent constrained scalar potential: substitute F22=s^3/det(A) into
// the NH energy before differentiating. The source E/nu and engineering stretch
// are written explicitly; no native material, derivative or inverse is used.
struct Oracle {
    internal_energy: f64,
    entropy: f64,
    mechanical_capacity: f64,
    pressure: f64,
    piola: [f64; 4],
    thickness: f64,
    jacobian: f64,
}

fn oracle(t: f64) -> Oracle {
    let e = 3000.0 - 12.0 * (t - 300.0);
    let n = 1.25 + 0.004 * (t - 300.0);
    let s = 1.0 + 0.001 * (t - 300.0);
    let mu = e / (2.0 * n);
    let mu_t = -12.0 / (2.0 * n) - mu * 0.004 / n;
    let mu_tt = -2.0 * 0.004 * mu_t / n;
    let d = A[0] * A[3] - A[1] * A[2];
    let cof = [A[3], -A[2], -A[1], A[0]];
    let i1: f64 = A.iter().map(|v| v * v).sum();
    let h = 0.5 * (s * i1 + s.powi(7) / d.powi(2) - 3.0 * s.powi(3));
    let h_t = 0.0005 * (i1 + 7.0 * s.powi(6) / d.powi(2) - 9.0 * s.powi(2));
    let h_tt = 0.0000005 * (42.0 * s.powi(5) / d.powi(2) - 18.0 * s);
    let entropy = -mu_t * h - mu * h_t;
    Oracle {
        internal_energy: mu * h + t * entropy,
        entropy,
        mechanical_capacity: -t * (mu_tt * h + 2.0 * mu_t * h_t + mu * h_tt),
        pressure: mu * (s.powi(4) / d.powi(2) - 1.0),
        piola: core::array::from_fn(|i| mu * (s * A[i] - s.powi(7) * cof[i] / d.powi(3))),
        thickness: s.powi(3) / d,
        jacobian: s.powi(3),
    }
}

fn plane_change(a: f64, b: f64) -> (f64, f64) {
    let first = oracle(a);
    let last = oracle(b);
    let (ha, sa) = base(a);
    let (hb, sb) = base(b);
    (
        2.0 * (hb - ha) + (last.internal_energy - first.internal_energy) / 3.0,
        2.0 * (sb - sa) + (last.entropy - first.entropy) / 3.0,
    )
}

#[test]
fn heat_changes_actual_thickness_pressure_and_stress_with_original_mass() {
    let card = card([-12.0, 0.004], 1.0);
    let curve = source_curve(&card, 0.0);
    let anchor = plane_state(&card, &curve);
    let initial = anchor.at_temperature(296.0, || false).unwrap();
    let (heat, entropy) = plane_change(296.0, 304.0);
    assert!(heat > 0.0);
    let transaction = initial.add_heat(heat, controls(), || false).unwrap();
    let final_state = transaction.final_state();
    let expected = oracle(304.0);
    near(final_state.temperature_k(), 304.0);
    near(transaction.internal_energy_change_j(), heat);
    near(transaction.entropy_change_j_per_k(), entropy);
    assert!(transaction.energy_residual_j().abs() <= controls().thermal.energy_tolerance_j);
    assert!(transaction.state_evaluations() > 9 && transaction.state_evaluations() <= 256);
    assert_eq!(
        final_state.storage_basis(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(final_state.reference_mass_kg(), 2.0);
    assert_eq!(
        final_state.reference_volume_m3(),
        anchor.reference_volume_m3()
    );
    near(final_state.current_volume_m3(), expected.jacobian / 3.0);
    near(final_state.current_density_kg_m3(), 6.0 / expected.jacobian);
    near(
        final_state.current_volume_m3() * final_state.current_density_kg_m3(),
        2.0,
    );
    let response = final_state.thermodynamics().response();
    near(response.deformation_gradient[8], expected.thickness);
    near(response.pressure_pa, expected.pressure);
    for (i, component) in [0, 1, 3, 4].into_iter().enumerate() {
        assert_eq!(response.deformation_gradient[component], A[i]);
        near(response.full_piola[component], expected.piola[i]);
    }
    for component in [2, 5, 8] {
        assert!(response.full_piola[component].abs() < 1e-9);
    }
    near(
        final_state.heat_capacity_j_per_k(),
        2.0 * 13.6 + expected.mechanical_capacity / 3.0,
    );
    let epsilon = 0.002;
    assert!(
        (final_state.heat_capacity_j_per_k()
            - plane_change(304.0 - epsilon, 304.0 + epsilon).0 / (2.0 * epsilon))
            .abs()
            < 2e-7
    );
    assert_ne!(
        response.deformation_gradient[8],
        initial.thermodynamics().response().deformation_gradient[8]
    );
    assert_ne!(
        response.full_piola,
        initial.thermodynamics().response().full_piola
    );
    let cooled = final_state.add_heat(-heat, controls(), || false).unwrap();
    near(cooled.final_state().temperature_k(), 296.0);
    assert!((transaction.entropy_change_j_per_k() + cooled.entropy_change_j_per_k()).abs() < 1e-9);
    assert_eq!(initial.temperature_k(), 296.0);
    assert_eq!(anchor.temperature_k(), 300.0);
}

#[test]
fn sub_ulp_caloric_heat_keeps_the_same_constrained_endpoint_and_original_sources() {
    let card = card([-12.0, 0.004], 1.0);
    let plain = plane_state(&card, &source_curve(&card, 0.0));
    let shifted = plane_state(&card, &source_curve(&card, 1e17));
    let target = plain.at_temperature(300.01, || false).unwrap();
    let heat = plain.internal_energy_change_to(&target).unwrap();
    assert!(heat > 0.0 && heat < 2.0 * (1e17_f64.next_up() - 1e17));
    let mut options = controls();
    options.thermal.bracket = [300.0, 300.1];
    options.maximum_temperature_step_k = 0.02;
    let a = plain.add_heat(heat, options, || false).unwrap();
    let b = shifted.add_heat(heat, options, || false).unwrap();
    assert_eq!(
        a.final_state().temperature_k(),
        b.final_state().temperature_k()
    );
    near(b.final_state().temperature_k(), 300.01);
    assert_eq!(a.entropy_change_j_per_k(), b.entropy_change_j_per_k());
    assert_eq!(a.energy_residual_j(), b.energy_residual_j());
    assert_eq!(
        a.final_state().thermodynamics().response().full_piola,
        b.final_state().thermodynamics().response().full_piola,
    );
    assert_eq!(b.final_state().coefficients(), shifted.coefficients());
    assert_eq!(b.final_state().material().law(), shifted.material().law());
    assert_ne!(
        plain.caloric_source().curve().identity(),
        shifted.caloric_source().curve().identity()
    );
}

#[test]
fn zero_expansion_and_constant_coefficients_reduce_to_original_caloric_integrals() {
    let card = card_with_expansion([0.0; 2], 1.0, [0.0; 3]);
    let curve = source_curve(&card, 0.0);
    let initial = plane_state(&card, &curve);
    let endpoint = initial.at_temperature(304.0, || false).unwrap();
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
    assert_eq!(
        endpoint.thermodynamics().response().deformation_gradient,
        initial.thermodynamics().response().deformation_gradient
    );
    assert_eq!(
        endpoint.thermodynamics().response().full_piola,
        initial.thermodynamics().response().full_piola
    );
    let heat = initial.internal_energy_change_to(&endpoint).unwrap();
    near(
        initial
            .add_heat(heat, controls(), || false)
            .unwrap()
            .final_state()
            .temperature_k(),
        304.0,
    );
}

#[test]
fn malformed_constraints_capacity_source_exits_and_late_cancellation_publish_nothing() {
    let card = card([-12.0, 0.004], 1.0);
    let curve = source_curve(&card, 0.0);
    let initial = plane_state(&card, &curve);
    let original_f = initial.thermodynamics().response().deformation_gradient;
    let original_p = initial.thermodynamics().response().full_piola;
    let (heat, _) = plane_change(300.0, 304.0);
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
    assert!(initial.at_temperature(311.0, || false).is_err());
    assert!(initial.add_heat(1e9, controls(), || false).is_err());
    assert!(
        initial
            .internal_energy_change_to(&plane_state(&card, &curve))
            .is_err()
    );
    let zero = initial.add_heat(0.0, controls(), || false).unwrap();
    assert_eq!(zero.state_evaluations(), 0);
    assert_eq!(
        zero.final_state()
            .thermodynamics()
            .response()
            .deformation_gradient,
        original_f
    );
    let material = material(&card);
    let coefficients = coefficients(&card, &material);
    assert!(
        material
            .incompressible_plane_stress_caloric_state(
                [1.0, 2.0, 2.0, 1.0],
                ReferenceDeformation::identity(),
                2.0,
                &coefficients,
                &curve,
                CaloricReference::StressFreeZeroPressure,
                PlaneStressSettings::default(),
                || false,
            )
            .is_err()
    );
    let cold_capacity = card_with_expansion([-12.0, 0.004], 1e-6, [-0.01, 0.0, 0.01]);
    let soft_material = super::super::material(&cold_capacity);
    let soft_coefficients = super::super::coefficients(&cold_capacity, &soft_material);
    assert!(
        soft_material
            .incompressible_plane_stress_caloric_state(
                A,
                ReferenceDeformation::identity(),
                2.0,
                &soft_coefficients,
                &source_curve(&cold_capacity, 0.0),
                CaloricReference::StressFreeZeroPressure,
                PlaneStressSettings::default(),
                || false,
            )
            .is_err()
    );
    assert_eq!(
        initial.thermodynamics().response().deformation_gradient,
        original_f
    );
    assert_eq!(initial.thermodynamics().response().full_piola, original_p);
}
