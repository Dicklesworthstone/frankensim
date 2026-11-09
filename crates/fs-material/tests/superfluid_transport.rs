//! G0/G1/G3: reciprocal bulk viscosity, independent constitutive work and heat.
//! Coefficients below are manufactured fixtures, not fitted helium data.

use fs_material::superfluid::transport::{LongitudinalTransport, SecondViscosities};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};

fn material(viscosity: f64, expansion: f64) -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: expansion,
        normal_viscosity_pa_s: viscosity,
    })
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 64.0 * f64::EPSILON * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g0_onsager_admission_checks_exact_products_without_overflow_or_tolerance() {
    assert_eq!(
        SecondViscosities::try_new(0.0, 0.0, 0.0).unwrap(),
        SecondViscosities::ZERO
    );
    for magnitude in [2.0_f64.powi(-900), 1.0, 2.0_f64.powi(900)] {
        for sign in [-1.0, 1.0] {
            let law = SecondViscosities::try_new(sign * magnitude, magnitude, magnitude).unwrap();
            assert_eq!(law.zeta1_m2_s(), sign * magnitude);
            assert!(
                SecondViscosities::try_new(sign * magnitude.next_up(), magnitude, magnitude)
                    .is_err()
            );
        }
    }
    // The rounded product is exactly one although the exact determinant is
    // negative by EPSILON². A multiplication-based test would accept it.
    let above = 1.0 + f64::EPSILON;
    let below = 1.0 - f64::EPSILON;
    assert_eq!(above * below, 1.0);
    assert!(SecondViscosities::try_new(1.0, above, below).is_err());
    let subnormal = f64::from_bits(1);
    assert!(SecondViscosities::try_new(subnormal, subnormal, subnormal).is_ok());
    assert!(SecondViscosities::try_new(subnormal.next_up(), subnormal, subnormal).is_err());
    assert!(SecondViscosities::try_new(subnormal, 0.0, 1.0).is_err());
    assert!(SecondViscosities::try_new(subnormal, 1.0, 0.0).is_err());
    for (cross, normal, relative) in [
        (0.0, -1.0, 1.0),
        (0.0, 1.0, -1.0),
        (f64::NAN, 1.0, 1.0),
        (0.0, f64::INFINITY, 1.0),
        (0.0, 1.0, f64::NAN),
    ] {
        assert!(SecondViscosities::try_new(cross, normal, relative).is_err());
    }
}

#[test]
fn g1_stress_chemical_potential_and_heat_close_independent_constitutive_work() {
    let second = SecondViscosities::try_new(0.2, 1.2, 0.8).unwrap();
    let state = material(0.3, 0.0);
    let transport = LongitudinalTransport::try_new(state, second, 5.0).unwrap();
    assert_eq!(transport.material(), &state);
    assert_eq!(transport.second_viscosities(), second);
    assert_eq!(transport.thermal_conductivity_w_m_k(), 5.0);
    let response = transport.response(-0.4, 1.2, 0.3).unwrap();
    close(response.longitudinal_viscous_stress_pa, -0.4);
    close(response.chemical_potential_correction_j_kg, -0.88);
    close(response.viscous_heat_w_m3, 1.216);
    close(
        response.viscous_heat_w_m3,
        response.longitudinal_viscous_stress_pa * -0.4
            - response.chemical_potential_correction_j_kg * 1.2,
    );
    close(response.conductive_heat_flux_w_m2, -1.5);
    close(response.conductive_entropy_flux_w_m2_k, -0.75);
    close(response.conductive_availability_loss_w_m3, 0.225);
    close(response.entropy_production_w_m3_k, (1.216 + 0.225) / 2.0);
}

#[test]
fn g1_energy_coordinate_diffusion_retains_both_densities_and_counterflow_sign() {
    let transport = LongitudinalTransport::try_new(
        material(0.3, 0.0),
        SecondViscosities::try_new(0.2, 1.2, 0.8).unwrap(),
        5.0,
    )
    .unwrap();
    let resistance = transport.longitudinal_resistance();
    close(resistance[0][0], 1.6);
    close(resistance[0][1], 0.6);
    assert_eq!(resistance[0][1], resistance[1][0]);
    close(resistance[1][1], 7.2);
    let diffusion = transport.kinetic_diffusivity_m2_s();
    close(diffusion[0][0], 0.4);
    close(diffusion[0][1], 0.6 / 3.0_f64.sqrt());
    assert_eq!(diffusion[0][1], diffusion[1][0]);
    close(diffusion[1][1], 9.6);
    // Arbitrary barycentric/counterflow gradients, evaluated independently in
    // normal velocity and relative MASS-current coordinates.
    for (v, w) in [(0.3, -0.7), (-0.2, 0.6), (1.0, 0.0), (0.0, 1.0)] {
        let response = transport.response(v + 0.75 * w, -3.0 * w, 0.0).unwrap();
        let from_resistance = v * (resistance[0][0] * v + resistance[0][1] * w)
            + w * (resistance[1][0] * v + resistance[1][1] * w);
        close(response.viscous_heat_w_m3, from_resistance);
        let y = [2.0 * v, 0.75_f64.sqrt() * w];
        close(
            response.viscous_heat_w_m3,
            y[0] * (diffusion[0][0] * y[0] + diffusion[0][1] * y[1])
                + y[1] * (diffusion[1][0] * y[0] + diffusion[1][1] * y[1]),
        );
    }
}

#[test]
fn g0_saturated_second_viscosity_retains_its_exact_nondissipative_direction() {
    for sign in [-1.0, 1.0] {
        let law = LongitudinalTransport::try_new(
            material(0.0, 0.0),
            SecondViscosities::try_new(sign * 2.0, 4.0, 1.0).unwrap(),
            0.0,
        )
        .unwrap();
        let null = law.response(1.0, -sign * 2.0, 0.0).unwrap();
        assert_eq!(null.viscous_heat_w_m3, 0.0);
        assert_eq!(null.longitudinal_viscous_stress_pa, 0.0);
        assert_eq!(null.chemical_potential_correction_j_kg, 0.0);
        let dissipative = law.response(1.0, sign * 2.0, 0.0).unwrap();
        assert_eq!(dissipative.viscous_heat_w_m3, 16.0);
    }
}

#[test]
fn g1_conduction_uses_cv_and_destroys_availability_without_creating_heat() {
    for expansion in [-0.2, 0.0, 0.2] {
        let state = material(0.0, expansion);
        let law = LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 3.0).unwrap();
        close(
            law.thermal_diffusivity_m2_s(),
            3.0 / (state.density_kg_m3() * state.specific_heat_cv_j_kg_k()),
        );
        let response = law.response(0.0, 0.0, -0.8).unwrap();
        assert_eq!(response.viscous_heat_w_m3, 0.0);
        assert!(response.conductive_heat_flux_w_m2 > 0.0);
        close(response.conductive_availability_loss_w_m3, 3.0 * 0.64 / 2.0);
        close(response.entropy_production_w_m3_k, 3.0 * 0.64 / 4.0);
        if expansion != 0.0 {
            assert!(law.thermal_diffusivity_m2_s() > 3.0 / 64.0);
        }
    }
}

#[test]
fn g3_cross_sign_covariance_and_zero_gradient_limits() {
    let state = material(0.1, 0.0);
    let positive = LongitudinalTransport::try_new(
        state,
        SecondViscosities::try_new(0.7, 1.0, 1.0).unwrap(),
        0.2,
    )
    .unwrap();
    let negative = LongitudinalTransport::try_new(
        state,
        SecondViscosities::try_new(-0.7, 1.0, 1.0).unwrap(),
        0.2,
    )
    .unwrap();
    let a = positive.response(-0.3, 0.7, 0.5).unwrap();
    let b = negative.response(-0.3, -0.7, -0.5).unwrap();
    assert_eq!(a.viscous_heat_w_m3, b.viscous_heat_w_m3);
    assert_eq!(
        a.longitudinal_viscous_stress_pa,
        b.longitudinal_viscous_stress_pa
    );
    assert_eq!(
        a.chemical_potential_correction_j_kg,
        -b.chemical_potential_correction_j_kg
    );
    assert_eq!(a.conductive_heat_flux_w_m2, -b.conductive_heat_flux_w_m2);
    assert_eq!(a.entropy_production_w_m3_k, b.entropy_production_w_m3_k);
    let zero = positive.response(0.0, 0.0, 0.0).unwrap();
    assert_eq!(zero.entropy_production_w_m3_k, 0.0);
}

#[test]
fn g0_nonfinite_or_unrepresentable_transport_refuses() {
    let state = material(0.3, 0.0);
    for kappa in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(LongitudinalTransport::try_new(state, SecondViscosities::ZERO, kappa).is_err());
    }
    let model = LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 1.0).unwrap();
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(model.response(invalid, 0.0, 0.0).is_err());
        assert!(model.response(0.0, invalid, 0.0).is_err());
        assert!(model.response(0.0, 0.0, invalid).is_err());
    }
    assert!(model.response(f64::MAX, 0.0, 0.0).is_err());
    assert!(model.response(0.0, 0.0, f64::MAX).is_err());
    let excessive = SecondViscosities::try_new(0.0, f64::MAX, f64::MAX).unwrap();
    assert!(LongitudinalTransport::try_new(state, excessive, 0.0).is_err());
}
