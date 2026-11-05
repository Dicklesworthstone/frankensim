//! G0/G1/G3: causal thermoelastic storage, positive diffusion modes and the
//! original Maxwell history. These are analytic implementation checks, not
//! experimental calibration or a spatial/large-strain validity certificate.

use fs_material::visco::slab_ted::{
    SlabThermoelasticBending, slab_relaxation_response, slab_relaxation_shape,
};
use fs_material::visco::{GeneralizedMaxwell, ViscoError};
use fs_math::c64::C64;

fn specimen() -> SlabThermoelasticBending {
    SlabThermoelasticBending::new_coupled(
        [8e9, 0.0, 0.0, 0.0, 4e9, 0.0, 0.0, 0.0, 2e9],
        [1e-5, 1e-5, -0.5e-5],
        300.0,
        3e6,
        10.0,
        0.002,
    )
    .unwrap()
}

// Independent Fourier eigenfunction expansion of the thermal field. This
// uses neither the production hyperbolic quotient nor emitted Prony weights.
fn thermal_modes(wt: f64, count: usize) -> C64 {
    let mut re = 0.0;
    let mut im = 0.0;
    for n in 0..count {
        let odd = (2 * n + 1) as f64;
        let x = wt / (odd * odd);
        let weight = 96.0 / (core::f64::consts::PI.powi(4) * odd.powi(4));
        re += weight * x * x / (1.0 + x * x);
        im += weight * x / (1.0 + x * x);
    }
    C64::new(re, im)
}

#[test]
fn g1_complete_response_matches_independent_heat_modes_without_cancellation() {
    for wt in [1e-150, 1e-12, 1e-6, 0.01, 0.2, 0.203, 1.0, 10.0, 100.0, 1e4] {
        let expected = thermal_modes(wt, 16_384);
        let actual = slab_relaxation_response(wt, 1.0);
        assert!(
            (actual.re / expected.re - 1.0).abs() < 3e-13,
            "storage at omega*tau={wt}: {actual:?}, {expected:?}"
        );
        assert!(
            (actual.im / expected.im - 1.0).abs() < 3e-10,
            "loss at omega*tau={wt}: {actual:?}, {expected:?}"
        );
        assert_eq!(actual.im, slab_relaxation_shape(wt, 1.0));
        assert_eq!(slab_relaxation_response(-wt, 1.0), actual.conj());
        assert!(actual.re > 0.0 && actual.re < 1.0 && actual.im > 0.0);
    }
    assert_eq!(slab_relaxation_response(0.0, 0.3), C64::ZERO);
    for (omega, tau) in [(1e100, 1.0), (1e200, 1e110)] {
        let actual = slab_relaxation_response(omega, tau);
        let tail = (12.0 / core::f64::consts::PI.powi(2) / omega) / tau;
        assert_eq!(actual.re, 1.0);
        assert!(actual.im > 0.0);
        assert!((actual.im / tail - 1.0).abs() < 1e-14);
    }
    for (omega, tau) in [(f64::NAN, 1.0), (1.0, 0.0), (f64::INFINITY, 1.0)] {
        let invalid = slab_relaxation_response(omega, tau);
        assert!(invalid.re.is_nan() && invalid.im.is_nan());
    }
}

#[test]
fn g1_coupled_capacity_retains_linear_thermodynamic_feedback() {
    let q = [8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0];
    let alpha = [0.5, 0.25, 0.125];
    let (temperature, capacity, conductivity, h) = (2.0, 10.0, 3.0, 0.2);
    // Independent diagonal constitutive relation C_sigma-C_epsilon=T alpha.Q.alpha.
    let decrement = 2.0 * (8.0 * 0.5 * 0.5 + 4.0 * 0.25 * 0.25 + 2.0 * 0.125 * 0.125);
    let fixed_strain_capacity = capacity - decrement;
    let weak =
        SlabThermoelasticBending::new(q, alpha, temperature, capacity, conductivity, h).unwrap();
    let coupled =
        SlabThermoelasticBending::new_coupled(q, alpha, temperature, capacity, conductivity, h)
            .unwrap();
    assert!(
        (coupled.first_relaxation_time_s() / weak.first_relaxation_time_s()
            - fixed_strain_capacity / capacity)
            .abs()
            < 1e-14
    );
    for (a, b) in coupled
        .relaxation_matrix()
        .iter()
        .zip(weak.relaxation_matrix())
    {
        assert!((*a / b - capacity / fixed_strain_capacity).abs() < 1e-13);
    }
    for invalid_capacity in [decrement - 0.25, decrement, decrement.next_up()] {
        assert!(
            SlabThermoelasticBending::new_coupled(
                q,
                alpha,
                temperature,
                invalid_capacity,
                conductivity,
                h
            )
            .is_err()
        );
    }
    let zero =
        SlabThermoelasticBending::new(q, [0.0; 3], temperature, capacity, conductivity, h).unwrap();
    assert_eq!(
        SlabThermoelasticBending::new_coupled(q, [0.0; 3], temperature, capacity, conductivity, h)
            .unwrap(),
        zero
    );
}

#[test]
fn g3_complete_tensor_retains_signed_cross_terms_and_zero_thermal_channel() {
    let law = specimen();
    let increment = law
        .complex_bending_increment(law.first_relaxation_time_s().recip())
        .unwrap();
    let factor = law.relaxation_factor();
    let shape = thermal_modes(1.0, 4096);
    for i in 0..3 {
        for j in 0..3 {
            let expected = shape.scale(factor[i] * factor[j]);
            assert!((increment[3 * i + j] - expected).abs() < 1e-12 * expected.abs());
        }
    }
    assert!(increment[2].re < 0.0 && increment[2].im < 0.0);
    assert_eq!(factor[0] - 2.0 * factor[1], 0.0);
    assert_eq!(law.complex_bending_increment(0.0).unwrap(), [C64::ZERO; 9]);
    let zero = SlabThermoelasticBending::new_coupled(
        [8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0],
        [0.0; 3],
        300.0,
        3e6,
        10.0,
        0.002,
    )
    .unwrap();
    assert_eq!(
        zero.complex_bending_increment(10.0).unwrap(),
        [C64::ZERO; 9]
    );
    assert!(law.complex_bending_increment(f64::NAN).is_err());
}

#[test]
fn g1_positive_spectrum_error_bound_encloses_independent_diffusion_ladders() {
    let law = specimen();
    let tau = law.first_relaxation_time_s();
    for maximum_wt in [1e-4, 1.0, 1e3, 1e8] {
        let mut previous_terms = 0;
        for tolerance in [1e-3, 1e-5, 1e-7] {
            let band = (0.0, maximum_wt / tau);
            let spectrum = law.lower_to_prony(band, tolerance, 64).unwrap();
            assert_eq!(spectrum.band_rad_s(), band);
            assert!(spectrum.terms().len() >= previous_terms);
            previous_terms = spectrum.terms().len();
            assert!(spectrum.normalized_error_bound() <= tolerance);
            assert!(spectrum.normalized_coefficient_bound() > 0.0);
            assert!(spectrum.normalized_error_bound() >= spectrum.normalized_tail_bound());
            let weight_sum: f64 = spectrum.terms().iter().map(|t| t.weight()).sum();
            assert!(
                weight_sum < 1.0,
                "the omitted thermal tail must not be reassigned"
            );
            for pair in spectrum.terms().windows(2) {
                assert!(pair[0].weight() > pair[1].weight() && pair[1].weight() > 0.0);
                assert!(pair[0].relaxation_time_s() > pair[1].relaxation_time_s());
            }
            for index in 0..65 {
                let wt = maximum_wt * (index as f64 / 64.0).powi(2);
                let expected = thermal_modes(wt, 4096);
                let mut actual = C64::ZERO;
                for term in spectrum.terms() {
                    let x = (wt / tau) * term.relaxation_time_s();
                    actual = actual
                        + C64::new(x * x / (1.0 + x * x), x / (1.0 + x * x)).scale(term.weight());
                }
                assert!(
                    (actual - expected).abs() <= spectrum.normalized_error_bound(),
                    "wt={wt}, tolerance={tolerance}, terms={}, error={}, bound={}",
                    spectrum.terms().len(),
                    (actual - expected).abs(),
                    spectrum.normalized_error_bound()
                );
            }
        }
    }
}

#[test]
fn g1_diffusion_terms_use_the_existing_maxwell_energy_and_relaxation_history() {
    let law = specimen();
    let tau = law.first_relaxation_time_s();
    let spectrum = law.lower_to_prony((0.0, 100.0 / tau), 1e-6, 64).unwrap();
    let material = GeneralizedMaxwell::new(
        1.0,
        spectrum
            .terms()
            .iter()
            .map(|t| (t.weight(), t.relaxation_time_s()))
            .collect(),
    )
    .unwrap();
    let mut state = material.state();
    let mut previous_heat = 0.0;
    for index in 0..128 {
        let strain = (0.17 * index as f64).sin() * 0.03;
        material.step(&mut state, strain, 0.02 * tau);
        let stored = 0.5 * strain * strain
            + 0.5
                * state
                    .q
                    .iter()
                    .zip(&material.terms)
                    .map(|(q, (e, _))| e * q * q)
                    .sum::<f64>();
        assert!(state.dissipated >= previous_heat);
        assert!((state.work - stored - state.dissipated).abs() < 1e-14);
        previous_heat = state.dissipated;
    }
    let held = state.strain;
    let branch_before: Vec<_> = state.q.clone();
    material.step(&mut state, held, tau);
    for (index, term) in spectrum.terms().iter().enumerate() {
        let expected = branch_before[index] * (-tau / term.relaxation_time_s()).exp();
        assert!((state.q[index] - expected).abs() <= 1e-13 * branch_before[index].abs());
    }
}

#[test]
fn g0_spectral_limits_and_unresolved_capacity_refuse_without_fabricated_loss() {
    let law = specimen();
    for band in [
        (-1.0, 1.0),
        (1.0, 1.0),
        (0.0, f64::INFINITY),
        (f64::NAN, 1.0),
    ] {
        assert!(law.lower_to_prony(band, 1e-6, 64).is_err());
    }
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(law.lower_to_prony((0.0, 1.0), tolerance, 64).is_err());
    }
    for count in [0, 65] {
        assert!(law.lower_to_prony((0.0, 1.0), 1e-6, count).is_err());
    }
    assert!(matches!(
        law.lower_to_prony((0.0, 1e12), 1e-10, 1),
        Err(ViscoError::TruncationTolerance { terms: 1, .. })
    ));
    assert!(matches!(
        law.lower_to_prony((0.0, 1e12), 1e-16, 64),
        Err(ViscoError::TruncationTolerance { terms: 64, .. })
    ));
}

#[test]
fn g1_dimensioned_increment_survives_an_unrepresentable_unit_shape() {
    let q = [1e150, 0.0, 0.0, 0.0, 1e150, 0.0, 0.0, 0.0, 1e150];
    let alpha = [1e-75, 0.0, 0.0];
    let high = SlabThermoelasticBending::new_coupled(q, alpha, 1.0, 100.0, 1e-100, 1e50).unwrap();
    let omega = 1e200;
    let tau = high.first_relaxation_time_s();
    let stiffness = high.relaxation_matrix()[0];
    assert_eq!(slab_relaxation_response(omega, tau).im, 0.0);
    let actual = high.complex_bending_increment(omega).unwrap()[0];
    let expected = ((stiffness / omega) / tau) * (12.0 / core::f64::consts::PI.powi(2));
    assert!(actual.im > 0.0);
    assert!((actual.im / expected - 1.0).abs() < 1e-14);

    let low = SlabThermoelasticBending::new_coupled(q, alpha, 1.0, 100.0, 1e308, 1e50).unwrap();
    let tau = low.first_relaxation_time_s();
    let stiffness = low.relaxation_matrix()[0];
    // Unit storage underflows, but the actual dimensional storage survives.
    let omega = 1e-100;
    assert_eq!(slab_relaxation_response(omega, tau).re, 0.0);
    let actual = low.complex_bending_increment(omega).unwrap()[0];
    let expected =
        ((((stiffness * (17.0 * core::f64::consts::PI.powi(4) / 1680.0)) * omega) * tau) * omega)
            * tau;
    assert!(actual.re > 0.0);
    assert!((actual.re - expected).abs() <= 4.0 * f64::from_bits(1));

    // Even omega*tau itself may underflow without losing physical loss.
    let omega = 1e-200;
    assert_eq!(omega * tau, 0.0);
    let actual = low.complex_bending_increment(omega).unwrap()[0];
    let expected = ((stiffness * (core::f64::consts::PI.powi(2) / 10.0)) * omega) * tau;
    assert!(actual.im > 0.0);
    assert!((actual.im / expected - 1.0).abs() < 1e-14);
}

#[test]
fn g1_coupled_capacity_and_weak_stress_channel_survive_cancellation() {
    // 110-digit exact-binary64 arithmetic gives this fixed-strain capacity.
    // Ordinary products/subtraction bias it by about two percent despite an
    // established positive interval lower bound.
    let expected_capacity = 2.2272570940340408e-11;
    let law = SlabThermoelasticBending::new_coupled(
        [70e9, 0.0, 0.0, 0.0, 70e9, 0.0, 0.0, 0.0, 70e9],
        [2e-5, 0.0, 0.0],
        300.0,
        8400.000000000024,
        100.0,
        0.002,
    )
    .unwrap();
    let recovered = law.first_relaxation_time_s() * (core::f64::consts::PI / 0.002).powi(2) * 100.0;
    assert!((recovered / expected_capacity - 1.0).abs() < 1e-13);

    // Q is SPD, yet its first thermal channel is 1e16 + 1 - 1e16.
    // A plain three-term dot product loses the represented positive channel.
    let q = [1e16, 1.0, 1e16, 1.0, 1e16, 0.0, 1e16, 0.0, 2e16];
    let law =
        SlabThermoelasticBending::new_coupled(q, [1.0, 1.0, -1.0], 1.0, 1e18, 100.0, 0.1).unwrap();
    let scale = (1.0 / 9.8e17_f64).sqrt() * 0.1 * (0.1_f64 / 12.0).sqrt();
    let factor = law.relaxation_factor();
    assert!(factor[0] > 0.0);
    assert!((factor[0] / scale - 1.0).abs() < 1e-13);
}
