use super::*;

fn close(a: f64, b: f64, relative: f64) {
    assert!((a-b).abs() <= relative*b.abs(), "{a:e} != {b:e}");
}

#[test]
fn photon_statistics_match_the_bose_entropy_not_a_scaled_blackbody() {
    for x in [0.001_f64, 0.1, 1.0, 10.0, 100.0] {
        for emissivity in [1.0, 0.5, 0.01] {
            let n = emissivity / x.exp_m1();
            let direct = n.ln_1p() + n*(1.0/n).ln_1p();
            close(det::exp(log_mode_entropy(log_occupation(x, det::ln(emissivity)))), direct, 2e-12);
        }
    }
    let black = GrayPhotonEmission::new(1000.0, 1.0).unwrap().spectrum_nm(3000.0).unwrap();
    let gray = GrayPhotonEmission::new(1000.0, 0.1).unwrap().spectrum_nm(3000.0).unwrap();
    close(gray.energy_w_m2_sr_nm, 0.1*black.energy_w_m2_sr_nm, 1e-13);
    assert!(gray.entropy_w_m2_sr_k_nm > 0.1*black.entropy_w_m2_sr_k_nm);
    assert!(gray.entropy_w_m2_sr_k_nm < black.entropy_w_m2_sr_k_nm);
}

#[test]
fn black_spectral_integral_recovers_four_thirds_stefan_boltzmann() {
    for temperature in [300.0, 1000.0, 5778.0] {
        let emitter = GrayPhotonEmission::new(temperature, 1.0).unwrap();
        let lo = (C2_NM_K/(100.0*temperature)).ln();
        let hi = (C2_NM_K*1e6/temperature).ln();
        let n = 8192;
        let step = (hi-lo)/n as f64;
        let (mut energy, mut entropy) = (0.0, 0.0);
        for i in 0..=n {
            let wavelength = (lo+i as f64*step).exp();
            let s = emitter.spectrum_nm(wavelength).unwrap();
            let weight = if i == 0 || i == n { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
            energy += weight*wavelength*s.energy_w_m2_sr_nm;
            entropy += weight*wavelength*s.entropy_w_m2_sr_k_nm;
        }
        let flux = emitter.bolometric(PhotonEntropyQuadrature::default(), || false).unwrap();
        close(PI*step*energy/3.0, flux.energy_w_m2, 2e-9);
        close(PI*step*entropy/3.0, flux.entropy_w_m2_k, 2e-9);
        close(flux.entropy_w_m2_k, (4.0/3.0)*flux.energy_w_m2/temperature, 2e-14);
        assert_eq!(flux.evaluations, 0);
    }
}

#[test]
fn gray_integral_matches_independent_high_precision_references_and_temperature_scaling() {
    // Independent 70-digit quadrature in photon-energy x, NOT the production map.
    let references = [
        (0.1, 13.932907136474484),
        (0.001, 24.964363186199822),
        (1e-100, 561.9245528209042),
    ];
    for (emissivity, integral) in references {
        let emitter = GrayPhotonEmission::new(500.0, emissivity).unwrap();
        let a = emitter.bolometric(PhotonEntropyQuadrature::default(), || false).unwrap();
        let b = GrayPhotonEmission::new(1000.0, emissivity).unwrap()
            .bolometric(PhotonEntropyQuadrature::default(), || false).unwrap();
        let expected = emissivity * (4.0/3.0)*crate::radiation::STEFAN_BOLTZMANN_W_M2_K4
            * 500.0_f64.powi(3)*integral/BLACK_INTEGRAL;
        close(a.entropy_w_m2_k, expected, 2e-9);
        close(b.entropy_w_m2_k, 8.0*a.entropy_w_m2_k, 5e-13);
        close(b.energy_w_m2, 16.0*a.energy_w_m2, 5e-13);
        assert!(a.estimated_entropy_error_w_m2_k <= 1e-9*a.entropy_w_m2_k);
        assert!(a.evaluations <= PhotonEntropyQuadrature::default().maximum_evaluations);
        assert!(a.entropy_w_m2_k > (4.0/3.0)*a.energy_w_m2/500.0);
    }
}

#[test]
fn wien_and_rayleigh_jeans_limits_do_not_lose_representable_entropy() {
    // The occupation exp(-800) is zero in binary64, but its scaled spectra are not.
    let wavelength = 0.001;
    let emitter = GrayPhotonEmission::new(C2_NM_K/(800.0*wavelength), 1.0).unwrap();
    let spectrum = emitter.spectrum_nm(wavelength).unwrap();
    assert_eq!(det::exp(-800.0), 0.0);
    assert!(spectrum.energy_w_m2_sr_nm > 0.0);
    assert!(spectrum.entropy_w_m2_sr_k_nm > 0.0);
    // g(n) tends to log(n)+1 without subtracting O(n log n) terms.
    for log_n in [40.0, 700.0, 1400.0] {
        close(det::exp(log_mode_entropy(log_n)), log_n+1.0, 1e-13);
    }
    let weak = GrayPhotonEmission::new(1000.0, f64::from_bits(1)).unwrap()
        .bolometric(PhotonEntropyQuadrature::default(), || false).unwrap();
    assert!(weak.energy_w_m2 > 0.0 && weak.entropy_w_m2_k > 0.0);
}

#[test]
fn exact_zero_invalid_coordinates_budgets_and_cancellation_are_distinct() {
    let none = GrayPhotonEmission::new(300.0, 0.0).unwrap();
    assert_eq!(none.spectrum_nm(1000.0).unwrap().entropy_w_m2_sr_k_nm, 0.0);
    assert_eq!(none.bolometric(PhotonEntropyQuadrature::default(), || false).unwrap().evaluations, 0);
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(GrayPhotonEmission::new(t, 0.5), Err(RadiationEntropyError::Temperature));
    }
    for e in [-0.1, 1.1, f64::NAN] {
        assert_eq!(GrayPhotonEmission::new(300.0, e), Err(RadiationEntropyError::Emissivity));
    }
    assert_eq!(none.spectrum_nm(0.0), Err(RadiationEntropyError::Wavelength));
    let emitter = GrayPhotonEmission::new(300.0, 0.2).unwrap();
    assert_eq!(emitter.bolometric(PhotonEntropyQuadrature::default(), || true), Err(RadiationEntropyError::Cancelled));
    let config = PhotonEntropyQuadrature { relative_tolerance: 1e-12, maximum_evaluations: 33 };
    assert!(matches!(emitter.bolometric(config, || false), Err(RadiationEntropyError::BudgetExceeded { evaluations: 33 })));
    assert!(matches!(emitter.bolometric(PhotonEntropyQuadrature { maximum_evaluations: 0, ..config }, || false), Err(RadiationEntropyError::Quadrature)));
    let mut polls = 0;
    assert_eq!(emitter.bolometric(PhotonEntropyQuadrature::default(), || { polls += 1; polls > 4 }), Err(RadiationEntropyError::Cancelled));
    assert!(emitter.bolometric(PhotonEntropyQuadrature::default(), || false).is_ok());
}
