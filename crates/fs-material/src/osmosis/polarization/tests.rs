use super::{FilmPermeance, PolarizedMembraneTransport};
use crate::osmosis::{
    HydraulicPermeance, MembraneTransport, R, ReflectionCoefficient, SolutePermeance,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-13 * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}"
    );
}

fn membrane(lp: f64, sigma: f64, p: f64) -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(lp).unwrap(),
        vec![(
            ReflectionCoefficient::new(sigma).unwrap(),
            SolutePermeance::new(p).unwrap(),
        )],
    )
    .unwrap()
}

fn films(sigma: f64, p: f64, left: Option<f64>, right: Option<f64>) -> PolarizedMembraneTransport {
    PolarizedMembraneTransport::new(
        membrane(1.0, sigma, p),
        vec![left.map(|k| FilmPermeance::new(k).unwrap())],
        vec![right.map(|k| FilmPermeance::new(k).unwrap())],
    )
    .unwrap()
}

// Independent integrated advection-diffusion oracle; std transcendental math
// appears only in tests, not in the prepared deterministic constitutive law.
fn layer_flux(p: f64, u: f64, left: f64, right: f64) -> f64 {
    if u == 0.0 {
        p * (left - right)
    } else if p == 0.0 {
        u * if u > 0.0 { left } else { right }
    } else if u > 0.0 {
        u * (left - (-u / p).exp() * right) / -(-u / p).exp_m1()
    } else {
        u * (right - (u / p).exp() * left) / -(u / p).exp_m1()
    }
}

#[test]
fn g0_film_permeance_preserves_units_and_refuses_invalid_geometry() {
    close(
        FilmPermeance::from_diffusion_layer(1e-9, 1e-4)
            .unwrap()
            .value(),
        1e-5,
    );
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(FilmPermeance::new(invalid).is_err());
        assert!(FilmPermeance::from_diffusion_layer(invalid, 1.0).is_err());
        assert!(FilmPermeance::from_diffusion_layer(1.0, invalid).is_err());
    }
    assert!(FilmPermeance::from_diffusion_layer(1e300, 1e-300).is_err());
    assert!(FilmPermeance::from_diffusion_layer(1e-300, 1e300).is_err());
}

#[test]
fn g0_preparation_and_evaluation_validate_layout_and_domain() {
    assert!(PolarizedMembraneTransport::new(membrane(1.0, 0.5, 0.2), vec![], vec![None]).is_err());
    assert!(PolarizedMembraneTransport::new(membrane(1.0, 0.5, 0.2), vec![None], vec![]).is_err());
    let law = films(0.5, 0.2, Some(0.4), Some(0.6));
    assert!(law.coefficients(1, 0.0).is_err());
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(law.coefficients(0, invalid).is_err());
    }
    let coefficients = law.coefficients(0, 0.1).unwrap();
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            coefficients
                .membrane_face_concentrations(invalid, 1.0)
                .is_err()
        );
        assert!(
            coefficients
                .membrane_face_concentrations(1.0, invalid)
                .is_err()
        );
    }
}

#[test]
fn g0_absent_films_are_the_exact_original_interface() {
    let original = membrane(0.4, 0.35, 0.7);
    let law = PolarizedMembraneTransport::new(original.clone(), vec![None], vec![None]).unwrap();
    assert_eq!(law.membrane(), &original);
    for j in [-3.0, -1e-10, 0.0, 1e-10, 3.0] {
        let coefficients = law.coefficients(0, j).unwrap();
        assert_eq!(
            coefficients.donor_coefficients(),
            original.solute_stream_coefficients(0, j).unwrap()
        );
        assert_eq!(
            coefficients.membrane_face_concentrations(3.0, 7.0).unwrap(),
            (3.0, 7.0)
        );
    }
}

#[test]
fn g1_zero_flow_adds_diffusion_resistances_and_resolves_each_face() {
    let law = films(0.8, 0.5, Some(0.25), Some(1.0));
    let coefficients = law.coefficients(0, 0.0).unwrap();
    let (kl, kr) = coefficients.donor_coefficients();
    let effective = 1.0 / (1.0 / 0.25 + 1.0 / 0.5 + 1.0);
    close(kl, effective);
    close(kr, effective);
    let flux = effective * (12.0 - 2.0);
    let (left, right) = coefficients
        .membrane_face_concentrations(12.0, 2.0)
        .unwrap();
    close(left, 12.0 - flux / 0.25);
    close(right, 2.0 + flux);
    close(0.5 * (left - right), flux);
}

#[test]
fn g1_leaking_membrane_carries_identical_flux_through_all_three_layers() {
    for sigma in [0.0, 0.4, 1.0] {
        let law = films(sigma, 0.3, Some(0.4), Some(1.1));
        for j in [-1.2, -0.03, 0.0, 0.03, 0.8] {
            let coefficients = law.coefficients(0, j).unwrap();
            let (kl, kr) = coefficients.donor_coefficients();
            let (left, right) = coefficients.membrane_face_concentrations(3.0, 8.0).unwrap();
            let flux = kl * 3.0 - kr * 8.0;
            close(layer_flux(0.4, j, 3.0, left), flux);
            close(layer_flux(0.3, (1.0 - sigma) * j, left, right), flux);
            close(layer_flux(1.1, j, right, 8.0), flux);
        }
    }
}

#[test]
fn g3_orientation_reversal_swaps_faces_and_donor_coefficients() {
    let forward = films(0.65, 0.2, Some(0.3), Some(0.9));
    let reverse = films(0.65, 0.2, Some(0.9), Some(0.3));
    for j in [-1.0, -0.1, 0.0, 0.1, 1.0] {
        let a = forward.coefficients(0, j).unwrap();
        let b = reverse.coefficients(0, -j).unwrap();
        let (al, ar) = a.donor_coefficients();
        let (bl, br) = b.donor_coefficients();
        close(al, br);
        close(ar, bl);
        let (afl, afr) = a.membrane_face_concentrations(2.0, 9.0).unwrap();
        let (bfl, bfr) = b.membrane_face_concentrations(9.0, 2.0).unwrap();
        close(afl, bfr);
        close(afr, bfl);
    }
}

#[test]
fn g1_nonselective_uniform_solution_has_no_polarization() {
    let law = films(0.0, 0.7, Some(0.2), Some(1.3));
    for j in [-30.0, -0.1, 0.0, 0.1, 30.0] {
        let coefficients = law.coefficients(0, j).unwrap();
        let (left, right) = coefficients.membrane_face_concentrations(4.0, 4.0).unwrap();
        close(left, 4.0);
        close(right, 4.0);
        let (kl, kr) = coefficients.donor_coefficients();
        assert!((4.0 * (kl - kr) - 4.0 * j).abs() <= 2e-13 * (1.0 + j.abs()));
    }
}

#[test]
fn g1_retained_solute_matches_both_exponential_faces_without_epsilon_leakage() {
    let law = films(1.0, 0.0, Some(0.7), Some(1.4));
    for j in [-2.0, -0.1, 0.0, 0.1, 1.3] {
        let coefficients = law.coefficients(0, j).unwrap();
        assert_eq!(coefficients.donor_coefficients(), (0.0, 0.0));
        let (left, right) = coefficients.membrane_face_concentrations(2.3, 7.1).unwrap();
        close(left, 2.3 * (j / 0.7).exp());
        close(right, 7.1 * (-j / 1.4).exp());
    }
    // A sealed zero-flow membrane retains its own two reservoirs, even when
    // the surrounding diffusion films have enormous conductance.
    let sealed = films(0.0, 0.0, Some(1e100), Some(1e100));
    let coefficients = sealed.coefficients(0, 0.0).unwrap();
    assert_eq!(coefficients.donor_coefficients(), (0.0, 0.0));
    assert_eq!(
        coefficients.membrane_face_concentrations(1.0, 9.0).unwrap(),
        (1.0, 9.0)
    );
}

#[test]
fn g1_absent_bulk_solute_has_finite_positive_faces_and_conservative_flux() {
    let law = films(0.6, 0.4, Some(0.5), Some(0.8));
    for (bulk_left, bulk_right) in [(0.0, 9.0), (9.0, 0.0), (0.0, 0.0)] {
        for j in [-0.7, 0.0, 0.7] {
            let coefficients = law.coefficients(0, j).unwrap();
            let (left, right) = coefficients
                .membrane_face_concentrations(bulk_left, bulk_right)
                .unwrap();
            assert!(left.is_finite() && left >= 0.0 && right.is_finite() && right >= 0.0);
            let (kl, kr) = coefficients.donor_coefficients();
            let flux = kl * bulk_left - kr * bulk_right;
            close(layer_flux(0.5, j, bulk_left, left), flux);
            close(layer_flux(0.8, j, right, bulk_right), flux);
        }
    }
}

#[test]
fn g1_high_peclet_nonselective_transport_reaches_finite_donor_limit() {
    let law = films(0.0, 1e-100, Some(1e-100), Some(1e-100));
    for j in [-1e100, 1e100] {
        let coefficients = law.coefficients(0, j).unwrap();
        let (kl, kr) = coefficients.donor_coefficients();
        close(kl, j.max(0.0));
        close(kr, (-j).max(0.0));
        let (left, right) = coefficients.membrane_face_concentrations(2.0, 7.0).unwrap();
        let donor = if j > 0.0 { 2.0 } else { 7.0 };
        close(left, donor);
        close(right, donor);
    }
    // Logarithmic weights recover a representable face even when exp(Pe)
    // overflows. Exact absence stays zero; actual concentration overflow fails.
    let retained = films(1.0, 0.0, Some(1.0), Some(1.0));
    for j in [-1000.0, 1000.0] {
        let coefficients = retained.coefficients(0, j).unwrap();
        let tiny_bulk = if j > 0.0 {
            (1e-300, 0.0)
        } else {
            (0.0, 1e-300)
        };
        let faces = coefficients
            .membrane_face_concentrations(tiny_bulk.0, tiny_bulk.1)
            .unwrap();
        let expected = (1e-300_f64.ln() + 1000.0).exp();
        close(if j > 0.0 { faces.0 } else { faces.1 }, expected);
        assert_eq!(
            coefficients.membrane_face_concentrations(0.0, 0.0).unwrap(),
            (0.0, 0.0)
        );
        assert!(coefficients.membrane_face_concentrations(1.0, 1.0).is_err());
    }
}

#[test]
fn g1_series_transport_has_nonnegative_bulk_isothermal_dissipation() {
    let sigma = 0.7;
    let law = films(sigma, 0.3, Some(0.6), Some(1.5));
    for (bulk_left, bulk_right) in [(1.0_f64, 9.0_f64), (9.0, 1.0), (4.0, 4.0)] {
        for j in [-3.0, -0.1, 0.0, 0.1, 3.0] {
            let coefficients = law.coefficients(0, j).unwrap();
            let (left, right) = coefficients
                .membrane_face_concentrations(bulk_left, bulk_right)
                .unwrap();
            let (kl, kr) = coefficients.donor_coefficients();
            let solute = kl * bulk_left - kr * bulk_right;
            // Reduced units RT=Lp=1. Choose pressure to support this candidate
            // flow; subtract the BULK solvent chemical-potential difference.
            let pressure = j + sigma * (left - right);
            let dissipation =
                j * (pressure - (bulk_left - bulk_right)) + solute * (bulk_left / bulk_right).ln();
            assert!(
                dissipation >= -1e-12,
                "negative dissipation {dissipation:e}"
            );
        }
    }
}

#[test]
fn g1_retained_water_oracle_resolves_polarization_limited_flow() {
    for (lp_rt, bulk_left, expected) in [
        (1.0, 1.0, 0.723_654_785_457_764_8),
        (10.0, 1.0, 0.900_237_151_872_033_6),
        (1e6, 1.0, 0.924_195_996_234_962_6),
        (1.0, 0.0, 1.705_211_004_027_451),
    ] {
        let law = PolarizedMembraneTransport::new(
            membrane(lp_rt / (R * 300.0), 1.0, 0.0),
            vec![Some(FilmPermeance::new(1.0).unwrap())],
            vec![Some(FilmPermeance::new(2.0).unwrap())],
        )
        .unwrap();
        let (mut lo, mut hi) = (0.0, 4.0);
        for _ in 0..100 {
            let j = 0.5 * (lo + hi);
            let coefficients = law.coefficients(0, j).unwrap();
            let (left, right) = coefficients
                .membrane_face_concentrations(bulk_left, 4.0)
                .unwrap();
            let water = law
                .membrane()
                .water_flux(300.0, 0.0, 0.0, &[left], &[right])
                .unwrap();
            if j > water.volume_flux_m_s {
                hi = j;
            } else {
                lo = j;
            }
        }
        let j = 0.5 * (lo + hi);
        close(j, expected);
        assert!(j < lp_rt * (4.0 - bulk_left));
        if bulk_left > 0.0 {
            assert!(j < 4.0_f64.ln() / 1.5);
        }
    }
}
