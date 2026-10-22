use super::*;
use core::f64::consts::PI;

fn surface(t: f64, e: f64) -> GrayDiffuseThermalSurface {
    GrayDiffuseThermalSurface::try_new(t, e, ContentHash([17; 32])).unwrap()
}

#[test]
fn g1_si_planck_spectrum_integrates_to_the_original_cooling_power() {
    // Integrate the actual wavelength API using x=hc/(lambda*k*T), including
    // its Jacobian and the Lambertian hemisphere factor pi. No gray-body
    // power formula enters the numerical integrand.
    for temperature in [300.0, 1500.0, 6000.0] {
        let s = surface(temperature, 0.37);
        let n = 16_384;
        let dx = 50.0 / f64::from(n);
        let mut power = 0.0;
        for i in 0..n {
            let x = (f64::from(i) + 0.5) * dx;
            let lambda = C2_NM_K / (temperature * x);
            power +=
                PI * s.radiance_w_m2_sr_nm(lambda).unwrap() * C2_NM_K / (temperature * x * x) * dx;
        }
        let expected = s.bolometric_exitance_w_m2();
        assert!(
            (power / expected - 1.0).abs() < 2.0e-9,
            "T={temperature}, {power} vs {expected}"
        );
    }
    // Independent direct SI expression at one in-band wavelength.
    let lambda_m: f64 = 650.0e-9;
    let t = 1500.0;
    let si = 2.0 * PLANCK_J_S * LIGHT_M_S * LIGHT_M_S
        / (lambda_m.powi(5) * (PLANCK_J_S * LIGHT_M_S / (lambda_m * BOLTZMANN_J_K * t)).exp_m1())
        * 1.0e-9;
    assert!((surface(t, 1.0).radiance_w_m2_sr_nm(650.0).unwrap() / si - 1.0).abs() < 2.0e-13);
}

#[test]
fn g3_temperature_and_emissivity_change_physical_power_and_color() {
    let cold = surface(1000.0, 0.25);
    let hot = surface(2000.0, 0.25);
    let brighter = surface(2000.0, 0.75);
    assert_eq!(
        hot.bolometric_exitance_w_m2(),
        16.0 * cold.bolometric_exitance_w_m2()
    );
    for lambda in [450.0, 650.0, 10_000.0] {
        assert!(
            hot.radiance_w_m2_sr_nm(lambda).unwrap() > cold.radiance_w_m2_sr_nm(lambda).unwrap()
        );
        assert!(
            (brighter.radiance_w_m2_sr_nm(lambda).unwrap()
                / hot.radiance_w_m2_sr_nm(lambda).unwrap()
                - 3.0)
                .abs()
                < 1.0e-14
        );
    }
    let ratio = |s: GrayDiffuseThermalSurface| {
        s.radiance_w_m2_sr_nm(450.0).unwrap() / s.radiance_w_m2_sr_nm(650.0).unwrap()
    };
    assert!(ratio(hot) > ratio(cold));
    assert_eq!(surface(300.0, 0.0).radiance_w_m2_sr_nm(500.0), Ok(0.0));
    assert_eq!(surface(1.0, 1.0).radiance_w_m2_sr_nm(500.0), Ok(0.0));
}

#[test]
fn g0_bad_thermodynamic_inputs_refuse_before_scene_construction() {
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(GrayDiffuseThermalSurface::try_new(t, 0.5, ContentHash([1; 32])).is_err());
    }
    for e in [-0.1, 1.1, f64::NAN] {
        assert!(GrayDiffuseThermalSurface::try_new(300.0, e, ContentHash([1; 32])).is_err());
    }
    assert!(GrayDiffuseThermalSurface::try_new(300.0, 0.5, ContentHash([0; 32])).is_err());
    assert!(surface(300.0, 0.5).radiance_w_m2_sr_nm(0.0).is_err());
    assert_ne!(
        surface(300.0, 0.5).identity(),
        surface(301.0, 0.5).identity()
    );
}
