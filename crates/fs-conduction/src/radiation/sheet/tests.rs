use super::*;

fn sheet(a: [f64; 2], t: f64) -> ReciprocalSheetOptics {
    ReciprocalSheetOptics::new(a, t, [1.0, 2_000.0], "analytic diffuse-sheet fixture").unwrap()
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn transparent_sheet_transmits_irradiation_without_heating() {
    let report = sheet([0.0, 0.0], 1.0)
        .evaluate_gray(400.0, [100.0, 900.0])
        .unwrap();
    assert_eq!(report.outgoing_w_m2, [900.0, 100.0]);
    assert_eq!(report.net_heat_w_m2, [0.0, 0.0]);
    assert_eq!(report.absorbed_w_m2, [0.0, 0.0]);
    assert_eq!(report.emitted_w_m2, [0.0, 0.0]);
    assert_eq!(report.net_heat_temperature_derivative_w_m2_k, [0.0, 0.0]);
}

#[test]
fn perfect_reflector_keeps_the_two_irradiations_on_their_original_sides() {
    let report = sheet([0.0, 0.0], 0.0)
        .evaluate_gray(400.0, [100.0, 900.0])
        .unwrap();
    assert_eq!(report.outgoing_w_m2, [100.0, 900.0]);
    assert_eq!(report.net_heat_w_m2, [0.0, 0.0]);
}

#[test]
fn asymmetric_absorption_conserves_the_paired_energy_accounts() {
    let optics = sheet([0.2, 0.5], 0.3);
    let input = [1_400.0, 70.0];
    let report = optics.evaluate_gray(325.0, input).unwrap();
    let a = optics.absorptance();
    let r = optics.reflectance();
    for i in 0..2 {
        close(
            report.net_heat_w_m2[i],
            a[i] * (input[i] - report.blackbody_w_m2),
            2e-15,
        );
        close(
            report.outgoing_w_m2[i],
            a[i] * report.blackbody_w_m2 + r[i] * input[i] + 0.3 * input[1 - i],
            2e-15,
        );
    }
    close(
        input.iter().sum::<f64>() - report.outgoing_w_m2.iter().sum::<f64>(),
        report.net_heat_w_m2.iter().sum::<f64>(),
        2e-14,
    );
    assert_ne!(
        input[0] - report.outgoing_w_m2[0],
        report.net_heat_w_m2[0],
        "one face's net optical flux includes transmitted power",
    );
}

#[test]
fn equilibrium_preserves_blackbody_radiance_with_zero_sheet_heat() {
    let optics = sheet([0.2, 0.6], 0.25);
    let b = optics
        .evaluate_gray(400.0, [0.0; 2])
        .unwrap()
        .blackbody_w_m2;
    let report = optics.evaluate_gray(400.0, [b; 2]).unwrap();
    assert_eq!(report.net_heat_w_m2, [0.0; 2]);
    for outgoing in report.outgoing_w_m2 {
        close(outgoing, b, 2e-15);
    }
}

#[test]
fn opaque_limit_is_two_independent_kirchhoff_surfaces() {
    let optics = sheet([0.2, 0.8], 0.0);
    let input = [1_200.0, 200.0];
    let report = optics.evaluate_gray(350.0, input).unwrap();
    assert_eq!(report.transmitted_w_m2, [0.0; 2]);
    for i in 0..2 {
        close(
            input[i] - report.outgoing_w_m2[i],
            report.net_heat_w_m2[i],
            2e-14,
        );
    }
}

#[test]
fn finite_band_retains_planck_coverage_and_analytic_tangent() {
    let optics = sheet([0.2, 0.5], 0.3);
    let band = BlackbodyBand::new(3e-6, 15e-6).unwrap();
    let input = [500.0, 200.0];
    let temperature = 700.0;
    let report = optics.evaluate_band(band, temperature, input).unwrap();
    let emission = band.evaluate(temperature).unwrap();
    assert_eq!(report.wavelengths_m, Some([3e-6, 15e-6]));
    assert_eq!(report.blackbody_w_m2, emission.power_w_m2);
    assert_eq!(
        report.omitted_blackbody_fraction,
        emission.short_wave_fraction + emission.long_wave_fraction,
    );
    assert!(report.omitted_blackbody_fraction > 0.01);
    let increment = 1e-3;
    let plus = optics
        .evaluate_band(band, temperature + increment, input)
        .unwrap();
    let minus = optics
        .evaluate_band(band, temperature - increment, input)
        .unwrap();
    for i in 0..2 {
        assert_eq!(
            report.net_heat_temperature_derivative_w_m2_k[i],
            -optics.absorptance()[i] * emission.temperature_derivative_w_m2_k
        );
        close(
            report.net_heat_temperature_derivative_w_m2_k[i],
            (plus.net_heat_w_m2[i] - minus.net_heat_w_m2[i]) / (2.0 * increment),
            2e-8,
        );
    }
}

#[test]
fn temperature_tangent_is_cooling_and_net_heat_keeps_tiny_departures() {
    let optics = sheet([0.125, 0.5], 0.25);
    let b = optics
        .evaluate_gray(300.0, [0.0; 2])
        .unwrap()
        .blackbody_w_m2;
    let input = [f64::from_bits(b.to_bits() + 1), b];
    let report = optics.evaluate_gray(300.0, input).unwrap();
    assert_eq!(report.net_heat_w_m2[0], 0.125 * (input[0] - b));
    assert!(report.net_heat_w_m2[0] > 0.0);
    assert_eq!(report.net_heat_w_m2[1], 0.0);
    for i in 0..2 {
        close(
            report.net_heat_temperature_derivative_w_m2_k[i],
            -4.0 * optics.absorptance()[i] * STEFAN_BOLTZMANN_W_M2_K4 * 300.0_f64.powi(3),
            2e-15,
        );
    }
}

#[test]
fn source_partition_and_support_refuse_invalid_physics() {
    for (a, t) in [
        ([-0.1, 0.0], 0.0),
        ([0.5, 0.8], 0.3),
        ([0.0, 0.0], 1.1),
        ([f64::NAN, 0.0], 0.0),
        ([0.0, 0.0], f64::INFINITY),
        // A rounded sum of 1 must not hide a strictly negative reflectance.
        ([f64::EPSILON / 4.0, 0.0], 1.0),
    ] {
        assert!(ReciprocalSheetOptics::new(a, t, [100.0, 1_000.0], "source").is_err());
    }
    assert!(ReciprocalSheetOptics::new([0.2; 2], 0.2, [0.0, 1_000.0], "source").is_err());
    assert!(ReciprocalSheetOptics::new([0.2; 2], 0.2, [100.0, 100.0], "source").is_err());
    assert!(ReciprocalSheetOptics::new([0.2; 2], 0.2, [100.0, 1_000.0], " ").is_err());
    let optics = sheet([0.2; 2], 0.2);
    for t in [0.0, 2_001.0, f64::NAN, f64::INFINITY] {
        assert!(optics.evaluate_gray(t, [0.0; 2]).is_err());
    }
    for input in [[-1.0, 0.0], [f64::NAN, 0.0], [f64::INFINITY, 0.0]] {
        assert!(optics.evaluate_gray(300.0, input).is_err());
    }
}
