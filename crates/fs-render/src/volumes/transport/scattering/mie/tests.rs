//! G0/G1/G2 sphere-scattering checks against independent identities and values.

use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

pub(super) fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

pub(super) fn sphere_for_x(x: f64, n: f64, k: f64) -> MieSphere {
    MieSphere::try_new(x * 550e-9 / core::f64::consts::TAU, n, k, 1.0, 1000.0).unwrap()
}

pub(super) fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn g2_published_mie_coefficients_and_efficiencies() {
    // Independent values from Prahl's Mie Algorithms, a_n/b_n and Mie
    // scattering calculations sections. exp(-i omega t) matches these a_n.
    with_cx(|cx, _| {
        for (x, n, k, a, b) in [
            (
                50.0,
                4.0 / 3.0,
                0.0,
                C64::new(0.531105889295, -0.499031485631),
                C64::new(0.791924475935, -0.405931152229),
            ),
            (
                2.0,
                1.5,
                1.0,
                C64::new(0.546520203397, -0.152373857258),
                C64::new(0.389714727888, 0.227896075256),
            ),
            (
                2.0,
                1.1,
                25.0,
                C64::new(0.322406907480, -0.465063542971),
                C64::new(0.575167279092, 0.492912495262),
            ),
        ] {
            let result = sphere_for_x(x, n, k)
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap();
            assert!(
                (result.coefficients[0].0 - a).abs() < 1e-11,
                "a1 at x={x}: {:?}",
                result.coefficients[0].0
            );
            assert!(
                (result.coefficients[0].1 - b).abs() < 1e-11,
                "b1 at x={x}: {:?}",
                result.coefficients[0].1
            );
        }
        let radius = 0.525e-6;
        let result = MieSphere::try_new(radius, 1.55, 0.0, 1.0, 1000.0)
            .unwrap()
            .at_wavelength(cx, 632.8, MieBudget::default())
            .unwrap();
        let c = result.cross_sections();
        close(
            c.scattering_m2 / (core::f64::consts::PI * radius * radius),
            3.1054255,
            6e-8,
        );
        close(c.mean_cosine, 0.6331368, 6e-8);
        assert_eq!(c.absorption_m2, 0.0);
        assert_eq!(c.extinction_m2(), c.scattering_m2);
    });
}

#[test]
fn g1_rayleigh_limit_has_correct_wavelength_size_and_phase_laws() {
    with_cx(|cx, _| {
        let sphere = sphere_for_x(0.001, 1.5, 0.0);
        let result = sphere
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        // Independent electrostatic polarizability: C_sca = (8 pi/3) k^4 r^6 |(m²-1)/(m²+2)|².
        let wave = core::f64::consts::TAU / 550e-9;
        let contrast = (1.5f64.powi(2) - 1.0) / (1.5f64.powi(2) + 2.0);
        let want = 8.0 * core::f64::consts::PI / 3.0
            * wave.powi(4)
            * sphere.radius_m.powi(6)
            * contrast.powi(2);
        close(result.cross_sections().scattering_m2 / want, 1.0, 1e-6);
        for mu in [-1.0, -0.5, 0.0, 0.4, 1.0] {
            let want = 3.0 * (1.0 + mu * mu) / (16.0 * core::f64::consts::PI);
            close(result.eval_cos(mu).unwrap(), want, 2e-7);
        }
        let blue = sphere
            .at_wavelength(cx, 400.0, MieBudget::default())
            .unwrap();
        let red = sphere
            .at_wavelength(cx, 800.0, MieBudget::default())
            .unwrap();
        close(
            blue.cross_sections().scattering_m2 / red.cross_sections().scattering_m2,
            16.0,
            2e-5,
        );
        let larger = MieSphere::try_new(2.0 * sphere.radius_m, 1.5, 0.0, 1.0, 1000.0)
            .unwrap()
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        close(
            larger.cross_sections().scattering_m2 / result.cross_sections().scattering_m2,
            64.0,
            3e-5,
        );
        close(
            larger.cross_sections().mass_scattering_m2_kg
                / result.cross_sections().mass_scattering_m2_kg,
            8.0,
            4e-6,
        );
    });
}

#[test]
fn g0_absorption_optical_theorem_and_index_matching() {
    with_cx(|cx, _| {
        for (x, n, k) in [
            (0.001, 1.5, 0.1),
            (0.1, 1.5, 0.1),
            (2.0, 1.5, 1.0),
            (50.0, 1.33, 0.01),
            (128.0, 1.33, 0.0),
        ] {
            let sphere = sphere_for_x(x, n, k);
            let result = sphere
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap();
            let c = result.cross_sections();
            let optical_sum: f64 = result
                .coefficients
                .iter()
                .enumerate()
                .map(|(i, (a, b))| (2 * i + 3) as f64 * (a.re + b.re))
                .sum();
            let optical_cross =
                2.0 * optical_sum / (x * x) * core::f64::consts::PI * sphere.radius_m.powi(2);
            close(c.extinction_m2() / optical_cross, 1.0, 1e-10);
            assert!(c.absorption_m2 >= 0.0 && c.scattering_m2 > 0.0);
            assert!(c.transport_m2() > 0.0);
            if k > 0.0 {
                assert!(c.absorption_m2 > 0.0);
            }
        }
        let matched = MieSphere::try_new(1e-7, 1.33, 0.0, 1.33, 1000.0)
            .unwrap()
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        assert_eq!(matched.cross_sections().extinction_m2(), 0.0);
        assert_eq!(matched.eval_cos(0.0), Err(MieError::NoScattering));
    });
}

#[test]
fn g1_weak_absorption_matches_independent_high_precision_bessel_values() {
    // Independent reference: direct Riccati--Bessel power series for psi_n(mx),
    // psi_n(x), their analytic derivatives and xi_n(x), then Re(a)-|a|² at
    // 150 and 200 decimal digits. Thirty extra orders beyond the production
    // cutoff agreed to over 27 digits. These are full-series reference values;
    // the tolerance includes the existing finite-order truncation at large x.
    with_cx(|cx, _| {
        for (x, n, k, expected) in [
            (0.001, 1.5, 1e-24, 1.993_080_701_856_626_6e-27),
            (0.1, 1.5, 1e-24, 2.004_237_769_078_977e-25),
            (1.0, 1.5, 1e-24, 2.901_262_980_870_932_4e-24),
            (2.0, 1.5, 1e-24, 9.169_658_552_356_753e-24),
            (10.0, 1.5, 1e-24, 5.479_637_066_985_330_2e-23),
            (50.0, 1.5, 1e-24, 2.038_180_437_356_593_7e-22),
            (128.0, 1.33, 1e-24, 4.827_614_348_366_774e-22),
            (2.0, 1.5, 1.0, 1.591_598_898_565_839_3),
            (2.0, 1.1, 25.0, 0.014_612_713_789_718_393),
        ] {
            let sphere = sphere_for_x(x, n, k);
            let result = sphere
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap();
            let area = core::f64::consts::PI * sphere.radius_m.powi(2);
            close(
                result.cross_sections().absorption_m2 / (area * expected),
                1.0,
                1e-8,
            );
        }
    });
}

#[test]
fn g3_weak_particle_heat_scales_with_loss_without_a_roundoff_floor() {
    with_cx(|cx, _| {
        for x in [0.001, 0.1, 1.0, 2.0, 10.0, 128.0] {
            let reference = sphere_for_x(x, 1.5, 1e-18)
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap()
                .cross_sections();
            for k in [1e-24, 1e-30] {
                let actual = sphere_for_x(x, 1.5, k)
                    .at_wavelength(cx, 550.0, MieBudget::default())
                    .unwrap()
                    .cross_sections();
                assert!(actual.absorption_m2 > 0.0);
                close(
                    actual.absorption_m2 / reference.absorption_m2,
                    k / 1e-18,
                    2e-12 * k / 1e-18,
                );
                close(actual.scattering_m2 / reference.scattering_m2, 1.0, 2e-12);
            }
            let lossless = sphere_for_x(x, 1.5, 0.0)
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap();
            assert_eq!(lossless.cross_sections().absorption_m2, 0.0);
        }
    });
}

#[test]
fn g1_weak_mie_absorption_reaches_native_spectral_beam_heating() {
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::scattering::ScatteringGrid;
    use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget};
    use fs_geom::{Point3, Vec3};

    with_cx(|cx, _| {
        for k in [0.0, 1e-24] {
            let sphere = sphere_for_x(2.0, 1.5, k);
            let spectrum = sphere
                .visible_spectrum(
                    cx,
                    &[380.0, 550.0, 780.0],
                    MieBudget::default(),
                    0.002,
                    4096,
                )
                .unwrap();
            let densities = [1e-4];
            let absorption = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &densities, [0.0; 3], [1.0; 3]),
                spectrum.absorption,
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    max_steps_per_segment: 8,
                },
            )
            .unwrap();
            let medium =
                ScatteringGrid::try_new(cx, absorption, spectrum.scattering, spectrum.phase)
                    .unwrap();
            let report = medium
                .deposit_beam_power(
                    cx,
                    Point3::new(0.5, 0.5, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    1.0,
                    &[550.0],
                    &[1e12],
                    1,
                )
                .unwrap();
            // Independent Q values above, converted to particle mass and then
            // integrated through the uniform cell. Elastic scattering can be
            // 10^23 times the heat; it is still not deposited as heat.
            let mass_factor = 0.75 / (1000.0 * sphere.radius_m);
            let ka = if k == 0.0 {
                0.0
            } else {
                9.169_658_552_356_753e-24 * mass_factor
            };
            let ks = 1.798_418_163_226_263_9 * mass_factor;
            let removed = 1e12 * -(-densities[0] * (ka + ks)).exp_m1();
            let expected_heat = removed * ka / (ka + ks);
            close(
                report.absorbed_power_w[0],
                expected_heat,
                2e-11 * expected_heat,
            );
            close(report.scattered_power_w[0], removed * ks / (ka + ks), 0.02);
            assert_eq!(report.cells.len(), 1);
            assert_eq!(
                report.cells[0].absorbed_power_w[0],
                report.absorbed_power_w[0]
            );
            close(
                report.transmitted_power_w[0]
                    + report.absorbed_power_w[0]
                    + report.scattered_power_w[0],
                1e12,
                0.02,
            );
        }
    });
}

#[test]
fn g1_angular_normalization_moment_and_resolved_phase_table() {
    with_cx(|cx, _| {
        for (x, n, k) in [
            (0.01, 1.5, 0.0),
            (2.0, 1.5, 1.0),
            (12.0, 1.33, 0.0),
            (128.0, 1.33, 0.0),
        ] {
            let solution = sphere_for_x(x, n, k)
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap();
            // Independent Gauss--Legendre quadrature resolves the angular
            // polynomial's degree, including x=128's narrow forward peak.
            // No normalized phase table is used to establish Mie normalization.
            let (mass, moment) = angular_integrals(&solution);
            close(mass, 1.0, 2e-11);
            close(moment, solution.cross_sections().mean_cosine, 2e-11);
            let report = solution.phase_table(cx, 0.002, MAX_PHASE_KNOTS).unwrap();
            assert!(report.max_sampled_relative_error <= 0.002);
            assert!(report.mean_cosine_error <= 0.002);
            assert!(report.normalization_error <= 0.002);
            // A separate angular grid checks held-out points after admission.
            for i in 0..1023 {
                let mu = -1.0 + 2.0 * (i as f64 + 0.37) / 1023.0;
                let actual = solution.eval_cos(mu).unwrap();
                let approx = report.phase.eval_cos(mu).unwrap();
                let scale = actual.max(approx).max(1.0 / (4.0 * core::f64::consts::PI));
                assert!((actual - approx).abs() / scale < 0.0025);
            }
        }
    });
}

pub(super) fn angular_integrals(solution: &impl MieAngularSolution) -> (f64, f64) {
    let count = solution.multipoles() + 1;
    let (mut mass, mut moment) = (0.0, 0.0);
    for index in 0..count {
        let mut root = (core::f64::consts::PI * (index as f64 + 0.75) / (count as f64 + 0.5)).cos();
        let mut derivative = 0.0;
        for _ in 0..32 {
            let (mut previous, mut current) = (1.0, root);
            for order in 2..=count {
                let next = ((2 * order - 1) as f64 * root * current
                    - (order - 1) as f64 * previous)
                    / order as f64;
                (previous, current) = (current, next);
            }
            derivative = count as f64 * (root * current - previous) / (root * root - 1.0);
            let correction = current / derivative;
            root -= correction;
            if correction.abs() < 2e-16 {
                break;
            }
        }
        let quadrature = 2.0 / ((1.0 - root * root) * derivative * derivative);
        let weight = core::f64::consts::TAU * quadrature * solution.eval_cos(root).unwrap();
        mass += weight;
        moment += root * weight;
    }
    (mass, moment)
}

#[test]
fn g1_visible_preparation_retains_spectral_physics_and_exposes_interpolation_error() {
    with_cx(|cx, _| {
        let sphere = sphere_for_x(0.15, 1.5, 0.02);
        let coarse = sphere
            .visible_spectrum(cx, &[380.0, 780.0], MieBudget::default(), 0.002, 4096)
            .unwrap();
        let wavelengths: Vec<f64> = (0..=15)
            .map(|i| 380.0 + i as f64 * (400.0 / 15.0))
            .collect();
        let fine = sphere
            .visible_spectrum(cx, &wavelengths, MieBudget::default(), 0.002, 4096)
            .unwrap();
        assert_eq!(fine.angular_reports.len(), wavelengths.len());
        assert!(
            fine.max_relative_spectral_midpoint_error
                < coarse.max_relative_spectral_midpoint_error / 20.0
        );
        assert!(fine.max_relative_spectral_midpoint_error < 0.015);
        assert!(fine.max_sampled_phase_midpoint_error < 0.003);
        assert!(fine.scattering.at(380.0).unwrap() > fine.scattering.at(780.0).unwrap());
        // Every declared knot is the original Mie cross section, not an HG or
        // Rayleigh cross-section fit; the phase really carries wavelength.
        for wavelength in wavelengths {
            let original = sphere
                .at_wavelength(cx, wavelength, MieBudget::default())
                .unwrap();
            assert_eq!(
                fine.scattering.at(wavelength).unwrap(),
                original.cross_sections().mass_scattering_m2_kg
            );
            assert_eq!(
                fine.absorption.at(wavelength).unwrap(),
                original.cross_sections().mass_absorption_m2_kg
            );
        }
        assert!(fine.phase.eval_cos(0.0).is_err());
        assert_ne!(
            fine.phase.eval_at(380.0, 1.0).unwrap(),
            fine.phase.eval_at(780.0, 1.0).unwrap()
        );
        assert_eq!(
            sphere.visible_spectrum(cx, &[400.0, 780.0], MieBudget::default(), 0.002, 4096),
            Err(MieError::Parameters)
        );
        assert_eq!(
            sphere.visible_spectrum(
                cx,
                &[380.0, 550.0, 550.0, 780.0],
                MieBudget::default(),
                0.002,
                4096
            ),
            Err(MieError::Parameters)
        );
    });
}

#[test]
fn g1_source_spectra_retain_infrared_particle_absorption_and_population_weights() {
    use fs_blake3::ContentHash;

    with_cx(|cx, _| {
        let source = ContentHash([53; 32]);
        let sphere = MieSphere::try_new(2e-9, 1.5, 0.02, 1.0, 1800.0).unwrap();
        let wavelengths = [8000.0, 10_000.0, 12_000.0];
        let spectrum = sphere
            .source_spectrum(cx, &wavelengths, source, MieBudget::default(), 0.002, 4096)
            .unwrap();
        assert_eq!(
            spectrum.absorption.wavelength_bounds_nm(),
            (8000.0, 12_000.0)
        );
        assert_eq!(spectrum.absorption.source_identity(), Some(source));
        assert_eq!(spectrum.scattering.source_identity(), Some(source));
        let super::super::PhaseFunction::SpectralTable(phase) = &spectrum.phase else {
            panic!("source phase was replaced by an achromatic approximation");
        };
        assert_eq!(phase.source_identity(), Some(source));
        // Independent electrostatic absorption of a small sphere:
        // C_abs/m = 6*pi*Im((m²-1)/(m²+2))/(lambda*rho).
        // This checks absolute m²/kg and vacuum-nanometre conversion in IR.
        let m = C64::new(1.5, 0.02);
        let contrast = (m * m - C64::ONE) / (m * m + C64::from_re(2.0));
        for wavelength in wavelengths {
            let expected = 6.0 * core::f64::consts::PI * contrast.im / (wavelength * 1e-9 * 1800.0);
            close(
                spectrum.absorption.at(wavelength).unwrap() / expected,
                1.0,
                2e-6,
            );
        }
        for wavelength in [550.0, 7999.0, 12_001.0] {
            assert!(spectrum.absorption.at(wavelength).is_err());
            assert!(spectrum.scattering.at(wavelength).is_err());
            assert!(spectrum.phase.eval_at(wavelength, 0.0).is_err());
        }
        let large = MieSphere::try_new(0.7e-6, 1.7, 0.3, 1.0, 3100.0).unwrap();
        let population = MiePopulation::try_new(&[(sphere, 2.0), (large, 1.0)]).unwrap();
        let prepared = population
            .source_spectrum(
                cx,
                &wavelengths,
                source,
                MiePopulationBudget::default(),
                0.002,
                4096,
            )
            .unwrap();
        for wavelength in wavelengths {
            let a = sphere
                .at_wavelength(cx, wavelength, MieBudget::default())
                .unwrap();
            let b = large
                .at_wavelength(cx, wavelength, MieBudget::default())
                .unwrap();
            let mass = 2.0 * sphere.particle_mass_kg().unwrap() + large.particle_mass_kg().unwrap();
            let expected =
                (2.0 * a.cross_sections().absorption_m2 + b.cross_sections().absorption_m2) / mass;
            close(
                prepared.absorption.at(wavelength).unwrap() / expected,
                1.0,
                2e-14,
            );
        }
    });
}

#[test]
fn g1_dispersive_source_is_queried_at_each_infrared_knot_and_midpoint() {
    use fs_blake3::ContentHash;

    with_cx(|cx, _| {
        let wavelengths = [8000.0, 10_000.0, 12_000.0];
        let index_at = |wavelength: f64| {
            C64::new(
                1.4 + 0.1 * (wavelength - 8000.0) / 4000.0,
                0.01 + 0.03 * (wavelength - 8000.0) / 4000.0,
            )
        };
        let mut sampled = Vec::new();
        let spectrum = MieOpticalSpectrum::from_population_source(
            cx,
            &wavelengths,
            ContentHash([54; 32]),
            MiePopulationBudget::default(),
            0.002,
            4096,
            |wavelength| {
                sampled.push(wavelength);
                let index = index_at(wavelength);
                MiePopulation::try_new(&[(
                    MieSphere::try_new(3e-9, index.re, index.im, 1.0, 1000.0)?,
                    1.0,
                )])
            },
        )
        .unwrap();
        assert_eq!(sampled, vec![8000.0, 10_000.0, 12_000.0, 9000.0, 11_000.0]);
        for wavelength in wavelengths {
            let m = index_at(wavelength);
            let contrast = (m * m - C64::ONE) / (m * m + C64::from_re(2.0));
            let expected = 6.0 * core::f64::consts::PI * contrast.im / (wavelength * 1e-9 * 1000.0);
            close(
                spectrum.absorption.at(wavelength).unwrap() / expected,
                1.0,
                3e-6,
            );
        }
        // The changing source index reverses the constant-index absorption trend.
        assert!(
            spectrum.absorption.at(12_000.0).unwrap() > spectrum.absorption.at(8000.0).unwrap()
        );
        assert!(spectrum.max_relative_spectral_midpoint_error > 0.0);
    });
}

#[test]
fn g0_source_spectrum_preserves_visible_results_and_refuses_unsupported_work() {
    use fs_blake3::ContentHash;

    with_cx(|cx, gate| {
        let sphere = sphere_for_x(0.2, 1.5, 0.02);
        let wavelengths = [380.0, 550.0, 780.0];
        let legacy = sphere
            .visible_spectrum(cx, &wavelengths, MieBudget::default(), 0.002, 4096)
            .unwrap();
        let source = ContentHash([55; 32]);
        let explicit = sphere
            .source_spectrum(cx, &wavelengths, source, MieBudget::default(), 0.002, 4096)
            .unwrap();
        for wavelength in [380.0, 465.0, 550.0, 665.0, 780.0] {
            assert_eq!(
                legacy.absorption.at(wavelength),
                explicit.absorption.at(wavelength)
            );
            assert_eq!(
                legacy.scattering.at(wavelength),
                explicit.scattering.at(wavelength)
            );
            assert_eq!(
                legacy.phase.eval_at(wavelength, 0.3),
                explicit.phase.eval_at(wavelength, 0.3)
            );
        }
        assert_eq!(legacy.angular_reports, explicit.angular_reports);
        assert_eq!(
            legacy.max_relative_spectral_midpoint_error,
            explicit.max_relative_spectral_midpoint_error
        );
        for invalid in [
            &[][..],
            &[8000.0][..],
            &[8000.0, 8000.0][..],
            &[0.0, 12_000.0][..],
            &[8000.0, f64::NAN][..],
        ] {
            assert_eq!(
                sphere.source_spectrum(cx, invalid, source, MieBudget::default(), 0.002, 4096),
                Err(MieError::Parameters)
            );
        }
        assert_eq!(
            sphere.source_spectrum(
                cx,
                &wavelengths,
                ContentHash([0; 32]),
                MieBudget::default(),
                0.002,
                4096
            ),
            Err(MieError::Parameters)
        );
        let limited = MieBudget {
            max_multipoles: 1,
            ..MieBudget::default()
        };
        assert_eq!(
            sphere.source_spectrum(cx, &wavelengths, source, limited, 0.002, 4096),
            Err(MieError::Budget)
        );
        gate.request();
        assert_eq!(
            sphere.source_spectrum(cx, &wavelengths, source, MieBudget::default(), 0.002, 4096),
            Err(MieError::Execution)
        );
    });
}

#[test]
fn g0_dispersive_preparation_preserves_source_refusals_and_fixed_particle_geometry() {
    use fs_blake3::ContentHash;

    #[derive(Debug, PartialEq)]
    enum SourceError {
        OutsideMeasuredSupport,
        Mie(MieError),
    }
    impl From<MieError> for SourceError {
        fn from(error: MieError) -> Self {
            Self::Mie(error)
        }
    }
    with_cx(|cx, _| {
        let population = MiePopulation::try_new(&[(
            MieSphere::try_new(1e-7, 1.5, 0.1, 1.0, 1000.0).unwrap(),
            1.0,
        )])
        .unwrap();
        let unavailable = MieOpticalSpectrum::from_population_source(
            cx,
            &[8000.0, 12_000.0],
            ContentHash([56; 32]),
            MiePopulationBudget::default(),
            0.002,
            4096,
            |wavelength| {
                if wavelength == 10_000.0 {
                    Err(SourceError::OutsideMeasuredSupport)
                } else {
                    Ok(population.clone())
                }
            },
        );
        assert_eq!(unavailable, Err(SourceError::OutsideMeasuredSupport));
        for changed in 0..3 {
            let result = MieOpticalSpectrum::from_population_source(
                cx,
                &[8000.0, 12_000.0],
                ContentHash([56; 32]),
                MiePopulationBudget::default(),
                0.002,
                4096,
                |wavelength| {
                    let later = wavelength != 8000.0;
                    MiePopulation::try_new(&[(
                        MieSphere::try_new(
                            if later && changed == 0 { 2e-7 } else { 1e-7 },
                            1.5,
                            0.1,
                            if later && changed == 1 { 1.1 } else { 1.0 },
                            if later && changed == 2 {
                                2000.0
                            } else {
                                1000.0
                            },
                        )?,
                        1.0,
                    )])
                },
            );
            assert_eq!(result, Err(MieError::Parameters));
        }
    });
}

#[test]
fn g0_medium_scaling_and_explicit_refusals() {
    with_cx(|cx, gate| {
        let original = MieSphere::try_new(1e-7, 1.5, 0.1, 1.0, 1000.0)
            .unwrap()
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        let transformed = MieSphere::try_new(1e-7, 3.0, 0.2, 2.0, 2000.0)
            .unwrap()
            .at_wavelength(cx, 1100.0, MieBudget::default())
            .unwrap();
        assert_eq!(
            original.cross_sections().scattering_m2,
            transformed.cross_sections().scattering_m2
        );
        close(
            original.cross_sections().mass_scattering_m2_kg
                / transformed.cross_sections().mass_scattering_m2_kg,
            2.0,
            1e-14,
        );
        assert!(MieSphere::try_new(1e-7, 1.5, -0.1, 1.0, 1000.0).is_err());
        assert!(MieSphere::try_new(f64::NAN, 1.5, 0.0, 1.0, 1000.0).is_err());
        assert_eq!(
            sphere_for_x(129.0, 1.5, 0.0).at_wavelength(cx, 550.0, MieBudget::default()),
            Err(MieError::Regime)
        );
        assert_eq!(
            sphere_for_x(2.0, 1.5, 0.0).at_wavelength(
                cx,
                550.0,
                MieBudget {
                    max_multipoles: 1,
                    max_recurrence_order: 560
                }
            ),
            Err(MieError::Budget)
        );
        assert_eq!(
            original.phase_table(cx, 1e-4, 8),
            Err(MieError::AngularResolution)
        );
        assert_eq!(
            original.phase_table(cx, 0.0, 4096),
            Err(MieError::Parameters)
        );
        gate.request();
        assert_eq!(
            sphere_for_x(2.0, 1.5, 0.0).at_wavelength(cx, 550.0, MieBudget::default()),
            Err(MieError::Execution)
        );
        assert_eq!(
            original.phase_table(cx, 0.01, 4096),
            Err(MieError::Execution)
        );
    });
}
