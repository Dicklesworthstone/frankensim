use super::super::tests::{close, sphere_for_x, with_cx};
use super::super::{MieBudget, MieOpticalSpectrum, MiePopulation, MiePopulationBudget, MieSphere};
use super::*;
use crate::polarization::Stokes;
use fs_math::c64::C64;

#[test]
fn mie_full_stokes_reduces_to_the_rayleigh_matrix() {
    with_cx(|cx, _| {
        for k in [0.0, 0.1] {
            let solution = sphere_for_x(2e-4, 1.5, k)
                .at_wavelength(cx, 550.0, MieBudget::default())
                .unwrap();
            for c in [-1.0, -0.3, 0.0, 0.4, 1.0] {
                let value = solution.polarized_value(c).unwrap();
                let m = value.mueller().0;
                let scale = 3.0 / (16.0 * core::f64::consts::PI);
                close(m[0][0], scale * (1.0 + c * c), 2e-8);
                close(m[0][1], scale * (1.0 - c * c), 2e-8);
                close(m[2][2], 2.0 * scale * c, 2e-8);
                close(m[3][2], 0.0, 2e-8);
                close(
                    value.unpolarized_per_sr(),
                    solution.eval_cos(c).unwrap(),
                    2e-15,
                );
            }
        }
    });
}

#[test]
fn actual_mie_phase_preserves_complex_jones_transport() {
    with_cx(|cx, _| {
        let solution = sphere_for_x(2.0, 1.5, 0.2)
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        for c in [-0.8, -0.1, 0.3, 0.95] {
            let [s, p] = solution.amplitudes(c).unwrap();
            let value = solution.polarized_value(c).unwrap();
            let input = [C64::new(0.3, 0.4), C64::new(-0.2, 0.7)];
            let cross = input[0] * input[1].conj();
            let state = Stokes {
                i: input[0].norm_sq() + input[1].norm_sq(),
                q: input[0].norm_sq() - input[1].norm_sq(),
                u: 2.0 * cross.re,
                v: 2.0 * cross.im,
            };
            let output = [s * input[0], p * input[1]];
            let z = output[0] * output[1].conj();
            let scale = 1.0 / (core::f64::consts::TAU * solution.scattering_sum);
            let expected = [
                (output[0].norm_sq() + output[1].norm_sq()) * scale,
                (output[0].norm_sq() - output[1].norm_sq()) * scale,
                2.0 * z.re * scale,
                2.0 * z.im * scale,
            ];
            let actual = value.apply(state).unwrap();
            for (a, b) in [actual.i, actual.q, actual.u, actual.v]
                .into_iter()
                .zip(expected)
            {
                close(a, b, 2e-14);
            }
        }
        assert!(
            solution
                .polarized_value(0.3)
                .unwrap()
                .cross_per_sr()
                .im
                .abs()
                > 1e-5
        );
    });
}

#[test]
fn independent_particle_population_uses_scattered_power_coherency() {
    with_cx(|cx, _| {
        let a = sphere_for_x(0.7, 1.33, 0.0);
        let b = sphere_for_x(2.0, 1.7, 0.2);
        let population = MiePopulation::try_new(&[(a, 1.0), (b, 2.0)]).unwrap();
        let mixed = population
            .at_wavelength(cx, 550.0, MiePopulationBudget::default())
            .unwrap();
        let x = a.at_wavelength(cx, 550.0, MieBudget::default()).unwrap();
        let y = b.at_wavelength(cx, 550.0, MieBudget::default()).unwrap();
        let total = x.cross_sections().scattering_m2 + 2.0 * y.cross_sections().scattering_m2;
        let weight = x.cross_sections().scattering_m2 / total;
        for c in [-1.0, -0.4, 0.0, 0.7, 1.0] {
            let expected = x
                .polarized_value(c)
                .unwrap()
                .interpolate(y.polarized_value(c).unwrap(), 1.0 - weight)
                .unwrap();
            let actual = mixed.polarized_value(c).unwrap();
            assert!(discrepancy(actual, expected) < 2e-14);
            close(
                actual.unpolarized_per_sr(),
                mixed.eval_cos(c).unwrap(),
                2e-15,
            );
        }
    });
}

#[test]
fn full_coherency_refinement_resolves_retardance_and_sampling() {
    with_cx(|cx, _| {
        let solution = sphere_for_x(2.0, 1.5, 0.2)
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        let report = solution.polarized_phase_table(cx, 0.003, 4096).unwrap();
        assert!(report.max_sampled_relative_error <= 0.003);
        let table = report.polarization.as_ref().unwrap();
        let mut worst = 0.0_f64;
        for j in 0..=2048 {
            let c = -1.0 + 2.0 * j as f64 / 2048.0;
            worst = worst.max(discrepancy(
                solution.polarized_value(c).unwrap(),
                table.evaluate(c).unwrap(),
            ));
        }
        assert!(worst < 0.0033, "{worst}");
        assert!(matches!(
            solution.polarized_phase_table(cx, 1e-5, 17),
            Err(MieError::AngularResolution)
        ));
        let matched = sphere_for_x(1.0, 1.0, 0.0)
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        assert_eq!(matched.amplitudes(0.4).unwrap(), [C64::ZERO; 2]);
        assert!(matches!(
            matched.polarized_phase_table(cx, 0.003, 4096),
            Err(MieError::NoScattering)
        ));
    });
}

#[test]
fn original_dispersive_population_prepares_polarized_spectral_transport_once() {
    with_cx(|cx, _| {
        let calls = core::cell::Cell::new(0);
        let spectrum = MieOpticalSpectrum::from_population_source_polarized::<MieError>(
            cx,
            &[380.0, 550.0, 780.0],
            fs_blake3::ContentHash([23; 32]),
            MiePopulationBudget::default(),
            0.005,
            4096,
            |w| {
                calls.set(calls.get() + 1);
                MiePopulation::try_new(&[(
                    MieSphere::try_new(1e-7, 1.4 + 50.0 / w, 0.02, 1.0, 1200.0)?,
                    1.0,
                )])
            },
        )
        .unwrap();
        assert_eq!(calls.get(), 5);
        assert!(
            spectrum
                .angular_reports
                .iter()
                .all(|r| r.polarization.is_some())
        );
        assert!(spectrum.max_sampled_phase_midpoint_error.is_finite());
        assert!(spectrum.phase.spherical_value_at(379.0, 0.3).is_err());
        for w in [380.0, 465.0, 550.0, 665.0, 780.0] {
            let value = spectrum.phase.spherical_value_at(w, 0.3).unwrap();
            close(
                value.unpolarized_per_sr(),
                spectrum.phase.eval_at(w, 0.3).unwrap(),
                2e-15,
            );
        }
    });
}
