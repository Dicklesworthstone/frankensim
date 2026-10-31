use super::*;
use crate::volumes::transport::scattering::tests::with_cx;
use core::f64::consts::PI;
use fs_blake3::ContentHash;

fn budget() -> OpticalMixtureBudget {
    OpticalMixtureBudget { max_components: 8 }
}

fn component(
    weight: f64,
    scattering: MassScatteringSpectrum,
    phase: PhaseFunction,
) -> ScatteringMixtureComponent {
    ScatteringMixtureComponent {
        mass_weight: weight,
        scattering,
        phase,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_mixture_retains_analytic_rayleigh_and_cross_section_weighted_phase() {
    with_cx(|cx, _| {
        let components = [
            component(
                0.25,
                MassScatteringSpectrum::rayleigh(550.0, 2.0).unwrap(),
                PhaseFunction::Rayleigh,
            ),
            component(
                0.75,
                MassScatteringSpectrum::constant(3.0).unwrap(),
                PhaseFunction::henyey_greenstein(0.6).unwrap(),
            ),
        ];
        let phase = PhaseFunction::scattering_mixture(cx, &components, budget()).unwrap();
        let PhaseFunction::SpectralMixture(retained) = &phase else {
            unreachable!()
        };
        assert_eq!(retained.components(), &components);
        assert_eq!(retained.wavelength_bounds_nm(), (380.0, 780.0));
        for wavelength in [380.0_f64, 477.3, 580.0, 713.2, 780.0] {
            let gas = 0.5 * (550.0 / wavelength).powi(4);
            let particle = 2.25;
            for cosine in [-0.9, -0.1, 0.4, 0.97] {
                let rayleigh = 3.0 * (1.0 + cosine * cosine) / (16.0 * PI);
                let hg = (1.0 - 0.6_f64.powi(2))
                    / (4.0 * PI * (1.0 + 0.6_f64.powi(2) - 1.2 * cosine).powf(1.5));
                let expected = (gas * rayleigh + particle * hg) / (gas + particle);
                close(phase.eval_at(wavelength, cosine).unwrap(), expected, 6e-15);
            }
        }
    });
}

#[test]
fn g1_packet_mixture_preserves_energy_and_dispersive_angular_moments() {
    with_cx(|cx, _| {
        let particle_phase = PhaseFunction::spectral_table(&[
            (380.0, PhaseFunction::henyey_greenstein(0.65).unwrap()),
            (780.0, PhaseFunction::henyey_greenstein(-0.45).unwrap()),
        ])
        .unwrap();
        let phase = PhaseFunction::scattering_mixture(
            cx,
            &[
                component(
                    0.15,
                    MassScatteringSpectrum::rayleigh(550.0, 5.0).unwrap(),
                    PhaseFunction::Rayleigh,
                ),
                component(
                    0.85,
                    MassScatteringSpectrum::table(&[(380.0, 0.5), (780.0, 4.0)]).unwrap(),
                    particle_phase,
                ),
            ],
            budget(),
        )
        .unwrap();
        let wavelengths = [380.0_f64, 537.0, 780.0];
        let incoming = Vec3::new(0.0, 0.0, 1.0);
        let mut mass = [0.0; 3];
        let mut first = [0.0; 3];
        let mut second = [0.0; 3];
        let count = 24_576;
        for index in 0..count {
            let sample = phase
                .sample_packet(
                    &wavelengths,
                    incoming,
                    (f64::from(index) + 0.5) / f64::from(count),
                    0.37,
                )
                .unwrap();
            let cosine = sample.direction.z;
            let evaluated = phase.eval_packet(&wavelengths, cosine).unwrap();
            assert_eq!(sample.pdf, evaluated.pdf);
            assert!((sample.direction.dot(sample.direction) - 1.0).abs() < 1e-14);
            for lane in 0..3 {
                assert_eq!(sample.weight[lane], evaluated.value[lane] / evaluated.pdf);
                let contribution = sample.weight[lane] / f64::from(count);
                mass[lane] += contribution;
                first[lane] += contribution * cosine;
                second[lane] += contribution * cosine * cosine;
            }
        }
        for lane in 0..3 {
            let t = (wavelengths[lane] - 380.0) / 400.0;
            let gas = 0.75 * (550.0 / wavelengths[lane]).powi(4);
            let particle = 0.85 * (0.5 + 3.5 * t);
            let expected_first = particle * ((1.0 - t) * 0.65 - t * 0.45) / (gas + particle);
            let particle_second = ((1.0 - t) * (1.0 + 2.0 * 0.65_f64.powi(2))
                + t * (1.0 + 2.0 * 0.45_f64.powi(2)))
                / 3.0;
            let expected_second = (0.4 * gas + particle * particle_second) / (gas + particle);
            assert!((mass[lane] - 1.0).abs() < 5e-4, "{mass:?}");
            assert!((first[lane] - expected_first).abs() < 5e-4, "{first:?}");
            assert!((second[lane] - expected_second).abs() < 5e-4, "{second:?}");
        }
    });
}

#[test]
fn g0_zero_scattering_lanes_keep_a_normalized_matched_packet_proposal() {
    with_cx(|cx, _| {
        let directional = PhaseFunction::henyey_greenstein(0.7).unwrap();
        let phase = PhaseFunction::scattering_mixture(
            cx,
            &[
                component(
                    1.0,
                    MassScatteringSpectrum::table(&[(380.0, 0.0), (780.0, 2.0)]).unwrap(),
                    directional.clone(),
                ),
                component(
                    0.0,
                    MassScatteringSpectrum::constant(30.0).unwrap(),
                    PhaseFunction::Rayleigh,
                ),
            ],
            budget(),
        )
        .unwrap();
        for cosine in [-1.0, -0.2, 0.9, 1.0] {
            assert_eq!(
                phase.eval_at(380.0, cosine).unwrap(),
                PhaseFunction::Isotropic.eval_cos(cosine).unwrap()
            );
            assert_eq!(
                phase.eval_at(780.0, cosine).unwrap(),
                directional.eval_cos(cosine).unwrap()
            );
        }
        for u in [0.0, 0.11, 0.49, 0.51, 0.83, fs_math::next_down(1.0)] {
            let sample = phase
                .sample_packet(&[380.0, 780.0], Vec3::new(0.0, 0.0, 1.0), u, 0.2)
                .unwrap();
            let values = [
                PhaseFunction::Isotropic
                    .eval_cos(sample.direction.z)
                    .unwrap(),
                directional.eval_cos(sample.direction.z).unwrap(),
            ];
            let expected_pdf = 0.5 * values[0] + 0.5 * values[1];
            assert_eq!(sample.pdf, expected_pdf);
            assert_eq!(sample.weight, values.map(|value| value / expected_pdf));
        }
        let zero = PhaseFunction::scattering_mixture(
            cx,
            &[
                component(
                    0.0,
                    MassScatteringSpectrum::constant(1.0).unwrap(),
                    directional,
                ),
                component(
                    2.0,
                    MassScatteringSpectrum::constant(0.0).unwrap(),
                    PhaseFunction::Rayleigh,
                ),
            ],
            budget(),
        )
        .unwrap();
        assert_eq!(zero.eval_at(527.0, 0.4).unwrap(), 1.0 / (4.0 * PI));
        assert_eq!(
            zero.sample_packet(&[400.0, 700.0], Vec3::new(0.0, 0.0, 1.0), 0.2, 0.3)
                .unwrap()
                .weight,
            [1.0; 2]
        );
    });
}

#[test]
fn g5_one_active_spectral_component_retains_its_original_phase_and_sampler() {
    with_cx(|cx, _| {
        let original = PhaseFunction::spectral_table(&[
            (380.0, PhaseFunction::Rayleigh),
            (780.0, PhaseFunction::henyey_greenstein(0.5).unwrap()),
        ])
        .unwrap();
        let mixture = PhaseFunction::scattering_mixture(
            cx,
            &[
                component(
                    0.0,
                    MassScatteringSpectrum::constant(200.0).unwrap(),
                    PhaseFunction::Isotropic,
                ),
                component(
                    0.3,
                    MassScatteringSpectrum::rayleigh(600.0, 4.0).unwrap(),
                    original.clone(),
                ),
            ],
            budget(),
        )
        .unwrap();
        let wavelengths = [380.0, 507.0, 780.0];
        for cosine in [-1.0, 0.0, 0.4, 1.0] {
            assert_eq!(
                mixture.eval_packet(&wavelengths, cosine).unwrap(),
                original.eval_packet(&wavelengths, cosine).unwrap()
            );
        }
        for u in [0.0, 0.21, 0.51, 0.87, fs_math::next_down(1.0)] {
            let incoming = Vec3::new(0.6, 0.0, 0.8);
            assert_eq!(
                mixture
                    .sample_packet(&wavelengths, incoming, u, 0.17)
                    .unwrap(),
                original
                    .sample_packet(&wavelengths, incoming, u, 0.17)
                    .unwrap()
            );
        }
    });
}

#[test]
fn g0_mixture_refuses_invalid_loading_domains_nesting_and_work_limits() {
    with_cx(|cx, gate| {
        let base = component(
            1.0,
            MassScatteringSpectrum::constant(1.0).unwrap(),
            PhaseFunction::Isotropic,
        );
        assert_eq!(
            PhaseFunction::scattering_mixture(cx, &[], budget()),
            Err(VolumeTransportError::Mixture)
        );
        for invalid in [-0.1, f64::NAN, f64::INFINITY] {
            let mut bad = base.clone();
            bad.mass_weight = invalid;
            assert_eq!(
                PhaseFunction::scattering_mixture(cx, &[bad], budget()),
                Err(VolumeTransportError::Mixture)
            );
        }
        let broad = OpticalMixtureBudget {
            max_components: 128,
        };
        assert!(
            PhaseFunction::scattering_mixture(
                cx,
                &vec![base.clone(); MAX_OPTICAL_MIXTURE_COMPONENTS],
                broad
            )
            .is_ok()
        );
        assert_eq!(
            PhaseFunction::scattering_mixture(
                cx,
                &vec![base.clone(); MAX_OPTICAL_MIXTURE_COMPONENTS + 1],
                broad
            ),
            Err(VolumeTransportError::Budget)
        );
        assert_eq!(
            PhaseFunction::scattering_mixture(
                cx,
                &[base.clone()],
                OpticalMixtureBudget { max_components: 0 }
            ),
            Err(VolumeTransportError::Budget)
        );
        let foreign = component(
            1.0,
            MassScatteringSpectrum::source_table(
                &[(250.0, 1.0), (780.0, 1.0)],
                ContentHash([9; 32]),
            )
            .unwrap(),
            PhaseFunction::Rayleigh,
        );
        assert_eq!(
            PhaseFunction::scattering_mixture(cx, &[base.clone(), foreign], budget()),
            Err(VolumeTransportError::Wavelength)
        );
        let narrow_phase = PhaseFunction::source_spectral_table(
            &[
                (400.0, PhaseFunction::Rayleigh),
                (700.0, PhaseFunction::Rayleigh),
            ],
            ContentHash([3; 32]),
        )
        .unwrap();
        assert_eq!(
            PhaseFunction::scattering_mixture(
                cx,
                &[component(1.0, base.scattering, narrow_phase)],
                budget()
            ),
            Err(VolumeTransportError::Wavelength)
        );
        let mixture = PhaseFunction::scattering_mixture(cx, &[base.clone()], budget()).unwrap();
        assert_eq!(
            PhaseFunction::scattering_mixture(
                cx,
                &[component(1.0, base.scattering, mixture.clone())],
                budget()
            ),
            Err(VolumeTransportError::Phase)
        );
        assert_eq!(
            PhaseFunction::spectral_table(&[
                (380.0, mixture.clone()),
                (780.0, PhaseFunction::Rayleigh)
            ]),
            Err(VolumeTransportError::Phase)
        );
        assert_eq!(mixture.eval_cos(0.0), Err(VolumeTransportError::Wavelength));
        assert_eq!(
            mixture.eval_at(800.0, 0.0),
            Err(VolumeTransportError::Wavelength)
        );
        assert_eq!(
            mixture.eval_packet(&[], 0.0),
            Err(VolumeTransportError::Sample)
        );
        gate.request();
        assert_eq!(
            PhaseFunction::scattering_mixture(cx, &[base], budget()),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
