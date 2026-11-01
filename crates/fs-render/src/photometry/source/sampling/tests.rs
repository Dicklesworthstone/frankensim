use super::*;
use crate::photometry::source::SourceSpectrum;
use core::f64::consts::{FRAC_PI_2, PI};

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 71,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance}"
    );
}

#[test]
fn g1_analytic_direction_distributions_recover_independent_cdfs_and_pdfs() {
    with_cx(|cx, _| {
        let laws = [
            (AngularDistribution::isotropic(), -1.0, 0.0),
            (
                AngularDistribution::uniform_cone(0.7).unwrap(),
                0.7_f64.cos(),
                0.0,
            ),
            (
                AngularDistribution::cosine_power(0.0, FRAC_PI_2).unwrap(),
                0.0,
                0.0,
            ),
            (
                AngularDistribution::cosine_power(3.5, 0.8).unwrap(),
                0.8_f64.cos(),
                3.5,
            ),
            (
                AngularDistribution::cosine_power(1_000_000.0, 1e-4).unwrap(),
                1e-4_f64.cos(),
                1_000_000.0,
            ),
        ];
        for (law, cutoff, exponent) in laws {
            let prepared = law.prepare_direction_sampler(cx).unwrap();
            for u in [0.001, 0.1, 0.25, 0.8, 0.999] {
                let sample = prepared.sample([u, 0.37]).unwrap();
                let mu = sample.direction.z;
                let recovered = if exponent == 0.0 {
                    (1.0 - mu) / (1.0 - cutoff)
                } else {
                    (1.0 - mu.powf(exponent + 1.0)) / (1.0 - cutoff.powf(exponent + 1.0))
                };
                close(recovered, u, if exponent > 1000.0 { 4e-8 } else { 4e-14 });
                close(sample.direction.dot(sample.direction), 1.0, 5e-15);
                let native_density = law.intensity_per_steradian(sample.direction).unwrap();
                close(
                    sample.pdf_per_sr * sample.angular_weight(),
                    native_density,
                    native_density * 3e-9 + 1e-15,
                );
                close(sample.angular_weight(), 1.0, 3e-15);
            }
        }
    });
}

fn asymmetric_table() -> AngularDistribution {
    AngularDistribution::tabulated(
        vec![0.0, FRAC_PI_2, PI],
        vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
        vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 2.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ],
    )
    .unwrap()
}

#[test]
fn g1_bilinear_table_cdf_pdf_and_independent_angular_moments_agree() {
    with_cx(|cx, _| {
        let distribution = asymmetric_table();
        let prepared = distribution.prepare_direction_sampler(cx).unwrap();
        let mut mu_mean = 0.0;
        let mut cos_phi_mean = 0.0;
        let n = 256;
        for i in 0..n {
            let u = (i as f64 + 0.5) / n as f64;
            for j in 0..n {
                let v = (j as f64 + 0.5) / n as f64;
                let sample = prepared.sample([u, v]).unwrap();
                let mu = sample.direction.z;
                let native = distribution
                    .intensity_per_steradian(sample.direction)
                    .unwrap();
                close(sample.pdf_per_sr * sample.angular_weight(), native, 2e-15);
                // Integrating the unnormalized polar marginal gives 5*pi.
                let cdf = if mu >= 0.0 {
                    let t = 1.0 - mu;
                    (2.0 * t + t * t) / 5.0
                } else {
                    let t = -mu;
                    (3.0 + 4.0 * t - 2.0 * t * t) / 5.0
                };
                close(cdf, u, 2e-15);
                mu_mean += mu;
                cos_phi_mean += sample.direction.x
                    / (sample.direction.x * sample.direction.x
                        + sample.direction.y * sample.direction.y)
                        .sqrt();
            }
        }
        close(mu_mean / (n * n) as f64, 2.0 / 15.0, 4e-5);
        close(cos_phi_mean / (n * n) as f64, -16.0 / (5.0 * PI * PI), 8e-5);
    });
}

#[test]
fn g0_zero_probability_cells_are_skipped_without_truncating_nonzero_support() {
    with_cx(|cx, _| {
        let distribution = AngularDistribution::tabulated(
            vec![0.0, FRAC_PI_2, PI],
            vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
            vec![
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        let prepared = distribution.prepare_direction_sampler(cx).unwrap();
        for i in 0..1000 {
            let sample = prepared.sample([(i as f64 + 0.5) / 1000.0, 0.267]).unwrap();
            assert!(sample.direction.x < 0.0);
            assert!(sample.pdf_per_sr > 0.0);
        }
        let beam = AngularDistribution::beam_field(0.2, 0.4, 0.8).unwrap();
        let prepared = beam.prepare_direction_sampler(cx).unwrap();
        for i in 0..1000 {
            let sample = prepared.sample([(i as f64 + 0.5) / 1000.0, 0.3]).unwrap();
            assert!(sample.direction.z >= 0.4_f64.cos());
        }
    });
}

#[test]
fn g3_world_rotation_preserves_probability_and_calibrated_packet_power() {
    with_cx(|cx, _| {
        let distribution = asymmetric_table();
        let transform = RigidTransform::try_new(
            [0.0, (0.41_f64).sin(), 0.0, (0.41_f64).cos()],
            [30.0, -10.0, 17.0],
        )
        .unwrap();
        let spectrum = SourceSpectrum::equal_energy();
        let source = PhotometricSource::try_new(
            transform,
            900.0,
            0.4,
            distribution.clone(),
            spectrum.clone(),
        )
        .unwrap();
        let local = distribution.prepare_direction_sampler(cx).unwrap();
        let world = source.prepare_direction_sampler(cx).unwrap();
        for uv in [[0.2, 0.8], [0.9, 0.5], [0.4, 0.06]] {
            let a = local.sample(uv).unwrap();
            let b = world.sample(uv).unwrap();
            let rotated = transform.transform_vector(a.direction);
            close(b.direction.x, rotated.x, 1e-15);
            close(b.direction.y, rotated.y, 1e-15);
            close(b.direction.z, rotated.z, 1e-15);
            assert_eq!(a.pdf_per_sr, b.pdf_per_sr);
            for wavelength in [382.5, 552.5, 777.5] {
                let packet = source.flux_lm()
                    * spectrum.value_w_per_lumen_nm(wavelength)
                    * b.angular_weight();
                let direct = source
                    .spectral_intensity_w_sr_nm(b.direction, wavelength)
                    .unwrap()
                    / b.pdf_per_sr;
                close(packet, direct, 2e-17);
            }
        }
    });
}

#[test]
fn g4_sampling_refuses_invalid_samples_lost_support_and_cancellation() {
    with_cx(|cx, gate| {
        let distribution = AngularDistribution::isotropic();
        let prepared = distribution.prepare_direction_sampler(cx).unwrap();
        for sample in [
            [0.0, 0.5],
            [0.5, 1.0],
            [f64::NAN, 0.5],
            [-0.1, 0.3],
            [0.2, f64::INFINITY],
        ] {
            assert_eq!(
                prepared.sample(sample),
                Err(SourceSamplingError::InvalidSample)
            );
        }
        let lost_support = AngularDistribution::tabulated(
            vec![0.0, FRAC_PI_2, PI],
            vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
            vec![
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1e-300, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            ],
        )
        .unwrap();
        assert!(matches!(
            lost_support.prepare_direction_sampler(cx),
            Err(SourceSamplingError::NumericalRange)
        ));
        gate.request();
        assert!(matches!(
            distribution.prepare_direction_sampler(cx),
            Err(SourceSamplingError::Cancelled(_))
        ));
    });
}
