//! Interpretation of sampling diagnostics and the calculation that produced them.
//! These checks do not assign confidence or physical-accuracy bounds.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::PHOTOMETRY_SEMANTICS_VERSION;
use fs_render::photometry::observer::OBSERVER_ID;
use fs_render::photometry::path::{
    IlluminanceSampling, PHOTOMETRIC_PATH_SEMANTICS_VERSION, PathCalculationContext,
    PathErrorBudget, PathIlluminanceJob, PhotometricQuantity, PointFixtureSampling,
    SamplingUncertainty,
};

fn check_context(context: PathCalculationContext, sampler: Sampler, seed: u64, depth: u32) {
    assert_eq!(context.sampler, sampler);
    assert_eq!(context.seed, seed);
    assert_eq!(context.max_depth, depth);
    assert_eq!(context.observer_id, OBSERVER_ID);
    assert_eq!(
        context.photometry_semantics_version,
        PHOTOMETRY_SEMANTICS_VERSION
    );
    assert_eq!(
        context.path_semantics_version,
        PHOTOMETRIC_PATH_SEMANTICS_VERSION
    );
}

fn check_unknown_physical_error(report: PathErrorBudget) {
    assert!(report.calibration_error_bound.is_none());
    assert!(report.spectral_error_bound.is_none());
    assert!(report.geometry_error_bound.is_none());
    assert!(report.physical_model_error_bound.is_none());
    assert!(report.total_error_bound.is_none());
    assert_eq!(report.numerical.kind, fs_evidence::NumericalKind::NoClaim);
}

fn check_sampling_report(
    reported: Option<SamplingUncertainty>,
    sampler: Sampler,
    samples: u32,
    legacy_error: Option<f64>,
) {
    match (sampler, samples) {
        (Sampler::Iid, 1) => {
            assert_eq!(
                reported,
                Some(SamplingUncertainty::InsufficientIidSamples { samples: 1 })
            );
            assert!(legacy_error.is_none());
        }
        (Sampler::Iid, _) => {
            let Some(SamplingUncertainty::EstimatedIidStandardError {
                standard_error,
                samples: count,
            }) = reported
            else {
                panic!("two or more IID paths must retain their estimated sampling error");
            };
            assert_eq!(count, samples);
            assert_eq!(standard_error.to_bits(), legacy_error.unwrap().to_bits());
        }
        (Sampler::OwenSobol | Sampler::OwenSobolFullPath, _) => {
            assert_eq!(
                reported,
                Some(SamplingUncertainty::CorrelatedSamples { samples })
            );
            assert!(legacy_error.is_none());
        }
    }
}

#[test]
fn g0_camera_reporting_distinguishes_missing_iid_and_correlated_samples() {
    with_cx(|cx, _| {
        for illuminated in [false, true] {
            let scene = receiver_scene();
            let sources = [isotropic(
                [0.0, 0.0, 3.0],
                90.0,
                SourceSpectrum::equal_energy(),
            )];
            let active_sources = if illuminated { &sources[..] } else { &[] };
            for (sampler, seed, depth) in [(Sampler::Iid, 91, 3), (Sampler::OwenSobol, 701, 5)] {
                let mut job = PathLuminanceJob::new(
                    &scene,
                    active_sources,
                    PathMedium::Vacuum,
                    TransportOptions::default(),
                    Settings {
                        sampler,
                        seed,
                        ..settings(depth, 4)
                    },
                    cx,
                )
                .unwrap();
                let context = job.calculation_context();
                check_context(context, sampler, seed, depth);
                assert!(job.sampling_uncertainty_cd_m2(0).is_none());
                assert!(job.sampling_uncertainty_cd_m2(1).is_none());
                assert!(job.error_budget_cd_m2(0).is_none());
                assert!(job.error_budget_cd_m2(1).is_none());
                let mut previous = 0;
                for samples in [1, 2, 4] {
                    assert_eq!(
                        job.advance(cx, samples - previous).unwrap(),
                        samples - previous
                    );
                    previous = samples;
                    let estimate = job.pixel_estimate(0).unwrap();
                    let report = job.error_budget_cd_m2(0).unwrap();
                    assert_eq!(report.quantity, PhotometricQuantity::Luminance);
                    assert_eq!(report.value.to_bits(), estimate.luminance_cd_m2.to_bits());
                    assert_eq!(Some(report.sampling), job.sampling_uncertainty_cd_m2(0));
                    assert_eq!(report.calculation, context);
                    assert_eq!(report.point_sampling, job.point_sampling());
                    assert_eq!(
                        report.diffuse_transport_remainder,
                        job.diffuse_transport_remainder()
                    );
                    assert!(report.transport_truncation_upper_bound().is_none());
                    check_unknown_physical_error(report);
                    assert_eq!(estimate.samples, samples as u32);
                    check_sampling_report(
                        job.sampling_uncertainty_cd_m2(0),
                        sampler,
                        samples as u32,
                        estimate.estimated_standard_error_cd_m2,
                    );
                    assert_eq!(job.calculation_context(), context);
                    assert!(job.sampling_uncertainty_cd_m2(1).is_none());
                    assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                    assert!(estimate.total_error_bound_cd_m2.is_none());
                    if !illuminated {
                        assert_eq!(estimate.luminance_cd_m2, 0.0);
                        if sampler == Sampler::Iid && samples >= 2 {
                            assert_eq!(estimate.estimated_standard_error_cd_m2, Some(0.0));
                        }
                    }
                }
                assert!(job.is_complete());
            }
        }
    });
}

#[test]
fn g0_probe_reporting_separates_zero_sampling_error_from_physical_accuracy() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                    .unwrap(),
            ];
        for (sampler, seed, depth) in [(Sampler::Iid, 41, 7), (Sampler::OwenSobol, 911, 2)] {
            let mut job = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 4,
                    max_depth: depth,
                    sampler,
                    seed,
                },
                cx,
            )
            .unwrap();
            let context = job.calculation_context();
            check_context(context, sampler, seed, depth);
            assert!(job.sampling_uncertainty_lux(0).is_none());
            assert!(job.sampling_uncertainty_lux(1).is_none());
            assert!(job.error_budget_lux(0).is_none());
            assert!(job.error_budget_lux(1).is_none());
            let mut previous = 0;
            for samples in [1, 2, 4] {
                assert_eq!(
                    job.advance(cx, samples - previous).unwrap(),
                    samples - previous
                );
                previous = samples;
                let estimate = job.probe_estimate(0).unwrap();
                let report = job.error_budget_lux(0).unwrap();
                assert_eq!(report.quantity, PhotometricQuantity::Illuminance);
                assert_eq!(report.value.to_bits(), estimate.illuminance_lux.to_bits());
                assert_eq!(Some(report.sampling), job.sampling_uncertainty_lux(0));
                assert_eq!(report.calculation, context);
                assert_eq!(report.point_sampling, job.point_sampling());
                assert!(report.transport_truncation_upper_bound().is_none());
                check_unknown_physical_error(report);
                assert_eq!(estimate.samples, samples as u32);
                check_sampling_report(
                    job.sampling_uncertainty_lux(0),
                    sampler,
                    samples as u32,
                    estimate.estimated_standard_error_lux,
                );
                assert_eq!(job.calculation_context(), context);
                assert!(job.sampling_uncertainty_lux(1).is_none());
                // Direct lux is deterministic and the continuous term is
                // identically zero. A zero IID sampling estimate is valid;
                // source calibration and arithmetic remain uncertified.
                close(
                    estimate.direct_lux,
                    10.0,
                    2.0e-13,
                    "direct reporting fixture",
                );
                assert_eq!(estimate.continuous_lux, 0.0);
                if sampler == Sampler::Iid && samples >= 2 {
                    assert_eq!(estimate.estimated_standard_error_lux, Some(0.0));
                }
                assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                assert!(estimate.total_error_bound_lux.is_none());
            }
            assert!(job.is_complete());
        }
    });
}

#[test]
fn g1_g5_planned_cavity_budget_reaches_results_without_promoting_accuracy() {
    use fs_render::photometry::area::{
        EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters,
    };
    use fs_render::spectral::SpectralReflectance;

    with_cx(|cx, _| {
        for rho in [0.0, 0.5] {
            let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
            let mut cube = super::area::diffuse_emitting_cube(0.0);
            cube.material = Material::SpectralDiffuse {
                reflectance: SpectralReflectance::constant(rho).unwrap(),
            };
            scene.primitives.push(cube);
            let surfaces = [
                PhotometricAreaEmitter::try_new(
                    0,
                    24.0 * PI * 8.0,
                    1.0,
                    SourceSpectrum::equal_energy(),
                    EmissionSidedness::Front,
                )
                .unwrap(),
            ];
            let emitters = PhotometricEmitters {
                points: &[],
                surfaces: &surfaces,
            };
            let probes = [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                    .unwrap(),
                IlluminanceProbe::try_new(Point3::new(0.2, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0))
                    .unwrap(),
            ];
            let mut sampling = IlluminanceSampling {
                samples_per_probe: 4,
                max_depth: 1,
                sampler: Sampler::Iid,
                seed: 173,
            };
            let policy = PointFixtureSampling::Importance {
                connections_per_vertex: 2,
            };
            let initial = PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                sampling,
                cx,
            )
            .unwrap();
            let plan = initial
                .diffuse_vacuum_remainder()
                .unwrap()
                .depth_for_illuminance_budget(0.05, 64)
                .unwrap();
            sampling.max_depth = plan.max_depth();
            let make = || {
                PathIlluminanceJob::new_with_emitters(
                    &scene,
                    emitters,
                    PathMedium::Vacuum,
                    &probes,
                    TransportOptions::default(),
                    sampling,
                    cx,
                )
                .unwrap()
                .with_point_sampling(policy)
                .unwrap()
            };
            let mut full = make();
            let mut split = make();
            full.advance(cx, usize::MAX).unwrap();
            split.advance(cx, 1).unwrap();
            assert!(split.error_budget_lux(0).is_some());
            assert!(split.error_budget_lux(1).is_none());
            split.advance(cx, usize::MAX).unwrap();
            for index in 0..probes.len() {
                let report = full.error_budget_lux(index).unwrap();
                assert_eq!(Some(report), split.error_budget_lux(index));
                assert_eq!(report.measurement_index, index);
                assert_eq!(report.point_sampling, policy);
                assert_eq!(report.calculation.max_depth, plan.max_depth());
                let actual = full.diffuse_transport_remainder().unwrap();
                assert_eq!(
                    report.transport_truncation_upper_bound(),
                    Some(actual.illuminance_bound_lux())
                );
                assert!(actual.illuminance_bound_lux() <= 0.05);
                close(
                    actual.source_luminance_ceiling_cd_m2(),
                    8.0,
                    3e-12,
                    "retained source ceiling",
                );
                check_unknown_physical_error(report);
                let mut csv = Vec::new();
                PathErrorBudget::write_csv_header(&mut csv).unwrap();
                report.write_csv_row(&mut csv).unwrap();
                assert!(
                    String::from_utf8(csv)
                        .unwrap()
                        .contains("diffuse-surface-volume-geometric-tail-v1")
                );
            }
            assert!(full.error_budget_lux(probes.len()).is_none());
            assert_eq!(initial.completed_path_samples(), 0);
            assert_eq!(initial.calculation_context().max_depth, 1);
            let mut image = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    sampler: Sampler::Iid,
                    ..settings(plan.max_depth(), 2)
                },
                cx,
            )
            .unwrap();
            image.advance(cx, usize::MAX).unwrap();
            let report = image.error_budget_cd_m2(0).unwrap();
            assert_eq!(report.quantity, PhotometricQuantity::Luminance);
            assert_eq!(
                report.transport_truncation_upper_bound(),
                Some(image.diffuse_vacuum_remainder().unwrap().luminance_bound_cd_m2())
            );
            check_unknown_physical_error(report);
        }
    });
}
