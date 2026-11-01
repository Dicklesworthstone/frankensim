//! Authored mathematical data, not measurements of a physical fixture.
use super::*;
use crate::instances::RigidTransform;
use crate::photometry::medium::HomogeneousMedium;
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::photometry::TransportOptions;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use std::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x4643_414c, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
fn id(byte: u8) -> ContentHash { ContentHash([byte; 32]) }
fn probe(z: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(0.0, 0.0, z), Vec3::new(0.0, 0.0, -1.0)).unwrap()
}
fn fixture(lumens: f64) -> PhotometricSource {
    PhotometricSource::try_new(RigidTransform::identity(), lumens, 1.0,
        AngularDistribution::isotropic(), SourceSpectrum::equal_energy()).unwrap()
}
fn measurement(z: f64, lo: f64, hi: f64, evidence: u8) -> IlluminanceMeasurement {
    IlluminanceMeasurement::try_new(probe(z), Interval::new(lo, hi), Interval::point(2.0), id(evidence)).unwrap()
}
fn observations() -> [IlluminanceMeasurement; 2] {
    // Nominal intensity 900 cd gives 100 lx at 3 m and 25 lx at 6 m.
    // True gain 4/5 plus background 2 lx gives 82 lx and 22 lx.
    [measurement(3.0, 81.0, 83.0, 1), measurement(6.0, 21.875, 22.125, 2)]
}

#[test]
fn inverse_constraints_keep_background_and_zero_support_distinct() {
    let fit = compatible_gain(Interval::point(100.0), Interval::new(81.0, 83.0), Interval::point(2.0)).unwrap();
    assert!(fit.contains(0.79) && fit.contains(0.81));
    assert!(fit.width() < 0.020000000001);
    assert!(compatible_gain(Interval::point(100.0), Interval::point(1.0), Interval::point(2.0)).is_none());
    let unbounded = Interval::new(0.0, f64::INFINITY);
    assert_eq!(compatible_gain(Interval::point(0.0), Interval::new(1.0, 3.0), Interval::point(2.0)), Some(unbounded));
    assert!(compatible_gain(Interval::point(0.0), Interval::point(3.0), Interval::point(2.0)).is_none());
    assert_eq!(compatible_gain(Interval::new(10.0, 20.0), Interval::point(2.0), Interval::point(2.0)), Some(Interval::point(0.0)));
    assert_eq!(compatible_gain(Interval::new(0.0, 20.0), Interval::point(2.0), Interval::point(2.0)), Some(unbounded));
    let ambiguous = compatible_gain(Interval::new(0.0, 100.0), Interval::new(81.0, 83.0), Interval::point(2.0)).unwrap();
    assert!(ambiguous.lo() <= 0.79 && ambiguous.lo() > 0.78);
    assert_eq!(ambiguous.hi(), f64::INFINITY);
}

#[test]
fn subnormal_and_overflow_constraints_never_manufacture_finite_certainty() {
    let tiny = f64::from_bits(1);
    let range = compatible_gain(Interval::point(f64::MAX), Interval::point(tiny), Interval::point(0.0)).unwrap();
    assert_eq!(range.lo(), 0.0);
    assert!(range.hi() > 0.0 && range.hi() < f64::MIN_POSITIVE);
    let range = compatible_gain(Interval::point(tiny), Interval::point(1.0), Interval::point(0.0)).unwrap();
    assert_eq!(range.hi(), f64::INFINITY);
    assert!(range.lo() > 1e300);
}

#[test]
fn two_distances_fit_one_multi_emitter_gain_and_feed_forward_flux_bounds() {
    with_cx(|cx, _| {
        let sources = [fixture(1800.0 * PI), fixture(1800.0 * PI)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let data = observations();
        let mut job = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
        assert!(job.result().is_none());
        assert_eq!(job.advance(cx, 1).unwrap(), 1);
        assert!(job.result().is_none());
        assert_eq!(job.advance(cx, usize::MAX).unwrap(), 1);
        let result = job.result().unwrap();
        assert_eq!(result.status(), FluxCalibrationStatus::NotExcluded);
        let gain = result.scale_enclosure().unwrap();
        assert!(gain.lo() <= 0.795 && gain.hi() >= 0.805);
        assert!(gain.lo() > 0.794999999 && gain.hi() < 0.805000001);
        let bounds = [result.flux_scale_bound(); 2];
        let report = scene.illuminance_with_flux_bounds(cx, &probe(4.0), &bounds).unwrap();
        // Independent inverse-square prediction: 900/16 * [0.795, 0.805].
        let predicted = report.flux_conditioned_lux.unwrap();
        assert!(predicted.contains(44.71875) && predicted.contains(45.28125));
        assert_eq!(report.nominal.estimate.lux.to_bits(), scene.illuminance(cx, &probe(4.0)).unwrap().lux.to_bits());
        assert!(report.nominal.estimate.total_error_bound_lux.is_none());
        assert_eq!(sources[0].flux_lm(), 1800.0 * PI);
        assert_eq!(result.processed_measurements(), 2);
    });
}

#[test]
fn additional_data_can_disprove_without_refitting_or_claiming_pass_on_overlap() {
    with_cx(|cx, _| {
        let sources = [fixture(3600.0 * PI)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let data = observations();
        let mut job = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
        let good = measurement(4.0, 46.9, 47.1, 4);
        assert!(matches!(job.check_measurement(cx, &good), Err(FluxCalibrationError::NoFiniteFit)));
        job.advance(cx, usize::MAX).unwrap();
        let before = job.result();
        let check = job.check_measurement(cx, &good).unwrap();
        assert_eq!(check.compatibility, MeasurementCompatibility::NotExcluded);
        assert!(check.predicted_measured_lux.contains(46.71875));
        assert!(check.predicted_measured_lux.contains(47.28125));
        let bad = measurement(4.0, 52.0, 53.0, 5);
        assert_eq!(job.check_measurement(cx, &bad).unwrap().compatibility, MeasurementCompatibility::Contradicted);
        assert_eq!(job.result(), before);
        assert_eq!(job.processed_measurements(), 2);
        assert!(check.nominal.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn contradictions_stop_before_later_unusable_measurements_and_export_no_gain() {
    with_cx(|cx, _| {
        let sources = [fixture(3600.0 * PI)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let data = [observations()[0], measurement(6.0, 27.0, 28.0, 6), measurement(0.0, 2.0, 2.0, 7)];
        let mut job = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
        assert_eq!(job.advance(cx, usize::MAX).unwrap(), 2);
        let result = job.result().unwrap();
        assert_eq!(result.status(), FluxCalibrationStatus::Inconsistent { measurement: 1 });
        assert!(result.scale_enclosure().is_none() && result.flux_scale_bound().is_none());
        assert_eq!(job.advance(cx, usize::MAX).unwrap(), 0);
    });
}

#[test]
fn dark_and_natively_occluded_experiments_do_not_identify_a_gain() {
    use crate::charts::TriMesh;
    use crate::tracer::{Material, Primitive, Shape};
    use crate::spectral::SpectralReflectance;
    with_cx(|cx, _| {
        let blocker = [Primitive {
            shape: Shape::Mesh(TriMesh::new(vec![[-1.0,-1.0,1.0],[1.0,-1.0,1.0],[0.0,1.0,1.0]], vec![[0,1,2]])),
            material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() }, emission: None,
        }];
        for blocked in [false, true] {
            let sources = [fixture(if blocked { 3600.0 * PI } else { 0.0 })];
            let geometry = if blocked { &blocker[..] } else { &[] };
            let scene = PhotometricScene::try_new(geometry, &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
            let data = [measurement(3.0, 2.0, 2.0, 1)];
            let mut job = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
            job.advance(cx, 1).unwrap();
            assert_eq!(job.result().unwrap().status(), FluxCalibrationStatus::Unbounded);
            assert!(job.result().unwrap().flux_scale_bound().is_none());
            let positive = [measurement(3.0, 3.0, 3.0, 2)];
            let mut job = FluxCalibrationJob::new(&scene, &positive, id(3), cx).unwrap();
            job.advance(cx, 1).unwrap();
            assert_eq!(job.result().unwrap().status(), FluxCalibrationStatus::Inconsistent { measurement: 0 });
        }
    });
}

#[test]
fn cancellation_replay_and_evidence_changes_preserve_constraint_meaning() {
    with_cx(|cx, gate| {
        let sources = [fixture(3600.0 * PI)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let data = observations();
        let mut full = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
        full.advance(cx, usize::MAX).unwrap();
        let mut split = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
        split.advance(cx, 1).unwrap();
        gate.request();
        assert!(matches!(split.advance(cx, 1), Err(FluxCalibrationError::Transport(PhotometryError::Cancelled))));
        assert_eq!(split.processed_measurements(), 1);
        assert!(split.result().is_none());
        with_cx(|resume, _| {
            split.advance(resume, 1).unwrap();
            assert_eq!(split.result(), full.result());
            let mut changed = data;
            changed[0].evidence = id(8);
            let mut other = FluxCalibrationJob::new(&scene, &changed, id(3), resume).unwrap();
            other.advance(resume, usize::MAX).unwrap();
            assert_eq!(other.result().unwrap().scale_enclosure(), full.result().unwrap().scale_enclosure());
            assert_ne!(other.result().unwrap().evidence_identity(), full.result().unwrap().evidence_identity());
        });
    });
}

#[test]
fn colored_native_grid_bounds_are_inverted_and_budget_refusal_is_not_skipped() {
    use crate::photometry::transport::PhotometricMedium;
    use crate::photometry::source::{SPECTRAL_BINS, wavelength_nm, photopic_weight_nm};
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
    with_cx(|cx, _| {
        let density = [1.0, 2.0];
        let grid = AbsorbingGrid::try_new(cx, VolumeGrid::new([1,1,2], &density, [-0.5,-0.5,1.0], [1.0;3]),
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,1.0)]).unwrap(), None,
            AbsorbingGridBudget { max_grid_cells: 2, max_steps_per_segment: 128 }).unwrap();
        let bin = 20;
        let mut bins = [0.0; SPECTRAL_BINS]; bins[bin] = 1.0;
        let sources = [PhotometricSource::try_new(RigidTransform::identity(), 3600.0 * PI, 1.0,
            AngularDistribution::isotropic(), SourceSpectrum::try_from_bins(bins).unwrap()).unwrap()];
        let scene = PhotometricScene::try_with_medium(&[], &sources, PhotometricMedium::Absorbing(&grid), TransportOptions::default(), cx).unwrap();
        // Independent two-cell column mass 3 kg/m2, retained monochromatic
        // source/observer and affine kappa(lambda). No renderer result is fitted.
        let transfer = 683.0 * sources[0].spectrum().bins_w_per_lumen_nm()[bin] * photopic_weight_nm(bin)
            * (-3.0 * (wavelength_nm(bin)-380.0)/400.0).exp();
        let readings = [100.0, 25.0].map(|nominal| 0.8 * nominal * transfer + 2.0);
        let data = [measurement(3.0, readings[0]-0.001, readings[0]+0.001, 1), measurement(6.0, readings[1]-0.001, readings[1]+0.001, 2)];
        let mut job = FluxCalibrationJob::new(&scene, &data, id(3), cx).unwrap();
        job.set_grid_enclosure_budget(GridEnclosureBudget { max_nodes_per_segment: 1 }).unwrap();
        assert!(matches!(job.advance(cx, 1), Err(FluxCalibrationError::NumericalUnavailable {
            measurement: 0, reason: EnclosureUnavailable::MediumBudgetExceeded,
        })));
        assert_eq!(job.processed_measurements(), 0);
        assert!(job.result().is_none());
        job.set_grid_enclosure_budget(GridEnclosureBudget::default()).unwrap();
        job.advance(cx, usize::MAX).unwrap();
        assert!(job.result().unwrap().scale_enclosure().unwrap().contains(0.8));
        assert!(job.result().unwrap().scale_enclosure().unwrap().width() < 0.001);
        assert!(job.set_grid_enclosure_budget(GridEnclosureBudget::default()).is_err());
    });
}

#[test]
fn invalid_declarations_and_empty_or_overbudget_jobs_refuse() {
    assert!(matches!(IlluminanceMeasurement::try_new(probe(3.0), Interval::new(-1.0, 1.0), Interval::point(0.0), id(1)), Err(FluxCalibrationError::InvalidMeasurement)));
    assert!(matches!(IlluminanceMeasurement::try_new(probe(3.0), Interval::point(1.0), Interval::point(f64::INFINITY), id(1)), Err(FluxCalibrationError::InvalidMeasurement)));
    assert!(matches!(IlluminanceMeasurement::try_new(probe(3.0), Interval::point(1.0), Interval::point(0.0), id(0)), Err(FluxCalibrationError::MissingEvidence)));
    with_cx(|cx, _| {
        let sources = [fixture(3600.0 * PI)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions { max_outputs: 1, ..TransportOptions::default() }, cx).unwrap();
        assert!(matches!(FluxCalibrationJob::new(&scene, &[], id(3), cx), Err(FluxCalibrationError::Transport(PhotometryError::InvalidOptions))));
        let data = observations();
        assert!(matches!(FluxCalibrationJob::new(&scene, &data, id(3), cx), Err(FluxCalibrationError::Transport(PhotometryError::OutputBudgetExceeded))));
        assert!(matches!(FluxCalibrationJob::new(&scene, &data[..1], id(0), cx), Err(FluxCalibrationError::MissingEvidence)));
    });
}
