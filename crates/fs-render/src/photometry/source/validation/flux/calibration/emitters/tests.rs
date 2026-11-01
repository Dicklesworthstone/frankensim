//! Authored deterministic observations, not measured fixture validation.
use super::*;
use crate::instances::RigidTransform;
use crate::photometry::{IlluminanceProbe, TransportOptions};
use crate::photometry::medium::HomogeneousMedium;
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::photometry::source::validation::flux::FluxEnclosureJob;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x454d_4341, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
fn evidence(value: u8) -> ContentHash { ContentHash([value; 32]) }
fn fixture(height: f64, intensity: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, height]).unwrap(),
        4.0 * std::f64::consts::PI * intensity, 1.0,
        AngularDistribution::isotropic(), SourceSpectrum::equal_energy()).unwrap()
}
fn probe() -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
}
fn measurement(lo: f64, hi: f64, background: f64, id: u8) -> IlluminanceMeasurement {
    IlluminanceMeasurement::try_new(probe(), Interval::new(lo, hi),
        Interval::point(background), evidence(id)).unwrap()
}
fn mixed_rows() -> [EmitterCalibrationObservation<'static>; 2] {
    [EmitterCalibrationObservation { measurement: measurement(111.9, 112.1, 2.0, 1), relative_output: &[1.0, 0.5] },
     EmitterCalibrationObservation { measurement: measurement(101.9, 102.1, 2.0, 2), relative_output: &[0.5, 1.0] }]
}
fn bound_contains(range: Interval, value: f64) {
    assert!(range.contains(value), "{range:?} excludes {value}");
}

#[test]
fn overlapping_emitters_are_separated_and_export_to_original_forward_flux() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
            TransportOptions::default(), cx).unwrap();
        let rows = mixed_rows();
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10),
            EmitterCalibrationBudget { max_sweeps: 64, ..Default::default() }, cx).unwrap();
        assert!(job.result().is_none());
        assert_eq!(job.advance(cx, 2).unwrap(), 2); // only prepare native rows
        assert!(job.result().is_none());
        assert_eq!(job.coefficient_count(), 4);
        job.advance(cx, usize::MAX).unwrap();
        let fit = job.result().unwrap();
        let scales = fit.scale_enclosures().unwrap();
        // Independent inverse: u=100*g0+25*g1, v=50*g0+50*g1.
        // g0=(u-v/2)/75 -> [.798,.802]; g1=(2*v-u)/75 -> [1.196,1.204].
        bound_contains(scales[0], 0.8); bound_contains(scales[1], 1.2);
        assert!(scales[0].lo() > 0.797_999_99 && scales[0].hi() < 0.802_000_01);
        assert!(scales[1].lo() > 1.195_999_99 && scales[1].hi() < 1.204_000_01);
        let identity = fit.evidence_identity;
        let holdout = EmitterCalibrationObservation {
            measurement: measurement(141.9, 142.1, 2.0, 3), relative_output: &[1.0, 1.0] };
        let check = job.check_observation(cx, holdout).unwrap();
        assert_eq!(check.compatibility, MeasurementCompatibility::NotExcluded);
        bound_contains(check.predicted_measured_lux, 142.0);
        assert_eq!(check.reference, scene.illuminance_with_enclosure(cx, &probe()).unwrap());
        assert!(check.reference.estimate.total_error_bound_lux.is_none());
        let bad = EmitterCalibrationObservation {
            measurement: measurement(160.0, 161.0, 2.0, 4), ..holdout };
        assert_eq!(job.check_observation(cx, bad).unwrap().compatibility, MeasurementCompatibility::Contradicted);
        assert_eq!(job.result().unwrap().evidence_identity, identity);
        let bounds = [fit.flux_scale_bound(0), fit.flux_scale_bound(1)];
        assert!(bounds.iter().all(Option::is_some));
        assert_ne!(bounds[0].unwrap().evidence_identity(), bounds[1].unwrap().evidence_identity());
        let probes = [probe()];
        let mut forward = FluxEnclosureJob::new(&scene, &probes, &bounds, cx).unwrap();
        forward.advance(cx, 1).unwrap();
        bound_contains(forward.results()[0].flux_conditioned_lux.unwrap(), 140.0);
        assert!(forward.results()[0].nominal.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn indistinguishable_emitters_are_not_given_a_unique_or_common_gain() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(3.0, 900.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = [EmitterCalibrationObservation {
            measurement: measurement(160.0, 160.0, 0.0, 1), relative_output: &[1.0, 1.0] }];
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        job.advance(cx, usize::MAX).unwrap();
        let fit = job.result().unwrap();
        assert_eq!(fit.termination, EmitterCalibrationTermination::Stable);
        for &range in fit.scale_enclosures().unwrap() {
            assert_eq!(range.lo(), 0.0);
            bound_contains(range, 1.6);
            assert!(range.hi() < 1.600_000_01);
        }
    });
}

#[test]
fn unobserved_emitters_stay_unbounded_and_known_off_predictions_remain_possible() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = [EmitterCalibrationObservation {
            measurement: measurement(60.0, 60.0, 0.0, 1), relative_output: &[0.0, 1.0] }];
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        job.advance(cx, usize::MAX).unwrap();
        let fit = job.result().unwrap();
        assert_eq!(job.coefficient_count(), 1);
        assert_eq!(fit.scale_enclosures().unwrap()[0], Interval::new(0.0, f64::INFINITY));
        assert!(fit.flux_scale_bound(0).is_none());
        assert!(fit.flux_scale_bound(1).is_some());
        assert_eq!(job.check_observation(cx, rows[0]).unwrap().compatibility, MeasurementCompatibility::NotExcluded);
        let active = EmitterCalibrationObservation { relative_output: &[1.0, 1.0], ..rows[0] };
        assert!(matches!(job.check_observation(cx, active), Err(FluxCalibrationError::PredictionRange)));
        let dark = [EmitterCalibrationObservation {
            measurement: measurement(3.0, 3.0, 2.0, 2), relative_output: &[0.0, 0.0] }];
        let mut conflict = EmitterCalibrationJob::new(&scene, &dark, evidence(10), Default::default(), cx).unwrap();
        conflict.advance(cx, usize::MAX).unwrap();
        let excluded = conflict.result().unwrap();
        assert_eq!(excluded.termination, EmitterCalibrationTermination::Inconsistent { measurement: 0 });
        assert!(excluded.scale_enclosures().is_none());
        assert!(excluded.flux_scale_bound(0).is_none());
    });
}

#[test]
fn sweep_cap_is_a_sound_outer_box_not_a_convergence_claim() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = mixed_rows();
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10),
            EmitterCalibrationBudget { max_sweeps: 1, ..Default::default() }, cx).unwrap();
        assert_eq!(job.advance(cx, usize::MAX).unwrap(), 4);
        let fit = job.result().unwrap();
        assert_eq!(fit.termination, EmitterCalibrationTermination::SweepLimit);
        assert_eq!(fit.completed_sweeps, 1);
        let gains = fit.scale_enclosures().unwrap();
        bound_contains(gains[0], 0.8); bound_contains(gains[1], 1.2);
        assert!(gains[0].width() > 0.01);
        assert!(fit.flux_scale_bound(0).is_some());
        assert_eq!(job.advance(cx, usize::MAX).unwrap(), 0);
    });
}

#[test]
fn cache_and_contraction_cancel_resume_are_partition_invariant() {
    with_cx(|cx, gate| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = mixed_rows();
        let make = |cx| EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        let mut full = make(cx); full.advance(cx, usize::MAX).unwrap();
        let mut split = make(cx); split.advance(cx, 1).unwrap();
        gate.request();
        assert!(split.advance(cx, 1).is_err());
        assert_eq!(split.prepared_measurements(), 1);
        assert!(split.result().is_none());
        with_cx(|resume, resume_gate| {
            split.advance(resume, 2).unwrap(); // finish cache, then one contracted row
            let saved = split.scales.clone();
            let cursor = split.next_row;
            resume_gate.request();
            assert!(split.advance(resume, 1).is_err());
            assert_eq!(split.scales, saved); assert_eq!(split.next_row, cursor);
            with_cx(|last, _| {
                while !split.is_complete() { split.advance(last, 3).unwrap(); }
            });
        });
        let a = full.result().unwrap(); let b = split.result().unwrap();
        assert_eq!(a.scale_enclosures(), b.scale_enclosures());
        assert_eq!(a.termination, b.termination);
        assert_eq!(a.completed_sweeps, b.completed_sweeps);
        assert_eq!(a.evidence_identity, b.evidence_identity);
    });
}

#[test]
fn sparse_storage_and_input_limits_refuse_without_publishing_partial_rows() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = [
            EmitterCalibrationObservation { measurement: measurement(80.0, 80.0, 0.0, 1), relative_output: &[1.0, 0.0] },
            EmitterCalibrationObservation { measurement: measurement(60.0, 60.0, 0.0, 2), relative_output: &[0.0, 1.0] },
        ];
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10),
            EmitterCalibrationBudget { max_coefficients: 1, ..Default::default() }, cx).unwrap();
        job.advance(cx, 1).unwrap();
        let identity = job.identity;
        assert!(matches!(job.advance(cx, 1), Err(FluxCalibrationError::EmitterBudgetExceeded)));
        assert_eq!(job.prepared_measurements(), 1); assert_eq!(job.coefficient_count(), 1);
        assert_eq!(job.identity, identity); assert!(job.result().is_none());
        for output in [&[f64::NAN, 1.0][..], &[1.1, 0.0], &[-0.1, 1.0], &[1.0], &[f64::INFINITY, 0.0]] {
            let invalid = [EmitterCalibrationObservation { relative_output: output, ..rows[0] }];
            assert!(EmitterCalibrationJob::new(&scene, &invalid, evidence(10), Default::default(), cx).is_err());
        }
        for budget in [
            EmitterCalibrationBudget { max_sources: 1, ..Default::default() },
            EmitterCalibrationBudget { max_measurements: 1, ..Default::default() },
            EmitterCalibrationBudget { max_sweeps: 0, ..Default::default() },
        ] { assert!(EmitterCalibrationJob::new(&scene, &rows, evidence(10), budget, cx).is_err()); }
        assert!(EmitterCalibrationJob::new(&scene, &rows, evidence(0), Default::default(), cx).is_err());
        assert!(job.set_grid_enclosure_budget(GridEnclosureBudget::default()).is_err());
    });
}

#[test]
fn positive_extremes_preserve_zero_support_and_do_not_subtract_infinities() {
    with_cx(|cx, _| {
        let unknown = Interval::new(0.0, f64::INFINITY);
        assert_eq!(positive_product(Interval::point(0.0), unknown), Interval::point(0.0));
        let tiny = f64::from_bits(1);
        let subnormal = positive_product(Interval::point(tiny), Interval::point(0.5));
        assert!(subnormal.lo() == 0.0 && subnormal.hi() > 0.0);
        let overflow = positive_product(Interval::point(f64::MAX), Interval::point(2.0));
        assert!(overflow.lo().is_finite() && overflow.hi().is_infinite());
        let entries = [Coefficient { source: 0, lux: Interval::point(2.0) },
                       Coefficient { source: 1, lux: Interval::point(1.0) }];
        let mut candidate = [Interval::point(0.0); 2];
        let mut suffix = [Interval::point(0.0); 3];
        assert_eq!(contract_row(&entries, measurement(4.0, 4.0, 0.0, 1),
            &[unknown; 2], &mut candidate, &mut suffix, cx).unwrap(), Some(true));
        bound_contains(candidate[0], 2.0); bound_contains(candidate[1], 4.0);
        assert!(candidate.iter().all(|g| finite(*g)));
        assert!(contract_row(&[], measurement(0.0, 0.0, 1.0, 2),
            &[unknown; 2], &mut candidate, &mut suffix, cx).unwrap().is_none());
    });
}

#[test]
fn controlled_outputs_are_bound_to_fit_and_holdout_identities() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = mixed_rows();
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        job.advance(cx, usize::MAX).unwrap();
        let a = job.check_observation(cx, rows[0]).unwrap();
        let b = job.check_observation(cx, EmitterCalibrationObservation {
            relative_output: &[1.0, 1.0], ..rows[0] }).unwrap();
        assert_ne!(a.observation_identity, b.observation_identity);
        assert_eq!(a.calibration_identity, b.calibration_identity);
        let changed_rows = [EmitterCalibrationObservation { relative_output: &[1.0, 0.25], ..rows[0] }, rows[1]];
        let mut changed = EmitterCalibrationJob::new(&scene, &changed_rows, evidence(10), Default::default(), cx).unwrap();
        changed.advance(cx, usize::MAX).unwrap();
        assert_ne!(job.result().unwrap().evidence_identity, changed.result().unwrap().evidence_identity);
    });
}

#[test]
fn independent_emitter_calibration_reuses_native_spectral_haze_and_receiving_floor() {
    use crate::charts::TriMesh;
    use crate::photometry::source::{SPECTRAL_BINS, photopic_weight_nm, wavelength_nm};
    use crate::photometry::transport::PhotometricMedium;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material, Primitive, Shape};
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
    with_cx(|cx, _| {
        let mut bins = [0.0; SPECTRAL_BINS]; bins[40] = 1.0;
        let spectrum = SourceSpectrum::try_from_bins(bins).unwrap();
        let make = |z, intensity| PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, z]).unwrap(),
            4.0 * std::f64::consts::PI * intensity, 1.0, AngularDistribution::isotropic(), spectrum.clone()).unwrap();
        let sources = [make(3.0, 900.0), make(6.0, 1800.0)];
        let density = [1.0, 2.0, 4.0];
        let grid = AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1, 1, 3], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 0.25]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (780.0, 0.4)]).unwrap(), None,
            AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 128 }).unwrap();
        let floor = [Primitive {
            shape: Shape::Mesh(TriMesh::new(vec![[-2.0,-2.0,0.0],[2.0,-2.0,0.0],[0.0,2.0,0.0]], vec![[0,1,2]])),
            material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.3).unwrap() }, emission: None,
        }];
        let scene = PhotometricScene::try_with_medium(&floor, &sources, PhotometricMedium::Absorbing(&grid), TransportOptions::default(), cx).unwrap();
        // Independent single-bin reference: 1.75 kg/m² through all layers.
        let transfer = 683.0 * spectrum.bins_w_per_lumen_nm()[40] * photopic_weight_nm(40)
            * (-1.75 * (wavelength_nm(40) - 380.0) / 1000.0).exp();
        let rows = [
            EmitterCalibrationObservation { measurement: measurement(80.0*transfer-1e-7, 80.0*transfer+1e-7, 0.0, 1), relative_output: &[1.0,0.0] },
            EmitterCalibrationObservation { measurement: measurement(60.0*transfer-1e-7, 60.0*transfer+1e-7, 0.0, 2), relative_output: &[0.0,1.0] },
        ];
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        job.advance(cx, usize::MAX).unwrap();
        let fit = job.result().unwrap();
        for (g, truth) in fit.scale_enclosures().unwrap().iter().zip([0.8,1.2]) {
            bound_contains(*g, truth); assert!(g.width() < 1e-6);
        }
        let check = job.check_observation(cx, rows[0]).unwrap();
        assert_eq!(check.compatibility, MeasurementCompatibility::NotExcluded);
        assert_eq!(check.reference.uncertain_visibility_sources, Some(0));
        assert!(check.reference.estimate.total_error_bound_lux.is_none());
        let mut refused = EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        refused.set_grid_enclosure_budget(GridEnclosureBudget { max_nodes_per_segment: 1 }).unwrap();
        assert!(matches!(refused.advance(cx, 1), Err(FluxCalibrationError::NumericalUnavailable { .. })));
        assert_eq!(refused.prepared_measurements(), 0); assert!(refused.result().is_none());
    });
}

#[test]
fn contradictory_positive_readings_are_not_silently_averaged_or_discarded() {
    with_cx(|cx, _| {
        let sources = [fixture(3.0, 900.0), fixture(6.0, 1800.0)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let rows = [
            EmitterCalibrationObservation { measurement: measurement(79.9, 80.1, 0.0, 1), relative_output: &[1.0, 0.0] },
            EmitterCalibrationObservation { measurement: measurement(99.9, 100.1, 0.0, 2), relative_output: &[1.0, 0.0] },
        ];
        let mut job = EmitterCalibrationJob::new(&scene, &rows, evidence(10), Default::default(), cx).unwrap();
        job.advance(cx, 3).unwrap(); // cache both rows, then accept only first contraction
        let accepted = job.scales.clone();
        job.advance(cx, 1).unwrap();
        assert_eq!(job.scales, accepted); // no partially contracted conflicting vector
        let fit = job.result().unwrap();
        assert_eq!(fit.termination, EmitterCalibrationTermination::Inconsistent { measurement: 1 });
        assert!(fit.scale_enclosures().is_none());
        assert!(matches!(job.check_observation(cx, rows[0]), Err(FluxCalibrationError::NoFiniteFit)));
    });
}
