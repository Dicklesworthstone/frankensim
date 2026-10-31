use super::*;
use asupersync::types::Budget;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::polarization::photometry::SpectralOptics;
use super::super::{SurfaceBinding, TraceSettings};
use super::super::super::transport::{PolarizationFrame, PolarizedBeam};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 31, kernel_id: 7, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)));
}
fn air() -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
}
fn plane(id: u64, z: f64, half_size: f64) -> GeometryInstance {
    let label = format!("{id}:{z}:{half_size}");
    GeometryInstance::try_new(id, hash_domain("polarized-meter-test", label.as_bytes()),
        SharedGeometry::mesh(TriMesh::new(vec![[-half_size, -half_size, z],
            [half_size, -half_size, z], [0.0, half_size, z]], vec![[0, 1, 2]])),
        RigidTransform::identity()).unwrap()
}
fn detector(id: u64) -> SurfaceBinding {
    SurfaceBinding { object_id: id,
        interaction: Interaction::Detector { reference_axis: Vec3::new(1.0, 0.0, 0.0) } }
}
fn settings() -> TraceSettings {
    TraceSettings { max_surface_events: 8, max_distance_m: 10.0,
        surface_epsilon_m: 1e-9, query_offset_m: 1e-7 }
}
fn ray(bin: usize, x: f64) -> PolarizedBeam {
    PolarizedBeam::try_new(Point3::new(x, 0.0, 0.0),
        PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        wavelength_nm(bin), air(), Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap()).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 2e-9, "{a} != {b}"); }

#[test]
fn g1_calibrated_source_optics_native_geometry_and_polarimeter_close_flux() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0, 10.0)]).unwrap();
        let media = [air()];
        let scene = PolarizedScene::try_new(&geometry, &media, vec![detector(7)]).unwrap();
        let mut meter = scene.power_meter(7).unwrap();
        let source = PhotometricSource::try_new(RigidTransform::identity(), 1000.0, 1.0,
            AngularDistribution::uniform_cone(0.1).unwrap(),
            SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap()).unwrap()
            .with_unpolarized_optics()
            .through(&SpectralOptics::linear_polarizer(0.0).unwrap()).unwrap();
        // The entire cone hits this receiver; its uniform source profile admits
        // this one-node angular quadrature for the direct, lossless fixture.
        let solid_angle = 2.0 * core::f64::consts::PI * (1.0 - det::cos(0.1));
        for bin in 0..SPECTRAL_BINS {
            let beam = source.emit_beam(Vec3::new(0.0, 0.0, 1.0), bin, solid_angle, air()).unwrap();
            let path = scene.trace(cx, beam, 0, settings(), 123, bin as u64).unwrap();
            meter.record(bin, &path).unwrap();
        }
        let reading = meter.reading().unwrap().require_complete().unwrap();
        near(reading.photopic_stokes_lm().i, 500.0);
        near(reading.linear_analyzer(0.0).unwrap().luminous_lm, 500.0);
        near(reading.linear_analyzer(core::f64::consts::FRAC_PI_2).unwrap().luminous_lm, 0.0);
        near(reading.circular_analyzer(true).unwrap().luminous_lm, 250.0);
        assert_eq!(reading.launched_paths(), SPECTRAL_BINS as u64);
        assert_eq!(reading.detected_paths(), reading.launched_paths());
    });
}

#[test]
fn g1_occluded_launches_remain_in_each_bin_normalizer() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0, 10.0), plane(8, 0.5, 0.1)]).unwrap();
        let media = [air()];
        let scene = PolarizedScene::try_new(&geometry, &media,
            vec![detector(7), SurfaceBinding { object_id: 8, interaction: Interaction::Opaque }]).unwrap();
        let mut meter = scene.power_meter(7).unwrap();
        for bin in 0..SPECTRAL_BINS {
            for x in [0.0, 1.0] {
                let path = scene.trace(cx, ray(bin, x), 0, settings(), 0, 0).unwrap();
                meter.record(bin, &path).unwrap();
            }
            near(meter.spectral_mean_w_nm(bin).unwrap().i, 0.5);
        }
        let reading = meter.reading().unwrap().require_complete().unwrap();
        near(reading.radiant_stokes_w().i, 0.5 * SPECTRAL_BINS as f64 * BIN_WIDTH_NM);
        assert_eq!(reading.detected_paths(), SPECTRAL_BINS as u64);
        assert_eq!(reading.launched_paths(), 2 * SPECTRAL_BINS as u64);
    });
}

#[test]
fn g0_truncated_paths_remain_visible_through_analysis_and_completeness() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0, 10.0)]).unwrap();
        let media = [air()];
        let scene = PolarizedScene::try_new(&geometry, &media, vec![detector(7)]).unwrap();
        let mut meter = scene.power_meter(7).unwrap();
        let mut limited = settings(); limited.max_distance_m = 0.25;
        for bin in 0..SPECTRAL_BINS {
            for config in [settings(), limited] {
                let path = scene.trace(cx, ray(bin, 0.0), 0, config, 0, 0).unwrap();
                meter.record(bin, &path).unwrap();
            }
        }
        let reading = meter.reading().unwrap();
        near(reading.radiant_stokes_w().i, 0.5 * SPECTRAL_BINS as f64 * BIN_WIDTH_NM);
        assert_eq!(reading.unresolved_paths(), SPECTRAL_BINS as u64);
        assert_eq!(reading.require_complete().unwrap_err(),
            MeterError::UnresolvedPaths { count: SPECTRAL_BINS as u64 });
        assert_eq!(reading.linear_analyzer(0.0).unwrap().unresolved_paths, SPECTRAL_BINS as u64);
    });
}

#[test]
fn g0_missing_bins_wrong_wavelength_and_changed_detector_refuse_atomically() {
    with_cx(|cx| {
        let media = [air()];
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0, 10.0)]).unwrap();
        let scene = PolarizedScene::try_new(&geometry, &media, vec![detector(7)]).unwrap();
        let mut meter = scene.power_meter(7).unwrap();
        let original = meter.clone();
        assert_eq!(meter.reading().unwrap_err(), MeterError::MissingSamples { bin: 0 });
        let path = scene.trace(cx, ray(0, 0.0), 0, settings(), 0, 0).unwrap();
        assert_eq!(meter.record(SPECTRAL_BINS, &path), Err(MeterError::InvalidBin));
        assert_eq!(meter.record(1, &path), Err(MeterError::WavelengthMismatch));
        assert_eq!(meter, original);
        let moved_geometry = InstanceScene::try_new(vec![plane(7, 2.0, 10.0)]).unwrap();
        let moved = PolarizedScene::try_new(&moved_geometry, &media, vec![detector(7)]).unwrap();
        let path = moved.trace(cx, ray(0, 0.0), 0, settings(), 0, 0).unwrap();
        assert_eq!(meter.record(0, &path), Err(MeterError::DetectorGeometryMismatch));
        assert_eq!(meter, original);
    });
}

#[test]
fn g0_signed_stokes_sums_retain_small_terms_and_refuse_overflow() {
    let mut sum = Sum::default();
    for value in [1e16, 1.0, -1e16] { sum.add(value).unwrap(); }
    assert_eq!(sum.total(), 1.0);
    let mut sum = Sum::default(); sum.add(f64::MAX).unwrap();
    let before = sum;
    assert_eq!(sum.add(f64::MAX), Err(MeterError::NumericalRange));
    assert_eq!(sum, before);
}

#[test]
fn g1_analyzer_pairs_close_power_including_large_physical_states() {
    for state in [Stokes::try_new(1.0, 0.3, 0.4, 0.5).unwrap(),
        Stokes::try_new(f64::MAX, f64::MAX, 0.0, 0.0).unwrap()] {
        let horizontal = project(state, [1.0, 0.0, 0.0]).unwrap();
        let vertical = project(state, [-1.0, 0.0, 0.0]).unwrap();
        assert!((horizontal/state.i + vertical/state.i - 1.0).abs() < 1e-14);
        let right = project(state, [0.0, 0.0, 1.0]).unwrap();
        let left = project(state, [0.0, 0.0, -1.0]).unwrap();
        assert!((right/state.i + left/state.i - 1.0).abs() < 1e-14);
    }
}
