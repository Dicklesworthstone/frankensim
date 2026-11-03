use super::*;
use crate::conductor::{CoatingIorSample, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorIorSample, ConductorSource};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::surface_detail::BrushingPattern;
use crate::tracer::{Material, opaque_bsdf_eval};

fn medium() -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(1.33).unwrap(), BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
}

fn optics() -> [ConductorOptics; 3] {
    let source = ConductorSource::try_new(fs_blake3::ContentHash([98; 32]), ConductorDataStatus::Representative).unwrap();
    let bare = ConductorOptics::try_new(core::array::from_fn(|i|
        ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()), source).unwrap();
    let transparent = bare.with_coating(CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(), 125e-9, source).unwrap());
    let absorbing = bare.with_absorbing_coating(CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.45, 0.07).unwrap()),
        100e-9, source).unwrap());
    [bare, transparent, absorbing]
}

fn beam(wo: Vec3, state: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(Point3::new(0.4, 0.0, 0.0),
        PolarizationFrame::try_new(wo.scale(-1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        550.0, medium(), state).unwrap()
}

fn source_frame(wi: Vec3) -> PolarizationFrame {
    PolarizationFrame::try_new(wi.scale(-1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-9 * a.abs().max(b.abs()).max(1e-20), "{a:.17e} != {b:.17e}");
}

#[test]
fn native_metal_connections_retain_sampling_and_coherent_polarized_power() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let mut accepted = 0;
    for law in optics() {
        for alpha in [0.05, 0.3, 1.0] {
            let surface = ConductorSurface::try_rough(alpha).unwrap();
            let material = Material::Conductor { optics: law, surface };
            for state in [Stokes::unpolarized(3.0), Stokes { i: 2.0, q: 0.4, u: 0.7, v: 0.9 }] {
                let camera = beam(Vec3::new(0.0, 0.6, 0.8), state);
                for draw in [[0.0, 0.0], [0.17, 0.31], [0.8, 0.7]] {
                    let Some(sample) = camera.sample_rough_conductor(normal, &law, surface, draw).unwrap() else { continue; };
                    accepted += 1;
                    let wi = sample.beam.frame().direction();
                    let response = camera.evaluate_rough_conductor_connection(normal, &law, surface,
                        source_frame(wi), Stokes::unpolarized(1.0)).unwrap();
                    near(response.pdf, sample.pdf);
                    near(response.factor * normal.dot(wi) / response.pdf,
                        sample.beam.power().i / camera.power().i);
                    assert_eq!(response.event, DielectricEvent::Reflection);
                    if state.q == 0.0 {
                        let scalar = opaque_bsdf_eval(&material, normal, camera.frame().direction().scale(-1.0),
                            wi, 550.0, Some(medium()), None).unwrap();
                        near(response.factor, scalar);
                    }
                }
            }
        }
    }
    assert!(accepted > 20);
}

#[test]
fn brushed_connections_use_geometric_irradiance_and_native_camera_radiance() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let transform = RigidTransform::identity();
    let point = Point3::new(0.4, 0.0, 0.0);
    let law = optics()[2];
    let mut admitted = 0;
    let mut distinct_adjoint = false;
    for tilt in [0.0, 0.2, 0.35] {
        let detail = SurfaceDetail::try_new(0.2, 3.0, BrushingPattern::Radial, tilt, None).unwrap();
        let material = Material::BrushedConductor { optics: law, detail };
        for cosine in [0.005_f64, 0.4, 1.0] {
            let wo = Vec3::new(0.0, (1.0 - cosine * cosine).sqrt(), cosine);
            let camera = beam(wo, Stokes::unpolarized(2.0));
            let frame = instanced_frame(camera, &detail, transform, point, normal).unwrap().unwrap();
            for draw in [[0.0, 0.0], [0.13, 0.27], [0.6, 0.9], [0.99, 0.4]] {
                let sample = camera.sample_instanced_brushed_conductor_for_transport(
                    &law, &detail, transform, point, normal, draw, true,
                ).unwrap();
                let Some(sample) = sample else { continue; };
                admitted += 1;
                let wi = sample.beam.frame().direction();
                let response = camera.evaluate_instanced_brushed_conductor_connection(&law, &detail,
                    transform, point, normal, source_frame(wi), Stokes::unpolarized(1.0)).unwrap();
                near(response.pdf, sample.pdf);
                near(response.factor * normal.dot(wi) / response.pdf, sample.beam.power().i / camera.power().i);
                let scalar = opaque_bsdf_eval(&material, normal, wo, wi, 550.0, Some(medium()), Some(&frame)).unwrap();
                near(response.factor, scalar);
                let forward = camera.sample_instanced_brushed_conductor(&law, &detail, transform, point, normal, draw)
                    .unwrap().unwrap();
                distinct_adjoint |= (forward.beam.power().i - sample.beam.power().i).abs() > 1e-4;
            }
        }
    }
    assert!(admitted > 0 && distinct_adjoint);
}

#[test]
fn brushed_poles_and_zero_detail_reuse_isotropic_and_delta_stays_discrete() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let camera = beam(Vec3::new(0.0, 0.6, 0.8), Stokes { i: 2.0, q: 0.4, u: 0.7, v: 0.9 });
    let law = optics()[1];
    let wi = Vec3::new(0.0, -0.6, 0.8);
    let frame = source_frame(wi);
    let state = Stokes { i: 1.0, q: 0.2, u: -0.3, v: 0.5 };
    let reference = camera.evaluate_rough_conductor_connection(normal, &law,
        ConductorSurface::try_rough(0.2).unwrap(), frame, state).unwrap();
    for (detail, point) in [
        (SurfaceDetail::none(0.2), camera.ray().origin),
        (SurfaceDetail::try_new(0.2, 3.0, BrushingPattern::Radial, 0.2, None).unwrap(), Point3::new(0.0, 0.0, 2.0)),
    ] {
        let result = camera.evaluate_instanced_brushed_conductor_connection(&law, &detail,
            RigidTransform::identity(), point, normal, frame, state).unwrap();
        near(reference.factor, result.factor);
        assert_eq!(reference.pdf, result.pdf);
    }
    let delta = camera.evaluate_rough_conductor_connection(normal, &law, ConductorSurface::SMOOTH,
        frame, state).unwrap();
    assert_eq!(delta.factor, 0.0);
    assert_eq!(delta.pdf, 0.0);
    let below = camera.evaluate_rough_conductor_connection(normal, &law, ConductorSurface::try_rough(0.2).unwrap(),
        source_frame(wi.scale(-1.0)), state).unwrap();
    assert_eq!(below.factor, 0.0);
    assert_eq!(below.pdf, 0.0);
}

#[test]
fn metal_connection_preserves_source_ellipticity_under_independent_frame_changes() {
    let camera = beam(Vec3::new(0.0, 0.6, 0.8), Stokes { i: 2.0, q: 0.4, u: 0.7, v: 0.9 });
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let law = optics()[2];
    let surface = ConductorSurface::try_rough(0.2).unwrap();
    let sample = camera.sample_rough_conductor(normal, &law, surface, [0.17, 0.31]).unwrap().unwrap();
    let frame = source_frame(sample.beam.frame().direction());
    let state = Stokes { i: 1.0, q: 0.1, u: -0.4, v: 0.7 };
    let reference = camera.evaluate_rough_conductor_connection(normal, &law, surface, frame, state).unwrap();
    let changed_camera = camera.with_reference_axis(Vec3::new(0.3, 0.7, 0.9)).unwrap();
    let changed_frame = PolarizationFrame::try_new(frame.direction(), Vec3::new(0.8, 0.4, 0.2)).unwrap();
    let changed_state = frame.express_in(state, changed_frame);
    let changed = changed_camera.evaluate_rough_conductor_connection(normal, &law, surface, changed_frame, changed_state).unwrap();
    near(reference.factor, changed.factor);
    assert_eq!(reference.pdf, changed.pdf);
    let intensity_only = camera.evaluate_rough_conductor_connection(normal, &law, surface,
        frame, Stokes::unpolarized(1.0)).unwrap();
    assert!((reference.factor - intensity_only.factor).abs() > reference.factor * 1e-3);
}
