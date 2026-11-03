use super::*;
use crate::conductor::{CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricCoating, DielectricGlass, GlassProvenance};
use fs_geom::Point3;

fn glass(index: f64) -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(index).unwrap(), BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
}

fn beam(direction: Vec3, medium: DielectricGlass, state: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        550.0, medium, state).unwrap()
}

fn source_frame(direction: Vec3) -> PolarizationFrame {
    PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

fn coated(alpha: f64) -> DielectricSurface {
    let source = ConductorSource::try_new(fs_blake3::ContentHash([97; 32]), ConductorDataStatus::Representative).unwrap();
    let cell = [
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.1).unwrap(), 0.04, 82e-9).unwrap(),
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.35).unwrap(), 0.01, 135e-9).unwrap(),
    ];
    let coating = DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cell, 2, source).unwrap()).unwrap();
    DielectricSurface::try_rough(alpha).unwrap().with_coating(coating)
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * a.abs().max(b.abs()).max(1e-20), "{a:.17e} != {b:.17e}");
}

#[test]
fn rough_connections_match_native_proposals_and_radiance_on_both_sides() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let mut reflections = 0;
    let mut transmissions = 0;
    for (eta_i, eta_t, sign) in [(1.0, 1.5, -1.0), (1.5, 1.0, 1.0)] {
        for surface in [DielectricSurface::try_rough(0.2).unwrap(), coated(0.2)] {
            for state in [Stokes::unpolarized(3.0), Stokes { i: 2.0, q: 0.4, u: 0.7, v: 0.9 }] {
                let camera = beam(Vec3::new(0.6 * sign, 0.0, 0.8 * sign), glass(eta_i), state);
                for draw in [[0.0, 0.0, 0.0], [0.17, 0.31, 0.0], [0.17, 0.31, 0.999], [0.4, 0.8, 0.999]] {
                    let Some(sample) = camera.sample_interface(normal, glass(eta_t), surface, draw).unwrap() else { continue; };
                    assert!(!sample.delta);
                    let wi = sample.beam.frame().direction();
                    let response = camera.evaluate_dielectric_connection(normal, glass(eta_t), surface,
                        source_frame(wi.scale(-1.0)), Stokes::unpolarized(1.0)).unwrap();
                    assert_eq!(response.event, sample.event);
                    near(response.pdf, sample.pdf);
                    let eta_radiance = if sample.event == DielectricEvent::Transmission {
                        transmissions += 1;
                        (eta_i / eta_t) * (eta_i / eta_t)
                    } else {
                        reflections += 1;
                        1.0
                    };
                    let expected = sample.beam.power().i / camera.power().i * sample.pdf
                        / normal.dot(wi).abs() * eta_radiance;
                    near(response.factor, expected);
                    // A camera using native power continuation and a terminal
                    // n_camera/n_light conversion must agree with evaluated
                    // radiance using n_camera/n_camera-arrival exactly once.
                    let n_camera = 1.33;
                    let n_light = if sample.event == DielectricEvent::Transmission { eta_t } else { eta_i };
                    let connection = response.factor * (n_camera / eta_i) * (n_camera / eta_i);
                    let continuation = sample.beam.power().i / camera.power().i * sample.pdf
                        / normal.dot(wi).abs() * (n_camera / n_light) * (n_camera / n_light);
                    near(connection, continuation);
                }
            }
        }
    }
    assert!(reflections > 0 && transmissions > 0);
}

#[test]
fn polarized_prefix_and_source_are_both_retained_at_brewster_reflection() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let cosine = 1.0 / 3.25_f64.sqrt();
    let sine = 1.5 * cosine;
    let camera = beam(Vec3::new(-sine, 0.0, -cosine), glass(1.0), Stokes { i: 1.0, q: 1.0, u: 0.0, v: 0.0 });
    let frame = source_frame(Vec3::new(sine, 0.0, -cosine));
    let surface = DielectricSurface::try_rough(0.3).unwrap();
    let s = camera.evaluate_dielectric_connection(normal, glass(1.5), surface, frame,
        Stokes { i: 1.0, q: 1.0, u: 0.0, v: 0.0 }).unwrap();
    let p = camera.evaluate_dielectric_connection(normal, glass(1.5), surface, frame,
        Stokes { i: 1.0, q: -1.0, u: 0.0, v: 0.0 }).unwrap();
    let unpolarized = beam(camera.frame().direction(), glass(1.0), Stokes::unpolarized(1.0))
        .evaluate_dielectric_connection(normal, glass(1.5), surface, frame, Stokes::unpolarized(1.0)).unwrap();
    assert!(s.factor > 0.0);
    near(s.factor, 4.0 * unpolarized.factor);
    assert!(p.factor <= s.factor * 1e-25);
    assert_eq!(s.pdf, p.pdf);
    assert!(p.pdf > 0.0, "physical polarization extinction must not alter the scalar proposal");
}

#[test]
fn coherent_connections_ignore_reference_axes_and_prefix_intensity() {
    let camera = beam(Vec3::new(-0.6, 0.0, -0.8), glass(1.0), Stokes { i: 2.0, q: 0.4, u: 0.7, v: 0.9 });
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let surface = coated(0.2);
    let sample = camera.sample_interface(normal, glass(1.5), surface, [0.17, 0.31, 0.999]).unwrap().unwrap();
    let source = source_frame(sample.beam.frame().direction().scale(-1.0));
    let state = Stokes { i: 1.0, q: 0.1, u: 0.4, v: -0.6 };
    let reference = camera.evaluate_dielectric_connection(normal, glass(1.5), surface, source, state).unwrap();
    let changed_camera = camera.with_reference_axis(Vec3::new(0.3, 0.7, 0.9)).unwrap();
    let changed_source = PolarizationFrame::try_new(source.direction(), Vec3::new(0.8, 0.4, 0.2)).unwrap();
    let changed_state = source.express_in(state, changed_source);
    let changed = changed_camera.evaluate_dielectric_connection(normal, glass(1.5), surface,
        changed_source, changed_state).unwrap();
    near(reference.factor, changed.factor);
    assert_eq!(reference.pdf, changed.pdf);
    let unit_camera = normalized_camera(camera, state).unwrap();
    let unit_response = unit_camera.evaluate_dielectric_connection(normal, glass(1.5), surface, source, state).unwrap();
    near(reference.factor, unit_response.factor);
}

#[test]
fn circular_reciprocity_keeps_v_and_reverses_only_u() {
    let camera = beam(Vec3::new(0.0, 0.0, 1.0), glass(1.0), Stokes { i: 1.0, q: 0.0, u: 0.0, v: 1.0 });
    let source = source_frame(Vec3::new(0.0, 0.0, -1.0));
    let same = Stokes { i: 1.0, q: 0.0, u: 0.0, v: 1.0 };
    let opposite = Stokes { v: -1.0, ..same };
    assert_eq!(reciprocal_contraction(camera, source, same).unwrap(), 2.0);
    assert_eq!(reciprocal_contraction(camera, source, opposite).unwrap(), 0.0);
    let diagonal = beam(Vec3::new(0.0, 0.0, 1.0), glass(1.0), Stokes { i: 1.0, q: 0.0, u: 1.0, v: 0.0 });
    assert_eq!(reciprocal_contraction(diagonal, source, diagonal.power()).unwrap(), 0.0);
}

#[test]
fn smooth_and_equal_index_nulls_have_no_continuous_connection() {
    let camera = beam(Vec3::new(0.0, 0.0, -1.0), glass(1.0), Stokes::unpolarized(1.0));
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let source = source_frame(Vec3::new(0.0, 0.0, 1.0));
    for (surface, other) in [(DielectricSurface::SMOOTH, glass(1.5)), (DielectricSurface::try_rough(0.2).unwrap(), glass(1.0))] {
        let response = camera.evaluate_dielectric_connection(normal, other, surface, source, Stokes::unpolarized(1.0)).unwrap();
        assert_eq!(response.factor, 0.0);
        assert_eq!(response.pdf, 0.0);
    }
    assert!(camera.evaluate_dielectric_connection(normal, glass(1.0), coated(0.2), source,
        Stokes::unpolarized(1.0)).is_err(), "coated equal-index rough mixture retains native refusal");
    assert!(camera.evaluate_dielectric_connection(normal, glass(1.5), DielectricSurface::try_rough(0.2).unwrap(), source,
        Stokes::unpolarized(2.0)).is_err());
}
