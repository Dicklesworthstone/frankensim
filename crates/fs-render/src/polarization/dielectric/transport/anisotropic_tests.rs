//! Polarized power and reciprocal connections through the authored GGX law.

use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoherentDielectricCoating, CoherentPeriodicCoating, ConductorDataStatus,
    ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricCoating, GlassProvenance};

fn glass(eta: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(eta).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn beam(direction: Vec3, eta: f64, power: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        550.0,
        glass(eta),
        power,
    )
    .unwrap()
}

fn rough() -> DielectricSurface {
    DielectricSurface::try_anisotropic_rough(0.17, 0.63, Vec3::new(0.6, 0.8, 0.2)).unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 8e-11 * scale.max(a.abs()).max(b.abs()).max(1e-20),
        "{a:.17e} != {b:.17e}"
    );
}

fn same_stokes(a: Stokes, b: Stokes) {
    for (a, b) in [a.i, a.q, a.u, a.v].into_iter().zip([b.i, b.q, b.u, b.v]) {
        close(a, b, 1.0);
    }
}

// Independent directional Walter D and Smith G1, evaluated in the authored
// tangent plane. This test does not call the production anisotropic kernels.
fn dg(normal: Vec3, h: Vec3, wo: Vec3, wi: Vec3) -> (f64, f64, f64) {
    let t = Vec3::new(0.6, 0.8, 0.0);
    let b = cross(normal, t);
    let d = 1.0
        / (core::f64::consts::PI
            * 0.17
            * 0.63
            * ((h.dot(t) / 0.17).powi(2) + (h.dot(b) / 0.63).powi(2) + h.dot(normal).powi(2))
                .powi(2));
    let g = |v: Vec3| {
        let z = v.dot(normal).abs();
        2.0 / (1.0
            + (1.0 + ((0.17 * v.dot(t)).powi(2) + (0.63 * v.dot(b)).powi(2)) / (z * z)).sqrt())
    };
    (d, g(wo), g(wi))
}

#[test]
fn anisotropic_beam_samples_match_independent_walter_density_and_full_stokes_power() {
    let geometric = Vec3::new(0.0, 0.0, 1.0);
    let mut count = [0; 2];
    for (eta_i, eta_t, side) in [(1.0, 1.5, 1.0), (1.5, 1.0, -1.0)] {
        let normal = geometric.scale(side);
        let wo = unit(Vec3::new(0.35, -0.15, side)).unwrap();
        for power in [
            Stokes::unpolarized(1.0),
            Stokes {
                i: 1.0,
                q: 0.3,
                u: -0.4,
                v: 0.5,
            },
        ] {
            let input = beam(wo.scale(-1.0), eta_i, power);
            for draw in [[0.01, 0.2, 0.0], [0.07, 0.39, 0.0], [0.07, 0.39, 0.999]] {
                let Some(sample) = input
                    .sample_interface(geometric, glass(eta_t), rough(), draw)
                    .unwrap()
                else {
                    continue;
                };
                assert!(!sample.delta);
                let wi = sample.beam.frame().direction();
                let transmission = sample.event == DielectricEvent::Transmission;
                count[usize::from(transmission)] += 1;
                let eta = if transmission { eta_t / eta_i } else { 1.0 };
                let mut h = unit(Vec3::new(
                    wo.x + eta * wi.x,
                    wo.y + eta * wi.y,
                    wo.z + eta * wi.z,
                ))
                .unwrap();
                if h.dot(normal) < 0.0 {
                    h = h.scale(-1.0);
                }
                let (d, go, gi) = dg(normal, h, wo, wi);
                let a = wo.dot(h);
                let b = wi.dot(h);
                let fresnel = fresnel_dielectric_polarized(a, eta_i, eta_t).unwrap();
                let (matrix, jacobian, geometry) = if transmission {
                    let denominator = (a + eta * b).powi(2);
                    (
                        fresnel.transmission,
                        eta * eta * b.abs() / denominator,
                        d * go * gi * eta * eta * (a * b).abs()
                            / (normal.dot(wo) * denominator * sample.pdf),
                    )
                } else {
                    (
                        fresnel.reflection,
                        1.0 / (4.0 * a.abs()),
                        d * go * gi / (4.0 * normal.dot(wo) * sample.pdf),
                    )
                };
                let density = d * go * a.abs() / normal.dot(wo) * jacobian * matrix.0[0][0];
                close(sample.pdf, density, 0.0);
                let s = unit(cross(input.frame().direction(), h)).unwrap();
                let incident = PolarizationFrame::try_new(input.frame().direction(), s).unwrap();
                let local = input.frame().express_in(input.power(), incident);
                let expected = matrix.apply(local);
                let expected = Stokes {
                    i: expected.i * geometry,
                    q: expected.q * geometry,
                    u: expected.u * geometry,
                    v: expected.v * geometry,
                };
                let outgoing = PolarizationFrame::try_new(wi, s).unwrap();
                same_stokes(
                    sample.beam.power(),
                    outgoing.express_in(expected, sample.beam.frame()),
                );
                assert_eq!(
                    sample.beam.medium(),
                    if transmission {
                        glass(eta_t).into()
                    } else {
                        input.medium()
                    }
                );
            }
        }
    }
    assert!(count.into_iter().all(|count| count > 0));
}

#[test]
fn normal_facet_keeps_analytic_polarization_and_anisotropic_probability() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let input = beam(
        n.scale(-1.0),
        1.0,
        Stokes {
            i: 1.0,
            q: 0.3,
            u: 0.4,
            v: 0.5,
        },
    );
    let surface =
        DielectricSurface::try_anisotropic_rough(0.11, 0.73, Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let reflected = input
        .sample_interface(n, glass(1.5), surface, [0.0, 0.0, 0.0])
        .unwrap()
        .unwrap();
    assert_eq!(reflected.event, DielectricEvent::Reflection);
    same_stokes(
        reflected.beam.power(),
        Stokes {
            i: 1.0,
            q: 0.3,
            u: -0.4,
            v: -0.5,
        },
    );
    close(
        reflected.pdf,
        0.04 / (4.0 * core::f64::consts::PI * 0.11 * 0.73),
        0.0,
    );
    let transmitted = input
        .sample_interface(n, glass(1.5), surface, [0.0, 0.0, 0.999])
        .unwrap()
        .unwrap();
    assert_eq!(transmitted.event, DielectricEvent::Transmission);
    same_stokes(transmitted.beam.power(), input.power());
    close(
        transmitted.pdf,
        0.96 * 9.0 / (core::f64::consts::PI * 0.11 * 0.73),
        0.0,
    );
}

#[test]
fn anisotropic_connections_keep_polarized_sources_and_reciprocal_radiance() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    for (eta_i, eta_t, side) in [(1.0, 1.5, -1.0), (1.5, 1.0, 1.0)] {
        let input = beam(
            Vec3::new(0.3, -0.2, side),
            eta_i,
            Stokes {
                i: 2.0,
                q: 0.4,
                u: 0.7,
                v: 0.9,
            },
        );
        for event in [0.0, 0.999] {
            let Some(sample) = input
                .sample_interface(n, glass(eta_t), rough(), [0.04, 0.3, event])
                .unwrap()
            else {
                continue;
            };
            let wi = sample.beam.frame().direction();
            let source =
                PolarizationFrame::try_new(wi.scale(-1.0), Vec3::new(0.2, 1.0, 0.3)).unwrap();
            let response = input
                .evaluate_dielectric_connection(
                    n,
                    glass(eta_t),
                    rough(),
                    source,
                    Stokes::unpolarized(1.0),
                )
                .unwrap();
            let radiance = if sample.event == DielectricEvent::Transmission {
                (eta_i / eta_t).powi(2)
            } else {
                1.0
            };
            close(
                response.factor,
                sample.beam.power().i / input.power().i * sample.pdf / n.dot(wi).abs() * radiance,
                0.0,
            );
            close(response.pdf, sample.pdf, 0.0);
            let state = Stokes {
                i: 1.0,
                q: 0.2,
                u: -0.3,
                v: 0.6,
            };
            let original = input
                .evaluate_dielectric_connection(n, glass(eta_t), rough(), source, state)
                .unwrap();
            let changed_source =
                PolarizationFrame::try_new(source.direction(), Vec3::new(0.7, 0.3, 0.2)).unwrap();
            let changed_input = input.with_reference_axis(Vec3::new(0.8, 0.3, 0.2)).unwrap();
            let changed = changed_input
                .evaluate_dielectric_connection(
                    n,
                    glass(eta_t),
                    rough(),
                    changed_source,
                    source.express_in(state, changed_source),
                )
                .unwrap();
            close(original.factor, changed.factor, 0.0);
            assert_eq!(original.pdf, changed.pdf);
        }
    }
    // The mirror facet at Brewster incidence passes only s polarization, even
    // though the anisotropic visible-normal proposal has positive p density.
    let c = 1.0 / 3.25_f64.sqrt();
    let input = beam(
        Vec3::new(-1.5 * c, 0.0, -c),
        1.0,
        Stokes {
            i: 1.0,
            q: 1.0,
            u: 0.0,
            v: 0.0,
        },
    );
    let source =
        PolarizationFrame::try_new(Vec3::new(1.5 * c, 0.0, -c), Vec3::new(0.0, 1.0, 0.0)).unwrap();
    let s = input
        .evaluate_dielectric_connection(
            n,
            glass(1.5),
            rough(),
            source,
            Stokes {
                i: 1.0,
                q: 1.0,
                u: 0.0,
                v: 0.0,
            },
        )
        .unwrap();
    let p = input
        .evaluate_dielectric_connection(
            n,
            glass(1.5),
            rough(),
            source,
            Stokes {
                i: 1.0,
                q: -1.0,
                u: 0.0,
                v: 0.0,
            },
        )
        .unwrap();
    assert!(s.factor > 0.0 && p.factor < s.factor * 1e-24);
    assert!(p.pdf > 0.0);
    assert_eq!(p.pdf, s.pdf);
}

#[test]
fn original_film_stack_keeps_anisotropy_order_phase_and_survival() {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([83; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let laws = [(1.9, 83e-9), (1.25, 127e-9)];
    let films = laws.map(|(eta, thickness)| {
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(eta).unwrap(), thickness, source)
            .unwrap()
    });
    let cells = laws.map(|(eta, thickness)| {
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(eta).unwrap(), 0.0, thickness).unwrap()
    });
    let compact = rough().with_coating(
        DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cells, 1, source).unwrap())
            .unwrap(),
    );
    let stack = DielectricCoatingStack::try_new(&films).unwrap();
    let n = Vec3::new(0.0, 0.0, 1.0);
    for (eta_i, eta_t, side) in [(1.0, 1.5, -1.0), (1.5, 1.0, 1.0)] {
        let input = beam(
            Vec3::new(0.25, -0.15, side),
            eta_i,
            Stokes {
                i: 1.0,
                q: 0.2,
                u: 0.3,
                v: 0.6,
            },
        );
        for event in [0.0, 0.999] {
            let draw = [0.03, 0.4, event];
            let actual = input
                .sample_stack_interface(n, glass(eta_t), rough(), stack, draw)
                .unwrap()
                .unwrap();
            let expected = input
                .sample_interface(n, glass(eta_t), compact, draw)
                .unwrap()
                .unwrap();
            same_stokes(actual.beam.power(), expected.beam.power());
            close(actual.pdf, expected.pdf, 0.0);
            close(
                actual
                    .beam
                    .frame()
                    .direction()
                    .dot(expected.beam.frame().direction()),
                1.0,
                1.0,
            );
            let frame = PolarizationFrame::try_new(
                actual.beam.frame().direction().scale(-1.0),
                Vec3::new(0.0, 1.0, 0.0),
            )
            .unwrap();
            let state = Stokes {
                i: 1.0,
                q: -0.1,
                u: 0.4,
                v: -0.7,
            };
            let a = input
                .evaluate_stack_dielectric_connection(n, glass(eta_t), rough(), stack, frame, state)
                .unwrap();
            let b = input
                .evaluate_dielectric_connection(n, glass(eta_t), compact, frame, state)
                .unwrap();
            close(a.factor, b.factor, 0.0);
            close(a.pdf, b.pdf, 0.0);
        }
    }
}

#[test]
fn isotropic_canonicalization_and_invalid_axis_keep_existing_measures() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let input = beam(Vec3::new(0.2, 0.1, -1.0), 1.0, Stokes::unpolarized(1.0));
    let iso = DielectricSurface::try_rough(0.25).unwrap();
    let equal =
        DielectricSurface::try_anisotropic_rough(0.25, 0.25, Vec3::new(0.3, 0.5, 0.8)).unwrap();
    for draw in [[0.04, 0.3, 0.0], [0.04, 0.3, 0.999], [0.99, 0.2, 0.0]] {
        assert_eq!(
            input.sample_interface(n, glass(1.5), iso, draw),
            input.sample_interface(n, glass(1.5), equal, draw)
        );
    }
    let invalid = DielectricSurface::try_anisotropic_rough(0.1, 0.7, n).unwrap();
    assert!(
        input
            .sample_interface(n, glass(1.5), invalid, [0.1, 0.2, 0.3])
            .is_err()
    );
    let frame =
        PolarizationFrame::try_new(Vec3::new(0.1, 0.0, -1.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
    assert!(
        input
            .evaluate_dielectric_connection(n, glass(1.5), invalid, frame, Stokes::unpolarized(1.0))
            .is_err()
    );
    let null = input
        .sample_interface(n, glass(1.0), rough(), [0.1, 0.2, 0.3])
        .unwrap()
        .unwrap();
    assert!(null.delta);
    assert_eq!(null.beam, input);
    assert_eq!(null.pdf, 1.0);
}
