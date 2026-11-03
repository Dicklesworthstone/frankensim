//! Optional native boundary-area proposals, in the chart's original frame.
//!
//! A chart opts in through all three [`crate::Chart`] methods: total area,
//! sampling and density evaluation. They describe the same immutable boundary
//! and area measure. Sampling need not be uniform, but the returned positive
//! finite density must be the actual proposal and must agree at native hits.
//! These numerical values carry no geometry conversion, visibility, material,
//! intersection, topology or certified differential authority. Consumers retain
//! the original chart owner and perform their native contact/departure checks.
//!
//! Each operation takes an execution context through the chart trait. Constant
//! work implementations checkpoint once; iterative implementations must retain
//! bounded cancellation and refuse exhausted work rather than return a partial
//! or success-conditioned proposal. Floating-point arithmetic is numerical,
//! with no cross-platform transcendental bit-stability claim.

use crate::{Point3, Vec3};
use fs_exec::Cancelled;

/// One sample of the original boundary, before any instance placement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceAreaSample {
    /// Finite position in the original chart frame, in metres.
    pub position_m: Point3,
    /// Finite unit outward normal of that boundary.
    pub normal: Vec3,
    /// Positive finite proposal density with respect to physical area.
    pub pdf_area_per_m2: f64,
}

/// Matching proposal evaluated at a caller's original native contact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceAreaDensity {
    /// Finite unit outward normal at the represented boundary location.
    pub normal: Vec3,
    /// Positive finite density of the same sampling measure.
    pub pdf_area_per_m2: f64,
}

/// Explicit refusal of native area sampling or density evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceAreaError {
    /// This exact chart implementation has not opted into the capability.
    Unsupported,
    /// The chart's defining geometry is malformed or degenerate.
    InvalidGeometry,
    /// A variate, query point or distance allowance is invalid.
    InvalidSample,
    /// The query is outside the caller's boundary-distance allowance.
    OffSurface,
    /// The query cannot select one regular oriented boundary location.
    Ambiguous,
    /// Required positive finite geometry, position or density is unrepresentable.
    NumericalRange,
    /// The operation observed cancellation before publishing a result.
    Cancelled(Cancelled),
}

impl From<Cancelled> for SurfaceAreaError {
    fn from(error: Cancelled) -> Self {
        Self::Cancelled(error)
    }
}

impl core::fmt::Display for SurfaceAreaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native surface-area proposal: {self:?}")
    }
}

impl core::error::Error for SurfaceAreaError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{BoxChart, SphereChart, TorusChart};
    use crate::{Aabb, Chart, ChartSample, TraceStepClaim};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

    fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
        let gate = CancelGate::new_clock_free();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(&gate, arena,
                StreamKey { seed: 71, kernel_id: 8, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            f(&cx, &gate);
        });
    }

    fn near(actual: f64, expected: f64, relative: f64) {
        assert!((actual - expected).abs() <= relative * expected.abs(),
            "{actual:.17e} != {expected:.17e}");
    }

    #[test]
    fn g1_sphere_area_sampler_has_uniform_moments_and_matching_density() {
        with_cx(|cx, _| {
            let sphere = SphereChart { center: Point3::new(3.0, -4.0, 7.0), radius: 2.5 };
            let chart: &dyn Chart = &sphere;
            let area = 25.0 * core::f64::consts::PI;
            near(chart.surface_area(cx).unwrap(), area, 2e-15);
            let mut first = [0.0_f64; 3];
            let mut second = [0.0_f64; 3];
            let count = 64 * 32;
            for i in 0..64 {
                for j in 0..32 {
                    let sample = chart.sample_surface_area(cx,
                        [(i as f64 + 0.5) / 64.0, (j as f64 + 0.5) / 32.0]).unwrap();
                    let density = chart.surface_area_pdf_at(cx, sample.position_m, 1e-13).unwrap();
                    near(sample.pdf_area_per_m2 * area, 1.0, 2e-15);
                    assert_eq!(density.pdf_area_per_m2, sample.pdf_area_per_m2);
                    near(sample.normal.dot(sample.normal), 1.0, 1e-14);
                    assert!(density.normal.dot(sample.normal) > 1.0 - 1e-14);
                    let radial = sample.position_m.delta_from(sphere.center).scale(1.0 / sphere.radius);
                    let coordinates = [radial.x, radial.y, radial.z];
                    for k in 0..3 {
                        first[k] += coordinates[k];
                        second[k] += coordinates[k] * coordinates[k];
                    }
                }
            }
            // Uniform physical area has E[n_i]=0 and E[n_i²]=1/3. This
            // independent moment check rejects uniform-latitude sampling.
            for k in 0..3 {
                assert!((first[k] / count as f64).abs() < 1e-14);
                assert!((second[k] / count as f64 - 1.0 / 3.0).abs() < 1e-4);
            }
            let pole = chart.sample_surface_area(cx, [0.0, 0.73]).unwrap();
            assert_eq!(pole.position_m, Point3::new(3.0, -4.0, 9.5));
            assert_eq!(pole.normal, Vec3::new(0.0, 0.0, 1.0));
            assert_eq!(chart.surface_area_pdf_at(cx, pole.position_m, 0.0).unwrap().normal,
                pole.normal);
        });
    }

    #[test]
    fn g3_sphere_area_density_scales_and_translates_without_overflowing_norms() {
        with_cx(|cx, _| {
            for radius in [1e-150, 0.75, 1e150] {
                let sphere = SphereChart {
                    center: Point3::new(3.0 * radius, -4.0 * radius, 7.0 * radius), radius,
                };
                let sample = sphere.sample_surface_area(cx, [0.173, 0.731]).unwrap();
                let area = sphere.surface_area(cx).unwrap();
                near(area / radius / radius, 4.0 * core::f64::consts::PI, 3e-15);
                near(sample.pdf_area_per_m2 * area, 1.0, 3e-15);
                let tolerance = 512.0 * f64::EPSILON * radius;
                let density = sphere.surface_area_pdf_at(cx, sample.position_m, tolerance).unwrap();
                assert_eq!(sample.pdf_area_per_m2, density.pdf_area_per_m2);
                assert!(sample.normal.dot(density.normal) > 1.0 - 1e-14);
                let nearby = sample.position_m.offset(sample.normal.scale(0.25 * tolerance));
                let at_contact = sphere.surface_area_pdf_at(cx, nearby, tolerance).unwrap();
                assert_eq!(at_contact.pdf_area_per_m2, sample.pdf_area_per_m2);
                assert!(at_contact.normal.dot(sample.normal) > 1.0 - 1e-14);
            }
        });
    }

    #[test]
    fn g0_sphere_area_refuses_malformed_geometry_and_unresolved_queries() {
        with_cx(|cx, _| {
            let sphere = SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 1.0 };
            for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert_eq!(SphereChart { radius, ..sphere }.surface_area(cx),
                    Err(SurfaceAreaError::InvalidGeometry));
            }
            assert_eq!(SphereChart { center: Point3::new(f64::NAN, 0.0, 0.0), ..sphere }
                .surface_area(cx), Err(SurfaceAreaError::InvalidGeometry));
            for radius in [1e-160, 1e154] {
                assert_eq!(SphereChart { radius, ..sphere }.surface_area(cx),
                    Err(SurfaceAreaError::NumericalRange));
            }
            assert_eq!(SphereChart { center: Point3::new(1e300, 0.0, 0.0), ..sphere }
                .surface_area(cx), Err(SurfaceAreaError::NumericalRange));
            let coarse = SphereChart { center: Point3::new(1e16, 1e16, 1e16), radius: 2.0 };
            assert!(coarse.surface_area(cx).is_ok());
            assert_eq!(coarse.sample_surface_area(cx, [0.25, 0.125]),
                Err(SurfaceAreaError::NumericalRange));
            for samples in [[1.0, 0.5], [-0.1, 0.5], [0.5, f64::NAN], [0.5, f64::INFINITY]] {
                assert_eq!(sphere.sample_surface_area(cx, samples), Err(SurfaceAreaError::InvalidSample));
            }
            assert_eq!(sphere.surface_area_pdf_at(cx, Point3::new(2.0, 0.0, 0.0), 0.1),
                Err(SurfaceAreaError::OffSurface));
            assert_eq!(sphere.surface_area_pdf_at(cx, sphere.center, 0.1),
                Err(SurfaceAreaError::OffSurface));
            assert_eq!(sphere.surface_area_pdf_at(cx, sphere.center, 1.0),
                Err(SurfaceAreaError::Ambiguous));
            for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
                assert_eq!(sphere.surface_area_pdf_at(cx, Point3::new(1.0, 0.0, 0.0), tolerance),
                    Err(SurfaceAreaError::InvalidSample));
            }
            assert_eq!(sphere.surface_area_pdf_at(cx, Point3::new(f64::NAN, 0.0, 0.0), 0.1),
                Err(SurfaceAreaError::InvalidSample));
        });
    }

    struct UnoptedSphere(SphereChart);

    impl Chart for UnoptedSphere {
        fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample { self.0.eval(point, cx) }
        fn support(&self) -> Aabb { self.0.support() }
        fn trace_step_claim(&self) -> TraceStepClaim { TraceStepClaim::ExactDistance }
        fn name(&self) -> &'static str { "fixture/sphere" }
    }

    #[test]
    fn g0_area_capability_is_not_inferred_from_names_or_distance_claims() {
        with_cx(|cx, _| {
            let unopted = UnoptedSphere(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0), radius: 1.0,
            });
            let box_chart = BoxChart {
                aabb: Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0)),
            };
            let torus = TorusChart { center: Point3::new(0.0, 0.0, 0.0), major: 2.0, minor: 0.5 };
            for chart in [&unopted as &dyn Chart, &box_chart, &torus] {
                assert_eq!(chart.surface_area(cx), Err(SurfaceAreaError::Unsupported));
                assert_eq!(chart.sample_surface_area(cx, [0.3, 0.7]), Err(SurfaceAreaError::Unsupported));
                assert_eq!(chart.surface_area_pdf_at(cx, Point3::new(1.0, 0.0, 0.0), 0.1),
                    Err(SurfaceAreaError::Unsupported));
            }
        });
    }

    #[test]
    fn g4_area_operations_preserve_cancellation_including_unsupported_charts() {
        with_cx(|cx, gate| {
            let sphere = SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 1.0 };
            let unopted = UnoptedSphere(sphere);
            gate.request();
            for chart in [&sphere as &dyn Chart, &unopted] {
                assert_eq!(chart.surface_area(cx), Err(SurfaceAreaError::Cancelled(Cancelled)));
                assert_eq!(chart.sample_surface_area(cx, [0.3, 0.7]),
                    Err(SurfaceAreaError::Cancelled(Cancelled)));
                assert_eq!(chart.surface_area_pdf_at(cx, Point3::new(1.0, 0.0, 0.0), 0.1),
                    Err(SurfaceAreaError::Cancelled(Cancelled)));
            }
        });
    }
}
