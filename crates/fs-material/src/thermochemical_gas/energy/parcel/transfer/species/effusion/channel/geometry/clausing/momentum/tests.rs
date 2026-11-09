//! Independent ray tracing and aperture-overlap quadrature for angular moments.
use super::*;

fn geometry(h: f64) -> CylindricalMolecularGeometry {
    CylindricalMolecularGeometry::try_new(
        1.0,
        h,
        100.0 * h.max(1.0),
        "synthetic diffuse circular channel",
    )
    .unwrap()
}
fn controls() -> ClausingControls {
    ClausingControls {
        initial_rings: 8,
        maximum_rings: 2048,
        maximum_iterations: 512,
        absolute_tolerance: 2e-7,
        relative_tolerance: 5e-5,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn aperture_limit_is_an_angular_moment_and_not_a_mass_probability() {
    let (g, moments) = geometry(0.0)
        .resolve_axial_momentum(controls(), || false)
        .unwrap();
    assert_eq!(g.transmission_probability(), 1.0);
    assert_eq!(moments.direct_axial_cosine(), 2.0 / 3.0);
    assert_eq!(moments.wall_transmitted_axial_cosine(), 0.0);
    assert_eq!(moments.reflected_axial_cosine(), 0.0);
    assert_eq!(moments.resolution().rings, 0);
    assert_eq!(moments.resolution().quadrature_evaluations, 0);
    assert_eq!(moments.resolution().direct_error_estimate, 0.0);
}

#[test]
fn shared_mass_and_moment_grids_obey_exit_reciprocity_and_joint_budgets() {
    for h in [1e-8, 0.01, 0.5, 2.0, 10.0] {
        let original = geometry(h);
        let mass_only = original
            .clone()
            .resolve_transmission(controls(), || false)
            .unwrap();
        let (g, m) = original
            .resolve_axial_momentum(controls(), || false)
            .unwrap();
        let report = m.resolution();
        let mass = g.resolution().unwrap();
        assert_eq!(report.rings, mass.rings);
        assert_eq!(report.iterations, mass.iterations);
        assert!(report.rings >= mass_only.resolution().unwrap().rings);
        assert!(report.angular_panels <= controls().maximum_rings / 2);
        assert!(report.refinements >= 2);
        for (moment, probability, error) in [
            (
                m.direct_axial_cosine(),
                g.direct_probability(),
                report.direct_error_estimate,
            ),
            (
                m.wall_transmitted_axial_cosine(),
                g.wall_transmitted_probability(),
                report.wall_transmitted_error_estimate,
            ),
            (
                m.reflected_axial_cosine(),
                g.reflected_probability(),
                report.reflected_error_estimate,
            ),
        ] {
            assert!(moment > 0.0 && moment <= probability);
            assert!(
                error <= controls().absolute_tolerance + controls().relative_tolerance * moment
            );
        }
        near(
            m.transmitted_axial_cosine() + m.reflected_axial_cosine(),
            2.0 / 3.0,
            2e-12,
        );
        assert!(m.transmitted_axial_cosine() / g.transmission_probability() > 2.0 / 3.0);
        near(
            g.transmission_probability(),
            mass_only.transmission_probability(),
            mass.estimated_absolute_error
                + mass_only.resolution().unwrap().estimated_absolute_error,
        );
    }
}

// Independent direct-flight formulation: integrate the overlap of entrance and
// exit discs shifted by L*tan(theta). Cosine-law incidence has density 2*mu,
// so the axial moment uses 2*mu^2 times the geometric overlap fraction.
fn overlap_direct_moment(h: f64) -> f64 {
    let lower = h / (1.0 + h * h).sqrt();
    let n = 16384;
    let delta = (1.0 - lower) / n as f64;
    let value = |i: usize| {
        let mu = lower + i as f64 * delta;
        let shift = (h * (1.0 - mu * mu).max(0.0).sqrt() / mu).min(1.0);
        let overlap = 2.0 / core::f64::consts::PI
            * (shift.acos() - shift * (1.0 - shift * shift).sqrt());
        2.0 * mu * mu * overlap
    };
    let interior = (1..n).map(|i| {
        if i % 2 == 0 {
            2.0 * value(i)
        } else {
            4.0 * value(i)
        }
    });
    delta / 3.0 * (value(0) + value(n) + interior.sum::<f64>())
}

#[test]
fn direct_angular_integral_matches_independent_shifted_disc_overlap() {
    for h in [0.125, 1.0, 8.0] {
        let (_, moments) = geometry(h)
            .resolve_axial_momentum(controls(), || false)
            .unwrap();
        near(
            moments.direct_axial_cosine(),
            overlap_direct_moment(h),
            3e-8,
        );
    }
}

// Fixed-seed test-only ray tracing, independent of the Clausing kernel and all
// production angular integrals. Samples uniform entry area, cosine incidence
// and fresh diffuse wall directions; counts actual first exit direction cosines.
fn ray_trace(h: f64, histories: usize) -> [f64; 3] {
    let mut seed = 0x4d595df4d0f33173_u64;
    let mut uniform = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        ((seed >> 32) as f64 + 0.5) / 4294967296.0
    };
    let mut totals = [0.0; 3];
    for _ in 0..histories {
        let radius = 0.5 * uniform().sqrt();
        let azimuth = 2.0 * core::f64::consts::PI * uniform();
        let (mut x, mut y, mut z) = (radius * azimuth.cos(), radius * azimuth.sin(), 0.0);
        let mut vz = uniform().sqrt();
        let transverse = (1.0 - vz * vz).sqrt();
        let azimuth = 2.0 * core::f64::consts::PI * uniform();
        let (mut vx, mut vy) = (transverse * azimuth.cos(), transverse * azimuth.sin());
        let mut wall_hit = false;
        let mut completed = false;
        for _ in 0..10000 {
            let transverse2 = vx * vx + vy * vy;
            let b = x * vx + y * vy;
            let wall_distance = if wall_hit {
                -2.0 * b / transverse2
            } else {
                let gap = 0.25 - x * x - y * y;
                let root = (b * b + transverse2 * gap).sqrt();
                if b > 0.0 {
                    gap / (b + root)
                } else {
                    (root - b) / transverse2
                }
            };
            let exit_distance = if vz > 0.0 { (h - z) / vz } else { -z / vz };
            if exit_distance < wall_distance {
                let index = if vz < 0.0 {
                    2
                } else if wall_hit {
                    1
                } else {
                    0
                };
                totals[index] += vz.abs();
                completed = true;
                break;
            }
            x += wall_distance * vx;
            y += wall_distance * vy;
            z += wall_distance * vz;
            let radius = (x * x + y * y).sqrt();
            x *= 0.5 / radius;
            y *= 0.5 / radius;
            let normal = uniform().sqrt();
            let transverse = (1.0 - normal * normal).sqrt();
            let azimuth = 2.0 * core::f64::consts::PI * uniform();
            let tangent = transverse * azimuth.cos();
            vz = transverse * azimuth.sin();
            vx = -2.0 * (normal * x + tangent * y);
            vy = 2.0 * (tangent * x - normal * y);
            wall_hit = true;
        }
        assert!(
            completed,
            "independent ray history exceeded its explicit budget"
        );
    }
    totals.map(|sum| sum / histories as f64)
}

#[test]
fn reflected_and_transmitted_moments_match_independent_diffuse_ray_histories() {
    let h = 2.0;
    let (_, moments) = geometry(h)
        .resolve_axial_momentum(controls(), || false)
        .unwrap();
    let reference = ray_trace(h, 160_000);
    for (actual, ray) in [
        moments.direct_axial_cosine(),
        moments.wall_transmitted_axial_cosine(),
        moments.reflected_axial_cosine(),
    ]
    .into_iter()
    .zip(reference)
    {
        near(actual, ray, 0.0035);
    }
}

#[test]
fn dimensionless_results_scale_deterministically_and_retain_tiny_wall_paths() {
    let reference = geometry(2.0)
        .resolve_axial_momentum(controls(), || false)
        .unwrap()
        .1;
    for scale in [1e-90, 1.0, 1e90] {
        let g = CylindricalMolecularGeometry::try_new(
            scale,
            2.0 * scale,
            200.0 * scale,
            "scaled synthetic diffuse tube",
        )
        .unwrap();
        let first = g
            .clone()
            .resolve_axial_momentum(controls(), || false)
            .unwrap()
            .1;
        let repeated = g
            .resolve_axial_momentum(controls(), || false)
            .unwrap()
            .1;
        assert_eq!(first, repeated);
        assert_eq!(first, reference);
    }
    let (g, tiny) = geometry(1e-12)
        .resolve_axial_momentum(controls(), || false)
        .unwrap();
    assert!(tiny.wall_transmitted_axial_cosine() > 0.0);
    assert!(tiny.reflected_axial_cosine() > 0.0);
    near(
        tiny.reflected_axial_cosine() / g.reflected_probability(),
        4.0 / (3.0 * core::f64::consts::PI),
        2e-10,
    );
}

#[test]
fn invalid_controls_refinement_exhaustion_and_cancellation_refuse_atomically() {
    let original = geometry(2.0);
    let mut invalid_controls = controls();
    invalid_controls.initial_rings = 3;
    assert!(
        original
            .clone()
            .resolve_axial_momentum(invalid_controls, || false)
            .is_err()
    );
    let mut exhausted = controls();
    exhausted.initial_rings = 4;
    exhausted.maximum_rings = 16;
    exhausted.absolute_tolerance = 1e-15;
    exhausted.relative_tolerance = 0.0;
    assert!(
        original
            .clone()
            .resolve_axial_momentum(exhausted, || false)
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        original
            .clone()
            .resolve_axial_momentum(controls(), || {
                polls += 1;
                polls > 70
            })
            .unwrap_err(),
        GasEnergyError::Cancelled
    );
    assert_eq!(polls, 71);
    assert_eq!(original, geometry(2.0));
    assert!(
        original
            .resolve_axial_momentum(controls(), || false)
            .is_ok()
    );
}
