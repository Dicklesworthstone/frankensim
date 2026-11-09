//! G1/G3: independent static Green function, physical equations and limits.
use super::*;
use crate::diffraction::{HalfPlaneDiffraction, HalfPlaneFace};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 340.0).unwrap()
}
fn frame() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}
fn limits(k: f64) -> HalfPlaneMonopoleLimits {
    HalfPlaneMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: k * 340.0 * 1.3 * 1e-12,
        absolute_velocity_tolerance_per_m2: 1e-11,
        relative_tolerance: 1e-10,
        maximum_evaluations: 100_000,
        maximum_depth: 42,
    }
}
fn prepare(
    source: [f64; 3],
    receiver: HalfPlaneReceiver,
    k: f64,
    boundary: IdealPlaneBoundary,
) -> PreparedHalfPlaneMonopole {
    PreparedHalfPlaneMonopole::new(
        medium(),
        k * 340.0 / TAU,
        frame(),
        boundary,
        source,
        receiver,
        limits(k),
        || false,
    )
    .unwrap_or_else(|error| panic!("source={source:?}, receiver={receiver:?}, k={k}: {error:?}"))
}
fn green(
    source: [f64; 3],
    point: [f64; 3],
    k: f64,
    boundary: IdealPlaneBoundary,
) -> (C64, [C64; 3]) {
    let owner = prepare(source, HalfPlaneReceiver::WorldPosition(point), k, boundary);
    let line = owner.evaluate_harmonic_line(C64::ONE).unwrap();
    (
        line.pressure_pa / C64::new(0.0, -k * 340.0 * 1.3),
        line.particle_velocity_m_s.map(|u| -u),
    )
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn static_half_plane_and_open_plane_source_have_independent_exact_limits() {
    let source = [-0.4_f64, 0.3, 0.2];
    let point = [0.3_f64, -0.2, 0.1];
    let rs = source[0].hypot(source[1]);
    let r = point[0].hypot(point[1]);
    let theta = point[1].atan2(point[0]).rem_euclid(TAU);
    let theta_s = source[1].atan2(source[0]).rem_euclid(TAU);
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let mut reference = 0.0;
        for (i, angle) in [theta - theta_s, theta + theta_s].into_iter().enumerate() {
            let q = 2.0 * (r * rs).sqrt() * (0.5 * angle).cos();
            let d = ((r - rs).powi(2)
                + 4.0 * r * rs * (0.5 * angle).sin().powi(2)
                + (point[2] - source[2]).powi(2))
            .sqrt();
            reference += (if i == 0 { 1.0 } else { boundary.sign() }) * (0.5 + (q / d).atan() / PI)
                / (4.0 * PI * d);
        }
        let (g, _) = green(source, point, 1e-7, boundary);
        assert!((g.re - reference).abs() < 2e-11);
        assert!(g.im.abs() < 2e-8);
    }
    // A rigid sheet cannot scatter a monopole on its open continuation:
    // the original spherical field already has zero normal velocity there.
    let source = [-0.4, 0.0, 0.2];
    for k in [0.03, 7.0, 120.0] {
        let (actual, gradient) = green(source, point, k, IdealPlaneBoundary::Rigid);
        let delta: [f64; 3] = core::array::from_fn(|i| point[i] - source[i]);
        let d = delta[0].hypot(delta[1]).hypot(delta[2]);
        let expected = C64::new((k * d).cos(), (k * d).sin()).scale(1.0 / (4.0 * PI * d));
        close(actual, expected, 3e-10);
        for i in 0..3 {
            close(
                gradient[i],
                expected * C64::new(-1.0 / d, k).scale(delta[i] / d),
                3e-10,
            );
        }
    }
}

#[test]
fn analytic_velocity_and_original_helmholtz_equation_hold_in_light_and_shadow() {
    let source = [-0.4, 0.3, 0.2];
    let k = 7.0;
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for point in [[0.3, 0.2, 0.1], [0.4, -0.2, -0.3], [-0.15, -0.25, 0.2]] {
            let (g, gradient) = green(source, point, k, boundary);
            let h = 1e-4;
            let mut laplacian = C64::ZERO;
            for axis in 0..3 {
                let sample = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
                    let mut at = point;
                    at[axis] += offset * h;
                    green(source, at, k, boundary).0
                });
                let derivative = (sample[0] - sample[1].scale(8.0) + sample[2].scale(8.0)
                    - sample[3])
                    .scale(1.0 / (12.0 * h));
                assert!((derivative - gradient[axis]).abs() < 3e-8);
                laplacian = laplacian
                    + (-sample[0] + sample[1].scale(16.0) - g.scale(30.0) + sample[2].scale(16.0)
                        - sample[3])
                        .scale(1.0 / (12.0 * h * h));
            }
            assert!((laplacian + g.scale(k * k)).abs() < 1e-5);
        }
    }
}

#[test]
fn reciprocity_both_sheet_faces_and_occluded_image_coincidence_remain_physical() {
    let source = [0.4, 0.3, 0.2];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for point in [[0.3, 0.2, 0.1], [0.4, -0.3, 0.2], [-0.15, -0.25, 0.2]] {
            let (forward, _) = green(source, point, 7.0, boundary);
            let (reverse, _) = green(point, source, 7.0, boundary);
            close(forward, reverse, 2e-9);
        }
        for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
            let owner = prepare(
                source,
                HalfPlaneReceiver::SheetFace {
                    distance_from_edge_m: 0.3,
                    along_edge_m: 0.1,
                    face,
                },
                7.0,
                boundary,
            );
            let line = owner.evaluate_harmonic_line(C64::ONE).unwrap();
            match boundary {
                IdealPlaneBoundary::Rigid => assert_eq!(line.particle_velocity_m_s[1], C64::ZERO),
                IdealPlaneBoundary::PressureRelease => assert_eq!(line.pressure_pa, C64::ZERO),
            }
        }
    }
}

#[test]
fn geometrical_shadow_boundary_has_a_finite_value_and_analytic_gradient() {
    let source = [0.4, 0.3, 0.2];
    let boundary_point = [-0.4, -0.3, 0.5];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let (center, gradient) = green(source, boundary_point, 7.0, boundary);
        for h in [1e-4, 1e-6, 1e-8] {
            let mut lo = boundary_point;
            lo[0] -= h;
            let mut hi = boundary_point;
            hi[0] += h;
            let left = green(source, lo, 7.0, boundary).0;
            let right = green(source, hi, 7.0, boundary).0;
            assert!((left - center).abs() < 2.0 * h && (right - center).abs() < 2.0 * h);
            assert!(((right - left).scale(0.5 / h) - gradient[0]).abs() < 3e-7);
        }
    }
}

#[test]
fn finite_source_converges_to_the_original_sommerfeld_plane_wave() {
    let direction = [-0.48, -0.64, 0.6];
    let point = [0.3, -0.2, 0.1];
    let k = 7.0;
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let reference = HalfPlaneDiffraction::new(medium(), k * 340.0 / TAU, direction, boundary)
            .unwrap()
            .evaluate(point, C64::ONE, || false)
            .unwrap();
        let mut errors = Vec::new();
        for distance in [100.0, 1000.0, 10000.0] {
            let source = direction.map(|d| -distance * d);
            let owner = prepare(source, HalfPlaneReceiver::WorldPosition(point), k, boundary);
            // Q chosen so the incident spherical pressure at the edge origin
            // is exactly one; no plane-wave transfer is used in preparation.
            let clock = C64::new((-k * distance).cos(), (-k * distance).sin());
            let flow = clock.scale(4.0 * PI * distance) / C64::new(0.0, -k * 340.0 * 1.3);
            let field = owner.evaluate(flow).unwrap();
            let error = (field.pressure_pa - reference.pressure_pa).abs();
            errors.push(error);
        }
        assert!(
            errors[1] < 0.11 * errors[0] && errors[2] < 0.11 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1e-5, "{errors:?}");
    }
}

#[test]
fn placed_source_velocity_and_source_phase_use_one_original_world_frame() {
    let source = [-0.4, 0.3, 0.2];
    let point = [0.3, -0.2, 0.1];
    let origin = [1.25, -2.0, 0.5];
    let world = |p: [f64; 3]| [origin[0] + p[2], origin[1] + p[0], origin[2] + p[1]];
    let placed = HalfPlaneFrame::new(origin, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let expected = prepare(
        source,
        HalfPlaneReceiver::WorldPosition(point),
        7.0,
        IdealPlaneBoundary::Rigid,
    );
    let actual = PreparedHalfPlaneMonopole::new(
        medium(),
        7.0 * 340.0 / TAU,
        placed,
        IdealPlaneBoundary::Rigid,
        world(source),
        HalfPlaneReceiver::WorldPosition(world(point)),
        limits(7.0),
        || false,
    )
    .unwrap();
    let flow = C64::new(0.3, -0.7);
    let a = actual.evaluate(flow).unwrap();
    let b = expected.evaluate(flow).unwrap();
    close(a.pressure_pa, b.pressure_pa, 3e-10);
    for (i, j) in [(0, 2), (1, 0), (2, 1)] {
        close(
            a.particle_velocity_m_s[i],
            b.particle_velocity_m_s[j],
            3e-10,
        );
    }
    let clock = C64::new(0.41_f64.cos(), 0.41_f64.sin());
    let phased = actual.evaluate(flow * clock).unwrap();
    close(phased.pressure_pa, a.pressure_pa * clock, 3e-15);
    for i in 0..3 {
        assert!((phased.active_intensity_w_m2[i] - a.active_intensity_w_m2[i]).abs() < 1e-11);
    }
    assert!(actual.diagnostics().evaluations > 0);
    assert_eq!(actual.source_position_m(), world(source));
    assert_eq!(actual.frame(), placed);
    assert_eq!(actual.boundary(), IdealPlaneBoundary::Rigid);
    assert_eq!(actual.acoustic(), medium());
}

#[test]
fn invalid_geometry_cancelled_or_unresolved_work_publishes_no_transfer() {
    let source = [0.4, 0.3, 0.2];
    let point = HalfPlaneReceiver::WorldPosition([0.3, -0.2, 0.1]);
    let build = |s, receiver, controls, cancelled| {
        PreparedHalfPlaneMonopole::new(
            medium(),
            7.0 * 340.0 / TAU,
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            receiver,
            controls,
            || cancelled,
        )
    };
    assert!(
        build(
            source,
            HalfPlaneReceiver::WorldPosition(source),
            limits(7.0),
            false
        )
        .is_err()
    );
    for s in [[0.0, 0.0, 0.0], [0.4, 0.0, 0.2], [f64::NAN, 0.2, 0.1]] {
        assert!(build(s, point, limits(7.0), false).is_err());
    }
    assert!(
        build(
            source,
            HalfPlaneReceiver::WorldPosition([0.4, 0.0, 0.2]),
            limits(7.0),
            false
        )
        .is_err()
    );
    assert!(matches!(
        build(source, point, limits(7.0), true),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut limited = limits(7.0);
    limited.maximum_evaluations = 1;
    assert!(build(source, point, limited, false).is_err());
    let mut polls = 0;
    assert!(matches!(
        PreparedHalfPlaneMonopole::new(
            medium(),
            7.0 * 340.0 / TAU,
            frame(),
            IdealPlaneBoundary::Rigid,
            source,
            point,
            limits(7.0),
            || {
                polls += 1;
                polls == 20
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let owner = build(source, point, limits(7.0), false).unwrap();
    assert!(owner.evaluate(C64::new(f64::NAN, 0.0)).is_err());
    assert_eq!(
        owner.evaluate_harmonic_line(C64::ZERO).unwrap().pressure_pa,
        C64::ZERO
    );
}

#[test]
fn physical_transfer_and_tail_resolution_survive_raw_contour_decay_underflow() {
    let omega_rho = 1e300;
    let k = 7.0;
    let mut kernel = PointKernel::new(
        [-0.4, 0.3, 0.2],
        [0.3, -0.2, 0.1],
        k,
        omega_rho,
        IdealPlaneBoundary::Rigid,
        frame().axes(),
        QuadratureLimits::from(limits(k)),
    )
    .unwrap();
    // Isolate the contour correction from the geometrical field; otherwise
    // adding a much larger geometrical term hides this finite transfer.
    kernel.base = [C64::ZERO; 4];
    kernel.absolute[0] = 1e-200;
    let region = kernel.count - 1;
    let delta = 220.0 - kernel.breaks[region];
    let x = delta / (delta + kernel.scale);
    let (w, jacobian) = kernel.coordinate(region, x);
    let s = C64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2).scale(w * kernel.distance);
    let radius = (C64::from_re(kernel.distance * kernel.distance) + s * s).sqrt();
    let attenuation = k * radius.im;
    assert_eq!(det::exp(-attenuation), 0.0);
    assert!(!kernel.tail_is_small(w));
    // Independent dimensional f'(s) ds, with two balanced half-decays.
    // Neither the reference nor production materializes the vanishing full
    // exponential before its original physical source coefficient.
    let half_decay = det::exp(-0.5 * attenuation);
    let exponential = cis(Dd::from_f64(k) * Dd::from_f64(radius.re)).scale(half_decay);
    let derivative =
        s * exponential * (C64::new(0.0, k) * radius - C64::ONE) / (radius * radius * radius);
    let angular = contour_arctangent(kernel.eta[0], w) + contour_arctangent(kernel.eta[1], w);
    let measure = C64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2).scale(kernel.distance * jacobian);
    let integral = derivative * angular * measure;
    let expected =
        C64::new(-integral.im, integral.re).scale((omega_rho * half_decay) / (4.0 * PI * PI));
    let actual = kernel.sample(region, x).unwrap()[0];
    assert!(expected.abs() > 0.0 && actual.abs() > 0.0);
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs(),
        "{actual:?} != balanced dimensional reference {expected:?}"
    );
}
