//! G1/G3: physical force units, analytic mixed derivatives, source/receiver
//! reciprocity, exact sheet traces and geometrical-boundary continuity.
use super::*;
use crate::diffraction::HalfPlaneFace;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 340.0).unwrap()
}
fn frame() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}
fn limits() -> HalfPlaneForceLimits {
    HalfPlaneForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-11,
        absolute_velocity_tolerance_m_s_per_n: 1e-13,
        relative_tolerance: 1e-10,
        maximum_evaluations: 100_000,
        maximum_depth: 42,
    }
}
fn prepare(s: [f64; 3], x: [f64; 3], k: f64, b: IdealPlaneBoundary) -> PreparedHalfPlaneForce {
    PreparedHalfPlaneForce::new(
        medium(),
        k * 340.0 / TAU,
        frame(),
        b,
        s,
        HalfPlaneReceiver::WorldPosition(x),
        limits(),
        || false,
    )
    .unwrap_or_else(|e| panic!("{s:?} -> {x:?}, k={k}, {b:?}: {e:?}"))
}
fn mono(s: [f64; 3], x: [f64; 3], k: f64, b: IdealPlaneBoundary) -> HarmonicFieldLine<3> {
    PreparedHalfPlaneMonopole::new(
        medium(),
        k * 340.0 / TAU,
        frame(),
        b,
        s,
        HalfPlaneReceiver::WorldPosition(x),
        HalfPlaneMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: k * 340.0 * 1.3 * 1e-12,
            absolute_velocity_tolerance_per_m2: 1e-11,
            relative_tolerance: 1e-11,
            maximum_evaluations: 100_000,
            maximum_depth: 42,
        },
        || false,
    )
    .unwrap()
    .evaluate_harmonic_line(C64::ONE)
    .unwrap()
}
fn values(line: HarmonicFieldLine<3>) -> [C64; 4] {
    [
        line.pressure_pa,
        line.particle_velocity_m_s[0],
        line.particle_velocity_m_s[1],
        line.particle_velocity_m_s[2],
    ]
}
fn close(actual: C64, expected: C64, scale: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * scale.max(expected.abs()).max(1e-30),
        "{actual:?} != {expected:?}, scaled error={}",
        (actual - expected).abs() / scale.max(expected.abs()).max(1e-30)
    );
}
fn original_force(s: [f64; 3], x: [f64; 3], k: f64, f: [C64; 3]) -> [C64; 4] {
    let delta: [f64; 3] = core::array::from_fn(|i| x[i] - s[i]);
    let d = delta.into_iter().fold(0.0_f64, f64::hypot);
    let n = delta.map(|v| v / d);
    let g = C64::new((k * d).cos(), (k * d).sin()).scale(1.0 / (4.0 * PI * d));
    let first = g * C64::new(-1.0 / d, k);
    let second = g * C64::new(2.0 / d / d - k * k, -2.0 * k / d);
    let mut p = C64::ZERO;
    let mut u = [C64::ZERO; 3];
    for j in 0..3 {
        p = p - first.scale(n[j]) * f[j];
        for i in 0..3 {
            let h = second.scale(n[i] * n[j])
                + first.scale(((if i == j { 1.0 } else { 0.0 }) - n[i] * n[j]) / d);
            u[i] = u[i] + C64::new(-h.im, h.re).scale(1.0 / (k * 340.0 * 1.3)) * f[j];
        }
    }
    [p, u[0], u[1], u[2]]
}

#[test]
fn physical_force_recovers_independent_free_spherical_dipole_fields() {
    // These source symmetries satisfy the full plane boundary exactly, so the
    // half-plane cannot scatter: tangential forces for rigid, normal for release.
    let source = [-0.4, 0.0, 0.2];
    let point = [0.3, -0.2, 0.1];
    for (boundary, axis) in [
        (IdealPlaneBoundary::Rigid, 0),
        (IdealPlaneBoundary::Rigid, 2),
        (IdealPlaneBoundary::PressureRelease, 1),
    ] {
        for k in [0.03, 7.0, 120.0] {
            let mut f = [C64::ZERO; 3];
            f[axis] = C64::new(0.7, -0.2);
            let owner = prepare(source, point, k, boundary);
            let a = values(owner.evaluate_harmonic_line(f).unwrap());
            let e = original_force(source, point, k, f);
            for i in 0..4 {
                close(a[i], e[i], e[i].abs(), 3e-9);
            }
            let observed = owner.evaluate(f).unwrap();
            for i in 0..3 {
                let flux = a[0].scale(0.5) * a[i + 1].conj();
                assert_eq!(observed.active_intensity_w_m2[i], flux.re);
                assert_eq!(observed.reactive_intensity_w_m2[i], flux.im);
            }
        }
    }
}

#[test]
fn analytic_source_derivatives_match_original_monopole_in_light_and_shadow() {
    let source = [-0.4, 0.3, 0.2];
    let k = 7.0;
    let h = 2e-4;
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for point in [[0.3, 0.2, 0.1], [0.4, -0.2, -0.3], [-0.15, -0.25, 0.2]] {
            let force = prepare(source, point, k, boundary);
            for axis in 0..3 {
                let samples = [-2.0, -1.0, 1.0, 2.0].map(|step| {
                    let mut at = source;
                    at[axis] += step * h;
                    values(mono(at, point, k, boundary))
                });
                let reference: [C64; 4] = core::array::from_fn(|i| {
                    (samples[0][i] - samples[1][i].scale(8.0) + samples[2][i].scale(8.0)
                        - samples[3][i])
                        .scale(1.0 / (12.0 * h))
                        / C64::new(0.0, -k * 340.0 * 1.3)
                });
                let mut f = [C64::ZERO; 3];
                f[axis] = C64::ONE;
                let actual = values(force.evaluate_harmonic_line(f).unwrap());
                close(actual[0], reference[0], 1.0, 5e-8);
                for i in 1..4 {
                    close(actual[i], reference[i], 1.0 / (k * 340.0 * 1.3), 5e-8);
                }
            }
        }
    }
}

#[test]
fn full_velocity_satisfies_euler_continuity_and_outgoing_energy_balance() {
    let source = [0.4, 0.3, 0.2];
    let k = 7.0;
    let h = 1e-4;
    let f = [C64::new(0.7, -0.2), C64::new(-0.1, 0.3), C64::new(0.2, 0.1)];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for point in [[0.3, 0.2, 0.1], [0.4, -0.2, -0.3]] {
            let owner = prepare(source, point, k, boundary);
            let field = values(owner.evaluate_harmonic_line(f).unwrap());
            let mut divergence = C64::ZERO;
            let mut divergence_flux = 0.0;
            for axis in 0..3 {
                let samples = [-2.0, -1.0, 1.0, 2.0].map(|step| {
                    let mut at = point;
                    at[axis] += step * h;
                    let p = prepare(source, at, k, boundary);
                    (
                        values(p.evaluate_harmonic_line(f).unwrap()),
                        p.evaluate(f).unwrap().active_intensity_w_m2,
                    )
                });
                let derivative = |component: usize| {
                    (samples[0].0[component] - samples[1].0[component].scale(8.0)
                        + samples[2].0[component].scale(8.0)
                        - samples[3].0[component])
                        .scale(1.0 / (12.0 * h))
                };
                close(
                    derivative(0),
                    C64::new(0.0, k * 340.0 * 1.3) * field[axis + 1],
                    1.0,
                    2e-6,
                );
                divergence = divergence + derivative(axis + 1);
                divergence_flux += (samples[0].1[axis] - 8.0 * samples[1].1[axis]
                    + 8.0 * samples[2].1[axis]
                    - samples[3].1[axis])
                    / (12.0 * h);
            }
            close(
                divergence,
                C64::new(0.0, k / (1.3 * 340.0)) * field[0],
                1e-3,
                3e-6,
            );
            assert!(divergence_flux.abs() < 3e-8);
            let energy = owner.evaluate_energy_transport(f).unwrap();
            assert!(energy.energy_density_j_m3 > 0.0);
        }
    }
}

#[test]
fn exact_sheet_faces_and_shadow_boundaries_keep_all_force_directions() {
    let source = [0.4, 0.3, 0.2];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
            let owner = PreparedHalfPlaneForce::new(
                medium(),
                7.0 * 340.0 / TAU,
                frame(),
                boundary,
                source,
                HalfPlaneReceiver::SheetFace {
                    distance_from_edge_m: 0.3,
                    along_edge_m: 0.1,
                    face,
                },
                limits(),
                || false,
            )
            .unwrap();
            match boundary {
                IdealPlaneBoundary::Rigid => assert_eq!(owner.velocity[1], [C64::ZERO; 3]),
                IdealPlaneBoundary::PressureRelease => assert_eq!(owner.pressure, [C64::ZERO; 3]),
            }
        }
        // An occluded image coincidence is finite, including its mixed tensor.
        let image = prepare(source, [0.4, -0.3, 0.2], 7.0, boundary);
        assert!(image.pressure.into_iter().all(finite));
        assert!(image.velocity.into_iter().flatten().all(finite));

        // q_direct=0 exactly; the half-weight Hessian must agree with both
        // one-sided physical limits, not jump by the visible ray's Hessian.
        let source = [0.0, 1.0, 0.2];
        let point = [0.0, -1.0, -0.1];
        let exact = prepare(source, point, 7.0, boundary);
        for dx in [-1e-7, 1e-7, -1e-10, 1e-10, -1e-13, 1e-13] {
            let mut at = point;
            at[0] += dx;
            let near = prepare(source, at, 7.0, boundary);
            for j in 0..3 {
                close(exact.pressure[j], near.pressure[j], 1.0, 3e-6);
                for i in 0..3 {
                    close(exact.velocity[i][j], near.velocity[i][j], 1e-3, 3e-6);
                }
            }
        }
        let h = 2e-4;
        for axis in 0..3 {
            let samples = [-2.0, -1.0, 1.0, 2.0].map(|step| {
                let mut at = source;
                at[axis] += step * h;
                values(mono(at, point, 7.0, boundary))
            });
            for i in 0..3 {
                let derivative = (samples[0][i + 1] - samples[1][i + 1].scale(8.0)
                    + samples[2][i + 1].scale(8.0)
                    - samples[3][i + 1])
                    .scale(1.0 / (12.0 * h))
                    / C64::new(0.0, -7.0 * 340.0 * 1.3);
                close(exact.velocity[i][axis], derivative, 1e-3, 5e-7);
            }
        }
    }
}

#[test]
fn original_source_clock_rotated_tensor_and_reciprocity_are_preserved() {
    let source = [-0.4, 0.3, 0.2];
    let point = [0.3, -0.2, 0.1];
    let k = 7.0;
    let axes = [[0.8, 0.0, -0.6], [0.0, 1.0, 0.0], [0.6, 0.0, 0.8]];
    let origin = [2.0, -1.0, 0.7];
    let rotation = |v: [f64; 3]| {
        core::array::from_fn::<_, 3, _>(|i| {
            axes[0][i] * v[0] + axes[1][i] * v[1] + axes[2][i] * v[2]
        })
    };
    let placed = |v: [f64; 3]| {
        let v = rotation(v);
        core::array::from_fn(|i| origin[i] + v[i])
    };
    let f = [C64::new(0.7, -0.2), C64::new(-0.1, 0.3), C64::new(0.2, 0.1)];
    let phase = C64::new(0.7_f64.cos(), 0.7_f64.sin());
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let a = prepare(source, point, k, boundary);
        let reverse = prepare(point, source, k, boundary);
        let m = mono(point, source, k, boundary);
        for j in 0..3 {
            close(a.pressure[j], -m.particle_velocity_m_s[j], 1.0, 2e-10);
            for i in 0..3 {
                close(a.velocity[i][j], reverse.velocity[j][i], 1e-3, 2e-10);
            }
        }
        let b = PreparedHalfPlaneForce::new(
            medium(),
            k * 340.0 / TAU,
            HalfPlaneFrame::new(origin, axes[2], axes[1]).unwrap(),
            boundary,
            placed(source),
            HalfPlaneReceiver::WorldPosition(placed(point)),
            limits(),
            || false,
        )
        .unwrap();
        let original = a.evaluate_harmonic_line(f).unwrap();
        let placed = b
            .evaluate_harmonic_line(rotate(axes, f).map(|v| v * phase))
            .unwrap();
        close(placed.pressure_pa, original.pressure_pa * phase, 1.0, 3e-10);
        for (a, b) in placed
            .particle_velocity_m_s
            .into_iter()
            .zip(rotate(axes, original.particle_velocity_m_s))
        {
            close(a, b * phase, 1e-3, 3e-10);
        }
        assert_eq!(b.frequency_hz(), a.frequency_hz());
        assert_eq!(b.acoustic(), medium());
        assert_eq!(
            values(a.evaluate_harmonic_line([C64::ZERO; 3]).unwrap()),
            [C64::ZERO; 4]
        );
    }
}

#[test]
fn force_geometry_physical_units_work_limits_and_cancellation_refuse_explicitly() {
    let source = [-0.4, 0.3, 0.2];
    let point = [0.3, -0.2, 0.1];
    let try_owner = |s, x, f, limits| {
        PreparedHalfPlaneForce::new(
            medium(),
            f,
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            HalfPlaneReceiver::WorldPosition(x),
            limits,
            || false,
        )
    };
    for s in [[0.0, 0.0, 0.2], [0.4, 0.0, 0.2], [f64::NAN, 0.1, 0.0]] {
        assert!(try_owner(s, point, 100.0, limits()).is_err());
    }
    for x in [source, [0.0, 0.0, 0.2], [0.4, 0.0, 0.1]] {
        assert!(try_owner(source, x, 100.0, limits()).is_err());
    }
    for frequency in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(try_owner(source, point, frequency, limits()).is_err());
    }
    let mut short = limits();
    short.maximum_evaluations = 15;
    assert!(try_owner(source, point, 100.0, short).is_err());
    let mut wrong = limits();
    wrong.absolute_velocity_tolerance_m_s_per_n = 0.0;
    assert!(try_owner(source, point, 100.0, wrong).is_err());
    let owner = prepare(source, point, 7.0, IdealPlaneBoundary::Rigid);
    assert!(owner.evaluate([C64::new(f64::NAN, 0.0); 3]).is_err());
    let mut polls = 0;
    let error = PreparedHalfPlaneForce::new(
        medium(),
        100.0,
        frame(),
        IdealPlaneBoundary::Rigid,
        source,
        HalfPlaneReceiver::WorldPosition(point),
        limits(),
        || {
            polls += 1;
            polls >= 4
        },
    );
    assert!(matches!(error, Err(AeroacError::Cancelled { .. })));
    assert_eq!(polls, 4);
}
