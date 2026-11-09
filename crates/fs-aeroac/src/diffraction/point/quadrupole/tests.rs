//! Independent physical oracles for all six stress entries and full Euler motion.
use super::*;
use crate::diffraction::HalfPlaneFace;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 340.0).unwrap()
}
fn frame() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}
fn limits() -> HalfPlaneQuadrupoleLimits {
    HalfPlaneQuadrupoleLimits {
        absolute_pressure_tolerance_pa_per_n_m: 1e-10,
        absolute_velocity_tolerance_m_s_per_n_m: 1e-12,
        relative_tolerance: 2e-10,
        maximum_evaluations: 200_000,
        maximum_depth: 42,
    }
}
fn prepare(s: [f64; 3], x: [f64; 3], k: f64, b: IdealPlaneBoundary) -> PreparedHalfPlaneQuadrupole {
    PreparedHalfPlaneQuadrupole::new(
        medium(),
        k * 340.0 / TAU,
        frame(),
        b,
        s,
        HalfPlaneReceiver::WorldPosition(x),
        limits(),
        || false,
    )
    .unwrap_or_else(|e| panic!("{s:?}->{x:?},k={k},{b:?}: {e:?}"))
}
fn values(line: HarmonicFieldLine<3>) -> [C64; 4] {
    [
        line.pressure_pa,
        line.particle_velocity_m_s[0],
        line.particle_velocity_m_s[1],
        line.particle_velocity_m_s[2],
    ]
}
fn close(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(b.abs()).max(1e-30),
        "{a:?} != {b:?}, scaled error={}",
        (a - b).abs() / scale.max(b.abs()).max(1e-30)
    );
}
fn stress() -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.7, 0.1),
            C64::new(-0.3, 0.2),
            C64::new(0.9, -0.4),
            C64::new(0.23, -0.15),
            C64::new(-0.11, 0.3),
            C64::new(0.4, 0.2),
        ],
    }
}
// Direct spherical Hessian and third derivative, independent of the contour
// and its local derivative algebra. Physical shear is counted twice.
fn free(s: [f64; 3], x: [f64; 3], k: f64, m: QuadrupolePhasor) -> [C64; 4] {
    let d = core::array::from_fn::<_, 3, _>(|i| x[i] - s[i]);
    let r = d.into_iter().fold(0.0_f64, f64::hypot);
    let n = d.map(|v| v / r);
    let g = C64::new((k * r).cos(), (k * r).sin()).scale(1.0 / (4.0 * PI * r));
    let a = C64::new(3.0 / (r * r) - k * k, -3.0 * k / r);
    let b = C64::new(-1.0 / (r * r), k / r);
    let c = C64::new(
        6.0 * k * k / r - 15.0 / (r * r * r),
        15.0 * k / (r * r) - k * k * k,
    );
    let mut out = [C64::ZERO; 4];
    for (entry, (i, j)) in derivatives::PAIRS.into_iter().enumerate() {
        let multiplicity = if i == j { 1.0 } else { 2.0 };
        let q = m.stress_moment_n_m[entry].scale(multiplicity);
        out[0] = out[0] + g * (a.scale(n[i] * n[j]) + b.scale(if i == j { 1.0 } else { 0.0 })) * q;
        for l in 0..3 {
            let angular = (if i == j { n[l] } else { 0.0 })
                + (if i == l { n[j] } else { 0.0 })
                + (if j == l { n[i] } else { 0.0 });
            let third = g * (c.scale(n[i] * n[j] * n[l]) + a.scale(angular / r));
            out[l + 1] = out[l + 1] + third / C64::new(0.0, 1.3 * 340.0 * k) * q;
        }
    }
    out
}
#[test]
fn all_physical_tensor_entries_recover_independent_free_fields() {
    let s = [-0.4, 0.0, 0.2];
    let x = [0.3, -0.2, 0.1];
    for k in [0.03, 7.0, 120.0] {
        for b in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let owner = prepare(s, x, k, b);
            eprintln!(
                "stress free k={k} {b:?}: {} evaluations",
                owner.diagnostics.evaluations
            );
            for (j, (i, l)) in derivatives::PAIRS.into_iter().enumerate() {
                let odd = (i == 1) ^ (l == 1);
                if odd != (b == IdealPlaneBoundary::PressureRelease) {
                    continue;
                }
                let mut m = QuadrupolePhasor::default();
                m.stress_moment_n_m[j] = C64::new(0.7, -0.2);
                let a = values(owner.evaluate_harmonic_line(m).unwrap());
                let e = free(s, x, k, m);
                for c in 0..4 {
                    close(a[c], e[c], e[c].abs().max(1e-10), 2e-8);
                }
            }
        }
    }
}
#[test]
fn original_helmholtz_trace_survives_light_shadow_and_exact_boundary() {
    for (s, x) in [
        ([0.4, 0.3, 0.2], [0.3, 0.2, 0.1]),
        ([0.4, 0.3, 0.2], [0.4, -0.2, -0.3]),
        ([0.0, 1.0, 0.2], [0.0, -1.0, -0.1]),
    ] {
        for b in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let k = 7.0;
            let owner = prepare(s, x, k, b);
            let mono = PreparedHalfPlaneMonopole::new(
                medium(),
                k * 340.0 / TAU,
                frame(),
                b,
                s,
                HalfPlaneReceiver::WorldPosition(x),
                HalfPlaneMonopoleLimits {
                    absolute_pressure_tolerance_pa_per_m3_s: 1e-8,
                    absolute_velocity_tolerance_per_m2: 1e-11,
                    relative_tolerance: 1e-11,
                    maximum_evaluations: 100_000,
                    maximum_depth: 42,
                },
                || false,
            )
            .unwrap();
            let moment = C64::new(0.7, -0.2);
            let m = QuadrupolePhasor {
                stress_moment_n_m: [moment, moment, moment, C64::ZERO, C64::ZERO, C64::ZERO],
            };
            let a = values(owner.evaluate_harmonic_line(m).unwrap());
            let factor = C64::new(0.0, -k / (1.3 * 340.0)) * moment;
            let e = values(mono.evaluate_harmonic_line(factor).unwrap());
            for c in 0..4 {
                close(a[c], e[c], if c == 0 { 1.0 } else { 1e-3 }, 3e-8);
            }
        }
    }
}

fn force_map(s: [f64; 3], x: [f64; 3], k: f64, b: IdealPlaneBoundary) -> [[C64; 3]; 4] {
    let owner = PreparedHalfPlaneForce::new(
        medium(),
        k * 340.0 / TAU,
        frame(),
        b,
        s,
        HalfPlaneReceiver::WorldPosition(x),
        HalfPlaneForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-12,
            absolute_velocity_tolerance_m_s_per_n: 1e-14,
            relative_tolerance: 1e-11,
            maximum_evaluations: 100_000,
            maximum_depth: 42,
        },
        || false,
    )
    .unwrap();
    [
        *owner.pressure_transfer_pa_per_n(),
        owner.velocity_transfer_m_s_per_n()[0],
        owner.velocity_transfer_m_s_per_n()[1],
        owner.velocity_transfer_m_s_per_n()[2],
    ]
}
#[test]
fn all_source_hessians_and_mixed_third_derivatives_match_original_force_fields() {
    let k = 7.0;
    let h = 2e-4;
    for (s, x) in [
        ([0.4, 0.3, 0.2], [0.3, 0.2, 0.1]),
        ([0.4, 0.3, 0.2], [0.4, -0.2, -0.3]),
        ([0.0, 1.0, 0.2], [0.0, -1.0, -0.1]),
    ] {
        for b in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let owner = prepare(s, x, k, b);
            let d = core::array::from_fn::<_, 3, _>(|axis| {
                let sample = [-2.0, -1.0, 1.0, 2.0].map(|step| {
                    let mut at = s;
                    at[axis] += step * h;
                    force_map(at, x, k, b)
                });
                core::array::from_fn::<_, 4, _>(|row| {
                    core::array::from_fn::<_, 3, _>(|col| {
                        (sample[0][row][col] - sample[1][row][col].scale(8.0)
                            + sample[2][row][col].scale(8.0)
                            - sample[3][row][col])
                            .scale(1.0 / (12.0 * h))
                    })
                })
            });
            for (col, (i, j)) in derivatives::PAIRS.into_iter().enumerate() {
                for row in 0..4 {
                    let reference = if i == j {
                        d[j][row][i]
                    } else {
                        d[j][row][i] + d[i][row][j]
                    };
                    close(
                        owner.transfer[row][col],
                        reference,
                        if row == 0 { 1.0 } else { 1e-3 },
                        3e-6,
                    );
                }
            }
        }
    }
}
#[test]
fn full_near_field_closes_momentum_compression_and_energy_flux() {
    let s = [0.4, 0.3, 0.2];
    let k = 7.0;
    let h = 1e-4;
    let m = stress();
    for (x, b) in [
        ([0.3, 0.2, 0.1], IdealPlaneBoundary::Rigid),
        ([0.4, -0.2, -0.3], IdealPlaneBoundary::PressureRelease),
    ] {
        let owner = prepare(s, x, k, b);
        let f = values(owner.evaluate_harmonic_line(m).unwrap());
        let flux = owner.evaluate(m).unwrap();
        let energy = owner.evaluate_energy_transport(m).unwrap();
        assert!(energy.energy_density_j_m3 > 0.0);
        let mut divergence = C64::ZERO;
        let mut flux_divergence = 0.0;
        for axis in 0..3 {
            let samples = [-2.0, -1.0, 1.0, 2.0].map(|step| {
                let mut at = x;
                at[axis] += step * h;
                let observer = prepare(s, at, k, b);
                (
                    values(observer.evaluate_harmonic_line(m).unwrap()),
                    observer.evaluate(m).unwrap().active_intensity_w_m2,
                )
            });
            let d = |row: usize| {
                (samples[0].0[row] - samples[1].0[row].scale(8.0) + samples[2].0[row].scale(8.0)
                    - samples[3].0[row])
                    .scale(1.0 / (12.0 * h))
            };
            close(
                d(0),
                C64::new(0.0, k * 340.0 * 1.3) * f[axis + 1],
                1.0,
                3e-6,
            );
            divergence = divergence + d(axis + 1);
            flux_divergence += (samples[0].1[axis] - 8.0 * samples[1].1[axis]
                + 8.0 * samples[2].1[axis]
                - samples[3].1[axis])
                / (12.0 * h);
            let product = f[0].scale(0.5) * f[axis + 1].conj();
            assert_eq!(flux.active_intensity_w_m2[axis], product.re);
            assert_eq!(flux.reactive_intensity_w_m2[axis], product.im);
        }
        close(
            divergence,
            C64::new(0.0, k / (1.3 * 340.0)) * f[0],
            1e-3,
            5e-6,
        );
        let flux_scale = flux
            .active_intensity_w_m2
            .into_iter()
            .map(f64::abs)
            .sum::<f64>();
        assert!(
            flux_divergence.abs() < 3e-6 * flux_scale.max(1e-8),
            "{flux_divergence}"
        );
    }
}
#[test]
fn sheet_faces_and_arbitrarily_narrow_shadow_transitions_are_regular() {
    let s = [0.4, 0.3, 0.2];
    let k = 7.0;
    for b in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
            let observer = PreparedHalfPlaneQuadrupole::new(
                medium(),
                k * 340.0 / TAU,
                frame(),
                b,
                s,
                HalfPlaneReceiver::SheetFace {
                    distance_from_edge_m: 0.3,
                    along_edge_m: 0.1,
                    face,
                },
                limits(),
                || false,
            )
            .unwrap();
            let row = if b == IdealPlaneBoundary::Rigid { 2 } else { 0 };
            for &coefficient in &observer.transfer[row] {
                assert!(coefficient.abs() < 2e-11, "{b:?} {face:?}: {coefficient:?}");
            }
        }
        let occluded = prepare(s, [0.4, -0.3, 0.2], k, b);
        assert!(occluded.transfer.into_iter().flatten().all(finite));
        let s = [0.0, 1.0, 0.2];
        let x = [0.0, -1.0, -0.1];
        let exact = prepare(s, x, k, b);
        for dx in [-1e-7, 1e-7, -1e-10, 1e-10, -1e-13, 1e-13] {
            let mut at = x;
            at[0] += dx;
            let near = prepare(s, at, k, b);
            for row in 0..4 {
                for col in 0..6 {
                    close(
                        near.transfer[row][col],
                        exact.transfer[row][col],
                        if row == 0 { 1.0 } else { 1e-3 },
                        4e-6,
                    );
                }
            }
        }
    }
}
#[test]
fn original_world_tensor_source_clock_and_coherent_null_are_retained() {
    let s = [0.4, 0.3, 0.2];
    let x = [0.4, -0.2, -0.3];
    let m = stress();
    let k = 7.0;
    let rotation = [
        [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    ];
    let origin = [1.3, -0.4, 0.7];
    let rotate_real = |v: [f64; 3]| rotation.map(|row| (0..3).map(|i| row[i] * v[i]).sum::<f64>());
    let place = |v: [f64; 3]| {
        let r = rotate_real(v);
        core::array::from_fn(|i| origin[i] + r[i])
    };
    let mut tensor = [[C64::ZERO; 3]; 3];
    for (entry, (i, j)) in derivatives::PAIRS.into_iter().enumerate() {
        tensor[i][j] = m.stress_moment_n_m[entry];
        tensor[j][i] = tensor[i][j];
    }
    let rotated = QuadrupolePhasor {
        stress_moment_n_m: derivatives::PAIRS.map(|(i, j)| {
            let mut v = C64::ZERO;
            for a in 0..3 {
                for b in 0..3 {
                    v = v + tensor[a][b].scale(rotation[i][a] * rotation[j][b]);
                }
            }
            v
        }),
    };
    let placed = HalfPlaneFrame::new(
        origin,
        rotate_real([0.0, 0.0, 1.0]),
        rotate_real([0.0, 1.0, 0.0]),
    )
    .unwrap();
    for b in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let original = prepare(s, x, k, b);
        let world = PreparedHalfPlaneQuadrupole::new(
            medium(),
            k * 340.0 / TAU,
            placed,
            b,
            place(s),
            HalfPlaneReceiver::WorldPosition(place(x)),
            limits(),
            || false,
        )
        .unwrap();
        let a = values(original.evaluate_harmonic_line(m).unwrap());
        let e = values(world.evaluate_harmonic_line(rotated).unwrap());
        close(e[0], a[0], 1.0, 3e-8);
        for axis in 0..3 {
            let expected = (0..3).fold(C64::ZERO, |v, j| v + a[j + 1].scale(rotation[axis][j]));
            close(e[axis + 1], expected, 1e-3, 3e-8);
        }
        let phase = C64::new(0.0, 1.0);
        let rotated_clock = QuadrupolePhasor {
            stress_moment_n_m: m.stress_moment_n_m.map(|v| v * phase),
        };
        let p = values(original.evaluate_harmonic_line(rotated_clock).unwrap());
        for i in 0..4 {
            close(p[i], a[i] * phase, if i == 0 { 1.0 } else { 1e-3 }, 1e-14);
        }
        let zero = values(
            original
                .evaluate_harmonic_line(QuadrupolePhasor::default())
                .unwrap(),
        );
        assert!(zero.into_iter().all(|v| v == C64::ZERO));
        assert_eq!(original.acoustic(), medium());
        assert_eq!(original.source_position_m(), s);
        match original.receiver() {
            HalfPlaneReceiver::WorldPosition(at) => assert_eq!(at, x),
            _ => panic!("original receiver changed"),
        }
    }
}
#[test]
fn physical_inputs_work_limits_and_bounded_cancellation_refuse() {
    let s = [0.4, 0.3, 0.2];
    let x = [0.4, -0.2, -0.3];
    let make = |source, receiver, frequency, control, cancel: fn() -> bool| {
        PreparedHalfPlaneQuadrupole::new(
            medium(),
            frequency,
            frame(),
            IdealPlaneBoundary::Rigid,
            source,
            HalfPlaneReceiver::WorldPosition(receiver),
            control,
            cancel,
        )
    };
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(make(s, x, frequency, limits(), || false).is_err());
    }
    for source in [[0.0; 3], [0.4, 0.0, 0.2], [f64::NAN, 0.3, 0.2], x] {
        assert!(make(source, x, 300.0, limits(), || false).is_err());
    }
    let mut budget = limits();
    budget.maximum_evaluations = 1;
    assert!(make(s, x, 300.0, budget, || false).is_err());
    let mut invalid_limits = limits();
    invalid_limits.relative_tolerance = 1.0;
    assert!(make(s, x, 300.0, invalid_limits, || false).is_err());
    assert!(matches!(
        make(s, x, 300.0, limits(), || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut calls = 0;
    assert!(matches!(
        PreparedHalfPlaneQuadrupole::new(
            medium(),
            300.0,
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            HalfPlaneReceiver::WorldPosition(x),
            limits(),
            || {
                calls += 1;
                calls >= 20
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(calls, 20);
    let observer = prepare(s, x, 7.0, IdealPlaneBoundary::Rigid);
    let mut bad = stress();
    bad.stress_moment_n_m[3] = C64::new(f64::NAN, 0.0);
    assert!(observer.evaluate(bad).is_err());
}
