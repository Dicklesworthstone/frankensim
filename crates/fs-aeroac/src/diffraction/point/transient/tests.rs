//! Independent real-time field, causality, refinement and reciprocal identities.
use super::*;
use crate::curle3d::moving::CompactSourceStrength;
use crate::diffraction::HalfPlaneFace;
fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.7, 1.3).unwrap()
}
fn frame() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}
fn limits() -> HalfPlaneTransientLimits {
    HalfPlaneTransientLimits {
        quadrature: HalfPlaneMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 2e-9,
            absolute_velocity_tolerance_per_m2: 2e-10,
            relative_tolerance: 2e-10,
            maximum_evaluations: 8000,
            maximum_depth: 42,
        },
        maximum_total_evaluations: 2_000_000,
        maximum_payload_bytes: 64_000_000,
    }
}
fn prepare(
    s: [f64; 3],
    x: HalfPlaneReceiver,
    b: IdealPlaneBoundary,
    dt: f64,
    n: usize,
) -> PreparedHalfPlaneTransient {
    PreparedHalfPlaneTransient::new(medium(), frame(), b, s, x, dt, n, limits(), || false).unwrap()
}
fn history(dt: f64, n: usize, start: f64, delay: f64, degree: i32) -> Vec<TimedCompactSource> {
    (0..n)
        .map(|i| {
            let time = (Dd::from_f64(start) + Dd::from_f64(dt) * Dd::from_f64(i as f64)).to_f64();
            let age = (i as f64 * dt - delay).max(0.0);
            TimedCompactSource {
                source_time_s: time,
                strength: CompactSourceStrength {
                    volume_flow_m3_s: age.powi(degree),
                    volume_flow_rate_m3_s2: degree as f64 * age.powi(degree - 1),
                    ..CompactSourceStrength::default()
                },
            }
        })
        .collect()
}
fn values(field: HalfPlaneTransientField) -> [f64; 4] {
    [
        field.pressure_pa,
        field.particle_velocity_m_s[0],
        field.particle_velocity_m_s[1],
        field.particle_velocity_m_s[2],
    ]
}
fn error(actual: [f64; 4], expected: [f64; 4]) -> f64 {
    let mut result = ((actual[0] - expected[0]).abs() / expected[0].abs().max(1e-13)).abs();
    let scale = expected[1..]
        .iter()
        .map(|v| v * v)
        .sum::<f64>()
        .sqrt()
        .max(1e-13);
    result = result.max(
        actual[1..]
            .iter()
            .zip(&expected[1..])
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
            / scale,
    );
    result
}
fn free(source: [f64; 3], point: [f64; 3], time: f64, delay: f64, degree: i32) -> [f64; 4] {
    let delta = core::array::from_fn::<_, 3, _>(|i| point[i] - source[i]);
    let r = delta.iter().map(|x| x * x).sum::<f64>().sqrt();
    let t = (time - delay - r / medium().sound_speed_m_s()).max(0.0);
    let q = t.powi(degree);
    let dq = degree as f64 * t.powi(degree - 1);
    let p = medium().density_kg_m3() * dq / (4.0 * PI * r);
    let radial = (q / r + dq / medium().sound_speed_m_s()) / (4.0 * PI * r);
    [
        p,
        radial * delta[0] / r,
        radial * delta[1] / r,
        radial * delta[2] / r,
    ]
}

#[test]
fn causal_free_spherical_field_cubic_history_and_fft_linear_convolution_agree() {
    let s = [-0.04, 0.0, 0.012];
    let x = [0.024, 0.038, -0.02];
    let dt = 0.0025;
    let n = 192;
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let source = history(dt, n, 0.0, 0.0, 3);
    let all = owner.observe_record(&source, || false).unwrap();
    for i in 0..n {
        let expected = free(s, x, i as f64 * dt, 0.0, 3);
        let direct = values(owner.observe_at(&source, i, || false).unwrap());
        if expected == [0.0; 4] {
            assert_eq!(values(all[i]), [0.0; 4]);
            assert_eq!(direct, [0.0; 4]);
        } else {
            assert!(
                error(direct, expected) < 2e-10,
                "direct {i}: {direct:?} != {expected:?}"
            );
            assert!(
                error(values(all[i]), expected) < 2e-8,
                "FFT {i}: {:?} != {expected:?}",
                values(all[i])
            );
        }
    }
    // A long quiescent prefix still has exact zeros before its delayed
    // physical arrival; circular FFT roundoff must not create a pre-echo.
    let delayed = history(dt, n, 0.0, 0.25, 3);
    let delayed_output = owner.observe_record(&delayed, || false).unwrap();
    let first = delayed
        .iter()
        .position(|sample| sample.strength.volume_flow_m3_s != 0.0)
        .unwrap();
    let start = delayed[first - 1].source_time_s;
    for (sample, field) in delayed.iter().zip(&delayed_output) {
        if sample.source_time_s - start < owner.diagnostics().minimum_arrival_delay_s {
            assert_eq!(values(*field), [0.0; 4]);
        } else {
            let expected = free(s, x, sample.source_time_s, 0.25, 3);
            assert!(error(values(*field), expected) < 2e-8);
        }
    }
    let zero = history(dt, n, 0.0, 1.0, 3);
    assert!(
        owner
            .observe_record(&zero, || false)
            .unwrap()
            .iter()
            .all(|f| values(*f) == [0.0; 4])
    );
    // A valid large absolute epoch rounds adjacent timestamps differently
    // from the declared relative clock; it must not erase a real arrival.
    let large_dt = 4.3;
    let large_source = [-1.0, 0.0, 0.0];
    let distance = medium().sound_speed_m_s() * 4.15;
    let large_point = [0.0, (distance * distance - 1.0).sqrt(), 0.0];
    let large = prepare(
        large_source,
        HalfPlaneReceiver::WorldPosition(large_point),
        IdealPlaneBoundary::Rigid,
        large_dt,
        2,
    );
    let large_history = history(large_dt, 2, (1u64 << 54) as f64, 0.0, 3);
    assert_eq!(
        large_history[1].source_time_s - large_history[0].source_time_s,
        4.0
    );
    let large_batch = large.observe_record(&large_history, || false).unwrap();
    let large_direct = large.observe_at(&large_history, 1, || false).unwrap();
    let expected = free(large_source, large_point, large_dt, 0.0, 3);
    assert!(large_direct.pressure_pa > 0.0);
    assert!(error(values(large_direct), expected) < 2e-10);
    assert!(error(values(large_batch[1]), expected) < 2e-8);
    let i = 151;
    let p = all[i].pressure_pa;
    let u = all[i].particle_velocity_m_s;
    assert_eq!(
        all[i].instantaneous_intensity_w_m2,
        u.map(|v| (Dd::from_f64(p) * Dd::from_f64(v)).to_f64())
    );
    let energy = p * p / (2.0 * medium().adiabatic_bulk_modulus_pa())
        + 0.5 * medium().density_kg_m3() * u.iter().map(|v| v * v).sum::<f64>();
    assert!((all[i].perturbation_energy_j_m3 - energy).abs() < 2e-15 * energy);
}

// Independent UNINTEGRATED-BY-PARTS Biot--Tolstoy real kernel. Geometry uses
// ordinary polar angles, and a composite Simpson integral supplies its own
// pressure/gradient reference; it uses no production taps, quadrature or FFT.
fn reference(
    source: [f64; 3],
    point: [f64; 3],
    time: f64,
    degree: i32,
    boundary: IdealPlaneBoundary,
) -> [f64; 4] {
    let rho = medium().density_kg_m3();
    let c = medium().sound_speed_m_s();
    let r = point[0].hypot(point[1]);
    let rs = source[0].hypot(source[1]);
    let dz = point[2] - source[2];
    let l = (r + rs).hypot(dz);
    let tr = point[1].atan2(point[0]).rem_euclid(2.0 * PI);
    let ts = source[1].atan2(source[0]).rem_euclid(2.0 * PI);
    let radial = [point[0] / r, point[1] / r, 0.0];
    let tangent = [-point[1] / r, point[0] / r, 0.0];
    let dl = [(r + rs) / l * radial[0], (r + rs) / l * radial[1], dz / l];
    let mut result = [0.0; 4];
    for j in 0..2 {
        let angle = if j == 0 { tr - ts } else { tr + ts };
        let q = 2.0 * (r * rs).sqrt() * (0.5 * angle).cos();
        let dq = core::array::from_fn::<_, 3, _>(|i| {
            (rs / r).sqrt() * ((0.5 * angle).cos() * radial[i] - (0.5 * angle).sin() * tangent[i])
        });
        let sign = if j == 0 { 1.0 } else { boundary.sign() };
        if q > 0.0 {
            let image = [
                source[0],
                if j == 0 { source[1] } else { -source[1] },
                source[2],
            ];
            let field = free(image, point, time, 0.0, degree);
            for i in 0..4 {
                result[i] += sign * field[i];
            }
        }
        if c * time <= l {
            continue;
        }
        let upper = ((c * time - l) * (c * time + l)).sqrt();
        let count = 32_768;
        let step = upper / count as f64;
        let mut integral = [0.0; 4];
        for k in 0..=count {
            let s = k as f64 * step;
            let radius = l.hypot(s);
            let age = (time - radius / c).max(0.0);
            let value = age.powi(degree);
            let first = degree as f64 * age.powi(degree - 1);
            let denominator = s * s + q * q;
            let pressure = -rho * q / denominator * first / radius;
            let gradient = -(value + radius / c * first) / (radius * radius);
            let weight = if k == 0 || k == count {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral[0] += weight * pressure;
            for axis in 0..3 {
                integral[axis + 1] += weight
                    * (q / denominator * gradient * l / radius * dl[axis]
                        + (s * s - q * q) / (denominator * denominator) * value / radius
                            * dq[axis]);
            }
        }
        for i in 0..4 {
            result[i] += sign * step / 3.0 * integral[i] / (4.0 * PI * PI);
        }
    }
    result
}

#[test]
fn complete_causal_edge_tail_matches_independent_real_biot_tolstoy_integral() {
    let source = [-0.04, 0.03, 0.01];
    let dt = 0.004;
    let n = 81;
    let history = history(dt, n, 0.0, 0.0, 3);
    for point in [[0.03, -0.02, -0.01], [-0.025, -0.07, 0.03]] {
        for b in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let owner = prepare(source, HalfPlaneReceiver::WorldPosition(point), b, dt, n);
            for index in [20, 40, 80] {
                let actual = values(owner.observe_at(&history, index, || false).unwrap());
                let expected = reference(source, point, index as f64 * dt, 3, b);
                assert!(
                    error(actual, expected) < 3e-8,
                    "{point:?} {b:?} t={index}: {actual:?} != {expected:?}, err={}",
                    error(actual, expected)
                );
            }
        }
    }
}

#[test]
fn temporal_hermite_refinement_converges_to_independent_continuous_source() {
    let s = [-0.04, 0.03, 0.01];
    let x = [0.03, -0.02, -0.01];
    let time = 0.32;
    let expected = reference(s, x, time, 4, IdealPlaneBoundary::Rigid);
    let mut errors = Vec::new();
    for dt in [0.01, 0.005, 0.0025] {
        let n = (time / dt).round() as usize + 1;
        let owner = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(x),
            IdealPlaneBoundary::Rigid,
            dt,
            n,
        );
        let source = history(dt, n, 0.0, 0.0, 4);
        let actual = values(owner.observe_at(&source, n - 1, || false).unwrap());
        errors.push(error(actual, expected));
    }
    eprintln!("causal source-interpolation ladder {errors:?}");
    assert!(errors[1] < errors[0] * 0.4 && errors[2] < errors[1] * 0.4);
    assert!(errors[2] < 2e-6);
}

#[test]
fn exact_sheet_faces_and_narrow_shadow_transitions_retain_causal_fields() {
    let source = [0.03, 0.04, 0.01];
    let dt = 0.004;
    let n = 81;
    let input = history(dt, n, 0.0, 0.0, 3);
    for b in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
            let owner = prepare(
                source,
                HalfPlaneReceiver::SheetFace {
                    distance_from_edge_m: 0.06,
                    along_edge_m: -0.02,
                    face,
                },
                b,
                dt,
                n,
            );
            let out = owner.observe_at(&input, 80, || false).unwrap();
            if b == IdealPlaneBoundary::Rigid {
                assert!(out.particle_velocity_m_s[1].abs() < 1e-12);
            } else {
                assert!(out.pressure_pa.abs() < 1e-12);
            }
        }
        let at = [-0.03, 0.04, -0.02];
        let exact = values(
            prepare(source, HalfPlaneReceiver::WorldPosition(at), b, dt, n)
                .observe_at(&input, 80, || false)
                .unwrap(),
        );
        for perturbation in [-1e-9, -1e-12, 1e-12, 1e-9] {
            let mut point = at;
            point[1] += perturbation;
            let actual = values(
                prepare(source, HalfPlaneReceiver::WorldPosition(point), b, dt, n)
                    .observe_at(&input, 80, || false)
                    .unwrap(),
            );
            assert!(
                error(actual, exact) < 2e-6,
                "shadow {b:?} {perturbation}: {actual:?} {exact:?}"
            );
        }
    }
}

#[test]
fn causal_fields_close_momentum_continuity_and_instantaneous_energy() {
    let s = [0.03, 0.04, -0.01];
    let x = [-0.02, -0.06, 0.04];
    let dt = 0.002;
    let n = 161;
    let time = (n - 1) as f64 * dt;
    let h = 0.00001;
    let start = 0.02;
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let central = owner
        .observe_at(&history(dt, n, 0.0, start, 3), n - 1, || false)
        .unwrap();
    let plus = owner
        .observe_at(&history(dt, n, 0.0, start - h, 3), n - 1, || false)
        .unwrap();
    let minus = owner
        .observe_at(&history(dt, n, 0.0, start + h, 3), n - 1, || false)
        .unwrap();
    let dpdt = (plus.pressure_pa - minus.pressure_pa) / (2.0 * h);
    let dudt = core::array::from_fn::<_, 3, _>(|i| {
        (plus.particle_velocity_m_s[i] - minus.particle_velocity_m_s[i]) / (2.0 * h)
    });
    let dedt = (plus.perturbation_energy_j_m3 - minus.perturbation_energy_j_m3) / (2.0 * h);
    let mut gradient = [0.0; 3];
    let mut divergence = 0.0;
    let mut power_divergence = 0.0;
    for axis in 0..3 {
        let mut p = x;
        let mut m = x;
        p[axis] += h;
        m[axis] -= h;
        let fp = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(p),
            IdealPlaneBoundary::Rigid,
            dt,
            n,
        )
        .observe_at(&history(dt, n, 0.0, start, 3), n - 1, || false)
        .unwrap();
        let fm = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(m),
            IdealPlaneBoundary::Rigid,
            dt,
            n,
        )
        .observe_at(&history(dt, n, 0.0, start, 3), n - 1, || false)
        .unwrap();
        gradient[axis] = (fp.pressure_pa - fm.pressure_pa) / (2.0 * h);
        divergence += (fp.particle_velocity_m_s[axis] - fm.particle_velocity_m_s[axis]) / (2.0 * h);
        power_divergence += (fp.instantaneous_intensity_w_m2[axis]
            - fm.instantaneous_intensity_w_m2[axis])
            / (2.0 * h);
    }
    let rho = medium().density_kg_m3();
    let k = medium().adiabatic_bulk_modulus_pa();
    let residual = gradient
        .iter()
        .zip(dudt)
        .map(|(g, u)| (g + rho * u).powi(2))
        .sum::<f64>()
        .sqrt();
    let scale = gradient.iter().map(|v| v * v).sum::<f64>().sqrt();
    eprintln!(
        "causal PDE t={time}, momentum={} continuity={} energy={}",
        residual / scale,
        (dpdt / k + divergence).abs() / (dpdt / k).abs(),
        (dedt + power_divergence).abs() / dedt.abs()
    );
    assert!(residual < scale * 2e-4);
    assert!((dpdt / k + divergence).abs() < 2e-4 * (dpdt / k).abs());
    assert!((dedt + power_divergence).abs() < 3e-4 * dedt.abs());
    assert!(central.perturbation_energy_j_m3 > 0.0);
}

#[test]
fn original_world_frame_clock_linearity_and_complete_fft_tail_are_preserved() {
    let s = [-0.04, 0.03, 0.01];
    let x = [0.03, -0.02, -0.01];
    let dt = 0.004;
    let n = 96;
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let transform = |v: [f64; 3]| [-v[1] + 0.2, v[2] - 0.3, -v[0] + 0.1];
    let placed = HalfPlaneFrame::new([0.2, -0.3, 0.1], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]).unwrap();
    let other = PreparedHalfPlaneTransient::new(
        medium(),
        placed,
        IdealPlaneBoundary::Rigid,
        transform(s),
        HalfPlaneReceiver::WorldPosition(transform(x)),
        dt,
        n,
        limits(),
        || false,
    )
    .unwrap();
    let input = history(dt, n, 2.0, 0.0, 3);
    let a = owner.observe_record(&input, || false).unwrap();
    let b = other.observe_record(&input, || false).unwrap();
    for index in [30, 60, 95] {
        let expect = [
            a[index].pressure_pa,
            -a[index].particle_velocity_m_s[1],
            a[index].particle_velocity_m_s[2],
            -a[index].particle_velocity_m_s[0],
        ];
        assert!(error(values(b[index]), expect) < 2e-10);
        assert_eq!(b[index].receiver_time_s, input[index].source_time_s);
        let direct = owner.observe_at(&input, index, || false).unwrap();
        assert!(error(values(a[index]), values(direct)) < 2e-10);
    }
    // Future accepted strengths cannot influence an earlier output.
    let mut changed = input.clone();
    for sample in &mut changed[61..] {
        sample.strength.volume_flow_m3_s *= 13.0;
        sample.strength.volume_flow_rate_m3_s2 *= -7.0;
    }
    let c = owner.observe_record(&changed, || false).unwrap();
    for index in [30, 45, 60] {
        assert!(error(values(c[index]), values(a[index])) < 2e-10);
    }
    // A C1 compact pulse leaves a physical diffraction tail after every direct
    // source image has returned to zero.
    let mut pulse = input.clone();
    let duration = 0.08;
    for (i, sample) in pulse.iter_mut().enumerate() {
        let t = i as f64 * dt;
        if t < duration {
            let z = t / duration;
            sample.strength.volume_flow_m3_s = z * z * (1.0 - z) * (1.0 - z);
            sample.strength.volume_flow_rate_m3_s2 =
                2.0 * z * (1.0 - z) * (1.0 - 2.0 * z) / duration;
        } else {
            sample.strength = CompactSourceStrength::default();
        }
    }
    let tail = owner.observe_record(&pulse, || false).unwrap();
    assert!(
        tail[95].pressure_pa.abs() > 1e-7,
        "physical edge tail was dropped"
    );
    assert!(
        error(
            values(tail[95]),
            values(owner.observe_at(&pulse, 95, || false).unwrap())
        ) < 2e-10
    );
    let mut opposite = pulse.clone();
    for sample in &mut opposite {
        sample.strength.volume_flow_m3_s *= -2.0;
        sample.strength.volume_flow_rate_m3_s2 *= -2.0;
    }
    let doubled = owner.observe_record(&opposite, || false).unwrap();
    for i in [30, 60, 95] {
        assert!(error(values(doubled[i]), values(tail[i]).map(|v| -2.0 * v)) < 2e-10);
        assert!(
            (doubled[i].perturbation_energy_j_m3 - 4.0 * tail[i].perturbation_energy_j_m3).abs()
                < 1e-10 * tail[i].perturbation_energy_j_m3
        );
    }
}

#[test]
fn declared_history_work_storage_clock_and_source_boundaries_refuse_transactionally() {
    let s = [-0.04, 0.03, 0.01];
    let x = HalfPlaneReceiver::WorldPosition([0.03, -0.02, -0.01]);
    let dt = 0.004;
    let n = 24;
    let owner = prepare(s, x, IdealPlaneBoundary::Rigid, dt, n);
    let source = history(dt, n, 0.0, 0.0, 3);
    let expected = owner.observe_record(&source, || false).unwrap();
    assert!(owner.observe_record(&[], || false).is_err());
    assert!(
        owner
            .observe_record(&history(dt, n + 1, 0.0, 0.0, 3), || false)
            .is_err()
    );
    assert!(owner.observe_at(&source, n, || false).is_err());
    for choice in 0..5 {
        let mut bad = source.clone();
        match choice {
            0 => bad[0].strength.volume_flow_m3_s = 1.0,
            1 => bad[10].source_time_s += dt * 0.1,
            2 => bad[10].strength.force_on_fluid_n[1] = 1.0,
            3 => bad[10].strength.volume_flow_rate_m3_s2 = f64::NAN,
            _ => bad[10].strength.force_rate_on_fluid_n_s[2] = f64::NAN,
        }
        assert!(owner.observe_record(&bad, || false).is_err());
    }
    let mut count = 0;
    assert!(matches!(
        owner.observe_record(&source, || {
            count += 1;
            count > 8
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        values(owner.observe_record(&source, || false).unwrap()[23]),
        values(expected[23])
    );
    let make = |budget: HalfPlaneTransientLimits, stop: bool| {
        PreparedHalfPlaneTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            x,
            dt,
            n,
            budget,
            || stop,
        )
    };
    let peak = owner.diagnostics().peak_payload_bytes;
    assert!(
        make(
            HalfPlaneTransientLimits {
                maximum_payload_bytes: peak,
                ..limits()
            },
            false
        )
        .is_ok()
    );
    assert!(
        make(
            HalfPlaneTransientLimits {
                maximum_payload_bytes: peak - 1,
                ..limits()
            },
            false
        )
        .is_err()
    );
    assert!(
        make(
            HalfPlaneTransientLimits {
                maximum_total_evaluations: 14,
                ..limits()
            },
            false
        )
        .is_err()
    );
    assert!(matches!(
        make(limits(), true),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut count = 0;
    assert!(matches!(
        PreparedHalfPlaneTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            x,
            dt,
            n,
            limits(),
            || {
                count += 1;
                count > 100
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedHalfPlaneTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            x,
            dt,
            MAX_SAMPLES + 1,
            limits(),
            || false
        )
        .is_err()
    );
    assert!(
        PreparedHalfPlaneTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            [0.0; 3],
            x,
            dt,
            n,
            limits(),
            || false
        )
        .is_err()
    );
}
