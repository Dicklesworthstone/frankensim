//! G1/G3: physical force, independent real Green derivatives and causal impulse.
use super::*;
use crate::curle3d::moving::CompactSourceStrength;
use crate::diffraction::HalfPlaneFace;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.7, 1.3).unwrap()
}
fn frame() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}
fn limits() -> HalfPlaneForceTransientLimits {
    HalfPlaneForceTransientLimits {
        quadrature: HalfPlaneForceLimits {
            absolute_pressure_tolerance_pa_per_n: 2e-8,
            absolute_velocity_tolerance_m_s_per_n: 2e-8,
            relative_tolerance: 3e-10,
            maximum_evaluations: 16000,
            maximum_depth: 44,
        },
        maximum_total_evaluations: 4_000_000,
        maximum_payload_bytes: 128_000_000,
    }
}
fn prepare(
    s: [f64; 3],
    x: HalfPlaneReceiver,
    b: IdealPlaneBoundary,
    dt: f64,
    n: usize,
) -> PreparedHalfPlaneForceTransient {
    PreparedHalfPlaneForceTransient::new(medium(), frame(), b, s, x, dt, n, limits(), || false)
        .unwrap()
}
fn record(
    dt: f64,
    n: usize,
    epoch: f64,
    f: impl Fn(f64) -> ([f64; 3], [f64; 3]),
) -> Vec<TimedCompactSource> {
    (0..n)
        .map(|i| {
            let (force, rate) = f(i as f64 * dt);
            TimedCompactSource {
                source_time_s: (Dd::from_f64(epoch) + Dd::from_f64(dt) * Dd::from_f64(i as f64))
                    .to_f64(),
                strength: CompactSourceStrength {
                    force_on_fluid_n: force,
                    force_rate_on_fluid_n_s: rate,
                    ..CompactSourceStrength::default()
                },
            }
        })
        .collect()
}
fn polynomial(
    dt: f64,
    n: usize,
    degree: i32,
    direction: [f64; 3],
    start: f64,
) -> Vec<TimedCompactSource> {
    record(dt, n, 0.0, |t| {
        let a = (t - start).max(0.0);
        (
            direction.map(|d| d * a.powi(degree)),
            direction.map(|d| d * degree as f64 * a.powi(degree - 1)),
        )
    })
}
fn values(x: HalfPlaneTransientField) -> [f64; 4] {
    [
        x.pressure_pa,
        x.particle_velocity_m_s[0],
        x.particle_velocity_m_s[1],
        x.particle_velocity_m_s[2],
    ]
}
fn error(a: [f64; 4], b: [f64; 4]) -> f64 {
    let p = (a[0] - b[0]).abs() / b[0].abs().max(1e-11);
    let u = a[1..]
        .iter()
        .zip(&b[1..])
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
        / b[1..].iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-11);
    p.max(u)
}
fn pulse(t: f64, span: f64) -> [f64; 3] {
    if t <= 0.0 {
        return [0.0; 3];
    }
    if t >= 2.0 * span {
        return [0.0, 0.0, span];
    }
    let u = if t <= span {
        t / span
    } else {
        (2.0 * span - t) / span
    };
    let f = u * u * (3.0 - 2.0 * u);
    let fd = 6.0 * u * (1.0 - u) / span * if t <= span { 1.0 } else { -1.0 };
    let primitive = span * u * u * u * (1.0 - 0.5 * u);
    [
        f,
        fd,
        if t <= span {
            primitive
        } else {
            span - primitive
        },
    ]
}
fn free(
    s: [f64; 3],
    x: [f64; 3],
    time: f64,
    direction: [f64; 3],
    history: impl Fn(f64) -> [f64; 3],
) -> [f64; 4] {
    let delta = core::array::from_fn::<_, 3, _>(|i| x[i] - s[i]);
    let r = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
    let n = delta.map(|v| v / r);
    let c = medium().sound_speed_m_s();
    let [f, fd, impulse] = history(time - r / c);
    let dot = (0..3).map(|i| direction[i] * n[i]).sum::<f64>();
    let mut out = [0.0; 4];
    out[0] = dot * (f / r + fd / c) / (4.0 * PI * r);
    for i in 0..3 {
        out[i + 1] = ((3.0 * n[i] * dot - direction[i]) * (impulse / r + f / c)
            + n[i] * dot * fd * r / c / c)
            / (4.0 * PI * medium().density_kg_m3() * r * r);
    }
    out
}

#[test]
fn original_force_impulse_and_fft_recover_independent_free_near_field() {
    let s = [-0.04, 0.0, 0.012];
    let x = [0.024, 0.038, -0.02];
    let dt = 0.003;
    let n = 161;
    let direction = [1.0, 0.0, -0.35];
    let span = 0.06;
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let input = record(dt, n, 12000.0, |t| {
        let [f, fd, _] = pulse(t, span);
        (direction.map(|d| d * f), direction.map(|d| d * fd))
    });
    let all = owner.observe_record(&input, || false).unwrap();
    let reference: Vec<_> = (0..n)
        .map(|i| free(s, x, i as f64 * dt, direction, |t| pulse(t, span)))
        .collect();
    let scales = core::array::from_fn::<_, 4, _>(|j| {
        reference.iter().map(|v| v[j].abs()).fold(0.0_f64, f64::max)
    });
    let mut worst = 0.0_f64;
    let mut fft_worst = 0.0_f64;
    for i in 0..n {
        let expected = reference[i];
        let actual = values(owner.observe_at(&input, i, || false).unwrap());
        if expected == [0.0; 4] {
            assert_eq!(values(all[i]), [0.0; 4]);
            assert_eq!(actual, [0.0; 4]);
        } else {
            assert!(
                error(actual, expected) < 2e-7,
                "direct {i}: {actual:?} != {expected:?}"
            );
            // FFT roundoff follows the complete signal scale, including at
            // exact pressure zeros after the force pulse. No relative null test.
            for j in 0..4 {
                let normalized = (values(all[i])[j] - expected[j]).abs() / scales[j];
                fft_worst = fft_worst.max(normalized);
                assert!(normalized < 2e-11, "FFT {i}/{j}: {normalized:e}");
            }
            worst = worst.max(error(actual, expected));
        }
    }
    assert_eq!(input[n - 1].strength.force_on_fluid_n, [0.0; 3]);
    assert!(
        all[n - 1]
            .particle_velocity_m_s
            .iter()
            .any(|v| v.abs() > 0.01)
    );
    assert!(all[n - 1].pressure_pa.abs() < 2e-8);
    eprintln!(
        "causal force independent impulse/free-field worst {worst:e}; FFT record-scaled {fft_worst:e}"
    );
}

// Independent UNINTEGRATED Biot--Tolstoy kernel, polar geometry and composite
// Simpson quadrature. This returns grad_receiver G[t_+^degree], using no taps,
// production geometry, transformed quadrature or numerical source differences.
fn original_gradient(
    s: [f64; 3],
    x: [f64; 3],
    time: f64,
    degree: i32,
    b: IdealPlaneBoundary,
) -> [f64; 3] {
    let r = x[0].hypot(x[1]);
    let rs = s[0].hypot(s[1]);
    let l = (r + rs).hypot(x[2] - s[2]);
    let theta = x[1].atan2(x[0]).rem_euclid(2.0 * PI);
    let theta_s = s[1].atan2(s[0]).rem_euclid(2.0 * PI);
    let radial = [x[0] / r, x[1] / r, 0.0];
    let tangent = [-radial[1], radial[0], 0.0];
    let dl = [
        (r + rs) / l * radial[0],
        (r + rs) / l * radial[1],
        (x[2] - s[2]) / l,
    ];
    let c = medium().sound_speed_m_s();
    let mut out = [0.0; 3];
    for term in 0..2 {
        let angle = if term == 0 {
            theta - theta_s
        } else {
            theta + theta_s
        };
        let q = 2.0 * (r * rs).sqrt() * (0.5 * angle).cos();
        let dq = core::array::from_fn::<_, 3, _>(|i| {
            (rs / r).sqrt() * ((0.5 * angle).cos() * radial[i] - (0.5 * angle).sin() * tangent[i])
        });
        let sign = if term == 0 { 1.0 } else { b.sign() };
        if q > 0.0 {
            let image = [s[0], if term == 0 { s[1] } else { -s[1] }, s[2]];
            let delta = core::array::from_fn::<_, 3, _>(|i| x[i] - image[i]);
            let d = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
            let t = (time - d / c).max(0.0);
            let derivative =
                -(t.powi(degree) + d / c * degree as f64 * t.powi(degree - 1)) / (4.0 * PI * d * d);
            for i in 0..3 {
                out[i] += sign * derivative * delta[i] / d;
            }
        }
        if c * time <= l {
            continue;
        }
        let upper = ((c * time - l) * (c * time + l)).sqrt();
        let count = 16384;
        let step = upper / count as f64;
        let mut sum = [0.0; 3];
        for j in 0..=count {
            let t = j as f64 * step;
            let radius = l.hypot(t);
            let age = (time - radius / c).max(0.0);
            let f = age.powi(degree);
            let first = degree as f64 * age.powi(degree - 1);
            let den = t * t + q * q;
            let derivative = -(f + radius / c * first) / (radius * radius);
            let weight = if j == 0 || j == count {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for i in 0..3 {
                sum[i] += weight
                    * (q / den * derivative * l / radius * dl[i]
                        + (t * t - q * q) / (den * den) * f / radius * dq[i]);
            }
        }
        for i in 0..3 {
            out[i] -= sign * step / 3.0 * sum[i] / (4.0 * PI * PI);
        }
    }
    out
}
fn reference(
    s: [f64; 3],
    x: [f64; 3],
    time: f64,
    degree: i32,
    direction: [f64; 3],
    b: IdealPlaneBoundary,
) -> [f64; 4] {
    let pressure = |point, order| {
        let g = original_gradient(point, s, time, order, b);
        (0..3).map(|j| direction[j] * g[j]).sum::<f64>()
    };
    let mut out = [0.0; 4];
    out[0] = pressure(x, degree);
    let h = 0.00002;
    for axis in 0..3 {
        let at = |multiple: f64| {
            let mut point = x;
            point[axis] += multiple * h;
            pressure(point, degree + 1)
        };
        out[axis + 1] = -(at(-2.0) - 8.0 * at(-1.0) + 8.0 * at(1.0) - at(2.0))
            / (12.0 * h * (degree + 1) as f64 * medium().density_kg_m3());
    }
    out
}

#[test]
fn all_original_force_directions_match_independent_real_green_derivatives() {
    let s = [-0.04, 0.03, 0.01];
    let dt = 0.004;
    let n = 81;
    let direction = [0.7, -0.4, 0.3];
    let input = polynomial(dt, n, 3, direction, 0.0);
    for point in [[0.03, -0.02, -0.01], [-0.025, -0.07, 0.03]] {
        for b in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let owner = prepare(s, HalfPlaneReceiver::WorldPosition(point), b, dt, n);
            for i in [40, 80] {
                let actual = values(owner.observe_at(&input, i, || false).unwrap());
                let expected = reference(s, point, i as f64 * dt, 3, direction, b);
                let e = error(actual, expected);
                eprintln!("force raw BTT {b:?} {point:?} t={i}: {e:e}");
                assert!(e < 2e-7, "{actual:?} != {expected:?}");
            }
        }
    }
}

#[test]
fn temporal_force_refinement_converges_to_original_continuous_loading() {
    let s = [-0.04, 0.03, 0.01];
    let x = [0.03, -0.02, -0.01];
    let direction = [0.7, -0.4, 0.3];
    let expected = reference(s, x, 0.32, 4, direction, IdealPlaneBoundary::Rigid);
    let mut errors = Vec::new();
    for dt in [0.02_f64, 0.01, 0.005] {
        let n = (0.32 / dt).round() as usize + 1;
        let owner = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(x),
            IdealPlaneBoundary::Rigid,
            dt,
            n,
        );
        errors.push(error(
            values(
                owner
                    .observe_at(&polynomial(dt, n, 4, direction, 0.0), n - 1, || false)
                    .unwrap(),
            ),
            expected,
        ));
    }
    eprintln!("causal force temporal refinement {errors:?}");
    assert!(errors[1] < 0.5 * errors[0] && errors[2] < 0.5 * errors[1]);
    assert!(errors[2] < 2e-5);
}

#[test]
fn physical_force_sheet_traces_and_narrow_shadow_hessian_are_continuous() {
    let s = [0.03, 0.04, 0.01];
    let dt = 0.004;
    let n = 81;
    let input = polynomial(dt, n, 3, [0.7, -0.4, 0.3], 0.0);
    for b in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
            let owner = prepare(
                s,
                HalfPlaneReceiver::SheetFace {
                    distance_from_edge_m: 0.06,
                    along_edge_m: -0.02,
                    face,
                },
                b,
                dt,
                n,
            );
            let field = owner.observe_at(&input, n - 1, || false).unwrap();
            if b == IdealPlaneBoundary::Rigid {
                assert!(field.particle_velocity_m_s[1].abs() < 2e-10);
            } else {
                assert!(field.pressure_pa.abs() < 2e-10);
            }
        }
        let x = [-0.03, 0.04, -0.02];
        let exact = values(
            prepare(s, HalfPlaneReceiver::WorldPosition(x), b, dt, n)
                .observe_at(&input, n - 1, || false)
                .unwrap(),
        );
        for perturbation in [-1e-9, -1e-12, 1e-12, 1e-9] {
            let mut y = x;
            y[1] += perturbation;
            let actual = values(
                prepare(s, HalfPlaneReceiver::WorldPosition(y), b, dt, n)
                    .observe_at(&input, n - 1, || false)
                    .unwrap(),
            );
            assert!(
                error(actual, exact) < 2e-6,
                "shadow {b:?} {perturbation:e}: {actual:?} != {exact:?}"
            );
        }
    }
}

#[test]
fn original_force_world_axes_clock_and_coherent_work_are_retained() {
    let s = [-0.04, 0.03, 0.01];
    let x = [0.03, -0.02, -0.01];
    let dt = 0.004;
    let n = 81;
    let direction = [0.7, -0.4, 0.3];
    let input = polynomial(dt, n, 3, direction, 0.0);
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let a = owner.observe_at(&input, n - 1, || false).unwrap();
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let origin = [1.0, -0.5, 0.25];
    let place = |v: [f64; 3]| {
        let v = rotate(v);
        core::array::from_fn(|i| v[i] + origin[i])
    };
    let placed = PreparedHalfPlaneForceTransient::new(
        medium(),
        HalfPlaneFrame::new(origin, rotate([0.0, 0.0, 1.0]), rotate([0.0, 1.0, 0.0])).unwrap(),
        IdealPlaneBoundary::Rigid,
        place(s),
        HalfPlaneReceiver::WorldPosition(place(x)),
        dt,
        n,
        limits(),
        || false,
    )
    .unwrap();
    let mut rotated = input.clone();
    for sample in &mut rotated {
        sample.strength.force_on_fluid_n = rotate(sample.strength.force_on_fluid_n);
        sample.strength.force_rate_on_fluid_n_s = rotate(sample.strength.force_rate_on_fluid_n_s);
    }
    let b = placed.observe_at(&rotated, n - 1, || false).unwrap();
    let v = rotate(a.particle_velocity_m_s);
    assert!(error(values(b), [a.pressure_pa, v[0], v[1], v[2]]) < 2e-9);
    for i in 0..3 {
        assert_eq!(
            a.instantaneous_intensity_w_m2[i],
            a.pressure_pa * a.particle_velocity_m_s[i]
        );
    }
    assert_eq!(
        a.instantaneous_energy_flux_w_m2,
        a.instantaneous_intensity_w_m2
    );
    let expected = a.pressure_pa * a.pressure_pa / (2.0 * medium().adiabatic_bulk_modulus_pa())
        + 0.5
            * medium().density_kg_m3()
            * a.particle_velocity_m_s.iter().map(|v| v * v).sum::<f64>();
    assert!((a.perturbation_energy_j_m3 - expected).abs() < 2e-15 * expected);
}

#[test]
fn force_fields_satisfy_original_momentum_continuity_and_energy_transport() {
    let s = [0.03, 0.04, -0.01];
    let x = [-0.02, -0.06, 0.04];
    let dt = 0.002;
    let n = 161;
    let direction = [0.7, -0.4, 0.3];
    let h = 0.00001;
    let start = 0.02;
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let at = |shift| {
        owner
            .observe_at(
                &polynomial(dt, n, 3, direction, start - shift),
                n - 1,
                || false,
            )
            .unwrap()
    };
    let plus = at(h);
    let minus = at(-h);
    let center = at(0.0);
    let dp = (plus.pressure_pa - minus.pressure_pa) / (2.0 * h);
    let du = core::array::from_fn::<_, 3, _>(|i| {
        (plus.particle_velocity_m_s[i] - minus.particle_velocity_m_s[i]) / (2.0 * h)
    });
    let de = (plus.perturbation_energy_j_m3 - minus.perturbation_energy_j_m3) / (2.0 * h);
    let mut grad = [0.0; 3];
    let mut div = 0.0;
    let mut flux = 0.0;
    let input = polynomial(dt, n, 3, direction, start);
    for axis in 0..3 {
        let mut p = x;
        let mut m = x;
        p[axis] += h;
        m[axis] -= h;
        let a = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(p),
            IdealPlaneBoundary::Rigid,
            dt,
            n,
        )
        .observe_at(&input, n - 1, || false)
        .unwrap();
        let b = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(m),
            IdealPlaneBoundary::Rigid,
            dt,
            n,
        )
        .observe_at(&input, n - 1, || false)
        .unwrap();
        grad[axis] = (a.pressure_pa - b.pressure_pa) / (2.0 * h);
        div += (a.particle_velocity_m_s[axis] - b.particle_velocity_m_s[axis]) / (2.0 * h);
        flux += (a.instantaneous_energy_flux_w_m2[axis] - b.instantaneous_energy_flux_w_m2[axis])
            / (2.0 * h);
    }
    let momentum = (0..3)
        .map(|i| (medium().density_kg_m3() * du[i] + grad[i]).abs())
        .sum::<f64>()
        / (grad.iter().map(|v| v.abs()).sum::<f64>() + 1e-30);
    let continuity = (dp + medium().adiabatic_bulk_modulus_pa() * div).abs() / dp.abs();
    let energy = (de + flux).abs() / de.abs();
    eprintln!("causal force momentum={momentum:e}, continuity={continuity:e}, energy={energy:e}");
    assert!(momentum < 2e-5 && continuity < 2e-5 && energy < 2e-5);
    assert!(center.perturbation_energy_j_m3 > 0.0);
}

#[test]
fn original_force_admission_causality_budgets_and_exact_retry_are_explicit() {
    let s = [-0.04, 0.03, 0.01];
    let x = HalfPlaneReceiver::WorldPosition([0.03, -0.02, -0.01]);
    let dt = 0.004;
    let n = 81;
    let owner = prepare(s, x, IdealPlaneBoundary::Rigid, dt, n);
    let input = polynomial(dt, n, 3, [0.7, -0.4, 0.3], 0.0);
    assert!(owner.observe_record(&[], || false).is_err());
    assert!(owner.observe_at(&input, n, || false).is_err());
    for changed in [0, 1, 2, 3] {
        let mut bad = input.clone();
        match changed {
            0 => bad[0].strength.force_on_fluid_n[0] = 1.0,
            1 => bad[2].strength.volume_flow_m3_s = 1.0,
            2 => bad[2].strength.force_rate_on_fluid_n_s[1] = f64::NAN,
            _ => bad[2].source_time_s += 1e-8,
        }
        assert!(owner.observe_record(&bad, || false).is_err());
    }
    let original = owner.observe_at(&input, 40, || false).unwrap();
    let mut future = input.clone();
    for v in &mut future[41..] {
        v.strength.force_on_fluid_n = [1e6; 3];
    }
    assert_eq!(
        values(original),
        values(owner.observe_at(&future, 40, || false).unwrap())
    );
    let all = owner.observe_record(&input, || false).unwrap();
    for limit in [0, 1, 5, 20] {
        let mut polls = 0;
        assert!(matches!(
            owner.observe_record(&input, || {
                polls += 1;
                polls > limit
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    assert_eq!(
        all.iter().map(|f| values(*f)).collect::<Vec<_>>(),
        owner
            .observe_record(&input, || false)
            .unwrap()
            .iter()
            .map(|f| values(*f))
            .collect::<Vec<_>>()
    );
    let mut small = limits();
    small.maximum_payload_bytes = 1;
    assert!(
        PreparedHalfPlaneForceTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            x,
            dt,
            n,
            small,
            || false
        )
        .is_err()
    );
    let mut tiny = limits();
    tiny.maximum_total_evaluations = 1;
    assert!(
        PreparedHalfPlaneForceTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            x,
            dt,
            n,
            tiny,
            || false
        )
        .is_err()
    );
    assert!(
        PreparedHalfPlaneForceTransient::new(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            x,
            dt,
            n,
            limits(),
            || true
        )
        .is_err()
    );
}
