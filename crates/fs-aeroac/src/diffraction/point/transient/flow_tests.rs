//! G1/G3: original-clock convection, independent retarded fields and transport.
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
    source: [f64; 3],
    receiver: HalfPlaneReceiver,
    flow: f64,
    boundary: IdealPlaneBoundary,
    dt: f64,
    n: usize,
) -> PreparedHalfPlaneTransient {
    PreparedHalfPlaneTransient::with_edge_parallel_flow(
        medium(),
        frame(),
        boundary,
        source,
        receiver,
        flow,
        dt,
        n,
        limits(),
        || false,
    )
    .unwrap()
}
fn history(dt: f64, n: usize, start: f64, degree: i32) -> Vec<TimedCompactSource> {
    (0..n)
        .map(|i| {
            let age = i as f64 * dt;
            TimedCompactSource {
                source_time_s: (Dd::from_f64(start) + Dd::from_f64(dt) * Dd::from_f64(i as f64))
                    .to_f64(),
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
fn relative(a: [f64; 4], b: [f64; 4]) -> f64 {
    (a[0] - b[0]).abs() / b[0].abs().max(1e-12)
        + a[1..]
            .iter()
            .zip(&b[1..])
            .map(|(x, y)| (x - y) * (x - y))
            .sum::<f64>()
            .sqrt()
            / b[1..].iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-12)
}

// Independent real retarded free-space potential Q(t-tau)/(4*pi*R_M).
// Its derivatives give the original linearized-Euler pressure and velocity.
fn free(source: [f64; 3], point: [f64; 3], time: f64, flow: f64, degree: i32) -> ([f64; 4], f64) {
    let d: [f64; 3] = core::array::from_fn(|i| point[i] - source[i]);
    let c = medium().sound_speed_m_s();
    let m = flow / c;
    let beta2 = 1.0 - m * m;
    let rm = (beta2 * (d[0] * d[0] + d[1] * d[1]) + d[2] * d[2]).sqrt();
    let distance = d.iter().map(|x| x * x).sum::<f64>().sqrt();
    let delay = if m * d[2] > 0.0 {
        distance / c * (distance / (rm + m * d[2]))
    } else {
        (rm - m * d[2]) / (beta2 * c)
    };
    let age = (time - delay).max(0.0);
    let q = age.powi(degree);
    let dq = degree as f64 * age.powi(degree - 1);
    let dr: [f64; 3] = core::array::from_fn(|i| if i == 2 { d[i] / rm } else { beta2 * d[i] / rm });
    let dtau: [f64; 3] =
        core::array::from_fn(|i| (dr[i] - if i == 2 { m } else { 0.0 }) / (beta2 * c));
    let amplitude = 1.0 / (4.0 * PI * rm);
    let u = core::array::from_fn::<_, 3, _>(|i| amplitude * (dr[i] / rm * q + dtau[i] * dq));
    let p = medium().density_kg_m3()
        * amplitude
        * ((1.0 - flow * dtau[2]) * dq - flow * dr[2] / rm * q);
    ([p, u[0], u[1], u[2]], delay)
}

#[test]
fn zero_mean_flow_retains_original_taps_clock_and_physical_material() {
    let (s, x) = ([0.04, 0.03, 0.02], [0.03, -0.025, -0.01]);
    let (dt, n) = (0.005, 96);
    let old = PreparedHalfPlaneTransient::new(
        medium(),
        frame(),
        IdealPlaneBoundary::Rigid,
        s,
        HalfPlaneReceiver::WorldPosition(x),
        dt,
        n,
        limits(),
        || false,
    )
    .unwrap();
    let zero = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        0.0,
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    assert_eq!(old.coefficients(), zero.coefficients());
    assert_eq!(
        old.diagnostics().evaluations,
        zero.diagnostics().evaluations
    );
    assert_eq!(zero.acoustic(), medium());
    assert_eq!(zero.mean_flow_velocity_m_s(), [0.0; 3]);
    let source = history(dt, n, 2.0, 3);
    for (a, b) in old
        .observe_record(&source, || false)
        .unwrap()
        .into_iter()
        .zip(zero.observe_record(&source, || false).unwrap())
    {
        assert_eq!(values(a), values(b));
        assert_eq!(
            b.instantaneous_intensity_w_m2,
            b.instantaneous_energy_flux_w_m2
        );
    }
}

#[test]
fn convected_causal_cubic_source_recovers_independent_free_pressure_velocity_and_arrival() {
    let (s, dt, n) = ([-0.04, 0.0, 0.01], 0.005, 128);
    for mach in [-0.75, 0.4, 0.85, 0.99] {
        let flow = mach * medium().sound_speed_m_s();
        for x in [[0.024, 0.038, -0.02], [-0.025, -0.03, 0.08]] {
            let owner = prepare(
                s,
                HalfPlaneReceiver::WorldPosition(x),
                flow,
                IdealPlaneBoundary::Rigid,
                dt,
                n,
            );
            assert_eq!(owner.acoustic(), medium());
            assert_eq!(owner.edge_flow_velocity_m_s(), flow);
            let source = history(dt, n, 0.0, 3);
            let batch = owner.observe_record(&source, || false).unwrap();
            let mut direct_error = 0.0_f64;
            let mut batch_error = 0.0_f64;
            for i in 0..n {
                let (expected, delay) = free(s, x, i as f64 * dt, flow, 3);
                let direct = owner.observe_at(&source, i, || false).unwrap();
                if i as f64 * dt < delay {
                    assert_eq!(values(batch[i]), [0.0; 4]);
                    assert_eq!(values(direct), [0.0; 4]);
                } else if i as f64 * dt - delay > 0.01 {
                    direct_error = direct_error.max(relative(values(direct), expected));
                    batch_error = batch_error.max(relative(values(batch[i]), expected));
                }
                let energy = batch[i].pressure_pa.powi(2)
                    / (2.0 * medium().adiabatic_bulk_modulus_pa())
                    + 0.5
                        * medium().density_kg_m3()
                        * batch[i]
                            .particle_velocity_m_s
                            .iter()
                            .map(|v| v * v)
                            .sum::<f64>();
                for axis in 0..3 {
                    let expected_flux = batch[i].instantaneous_intensity_w_m2[axis]
                        + if axis == 2 { flow * energy } else { 0.0 };
                    assert!(
                        (batch[i].instantaneous_energy_flux_w_m2[axis] - expected_flux).abs()
                            < 2e-13 * expected_flux.abs().max(1e-12)
                    );
                }
            }
            eprintln!(
                "causal convected free M={mach}: direct={direct_error:e}, FFT={batch_error:e}"
            );
            assert!(direct_error < 2e-8);
            assert!(batch_error < 2e-7);
        }
    }
}

// Independent UNTRANSFORMED-BY-PARTS real BTT potential. Its source/receiver
// derivatives below are finite differences of the whole potential, not the
// production analytic field map, Hermite taps, quadrature or FFT.
fn reference_potential(
    s: [f64; 3],
    x: [f64; 3],
    time: f64,
    flow: f64,
    boundary: IdealPlaneBoundary,
    degree: i32,
) -> f64 {
    let c = medium().sound_speed_m_s();
    let m = flow / c;
    let beta2 = 1.0 - m * m;
    let beta = beta2.sqrt();
    let dz = x[2] - s[2];
    let r = x[0].hypot(x[1]);
    let rs = s[0].hypot(s[1]);
    let l = (r + rs).hypot(dz / beta);
    let rm = (beta * (r + rs)).hypot(dz);
    let physical = (r + rs).hypot(dz);
    let edge_delay = if m * dz > 0.0 {
        physical / c * (physical / (rm + m * dz))
    } else {
        (rm - m * dz) / (beta2 * c)
    };
    let theta = x[1].atan2(x[0]).rem_euclid(2.0 * PI);
    let theta_s = s[1].atan2(s[0]).rem_euclid(2.0 * PI);
    let mut result = 0.0;
    for j in 0..2 {
        let angle = if j == 0 {
            theta - theta_s
        } else {
            theta + theta_s
        };
        let q = 2.0 * (r * rs).sqrt() * (0.5 * angle).cos();
        let sign = if j == 0 { 1.0 } else { boundary.sign() };
        if q > 0.0 {
            let dy = x[1] - if j == 0 { s[1] } else { -s[1] };
            let radius = (beta2 * ((x[0] - s[0]).powi(2) + dy * dy) + dz * dz).sqrt();
            let delay = (radius - m * dz) / (beta2 * c);
            result += sign * (time - delay).max(0.0).powi(degree) / (4.0 * PI * radius);
        }
        if time <= edge_delay {
            continue;
        }
        let delta = beta * c * (time - edge_delay);
        let upper = delta.sqrt() * (2.0 * l + delta).sqrt();
        let count = 16_384;
        let step = upper / count as f64;
        let mut integral = 0.0;
        for i in 0..=count {
            let w = i as f64 * step;
            let radius = l.hypot(w);
            let extra = (w / (radius + l)) * w / (beta * c);
            let age = (time - edge_delay - extra).max(0.0);
            let coefficient = if i == 0 || i == count {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += coefficient * age.powi(degree) / radius * q / (w * w + q * q);
        }
        result -= sign * integral * step / (12.0 * PI * PI * beta);
    }
    result
}
fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}
fn reference(
    s: [f64; 3],
    x: [f64; 3],
    time: f64,
    flow: f64,
    b: IdealPlaneBoundary,
    degree: i32,
) -> [f64; 4] {
    let dt = derivative(
        |h| reference_potential(s, x, time + h, flow, b, degree),
        2e-5,
    );
    let gradient = core::array::from_fn::<_, 3, _>(|i| {
        derivative(
            |h| {
                let mut p = x;
                p[i] += h;
                reference_potential(s, p, time, flow, b, degree)
            },
            2e-5,
        )
    });
    [
        medium().density_kg_m3() * (dt + flow * gradient[2]),
        -gradient[0],
        -gradient[1],
        -gradient[2],
    ]
}

#[test]
fn convected_edge_tail_matches_independent_real_potential_and_temporal_refinement() {
    let (s, x) = ([0.04, 0.03, 0.02], [0.03, -0.025, -0.01]);
    let time = 0.32;
    for mach in [-0.6, 0.5] {
        let flow = mach * medium().sound_speed_m_s();
        for b in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let expected = reference(s, x, time, flow, b, 3);
            let owner = prepare(s, HalfPlaneReceiver::WorldPosition(x), flow, b, 0.005, 65);
            let actual = owner
                .observe_at(&history(0.005, 65, 0.0, 3), 64, || false)
                .unwrap();
            let error = relative(values(actual), expected);
            eprintln!("causal convected BTT M={mach}, {b:?}: error={error:e}");
            assert!(error < 2e-7);
        }
    }
    let flow = 0.5 * medium().sound_speed_m_s();
    let b = IdealPlaneBoundary::Rigid;
    let expected = reference(s, x, time, flow, b, 4);
    let mut errors = Vec::new();
    for count in [32, 64, 128] {
        let dt = time / count as f64;
        let owner = prepare(
            s,
            HalfPlaneReceiver::WorldPosition(x),
            flow,
            b,
            dt,
            count + 1,
        );
        let actual = owner
            .observe_at(&history(dt, count + 1, 0.0, 4), count, || false)
            .unwrap();
        errors.push(relative(values(actual), expected));
    }
    eprintln!("convected Hermite temporal refinement: {errors:?}");
    assert!(errors[1] < errors[0] && errors[2] < errors[1]);
    assert!(errors[2] < 2e-6);
}

fn observed(
    s: [f64; 3],
    x: [f64; 3],
    time: f64,
    flow: f64,
    b: IdealPlaneBoundary,
) -> HalfPlaneTransientField {
    let count = 64;
    let dt = time / count as f64;
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        flow,
        b,
        dt,
        count + 1,
    );
    owner
        .observe_at(&history(dt, count + 1, 0.0, 3), count, || false)
        .unwrap()
}
#[test]
fn real_convected_fields_close_momentum_continuity_and_advected_energy() {
    let (s, x, time) = ([0.04, 0.03, 0.02], [0.03, -0.025, -0.01], 0.32);
    let (flow, b) = (0.5 * medium().sound_speed_m_s(), IdealPlaneBoundary::Rigid);
    let h = 5e-5;
    let time_samples =
        [-2.0, -1.0, 1.0, 2.0].map(|offset| observed(s, x, time + offset * h, flow, b));
    let stencil = |f: [f64; 4]| (f[0] - 8.0 * f[1] + 8.0 * f[2] - f[3]) / (12.0 * h);
    let time_rate =
        core::array::from_fn::<_, 4, _>(|i| stencil(time_samples.map(|f| values(f)[i])));
    let energy_rate = stencil(time_samples.map(|f| f.perturbation_energy_j_m3));
    let mut gradient = [[0.0; 4]; 3];
    let mut divergence_flux = 0.0;
    for axis in 0..3 {
        let samples = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
            let mut at = x;
            at[axis] += offset * h;
            observed(s, at, time, flow, b)
        });
        gradient[axis] = core::array::from_fn(|i| stencil(samples.map(|f| values(f)[i])));
        divergence_flux += stencil(samples.map(|f| f.instantaneous_energy_flux_w_m2[axis]));
    }
    let rho = medium().density_kg_m3();
    let bulk = medium().adiabatic_bulk_modulus_pa();
    let mut momentum = 0.0_f64;
    for axis in 0..3 {
        let material = rho * (time_rate[axis + 1] + flow * gradient[2][axis + 1]);
        let pressure = gradient[axis][0];
        momentum = momentum
            .max((material + pressure).abs() / (material.abs() + pressure.abs()).max(1e-12));
    }
    let material = time_rate[0] + flow * gradient[2][0];
    let compression = bulk * (gradient[0][1] + gradient[1][2] + gradient[2][3]);
    let continuity =
        (material + compression).abs() / (material.abs() + compression.abs()).max(1e-12);
    let energy = (energy_rate + divergence_flux).abs()
        / (energy_rate.abs() + divergence_flux.abs()).max(1e-12);
    eprintln!(
        "convected time PDE: momentum={momentum:e} continuity={continuity:e} energy={energy:e}"
    );
    assert!(momentum < 3e-6 && continuity < 3e-6 && energy < 3e-6);
}

#[test]
fn both_sheet_faces_and_shadow_transitions_preserve_physical_flow_conditions() {
    let source = [0.04, 0.03, 0.02];
    for flow in [-0.65, 0.78] {
        for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
            for b in [
                IdealPlaneBoundary::Rigid,
                IdealPlaneBoundary::PressureRelease,
            ] {
                let owner = prepare(
                    source,
                    HalfPlaneReceiver::SheetFace {
                        distance_from_edge_m: 0.03,
                        along_edge_m: -0.01,
                        face,
                    },
                    flow,
                    b,
                    0.005,
                    65,
                );
                for field in owner
                    .observe_record(&history(0.005, 65, 0.0, 3), || false)
                    .unwrap()
                {
                    match b {
                        IdealPlaneBoundary::Rigid => {
                            assert_eq!(field.particle_velocity_m_s[1], 0.0)
                        }
                        IdealPlaneBoundary::PressureRelease => assert_eq!(field.pressure_pa, 0.0),
                    }
                }
            }
        }
        // theta_x = theta_s + pi is a direct-shadow boundary.
        let x = [-0.04, -0.03, -0.01];
        let center = observed(source, x, 0.32, flow, IdealPlaneBoundary::Rigid);
        for epsilon in [1e-9, 1e-12] {
            for sign in [-1.0, 1.0] {
                let mut point = x;
                point[1] += sign * epsilon;
                let field = observed(source, point, 0.32, flow, IdealPlaneBoundary::Rigid);
                assert!(relative(values(field), values(center)) < 1e-6);
            }
        }
    }
}

#[test]
fn placed_flow_original_clock_and_post_emission_memory_survive_cached_convolution() {
    let (s, x, flow, dt, n) = ([0.04, 0.03, 0.02], [0.03, -0.025, -0.01], 0.65, 0.005, 161);
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        flow,
        IdealPlaneBoundary::Rigid,
        dt,
        n,
    );
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let origin = [1.0, -2.0, 0.5];
    let placed = |p: [f64; 3]| {
        let r = rotate(p);
        core::array::from_fn(|i| r[i] + origin[i])
    };
    let world_frame =
        HalfPlaneFrame::new(origin, rotate([0.0, 0.0, 1.0]), rotate([0.0, 1.0, 0.0])).unwrap();
    let world = PreparedHalfPlaneTransient::with_edge_parallel_flow(
        medium(),
        world_frame,
        IdealPlaneBoundary::Rigid,
        placed(s),
        HalfPlaneReceiver::WorldPosition(placed(x)),
        flow,
        dt,
        n,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(world.mean_flow_velocity_m_s(), rotate([0.0, 0.0, flow]));
    let mut source = history(dt, n, 2.0, 3);
    let duration = 0.08;
    for (i, sample) in source.iter_mut().enumerate() {
        let time = i as f64 * dt;
        let z = time / duration;
        if z < 1.0 {
            sample.strength.volume_flow_m3_s = time.powi(3) * (1.0 - z).powi(3);
            sample.strength.volume_flow_rate_m3_s2 =
                3.0 * time.powi(2) * (1.0 - z).powi(2) * (1.0 - 2.0 * z);
        } else {
            sample.strength = CompactSourceStrength::default()
        }
    }
    let fields = owner.observe_record(&source, || false).unwrap();
    let rotated = world.observe_record(&source, || false).unwrap();
    for i in [50, 90, 150] {
        let direct = owner.observe_at(&source, i, || false).unwrap();
        assert!(relative(values(fields[i]), values(direct)) < 2e-6);
        let expected = [
            fields[i].pressure_pa,
            fields[i].particle_velocity_m_s[2],
            fields[i].particle_velocity_m_s[0],
            fields[i].particle_velocity_m_s[1],
        ];
        assert!(relative(values(rotated[i]), expected) < 2e-6);
    }
    assert!(
        fields[90].pressure_pa.abs() > 1e-13,
        "physical diffracted tail after source stopped"
    );
    let mut changed = source.clone();
    for sample in &mut changed[100..] {
        sample.strength.volume_flow_m3_s = 1e-3;
    }
    let changed = owner.observe_record(&changed, || false).unwrap();
    for i in 0..90 {
        for (a, b) in values(fields[i]).into_iter().zip(values(changed[i])) {
            assert!((a - b).abs() < 2e-12);
        }
    }
}

#[test]
fn flow_admission_quiescent_clock_and_cancellation_are_transactional() {
    let (s, x) = ([0.04, 0.03, 0.02], [0.03, -0.025, -0.01]);
    for flow in [f64::NAN, f64::INFINITY, -f64::INFINITY, 1.3, -1.3, 1.4] {
        assert!(
            PreparedHalfPlaneTransient::with_edge_parallel_flow(
                medium(),
                frame(),
                IdealPlaneBoundary::Rigid,
                s,
                HalfPlaneReceiver::WorldPosition(x),
                flow,
                0.005,
                65,
                limits(),
                || false
            )
            .is_err()
        );
    }
    let mut count = 0;
    assert!(
        PreparedHalfPlaneTransient::with_edge_parallel_flow(
            medium(),
            frame(),
            IdealPlaneBoundary::Rigid,
            s,
            HalfPlaneReceiver::WorldPosition(x),
            0.65,
            0.005,
            65,
            limits(),
            || {
                count += 1;
                count > 10
            }
        )
        .is_err()
    );
    let owner = prepare(
        s,
        HalfPlaneReceiver::WorldPosition(x),
        0.65,
        IdealPlaneBoundary::Rigid,
        0.005,
        65,
    );
    let source = history(0.005, 65, 2.0, 3);
    let expected = owner.observe_record(&source, || false).unwrap();
    assert!(owner.observe_record(&source, || true).is_err());
    let again = owner.observe_record(&source, || false).unwrap();
    for (a, b) in expected.into_iter().zip(again) {
        assert_eq!(values(a), values(b));
    }
    let mut bad = source.clone();
    bad[1].source_time_s += 1e-6;
    assert!(owner.observe_record(&bad, || false).is_err());
    let mut bad = source.clone();
    bad[8].strength.force_on_fluid_n[0] = 1.0;
    assert!(owner.observe_record(&bad, || false).is_err());
    let mut bad = source.clone();
    bad[0].strength.volume_flow_m3_s = 1.0;
    assert!(owner.observe_record(&bad, || false).is_err());
    // Large original epochs still use the prepared dt for causal arrival.
    let source = [-1.0, 0.0, 0.0];
    let point = [0.0, 4.0, 2.0];
    let flow = 0.5;
    let (_, delay) = free(source, point, 8.0, flow, 3);
    let dt = delay + 0.1;
    let model = prepare(
        source,
        HalfPlaneReceiver::WorldPosition(point),
        flow,
        IdealPlaneBoundary::Rigid,
        dt,
        2,
    );
    let knots = history(dt, 2, (1u64 << 53) as f64, 3);
    let a = model.observe_record(&knots, || false).unwrap();
    let b = model.observe_at(&knots, 1, || false).unwrap();
    assert!(b.pressure_pa != 0.0);
    assert!(relative(values(a[1]), values(b)) < 2e-8);
}
