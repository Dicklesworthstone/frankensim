//! G1/G3: analytic axes, original Euler equations and coherent interference.

use super::*;
use core::f64::consts::{PI, TAU};

fn limits() -> ConvectedMonopoleLimits {
    ConvectedMonopoleLimits {
        maximum_mach: 0.95,
        minimum_distance_m: 1e-5,
    }
}
fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.4, 320.0).unwrap()
}
fn phase(angle: f64) -> C64 {
    let (sin, cos) = angle.sin_cos();
    C64::new(cos, sin)
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn observer(
    flow: [f64; 3],
    position: [f64; 3],
    sources: &[[f64; 3]],
) -> PreparedConvectedMonopoleObserver {
    PreparedConvectedMonopoleObserver::new(
        medium(),
        flow,
        83.0,
        position,
        sources,
        sources.len(),
        limits(),
    )
    .unwrap()
}

#[test]
fn zero_flow_recovers_outgoing_spherical_pressure_velocity_and_work() {
    let sources = [[0.1, -0.4, 0.2], [-0.3, 0.2, -0.1]];
    let point = [1.4, -0.7, 0.9];
    let flows = [C64::new(0.003, -0.002), C64::new(-0.001, 0.004)];
    let prepared = observer([0.0; 3], point, &sources);
    let result = prepared.evaluate(&flows).unwrap();
    let omega = TAU * 83.0;
    let mut p = C64::ZERO;
    let mut u = [C64::ZERO; 3];
    for (j, (&source, &q)) in sources.iter().zip(&flows).enumerate() {
        let d = core::array::from_fn::<_, 3, _>(|i| point[i] - source[i]);
        let r = d.iter().map(|v| v * v).sum::<f64>().sqrt();
        let gq = phase(omega * r / 320.0).scale(1.0 / (4.0 * PI * r)) * q;
        p = p + C64::new(0.0, -1.4 * omega) * gq;
        for axis in 0..3 {
            u[axis] = u[axis] + gq * C64::new(1.0 / r, -omega / 320.0).scale(d[axis] / r);
        }
        assert!((prepared.source_delay_s(j).unwrap() - r / 320.0).abs() < 1e-17);
        assert_eq!(prepared.source_position_m(j), Some(source));
    }
    close(result.field.pressure_pa, p, 3e-14);
    for axis in 0..3 {
        close(result.field.particle_velocity_m_s[axis], u[axis], 3e-14);
        let work = (p * u[axis].conj()).re * 0.5;
        assert!((result.mean_perturbation_energy_flux_w_m2[axis] - work).abs() < 1e-16);
        assert_eq!(
            result.mean_perturbation_energy_flux_w_m2[axis],
            result.field.active_intensity_w_m2[axis]
        );
    }
    assert_eq!(prepared.acoustic(), medium());
    assert_eq!(prepared.flow_velocity_m_s(), [0.0; 3]);
    assert_eq!(prepared.observer_position_m(), point);
    assert_eq!(prepared.frequency_hz(), 83.0);
}

#[test]
fn original_galilean_owner_supplies_every_pressure_path_and_delay() {
    let positions = [[0.2, -0.3, 0.1], [-0.4, 0.5, -0.2]];
    let point = [1.4, 0.7, -0.3];
    let flows = [C64::new(0.002, -0.001), C64::new(-0.003, 0.004)];
    for velocity in [[0.0; 3], [100.0, -53.0, 29.0], [-200.0, 31.0, 9.0]] {
        let prepared = observer(velocity, point, &positions);
        let mut expected = ComplexSum::default();
        for (index, (&source, &flow)) in positions.iter().zip(&flows).enumerate() {
            let original = UniformFlowCompactObserver::new(
                medium(),
                velocity,
                point,
                source,
                MovingCompactLimits {
                    maximum_mach: limits().maximum_mach,
                    minimum_distance_m: limits().minimum_distance_m,
                },
            )
            .unwrap();
            assert_eq!(
                prepared.source_delay_s(index),
                Some(original.propagation_delay_s())
            );
            expected
                .add(
                    original
                        .pressure_phasor(83.0, CompactSourcePhasor::monopole(flow))
                        .unwrap(),
                )
                .unwrap();
        }
        close(
            prepared.evaluate(&flows).unwrap().field.pressure_pa,
            expected.value(),
            1e-14,
        );
    }
}

#[test]
fn upstream_downstream_and_transverse_keep_convective_amplitude_and_travel_time() {
    let c = 320.0;
    let m = 0.6;
    let omega = TAU * 83.0;
    let r = 2.3;
    let q = C64::new(0.002, -0.001);
    for direction in [-1.0, 1.0] {
        let prepared = observer([c * m, 0.0, 0.0], [direction * r, 0.0, 0.0], &[[0.0; 3]]);
        let delay = r / (c * (1.0 + direction * m));
        assert!((prepared.source_delay_s(0).unwrap() - delay).abs() < 2e-17);
        let field = prepared.evaluate(&[q]).unwrap().field;
        let gq = q * phase(omega * delay).scale(1.0 / (4.0 * PI * r));
        close(
            field.pressure_pa,
            gq * C64::new(
                -1.4 * c * m * direction / r,
                -1.4 * omega / (1.0 + direction * m),
            ),
            2e-14,
        );
        close(
            field.particle_velocity_m_s[0],
            gq * C64::new(
                direction / r,
                -direction * omega / (c * (1.0 + direction * m)),
            ),
            2e-14,
        );
        assert_eq!(field.particle_velocity_m_s[1..], [C64::ZERO; 2]);
    }
    let transverse = observer([c * m, 0.0, 0.0], [0.0, r, 0.0], &[[0.0; 3]]);
    let beta = (1.0 - m * m).sqrt();
    let delay = r / (c * beta);
    assert!((transverse.source_delay_s(0).unwrap() - delay).abs() < 2e-17);
    let field = transverse.evaluate(&[q]).unwrap();
    let gq = q * phase(omega * delay).scale(1.0 / (4.0 * PI * r * beta));
    close(
        field.field.pressure_pa,
        gq * C64::new(0.0, -1.4 * omega / (beta * beta)),
        3e-14,
    );
    close(
        field.field.particle_velocity_m_s[0],
        gq * C64::new(0.0, omega * m / (c * beta * beta)),
        3e-14,
    );
    assert!(field.mean_perturbation_energy_flux_w_m2[0].abs() > 0.0);
    assert_ne!(
        field.mean_perturbation_energy_flux_w_m2[0],
        field.field.active_intensity_w_m2[0]
    );
}

#[test]
fn actual_linearized_euler_continuity_and_energy_transport_close_away_from_sources() {
    // Independent finite differences of the published p/u and flux, not the
    // analytic kernel's internal gradients. Include two interfering sources.
    let flow = [77.0, -31.0, 42.0];
    let sources = [[-0.2, 0.3, 0.1], [0.4, -0.1, -0.3]];
    let flows = [C64::new(0.007, -0.002), C64::new(-0.004, 0.003)];
    let x = [1.3, 0.7, -0.5];
    let center = observer(flow, x, &sources).evaluate(&flows).unwrap();
    let h = 1e-4;
    let mut pressure_gradient = [C64::ZERO; 3];
    let mut velocity_gradient = [[C64::ZERO; 3]; 3];
    let mut energy_divergence = 0.0;
    for axis in 0..3 {
        let offset = |distance: f64| {
            let mut point = x;
            point[axis] += distance;
            observer(flow, point, &sources).evaluate(&flows).unwrap()
        };
        let values = [offset(-2.0 * h), offset(-h), offset(h), offset(2.0 * h)];
        for (value, weight) in values.iter().zip([1.0, -8.0, 8.0, -1.0]) {
            pressure_gradient[axis] =
                pressure_gradient[axis] + value.field.pressure_pa.scale(weight / (12.0 * h));
            for component in 0..3 {
                velocity_gradient[axis][component] = velocity_gradient[axis][component]
                    + value.field.particle_velocity_m_s[component].scale(weight / (12.0 * h));
            }
            energy_divergence +=
                value.mean_perturbation_energy_flux_w_m2[axis] * weight / (12.0 * h);
        }
    }
    let minus_i_omega = C64::new(0.0, -TAU * 83.0);
    let mut continuity = center.field.pressure_pa * minus_i_omega;
    for axis in 0..3 {
        continuity = continuity
            + pressure_gradient[axis].scale(flow[axis])
            + velocity_gradient[axis][axis].scale(1.4 * 320.0 * 320.0);
    }
    assert!(continuity.abs() < 1e-7 * center.field.pressure_pa.abs());
    for component in 0..3 {
        let mut momentum = center.field.particle_velocity_m_s[component] * minus_i_omega;
        for axis in 0..3 {
            momentum = momentum + velocity_gradient[axis][component].scale(flow[axis]);
        }
        momentum = momentum.scale(1.4) + pressure_gradient[component];
        assert!(momentum.abs() < 2e-10 * center.field.pressure_pa.abs());
    }
    assert!(
        energy_divergence.abs()
            < 2e-9
                * center
                    .mean_perturbation_energy_flux_w_m2
                    .iter()
                    .map(|v| v.abs())
                    .sum::<f64>()
    );
}

#[test]
fn source_interference_and_partial_coherence_survive_before_energy_is_formed() {
    let prepared = observer([110.0, 20.0, -10.0], [1.4, 0.2, 0.1], &[[0.0; 3], [0.0; 3]]);
    let q = C64::new(0.003, -0.001);
    let dark = prepared.evaluate(&[q, -q]).unwrap();
    assert_eq!(dark.field.pressure_pa, C64::ZERO);
    assert_eq!(dark.field.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(dark.mean_perturbation_energy_density_j_m3, 0.0);
    let one = prepared.evaluate(&[q, C64::ZERO]).unwrap();
    let bright = prepared.evaluate(&[q, q]).unwrap();
    assert!(
        (bright.mean_perturbation_energy_density_j_m3 / one.mean_perturbation_energy_density_j_m3
            - 4.0)
            .abs()
            < 1e-14
    );
    let factor_a = [q, q * C64::new(0.3, -0.4)];
    let factor_b = [C64::ZERO, q * C64::from_re(0.75_f64.sqrt())];
    let stats = prepared
        .evaluate_coherence_modes(
            &[&factor_a, &factor_b],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 4,
            },
            || false,
        )
        .unwrap();
    // Explicit four equally likely latent phase choices reproduce this second
    // moment, including the complex correlation, without using factor algebra.
    let a = prepared.evaluate(&factor_a).unwrap();
    let b = prepared.evaluate(&factor_b).unwrap();
    let mut energy = 0.0;
    let mut transport = [0.0; 3];
    for z in [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)] {
        let sample = prepared
            .evaluate(&[factor_a[0] + z * factor_b[0], factor_a[1] + z * factor_b[1]])
            .unwrap();
        energy += sample.mean_perturbation_energy_density_j_m3 * 0.25;
        for axis in 0..3 {
            transport[axis] += sample.mean_perturbation_energy_flux_w_m2[axis] * 0.25;
        }
    }
    assert!((stats.mean_perturbation_energy_density_j_m3 - energy).abs() < 1e-14 * energy);
    assert!(
        (energy
            - a.mean_perturbation_energy_density_j_m3
            - b.mean_perturbation_energy_density_j_m3)
            .abs()
            < 1e-14 * energy
    );
    for axis in 0..3 {
        assert!(
            (stats.mean_perturbation_energy_flux_w_m2[axis] - transport[axis]).abs()
                < 1e-14 * transport[axis].abs()
        );
    }
}

#[test]
fn coherent_material_scaling_and_rigid_coordinate_changes_preserve_the_same_field() {
    let x = [1.7, -0.6, 0.8];
    let sources = [[0.2, -0.1, 0.3], [-0.2, 0.3, -0.4]];
    let u = [80.0, 40.0, -20.0];
    let q = [C64::new(0.004, 0.003), C64::new(-0.001, 0.002)];
    let original = observer(u, x, &sources).evaluate(&q).unwrap();
    let rotate = |v: [f64; 3]| [-v[1], v[2], -v[0]];
    let translate = |v: [f64; 3]| {
        let rotated = rotate(v);
        [rotated[0] + 9.0, rotated[1] - 3.0, rotated[2] + 2.0]
    };
    let changed = observer(rotate(u), translate(x), &sources.map(translate))
        .evaluate(&q)
        .unwrap();
    close(changed.field.pressure_pa, original.field.pressure_pa, 2e-14);
    let velocity = original.field.particle_velocity_m_s;
    for (actual, expected) in changed.field.particle_velocity_m_s.into_iter().zip([
        -velocity[1],
        velocity[2],
        -velocity[0],
    ]) {
        close(actual, expected, 2e-14);
    }
    let flux = rotate(original.mean_perturbation_energy_flux_w_m2);
    for (actual, expected) in changed
        .mean_perturbation_energy_flux_w_m2
        .into_iter()
        .zip(flux)
    {
        assert!((actual - expected).abs() < 3e-14 * expected.abs());
    }
    // At fixed c and U, physical Q fixes velocity; changing rho changes p and
    // energy once. No hard-coded air density or separate c refresh is used.
    let doubled = PreparedConvectedMonopoleObserver::new(
        FluidAcousticState::try_new(2.8, 320.0).unwrap(),
        u,
        83.0,
        x,
        &sources,
        2,
        limits(),
    )
    .unwrap()
    .evaluate(&q)
    .unwrap();
    close(
        doubled.field.pressure_pa,
        original.field.pressure_pa.scale(2.0),
        1e-14,
    );
    for axis in 0..3 {
        close(
            doubled.field.particle_velocity_m_s[axis],
            velocity[axis],
            1e-14,
        );
    }
    assert!(
        (doubled.mean_perturbation_energy_density_j_m3
            / original.mean_perturbation_energy_density_j_m3
            - 2.0)
            .abs()
            < 1e-14
    );
}

#[test]
fn flow_reversal_preserves_green_reciprocity_but_not_ordinary_pressure_reciprocity() {
    let source = [0.2, -0.4, 0.1];
    let receiver = [1.4, 0.3, -0.7];
    let flow = [91.0, -23.0, 31.0];
    let a = observer(flow, receiver, &[source]);
    let b = observer(flow.map(|v| -v), source, &[receiver]);
    assert_eq!(a.source_delay_s(0), b.source_delay_s(0));
    let q = [C64::new(0.002, -0.001)];
    close(
        a.evaluate(&q).unwrap().field.pressure_pa,
        b.evaluate(&q).unwrap().field.pressure_pa,
        2e-14,
    );
    let unflipped = observer(flow, source, &[receiver]);
    assert!((a.source_delay_s(0).unwrap() - unflipped.source_delay_s(0).unwrap()).abs() > 1e-3);
}

#[test]
fn nearly_sonic_admitted_flow_keeps_both_axial_travel_times() {
    let acoustic = FluidAcousticState::try_new(1.0, 1.0).unwrap();
    let mach = 1.0 - 2.0_f64.powi(-40);
    for direction in [-1.0, 1.0] {
        let prepared = PreparedConvectedMonopoleObserver::new(
            acoustic,
            [mach, 0.0, 0.0],
            0.1,
            [direction, 0.0, 0.0],
            &[[0.0; 3]],
            1,
            ConvectedMonopoleLimits {
                maximum_mach: mach,
                minimum_distance_m: 0.1,
            },
        )
        .unwrap();
        let delay = prepared.source_delay_s(0).unwrap();
        let exact = 1.0 / (1.0 + direction * mach);
        assert!((delay / exact - 1.0).abs() < 3e-16);
        assert!(prepared.evaluate(&[C64::from_re(1e-10)]).is_ok());
    }
}

#[test]
fn subsonic_domain_work_and_cancellation_refuse_without_mutating_prepared_paths() {
    let positions = [[0.0; 3], [0.1, 0.0, 0.0]];
    let prepared = observer([64.0, 0.0, 0.0], [1.0, 0.0, 0.0], &positions);
    let q = [C64::ONE, C64::new(-0.2, 0.1)];
    let baseline = prepared.evaluate(&q).unwrap();
    assert!(
        PreparedConvectedMonopoleObserver::new(
            medium(),
            [320.0, 0.0, 0.0],
            83.0,
            [1.0, 0.0, 0.0],
            &positions,
            2,
            limits()
        )
        .is_err()
    );
    assert!(
        PreparedConvectedMonopoleObserver::new(
            medium(),
            [0.0; 3],
            83.0,
            [1.0, 0.0, 0.0],
            &positions,
            1,
            limits()
        )
        .is_err()
    );
    assert!(
        PreparedConvectedMonopoleObserver::new(
            medium(),
            [0.0; 3],
            83.0,
            [0.0; 3],
            &positions,
            2,
            limits()
        )
        .is_err()
    );
    assert!(
        PreparedConvectedMonopoleObserver::new(
            medium(),
            [f64::NAN, 0.0, 0.0],
            83.0,
            [1.0, 0.0, 0.0],
            &positions,
            2,
            limits()
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        prepared.evaluate_cancellable(&q, || {
            polls += 1;
            polls == 3
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(prepared.evaluate(&q[..1]).is_err());
    assert!(
        prepared
            .evaluate(&[C64::new(f64::NAN, 0.0), C64::ZERO])
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&q],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 1
                },
                || false
            )
            .is_err()
    );
    assert!(matches!(
        PreparedConvectedMonopoleObserver::new_cancellable(
            medium(),
            [0.0; 3],
            83.0,
            [1.0, 0.0, 0.0],
            &positions,
            2,
            limits(),
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        prepared.evaluate(&q).unwrap().field.pressure_pa,
        baseline.field.pressure_pa
    );
    let line = prepared.evaluate_harmonic_line(&q).unwrap();
    assert_eq!(line.pressure_pa, baseline.field.pressure_pa);
    assert_eq!(
        line.particle_velocity_m_s,
        baseline.field.particle_velocity_m_s
    );
    let silent = prepared
        .evaluate_coherence_modes(
            &[],
            CoherenceBudget {
                max_modes: 0,
                max_source_mode_pairs: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(silent.mean_perturbation_energy_density_j_m3, 0.0);
    assert_eq!(silent.mean_perturbation_energy_flux_w_m2, [0.0; 3]);
}
