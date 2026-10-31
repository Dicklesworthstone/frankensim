//! Independent SI 3D Green-function, momentum, energy and coherence regressions.

use fs_aeroac::curle3d::{
    CoherenceBudget, CompactForceObserver, CompactSourcePhasor, PreparedCompactSourceObserver,
};
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}

fn phase(theta: f64) -> C64 {
    let (sin, cos) = theta.sin_cos();
    C64::new(cos, sin)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-8),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn close_complex(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 5e-12 * expected.abs().max(1e-8),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn observe(
    rho: f64,
    c: f64,
    frequency: f64,
    observer: [f64; 3],
    positions: &[[f64; 3]],
) -> PreparedCompactSourceObserver {
    PreparedCompactSourceObserver::new(
        medium(rho, c),
        frequency,
        observer,
        positions,
        positions.len(),
    )
    .unwrap()
}

fn reference_pressure(
    rho: f64,
    c: f64,
    frequency: f64,
    observer: [f64; 3],
    positions: &[[f64; 3]],
    sources: &[CompactSourcePhasor],
) -> C64 {
    let omega = TAU * frequency;
    let mut result = C64::ZERO;
    for (position, source) in positions.iter().zip(sources) {
        let displacement = core::array::from_fn::<_, 3, _>(|i| observer[i] - position[i]);
        let range = displacement.iter().map(|v| v * v).sum::<f64>().sqrt();
        let radial_force = (0..3).fold(C64::ZERO, |sum, axis| {
            sum + source.force_on_fluid_n[axis].scale(displacement[axis] / range)
        });
        // Independently assemble the harmonic form of the time-domain owners.
        result = result
            + phase(omega * range / c)
                * (source.volume_flow_m3_s * C64::new(0.0, -rho * omega / (4.0 * PI * range))
                    + radial_force
                        * C64::new(
                            1.0 / (4.0 * PI * range * range),
                            -omega / (4.0 * PI * c * range),
                        ));
    }
    result
}

#[test]
fn harmonic_bridge_matches_retarded_force_at_multiple_times() {
    let rho = 1.3;
    let c = 317.0;
    let frequency = 173.0;
    let omega = TAU * frequency;
    let source = [0.4, -0.3, 0.2];
    let observer = [2.1, 3.2, -0.7];
    let force = [C64::new(2.0, -0.7), C64::new(-0.4, 0.9), C64::new(1.1, 0.3)];
    let retarded = CompactForceObserver::new(medium(rho, c), observer, source).unwrap();
    let prepared = retarded.prepare_harmonic(frequency).unwrap();
    assert_eq!(prepared.observer_position_m(), observer);
    assert_eq!(prepared.source_position_m(0), Some(source));
    assert_eq!(prepared.source_position_m(1), None);
    assert_eq!(prepared.source_delay_s(0), Some(retarded.delay_s()));
    let field = prepared
        .evaluate(&[CompactSourcePhasor::force_dipole(force)])
        .unwrap();
    for time in [0.0, 0.001, 0.023, 0.1] {
        let source_phase = phase(-omega * (time - retarded.delay_s()));
        let f = force.map(|value| (value * source_phase).re);
        let df = force.map(|value| (value * C64::new(0.0, -omega) * source_phase).re);
        close(
            (field.pressure_pa * phase(-omega * time)).re,
            retarded.pressure_at_retarded(f, df).unwrap(),
        );
    }
}

#[test]
fn harmonic_monopole_matches_the_causal_volume_acceleration_owner() {
    let frequency = 480.0;
    let omega = TAU * frequency;
    let flow = C64::new(1e-5, -2e-5);
    let distance = 2.4;
    let acoustic = medium(997.0, 1490.0);
    let retarded = CompactMonopoleObserver::new(
        acoustic,
        distance,
        1000.0,
        MonopoleLimits {
            maximum_ka: 0.1,
            minimum_distance_radii: 10.0,
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap();
    let prepared = PreparedCompactSourceObserver::new(
        acoustic,
        frequency,
        [0.0, 0.0, distance],
        &[[0.0; 3]],
        1,
    )
    .unwrap();
    let field = prepared
        .evaluate(&[CompactSourcePhasor::monopole(flow)])
        .unwrap();
    for time in [0.0, 0.0021, 0.0047] {
        let acceleration = (flow * C64::new(0.0, -omega) * phase(-omega * time)).re;
        let sample = retarded.sample(time, 0.001, acceleration).unwrap();
        close(
            (field.pressure_pa * phase(-omega * sample.arrival_time_s)).re,
            sample.pressure_pa,
        );
    }
}

#[test]
fn monopole_spherical_power_is_radius_invariant_with_reactive_near_field() {
    let rho = 1.25;
    let c = 340.0;
    let frequency = 1200.0;
    let omega = TAU * frequency;
    let flow = C64::new(1e-4, -2e-4);
    let power = rho * omega * omega * flow.norm_sq() / (8.0 * PI * c);
    for radius in [0.02, 0.7, 13.0] {
        let field = observe(rho, c, frequency, [0.0, radius, 0.0], &[[0.0; 3]])
            .evaluate(&[CompactSourcePhasor::monopole(flow)])
            .unwrap();
        close(
            field.active_intensity_w_m2[1] * 4.0 * PI * radius * radius,
            power,
        );
        close(
            field.reactive_intensity_w_m2[1],
            -field.active_intensity_w_m2[1] * c / (omega * radius),
        );
        assert_eq!(field.particle_velocity_m_s[0], C64::ZERO);
        assert_eq!(field.particle_velocity_m_s[2], C64::ZERO);
    }
}

#[test]
fn dipole_has_tangential_velocity_at_its_pressure_null_and_correct_total_power() {
    let rho = 1.2;
    let c = 300.0;
    let frequency = 170.0;
    let force = [C64::new(2.0, 0.7), C64::ZERO, C64::ZERO];
    let null = observe(rho, c, frequency, [0.0, 2.0, 0.0], &[[0.0; 3]])
        .evaluate(&[CompactSourcePhasor::force_dipole(force)])
        .unwrap();
    assert_eq!(null.pressure_pa, C64::ZERO);
    assert!(null.particle_velocity_m_s[0].abs() > 0.0);
    assert_eq!(null.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(null.active_intensity_w_m2, [0.0; 3]);

    let force = [force[0], C64::new(-0.3, 1.1), C64::new(0.8, 0.2)];
    let exact_power = (TAU * frequency).powi(2) * force.iter().map(|f| f.norm_sq()).sum::<f64>()
        / (24.0 * PI * rho * c.powi(3));
    for radius in [0.03, 3.0] {
        // Six-axis spherical quadrature exactly integrates the quadratic dipole pattern.
        let mut flux = 0.0;
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut observer = [0.0; 3];
                observer[axis] = sign * radius;
                let field = observe(rho, c, frequency, observer, &[[0.0; 3]])
                    .evaluate(&[CompactSourcePhasor::force_dipole(force)])
                    .unwrap();
                flux += sign * field.active_intensity_w_m2[axis];
            }
        }
        close(flux * 4.0 * PI * radius * radius / 6.0, exact_power);
    }
}

#[test]
fn all_velocity_components_obey_independently_differentiated_pressure() {
    let positions = [[0.2, -0.1, 0.3], [-0.5, 0.4, -0.2]];
    let sources = [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(3e-4, 2e-4),
            force_on_fluid_n: [C64::new(0.2, 0.4), C64::new(-0.5, 0.1), C64::new(0.8, -0.3)],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-1e-4, 4e-4),
            force_on_fluid_n: [C64::new(0.6, 0.2), C64::new(0.1, -0.4), C64::new(-0.3, 0.2)],
        },
    ];
    let observer = [1.2, 0.8, 1.7];
    let (rho, c, frequency) = (1.4, 329.0, 211.0);
    let field = observe(rho, c, frequency, observer, &positions)
        .evaluate(&sources)
        .unwrap();
    close_complex(
        field.pressure_pa,
        reference_pressure(rho, c, frequency, observer, &positions, &sources),
    );
    for axis in 0..3 {
        let mut plus = observer;
        let mut minus = observer;
        let h = 2e-5;
        plus[axis] += h;
        minus[axis] -= h;
        let gradient = (reference_pressure(rho, c, frequency, plus, &positions, &sources)
            - reference_pressure(rho, c, frequency, minus, &positions, &sources))
        .scale(0.5 / h);
        let expected = gradient / C64::new(0.0, TAU * frequency * rho);
        assert!((field.particle_velocity_m_s[axis] - expected).abs() < 2e-8 * expected.abs());
    }
}

#[test]
fn separated_sources_keep_path_phase_and_standing_waves_keep_zero_flux() {
    let (c, frequency) = (300.0, 300.0);
    let positions = [[0.0; 3], [0.5, 0.0, 0.0]];
    let sources = [
        CompactSourcePhasor::monopole(C64::from_re(1e-3)),
        CompactSourcePhasor::monopole(C64::from_re(0.95e-3)),
    ];
    let separated = observe(1.2, c, frequency, [10.0, 0.0, 0.0], &positions);
    let coherent = separated.evaluate(&sources).unwrap();
    let independent = separated.evaluate_independent_sources(&sources).unwrap();
    assert!(coherent.pressure_pa.norm_sq() < 1e-24 * independent.mean_square_pressure_pa2);
    let collapsed = observe(1.2, c, frequency, [10.0, 0.0, 0.0], &[[0.0; 3]; 2])
        .evaluate(&sources)
        .unwrap();
    assert!(collapsed.pressure_pa.norm_sq() > independent.mean_square_pressure_pa2);

    let standing = observe(
        1.2,
        c,
        frequency,
        [0.0; 3],
        &[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    );
    let q = CompactSourcePhasor::monopole(C64::new(1e-3, 2e-3));
    let same = standing.evaluate(&[q, q]).unwrap();
    assert!(same.pressure_pa.abs() > 0.0);
    assert_eq!(same.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(same.active_intensity_w_m2, [0.0; 3]);
    let opposite = standing
        .evaluate(&[q, CompactSourcePhasor::monopole(-q.volume_flow_m3_s)])
        .unwrap();
    assert_eq!(opposite.pressure_pa, C64::ZERO);
    assert!(opposite.particle_velocity_m_s[0].abs() > 0.0);
    let independent = standing.evaluate_independent_sources(&[q, q]).unwrap();
    close(
        independent.mean_square_pressure_pa2,
        0.25 * same.pressure_pa.norm_sq(),
    );
    assert!(independent.velocity_covariance_m2_s2[0][0].re > 0.0);
}

#[test]
fn volume_and_force_cross_terms_are_retained_even_for_independent_source_vectors() {
    let (rho, c, frequency, radius) = (1.2, 300.0, 240.0, 1.3);
    let omega = TAU * frequency;
    let flow = C64::new(1e-3, 2e-3);
    let force = flow * C64::new(0.0, rho * omega) / C64::new(1.0 / radius, -omega / c);
    let mixed = CompactSourcePhasor {
        volume_flow_m3_s: flow,
        force_on_fluid_n: [force, C64::ZERO, C64::ZERO],
    };
    let prepared = observe(rho, c, frequency, [radius, 0.0, 0.0], &[[0.0; 3]]);
    let field = prepared.evaluate(&[mixed]).unwrap();
    let statistics = prepared.evaluate_independent_sources(&[mixed]).unwrap();
    assert!(field.pressure_pa.abs() < 1e-14);
    assert!(statistics.mean_square_pressure_pa2 < 1e-28);
    assert!(field.particle_velocity_m_s[0].abs() > 0.0);
    let independent_mechanisms = observe(rho, c, frequency, [radius, 0.0, 0.0], &[[0.0; 3]; 2])
        .evaluate_independent_sources(&[
            CompactSourcePhasor::monopole(flow),
            CompactSourcePhasor::force_dipole(mixed.force_on_fluid_n),
        ])
        .unwrap();
    assert!(independent_mechanisms.mean_square_pressure_pa2 > 0.0);
}

#[test]
fn partial_coherence_matches_dense_complex_cross_spectral_algebra() {
    let positions = [[0.0, 0.0, 0.0], [0.3, -0.2, 0.4]];
    let sources = [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(1e-3, 0.5e-3),
            force_on_fluid_n: [C64::from_re(0.3), C64::new(0.2, -0.1), C64::ZERO],
        },
        CompactSourcePhasor::force_dipole([C64::new(-0.2, 0.1), C64::ZERO, C64::new(0.4, 0.3)]),
    ];
    let prepared = observe(1.2, 340.0, 450.0, [1.2, 1.7, -0.8], &positions);
    let silence = CompactSourcePhasor::default();
    let a = prepared.evaluate(&[sources[0], silence]).unwrap();
    let b = prepared.evaluate(&[silence, sources[1]]).unwrap();
    let fields = [
        [
            a.pressure_pa,
            a.particle_velocity_m_s[0],
            a.particle_velocity_m_s[1],
            a.particle_velocity_m_s[2],
        ],
        [
            b.pressure_pa,
            b.particle_velocity_m_s[0],
            b.particle_velocity_m_s[1],
            b.particle_velocity_m_s[2],
        ],
    ];
    let scaled = |source: CompactSourcePhasor, value: C64| CompactSourcePhasor {
        volume_flow_m3_s: source.volume_flow_m3_s * value,
        force_on_fluid_n: source.force_on_fluid_n.map(|f| f * value),
    };
    for gamma in [C64::ONE, -C64::ONE, C64::ZERO, C64::new(0.3, 0.4)] {
        let mode0 = [sources[0], scaled(sources[1], gamma.conj())];
        let mode1 = [
            silence,
            scaled(sources[1], C64::from_re((1.0 - gamma.norm_sq()).sqrt())),
        ];
        let statistics = prepared
            .evaluate_coherence_modes(
                &[&mode0, &mode1],
                CoherenceBudget {
                    max_modes: 2,
                    max_source_mode_pairs: 4,
                },
            )
            .unwrap();
        let moment = |i: usize, j: usize| {
            (fields[0][i] * fields[0][j].conj()
                + fields[1][i] * fields[1][j].conj()
                + fields[0][i] * gamma * fields[1][j].conj()
                + fields[1][i] * gamma.conj() * fields[0][j].conj())
            .scale(0.5)
        };
        close(statistics.mean_square_pressure_pa2, moment(0, 0).re);
        for i in 0..3 {
            close(statistics.active_intensity_w_m2[i], moment(0, i + 1).re);
            close(statistics.reactive_intensity_w_m2[i], moment(0, i + 1).im);
            for j in 0..3 {
                close_complex(
                    statistics.velocity_covariance_m2_s2[i][j],
                    moment(i + 1, j + 1),
                );
                assert_eq!(
                    statistics.velocity_covariance_m2_s2[i][j],
                    statistics.velocity_covariance_m2_s2[j][i].conj()
                );
            }
        }
        if gamma == C64::ZERO {
            let independent = prepared.evaluate_independent_sources(&sources).unwrap();
            close(
                independent.mean_square_pressure_pa2,
                statistics.mean_square_pressure_pa2,
            );
            assert_eq!(
                independent.velocity_covariance_m2_s2,
                statistics.velocity_covariance_m2_s2
            );
        }
    }
}

#[test]
fn prescribed_force_and_volume_flow_have_distinct_material_density_scalings() {
    let observe_source = |rho, source| {
        observe(rho, 300.0, 700.0, [0.8, 1.2, -0.3], &[[0.0; 3]])
            .evaluate(&[source])
            .unwrap()
    };
    for source in [
        CompactSourcePhasor::monopole(C64::new(1e-3, 2e-3)),
        CompactSourcePhasor::force_dipole([C64::from_re(0.5), C64::new(0.3, -0.2), C64::ZERO]),
    ] {
        let base = observe_source(1.0, source);
        for density in [0.1, 10.0, 1000.0] {
            let field = observe_source(density, source);
            let is_volume = source.volume_flow_m3_s != C64::ZERO;
            close_complex(
                field.pressure_pa,
                base.pressure_pa
                    .scale(if is_volume { density } else { 1.0 }),
            );
            for axis in 0..3 {
                close_complex(
                    field.particle_velocity_m_s[axis],
                    base.particle_velocity_m_s[axis].scale(if is_volume {
                        1.0
                    } else {
                        1.0 / density
                    }),
                );
                close(
                    field.active_intensity_w_m2[axis],
                    base.active_intensity_w_m2[axis]
                        * if is_volume { density } else { 1.0 / density },
                );
            }
        }
    }
}

#[test]
fn strong_cancellation_retains_weak_fields_and_common_clock_rotation() {
    let prepared = observe(1.2, 300.0, 500.0, [1.2, 0.7, -0.2], &[[0.0; 3]; 3]);
    let q = |value| CompactSourcePhasor::monopole(value);
    let field = prepared
        .evaluate(&[q(C64::from_re(1e16)), q(C64::ONE), q(C64::from_re(-1e16))])
        .unwrap();
    let weak = prepared
        .evaluate(&[q(C64::ZERO), q(C64::ONE), q(C64::ZERO)])
        .unwrap();
    assert_eq!(field.pressure_pa, weak.pressure_pa);
    assert_eq!(field.particle_velocity_m_s, weak.particle_velocity_m_s);
    let rotated = prepared
        .evaluate(&[q(C64::ZERO), q(phase(0.7)), q(C64::ZERO)])
        .unwrap();
    close_complex(rotated.pressure_pa, weak.pressure_pa * phase(0.7));
    for axis in 0..3 {
        close(
            rotated.active_intensity_w_m2[axis],
            weak.active_intensity_w_m2[axis],
        );
        close(
            rotated.reactive_intensity_w_m2[axis],
            weak.reactive_intensity_w_m2[axis],
        );
    }
}

#[test]
fn source_admission_and_coherence_work_are_explicit_and_bounded() {
    let acoustic = medium(1.2, 300.0);
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedCompactSourceObserver::new(
                acoustic,
                frequency,
                [1.0, 0.0, 0.0],
                &[[0.0; 3]],
                1
            )
            .is_err()
        );
    }
    for position in [[0.0; 3], [f64::NAN, 0.0, 0.0], [1e-200, 0.0, 0.0]] {
        assert!(
            PreparedCompactSourceObserver::new(acoustic, 500.0, position, &[[0.0; 3]], 1).is_err()
        );
    }
    assert!(
        PreparedCompactSourceObserver::new(acoustic, 500.0, [1.0, 0.0, 0.0], &[[0.0; 3]], 0)
            .is_err()
    );
    assert!(
        PreparedCompactSourceObserver::new(acoustic, 1e-300, [1e300, 0.0, 0.0], &[[0.0; 3]], 1)
            .is_err()
    );
    let prepared = observe(1.2, 300.0, 500.0, [1.0, 0.0, 0.0], &[[0.0; 3]; 2]);
    let mode = [CompactSourcePhasor::default(); 2];
    assert!(prepared.evaluate(&mode[..1]).is_err());
    assert!(
        prepared
            .evaluate(&[CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0)); 2])
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&mode, &mode],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 4
                }
            )
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&mode, &mode],
                CoherenceBudget {
                    max_modes: 2,
                    max_source_mode_pairs: 3
                }
            )
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&mode[..1]],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 2
                }
            )
            .is_err()
    );
    let zero = prepared
        .evaluate_coherence_modes(
            &[],
            CoherenceBudget {
                max_modes: 0,
                max_source_mode_pairs: 0,
            },
        )
        .unwrap();
    assert_eq!(zero.mean_square_pressure_pa2, 0.0);
    assert_eq!(zero.velocity_covariance_m2_s2, [[C64::ZERO; 3]; 3]);
    let empty = observe(1.2, 300.0, 500.0, [0.0; 3], &[]);
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
}
