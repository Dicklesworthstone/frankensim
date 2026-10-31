//! G0/G1/G3: coherent line sources, spatial phase, fluid momentum and energy.
use fs_aeroac::SCOPE_STATEMENT;
use fs_aeroac::curle2d::{
    CoherenceBudget, DipoleFieldStatistics, PreparedDipoleObserver, dipole_pressure,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, det};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 343.0).unwrap()
}
fn scaled(force: [C64; 2], scale: C64) -> [C64; 2] {
    force.map(|v| v * scale)
}
fn close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-12),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn g0_constructive_destructive_and_quadrature_fields_have_cross_terms() {
    let force = [C64::new(1.0, 0.3), C64::new(-0.2, 0.4)];
    let one = PreparedDipoleObserver::new(medium(), 700.0, [2.0, 0.5], &[[0.0; 2]], 1).unwrap();
    let pair = PreparedDipoleObserver::new(medium(), 700.0, [2.0, 0.5], &[[0.0; 2]; 2], 2).unwrap();
    let single = one.evaluate(&[force]).unwrap();
    let constructive = pair.evaluate(&[force, force]).unwrap();
    close(
        constructive.pressure_pa,
        single.pressure_pa.scale(2.0),
        2.0e-15,
    );
    for axis in 0..2 {
        close(
            constructive.particle_velocity_m_s[axis],
            single.particle_velocity_m_s[axis].scale(2.0),
            2.0e-15,
        );
        assert!(
            (constructive.active_intensity_w_m2[axis] / single.active_intensity_w_m2[axis] - 4.0)
                .abs()
                < 3.0e-15
        );
    }
    let destructive = pair
        .evaluate(&[force, scaled(force, C64::from_re(-1.0))])
        .unwrap();
    assert_eq!(destructive.pressure_pa, C64::ZERO);
    assert_eq!(destructive.particle_velocity_m_s, [C64::ZERO; 2]);
    assert_eq!(destructive.active_intensity_w_m2, [0.0; 2]);
    let quadrature = pair
        .evaluate(&[force, scaled(force, C64::new(0.0, 1.0))])
        .unwrap();
    close(
        quadrature.pressure_pa,
        single.pressure_pa * C64::new(1.0, 1.0),
        2.0e-15,
    );
    assert!(
        (quadrature.active_intensity_w_m2[0] / single.active_intensity_w_m2[0] - 2.0).abs()
            < 3.0e-15
    );
    assert_eq!(quadrature.scope, SCOPE_STATEMENT);
}

#[test]
fn g1_pressure_antinode_has_zero_net_energy_flux() {
    // Opposing equal dipoles at +/-a produce equal pressures and opposite
    // velocities at the center. Squaring pressure would report false power.
    let a = 0.3;
    let observer =
        PreparedDipoleObserver::new(medium(), 430.0, [0.0; 2], &[[-a, 0.0], [a, 0.0]], 2).unwrap();
    let field = observer
        .evaluate(&[[C64::ONE, C64::ZERO], [-C64::ONE, C64::ZERO]])
        .unwrap();
    assert!(field.pressure_pa.abs() > 1.0);
    assert_eq!(field.particle_velocity_m_s, [C64::ZERO; 2]);
    assert_eq!(field.active_intensity_w_m2, [0.0; 2]);
    let node = observer.evaluate(&[[C64::ONE, C64::ZERO]; 2]).unwrap();
    assert_eq!(node.pressure_pa, C64::ZERO);
    assert!(node.particle_velocity_m_s[0].abs() > 0.0);
}

#[test]
fn g1_velocity_matches_independent_pressure_gradient_near_and_far() {
    let frequency = 310.0;
    let omega = core::f64::consts::TAU * frequency;
    let k = omega / medium().sound_speed_m_s();
    let force = [C64::new(0.7, -0.8), C64::new(-0.3, 0.6)];
    for observer in [[0.007, 0.013], [0.35, -0.27], [7.0, 4.0]] {
        let prepared =
            PreparedDipoleObserver::new(medium(), frequency, observer, &[[0.0; 2]], 1).unwrap();
        let field = prepared.evaluate(&[force]).unwrap();
        let reference = dipole_pressure(force, k, observer, [0.0; 2])
            .unwrap()
            .pressure;
        close(field.pressure_pa, reference, 3.0e-14);
        for axis in 0..2 {
            let h = 1.0e-6 * det::hypot(observer[0], observer[1]);
            let mut plus = observer;
            let mut minus = observer;
            plus[axis] += h;
            minus[axis] -= h;
            let gradient = (dipole_pressure(force, k, plus, [0.0; 2]).unwrap().pressure
                - dipole_pressure(force, k, minus, [0.0; 2]).unwrap().pressure)
                .scale(0.5 / h);
            let velocity =
                C64::new(gradient.im, -gradient.re).scale(1.0 / (medium().density_kg_m3() * omega));
            close(field.particle_velocity_m_s[axis], velocity, 2.0e-8);
        }
    }
}

#[test]
fn g1_half_wavelength_path_difference_suppresses_far_field() {
    let frequency = 900.0;
    let wavelength = medium().sound_speed_m_s() / frequency;
    let distance = 1000.0;
    let sources = [[-0.25 * wavelength, 0.0], [0.25 * wavelength, 0.0]];
    let pair =
        PreparedDipoleObserver::new(medium(), frequency, [distance, 0.0], &sources, 2).unwrap();
    let first = PreparedDipoleObserver::new(medium(), frequency, [distance, 0.0], &sources[..1], 1)
        .unwrap();
    let force = [C64::ONE, C64::ZERO];
    let one = first.evaluate(&[force]).unwrap().pressure_pa;
    let null = pair.evaluate(&[force, force]).unwrap().pressure_pa;
    // Finite-distance cylindrical spreading prevents an exactly zero null.
    assert!(null.abs() / one.abs() < 1.0e-4);
    let shifted = pair
        .evaluate(&[force, scaled(force, -C64::ONE)])
        .unwrap()
        .pressure_pa;
    assert!((shifted.abs() / one.abs() - 2.0).abs() < 2.0e-4);
}

#[test]
fn g3_global_phase_rigid_motion_and_source_permutation_preserve_field() {
    let sources = [[0.1, 0.2], [-0.3, -0.1], [0.0, 0.4]];
    let forces = [
        [C64::new(1.0, 0.3), C64::new(0.2, -0.5)],
        [C64::new(-0.7, 0.4), C64::new(0.3, 0.8)],
        [C64::ZERO, C64::ONE],
    ];
    let observer = [1.3, -0.5];
    let prepared = PreparedDipoleObserver::new(medium(), 330.0, observer, &sources, 3).unwrap();
    let baseline = prepared.evaluate(&forces).unwrap();
    let phase = C64::new(det::cos(0.7), det::sin(0.7));
    let shifted = prepared
        .evaluate(&forces.map(|force| scaled(force, phase)))
        .unwrap();
    close(shifted.pressure_pa, baseline.pressure_pa * phase, 4.0e-15);
    for axis in 0..2 {
        assert!(
            (shifted.active_intensity_w_m2[axis] - baseline.active_intensity_w_m2[axis]).abs()
                < 1.0e-13 * baseline.active_intensity_w_m2[axis].abs()
        );
    }
    // Quarter turn and translation, including force orientation.
    let move_point = |p: [f64; 2]| [2.0 - p[1], -3.0 + p[0]];
    let rotated = PreparedDipoleObserver::new(
        medium(),
        330.0,
        move_point(observer),
        &sources.map(move_point),
        3,
    )
    .unwrap()
    .evaluate(&forces.map(|f| [-f[1], f[0]]))
    .unwrap();
    close(rotated.pressure_pa, baseline.pressure_pa, 1.0e-13);
    close(
        rotated.particle_velocity_m_s[0],
        -baseline.particle_velocity_m_s[1],
        1.0e-13,
    );
    close(
        rotated.particle_velocity_m_s[1],
        baseline.particle_velocity_m_s[0],
        1.0e-13,
    );
    let permuted = PreparedDipoleObserver::new(
        medium(),
        330.0,
        observer,
        &[sources[2], sources[0], sources[1]],
        3,
    )
    .unwrap()
    .evaluate(&[forces[2], forces[0], forces[1]])
    .unwrap();
    close(permuted.pressure_pa, baseline.pressure_pa, 2.0e-15);
    assert_eq!(
        prepared.evaluate(&forces).unwrap().pressure_pa,
        baseline.pressure_pa
    );
}

#[test]
fn g0_small_residual_survives_strong_cancelling_sources() {
    let prepared =
        PreparedDipoleObserver::new(medium(), 343.0, [1.0, 0.0], &[[0.0; 2]; 3], 3).unwrap();
    let force = [C64::ONE, C64::ZERO];
    let expected = prepared
        .evaluate(&[[C64::ZERO; 2], force, [C64::ZERO; 2]])
        .unwrap();
    let result = prepared
        .evaluate(&[
            scaled(force, C64::from_re(1.0e18)),
            force,
            scaled(force, C64::from_re(-1.0e18)),
        ])
        .unwrap();
    close(result.pressure_pa, expected.pressure_pa, 2.0e-15);
    close(
        result.particle_velocity_m_s[0],
        expected.particle_velocity_m_s[0],
        2.0e-15,
    );
}

#[test]
fn g0_budget_geometry_frequency_and_force_refusals_are_transactional() {
    assert!(PreparedDipoleObserver::new(medium(), 100.0, [1.0, 0.0], &[[0.0; 2]; 2], 1).is_err());
    assert!(PreparedDipoleObserver::new(medium(), 100.0, [0.0; 2], &[[0.0; 2]], 1).is_err());
    assert!(
        PreparedDipoleObserver::new(medium(), 100.0, [f64::MAX, 0.0], &[[-f64::MAX, 0.0]], 1)
            .is_err()
    );
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            PreparedDipoleObserver::new(medium(), frequency, [1.0, 0.0], &[[0.0; 2]], 1).is_err()
        );
    }
    let empty = PreparedDipoleObserver::new(medium(), 100.0, [0.0; 2], &[], 0)
        .unwrap()
        .evaluate(&[])
        .unwrap();
    assert_eq!(empty.pressure_pa, C64::ZERO);
    let prepared =
        PreparedDipoleObserver::new(medium(), 100.0, [1.0, 0.0], &[[0.0; 2]], 1).unwrap();
    assert!(prepared.evaluate(&[]).is_err());
    assert!(
        prepared
            .evaluate(&[[C64::new(f64::NAN, 0.0), C64::ZERO]])
            .is_err()
    );
    assert!(prepared.evaluate(&[[C64::from_re(f64::MAX); 2]]).is_err());
    assert!(prepared.evaluate(&[[C64::ONE, C64::ZERO]]).is_ok());
}

fn moments_close(actual: &DipoleFieldStatistics, expected: &DipoleFieldStatistics) {
    close(
        C64::from_re(actual.mean_square_pressure_pa2),
        C64::from_re(expected.mean_square_pressure_pa2),
        2.0e-13,
    );
    for axis in 0..2 {
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            C64::new(
                expected.active_intensity_w_m2[axis],
                expected.reactive_intensity_w_m2[axis],
            ),
            2.0e-13,
        );
        for other in 0..2 {
            close(
                actual.velocity_covariance_m2_s2[axis][other],
                expected.velocity_covariance_m2_s2[axis][other],
                2.0e-13,
            );
        }
    }
}

const TWO_MODES: CoherenceBudget = CoherenceBudget {
    max_modes: 2,
    max_source_mode_pairs: 4,
};

#[test]
fn g1_partial_coherence_matches_complex_cross_spectral_propagation() {
    let positions = [[-0.17, 0.07], [0.29, -0.13]];
    let forces = [
        [C64::new(1.0, 0.3), C64::new(-0.2, 0.5)],
        [C64::new(-0.7, 0.4), C64::new(0.3, -0.2)],
    ];
    let observer = [0.61, 0.38];
    let prepared = PreparedDipoleObserver::new(medium(), 530.0, observer, &positions, 2).unwrap();
    let fields = core::array::from_fn::<_, 2, _>(|source| {
        PreparedDipoleObserver::new(medium(), 530.0, observer, &positions[source..source + 1], 1)
            .unwrap()
            .evaluate(&forces[source..source + 1])
            .unwrap()
    });
    // E[xi_1 conj(xi_2)] = gamma. A non-real gamma catches reversed
    // conjugations that real in-phase/antiphase comparisons cannot reveal.
    for gamma in [-C64::ONE, C64::ZERO, C64::ONE, C64::new(0.3, 0.7)] {
        let independent_weight = det::sqrt((1.0 - gamma.norm_sq()).max(0.0));
        let first = [forces[0], scaled(forces[1], gamma.conj())];
        let second = [
            [C64::ZERO; 2],
            scaled(forces[1], C64::from_re(independent_weight)),
        ];
        let statistics = prepared
            .evaluate_coherence_modes(&[&first, &second], TWO_MODES)
            .unwrap();
        let covariance = [[C64::ONE, gamma], [gamma.conj(), C64::ONE]];
        // Independently assemble H C H*: do not use the factorization for the
        // reference. Source positions and all source phases remain explicit.
        let moment = |left: [C64; 2], right: [C64; 2]| {
            let mut result = C64::ZERO;
            for i in 0..2 {
                for j in 0..2 {
                    result = result + left[i] * covariance[i][j] * right[j].conj();
                }
            }
            result.scale(0.5)
        };
        let pressure = fields.each_ref().map(|field| field.pressure_pa);
        close(
            C64::from_re(statistics.mean_square_pressure_pa2),
            moment(pressure, pressure),
            4.0e-14,
        );
        for axis in 0..2 {
            let velocity = fields
                .each_ref()
                .map(|field| field.particle_velocity_m_s[axis]);
            close(
                C64::new(
                    statistics.active_intensity_w_m2[axis],
                    statistics.reactive_intensity_w_m2[axis],
                ),
                moment(pressure, velocity),
                4.0e-14,
            );
            for other in 0..2 {
                close(
                    statistics.velocity_covariance_m2_s2[axis][other],
                    moment(
                        velocity,
                        fields
                            .each_ref()
                            .map(|field| field.particle_velocity_m_s[other]),
                    ),
                    4.0e-14,
                );
            }
        }
    }
}

#[test]
fn g1_independent_sources_add_power_while_coherent_sources_can_cancel() {
    let prepared =
        PreparedDipoleObserver::new(medium(), 700.0, [0.8, 0.3], &[[0.0; 2]; 2], 2).unwrap();
    let force = [C64::new(1.0, 0.2), C64::new(-0.1, 0.3)];
    let forces = [force, scaled(force, -C64::ONE)];
    let field = prepared.evaluate(&forces).unwrap();
    assert_eq!(field.pressure_pa, C64::ZERO);
    let independent = prepared.evaluate_independent_sources(&forces).unwrap();
    let first = [forces[0], [C64::ZERO; 2]];
    let second = [[C64::ZERO; 2], forces[1]];
    let factored = prepared
        .evaluate_coherence_modes(&[&first, &second], TWO_MODES)
        .unwrap();
    moments_close(&independent, &factored);
    let single = prepared.evaluate(&first).unwrap();
    close(
        C64::from_re(independent.mean_square_pressure_pa2),
        C64::from_re(single.pressure_pa.norm_sq()),
        4.0e-15,
    );
    assert!(independent.mean_square_pressure_pa2 > 0.0);
    assert_eq!(independent.scope, SCOPE_STATEMENT);
}

#[test]
fn g1_statistical_standing_wave_separates_pressure_from_flux() {
    let prepared =
        PreparedDipoleObserver::new(medium(), 430.0, [0.0; 2], &[[-0.3, 0.0], [0.3, 0.0]], 2)
            .unwrap();
    let forces = [[C64::ONE, C64::ZERO], [-C64::ONE, C64::ZERO]];
    let coherent = prepared
        .evaluate_coherence_modes(&[&forces], TWO_MODES)
        .unwrap();
    let independent = prepared.evaluate_independent_sources(&forces).unwrap();
    assert!(coherent.mean_square_pressure_pa2 > 0.0);
    assert_eq!(coherent.velocity_covariance_m2_s2, [[C64::ZERO; 2]; 2]);
    assert_eq!(coherent.active_intensity_w_m2, [0.0; 2]);
    assert_eq!(independent.active_intensity_w_m2, [0.0; 2]);
    assert!(independent.velocity_covariance_m2_s2[0][0].re > 0.0);
    close(
        C64::from_re(coherent.mean_square_pressure_pa2),
        C64::from_re(2.0 * independent.mean_square_pressure_pa2),
        4.0e-15,
    );
}

#[test]
fn g3_coherence_factor_basis_and_clock_do_not_change_statistics() {
    let prepared =
        PreparedDipoleObserver::new(medium(), 343.0, [0.4, -0.2], &[[0.0, 0.0], [0.2, 0.1]], 2)
            .unwrap();
    let first = [
        [C64::ONE, C64::new(0.2, -0.7)],
        [C64::new(0.3, 0.4), C64::ZERO],
    ];
    let second = [
        [C64::new(0.3, -0.1), C64::ONE],
        [C64::ZERO, C64::new(-0.7, 0.2)],
    ];
    let baseline = prepared
        .evaluate_coherence_modes(&[&first, &second], TWO_MODES)
        .unwrap();
    // A unitary rotation of the independent factor basis leaves B B* fixed.
    let i = C64::new(0.0, 1.0);
    let factor = 1.0 / det::sqrt(2.0);
    let rotated_first = core::array::from_fn::<_, 2, _>(|source| {
        core::array::from_fn(|axis| (first[source][axis] + i * second[source][axis]).scale(factor))
    });
    let rotated_second = core::array::from_fn::<_, 2, _>(|source| {
        core::array::from_fn(|axis| (i * first[source][axis] + second[source][axis]).scale(factor))
    });
    let rotated = prepared
        .evaluate_coherence_modes(&[&rotated_first, &rotated_second], TWO_MODES)
        .unwrap();
    moments_close(&rotated, &baseline);
    let phase = C64::new(det::cos(0.7), det::sin(0.7));
    let shifted_first = first.map(|force| scaled(force, phase));
    let shifted_second = second.map(|force| scaled(force, phase));
    let shifted = prepared
        .evaluate_coherence_modes(&[&shifted_first, &shifted_second], TWO_MODES)
        .unwrap();
    moments_close(&shifted, &baseline);
}

#[test]
fn g0_statistics_preserve_quiet_coherent_residuals_and_refuse_invalid_factors() {
    let prepared =
        PreparedDipoleObserver::new(medium(), 343.0, [1.0, 0.0], &[[0.0; 2]; 3], 3).unwrap();
    let budget = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 3,
    };
    let force = [C64::ONE, C64::ZERO];
    let quiet = [[C64::ZERO; 2], force, [C64::ZERO; 2]];
    let strong = [
        scaled(force, C64::from_re(1.0e18)),
        force,
        scaled(force, C64::from_re(-1.0e18)),
    ];
    let expected = prepared
        .evaluate_coherence_modes(&[&quiet], budget)
        .unwrap();
    let actual = prepared
        .evaluate_coherence_modes(&[&strong], budget)
        .unwrap();
    moments_close(&actual, &expected);
    assert!(
        prepared
            .evaluate_coherence_modes(&[&quiet, &quiet], budget)
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&quiet],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 2,
                }
            )
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(&[&quiet[..2]], budget)
            .is_err()
    );
    assert!(prepared.evaluate_independent_sources(&quiet[..2]).is_err());
    for bad in [f64::NAN, f64::INFINITY, 1.0e200] {
        let bad_forces = [[C64::from_re(bad); 2]; 3];
        assert!(
            prepared
                .evaluate_coherence_modes(&[&bad_forces], budget)
                .is_err()
        );
        assert!(prepared.evaluate_independent_sources(&bad_forces).is_err());
    }
    let silence = prepared
        .evaluate_coherence_modes(
            &[],
            CoherenceBudget {
                max_modes: 0,
                max_source_mode_pairs: 0,
            },
        )
        .unwrap();
    assert_eq!(silence.mean_square_pressure_pa2, 0.0);
    let empty = PreparedDipoleObserver::new(medium(), 343.0, [0.0; 2], &[], 0).unwrap();
    assert_eq!(
        empty
            .evaluate_independent_sources(&[])
            .unwrap()
            .mean_square_pressure_pa2,
        0.0
    );
    assert!(
        empty
            .evaluate_coherence_modes(
                &[&[]],
                CoherenceBudget {
                    max_modes: 0,
                    max_source_mode_pairs: 0,
                }
            )
            .is_err()
    );
    // A failed observation cannot mutate the prepared operator.
    moments_close(
        &prepared
            .evaluate_coherence_modes(&[&quiet], budget)
            .unwrap(),
        &expected,
    );
}
