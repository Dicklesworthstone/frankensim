//! G1/G3/G4: original-source directional interference and outgoing wave momentum.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{
    CompactRadiatedMomentum, CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary,
    PreparedCompactRadiatedMomentum, PreparedCompactRadiatedPower, PreparedCompactSourceObserver,
    RadiatedPowerBudget, RadiatedPowerModeBudget,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 331.0).unwrap()
}

fn budget(n: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget {
        max_sources: n,
        max_source_pairs: n * n,
    }
}

fn mode_budget(n: usize, rank: usize) -> RadiatedPowerModeBudget {
    RadiatedPowerModeBudget {
        max_modes: rank,
        max_source_pair_mode_products: n * n * rank,
    }
}

fn phase(angle: f64) -> C64 {
    let (sine, cosine) = angle.sin_cos();
    C64::new(cosine, sine)
}

fn scale(source: CompactSourcePhasor, factor: C64) -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: source.volume_flow_m3_s * factor,
        force_on_fluid_n: source.force_on_fluid_n.map(|v| v * factor),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn assert_physical(result: &CompactRadiatedMomentum, acoustic: FluidAcousticState) {
    let first = result.power_first_moment_w;
    let norm = first[0].hypot(first[1]).hypot(first[2]);
    assert!(result.power_w >= 0.0);
    assert!(norm <= result.power_w * (1.0 + 2e-13), "{result:?}");
    for axis in 0..3 {
        close(
            result.wave_momentum_flux_n[axis] * acoustic.sound_speed_m_s(),
            first[axis],
            result.power_w * 2e-15,
        );
    }
}

fn cardioid() -> CompactSourcePhasor {
    // Z Q = Fx: the far-pressure pattern is proportional to 1+n_x.
    CompactSourcePhasor {
        volume_flow_m3_s: C64::from_re(1.0 / medium().characteristic_impedance_pa_s_m()),
        force_on_fluid_n: [C64::from_re(1.0), C64::ZERO, C64::ZERO],
    }
}

#[test]
fn colocated_volume_force_phase_controls_direction_without_changing_total_power() {
    let prepared =
        PreparedCompactRadiatedMomentum::new(medium(), 137.0, &[[0.0; 3]], budget(1)).unwrap();
    let base = cardioid();
    let rho_c = medium().characteristic_impedance_pa_s_m();
    let k = TAU * 137.0 / medium().sound_speed_m_s();
    let gain_squared = k * k / (8.0 * PI * rho_c);
    for angle in [0.0, PI / 3.0, PI / 2.0, PI, -PI / 3.0] {
        let mut source = base;
        source.force_on_fluid_n[0] = phase(angle);
        let result = prepared.evaluate(&[source]).unwrap();
        close(result.power_w, (4.0 / 3.0) * gain_squared, gain_squared * 2e-14);
        close(
            result.power_first_moment_w[0],
            (2.0 / 3.0) * gain_squared * angle.cos(),
            gain_squared * 2e-14,
        );
        assert_eq!(result.power_first_moment_w[1..], [0.0; 2]);
        assert_physical(&result, medium());
    }
    for source in [
        CompactSourcePhasor::monopole(base.volume_flow_m3_s),
        CompactSourcePhasor::force_dipole(base.force_on_fluid_n),
    ] {
        let result = prepared.evaluate(&[source]).unwrap();
        assert_eq!(result.power_first_moment_w, [0.0; 3]);
        assert!(result.power_w > 0.0);
    }
}

fn fixture() -> ([[f64; 3]; 3], [CompactSourcePhasor; 3]) {
    let positions = [
        [-0.25, 0.125, 0.0],
        [0.375, -0.25, 0.125],
        [0.0, 0.25, -0.375],
    ];
    let sources = [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(4e-4, 2e-4),
            force_on_fluid_n: [
                C64::new(0.3, 0.1),
                C64::new(-0.2, 0.4),
                C64::new(0.1, -0.3),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-1e-4, 3e-4),
            force_on_fluid_n: [
                C64::new(-0.1, 0.5),
                C64::new(0.2, -0.3),
                C64::new(0.4, 0.1),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(2e-4, -1e-4),
            force_on_fluid_n: [
                C64::new(0.2, -0.1),
                C64::new(-0.3, 0.2),
                C64::new(0.1, 0.3),
            ],
        },
    ];
    (positions, sources)
}

// Independent Gauss-Legendre / periodic-azimuth integration of the actual
// physical source pattern. No production angular moments or field transfers.
fn quadrature(
    hz: f64,
    positions: &[[f64; 3]],
    sources: &[CompactSourcePhasor],
) -> [f64; 4] {
    fn legendre(order: usize, x: f64) -> (f64, f64) {
        let (mut previous, mut value) = (0.0, 1.0);
        for degree in 1..=order {
            let next = ((2 * degree - 1) as f64 * x * value
                - (degree - 1) as f64 * previous)
                / degree as f64;
            previous = value;
            value = next;
        }
        (value, previous)
    }
    let order = 64;
    let azimuths = 128;
    let mut result = [0.0; 4];
    let z = medium().characteristic_impedance_pa_s_m();
    let k = TAU * hz / medium().sound_speed_m_s();
    for index in 0..order {
        let mut mu = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..20 {
            let (p, previous) = legendre(order, mu);
            let derivative = order as f64 * (mu * p - previous) / (mu * mu - 1.0);
            let correction = p / derivative;
            mu -= correction;
            if correction.abs() < 1e-15 {
                break;
            }
        }
        let (p, previous) = legendre(order, mu);
        let derivative = order as f64 * (mu * p - previous) / (mu * mu - 1.0);
        // Weight includes the sphere AVERAGE, rather than a solid-angle integral.
        let weight =
            1.0 / ((1.0 - mu * mu) * derivative * derivative * azimuths as f64);
        let radial = (1.0 - mu * mu).sqrt();
        for azimuth in 0..azimuths {
            let angle = TAU * azimuth as f64 / azimuths as f64;
            let direction = [radial * angle.cos(), radial * angle.sin(), mu];
            let mut amplitude = C64::ZERO;
            for (position, source) in positions.iter().zip(sources) {
                let projected = (0..3).fold(C64::ZERO, |sum, axis| {
                    sum + source.force_on_fluid_n[axis].scale(direction[axis])
                });
                let path_phase = -k
                    * (0..3)
                        .map(|axis| position[axis] * direction[axis])
                        .sum::<f64>();
                amplitude = amplitude
                    + (source.volume_flow_m3_s.scale(z) + projected) * phase(path_phase);
            }
            let power = weight * k * k * amplitude.norm_sq() / (8.0 * PI * z);
            result[0] += power;
            for axis in 0..3 {
                result[axis + 1] += power * direction[axis];
            }
        }
    }
    result
}

#[test]
fn separated_full_vector_interference_matches_independent_angular_quadrature() {
    let (positions, original) = fixture();
    // Pure force isolates the rank-three terms; mixed sources also cover all
    // second-rank Q/F terms. Frequencies straddle the small-x series switch.
    for force_only in [false, true] {
        let sources = original.map(|mut source| {
            if force_only {
                source.volume_flow_m3_s = C64::ZERO;
            }
            source
        });
        for hz in [13.0, 217.0, 1800.0] {
            let prepared =
                PreparedCompactRadiatedMomentum::new(medium(), hz, &positions, budget(3))
                    .unwrap();
            let result = prepared.evaluate(&sources).unwrap();
            let reference = quadrature(hz, &positions, &sources);
            close(result.power_w, reference[0], reference[0] * 3e-12);
            for axis in 0..3 {
                close(
                    result.power_first_moment_w[axis],
                    reference[axis + 1],
                    reference[0] * 3e-12,
                );
            }
            assert_eq!(
                result.power_w,
                prepared.power_operator().evaluate(&sources).unwrap().power_w,
            );
            assert_physical(&result, medium());
        }
    }
}

#[test]
fn quadrature_phase_reversal_steers_separated_monopoles_in_the_expected_direction() {
    let hz = 101.0;
    let d = 0.25;
    let k = TAU * hz / medium().sound_speed_m_s();
    let x = k * d;
    let j1 = (x.sin() / x - x.cos()) / x;
    let q = 0.001;
    let isolated =
        k * k * medium().characteristic_impedance_pa_s_m() * q * q / (8.0 * PI);
    let prepared = PreparedCompactRadiatedMomentum::new(
        medium(),
        hz,
        &[[0.0; 3], [d, 0.0, 0.0]],
        budget(2),
    )
    .unwrap();
    for sign in [-1.0, 1.0] {
        let sources = [
            CompactSourcePhasor::monopole(C64::from_re(q)),
            CompactSourcePhasor::monopole(C64::new(0.0, sign * q)),
        ];
        let result = prepared.evaluate(&sources).unwrap();
        close(result.power_w, 2.0 * isolated, isolated * 2e-14);
        close(
            result.power_first_moment_w[0],
            sign * 2.0 * j1 * isolated,
            isolated * 2e-14,
        );
        assert_eq!(result.power_first_moment_w[1..], [0.0; 2]);
        assert_physical(&result, medium());
    }
}

#[test]
fn weak_interference_survives_colocated_and_sub_ulp_separation_cancellation() {
    let base = cardioid();
    let sources = [
        scale(base, C64::from_re(1e16)),
        base,
        scale(base, C64::from_re(-1e16)),
    ];
    let prepared =
        PreparedCompactRadiatedMomentum::new(medium(), 100.0, &[[0.0; 3]; 3], budget(3))
            .unwrap();
    let single =
        PreparedCompactRadiatedMomentum::new(medium(), 100.0, &[[0.0; 3]], budget(1))
            .unwrap()
            .evaluate(&[base])
            .unwrap();
    assert_eq!(prepared.evaluate(&sources).unwrap(), single);

    let d = 1e-10;
    let hz = medium().sound_speed_m_s() / TAU;
    let k = TAU * hz / medium().sound_speed_m_s();
    let x = k * d;
    let delta = x / 3.0_f64.sqrt();
    let prepared = PreparedCompactRadiatedMomentum::new(
        medium(),
        hz,
        &[[0.0; 3], [d, 0.0, 0.0]],
        budget(2),
    )
    .unwrap();
    let result = prepared
        .evaluate(&[
            CompactSourcePhasor::monopole(C64::from_re(1.0)),
            CompactSourcePhasor::monopole(C64::new(-1.0, -delta)),
        ])
        .unwrap();
    let gain = k * k * medium().characteristic_impedance_pa_s_m() / (8.0 * PI);
    let expected_power = gain * (delta * delta + x * x / 3.0);
    let expected_moment = -gain * 2.0 * x * delta / 3.0;
    close(result.power_w, expected_power, expected_power * 3e-14);
    close(result.power_first_moment_w[0], expected_moment, expected_power * 3e-14);
    assert!(result.power_first_moment_w[0] < 0.0);
    assert_physical(&result, medium());
}

#[test]
fn source_frame_phase_and_original_material_are_preserved() {
    let (positions, sources) = fixture();
    let prepared =
        PreparedCompactRadiatedMomentum::new(medium(), 417.0, &positions, budget(3)).unwrap();
    let base = prepared.evaluate(&sources).unwrap();
    let shifted = positions.map(|p| [p[1] + 8.0, p[2] - 4.0, p[0] + 2.0]);
    let transformed = sources.map(|s| {
        scale(
            CompactSourcePhasor {
                volume_flow_m3_s: s.volume_flow_m3_s,
                force_on_fluid_n: [
                    s.force_on_fluid_n[1],
                    s.force_on_fluid_n[2],
                    s.force_on_fluid_n[0],
                ],
            },
            phase(0.73),
        )
    });
    let rotated =
        PreparedCompactRadiatedMomentum::new(medium(), 417.0, &shifted, budget(3))
            .unwrap()
            .evaluate(&transformed)
            .unwrap();
    close(rotated.power_w, base.power_w, base.power_w * 2e-14);
    for axis in 0..3 {
        close(
            rotated.power_first_moment_w[axis],
            base.power_first_moment_w[(axis + 1) % 3],
            base.power_w * 2e-14,
        );
    }
    let reordered = PreparedCompactRadiatedMomentum::new(
        medium(),
        417.0,
        &[positions[2], positions[0], positions[1]],
        budget(3),
    )
    .unwrap()
    .evaluate(&[sources[2], sources[0], sources[1]])
    .unwrap();
    close(reordered.power_w, base.power_w, base.power_w * 2e-14);
    for axis in 0..3 {
        close(
            reordered.power_first_moment_w[axis],
            base.power_first_moment_w[axis],
            base.power_w * 2e-14,
        );
    }
    let observer =
        PreparedCompactSourceObserver::new(medium(), 417.0, [2.0; 3], &positions, 3)
            .unwrap();
    let from_observer = PreparedCompactRadiatedPower::from_observer(&observer, budget(3))
        .unwrap()
        .with_momentum()
        .unwrap();
    assert_eq!(from_observer.power_operator().acoustic_state(), medium());
    assert_eq!(from_observer.evaluate(&sources).unwrap(), base);
    // A second original medium tests the final 1/c separately from the power
    // normalization; equal Z and k give equal patterns but half the flux.
    let doubled_c = FluidAcousticState::try_new(0.65, 662.0).unwrap();
    let changed =
        PreparedCompactRadiatedMomentum::new(doubled_c, 834.0, &positions, budget(3))
            .unwrap()
            .evaluate(&sources)
            .unwrap();
    assert_eq!(changed.power_w, base.power_w);
    for axis in 0..3 {
        assert_eq!(
            2.0 * changed.wave_momentum_flux_n[axis],
            base.wave_momentum_flux_n[axis],
        );
    }
}

#[test]
fn explicit_coherence_keeps_intrinsic_source_direction_without_inventing_cross_terms() {
    let base = cardioid();
    let prepared =
        PreparedCompactRadiatedMomentum::new(medium(), 170.0, &[[0.0; 3]; 2], budget(2))
            .unwrap();
    let cancel = [base, scale(base, C64::from_re(-1.0))];
    let zero = prepared.evaluate(&cancel).unwrap();
    assert_eq!(zero.power_w, 0.0);
    assert_eq!(zero.wave_momentum_flux_n, [0.0; 3]);
    let independent = prepared.evaluate_independent_sources(&cancel).unwrap();
    assert!(independent.wave_momentum_flux_n[0] > 0.0);
    close(
        independent.power_first_moment_w[0],
        0.5 * independent.power_w,
        independent.power_w * 2e-14,
    );
    let left = [base, CompactSourcePhasor::default()];
    let right = [CompactSourcePhasor::default(), cancel[1]];
    assert_eq!(
        independent,
        prepared
            .evaluate_coherence_modes(&[&left, &right], mode_budget(2, 2))
            .unwrap(),
    );
    // A unitary factor mixing must leave the physical second moments unchanged.
    let s = C64::from_re(1.0 / 2.0_f64.sqrt());
    let mixed_a = [scale(base, s), scale(base, -s)];
    let mixed_b = [scale(base, s), scale(base, s)];
    let mixed = prepared
        .evaluate_coherence_modes(&[&mixed_a, &mixed_b], mode_budget(2, 2))
        .unwrap();
    close(mixed.power_w, independent.power_w, independent.power_w * 2e-14);
    close(
        mixed.power_first_moment_w[0],
        independent.power_first_moment_w[0],
        independent.power_w * 2e-14,
    );
    assert_physical(&mixed, medium());
}

#[test]
fn original_budgets_wall_domain_refusal_and_cancellation_do_not_publish_partial_results() {
    let (positions, sources) = fixture();
    assert!(
        PreparedCompactRadiatedMomentum::new(medium(), 100.0, &positions, budget(2)).is_err()
    );
    let base =
        PreparedCompactRadiatedPower::new(medium(), 100.0, &positions, budget(3)).unwrap();
    let mut preparation_polls = 0;
    let prepared = base
        .clone()
        .with_momentum_cancellable(|| {
            preparation_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut polls = 0;
        assert!(matches!(
            base.clone().with_momentum_cancellable(|| {
                polls += 1;
                polls == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let mut application_polls = 0;
    let expected = prepared
        .evaluate_cancellable(&sources, || {
            application_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, application_polls / 2, application_polls] {
        let mut polls = 0;
        assert!(matches!(
            prepared.evaluate_cancellable(&sources, || {
                polls += 1;
                polls == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
        assert_eq!(prepared.evaluate(&sources).unwrap(), expected);
    }
    assert!(prepared.evaluate(&sources[..2]).is_err());
    let mut nonfinite = sources;
    nonfinite[1].force_on_fluid_n[0].im = f64::NAN;
    assert!(matches!(
        prepared.evaluate(&nonfinite),
        Err(AeroacError::NonFinite { .. })
    ));
    assert!(
        prepared
            .evaluate_coherence_modes(&[&sources], mode_budget(3, 0))
            .is_err()
    );
    assert!(
        prepared
            .evaluate_coherence_modes(&[&sources[..2]], mode_budget(3, 1))
            .is_err()
    );
    assert!(matches!(
        prepared.evaluate_coherence_modes_cancellable(
            &[&sources],
            mode_budget(3, 1),
            || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        prepared.evaluate_independent_sources_cancellable(&sources, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let silence = prepared
        .evaluate_coherence_modes(&[], mode_budget(3, 0))
        .unwrap();
    assert_eq!(silence.power_w, 0.0);
    assert_eq!(silence.wave_momentum_flux_n, [0.0; 3]);
    let empty =
        PreparedCompactRadiatedMomentum::new(medium(), 100.0, &[], budget(0)).unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap(), silence);
    let plane = IdealAcousticPlane::new(
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
        IdealPlaneBoundary::Rigid,
    )
    .unwrap();
    let observer =
        PreparedCompactSourceObserver::new(medium(), 100.0, [0.0, 0.0, 2.0], &positions, 3)
            .unwrap()
            .with_plane(plane, 6)
            .unwrap();
    let wall_power =
        PreparedCompactRadiatedPower::from_half_space(&observer, budget(6)).unwrap();
    assert!(matches!(
        wall_power.with_momentum(),
        Err(AeroacError::InvalidParameter { .. })
    ));
}
