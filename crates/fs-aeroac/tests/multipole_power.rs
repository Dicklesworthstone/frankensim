//! Physical power of volume stresses and their coherent compact-source mixtures.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::power::{CompactMultipolePhasor, PreparedMultipoleRadiatedPower};
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::{
    CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary, PreparedCompactRadiatedPower,
    PreparedCompactSourceObserver, RadiatedPowerBudget, RadiatedPowerModeBudget,
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

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-11 * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn fixture() -> ([[f64; 3]; 3], [CompactMultipolePhasor; 3]) {
    let positions = [
        [-0.25, 0.125, 0.0],
        [0.375, -0.25, 0.125],
        [0.0, 0.25, -0.375],
    ];
    let sources = core::array::from_fn(|i| {
        let t = (i + 1) as f64;
        CompactMultipolePhasor {
            compact: CompactSourcePhasor {
                volume_flow_m3_s: C64::new((t - 1.5) * 3e-4, (2.5 - t) * 2e-4),
                force_on_fluid_n: [
                    C64::new(0.3 * t, 0.1),
                    C64::new(-0.2, 0.4 * t),
                    C64::new(0.1 / t, -0.3),
                ],
            },
            quadrupole: QuadrupolePhasor {
                stress_moment_n_m: [
                    C64::new(0.11 * t, 0.05),
                    C64::new(-0.07, 0.13 * t),
                    C64::new(0.04 / t, -0.12),
                    C64::new(0.15 * t, -0.06),
                    C64::new(-0.04, -0.09 * t),
                    C64::new(0.07 * t, 0.03 / t),
                ],
            },
        }
    });
    (positions, sources)
}

fn tensor(source: QuadrupolePhasor) -> [[C64; 3]; 3] {
    let t = source.stress_moment_n_m;
    [[t[0], t[3], t[4]], [t[3], t[1], t[5]], [t[4], t[5], t[2]]]
}

fn scale(source: CompactMultipolePhasor, factor: C64) -> CompactMultipolePhasor {
    CompactMultipolePhasor {
        compact: CompactSourcePhasor {
            volume_flow_m3_s: source.compact.volume_flow_m3_s * factor,
            force_on_fluid_n: source.compact.force_on_fluid_n.map(|value| value * factor),
        },
        quadrupole: QuadrupolePhasor {
            stress_moment_n_m: source
                .quadrupole
                .stress_moment_n_m
                .map(|value| value * factor),
        },
    }
}

fn legendre(order: usize, x: f64) -> (f64, f64) {
    let (mut previous, mut value) = (0.0, 1.0);
    for degree in 1..=order {
        let next =
            ((2 * degree - 1) as f64 * x * value - (degree - 1) as f64 * previous) / degree as f64;
        previous = value;
        value = next;
    }
    (value, previous)
}

// Independent Gauss-Legendre/periodic angular quadrature. No production
// spherical moments or tensor-contraction coefficients enter this oracle.
fn sphere_integral(order: usize, mut field: impl FnMut([f64; 3]) -> f64) -> f64 {
    let azimuths = 2 * order;
    let mut integral = 0.0;
    for index in 0..order {
        let mut mu = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..20 {
            let (p, previous) = legendre(order, mu);
            let derivative = order as f64 * (mu * p - previous) / (mu * mu - 1.0);
            let step = p / derivative;
            mu -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        let (p, previous) = legendre(order, mu);
        let derivative = order as f64 * (mu * p - previous) / (mu * mu - 1.0);
        let weight = 2.0 / ((1.0 - mu * mu) * derivative * derivative);
        let radial = (1.0 - mu * mu).sqrt();
        for azimuth in 0..azimuths {
            let angle = TAU * azimuth as f64 / azimuths as f64;
            integral += weight * TAU / azimuths as f64
                * field([radial * angle.cos(), radial * angle.sin(), mu]);
        }
    }
    integral
}

fn far_power(
    acoustic: FluidAcousticState,
    frequency: f64,
    positions: &[[f64; 3]],
    sources: &[CompactMultipolePhasor],
) -> f64 {
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let impedance = acoustic.characteristic_impedance_pa_s_m();
    let integral = sphere_integral(48, |n| {
        let mut field = C64::ZERO;
        for (position, source) in positions.iter().zip(sources) {
            let phase = -k * n.iter().zip(position).map(|(a, b)| a * b).sum::<f64>();
            let direction = C64::new(phase.cos(), phase.sin());
            let mut amplitude = source.compact.volume_flow_m3_s.scale(impedance);
            let stress = tensor(source.quadrupole);
            for i in 0..3 {
                amplitude = amplitude + source.compact.force_on_fluid_n[i].scale(n[i]);
                for j in 0..3 {
                    amplitude = amplitude + C64::new(0.0, -k) * stress[i][j].scale(n[i] * n[j]);
                }
            }
            field = field + amplitude * direction;
        }
        field.norm_sq()
    });
    k * k * integral / (32.0 * PI * PI * impedance)
}

#[test]
fn isolated_full_stress_tensor_has_the_physical_si_power_and_frequency_scaling() {
    let (_, sources) = fixture();
    let acoustic = medium();
    let frequency = 197.0;
    let source = sources[0].quadrupole;
    let stress = tensor(source);
    let trace = stress[0][0] + stress[1][1] + stress[2][2];
    let frobenius = stress
        .iter()
        .flatten()
        .map(|value| value.norm_sq())
        .sum::<f64>();
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let expected = k.powi(4) * (trace.norm_sq() + 2.0 * frobenius)
        / (120.0 * PI * acoustic.characteristic_impedance_pa_s_m());
    let mut powers = Vec::new();
    for factor in [1.0, 2.0] {
        let prepared = PreparedMultipoleRadiatedPower::new(
            acoustic,
            frequency * factor,
            &[[0.25, -0.5, 0.125]],
            budget(1),
        )
        .unwrap();
        let actual = prepared.evaluate(&[source.into()]).unwrap();
        close(actual.power_w, expected * factor.powi(4));
        powers.push(actual.power_w);
        assert!(actual.scope.contains("stress"));
        assert_eq!(prepared.source_count(), 1);
    }
    close(powers[1], 16.0 * powers[0]);
}

#[test]
fn coherent_separated_mixtures_match_independent_angular_integration() {
    let (positions, sources) = fixture();
    for frequency in [7.0, 145.0, 700.0] {
        let prepared =
            PreparedMultipoleRadiatedPower::new(medium(), frequency, &positions, budget(3))
                .unwrap();
        close(
            prepared.evaluate(&sources).unwrap().power_w,
            far_power(medium(), frequency, &positions, &sources),
        );
        let stresses = sources.map(|s| CompactMultipolePhasor::from(s.quadrupole));
        close(
            prepared.evaluate(&stresses).unwrap().power_w,
            far_power(medium(), frequency, &positions, &stresses),
        );
        // Removing phases/cross terms would give the same independent result
        // for this rotation, but the coherent mixed power must change.
        let mut rotated = sources;
        rotated[1] = scale(rotated[1], C64::new(0.0, 1.0));
        let changed = prepared.evaluate(&rotated).unwrap().power_w;
        close(
            changed,
            far_power(medium(), frequency, &positions, &rotated),
        );
        assert!((changed - prepared.evaluate(&sources).unwrap().power_w).abs() > 1e-5 * changed);
    }
}

#[test]
fn stress_trace_interferes_with_monopoles_in_the_correct_time_phase() {
    let acoustic = FluidAcousticState::try_new(1.0, 1.0).unwrap();
    let prepared =
        PreparedMultipoleRadiatedPower::new(acoustic, 1.0 / TAU, &[[0.0; 3]], budget(1)).unwrap();
    let stress = QuadrupolePhasor {
        stress_moment_n_m: [
            C64::ONE,
            C64::ONE,
            C64::ONE,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ],
    };
    let mut source = CompactMultipolePhasor {
        compact: CompactSourcePhasor::monopole(C64::new(0.0, 1.0)),
        quadrupole: stress,
    };
    let stress_only = prepared.evaluate(&[stress.into()]).unwrap().power_w;
    close(stress_only, 1.0 / (8.0 * PI));
    assert!(prepared.evaluate(&[source]).unwrap().power_w < 1e-60);
    source.compact.volume_flow_m3_s.im = -1.0;
    close(
        prepared.evaluate(&[source]).unwrap().power_w,
        4.0 * stress_only,
    );

    let amplitude = 1e12;
    source.quadrupole.stress_moment_n_m[..3].fill(C64::from_re(amplitude));
    source.compact.volume_flow_m3_s.im = amplitude + 0.001;
    let residual = source.compact.volume_flow_m3_s.im - amplitude;
    close(
        prepared.evaluate(&[source]).unwrap().power_w,
        residual * residual / (8.0 * PI),
    );
}

#[test]
fn a_cancelling_subwavelength_tensor_pair_retains_its_weak_radiation() {
    let acoustic = FluidAcousticState::try_new(1.0, 1.0).unwrap();
    let mut source = CompactMultipolePhasor::default();
    source.quadrupole.stress_moment_n_m[2] = C64::ONE;
    for separation in [1.0 / 1_073_741_824.0, 1e-20, 1e-100] {
        let prepared = PreparedMultipoleRadiatedPower::new(
            acoustic,
            1.0 / TAU,
            &[[0.0; 3], [0.0, 0.0, separation]],
            budget(2),
        )
        .unwrap();
        let actual = prepared
            .evaluate(&[source, scale(source, C64::from_re(-1.0))])
            .unwrap()
            .power_w;
        // 2 <n_z^4 (1-cos(x n_z))> = x^2/7 + O(x^4).
        let expected = separation * separation / (56.0 * PI);
        assert!(actual > 0.0);
        assert!(
            (actual / expected - 1.0).abs() < 2e-13,
            "separation={separation:e}, actual={actual:e}, expected={expected:e}"
        );
    }
}

#[test]
fn the_extension_preserves_the_original_low_order_results() {
    let (positions, sources) = fixture();
    let compact_sources = sources.map(|source| source.compact);
    let multipoles = compact_sources.map(CompactMultipolePhasor::from);
    let original =
        PreparedCompactRadiatedPower::new(medium(), 383.0, &positions, budget(3)).unwrap();
    let extended = original.clone().with_quadrupoles().unwrap();
    assert_eq!(
        original.evaluate(&compact_sources).unwrap().power_w,
        extended.evaluate(&multipoles).unwrap().power_w
    );
    assert_eq!(
        original
            .evaluate_independent_sources(&compact_sources)
            .unwrap()
            .power_w,
        extended
            .evaluate_independent_sources(&multipoles)
            .unwrap()
            .power_w
    );
    assert_eq!(
        extended.compact_power().frequency_hz(),
        original.frequency_hz()
    );
}

#[test]
fn independent_sources_match_explicit_factors_without_losing_internal_cross_terms() {
    let (positions, sources) = fixture();
    let prepared =
        PreparedMultipoleRadiatedPower::new(medium(), 157.0, &positions, budget(3)).unwrap();
    let expected = sources
        .iter()
        .map(|source| far_power(medium(), 157.0, &[[0.0; 3]], &[*source]))
        .sum::<f64>();
    close(
        prepared
            .evaluate_independent_sources(&sources)
            .unwrap()
            .power_w,
        expected,
    );
    let mut factors = [[CompactMultipolePhasor::default(); 3]; 3];
    for i in 0..3 {
        factors[i][i] = sources[i];
    }
    let modes: Vec<&[CompactMultipolePhasor]> = factors.iter().map(|row| row.as_slice()).collect();
    close(
        prepared
            .evaluate_coherence_modes(
                &modes,
                RadiatedPowerModeBudget {
                    max_modes: 3,
                    max_source_pair_mode_products: 27,
                },
            )
            .unwrap()
            .power_w,
        expected,
    );

    let colocated =
        PreparedMultipoleRadiatedPower::new(medium(), 157.0, &[[0.0; 3]; 2], budget(2)).unwrap();
    let opposing = [sources[0], scale(sources[0], C64::from_re(-1.0))];
    assert_eq!(colocated.evaluate(&opposing).unwrap().power_w, 0.0);
    assert!(
        colocated
            .evaluate_independent_sources(&opposing)
            .unwrap()
            .power_w
            > 0.0
    );
}

#[test]
fn physical_near_field_flux_integrates_to_the_same_radius_independent_power() {
    let (_, sources) = fixture();
    let source = sources[1];
    let frequency = 83.0;
    let prepared =
        PreparedMultipoleRadiatedPower::new(medium(), frequency, &[[0.0; 3]], budget(1)).unwrap();
    let expected = prepared.evaluate(&[source]).unwrap().power_w;
    for radius in [0.15, 0.7, 5.0] {
        let flux = sphere_integral(12, |n| {
            let receiver = n.map(|value| radius * value);
            let compact =
                PreparedCompactSourceObserver::new(medium(), frequency, receiver, &[[0.0; 3]], 1)
                    .unwrap()
                    .evaluate(&[source.compact])
                    .unwrap();
            let stress =
                PreparedQuadrupoleObserver::new(medium(), frequency, receiver, &[[0.0; 3]], 1)
                    .unwrap()
                    .evaluate(&[source.quadrupole])
                    .unwrap();
            let pressure = compact.pressure_pa + stress.pressure_pa;
            let velocity = (0..3).fold(C64::ZERO, |sum, axis| {
                sum + (compact.particle_velocity_m_s[axis] + stress.particle_velocity_m_s[axis])
                    .scale(n[axis])
            });
            0.5 * (pressure * velocity.conj()).re * radius * radius
        });
        close(flux, expected);
    }
}

fn reflected(
    source: CompactMultipolePhasor,
    normal: [f64; 3],
    sign: f64,
) -> CompactMultipolePhasor {
    let rotation: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| (if i == j { 1.0 } else { 0.0 }) - 2.0 * normal[i] * normal[j])
    });
    let stress = tensor(source.quadrupole);
    let reflected_tensor: [[C64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .flat_map(|a| (0..3).map(move |b| (a, b)))
                .fold(C64::ZERO, |sum, (a, b)| {
                    sum + stress[a][b].scale(sign * rotation[i][a] * rotation[j][b])
                })
        })
    });
    CompactMultipolePhasor {
        compact: CompactSourcePhasor {
            volume_flow_m3_s: source.compact.volume_flow_m3_s.scale(sign),
            force_on_fluid_n: core::array::from_fn(|i| {
                (0..3).fold(C64::ZERO, |sum, a| {
                    sum + source.compact.force_on_fluid_n[a].scale(sign * rotation[i][a])
                })
            }),
        },
        quadrupole: QuadrupolePhasor {
            stress_moment_n_m: [
                reflected_tensor[0][0],
                reflected_tensor[1][1],
                reflected_tensor[2][2],
                reflected_tensor[0][1],
                reflected_tensor[0][2],
                reflected_tensor[1][2],
            ],
        },
    }
}

#[test]
fn oblique_ideal_walls_transform_the_tensor_and_retain_image_coherence() {
    let (mut positions, sources) = fixture();
    let normal = [0.6, 0.8, 0.0];
    let origin = [0.2, -0.1, 0.3];
    for position in &mut positions {
        for i in 0..3 {
            position[i] += origin[i] + normal[i];
        }
    }
    let frequency = 105.0;
    for (boundary, sign) in [
        (IdealPlaneBoundary::Rigid, 1.0),
        (IdealPlaneBoundary::PressureRelease, -1.0),
    ] {
        let plane = IdealAcousticPlane::new(origin, normal, boundary).unwrap();
        let field =
            PreparedCompactSourceObserver::new(medium(), frequency, [2.0, 3.0, 4.0], &positions, 3)
                .unwrap()
                .with_plane(plane, 6)
                .unwrap();
        let prepared = PreparedCompactRadiatedPower::from_half_space(&field, budget(6))
            .unwrap()
            .with_quadrupoles()
            .unwrap();
        assert_eq!(prepared.source_count(), 3);
        assert_eq!(prepared.compact_power().source_path_count(), 6);
        let mut image_positions = Vec::new();
        let mut image_sources = Vec::new();
        let mut independent = 0.0;
        for (position, source) in positions.iter().zip(sources) {
            let distance = (0..3)
                .map(|i| normal[i] * (position[i] - origin[i]))
                .sum::<f64>();
            let image = core::array::from_fn(|i| position[i] - 2.0 * normal[i] * distance);
            let image_source = reflected(source, normal, sign);
            image_positions.extend([*position, image]);
            image_sources.extend([source, image_source]);
            independent += far_power(
                medium(),
                frequency,
                &[*position, image],
                &[source, image_source],
            ) / 2.0;
        }
        close(
            prepared.evaluate(&sources).unwrap().power_w,
            far_power(medium(), frequency, &image_positions, &image_sources) / 2.0,
        );
        close(
            prepared
                .evaluate_independent_sources(&sources)
                .unwrap()
                .power_w,
            independent,
        );
        assert!(prepared.scope().contains("physical") || prepared.scope().contains("fluid domain"));
    }
}

#[test]
fn budgets_nonfinite_sources_and_cancellation_refuse_without_mutating_the_owner() {
    let (positions, mut sources) = fixture();
    assert!(PreparedMultipoleRadiatedPower::new(medium(), 200.0, &positions, budget(2)).is_err());
    assert!(matches!(
        PreparedMultipoleRadiatedPower::new_cancellable(
            medium(),
            200.0,
            &positions,
            budget(3),
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let original =
        PreparedCompactRadiatedPower::new(medium(), 200.0, &positions, budget(3)).unwrap();
    let mut polls = 0;
    assert!(matches!(
        original.clone().with_quadrupoles_cancellable(|| {
            polls += 1;
            polls == 3
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let prepared = original.with_quadrupoles().unwrap();
    let expected = prepared.evaluate(&sources).unwrap().power_w;
    assert!(prepared.evaluate(&sources[..2]).is_err());
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&sources],
                RadiatedPowerModeBudget {
                    max_modes: 1,
                    max_source_pair_mode_products: 8,
                }
            )
            .is_err()
    );
    for stop in [1, 4, 8] {
        let mut polls = 0;
        assert!(matches!(
            prepared.evaluate_cancellable(&sources, || {
                polls += 1;
                polls == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    assert!(matches!(
        prepared.evaluate_independent_sources_cancellable(&sources, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        prepared.evaluate_coherence_modes_cancellable(
            &[&sources],
            RadiatedPowerModeBudget {
                max_modes: 1,
                max_source_pair_mode_products: 9
            },
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(prepared.evaluate(&sources).unwrap().power_w, expected);
    sources[1].quadrupole.stress_moment_n_m[4].re = f64::NAN;
    assert!(matches!(
        prepared.evaluate(&sources),
        Err(AeroacError::NonFinite { .. })
    ));
}
