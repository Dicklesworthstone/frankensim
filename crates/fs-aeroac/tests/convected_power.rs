//! Convected analytic power against enclosing physical-field surface integrals.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{
    CompactSourcePhasor, PreparedCompactRadiatedPower, PreparedConvectedMonopolePower,
    RadiatedPowerBudget, RadiatedPowerModeBudget,
    array::convected::{ConvectedMonopoleLimits, PreparedConvectedMonopoleObserver},
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
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn unit_power(f: f64, mach: f64) -> f64 {
    medium().density_kg_m3() * (TAU * f).powi(2)
        / (8.0 * PI * medium().sound_speed_m_s() * (1.0 - mach * mach).powi(2))
}
fn fixture() -> ([[f64; 3]; 3], [C64; 3]) {
    (
        [
            [-0.25, 0.125, 0.0],
            [0.375, -0.25, 0.125],
            [0.0, 0.25, -0.375],
        ],
        [
            C64::new(4e-4, 2e-4),
            C64::new(-1e-4, 3e-4),
            C64::new(2e-4, -1e-4),
        ],
    )
}

// Independent Gauss-Legendre x periodic trapezoid integration of actual
// near-field pressure, velocity and convective energy flux on a sphere.
// This does not call the production pair moments or transformed far field.
fn sphere_integral(mut value: impl FnMut([f64; 3]) -> f64) -> f64 {
    let order = 48;
    let azimuths = 96;
    let legendre = |x: f64| {
        let (mut previous, mut p) = (0.0, 1.0);
        for degree in 1..=order {
            let next =
                ((2 * degree - 1) as f64 * x * p - (degree - 1) as f64 * previous) / degree as f64;
            previous = p;
            p = next;
        }
        (p, order as f64 * (x * p - previous) / (x * x - 1.0))
    };
    let mut result = 0.0;
    for index in 0..order {
        let mut mu = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..24 {
            let (p, derivative) = legendre(mu);
            let change = p / derivative;
            mu -= change;
            if change.abs() < 1e-15 {
                break;
            }
        }
        let (_, derivative) = legendre(mu);
        let weight = 2.0 / ((1.0 - mu * mu) * derivative * derivative);
        let radial = (1.0 - mu * mu).sqrt();
        for azimuth in 0..azimuths {
            let angle = TAU * azimuth as f64 / azimuths as f64;
            result += weight * TAU / azimuths as f64
                * value([radial * angle.cos(), radial * angle.sin(), mu]);
        }
    }
    result
}

#[test]
fn isolated_si_power_and_mach_amplification() {
    let c = medium().sound_speed_m_s();
    let q = C64::new(2e-4, -3e-4);
    for mach in [0.0, 0.2, 0.6, 0.9] {
        let prepared = PreparedConvectedMonopolePower::new(
            medium(),
            [mach * c, 0.0, 0.0],
            413.0,
            &[[0.125, -0.25, 0.375]],
            0.95,
            budget(1),
        )
        .unwrap();
        close(
            prepared.evaluate(&[q]).unwrap().power_w,
            unit_power(413.0, mach) * q.norm_sq(),
            3e-15,
        );
    }
}

#[test]
fn quiescent_limit_preserves_full_array_interference() {
    let (positions, q) = fixture();
    let convected =
        PreparedConvectedMonopolePower::new(medium(), [0.0; 3], 230.0, &positions, 0.0, budget(3))
            .unwrap();
    let still = PreparedCompactRadiatedPower::new(medium(), 230.0, &positions, budget(3)).unwrap();
    close(
        convected.evaluate(&q).unwrap().power_w,
        still
            .evaluate(&q.map(CompactSourcePhasor::monopole))
            .unwrap()
            .power_w,
        2e-14,
    );
}

#[test]
fn total_transport_matches_near_field_on_different_enclosing_surfaces() {
    let (positions, q) = fixture();
    let flow = [0.36, -0.48, 0.0].map(|m| m * medium().sound_speed_m_s());
    let prepared =
        PreparedConvectedMonopolePower::new(medium(), flow, 180.0, &positions, 0.7, budget(3))
            .unwrap();
    let expected = prepared.evaluate(&q).unwrap().power_w;
    for radius in [1.4, 3.7] {
        let power = sphere_integral(|normal| {
            let observer = PreparedConvectedMonopoleObserver::new(
                medium(),
                flow,
                180.0,
                normal.map(|n| radius * n),
                &positions,
                3,
                ConvectedMonopoleLimits {
                    maximum_mach: 0.7,
                    minimum_distance_m: 0.01,
                },
            )
            .unwrap();
            let field = observer.evaluate(&q).unwrap();
            radius
                * radius
                * (0..3)
                    .map(|axis| normal[axis] * field.mean_perturbation_energy_flux_w_m2[axis])
                    .sum::<f64>()
        });
        close(power, expected, 8e-11);
    }
}

#[test]
fn subwavelength_opposed_sources_preserve_directional_radiation() {
    let mach = 0.6;
    let beta_squared = 1.0 - mach * mach;
    let c = medium().sound_speed_m_s();
    let f = 180.0;
    let k = TAU * f / c;
    for x in [1e-6, 1e-10, 1e-14] {
        for along_flow in [false, true] {
            let displacement = if along_flow {
                [x * beta_squared / k, 0.0, 0.0]
            } else {
                [0.0, x * beta_squared.sqrt() / k, 0.0]
            };
            let prepared = PreparedConvectedMonopolePower::new(
                medium(),
                [mach * c, 0.0, 0.0],
                f,
                &[[0.0; 3], displacement],
                0.7,
                budget(2),
            )
            .unwrap();
            let power = prepared.evaluate(&[C64::ONE, -C64::ONE]).unwrap().power_w;
            // Independent small-separation angular moment:
            // <(1-M.n)[(n-M).d_hat]^2> = (1+5(M.d_hat)^2)/3.
            let directional = if along_flow {
                1.0 + 5.0 * mach * mach
            } else {
                1.0
            };
            close(
                power,
                unit_power(f, mach) * x * x * directional / 3.0,
                2e-12,
            );
        }
    }
}

#[test]
fn coincident_cancellation_and_explicit_independence() {
    let flow = [120.0, -60.0, 30.0];
    let position = [[0.25, -0.125, 0.375]; 3];
    let prepared =
        PreparedConvectedMonopolePower::new(medium(), flow, 180.0, &position, 0.6, budget(3))
            .unwrap();
    let small = C64::new(2e-4, -3e-4);
    let large = small.scale(1e16);
    let single =
        PreparedConvectedMonopolePower::new(medium(), flow, 180.0, &position[..1], 0.6, budget(1))
            .unwrap();
    assert_eq!(prepared.position_count(), 1);
    close(
        prepared.evaluate(&[large, small, -large]).unwrap().power_w,
        single.evaluate(&[small]).unwrap().power_w,
        3e-15,
    );
    let pair =
        PreparedConvectedMonopolePower::new(medium(), flow, 180.0, &position[..2], 0.6, budget(2))
            .unwrap();
    assert_eq!(pair.evaluate(&[small, -small]).unwrap().power_w, 0.0);
    let a = [small, C64::ZERO];
    let b = [C64::ZERO, -small];
    let independent = pair
        .evaluate_coherence_modes(
            &[&a, &b],
            RadiatedPowerModeBudget {
                max_modes: 2,
                max_source_pair_mode_products: 8,
            },
        )
        .unwrap()
        .power_w;
    close(
        independent,
        2.0 * single.evaluate(&[small]).unwrap().power_w,
        3e-15,
    );
}

#[test]
fn source_owner_translation_rotation_and_input_order_invariance() {
    let (positions, q) = fixture();
    let flow = [120.0, -60.0, 30.0];
    let observer = PreparedConvectedMonopoleObserver::new(
        medium(),
        flow,
        180.0,
        [0.0, 0.0, 2.0],
        &positions,
        3,
        ConvectedMonopoleLimits {
            maximum_mach: 0.6,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let prepared =
        PreparedConvectedMonopolePower::from_observer(&observer, 0.6, budget(3)).unwrap();
    assert_eq!(prepared.source_count(), 3);
    assert_eq!(prepared.flow_velocity_m_s(), flow);
    assert_eq!(prepared.frequency_hz(), 180.0);
    let expected = prepared.evaluate(&q).unwrap().power_w;
    let rotate = |v: [f64; 3]| [v[2], -v[0], -v[1]];
    let rotated = positions.map(|p| rotate(p).map(|v| v + 128.0));
    let other = PreparedConvectedMonopolePower::new(
        medium(),
        rotate(flow),
        180.0,
        &[rotated[2], rotated[0], rotated[1]],
        0.6,
        budget(3),
    )
    .unwrap();
    close(
        other.evaluate(&[q[2], q[0], q[1]]).unwrap().power_w,
        expected,
        2e-14,
    );
}

#[test]
fn admission_cancellation_and_nonfinite_refusals() {
    let position = [[0.0; 3], [0.1, 0.2, 0.3]];
    for mach in [1.0, f64::NAN, -0.1] {
        assert!(
            PreparedConvectedMonopolePower::new(
                medium(),
                [0.0; 3],
                180.0,
                &position,
                mach,
                budget(2)
            )
            .is_err()
        );
    }
    assert!(
        PreparedConvectedMonopolePower::new(
            medium(),
            [200.0, 0.0, 0.0],
            180.0,
            &position,
            0.5,
            budget(2)
        )
        .is_err()
    );
    assert!(
        PreparedConvectedMonopolePower::new(
            medium(),
            [0.0; 3],
            180.0,
            &position,
            0.5,
            RadiatedPowerBudget {
                max_sources: 2,
                max_source_pairs: 3
            }
        )
        .is_err()
    );
    assert!(matches!(
        PreparedConvectedMonopolePower::new_cancellable(
            medium(),
            [0.0; 3],
            180.0,
            &position,
            0.5,
            budget(2),
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let prepared = PreparedConvectedMonopolePower::new(
        medium(),
        [100.0, 0.0, 0.0],
        180.0,
        &position,
        0.5,
        budget(2),
    )
    .unwrap();
    assert!(prepared.evaluate(&[C64::ONE]).is_err());
    assert!(
        prepared
            .evaluate(&[C64::ONE, C64::new(f64::NAN, 0.0)])
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        prepared.evaluate_cancellable(&[C64::ONE; 2], || {
            polls += 1;
            polls >= 5
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&[C64::ONE; 2]],
                RadiatedPowerModeBudget {
                    max_modes: 1,
                    max_source_pair_mode_products: 3
                }
            )
            .is_err()
    );
}
