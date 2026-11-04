//! G1/G3: physical-domain image interference, including independent sources.

use super::*;
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};

fn wall(boundary: IdealPlaneBoundary) -> IdealAcousticPlane {
    IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundary).unwrap()
}

fn observer(frequency: f64, positions: &[[f64; 3]]) -> PreparedCompactSourceObserver {
    PreparedCompactSourceObserver::new(
        medium(),
        frequency,
        [1.0, 1.0, 2.0],
        positions,
        positions.len(),
    )
    .unwrap()
}

#[test]
fn plane_power_retains_monopole_and_both_dipole_image_parities() {
    let acoustic = medium();
    let frequency = 180.0;
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let x = 1.3_f64;
    let positions = [[0.0, 0.0, x / (2.0 * k)]];
    let j0 = x.sin() / x;
    let j1_over_x = (x.sin() / x - x.cos()) / (x * x);
    let sources = [
        CompactSourcePhasor::monopole(C64::new(2e-4, 1e-4)),
        CompactSourcePhasor::force_dipole([C64::ONE, C64::ZERO, C64::ZERO]),
        CompactSourcePhasor::force_dipole([C64::ZERO, C64::ZERO, C64::ONE]),
    ];
    for (boundary, sign) in [
        (IdealPlaneBoundary::Rigid, 1.0),
        (IdealPlaneBoundary::PressureRelease, -1.0),
    ] {
        let plane = wall(boundary);
        let direct = observer(frequency, &positions);
        let separate = direct.clone().with_plane(plane, 2).unwrap();
        let folded = direct.with_orthogonal_planes(&[plane], 2).unwrap();
        let power = PreparedCompactRadiatedPower::from_observer(&folded, budget(2)).unwrap();
        let half_space =
            PreparedCompactRadiatedPower::from_half_space(&separate, budget(2)).unwrap();
        assert_eq!(power.source_count(), 1);
        assert_eq!(power.source_path_count(), 2);
        assert_eq!(power.position_count(), 2);
        let factors = [
            1.0 + sign * j0,
            1.0 + sign * 3.0 * j1_over_x,
            1.0 - sign * 3.0 * (j0 - 2.0 * j1_over_x),
        ];
        for (source, factor) in sources.into_iter().zip(factors) {
            let actual = power.evaluate(&[source]).unwrap();
            close(
                actual.power_w,
                isolated(acoustic, frequency, source) * factor,
            );
            assert_eq!(actual, half_space.evaluate(&[source]).unwrap());
            close(
                power
                    .evaluate_independent_sources(&[source])
                    .unwrap()
                    .power_w,
                actual.power_w,
            );
            assert_eq!(actual.scope, power.scope());
        }
    }
}

#[test]
fn tiny_pressure_release_power_survives_direct_image_cancellation() {
    let frequency = 180.0;
    let k = TAU * frequency / medium().sound_speed_m_s();
    let source = CompactSourcePhasor::monopole(C64::new(3e-4, -2e-4));
    let p0 = isolated(medium(), frequency, source);
    for x in [2.0_f64.powi(-40), 2.0_f64.powi(-80)] {
        let field = observer(frequency, &[[0.0, 0.0, x / (2.0 * k)]])
            .with_orthogonal_planes(&[wall(IdealPlaneBoundary::PressureRelease)], 2)
            .unwrap();
        let power = PreparedCompactRadiatedPower::from_observer(&field, budget(2)).unwrap();
        let expected = p0 * x * x / 6.0;
        for actual in [
            power.evaluate(&[source]).unwrap().power_w,
            power
                .evaluate_independent_sources(&[source])
                .unwrap()
                .power_w,
        ] {
            assert!(
                (actual / expected - 1.0).abs() < 2e-11,
                "x={x:e}, actual={actual:e}, expected={expected:e}"
            );
        }
    }
}

// Map the original independent Gauss/azimuth rule to the positive hemisphere.
fn hemisphere_integral(mut value: impl FnMut([f64; 3]) -> f64) -> f64 {
    sphere_integral(|n| {
        let mu = 0.5 * (n[2] + 1.0);
        let scale = ((1.0 - mu * mu) / (1.0 - n[2] * n[2])).sqrt();
        0.5 * value([n[0] * scale, n[1] * scale, mu])
    })
}

fn wedge_integral(angle: f64, mut value: impl FnMut([f64; 3]) -> f64) -> f64 {
    let order = 32;
    let mut nodes = Vec::new();
    for index in 0..order {
        let mut x = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..20 {
            let (p, previous) = legendre(order, x);
            let derivative = order as f64 * (x * p - previous) / (x * x - 1.0);
            let change = p / derivative;
            x -= change;
            if change.abs() < 1e-15 {
                break;
            }
        }
        let (p, previous) = legendre(order, x);
        let derivative = order as f64 * (x * p - previous) / (x * x - 1.0);
        nodes.push((x, 2.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    let mut integral = 0.0;
    for &(mu, weight) in &nodes {
        let radial = (1.0 - mu * mu).sqrt();
        for &(x, azimuth_weight) in &nodes {
            let azimuth = 0.5 * angle * (x + 1.0);
            integral += weight
                * azimuth_weight
                * 0.5
                * angle
                * value([radial * azimuth.cos(), radial * azimuth.sin(), mu]);
        }
    }
    integral
}

#[test]
fn physical_hemisphere_flux_and_rigid_frame_changes_match_total_power() {
    let frequency = 180.0;
    let (positions, sources) = fixture();
    let positions = positions.map(|p| [p[0], p[1], p[2].abs() + 0.3]);
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = wall(boundary);
        let field = observer(frequency, &positions)
            .with_orthogonal_planes(&[plane], 6)
            .unwrap();
        let power = PreparedCompactRadiatedPower::from_observer(&field, budget(6))
            .unwrap()
            .evaluate(&sources)
            .unwrap()
            .power_w;
        for radius in [1.8, 3.5] {
            let flux = hemisphere_integral(|normal| {
                let receiver = normal.map(|v| v * radius);
                let field = PreparedCompactSourceObserver::new(
                    medium(),
                    frequency,
                    receiver,
                    &positions,
                    3,
                )
                .unwrap()
                .with_plane(plane, 6)
                .unwrap()
                .evaluate(&sources)
                .unwrap();
                radius
                    * radius
                    * (0..3)
                        .map(|axis| normal[axis] * field.active_intensity_w_m2[axis])
                        .sum::<f64>()
            });
            close(power, flux);
        }
        let moved_positions = positions.map(|p| [p[2] + 2.0, p[0] - 0.5, p[1] + 1.0]);
        let moved_sources = sources.map(|s| CompactSourcePhasor {
            force_on_fluid_n: [
                s.force_on_fluid_n[2],
                s.force_on_fluid_n[0],
                s.force_on_fluid_n[1],
            ],
            ..s
        });
        let moved_plane =
            IdealAcousticPlane::new([2.0, -0.5, 1.0], [1.0, 0.0, 0.0], boundary).unwrap();
        let moved = PreparedCompactSourceObserver::new(
            medium(),
            frequency,
            [4.0, 0.5, 2.0],
            &moved_positions,
            3,
        )
        .unwrap()
        .with_orthogonal_planes(&[moved_plane], 6)
        .unwrap();
        close(
            PreparedCompactRadiatedPower::from_observer(&moved, budget(6))
                .unwrap()
                .evaluate(&moved_sources)
                .unwrap()
                .power_w,
            power,
        );
    }
}

#[test]
fn wedge_solid_angle_and_perpendicular_corner_share_the_original_image_group() {
    let frequency = 180.0;
    let (_, fixture_sources) = fixture();
    let sources = [fixture_sources[0], fixture_sources[1]];
    for (order, boundaries) in [
        (
            2,
            [
                IdealPlaneBoundary::Rigid,
                IdealPlaneBoundary::PressureRelease,
            ],
        ),
        (3, [IdealPlaneBoundary::Rigid; 2]),
        (
            4,
            [
                IdealPlaneBoundary::Rigid,
                IdealPlaneBoundary::PressureRelease,
            ],
        ),
    ] {
        let angle = PI / order as f64;
        let planes = [
            IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], boundaries[0]).unwrap(),
            IdealAcousticPlane::new([0.0; 3], [angle.sin(), -angle.cos(), 0.0], boundaries[1])
                .unwrap(),
        ];
        let positions = [
            [0.4 * (0.3 * angle).cos(), 0.4 * (0.3 * angle).sin(), 0.1],
            [0.3 * (0.7 * angle).cos(), 0.3 * (0.7 * angle).sin(), -0.2],
        ];
        let receiver = [angle.cos() + 2.0, angle.sin(), 1.0];
        let direct =
            PreparedCompactSourceObserver::new(medium(), frequency, receiver, &positions, 2)
                .unwrap();
        let wedge = direct.clone().with_wedge(planes, order, 4 * order).unwrap();
        let prepared =
            PreparedCompactRadiatedPower::from_observer(&wedge, budget(4 * order)).unwrap();
        let actual = prepared.evaluate(&sources).unwrap().power_w;
        let radius = 2.4;
        let flux = wedge_integral(angle, |n| {
            let field = PreparedCompactSourceObserver::new(
                medium(),
                frequency,
                n.map(|v| v * radius),
                &positions,
                2,
            )
            .unwrap()
            .with_wedge(planes, order, 4 * order)
            .unwrap()
            .evaluate(&sources)
            .unwrap();
            radius
                * radius
                * (0..3)
                    .map(|axis| n[axis] * field.active_intensity_w_m2[axis])
                    .sum::<f64>()
        });
        close(actual, flux);
        if order == 2 {
            let corner = direct.with_orthogonal_planes(&planes, 8).unwrap();
            close(
                PreparedCompactRadiatedPower::from_observer(&corner, budget(8))
                    .unwrap()
                    .evaluate(&sources)
                    .unwrap()
                    .power_w,
                actual,
            );
        }
    }
}

#[test]
fn independent_physical_sources_keep_their_images_coherent_and_budgets_count_paths() {
    let frequency = 180.0;
    let position = [0.2, -0.1, 0.3];
    let (_, sources) = fixture();
    let source = sources[0];
    let opposite = scale(source, -C64::ONE);
    let field = observer(frequency, &[position, position])
        .with_orthogonal_planes(&[wall(IdealPlaneBoundary::Rigid)], 4)
        .unwrap();
    for insufficient in [
        budget(2),
        RadiatedPowerBudget {
            max_source_pairs: 15,
            ..budget(4)
        },
    ] {
        assert!(PreparedCompactRadiatedPower::from_observer(&field, insufficient).is_err());
    }
    let mut polls = 0;
    let prepared =
        PreparedCompactRadiatedPower::from_observer_cancellable(&field, budget(4), || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(prepared.source_count(), 2);
    assert_eq!(prepared.source_path_count(), 4);
    assert_eq!(prepared.position_count(), 2);
    assert_eq!(prepared.evaluate(&[source, opposite]).unwrap().power_w, 0.0);
    let a = [source, CompactSourcePhasor::default()];
    let b = [CompactSourcePhasor::default(), opposite];
    let independent = prepared
        .evaluate_independent_sources(&[source, opposite])
        .unwrap()
        .power_w;
    close(independent, 2.0 * prepared.evaluate(&a).unwrap().power_w);
    let modes = RadiatedPowerModeBudget {
        max_modes: 2,
        max_source_pair_mode_products: 32,
    };
    close(
        prepared
            .evaluate_coherence_modes(&[&a, &b], modes)
            .unwrap()
            .power_w,
        independent,
    );
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&a, &b],
                RadiatedPowerModeBudget {
                    max_source_pair_mode_products: 31,
                    ..modes
                }
            )
            .is_err()
    );
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            PreparedCompactRadiatedPower::from_observer_cancellable(&field, budget(4), || {
                count += 1;
                count == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let mut queries = 0;
    let result = prepared
        .evaluate_independent_sources_cancellable(&[source, opposite], || {
            queries += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        prepared.evaluate_independent_sources_cancellable(&[source, opposite], || {
            count += 1;
            count == queries
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        prepared
            .evaluate_independent_sources(&[source, opposite])
            .unwrap(),
        result
    );
    assert!(prepared.evaluate(&[source; 4]).is_err());
}
