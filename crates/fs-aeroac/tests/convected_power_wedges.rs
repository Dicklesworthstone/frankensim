//! G0/G1/G3/G4: complete dihedral radiation and actual physical-sector flux.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{
    CompactSourcePhasor, PreparedCompactRadiatedPower, PreparedConvectedMonopolePower,
    RadiatedPowerBudget, RadiatedPowerModeBudget,
    array::convected::{ConvectedMonopoleLimits, PreparedConvectedMonopoleObserver},
    array::{IdealAcousticPlane, IdealPlaneBoundary},
    power::CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 331.0).unwrap()
}

fn budget(paths: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget {
        max_sources: paths,
        max_source_pairs: paths * paths,
    }
}

fn sign(boundary: IdealPlaneBoundary) -> f64 {
    match boundary {
        IdealPlaneBoundary::Rigid => 1.0,
        IdealPlaneBoundary::PressureRelease => -1.0,
    }
}

fn polar(radius: f64, angle: f64, z: f64) -> [f64; 3] {
    [radius * angle.cos(), radius * angle.sin(), z]
}

// The physical fluid sector is 0 < azimuth < pi/order, with an unchanged z edge.
fn planes(order: usize, walls: [IdealPlaneBoundary; 2]) -> [IdealAcousticPlane; 2] {
    let angle = PI / order as f64;
    [
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], walls[0]).unwrap(),
        IdealAcousticPlane::new([0.0; 3], [angle.sin(), -angle.cos(), 0.0], walls[1]).unwrap(),
    ]
}

fn observer(
    frequency: f64,
    flow: [f64; 3],
    point: [f64; 3],
    positions: &[[f64; 3]],
    walls: [IdealAcousticPlane; 2],
    order: usize,
) -> PreparedConvectedMonopoleObserver {
    PreparedConvectedMonopoleObserver::new(
        medium(),
        flow,
        frequency,
        point,
        positions,
        positions.len(),
        ConvectedMonopoleLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 1e-6,
        },
    )
    .unwrap()
    .with_wedge(walls, order, positions.len() * 2 * order)
    .unwrap()
}

fn fixture(order: usize) -> ([[f64; 3]; 3], [C64; 3]) {
    let angle = PI / order as f64;
    (
        [
            polar(0.22, 0.19 * angle, -0.12),
            polar(0.37, 0.63 * angle, 0.17),
            polar(0.16, 0.81 * angle, 0.29),
        ],
        [
            C64::new(4e-4, 2e-4),
            C64::new(-1e-4, 3e-4),
            C64::new(2e-4, -1e-4),
        ],
    )
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn sinc(value: f64) -> f64 {
    if value == 0.0 {
        1.0
    } else {
        value.sin() / value
    }
}

#[test]
fn single_source_matches_independent_dihedral_chord_sum() {
    let frequency = 413.0;
    let radius = 0.31;
    let q = C64::new(2e-4, -3e-4);
    for order in [2, 3, 4, 5] {
        let opening = PI / order as f64;
        let azimuth = 0.37 * opening;
        for walls in [
            [IdealPlaneBoundary::Rigid; 2],
            [IdealPlaneBoundary::PressureRelease; 2],
            [
                IdealPlaneBoundary::Rigid,
                IdealPlaneBoundary::PressureRelease,
            ],
            [
                IdealPlaneBoundary::PressureRelease,
                IdealPlaneBoundary::Rigid,
            ],
        ] {
            if walls[0] != walls[1] && order % 2 != 0 {
                continue;
            }
            for mach in [0.0_f64, 0.6] {
                let beta = (1.0 - mach * mach).sqrt();
                let k = TAU * frequency / medium().sound_speed_m_s();
                let field = observer(
                    frequency,
                    [0.0, 0.0, mach * medium().sound_speed_m_s()],
                    polar(0.9, opening / 2.0, 0.8),
                    &[polar(radius, azimuth, 0.17)],
                    planes(order, walls),
                    order,
                );
                let power =
                    PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(2 * order))
                        .unwrap();
                // Each group element contributes once after division by 2m.
                // For axial flow, every single-source image separation is
                // transverse, so its exact moment is sinc(k*chord/beta).
                // Direct trigonometric chords do not use the native images.
                let mut factor = Dd::ZERO;
                let mut rotation_sign = 1.0;
                for index in 0..order {
                    let angle = index as f64 * opening;
                    factor = factor
                        + Dd::from_f64(rotation_sign)
                            * (Dd::from_f64(sinc(2.0 * k * radius * angle.sin() / beta))
                                + Dd::from_f64(
                                    sign(walls[0])
                                        * sinc(2.0 * k * radius * (angle - azimuth).sin() / beta),
                                ));
                    rotation_sign *= sign(walls[0]) * sign(walls[1]);
                }
                let isolated = medium().density_kg_m3() * (TAU * frequency).powi(2) * q.norm_sq()
                    / (8.0 * PI * medium().sound_speed_m_s() * (1.0 - mach * mach).powi(2));
                let result = power.evaluate(&[q]).unwrap();
                close(result.power_w, isolated * factor.to_f64(), 3e-10);
                assert_eq!(result.scope, CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE);
                assert_eq!(power.source_count(), 1);
                assert_eq!(power.position_count(), 2 * order);
            }
        }
    }
}

fn gauss_legendre(order: usize) -> Vec<(f64, f64)> {
    let polynomial = |x: f64| {
        let (mut previous, mut current) = (0.0, 1.0);
        for degree in 1..=order {
            let next = ((2 * degree - 1) as f64 * x * current - (degree - 1) as f64 * previous)
                / degree as f64;
            previous = current;
            current = next;
        }
        (
            current,
            order as f64 * (x * current - previous) / (x * x - 1.0),
        )
    };
    (0..order)
        .map(|index| {
            let mut x = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
            for _ in 0..24 {
                let (value, derivative) = polynomial(x);
                let change = value / derivative;
                x -= change;
                if change.abs() < 1e-15 {
                    break;
                }
            }
            let (_, derivative) = polynomial(x);
            (x, 2.0 / ((1.0 - x * x) * derivative * derivative))
        })
        .collect()
}

#[test]
fn actual_near_field_energy_flux_integrates_over_the_original_wedge_sector() {
    let nodes = gauss_legendre(48);
    let flow = [0.0, 0.0, 0.55 * medium().sound_speed_m_s()];
    let frequency = 230.0;
    for (order, walls) in [
        (3, [IdealPlaneBoundary::Rigid; 2]),
        (3, [IdealPlaneBoundary::PressureRelease; 2]),
        (
            4,
            [
                IdealPlaneBoundary::Rigid,
                IdealPlaneBoundary::PressureRelease,
            ],
        ),
    ] {
        let opening = PI / order as f64;
        let (positions, q) = fixture(order);
        let walls = planes(order, walls);
        let field = observer(
            frequency,
            flow,
            polar(0.9, opening / 2.0, 0.8),
            &positions,
            walls,
            order,
        );
        let power = PreparedConvectedMonopolePower::from_observer(
            &field,
            0.9,
            budget(positions.len() * 2 * order),
        )
        .unwrap();
        let expected = power.evaluate(&q).unwrap().power_w;
        for radius in [1.4, 3.7] {
            let mut integral = 0.0;
            for &(mu, mu_weight) in &nodes {
                let radial = (1.0 - mu * mu).sqrt();
                for &(phi_node, phi_weight) in &nodes {
                    let phi = 0.5 * opening * (1.0 + phi_node);
                    let normal = [radial * phi.cos(), radial * phi.sin(), mu];
                    let field = observer(
                        frequency,
                        flow,
                        normal.map(|n| radius * n),
                        &positions,
                        walls,
                        order,
                    )
                    .evaluate(&q)
                    .unwrap();
                    let flux: f64 = (0..3)
                        .map(|axis| normal[axis] * field.mean_perturbation_energy_flux_w_m2[axis])
                        .sum();
                    integral += mu_weight * 0.5 * opening * phi_weight * radius * radius * flux;
                }
            }
            close(integral, expected, 3e-10);
        }
    }
}

#[test]
fn stationary_limit_matches_independently_rotated_signed_images() {
    for (order, walls) in [
        (3, [IdealPlaneBoundary::Rigid; 2]),
        (3, [IdealPlaneBoundary::PressureRelease; 2]),
        (
            4,
            [
                IdealPlaneBoundary::Rigid,
                IdealPlaneBoundary::PressureRelease,
            ],
        ),
        (5, [IdealPlaneBoundary::Rigid; 2]),
    ] {
        let (positions, q) = fixture(order);
        let opening = PI / order as f64;
        let count = 2 * order * positions.len();
        let field = observer(
            230.0,
            [0.0; 3],
            polar(0.9, opening / 2.0, 0.8),
            &positions,
            planes(order, walls),
            order,
        );
        let power =
            PreparedConvectedMonopolePower::from_observer(&field, 0.0, budget(count)).unwrap();
        let mut image_positions = Vec::new();
        let mut image_flows = Vec::new();
        for (&position, &flow) in positions.iter().zip(&q) {
            let mut rotation_sign = 1.0;
            for index in 0..order {
                let (sin, cos) = (2.0 * index as f64 * opening).sin_cos();
                for reflected in [false, true] {
                    let y = if reflected { -position[1] } else { position[1] };
                    image_positions.push([
                        cos * position[0] - sin * y,
                        sin * position[0] + cos * y,
                        position[2],
                    ]);
                    let image_sign = rotation_sign * if reflected { sign(walls[0]) } else { 1.0 };
                    image_flows.push(CompactSourcePhasor::monopole(flow.scale(image_sign)));
                }
                rotation_sign *= sign(walls[0]) * sign(walls[1]);
            }
        }
        let stationary =
            PreparedCompactRadiatedPower::new(medium(), 230.0, &image_positions, budget(count))
                .unwrap();
        close(
            power.evaluate(&q).unwrap().power_w,
            stationary.evaluate(&image_flows).unwrap().power_w / (2 * order) as f64,
            5e-12,
        );
    }
}

#[test]
fn affine_wall_anchors_rotation_and_translation_preserve_the_same_wedge() {
    let order = 4;
    let opening = PI / order as f64;
    let (positions, q) = fixture(order);
    let boundaries = [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ];
    let original_planes = planes(order, boundaries);
    let flow = [0.0, 0.0, 0.37 * medium().sound_speed_m_s()];
    let point = polar(0.9, opening / 2.0, 0.8);
    let original = observer(230.0, flow, point, &positions, original_planes, order);
    let expected = PreparedConvectedMonopolePower::from_observer(&original, 0.9, budget(24))
        .unwrap()
        .evaluate(&q)
        .unwrap()
        .power_w;
    let rotate = |v: [f64; 3]| [v[2], -v[0], -v[1]];
    let translate = |v: [f64; 3]| [v[0] + 32.0, v[1] - 16.0, v[2] + 8.0];
    let anchors = [[0.25, 0.0, -0.5], polar(0.625, opening, 0.375)];
    let moved_planes = core::array::from_fn(|index| {
        IdealAcousticPlane::new(
            translate(rotate(anchors[index])),
            rotate(original_planes[index].inward_normal()),
            boundaries[index],
        )
        .unwrap()
    });
    let moved = observer(
        230.0,
        rotate(flow),
        translate(rotate(point)),
        &positions.map(|p| translate(rotate(p))),
        moved_planes,
        order,
    );
    let power = PreparedConvectedMonopolePower::from_observer(&moved, 0.9, budget(24)).unwrap();
    close(power.evaluate(&q).unwrap().power_w, expected, 5e-12);
    assert_eq!(power.acoustic_state(), medium());
    assert_eq!(power.flow_velocity_m_s(), rotate(flow));
    assert_eq!(power.frequency_hz(), 230.0);
    assert_eq!(power.source_count(), positions.len());
}

#[test]
fn complete_wedge_paths_bound_work_and_keep_images_coherent_within_factors() {
    let order = 3;
    let point = polar(0.9, PI / (2 * order) as f64, 0.8);
    let position = polar(0.31, 0.37 * PI / order as f64, 0.17);
    let walls = planes(order, [IdealPlaneBoundary::Rigid; 2]);
    let flow = [0.0, 0.0, 0.55 * medium().sound_speed_m_s()];
    let field = observer(230.0, flow, point, &[position; 2], walls, order);
    for limits in [
        budget(11),
        RadiatedPowerBudget {
            max_sources: 12,
            max_source_pairs: 143,
        },
    ] {
        assert!(PreparedConvectedMonopolePower::from_observer(&field, 0.9, limits).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        PreparedConvectedMonopolePower::from_observer_cancellable(&field, 0.9, budget(12), || {
            polls += 1;
            polls == 8
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let power = PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(12)).unwrap();
    let single_field = observer(230.0, flow, point, &[position], walls, order);
    let single =
        PreparedConvectedMonopolePower::from_observer(&single_field, 0.9, budget(6)).unwrap();
    let q = C64::new(4e-4, 2e-4);
    assert_eq!(power.evaluate(&[q, -q]).unwrap().power_w, 0.0);
    let first = [q, C64::ZERO];
    let second = [C64::ZERO, -q];
    let modes = [&first[..], &second[..]];
    let result = power
        .evaluate_coherence_modes(
            &modes,
            RadiatedPowerModeBudget {
                max_modes: 2,
                max_source_pair_mode_products: 288,
            },
        )
        .unwrap();
    close(
        result.power_w,
        2.0 * single.evaluate(&[q]).unwrap().power_w,
        2e-14,
    );
    assert!(
        power
            .evaluate_coherence_modes(
                &modes,
                RadiatedPowerModeBudget {
                    max_modes: 2,
                    max_source_pair_mode_products: 287
                }
            )
            .is_err()
    );
    assert!(power.evaluate(&[q; 12]).is_err());

    // Exercise the original upper order, including every one of its 256 paths.
    let order = 128;
    let field = observer(
        230.0,
        flow,
        polar(0.9, PI / (2 * order) as f64, 0.8),
        &[polar(0.31, 0.37 * PI / order as f64, 0.17)],
        planes(order, [IdealPlaneBoundary::Rigid; 2]),
        order,
    );
    let power = PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(256)).unwrap();
    assert_eq!(power.position_count(), 256);
    assert!(power.evaluate(&[q]).unwrap().power_w > 0.0);

    let empty_field = observer(230.0, flow, point, &[], walls, 3);
    let empty =
        PreparedConvectedMonopolePower::from_observer(&empty_field, 0.9, budget(0)).unwrap();
    let result = empty.evaluate(&[]).unwrap();
    assert_eq!(result.power_w, 0.0);
    assert_eq!(result.scope, CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE);
}
