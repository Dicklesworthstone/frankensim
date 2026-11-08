//! Power into actual reflecting sectors, including mean-flow energy transport.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{
    CompactSourcePhasor, PreparedCompactRadiatedPower, PreparedConvectedMonopolePower,
    RadiatedPowerBudget, RadiatedPowerModeBudget,
    array::{IdealAcousticPlane, IdealPlaneBoundary},
    array::convected::{ConvectedMonopoleLimits, PreparedConvectedMonopoleObserver},
    power::CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 331.0).unwrap()
}
fn budget(paths: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget { max_sources: paths, max_source_pairs: paths * paths }
}
fn plane(axis: usize, boundary: IdealPlaneBoundary) -> IdealAcousticPlane {
    let mut normal = [0.0; 3];
    normal[axis] = 1.0;
    IdealAcousticPlane::new([0.0; 3], normal, boundary).unwrap()
}
fn sign(boundary: IdealPlaneBoundary) -> f64 {
    match boundary {
        IdealPlaneBoundary::Rigid => 1.0,
        IdealPlaneBoundary::PressureRelease => -1.0,
    }
}
fn observer(
    flow: [f64; 3],
    frequency: f64,
    point: [f64; 3],
    positions: &[[f64; 3]],
    planes: &[IdealAcousticPlane],
) -> PreparedConvectedMonopoleObserver {
    PreparedConvectedMonopoleObserver::new(
        medium(), flow, frequency, point, positions, positions.len(),
        ConvectedMonopoleLimits { maximum_mach: 0.9, minimum_distance_m: 1e-6 },
    ).unwrap().with_orthogonal_planes(planes, positions.len() * (1 << planes.len())).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}");
}
fn unit_power(frequency: f64, mach: f64) -> f64 {
    medium().density_kg_m3() * (TAU * frequency).powi(2)
        / (8.0 * PI * medium().sound_speed_m_s() * (1.0 - mach * mach).powi(2))
}
fn sinc(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

#[test]
fn one_wall_retains_height_interference_and_physical_hemisphere_normalization() {
    let q = C64::new(2e-4, -3e-4);
    for boundary in [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease] {
        for mach in [0.0_f64, 0.4, 0.8] {
            let beta = (1.0 - mach * mach).sqrt();
            let flow = [mach * medium().sound_speed_m_s(), 0.0, 0.0];
            for frequency in [83.0, 317.0] {
                let k = TAU * frequency / medium().sound_speed_m_s();
                for height in [0.031, 0.17, 0.5] {
                    let field = observer(flow, frequency, [0.8, 0.3, 1.0],
                        &[[0.2, -0.1, height]], &[plane(2, boundary)]);
                    let power = PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(2)).unwrap();
                    let result = power.evaluate(&[q]).unwrap();
                    let expected = unit_power(frequency, mach) * q.norm_sq()
                        * (1.0 + sign(boundary) * sinc(2.0 * k * height / beta));
                    close(result.power_w, expected, 5e-13);
                    assert_eq!(result.scope, CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE);
                    assert_eq!(power.source_count(), 1);
                }
            }
        }
    }
}

#[test]
fn two_perpendicular_walls_keep_all_four_image_interference_terms() {
    let mach = 0.55_f64;
    let beta = (1.0 - mach * mach).sqrt();
    let frequency = 184.0;
    let k = TAU * frequency / medium().sound_speed_m_s();
    let heights = [0.11_f64, 0.27];
    let q = C64::new(3e-4, 1e-4);
    for first in [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease] {
        for second in [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease] {
            let planes = [plane(0, first), plane(1, second)];
            let field = observer([0.0, 0.0, mach * medium().sound_speed_m_s()], frequency,
                [0.8, 0.9, 0.7], &[[heights[0], heights[1], -0.2]], &planes);
            let power = PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(4)).unwrap();
            let factor = 1.0
                + sign(first) * sinc(2.0 * k * heights[0] / beta)
                + sign(second) * sinc(2.0 * k * heights[1] / beta)
                + sign(first) * sign(second) * sinc(2.0 * k * heights[0].hypot(heights[1]) / beta);
            close(power.evaluate(&[q]).unwrap().power_w,
                unit_power(frequency, mach) * q.norm_sq() * factor, 1e-12);
        }
    }
}

// Independent quadrature nodes, with no production spherical-moment kernel.
fn gauss_legendre(order: usize) -> Vec<(f64, f64)> {
    let polynomial = |x: f64| {
        let (mut previous, mut current) = (0.0, 1.0);
        for degree in 1..=order {
            let next = ((2 * degree - 1) as f64 * x * current
                - (degree - 1) as f64 * previous) / degree as f64;
            previous = current;
            current = next;
        }
        (current, order as f64 * (x * current - previous) / (x * x - 1.0))
    };
    (0..order).map(|index| {
        let mut x = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..24 {
            let (value, derivative) = polynomial(x);
            let step = value / derivative;
            x -= step;
            if step.abs() < 1e-15 { break; }
        }
        let (_, derivative) = polynomial(x);
        (x, 2.0 / ((1.0 - x * x) * derivative * derivative))
    }).collect()
}

#[test]
fn actual_near_field_flux_integrates_over_the_hemisphere_and_right_angle_sector() {
    let positions = [[0.13, 0.21, 0.19], [0.37, -0.2, 0.28], [0.08, 0.15, 0.35]];
    let q = [C64::new(4e-4, 2e-4), C64::new(-1e-4, 3e-4), C64::new(2e-4, -1e-4)];
    let frequency = 230.0;
    let polar = gauss_legendre(48);
    let azimuth = gauss_legendre(96);
    for corner in [false, true] {
        let (planes, flow, phi_mid, phi_half) = if corner {
            (vec![plane(2, IdealPlaneBoundary::PressureRelease), plane(0, IdealPlaneBoundary::Rigid)],
                [0.0, 0.6 * medium().sound_speed_m_s(), 0.0], 0.0, PI / 2.0)
        } else {
            (vec![plane(2, IdealPlaneBoundary::Rigid)],
                [0.3 * medium().sound_speed_m_s(), 0.4 * medium().sound_speed_m_s(), 0.0], PI, PI)
        };
        let field = observer(flow, frequency, [0.8, 0.9, 0.7], &positions, &planes);
        let power = PreparedConvectedMonopolePower::from_observer(
            &field, 0.9, budget(positions.len() * (1 << planes.len())),
        ).unwrap();
        let expected = power.evaluate(&q).unwrap().power_w;
        for radius in [1.4, 3.7] {
            let mut integral = 0.0;
            for &(mu_node, mu_weight) in &polar {
                let mu = 0.5 * (1.0 + mu_node);
                let radial = (1.0 - mu * mu).sqrt();
                for &(phi_node, phi_weight) in &azimuth {
                    let phi = phi_mid + phi_half * phi_node;
                    let normal = [radial * phi.cos(), radial * phi.sin(), mu];
                    let field = observer(flow, frequency, normal.map(|n| radius * n), &positions, &planes)
                        .evaluate(&q).unwrap();
                    let flux: f64 = (0..3)
                        .map(|axis| normal[axis] * field.mean_perturbation_energy_flux_w_m2[axis]).sum();
                    integral += 0.5 * mu_weight * phi_half * phi_weight * radius * radius * flux;
                }
            }
            close(integral, expected, 2e-10);
        }
    }
}

#[test]
fn zero_flow_matches_stationary_power_of_the_independently_reflected_array() {
    let positions = [[0.13, 0.21, 0.19], [0.37, 0.12, 0.28]];
    let q = [C64::new(4e-4, 2e-4), C64::new(-1e-4, 3e-4)];
    let boundaries = [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease, IdealPlaneBoundary::Rigid];
    for count in 1..=3 {
        let planes: Vec<_> = (0..count).map(|axis| plane(axis, boundaries[axis])).collect();
        let field = observer([0.0; 3], 230.0, [0.8, 0.9, 0.7], &positions, &planes);
        let groups = 1 << count;
        let power = PreparedConvectedMonopolePower::from_observer(&field, 0.0, budget(2 * groups)).unwrap();
        // Sign changes of Cartesian coordinates are independent of the native
        // plane reflection helper. The stationary owner retains all Q phases.
        let mut image_positions = Vec::new();
        let mut image_q = Vec::new();
        for (&position, &flow) in positions.iter().zip(&q) {
            for mask in 0..groups {
                let mut image = position;
                let mut image_sign = 1.0;
                for axis in 0..count {
                    if mask & (1 << axis) != 0 {
                        image[axis] = -image[axis];
                        image_sign *= sign(boundaries[axis]);
                    }
                }
                image_positions.push(image);
                image_q.push(CompactSourcePhasor::monopole(flow.scale(image_sign)));
            }
        }
        let stationary = PreparedCompactRadiatedPower::new(medium(), 230.0, &image_positions, budget(2 * groups)).unwrap();
        close(power.evaluate(&q).unwrap().power_w,
            stationary.evaluate(&image_q).unwrap().power_w / groups as f64, 2e-13);
    }
}

#[test]
fn original_sources_keep_coherent_images_within_each_statistical_factor() {
    let flow = [0.0, 0.0, 0.55 * medium().sound_speed_m_s()];
    let position = [0.13, 0.21, 0.19];
    let planes = [plane(0, IdealPlaneBoundary::Rigid), plane(1, IdealPlaneBoundary::PressureRelease)];
    let q = C64::new(4e-4, 2e-4);
    let field = observer(flow, 230.0, [0.8, 0.9, 0.7], &[position; 2], &planes);
    let power = PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(8)).unwrap();
    let single_field = observer(flow, 230.0, [0.8, 0.9, 0.7], &[position], &planes);
    let single = PreparedConvectedMonopolePower::from_observer(&single_field, 0.9, budget(4)).unwrap();
    assert_eq!(power.evaluate(&[q, -q]).unwrap().power_w, 0.0);
    let first = [q, C64::ZERO];
    let second = [C64::ZERO, -q];
    let modes = power.evaluate_coherence_modes(&[&first, &second],
        RadiatedPowerModeBudget { max_modes: 2, max_source_pair_mode_products: 128 }).unwrap();
    close(modes.power_w, 2.0 * single.evaluate(&[q]).unwrap().power_w, 2e-14);
    assert_eq!(modes.scope, CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE);
    assert!(power.evaluate_coherence_modes(&[&first, &second],
        RadiatedPowerModeBudget { max_modes: 2, max_source_pair_mode_products: 127 }).is_err());
    assert!(power.evaluate(&[q; 8]).is_err());
}

#[test]
fn complete_image_budgets_cancellation_and_empty_domain_scope_are_retained() {
    let planes = [plane(0, IdealPlaneBoundary::Rigid), plane(1, IdealPlaneBoundary::PressureRelease)];
    let flow = [0.0, 0.0, 0.55 * medium().sound_speed_m_s()];
    let field = observer(flow, 230.0, [0.8, 0.9, 0.7], &[[0.13, 0.21, 0.19]; 2], &planes);
    for limits in [budget(7), RadiatedPowerBudget { max_sources: 8, max_source_pairs: 63 }] {
        assert!(PreparedConvectedMonopolePower::from_observer(&field, 0.9, limits).is_err());
    }
    let mut polls = 0;
    assert!(matches!(PreparedConvectedMonopolePower::from_observer_cancellable(&field, 0.9, budget(8), || {
        polls += 1;
        polls == 6
    }), Err(AeroacError::Cancelled { .. })));
    let power = PreparedConvectedMonopolePower::from_observer(&field, 0.9, budget(8)).unwrap();
    assert!(matches!(power.evaluate_cancellable(&[C64::ONE; 2], || true), Err(AeroacError::Cancelled { .. })));
    assert!(power.evaluate(&[C64::new(f64::NAN, 0.0); 2]).is_err());
    assert!(power.evaluate(&[C64::ONE; 2]).unwrap().power_w > 0.0);

    let empty_field = observer(flow, 230.0, [0.8, 0.9, 0.7], &[], &planes);
    let empty = PreparedConvectedMonopolePower::from_observer(&empty_field, 0.9, budget(0)).unwrap();
    let result = empty.evaluate(&[]).unwrap();
    assert_eq!(result.power_w, 0.0);
    assert_eq!(result.scope, CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE);
    let modes = empty.evaluate_coherence_modes(&[],
        RadiatedPowerModeBudget { max_modes: 0, max_source_pair_mode_products: 0 }).unwrap();
    assert_eq!(modes, result);
}
