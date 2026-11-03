use super::*;

const RADIUS: f64 = 2.0;
const GAP: f64 = 1.0 / 128.0;
const CUTOFF: f64 = 1.0 / 32.0;
const FLOOR: f64 = 1.0 / 512.0;
const MU: f64 = 0.007;

fn plane(velocity: [f64; 2]) -> LubricationPlane2 {
    LubricationPlane2::new([0.0; 2], [0.0, 1.0], velocity, 0.125, 1.0 / 1024.0)
        .unwrap()
        .with_tangential(CUTOFF, FLOOR)
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

#[test]
fn cylinder_slide_spin_and_rolling_have_the_independent_diagonal_stokes_law() {
    let center = [0.0, RADIUS + GAP];
    let origin = [-11.0, 7.0];
    // Independent leading-order closed form, without the production
    // rationalized cutoff arithmetic or its force/couple construction.
    let zeta =
        core::f64::consts::PI * MU * (8.0 * RADIUS).sqrt() * (GAP.powf(-0.5) - CUTOFF.powf(-0.5));
    for (u, omega, wall_u) in [
        (0.004, 0.0, 0.0),    // Pure sliding: no center torque.
        (0.0, 0.002, 0.0),    // Pure spinning: no net force.
        (-0.004, 0.002, 0.0), // Zero contact-point slip still dissipates.
        (0.006, 0.0, 0.006),  // Common translation: no relative film motion.
    ] {
        let prepared = PreparedPlane2::new(plane([wall_u, 0.0]), center, RADIUS, MU).unwrap();
        let receipt = prepared.exchange([u, 0.0], omega, origin).unwrap();
        let force = -zeta * (u - wall_u);
        let torque = -zeta * RADIUS * RADIUS * omega;
        close(receipt.tangential_resistance_increment, zeta);
        close(receipt.body_impulse[0], force);
        assert_eq!(receipt.body_impulse[1], 0.0);
        close(receipt.body_center_couple, torque);
        close(
            receipt.body_impulse[2],
            -(center[1] - origin[1]) * force + torque,
        );
        for k in 0..3 {
            assert_eq!(receipt.plane_impulse[k], -receipt.body_impulse[k]);
        }
        close(
            receipt.plane_impulse[2],
            cross(
                subtract(receipt.plane_footpoint, origin),
                [receipt.plane_impulse[0], receipt.plane_impulse[1]],
            ) + receipt.plane_couple,
        );
        close(receipt.body_work, force * u + torque * omega);
        close(receipt.plane_work, -force * wall_u);
        close(
            receipt.dissipated_energy,
            zeta * ((u - wall_u).powi(2) + (RADIUS * omega).powi(2)),
        );
        close(
            receipt.body_work + receipt.plane_work,
            -receipt.dissipated_energy,
        );
        if omega == 0.0 {
            assert_eq!(receipt.body_center_couple, 0.0);
        } else if u == 0.0 {
            assert_eq!(receipt.body_impulse[0], 0.0);
        }
        if u + RADIUS * omega == wall_u && omega != 0.0 {
            assert!(
                receipt.dissipated_energy > 0.0,
                "rolling is not a zero-resistance mode"
            );
        }
        if u == wall_u && omega == 0.0 {
            assert_eq!(receipt.body_impulse, [0.0; 3]);
            assert_eq!(receipt.dissipated_energy, 0.0);
        }
    }
}

#[test]
fn film_force_work_and_reaction_couple_transform_with_origin_and_common_velocity() {
    let center = [0.0, RADIUS + GAP];
    let origin = [-11.0, 7.0];
    let offset = [13.0, -19.0];
    let prepared = PreparedPlane2::new(plane([0.003, 0.0]), center, RADIUS, MU).unwrap();
    let first = prepared.exchange([0.007, 0.0], 0.002, origin).unwrap();
    let moved = prepared
        .exchange(
            [0.007, 0.0],
            0.002,
            [origin[0] + offset[0], origin[1] + offset[1]],
        )
        .unwrap();
    let boosted = PreparedPlane2::new(plane([0.013, 0.0]), center, RADIUS, MU)
        .unwrap()
        .exchange([0.017, 0.0], 0.002, origin)
        .unwrap();
    close(
        moved.body_impulse[2],
        first.body_impulse[2] - cross(offset, [first.body_impulse[0], first.body_impulse[1]]),
    );
    close(
        moved.plane_impulse[2],
        first.plane_impulse[2] - cross(offset, [first.plane_impulse[0], first.plane_impulse[1]]),
    );
    assert_eq!(moved.body_center_couple, first.body_center_couple);
    assert_eq!(moved.plane_couple, first.plane_couple);
    assert_eq!(moved.body_work, first.body_work);
    for k in 0..3 {
        close(boosted.body_impulse[k], first.body_impulse[k]);
    }
    close(
        boosted.body_work - first.body_work,
        0.01 * first.body_impulse[0],
    );
    close(
        boosted.plane_work - first.plane_work,
        0.01 * first.plane_impulse[0],
    );
    close(boosted.dissipated_energy, first.dissipated_energy);
}

#[test]
fn tangent_cutoff_floor_and_refusals_are_separate_from_the_normal_law() {
    let wall = plane([0.0; 2]);
    let below_floor = [1.0 / 1024.0, 3.0 / 2048.0];
    let first = PreparedPlane2::new(wall, [0.0, RADIUS + below_floor[0]], RADIUS, MU).unwrap();
    let second = PreparedPlane2::new(wall, [0.0, RADIUS + below_floor[1]], RADIUS, MU).unwrap();
    assert_eq!(first.tangential_resistance, second.tangential_resistance);
    assert_ne!(first.pair.gap(), second.pair.gap());
    for gap in [CUTOFF, 2.0 * CUTOFF] {
        let prepared = PreparedPlane2::new(wall, [0.0, RADIUS + gap], RADIUS, MU).unwrap();
        let receipt = prepared.exchange([0.004, 0.001], 0.002, [0.0; 2]).unwrap();
        assert_eq!(receipt.tangential_resistance_increment, 0.0);
        assert_eq!(receipt.body_center_couple, 0.0);
        assert!(receipt.native.resistance_increment > 0.0);
        assert_eq!(
            receipt.body_impulse,
            [
                receipt.native.first_solid_impulse[0],
                receipt.native.first_solid_impulse[1],
                receipt.native.first_solid_angular_impulse,
            ]
        );
    }
    let just_inside = f64::from_bits(CUTOFF.to_bits() - 1);
    assert!(wall.tangential_resistance(RADIUS, just_inside, MU).unwrap() > 0.0);
    for (cutoff, floor) in [(f64::NAN, FLOOR), (CUTOFF, 0.0), (CUTOFF, 2.0 * CUTOFF)] {
        assert!(wall.with_tangential(cutoff, floor).is_err());
    }
    assert!(
        PreparedPlane2::new(
            wall.with_tangential(0.1, FLOOR).unwrap(),
            [0.0, RADIUS + GAP],
            RADIUS,
            MU,
        )
        .is_err()
    );
    assert!(PreparedPlane2::new(wall, [0.0, RADIUS + GAP], RADIUS, f64::MAX).is_err());
    assert!(first.exchange([0.0; 2], f64::NAN, [0.0; 2]).is_err());
    assert!(first.exchange([0.0; 2], f64::MAX, [0.0; 2]).is_err());
}
