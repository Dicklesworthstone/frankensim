//! Physical continuous intake on the original native variable-mass state.

use fs_mbd::Vec3;
use fs_mbd::rocket::exhaust::{ParticleCaptureForcing, ParticleCaptureStep};
use fs_mbd::rocket::{RocketError, RocketState};

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 4.0e-13 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn close_vec(actual: Vec3, expected: Vec3, scale: f64) {
    close(actual.x, expected.x, scale);
    close(actual.y, expected.y, scale);
    close(actual.z, expected.z, scale);
}

fn state() -> RocketState {
    RocketState {
        position_m: Vec3::new(1.0, 2.0, -3.0),
        velocity_m_s: Vec3::new(12.0, -6.0, 3.0),
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1.0,
    }
}

fn stream() -> ParticleCaptureForcing {
    ParticleCaptureForcing {
        mass_flow_kg_s: 2.0,
        velocity_world_m_s: Vec3::new(-2.0, 4.0, 1.0),
        internal_energy_rate_w: -7.0,
    }
}

fn balances(step: ParticleCaptureStep, scale: f64) {
    close(step.mass_residual_kg, 0.0, scale);
    close_vec(step.momentum_residual_kg_m_s, Vec3::ZERO, scale);
    close_vec(step.angular_momentum_residual_kg_m2_s, Vec3::ZERO, scale);
    close(step.kinetic_energy_residual_j, 0.0, scale);
}

#[test]
fn stationary_capture_preserves_momentum_and_thermalizes_exact_ram_drag() {
    let initial = state();
    let result = initial
        .advance_with_inertial_capture(
            ParticleCaptureForcing {
                velocity_world_m_s: Vec3::ZERO,
                ..stream()
            },
            1.5,
            || false,
        )
        .unwrap();
    assert_eq!(result.captured_mass_kg, 3.0);
    assert_eq!(result.state_after.propellant_mass_kg, 4.0);
    close_vec(
        result.state_after.velocity_m_s,
        initial.velocity_m_s.scale(2.0 / 5.0),
        12.0,
    );
    // Integral of v0*M0/(M0+q*t), independently in terms of physical masses.
    close_vec(
        result.state_after.position_m,
        initial
            .position_m
            .add(initial.velocity_m_s.scale(2.5_f64.ln())),
        20.0,
    );
    let heat = 0.5 * (2.0 * 3.0 / 5.0) * initial.velocity_m_s.norm_squared();
    close(result.thermalization_energy_j, heat, heat);
    close(result.kinetic_energy_change_j, -heat, heat);
    assert_eq!(result.incoming_momentum_world_kg_m_s, Vec3::ZERO);
    assert_eq!(result.incoming_angular_momentum_world_kg_m2_s, Vec3::ZERO);
    assert_eq!(result.incoming_energy_world_j, -10.5);
    close(result.retained_internal_energy_j, heat - 10.5, heat);
    balances(result, 200.0);
    assert_eq!(initial, state());
}

#[test]
fn transverse_intake_carries_orbital_momentum_at_the_moving_capture_point() {
    let initial = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::new(1.0, 0.0, 0.0),
        ..state()
    };
    let result = initial
        .advance_with_inertial_capture(
            ParticleCaptureForcing {
                mass_flow_kg_s: 1.0,
                velocity_world_m_s: Vec3::new(0.0, 1.0, 0.0),
                internal_energy_rate_w: 0.0,
            },
            2.0,
            || false,
        )
        .unwrap();
    // x(t)=2*ln(1+t/2); q*integral_0^2 x(t) dt = 4*(2*ln(2)-1).
    // An endpoint-only flux would instead be either zero or 4*ln(2).
    let angular = 4.0 * (2.0 * 2.0_f64.ln() - 1.0);
    close_vec(
        result.incoming_angular_momentum_world_kg_m2_s,
        Vec3::new(0.0, 0.0, angular),
        4.0,
    );
    close_vec(
        result.state_after.velocity_m_s,
        Vec3::new(0.5, 0.5, 0.0),
        1.0,
    );
    close(result.thermalization_energy_j, 1.0, 1.0);
    balances(result, 4.0);
}

#[test]
fn exact_capture_composes_over_time_with_the_same_physical_fluxes() {
    let initial = state();
    let whole = initial
        .advance_with_inertial_capture(stream(), 3.0, || false)
        .unwrap();
    let mut next = initial;
    let mut heat = 0.0;
    let mut incoming_energy = 0.0;
    let mut incoming_angular = Vec3::ZERO;
    for _ in 0..12 {
        let part = next
            .advance_with_inertial_capture(stream(), 0.25, || false)
            .unwrap();
        next = part.state_after;
        heat += part.thermalization_energy_j;
        incoming_energy += part.incoming_energy_world_j;
        incoming_angular = incoming_angular.add(part.incoming_angular_momentum_world_kg_m2_s);
        balances(part, 400.0);
    }
    close_vec(next.position_m, whole.state_after.position_m, 30.0);
    close_vec(next.velocity_m_s, whole.state_after.velocity_m_s, 12.0);
    assert_eq!(
        next.propellant_mass_kg,
        whole.state_after.propellant_mass_kg
    );
    close(heat, whole.thermalization_energy_j, 400.0);
    close(incoming_energy, whole.incoming_energy_world_j, 400.0);
    close_vec(
        incoming_angular,
        whole.incoming_angular_momentum_world_kg_m2_s,
        400.0,
    );
    balances(whole, 400.0);
}

#[test]
fn capture_is_covariant_under_boost_rotation_and_origin_shift() {
    let initial = state();
    let forcing = stream();
    let h = 1.5;
    let original = initial
        .advance_with_inertial_capture(forcing, h, || false)
        .unwrap();
    // An exact proper rotation permits a strict independent frame test.
    let rotate = |v: Vec3| Vec3::new(v.z, v.x, v.y);
    let boost = Vec3::new(7.0, -3.0, 2.0);
    let origin = Vec3::new(-4.0, 1.0, 5.0);
    let transformed = RocketState {
        position_m: rotate(initial.position_m).add(origin),
        velocity_m_s: rotate(initial.velocity_m_s).add(boost),
        ..initial
    }
    .advance_with_inertial_capture(
        ParticleCaptureForcing {
            velocity_world_m_s: rotate(forcing.velocity_world_m_s).add(boost),
            ..forcing
        },
        h,
        || false,
    )
    .unwrap();
    close_vec(
        transformed.state_after.position_m,
        rotate(original.state_after.position_m)
            .add(origin)
            .add(boost.scale(h)),
        40.0,
    );
    close_vec(
        transformed.state_after.velocity_m_s,
        rotate(original.state_after.velocity_m_s).add(boost),
        20.0,
    );
    close_vec(
        transformed.incoming_momentum_world_kg_m_s,
        rotate(original.incoming_momentum_world_kg_m_s).add(boost.scale(original.captured_mass_kg)),
        40.0,
    );
    close(
        transformed.incoming_energy_world_j,
        original.incoming_energy_world_j
            + boost.dot(rotate(original.incoming_momentum_world_kg_m_s))
            + 0.5 * original.captured_mass_kg * boost.norm_squared(),
        400.0,
    );
    close(
        transformed.thermalization_energy_j,
        original.thermalization_energy_j,
        400.0,
    );
    close(
        transformed.retained_internal_energy_j,
        original.retained_internal_energy_j,
        400.0,
    );
    balances(transformed, 400.0);
}

#[test]
fn small_and_dominant_capture_retain_their_physical_limits() {
    let initial = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
        dry_mass_kg: 2.0,
        propellant_mass_kg: 0.0,
    };
    let dm = 2.0_f64.powi(-40);
    let small = initial
        .advance_with_inertial_capture(
            ParticleCaptureForcing {
                mass_flow_kg_s: dm,
                velocity_world_m_s: Vec3::new(2.0, 0.0, 0.0),
                internal_energy_rate_w: 0.0,
            },
            1.0,
            || false,
        )
        .unwrap();
    let fraction = dm / (2.0 + dm);
    close(small.state_after.velocity_m_s.x, 2.0 * fraction, dm);
    // C=1/2+x/6+O(x²); an uncorrected average-mass step misses this term.
    close(
        small.state_after.position_m.x,
        fraction * (1.0 + fraction / 3.0),
        dm,
    );
    close(small.thermalization_energy_j, 4.0 * fraction, dm);
    balances(small, dm);

    let moving = RocketState {
        velocity_m_s: Vec3::new(1.25, 0.0, 0.0),
        ..initial
    };
    let dominant = moving
        .advance_with_inertial_capture(
            ParticleCaptureForcing {
                mass_flow_kg_s: 1.0e20,
                velocity_world_m_s: Vec3::ZERO,
                internal_energy_rate_w: 0.0,
            },
            1.0,
            || false,
        )
        .unwrap();
    close(dominant.state_after.velocity_m_s.x, 2.5e-20, 2.5e-20);
    close(
        dominant.state_after.position_m.x,
        2.5e-20 * (5.0e19_f64).ln(),
        1.0e-18,
    );
    close(dominant.thermalization_energy_j, 1.5625, 2.0);
    close(dominant.kinetic_energy_change_j, -1.5625, 2.0);
    balances(dominant, 2.0);
}

#[test]
fn zero_rate_coasts_zero_time_is_identity_and_lost_recoil_stays_visible() {
    let initial = state();
    let zero = initial
        .advance_with_inertial_capture(stream(), 0.0, || false)
        .unwrap();
    assert_eq!(zero.state_after, initial);
    assert_eq!(zero.thermalization_energy_j, 0.0);
    let coast = initial
        .advance_with_inertial_capture(
            ParticleCaptureForcing {
                mass_flow_kg_s: 0.0,
                internal_energy_rate_w: 0.0,
                ..stream()
            },
            2.0,
            || false,
        )
        .unwrap();
    assert_eq!(coast.state_after.velocity_m_s, initial.velocity_m_s);
    assert_eq!(
        coast.state_after.position_m,
        initial.position_m.add(initial.velocity_m_s.scale(2.0))
    );
    balances(coast, 1.0);

    let boosted = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::new(1.0e16, 0.0, 0.0),
        ..initial
    };
    let unresolved = boosted
        .advance_with_inertial_capture(
            ParticleCaptureForcing {
                mass_flow_kg_s: 1.0,
                velocity_world_m_s: Vec3::new(1.0e16 + 2.0, 0.0, 0.0),
                internal_energy_rate_w: 0.0,
            },
            1.0,
            || false,
        )
        .unwrap();
    assert_eq!(unresolved.state_after.velocity_m_s, boosted.velocity_m_s);
    assert_eq!(unresolved.momentum_residual_kg_m_s.x, -2.0);
}

#[test]
fn invalid_capture_and_late_cancellation_publish_no_endpoint() {
    let initial = state();
    for forcing in [
        ParticleCaptureForcing {
            mass_flow_kg_s: -1.0,
            ..stream()
        },
        ParticleCaptureForcing {
            mass_flow_kg_s: 0.0,
            ..stream()
        },
        ParticleCaptureForcing {
            velocity_world_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
            ..stream()
        },
    ] {
        assert!(
            initial
                .advance_with_inertial_capture(forcing, 1.0, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    initial
        .advance_with_inertial_capture(stream(), 1.0, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut observed = 0;
    assert_eq!(
        initial.advance_with_inertial_capture(stream(), 1.0, || {
            observed += 1;
            observed == polls
        }),
        Err(RocketError::Cancelled)
    );
    assert_eq!(initial, state());
}
