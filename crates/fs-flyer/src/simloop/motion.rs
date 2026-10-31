//! Longitudinal momentum in the inertial flight-line/up frame.
//!
//! Body air velocity is an observation of ground momentum, attitude and wind.
//! Integrating its rotating-frame transport with Euler creates kinetic energy
//! even at zero force. Advance ground velocity instead, then rotate back.

use fs_math::det;

pub(super) const GRAVITY_MPS2: f64 = 9.80665;

#[derive(Clone, Copy, Debug)]
pub(super) struct LongitudinalMotion {
    pub position_m: [f64; 2],
    pub ground_mps: [f64; 2],
    pub theta_rad: f64,
    pub q_rad_s: f64,
}

/// The coupled airborne degrees of freedom. Wind and controls are sampled
/// inputs; a stage never advances their clocks or mutates accepted state.
#[derive(Clone, Copy, Debug)]
pub(super) struct FlightMotion {
    pub longitudinal: LongitudinalMotion,
    pub rotor_rad_s: f64,
}

impl FlightMotion {
    pub fn is_finite(self) -> bool {
        self.longitudinal
            .position_m
            .iter()
            .chain(&self.longitudinal.ground_mps)
            .chain(&[
                self.longitudinal.theta_rad,
                self.longitudinal.q_rad_s,
                self.rotor_rad_s,
            ])
            .all(|value| value.is_finite())
    }
}

/// Physical derivatives obtained from one native aerodynamic/rotor solve.
#[derive(Clone, Copy, Debug)]
pub(super) struct FlightAcceleration {
    pub world_mps2: [f64; 2],
    pub pitch_rad_s2: f64,
    pub rotor_rad_s2: f64,
}

/// Drift to the force point, then predict velocity-dependent loads there.
/// Pairing this half drift with a half drift at the accepted velocity makes
/// the orbital angular impulse exactly `dt * x_mid cross F_mid`, to roundoff.
/// The pitch predictor retains its quadratic constant-acceleration term.
pub(super) fn midpoint_stage(
    start: FlightMotion,
    acceleration: FlightAcceleration,
    dt: f64,
) -> FlightMotion {
    let half = 0.5 * dt;
    let mut stage = start;
    for axis in 0..2 {
        stage.longitudinal.position_m[axis] += half * start.longitudinal.ground_mps[axis];
        stage.longitudinal.ground_mps[axis] += half * acceleration.world_mps2[axis];
    }
    stage.longitudinal.theta_rad +=
        half * (start.longitudinal.q_rad_s + 0.5 * half * acceleration.pitch_rad_s2);
    stage.longitudinal.q_rad_s += half * acceleration.pitch_rad_s2;
    stage.rotor_rad_s += half * acceleration.rotor_rad_s2;
    stage
}

/// Accept the physical midpoint impulse and midpoint kinematics. The force
/// callback is owned by the lifecycle; prediction itself publishes nothing.
/// Smooth longitudinal/rotor dynamics are second order with controls and the
/// temporal wind realization held over the tick. Discrete control/OU updates,
/// envelope clipping and terminal contact events have no smooth-order claim.
pub(super) fn finish_midpoint(
    start: FlightMotion,
    stage: FlightMotion,
    acceleration: FlightAcceleration,
    dt: f64,
) -> FlightMotion {
    let mut next = start;
    for axis in 0..2 {
        next.longitudinal.ground_mps[axis] += dt * acceleration.world_mps2[axis];
        next.longitudinal.position_m[axis] =
            stage.longitudinal.position_m[axis] + 0.5 * dt * next.longitudinal.ground_mps[axis];
    }
    next.longitudinal.q_rad_s += dt * acceleration.pitch_rad_s2;
    next.longitudinal.theta_rad +=
        0.5 * dt * (start.longitudinal.q_rad_s + next.longitudinal.q_rad_s);
    next.rotor_rad_s += dt * acceleration.rotor_rad_s2;
    next
}

/// Body FRD air velocity to inertial [forward, upward] ground velocity.
pub(super) fn ground_velocity([u, w]: [f64; 2], theta: f64, wind_mps: [f64; 2]) -> [f64; 2] {
    let (s, c) = (det::sin(theta), det::cos(theta));
    [c * u + s * w + wind_mps[0], s * u - c * w + wind_mps[1]]
}

/// Inertial ground velocity to body FRD air velocity at the sampled wind.
pub(crate) fn body_air_velocity(ground_mps: [f64; 2], theta: f64, wind_mps: [f64; 2]) -> [f64; 2] {
    let vx = ground_mps[0] - wind_mps[0];
    let vh = ground_mps[1] - wind_mps[1];
    let (s, c) = (det::sin(theta), det::cos(theta));
    [c * vx + s * vh, s * vx - c * vh]
}

/// Remove force_buildup's level-flight (theta = alpha) gravity contribution.
/// The inertial integrators apply actual world gravity once, independently of
/// attitude, angle of attack and the scalar gust's alpha increment.
pub(super) fn nongravity_acceleration(
    force_body_n: [f64; 2],
    alpha: f64,
    mass_kg: f64,
) -> [f64; 2] {
    let weight = mass_kg * GRAVITY_MPS2;
    [
        (force_body_n[0] + weight * det::sin(alpha)) / mass_kg,
        (force_body_n[1] - weight * det::cos(alpha)) / mass_kg,
    ]
}

pub(super) fn world_acceleration(
    acceleration_body: [f64; 2],
    theta: f64,
    gravity: f64,
) -> [f64; 2] {
    let (s, c) = (det::sin(theta), det::cos(theta));
    [
        c * acceleration_body[0] + s * acceleration_body[1],
        s * acceleration_body[0] - c * acceleration_body[1] - gravity,
    ]
}

/// Advance the accepted body force and pitch moment held over this tick.
/// Constant pitch acceleration gives the exact midpoint/end attitudes. Rotate
/// the held non-gravity body acceleration at that midpoint, apply its world
/// impulse and drift the centre of mass with the mean ground velocity.
///
/// This is the held-load terminal-event predictor and analytic baseline.
/// Smooth airborne evolution instead resamples loads at `midpoint_stage`.
pub(super) fn held_load_step(
    state: LongitudinalMotion,
    acceleration_body: [f64; 2],
    pitch_accel: f64,
    gravity: f64,
    dt: f64,
) -> LongitudinalMotion {
    let half_dt = 0.5 * dt;
    let q_mid = state.q_rad_s + half_dt * pitch_accel;
    let theta_mid = state.theta_rad + half_dt * (state.q_rad_s + 0.25 * dt * pitch_accel);
    let acceleration_world = world_acceleration(acceleration_body, theta_mid, gravity);
    let mut next = state;
    for (axis, acceleration) in acceleration_world.into_iter().enumerate() {
        next.ground_mps[axis] = state.ground_mps[axis] + dt * acceleration;
        next.position_m[axis] =
            state.position_m[axis] + dt * (state.ground_mps[axis] + half_dt * acceleration);
    }
    next.theta_rad = state.theta_rad + dt * q_mid;
    next.q_rad_s = state.q_rad_s + dt * pitch_accel;
    next
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, error={:.3e}",
            actual - expected
        );
    }

    fn initial() -> LongitudinalMotion {
        LongitudinalMotion {
            position_m: [1.3, 30.0],
            ground_mps: [12.4, -1.8],
            theta_rad: 0.37,
            q_rad_s: 1.4,
        }
    }

    #[test]
    fn rotating_force_free_body_preserves_world_momentum_and_energy() {
        let start = initial();
        let wind = [-5.0, 0.0];
        let mut state = start;
        let mut air = body_air_velocity(state.ground_mps, state.theta_rad, wind);
        let air_energy = air[0] * air[0] + air[1] * air[1];
        let dt = 0.0025;
        for _ in 0..4_000 {
            // Exercise the ordinary lifecycle's air-state storage round trip.
            state.ground_mps = ground_velocity(air, state.theta_rad, wind);
            state = held_load_step(state, [0.0; 2], 0.0, 0.0, dt);
            air = body_air_velocity(state.ground_mps, state.theta_rad, wind);
        }
        for axis in 0..2 {
            close(state.ground_mps[axis], start.ground_mps[axis], 2e-10);
            close(
                state.position_m[axis],
                start.position_m[axis] + 10.0 * start.ground_mps[axis],
                2e-9,
            );
        }
        close(air[0] * air[0] + air[1] * air[1], air_energy, 2e-9);
        // Independent analytic rotation, with platform trig only in the oracle.
        let (s, c) = (start.theta_rad + 10.0 * start.q_rad_s).sin_cos();
        close(
            air[0],
            c * (start.ground_mps[0] - wind[0]) + s * start.ground_mps[1],
            2e-10,
        );
        close(
            air[1],
            s * (start.ground_mps[0] - wind[0]) - c * start.ground_mps[1],
            2e-10,
        );
    }

    #[test]
    fn gravity_is_ballistic_when_attitude_and_incidence_differ() {
        let mass = 341.0;
        for alpha in [-0.4, 0.21] {
            let weight = mass * GRAVITY_MPS2;
            let acceleration = nongravity_acceleration(
                [-weight * det::sin(alpha), weight * det::cos(alpha)],
                alpha,
                mass,
            );
            assert_eq!(acceleration, [0.0; 2]);
            let start = initial();
            let dt = 0.3;
            let next = held_load_step(start, acceleration, -0.17, GRAVITY_MPS2, dt);
            assert_eq!(next.ground_mps[0], start.ground_mps[0]);
            close(
                next.ground_mps[1],
                start.ground_mps[1] - GRAVITY_MPS2 * dt,
                1e-14,
            );
            close(
                next.position_m[0],
                start.position_m[0] + start.ground_mps[0] * dt,
                1e-14,
            );
            close(
                next.position_m[1],
                start.position_m[1] + start.ground_mps[1] * dt - 0.5 * GRAVITY_MPS2 * dt * dt,
                1e-14,
            );
        }
    }

    #[test]
    fn sampled_force_delivers_world_impulse_and_midpoint_com_drift() {
        let mass = 257.0;
        let alpha = -0.2;
        let force = [761.0, -413.0];
        let weight = mass * GRAVITY_MPS2;
        let acceleration = nongravity_acceleration(
            [
                force[0] - weight * det::sin(alpha),
                force[1] + weight * det::cos(alpha),
            ],
            alpha,
            mass,
        );
        let start = initial();
        let dt = 0.18;
        let pitch_accel = 0.53;
        let next = held_load_step(start, acceleration, pitch_accel, GRAVITY_MPS2, dt);
        let middle_angle = start.theta_rad + start.q_rad_s * dt / 2.0 + pitch_accel * dt * dt / 8.0;
        let (s, c) = middle_angle.sin_cos();
        let world_force = [
            c * force[0] + s * force[1],
            s * force[0] - c * force[1] - weight,
        ];
        for axis in 0..2 {
            close(
                mass * (next.ground_mps[axis] - start.ground_mps[axis]),
                dt * world_force[axis],
                1e-11,
            );
            close(
                next.position_m[axis],
                start.position_m[axis]
                    + dt * start.ground_mps[axis]
                    + 0.5 * dt * dt * world_force[axis] / mass,
                1e-14,
            );
        }
    }

    #[test]
    fn constant_pitch_acceleration_has_exact_angle_and_rate() {
        let start = initial();
        let acceleration = -1.25;
        for steps in [1, 200] {
            let dt = 2.0 / f64::from(steps);
            let mut state = start;
            for _ in 0..steps {
                state = held_load_step(state, [0.0; 2], acceleration, 0.0, dt);
            }
            close(state.q_rad_s, start.q_rad_s + 2.0 * acceleration, 1e-13);
            close(
                state.theta_rad,
                start.theta_rad + 2.0 * start.q_rad_s + 2.0 * acceleration,
                1e-13,
            );
        }
    }

    #[test]
    fn uniform_wind_is_galilean_and_allows_backward_ground_motion() {
        let wind = [-13.0, 0.0];
        let boost = [4.5, -1.3];
        let boosted_wind = [wind[0] + boost[0], wind[1] + boost[1]];
        let air = [8.1, 0.7];
        let mut state = initial();
        state.ground_mps = ground_velocity(air, state.theta_rad, wind);
        let mut boosted = state;
        boosted.ground_mps = ground_velocity(air, state.theta_rad, boosted_wind);
        let mut body = air;
        let mut boosted_body = air;
        let dt = 0.01;
        for _ in 0..300 {
            state.ground_mps = ground_velocity(body, state.theta_rad, wind);
            boosted.ground_mps = ground_velocity(boosted_body, boosted.theta_rad, boosted_wind);
            state = held_load_step(state, [1.2, -0.8], 0.12, GRAVITY_MPS2, dt);
            boosted = held_load_step(boosted, [1.2, -0.8], 0.12, GRAVITY_MPS2, dt);
            body = body_air_velocity(state.ground_mps, state.theta_rad, wind);
            boosted_body = body_air_velocity(boosted.ground_mps, boosted.theta_rad, boosted_wind);
        }
        for axis in 0..2 {
            close(
                boosted.ground_mps[axis] - state.ground_mps[axis],
                boost[axis],
                1e-11,
            );
            close(
                boosted.position_m[axis] - state.position_m[axis],
                3.0 * boost[axis],
                1e-11,
            );
            close(boosted_body[axis], body[axis], 1e-11);
        }
        assert!(
            state.position_m[0] < initial().position_m[0],
            "headwind stronger than airspeed must permit backward flight"
        );
    }

    #[test]
    fn rotating_held_body_force_converges_to_analytic_world_motion() {
        let start = initial();
        let acceleration = [1.3, -0.9];
        let duration = 1.3;
        let q = start.q_rad_s;
        let (s0, c0) = start.theta_rad.sin_cos();
        let (s1, c1) = (start.theta_rad + q * duration).sin_cos();
        let velocity_integrals = [(s1 - s0) / q, (c0 - c1) / q];
        let position_integrals = [
            (c0 - c1) / (q * q) - duration * s0 / q,
            duration * c0 / q - (s1 - s0) / (q * q),
        ];
        let mut errors = Vec::new();
        for steps in [16, 32, 64] {
            let mut state = start;
            for _ in 0..steps {
                state = held_load_step(state, acceleration, 0.0, 0.0, duration / f64::from(steps));
            }
            let mut error = 0.0;
            for (actual, integrals, base) in [
                (state.ground_mps, velocity_integrals, start.ground_mps),
                (
                    state.position_m,
                    position_integrals,
                    [
                        start.position_m[0] + duration * start.ground_mps[0],
                        start.position_m[1] + duration * start.ground_mps[1],
                    ],
                ),
            ] {
                let exact = [
                    base[0] + acceleration[0] * integrals[0] + acceleration[1] * integrals[1],
                    base[1] + acceleration[0] * integrals[1] - acceleration[1] * integrals[0],
                ];
                error += (actual[0] - exact[0]).abs() + (actual[1] - exact[1]).abs();
            }
            errors.push(error);
        }
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.95..4.05).contains(&ratio),
                "held rotating load error must be second order: {errors:?}, ratio={ratio}"
            );
        }
    }

    fn coupled_step(
        state: FlightMotion,
        dt: f64,
        mut loads: impl FnMut(FlightMotion) -> FlightAcceleration,
    ) -> FlightMotion {
        let stage = midpoint_stage(state, loads(state), dt);
        finish_midpoint(state, stage, loads(stage), dt)
    }

    #[test]
    fn velocity_pitch_and_rotor_dependent_loads_converge_at_second_order() {
        // Independent triangular linear system: rotor torque damps omega,
        // omega drives thrust, translation and pitch have different drag rates.
        let start = FlightMotion {
            longitudinal: initial(),
            rotor_rad_s: 43.0,
        };
        let duration = 1.3_f64;
        let (drag, rotor_drag, pitch_drag, heave_drag, thrust) = (1.1, 0.6, 0.8, 0.9, 0.4);
        let e = |rate: f64| (-rate * duration).exp();
        let integral = |rate: f64| (1.0 - e(rate)) / rate;
        let exact = [
            start.longitudinal.position_m[0]
                + start.longitudinal.ground_mps[0] * integral(drag)
                + thrust * start.rotor_rad_s / (drag - rotor_drag)
                    * (integral(rotor_drag) - integral(drag)),
            start.longitudinal.position_m[1]
                + start.longitudinal.ground_mps[1] * integral(heave_drag),
            start.longitudinal.ground_mps[0] * e(drag)
                + thrust * start.rotor_rad_s / (drag - rotor_drag) * (e(rotor_drag) - e(drag)),
            start.longitudinal.ground_mps[1] * e(heave_drag),
            start.longitudinal.theta_rad + start.longitudinal.q_rad_s * integral(pitch_drag),
            start.longitudinal.q_rad_s * e(pitch_drag),
            start.rotor_rad_s * e(rotor_drag),
        ];
        let mut errors = Vec::new();
        for steps in [20, 40, 80] {
            let mut state = start;
            for _ in 0..steps {
                state = coupled_step(state, duration / f64::from(steps), |s| FlightAcceleration {
                    world_mps2: [
                        -drag * s.longitudinal.ground_mps[0] + thrust * s.rotor_rad_s,
                        -heave_drag * s.longitudinal.ground_mps[1],
                    ],
                    pitch_rad_s2: -pitch_drag * s.longitudinal.q_rad_s,
                    rotor_rad_s2: -rotor_drag * s.rotor_rad_s,
                });
            }
            let actual = [
                state.longitudinal.position_m[0],
                state.longitudinal.position_m[1],
                state.longitudinal.ground_mps[0],
                state.longitudinal.ground_mps[1],
                state.longitudinal.theta_rad,
                state.longitudinal.q_rad_s,
                state.rotor_rad_s,
            ];
            errors.push(
                actual
                    .iter()
                    .zip(exact)
                    .map(|(a, e)| (a - e).abs())
                    .sum::<f64>(),
            );
        }
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.8..4.3).contains(&ratio),
                "coupled errors {errors:?}, ratio {ratio}"
            );
        }
    }

    #[test]
    fn midpoint_geometry_resolves_position_dependent_force() {
        let start = FlightMotion {
            longitudinal: initial(),
            rotor_rad_s: 43.0,
        };
        let duration = 0.9_f64;
        let frequency = 1.7_f64;
        let angle = frequency * duration;
        let exact_x = start.longitudinal.position_m[0] * angle.cos()
            + start.longitudinal.ground_mps[0] * angle.sin() / frequency;
        let exact_v = start.longitudinal.ground_mps[0] * angle.cos()
            - frequency * start.longitudinal.position_m[0] * angle.sin();
        let mut errors = Vec::new();
        for steps in [12, 24, 48] {
            let mut state = start;
            for _ in 0..steps {
                state = coupled_step(state, duration / f64::from(steps), |s| FlightAcceleration {
                    world_mps2: [-frequency * frequency * s.longitudinal.position_m[0], 0.0],
                    pitch_rad_s2: 0.0,
                    rotor_rad_s2: 0.0,
                });
            }
            errors.push(
                (state.longitudinal.position_m[0] - exact_x).abs()
                    + (state.longitudinal.ground_mps[0] - exact_v).abs(),
            );
            assert_eq!(state.rotor_rad_s, start.rotor_rad_s);
        }
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.8..4.3).contains(&ratio),
                "spatial-force errors {errors:?}, ratio {ratio}"
            );
        }
    }

    #[test]
    fn midpoint_impulse_ballistics_and_constant_angular_accelerations_are_exact() {
        let start = FlightMotion {
            longitudinal: initial(),
            rotor_rad_s: 43.0,
        };
        let rates = FlightAcceleration {
            world_mps2: [1.7, -GRAVITY_MPS2],
            pitch_rad_s2: -0.4,
            rotor_rad_s2: -1.3,
        };
        let dt = 0.37;
        let next = coupled_step(start, dt, |_| rates);
        for axis in 0..2 {
            close(
                next.longitudinal.ground_mps[axis],
                start.longitudinal.ground_mps[axis] + dt * rates.world_mps2[axis],
                1e-14,
            );
            close(
                next.longitudinal.position_m[axis],
                start.longitudinal.position_m[axis]
                    + dt * start.longitudinal.ground_mps[axis]
                    + 0.5 * dt * dt * rates.world_mps2[axis],
                1e-14,
            );
        }
        close(
            next.longitudinal.theta_rad,
            start.longitudinal.theta_rad
                + dt * start.longitudinal.q_rad_s
                + 0.5 * dt * dt * rates.pitch_rad_s2,
            1e-14,
        );
        close(
            next.longitudinal.q_rad_s,
            start.longitudinal.q_rad_s + dt * rates.pitch_rad_s2,
            1e-14,
        );
        close(
            next.rotor_rad_s,
            start.rotor_rad_s + dt * rates.rotor_rad_s2,
            1e-14,
        );
    }

    #[test]
    fn noncollinear_coupled_load_closes_orbital_and_spin_angular_impulse() {
        let start = FlightMotion {
            longitudinal: initial(),
            rotor_rad_s: 43.0,
        };
        let mass = 257.0;
        let pitch_inertia = 369.0;
        let origin = [-2.1, 4.7];
        let dt = 0.17;
        let force = |s: FlightMotion| {
            [
                431.0 + 23.0 * s.longitudinal.position_m[1] - 19.0 * s.longitudinal.ground_mps[0],
                -287.0 - 31.0 * s.longitudinal.position_m[0] + 11.0 * s.longitudinal.ground_mps[1],
            ]
        };
        let torque = |s: FlightMotion| 71.0 - 43.0 * s.longitudinal.q_rad_s + 0.2 * s.rotor_rad_s;
        let rates = |s: FlightMotion| FlightAcceleration {
            world_mps2: force(s).map(|value| value / mass),
            pitch_rad_s2: torque(s) / pitch_inertia,
            rotor_rad_s2: -0.3 * s.rotor_rad_s,
        };
        let stage = midpoint_stage(start, rates(start), dt);
        let next = finish_midpoint(start, stage, rates(stage), dt);
        let angular = |s: FlightMotion| {
            mass * ((s.longitudinal.position_m[0] - origin[0]) * s.longitudinal.ground_mps[1]
                - (s.longitudinal.position_m[1] - origin[1]) * s.longitudinal.ground_mps[0])
                + pitch_inertia * s.longitudinal.q_rad_s
        };
        let f = force(stage);
        let expected = dt
            * ((stage.longitudinal.position_m[0] - origin[0]) * f[1]
                - (stage.longitudinal.position_m[1] - origin[1]) * f[0]
                + torque(stage));
        close(angular(next) - angular(start), expected, 5e-11);
        assert!(
            (start.longitudinal.ground_mps[0] * f[1] - start.longitudinal.ground_mps[1] * f[0])
                .abs()
                > 1000.0,
            "the fixture must expose a spurious orbital couple from the drift"
        );
    }
}
