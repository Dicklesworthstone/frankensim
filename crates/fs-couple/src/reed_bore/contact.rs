//! Implicit moving-aperture and nonadhesive lay contact in the reed step.
//!
//! Opening, velocity, pressure and contact work share one midpoint state.
//! The contact potential remains fs-dcontact-owned. Incoming characteristic
//! pressure and external body flow are held inputs, not an implicit whole-bore
//! or plate solve. Dynamic jets use a continuous, algebraically solved
//! Bernoulli junction. Second-order accuracy is limited to smooth branches.

use super::{
    AcousticRealizeError, BeatingReed, FastSolveStats, Obstacle, ReedSolverMode,
    reed_pressure_face, reed_structural, reed_swept_flow,
};
use crate::unilateral_contact::distributed::ApertureContactStep;
use crate::bernoulli_aperture::moving::characteristic_pressure as solve_moving_aperture_wave;

#[derive(Clone, Copy)]
struct Trial {
    residual: f64,
    scale: f64,
    state: (f64, f64, f64),
}

fn resolved(trial: Trial) -> bool {
    trial.residual.abs() <= 1e-13 * trial.scale.max(f64::MIN_POSITIVE)
}

fn finish(trial: Trial) -> Result<(f64, f64, f64), AcousticRealizeError> {
    // Adjacent floating-point endpoints alone cannot establish convergence.
    // The equation and its scale are impulses, not pressure-scaled flows.
    if trial.residual.abs() > 2e-11 * trial.scale.max(f64::MIN_POSITIVE) {
        return Err(AcousticRealizeError::Reed {
            what: "implicit reed contact did not resolve the momentum balance",
        });
    }
    Ok(trial.state)
}

// Ordered binary64 keys bound refinement even for subnormal velocities.
fn velocity_key(value: f64) -> u64 {
    let bits = value.to_bits();
    if bits >> 63 == 0 {
        bits | (1_u64 << 63)
    } else {
        !bits
    }
}

fn from_velocity_key(key: u64) -> f64 {
    f64::from_bits(if key >> 63 == 0 {
        !key
    } else {
        key & !(1_u64 << 63)
    })
}

#[allow(clippy::too_many_arguments)] // one coherent mechanical/junction step
pub(super) fn step(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    p_minus: f64,
    p_m: f64,
    y: f64,
    v: f64,
    dt: f64,
    u_body: f64,
    lay: &Obstacle,
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    step_with_opening(reed,rho,zc,p_minus,p_m,y,v,dt,u_body,lay,
        |opening| Ok(opening.max(0.0)))
}

// Same nonlinear equation and solver; only the actual midpoint slit geometry
// differs. A supplied callback returns nonnegative flow-equivalent height.
#[allow(clippy::too_many_arguments)]
pub(crate) fn step_with_opening(
    reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64,
    y: f64, v: f64, dt: f64, u_body: f64, lay: &Obstacle,
    flow_opening: impl Fn(f64) -> Result<f64, AcousticRealizeError>,
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    if ![
        rho,
        zc,
        p_minus,
        p_m,
        y,
        v,
        dt,
        u_body,
        reed.mass_kg,
        reed.rest_opening_m,
        reed.width_m,
        reed.closing_pressure_pa,
        reed.stiffness_n_m,
        reed.damping_ratio,
    ]
    .iter()
    .all(|value| value.is_finite())
        || rho <= 0.0
        || zc <= 0.0
        || dt <= 0.0
        || reed.mass_kg <= 0.0
        || reed.rest_opening_m <= 0.0
        || reed.width_m <= 0.0
        || reed.closing_pressure_pa <= 0.0
        || reed.stiffness_n_m < 0.0
        || reed.damping_ratio < 0.0
    {
        return Err(AcousticRealizeError::Reed {
            what: "implicit reed contact requires finite admitted mechanics and positive time step",
        });
    }
    let contact = ApertureContactStep::new(lay, y)
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
    // The held-aperture solve is a predictor only. Even free flight must solve
    // again: changing opening changes jet flow and hence the pressure load.
    let free = super::step_massive_reed(
        reed,
        rho,
        zc,
        p_minus,
        p_m,
        y,
        v,
        dt,
        u_body,
        None,
        ReedSolverMode::Strict,
        &mut FastSolveStats::default(),
    )?;
    let face = reed_pressure_face(reed);
    let (stiffness, damping) = reed_structural(reed);
    let half_dt = 0.5 * dt;
    let mass_mid = reed.mass_kg + half_dt * damping + half_dt * half_dt * stiffness;
    let rhs = reed.mass_kg * v - half_dt * stiffness * (y - reed.rest_opening_m);
    if !mass_mid.is_finite() || mass_mid <= 0.0 || !rhs.is_finite() || !face.is_finite() {
        return Err(AcousticRealizeError::Reed {
            what: "implicit reed contact mechanics left the finite set",
        });
    }
    let evaluate = |velocity: f64| -> Result<Trial, AcousticRealizeError> {
        let next_y = y + dt * velocity;
        let force = contact.response(next_y, velocity)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?.force;
        // Clamp the MIDPOINT opening, not each endpoint separately. This is
        // the same configuration used by the midpoint spring and velocity.
        let opening = flow_opening(f64::midpoint(y, next_y))?;
        let outgoing = solve_moving_aperture_wave(
            reed.width_m,
            opening,
            rho,
            zc,
            p_minus,
            p_m,
            u_body + reed_swept_flow(face, velocity),
        )?;
        let pressure_impulse = half_dt * face * (p_m - (outgoing + p_minus));
        let contact_impulse = half_dt * force;
        let inertial = mass_mid * velocity;
        let residual = inertial - rhs + pressure_impulse - contact_impulse;
        let scale = inertial.abs() + rhs.abs() + pressure_impulse.abs() + contact_impulse.abs();
        let next_v = 2.0 * velocity - v;
        if ![residual, scale, outgoing, next_y, next_v]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(AcousticRealizeError::Reed {
                what: "implicit reed contact trial left the finite set",
            });
        }
        Ok(Trial {
            residual,
            scale,
            state: (outgoing, next_y, next_v),
        })
    };
    let mut velocity = f64::midpoint(v, free.2);
    let mut current = evaluate(velocity)?;
    if resolved(current) {
        return finish(current);
    }
    let mut span = v.abs().max((rhs / mass_mid).abs()).max(1.0);
    let (mut lo, mut hi) = (-span, span);
    let (mut left, mut right) = (evaluate(lo)?, evaluate(hi)?);
    for _ in 0..64 {
        if left.residual <= 0.0 && right.residual >= 0.0 {
            break;
        }
        span *= 2.0;
        if left.residual > 0.0 {
            lo -= span;
            left = evaluate(lo)?;
        }
        if right.residual < 0.0 {
            hi += span;
            right = evaluate(hi)?;
        }
    }
    if left.residual > 0.0 || right.residual < 0.0 {
        return Err(AcousticRealizeError::Reed {
            what: "implicit reed contact could not bracket the momentum balance",
        });
    }
    // Bounded, safeguarded Newton avoids exhaustive nested bisection on each
    // smooth audio sample. The aperture feedback need not be globally monotone:
    // every accepted trial still solves the original equation, and a failed
    // acceleration falls back to the retained sign-changing bracket.
    for _ in 0..8 {
        if velocity <= lo || velocity >= hi {
            break;
        }
        if resolved(current) {
            return finish(current);
        }
        if current.residual < 0.0 {
            lo = velocity;
            left = current;
        } else {
            hi = velocity;
            right = current;
        }
        let increment = 1e-6 * (1.0 + velocity.abs());
        let lower = (velocity - increment).max(lo);
        let upper = (velocity + increment).min(hi);
        if upper <= lower {
            break;
        }
        let slope = (evaluate(upper)?.residual - evaluate(lower)?.residual) / (upper - lower);
        if !slope.is_finite() || slope <= 0.0 {
            break;
        }
        let candidate = velocity - current.residual / slope;
        if !candidate.is_finite() || candidate <= lo || candidate >= hi {
            break;
        }
        let next = evaluate(candidate)?;
        if next.residual.abs() >= current.residual.abs() {
            break;
        }
        velocity = candidate;
        current = next;
    }
    // No uniqueness or branch-independent result is inferred from bisection.
    for _ in 0..=64 {
        let (low_key, high_key) = (velocity_key(lo), velocity_key(hi));
        if resolved(left) || resolved(right) || high_key - low_key <= 1 {
            return finish(if left.residual.abs() <= right.residual.abs() {
                left
            } else {
                right
            });
        }
        let mid = from_velocity_key(low_key + (high_key - low_key) / 2);
        let trial = evaluate(mid)?;
        if resolved(trial) {
            return finish(trial);
        }
        if trial.residual < 0.0 {
            lo = mid;
            left = trial;
        } else {
            hi = mid;
            right = trial;
        }
    }
    Err(AcousticRealizeError::Reed {
        what: "implicit reed contact exceeded the binary64 refinement budget",
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unilateral_contact::slit_lay;

    fn reed() -> BeatingReed {
        BeatingReed {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
            blowing_pressure_pa: 2800.0,
            attack_s: 0.008,
            mass_kg: 1e-5,
            stiffness_n_m: 500.0,
            damping_ratio: 0.35,
        }
    }

    #[test]
    fn g1_reed_contact_step_balances_impact_release_and_bore_work() {
        let reed = reed();
        let (rho, zc, incoming, body) = (1.2, 2.7e7, 75.0, 2e-7);
        let (_, damping) = reed_structural(reed);
        let mut transitions = [false; 4];
        for alpha in [1.5, 2.0, 2.3] {
            for internal_loss in [0.0, 5.0] {
                let lay = slit_lay(1e8, alpha)
                    .unwrap()
                    .with_internal_loss(internal_loss)
                    .unwrap();
                // Independent analytical potential, not the solver's gradient.
                let potential = |opening: f64| {
                    1e8 * (-opening).max(0.0).powf(alpha + 1.0) / (alpha + 1.0)
                };
                let energy = |opening: f64, velocity: f64| {
                    0.5 * reed.mass_kg * velocity * velocity
                        + 0.5 * reed.stiffness_n_m * (opening - reed.rest_opening_m).powi(2)
                        + potential(opening)
                };
                for y in [-1e-4, 0.0, 1e-5, 4e-4] {
                    for v in [-2.0, -0.2, 0.1, 1.0] {
                        for dt in [1e-6, 1e-5, 1e-4, 5e-4] {
                            for mouth in [-500.0, 0.0, 500.0, 10000.0] {
                                let (pressure, y1, v1) = super::super::step_massive_reed(
                                    reed, rho, zc, incoming, mouth, y, v, dt, body,
                                    Some(&lay), ReedSolverMode::Strict,
                                    &mut FastSolveStats::default(),
                                ).unwrap();
                                let vm = 0.5 * (v + v1);
                                let dp = mouth - (pressure + incoming);
                                let opening = (0.5 * (y + y1)).max(0.0);
                                let jet = reed.width_m * opening * dp.signum()
                                    * (2.0 * dp.abs() / rho).sqrt();
                                let wave = (pressure - incoming) / zc;
                                let swept = -reed_pressure_face(reed) * vm;
                                let flow_scale = jet.abs() + body.abs() + swept.abs() + wave.abs();
                                let roundoff = 64.0 * f64::EPSILON
                                    * (pressure.abs() + incoming.abs()) / zc;
                                assert!((jet + body + swept - wave).abs()
                                    <= 2e-10 * flow_scale + roundoff);
                                let before = energy(y, v);
                                let after = energy(y1, v1);
                                let supplied = dt * dp * (wave - body);
                                let inertia = reed.mass_kg * (v1 - v) / dt;
                                let spring = reed.stiffness_n_m
                                    * (0.5 * (y + y1) - reed.rest_opening_m);
                                let viscous = damping * vm;
                                let pressure_force = reed_pressure_face(reed) * dp;
                                let force = inertia + spring + viscous + pressure_force;
                                let force_scale = inertia.abs() + spring.abs()
                                    + viscous.abs() + pressure_force.abs();
                                let force_tolerance = 2e-11 * force_scale.max(f64::MIN_POSITIVE);
                                let contact_loss = potential(y) - potential(y1) - dt * force * vm;
                                let loss = dt * (damping * vm * vm + dp * jet) + contact_loss;
                                let scale = before.abs() + after.abs() + supplied.abs() + loss.abs();
                                let tolerance = 2e-11 * scale.max(f64::MIN_POSITIVE);
                                assert!((after - before + loss - supplied).abs() <= tolerance);
                                assert!(contact_loss >= -tolerance, "contact created energy");
                                assert!(force >= -force_tolerance, "adhesive lay");
                                if internal_loss == 0.0 {
                                    assert!(contact_loss.abs() <= tolerance, "elastic contact lost energy");
                                }
                                if (y1 - y).abs() > 1e-10 * y.abs().max(y1.abs()) {
                                    let elastic = (potential(y) - potential(y1)) / (y1 - y);
                                    let expected = (elastic * (1.0 - internal_loss * vm)).max(0.0);
                                    assert!((force - expected).abs() <= force_tolerance);
                                }
                                transitions[usize::from(y >= 0.0) * 2 + usize::from(y1 >= 0.0)] = true;
                            }
                        }
                    }
                }
            }
        }
        assert!(transitions.iter().all(|seen| *seen));
    }

    #[test]
    fn g1_midpoint_aperture_matches_an_independently_manufactured_step() {
        let reed = reed();
        let lay = slit_lay(1e8, 2.0).unwrap();
        let (rho, zc, incoming, body) = (1.2, 1e6, 75.0, 2e-7);
        let (y, y1, v, dt) = (4e-4, 3e-4, 0.2, 1e-4);
        let vm = (y1 - y) / dt;
        let v1 = 2.0 * vm - v;
        let face = reed.stiffness_n_m * reed.rest_opening_m / reed.closing_pressure_pa;
        let damping = 2.0 * reed.damping_ratio * (reed.stiffness_n_m * reed.mass_kg).sqrt();
        let dp = -(reed.mass_kg * (v1 - v) / dt
            + reed.stiffness_n_m * (0.5 * (y + y1) - reed.rest_opening_m)
            + damping * vm) / face;
        let jet = reed.width_m * (0.5 * (y + y1)) * dp.signum()
            * (2.0 * dp.abs() / rho).sqrt();
        let expected_pressure = incoming + zc * (jet + body - face * vm);
        let mouth = dp + expected_pressure + incoming;
        let actual = step(reed, rho, zc, incoming, mouth, y, v, dt, body, &lay).unwrap();
        assert!((actual.0 - expected_pressure).abs() < 1e-8);
        assert!((actual.1 - y1).abs() < 1e-13);
        assert!((actual.2 - v1).abs() < 1e-9);
        let old = super::super::step_massive_reed(
            reed, rho, zc, incoming, mouth, y, v, dt, body, None,
            ReedSolverMode::Strict, &mut FastSolveStats::default(),
        ).unwrap();
        assert!((old.1 - y1).abs() > 1e-7, "held-opening negative control must fail");
    }

    // An independent continuous ODE oracle: invert the Bernoulli quadratic
    // algebraically, not through the characteristic solver used in production.
    fn derivative(state: [f64; 2]) -> [f64; 2] {
        let reed = reed();
        let (rho, zc, incoming, body, mouth) = (1.2, 1e6, 75.0, 2e-7, 1000.0);
        let face = reed.stiffness_n_m * reed.rest_opening_m / reed.closing_pressure_pa;
        let damping = 2.0 * reed.damping_ratio * (reed.stiffness_n_m * reed.mass_kg).sqrt();
        let [y, v] = state;
        assert!(y > 0.0, "smooth oracle must remain outside contact");
        let drive = mouth - 2.0 * incoming - zc * (body - face * v);
        let b = zc * reed.width_m * y * (2.0_f64 / rho).sqrt();
        let root = 2.0 * drive.abs() / (b.hypot(2.0 * drive.abs().sqrt()) + b);
        let dp = drive.signum() * root * root;
        [v, (-reed.stiffness_n_m * (y - reed.rest_opening_m) - damping * v - face * dp)
            / reed.mass_kg]
    }

    fn rk4_reference(mut state: [f64; 2], duration: f64) -> [f64; 2] {
        let h = duration / 4096.0;
        for _ in 0..4096 {
            let a = derivative(state);
            let b = derivative([state[0] + 0.5 * h * a[0], state[1] + 0.5 * h * a[1]]);
            let c = derivative([state[0] + 0.5 * h * b[0], state[1] + 0.5 * h * b[1]]);
            let d = derivative([state[0] + h * c[0], state[1] + h * c[1]]);
            for i in 0..2 {
                state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
        state
    }

    #[test]
    fn g1_moving_aperture_refines_at_second_order_against_independent_ode() {
        let reed = reed();
        let lay = slit_lay(1e8, 2.0).unwrap();
        let initial = [1.1 * reed.rest_opening_m, 0.25];
        let duration = 2e-4;
        let reference = rk4_reference(initial, duration);
        let mut errors = Vec::new();
        let mut frozen_errors = Vec::new();
        for count in [64, 128, 256] {
            let dt = duration / f64::from(count);
            for (obstacle, output) in [
                (Some(&lay), &mut errors),
                (None, &mut frozen_errors),
            ] {
                let [mut y, mut v] = initial;
                for _ in 0..count {
                    let (_, next_y, next_v) = super::super::step_massive_reed(
                        reed, 1.2, 1e6, 75.0, 1000.0, y, v, dt, 2e-7, obstacle,
                        ReedSolverMode::Strict, &mut FastSolveStats::default(),
                    ).unwrap();
                    (y, v) = (next_y, next_v);
                }
                let omega = (reed.stiffness_n_m / reed.mass_kg).sqrt();
                output.push(((y - reference[0]) / reed.rest_opening_m)
                    .hypot((v - reference[1]) / (reed.rest_opening_m * omega)));
            }
        }
        for pair in errors.windows(2) {
            assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2, "{errors:?}");
        }
        let frozen_ratio = frozen_errors[1] / frozen_errors[2];
        assert!(frozen_ratio > 1.8 && frozen_ratio < 2.3, "{frozen_errors:?}");
        assert!(errors[2] < 0.2 * frozen_errors[2]);
    }

    #[test]
    fn g0_stationary_zero_load_has_no_spurious_aperture_motion() {
        let reed = reed();
        let lay = slit_lay(1e8, 2.0).unwrap();
        let result = step(reed, 1.2, 1e6, 0.0, 0.0, reed.rest_opening_m,
            0.0, 1e-5, 0.0, &lay).unwrap();
        assert_eq!(result, (0.0, reed.rest_opening_m, 0.0));
    }

    #[test]
    fn g0_contact_step_refuses_invalid_time_without_publishing_state() {
        let lay = slit_lay(1e8, 2.0).unwrap();
        for dt in [0.0, -1e-6, f64::NAN, f64::INFINITY] {
            assert!(step(reed(), 1.2, 2.7e7, 0.0, 0.0, -1e-4, -0.2, dt, 0.0, &lay).is_err());
        }
    }
}
