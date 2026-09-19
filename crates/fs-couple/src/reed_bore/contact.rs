//! Implicit, nonadhesive lay contact in the characteristic reed step.
//!
//! Contact uses the existing PHS discrete gradient of fs-dcontact storage.
//! Aperture geometry is still held over a step: this closes contact work,
//! not the whole-bore energy or material-validation obligations of MR68.

use super::{
    AcousticRealizeError, BeatingReed, FastSolveStats, Obstacle, ReedSolverMode,
    reed_pressure_face, reed_structural, reed_swept_flow, solve_moving_aperture_wave,
};
use crate::unilateral_contact::SlitContactStep;

#[derive(Clone, Copy)]
struct Trial {
    residual: f64,
    scale: f64,
    state: (f64, f64, f64),
}

fn finish(trial: Trial) -> Result<(f64, f64, f64), AcousticRealizeError> {
    // A narrow bracket is not permission to accept a discontinuity or an
    // unresolved equation. This is an impulse residual, not a pressure-scaled
    // flow tolerance. The scale is the sum of the actual equation terms.
    if trial.residual.abs() > 2e-11 * trial.scale.max(f64::MIN_POSITIVE) {
        return Err(AcousticRealizeError::Reed {
            what: "implicit reed contact did not resolve the momentum balance",
        });
    }
    Ok(trial.state)
}

// Ordered binary64 keys give a bounded refinement even at subnormal roots.
// For finite lo < hi, their midpoint key is finite and inside the bracket.
// Keep -0 and +0 adjacent; neither mapping relies on unsafe code.
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
    let contact = SlitContactStep::new(lay, y)
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
    // Preserve the established no-contact trajectory bit for bit when its
    // entire linear-in-time opening segment has zero contact potential.
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
    if contact
        .coefficients(free.1)
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?
        .0
        == 0.0
    {
        return Ok(free);
    }
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
        let (elastic, contact_loss) = contact
            .coefficients(next_y)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        // For a convex obstacle potential, elastic decreases with velocity.
        // This nonadhesive force and the passive bore load therefore leave
        // a monotone scalar momentum equation, including contact crossings.
        let force = (elastic - contact_loss * velocity).max(0.0);
        let outgoing = solve_moving_aperture_wave(
            reed.width_m,
            y.max(0.0),
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
    for _ in 0..=64 {
        let (low_key, high_key) = (velocity_key(lo), velocity_key(hi));
        if left.residual == 0.0 || right.residual == 0.0 || high_key - low_key <= 1 {
            return finish(if left.residual.abs() <= right.residual.abs() {
                left
            } else {
                right
            });
        }
        let mid = from_velocity_key(low_key + (high_key - low_key) / 2);
        let trial = evaluate(mid)?;
        if trial.residual == 0.0 {
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
                // Independent analytical power-law energy, not the adapter's
                // gradient or the solver's own acceptance residual.
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
                                    reed,
                                    rho,
                                    zc,
                                    incoming,
                                    mouth,
                                    y,
                                    v,
                                    dt,
                                    body,
                                    Some(&lay),
                                    ReedSolverMode::Strict,
                                    &mut FastSolveStats::default(),
                                )
                                .unwrap();
                                let vm = 0.5 * (v + v1);
                                let dp = mouth - (pressure + incoming);
                                let jet = reed.width_m
                                    * y.max(0.0)
                                    * dp.signum()
                                    * (2.0 * dp.abs() / rho).sqrt();
                                let before = energy(y, v);
                                let after = energy(y1, v1);
                                let supplied = dt * dp * ((pressure - incoming) / zc - body);
                                // Reconstruct the actual contact force independently.
                                let inertia = reed.mass_kg * (v1 - v) / dt;
                                let spring = reed.stiffness_n_m
                                    * (0.5 * (y + y1) - reed.rest_opening_m);
                                let viscous = damping * vm;
                                let pressure_force = reed_pressure_face(reed) * dp;
                                let force = inertia + spring + viscous + pressure_force;
                                let force_scale = inertia.abs()
                                    + spring.abs()
                                    + viscous.abs()
                                    + pressure_force.abs();
                                let force_tolerance =
                                    2e-11 * force_scale.max(f64::MIN_POSITIVE);
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
                                // Independent potential secant checks the force law,
                                // rather than merely rearranging the energy balance.
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
        assert!(
            transitions.iter().all(|seen| *seen),
            "must exercise entry, exit and both persistent regimes"
        );
    }

    #[test]
    fn g0_free_flight_with_a_lay_preserves_the_existing_step_bits() {
        let reed = reed();
        let lay = slit_lay(1e8, 2.0).unwrap();
        let run = |obstacle| {
            super::super::step_massive_reed(
                reed,
                1.2,
                2.7e7,
                0.0,
                0.0,
                reed.rest_opening_m,
                0.01,
                1e-6,
                0.0,
                obstacle,
                ReedSolverMode::Strict,
                &mut FastSolveStats::default(),
            )
            .unwrap()
        };
        let free = run(None);
        let contact = run(Some(&lay));
        assert_eq!(
            (free.0.to_bits(), free.1.to_bits(), free.2.to_bits()),
            (contact.0.to_bits(), contact.1.to_bits(), contact.2.to_bits())
        );
    }

    #[test]
    fn g0_contact_step_refuses_invalid_time_without_publishing_state() {
        let lay = slit_lay(1e8, 2.0).unwrap();
        for dt in [0.0, -1e-6, f64::NAN, f64::INFINITY] {
            assert!(step(reed(), 1.2, 2.7e7, 0.0, 0.0, -1e-4, -0.2, dt, 0.0, &lay).is_err());
        }
    }
}
