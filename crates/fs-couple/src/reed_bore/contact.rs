//! Implicit moving-aperture and nonadhesive lay contact in the reed step.
//!
//! Opening, velocity, pressure and contact work share one midpoint state.
//! The contact potential remains fs-dcontact-owned. Incoming characteristic
//! pressure and external body flow are held inputs, not an implicit whole-bore
//! or plate solve. The selected continuous slit law enters every candidate
//! characteristic junction. Second-order accuracy is limited to smooth branches.

use super::{
    AcousticRealizeError, BeatingReed, Obstacle,
    reed_pressure_face, reed_structural, reed_swept_flow,
};
use crate::unilateral_contact::distributed::ApertureContactStep;
use crate::bernoulli_aperture::moving::characteristic_pressure as solve_moving_aperture_wave;
use crate::bernoulli_aperture::paired::{BladeState, ExactBlade};

/// Linear propagation of the opening coordinate over one step.
#[derive(Clone, Copy)]
pub(crate) enum OpeningPropagator<'a> {
    /// The original midpoint mechanics, unchanged.
    Midpoint,
    /// The midpoint mechanics WITHOUT the rigid closure stop (`closure_stop off` on a card): a reed that the free
    /// solve carries below zero opening keeps that opening, and the lay contact alone does the closing. Every other
    /// number is `Midpoint`'s, bit for bit; a step whose free solution stays at or above zero is `Midpoint`'s step.
    MidpointFree,
    /// Exact zero-order hold of the mass/damping/stiffness block, with the
    /// contact, pressure and this external force held over the step.
    Exact {
        /// Prepared exact maps of the opening coordinate and its Maxwell arms.
        blade: &'a ExactBlade,
        /// Held external force on the opening coordinate [N].
        external_n: f64,
        /// Accepted energy-normalized arm memories at the start of the step.
        memory: &'a [f64],
    },
}

/// What a characteristic junction hands the momentum balance: the outgoing wave and, when the blade is not driven by the
/// total mouth-to-bore drop, the drop that does drive it (the staple recovery's constriction drop, Almeida, Vergez and
/// Caussé 2007 eq. 16). A junction that returns only the outgoing wave (an `f64`) drives the blade by the total drop,
/// `p_m - (outgoing + p_minus)`, bit for bit as before.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct JunctionPressure {
    pub(crate) outgoing: f64,
    pub(crate) blade_drop_pa: Option<f64>,
}

impl From<f64> for JunctionPressure {
    fn from(outgoing: f64) -> Self {
        Self { outgoing, blade_drop_pa: None }
    }
}

impl OpeningPropagator<'_> {
    // A selected constitutive junction shares the existing mechanics/contact
    // solve. Keep the original free function for inviscid callers.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_with_junction<J: Into<JunctionPressure>>(
        self, reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64,
        y: f64, v: f64, dt: f64, u_body: f64, lay: &Obstacle,
        flow_opening: impl Fn(f64) -> Result<f64, AcousticRealizeError>,
        material_force: Option<&dyn Fn(f64) -> Result<f64, AcousticRealizeError>>,
        junction_pressure: impl Fn(f64, f64) -> Result<J, AcousticRealizeError>,
    ) -> Result<(f64, f64, f64), AcousticRealizeError> {
        step_with_opening_and_junction(reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
            lay, flow_opening, material_force, self, junction_pressure)
    }
}

#[derive(Clone, Copy)]
pub(super) struct OpeningStep {
    pub(super) outgoing: f64,
    pub(super) opening: f64,
    pub(super) velocity: f64,
    pub(super) mean_velocity: f64,
}

#[derive(Clone, Copy)]
struct Trial {
    residual: f64,
    scale: f64,
    /// Rounding floor of the evaluated residual: the pressure impulse is a
    /// difference of junction pressures, so its error follows their
    /// magnitudes, not the (possibly cancelled) difference.
    floor: f64,
    state: OpeningStep,
}

fn resolved(trial: Trial) -> bool {
    trial.residual.abs() <= 1e-13 * trial.scale.max(f64::MIN_POSITIVE)
}

fn finish(trial: Trial) -> Result<OpeningStep, AcousticRealizeError> {
    // Adjacent floating-point endpoints alone cannot establish convergence.
    // The equation and its scale are impulses, not pressure-scaled flows.
    // Near rest the pressures nearly cancel; a residual inside the computed
    // rounding floor is as resolved as the arithmetic allows.
    if trial.residual.abs() > (2e-11 * trial.scale).max(trial.floor).max(f64::MIN_POSITIVE) {
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
        |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint)
}

// Same nonlinear equation and solver; only the actual midpoint slit geometry
// differs. A supplied callback returns nonnegative flow-equivalent height.
#[allow(clippy::too_many_arguments)]
pub(crate) fn step_with_opening(
    reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64,
    y: f64, v: f64, dt: f64, u_body: f64, lay: &Obstacle,
    flow_opening: impl Fn(f64) -> Result<f64, AcousticRealizeError>,
    material_force: Option<&dyn Fn(f64) -> Result<f64, AcousticRealizeError>>,
    propagator: OpeningPropagator<'_>,
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    step_with_opening_and_junction(reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
        lay, flow_opening, material_force, propagator,
        |opening, moving| solve_moving_aperture_wave(reed.width_m, opening, rho,
            zc, p_minus, p_m, moving))
}

// The constitutive junction is evaluated at each trial gap and swept flow,
// not applied as an after-the-fact pressure or energy correction.
#[allow(clippy::too_many_arguments)]
fn step_with_opening_and_junction<J: Into<JunctionPressure>>(
    reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64,
    y: f64, v: f64, dt: f64, u_body: f64, lay: &Obstacle,
    flow_opening: impl Fn(f64) -> Result<f64, AcousticRealizeError>,
    material_force: Option<&dyn Fn(f64) -> Result<f64, AcousticRealizeError>>,
    propagator: OpeningPropagator<'_>,
    junction_pressure: impl Fn(f64, f64) -> Result<J, AcousticRealizeError>,
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    let contact = ApertureContactStep::new(lay, y)
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
    step_with_contact_response(reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
        flow_opening, material_force, propagator, junction_pressure,
        |after, velocity| contact.response(after, velocity)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string())))
        .map(|step| (step.outgoing, step.opening, step.velocity))
}

/// Reuse the same mechanics and scalar root with an already admitted contact
/// response. The full-body held-flow adapter uses the allocation-free unit-slit
/// owner; existing characteristic callers retain their original contact path.
/// The reed cannot pass through its mouthpiece: see [`step_with_contact_response_and_stop`].
#[allow(clippy::too_many_arguments)]
pub(super) fn step_with_contact_response<J: Into<JunctionPressure>>(
    reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64,
    y: f64, v: f64, dt: f64, u_body: f64,
    flow_opening: impl Fn(f64) -> Result<f64, AcousticRealizeError>,
    material_force: Option<&dyn Fn(f64) -> Result<f64, AcousticRealizeError>>,
    propagator: OpeningPropagator<'_>,
    junction_pressure: impl Fn(f64, f64) -> Result<J, AcousticRealizeError>,
    contact_response: impl Fn(f64, f64) -> Result<fs_dcontact::AffineContactResponse, AcousticRealizeError>,
) -> Result<OpeningStep, AcousticRealizeError> {
    step_with_contact_response_and_stop(reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
        flow_opening, material_force, propagator, junction_pressure, contact_response, true)
}

/// The step of [`step_with_contact_response`] with the rigid closure stop selectable.
///
/// The stop is the unilateral constraint `0 <= opening` with a non-negative multiplier. The implicit
/// midpoint solve below is unchanged and is run first. If the reed starts the step at or above zero
/// opening and that solve ends it below zero, the step is replaced by the constrained one: the reed
/// ends at exactly 0 with mean velocity `-y/dt` and end velocity 0 (rigid and inelastic: no printed
/// restitution exists for a reed on its lay), the junction and the lay contact being evaluated at that
/// landing velocity by the same code. The momentum residual rises with the mean velocity, so at the
/// landing velocity, above the penetrating root, it is not negative; that residual is the stop's own
/// push, and a negative one would mean the stop pulls the reed shut, which refuses. A reed that
/// enters the step already below zero is left to the free solve. With `closure_stop` false, or when
/// the free solve ends at or above zero, the result is the original one bit for bit.
#[allow(clippy::too_many_arguments)]
pub(super) fn step_with_contact_response_and_stop<J: Into<JunctionPressure>>(
    reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64,
    y: f64, v: f64, dt: f64, u_body: f64,
    flow_opening: impl Fn(f64) -> Result<f64, AcousticRealizeError>,
    material_force: Option<&dyn Fn(f64) -> Result<f64, AcousticRealizeError>>,
    propagator: OpeningPropagator<'_>,
    junction_pressure: impl Fn(f64, f64) -> Result<J, AcousticRealizeError>,
    contact_response: impl Fn(f64, f64) -> Result<fs_dcontact::AffineContactResponse, AcousticRealizeError>,
    closure_stop: bool,
) -> Result<OpeningStep, AcousticRealizeError> {
    if material_force.is_some() && matches!(propagator, OpeningPropagator::Exact { .. }) {
        // The exact block already carries its Maxwell arms; a midpoint material
        // force on top would count the same memory twice.
        return Err(AcousticRealizeError::Reed {
            what: "the exact opening block takes material memory as arms, not a midpoint force",
        });
    }
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
    let face = reed_pressure_face(reed);
    let (stiffness, damping) = reed_structural(reed);
    let half_dt = 0.5 * dt;
    let start = BladeState { displacement_m: y - reed.rest_opening_m, velocity_m_s: v };
    // Both propagators share the impulse form mass_mid*v_mean - rhs = h*(Fc - A dp):
    // midpoint's coefficients, or the exact map's mean velocity free + mobility*F
    // rewritten with mass_mid = h/mobility and rhs = mass_mid*free.
    let (mass_mid, rhs, predicted) = match propagator {
        OpeningPropagator::Midpoint | OpeningPropagator::MidpointFree => {
            // The held-aperture solve is a predictor only. Even free flight must
            // solve again: changing opening changes jet flow and the pressure load.
            let free = super::massive_reed_predictor(
                reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
            )?;
            (reed.mass_kg + half_dt * damping + half_dt * half_dt * stiffness,
                reed.mass_kg * v - half_dt * stiffness * (y - reed.rest_opening_m),
                f64::midpoint(v, free.2))
        }
        OpeningPropagator::Exact { blade, external_n, memory } => {
            let (free, mobility) = blade.free_and_mobility(start, external_n, dt, memory);
            let mass_mid = half_dt / mobility;
            (mass_mid, mass_mid * free, free)
        }
    };
    if !mass_mid.is_finite() || mass_mid <= 0.0 || !rhs.is_finite() || !face.is_finite() {
        return Err(AcousticRealizeError::Reed {
            what: "implicit reed contact mechanics left the finite set",
        });
    }
    let evaluate = |velocity: f64| -> Result<Trial, AcousticRealizeError> {
        let next_y = y + dt * velocity;
        let force = contact_response(next_y, velocity)?.force;
        // Clamp the MIDPOINT opening, not each endpoint separately. This is
        // the same configuration used by the midpoint spring and velocity.
        let opening = flow_opening(f64::midpoint(y, next_y))?;
        let junction: JunctionPressure =
            junction_pressure(opening, u_body + reed_swept_flow(face, velocity))?.into();
        let outgoing = junction.outgoing;
        // The blade is driven by the total mouth-to-bore drop unless the junction names the drop that drives it.
        let blade_drop = junction.blade_drop_pa.unwrap_or(p_m - (outgoing + p_minus));
        let pressure_impulse = half_dt * face * blade_drop;
        let contact_impulse = half_dt * force;
        let inertial = mass_mid * velocity;
        let mut residual = inertial - rhs + pressure_impulse - contact_impulse;
        let mut scale = inertial.abs() + rhs.abs() + pressure_impulse.abs() + contact_impulse.abs();
        // The closed-form junction and the pressure difference each cost a few
        // roundings of their operands (8 eps of their magnitudes, bounded).
        let pressure_magnitude = half_dt * face.abs() * (p_m.abs() + p_minus.abs() + outgoing.abs()
            + (zc * (u_body + reed_swept_flow(face, velocity))).abs());
        let mut floor = 8.0 * f64::EPSILON * (pressure_magnitude + scale);
        if let Some(restoring) = material_force {
            // Condensed material force uses this very candidate. No lagged
            // post-step stress correction and no change to the original path.
            let impulse = half_dt * restoring(next_y)?;
            residual += impulse;
            scale += impulse.abs();
            floor += 8.0 * f64::EPSILON * impulse.abs();
        }
        let next_v = match propagator {
            OpeningPropagator::Midpoint | OpeningPropagator::MidpointFree => 2.0 * velocity - v,
            OpeningPropagator::Exact { blade, external_n, memory } => {
                blade.end_velocity(start, external_n + force - face * blade_drop, memory)
            }
        };
        if ![residual, scale, floor, outgoing, blade_drop, next_y, next_v]
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
            floor,
            state: OpeningStep { outgoing, opening: next_y, velocity: next_v,
                mean_velocity: velocity },
        })
    };
    // The original solve, unchanged, as one closure so the closure stop can follow it.
    let solve = || -> Result<OpeningStep, AcousticRealizeError> {
    let mut velocity = predicted;
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
    };
    let free = solve()?;
    if !closure_stop || matches!(propagator, OpeningPropagator::MidpointFree) || y < 0.0 || free.opening >= 0.0 {
        return Ok(free);
    }
    let landing = -y / dt;
    let trial = evaluate(landing)?;
    if trial.residual < -trial.floor.max(2e-11 * trial.scale) {
        return Err(AcousticRealizeError::Reed {
            what: "closure stop would pull the reed shut (non-monotone momentum balance)",
        });
    }
    Ok(OpeningStep {
        outgoing: trial.state.outgoing,
        opening: 0.0,
        velocity: 0.0,
        mean_velocity: landing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reed_bore::{FastSolveStats, ReedSolverMode};
    use crate::unilateral_contact::slit_lay;
    use crate::bernoulli_aperture::viscous::{StapleRecovery, ViscousSlit};

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
                                let (pressure, y1, v1) = step_through_the_lay(
                                    reed, rho, zc, incoming, mouth, y, v, dt, body, &lay);
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
        let actual_free = super::super::step_massive_reed(
            reed, rho, zc, incoming, mouth, y, v, dt, body, None,
            ReedSolverMode::Strict, &mut FastSolveStats::default(),
        ).unwrap();
        assert!((actual_free.0 - expected_pressure).abs() < 1e-8);
        assert!((actual_free.1 - y1).abs() < 1e-13);
        assert!((actual_free.2 - v1).abs() < 1e-9);
        let old = super::super::massive_reed_predictor(
            reed, rho, zc, incoming, mouth, y, v, dt, body,
        ).unwrap();
        assert!((old.1 - y1).abs() > 1e-7, "frozen-opening predictor is not a physical step");
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
        let mut contact_errors = Vec::new();
        let mut free_errors = Vec::new();
        let mut predictor_errors = Vec::new();
        for count in [64, 128, 256] {
            let dt = duration / f64::from(count);
            for (obstacle, output) in [
                (Some(&lay), &mut contact_errors),
                (None, &mut free_errors),
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
            let [mut y, mut v] = initial;
            for _ in 0..count {
                let (_, next_y, next_v) = super::super::massive_reed_predictor(
                    reed, 1.2, 1e6, 75.0, 1000.0, y, v, dt, 2e-7,
                ).unwrap();
                (y, v) = (next_y, next_v);
            }
            let omega = (reed.stiffness_n_m / reed.mass_kg).sqrt();
            predictor_errors.push(((y - reference[0]) / reed.rest_opening_m)
                .hypot((v - reference[1]) / (reed.rest_opening_m * omega)));
        }
        for errors in [&contact_errors, &free_errors] {
            for pair in errors.windows(2) {
                assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2, "{errors:?}");
            }
        }
        for (a, b) in contact_errors.iter().zip(&free_errors) {
            assert!((a - b).abs() <= 1e-12 * a.max(*b).max(f64::MIN_POSITIVE));
        }
        let predictor_ratio = predictor_errors[1] / predictor_errors[2];
        assert!(predictor_ratio > 1.8 && predictor_ratio < 2.3, "{predictor_errors:?}");
        assert!(free_errors[2] < 0.2 * predictor_errors[2]);
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

    /// The lay step WITHOUT the closure stop: the free implicit midpoint scheme with the lay's own penetration branch, whose
    /// impact, release and work balance `g1_reed_contact_step_balances_impact_release_and_bore_work` asserts. With the stop
    /// on, a step that would end below zero opening ends shut with end velocity 0 and no longer satisfies those identities
    /// by construction; the closure-stop tests below cover it.
    #[allow(clippy::too_many_arguments)]
    fn step_through_the_lay(
        reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64, y: f64, v: f64, dt: f64, u_body: f64,
        lay: &Obstacle,
    ) -> (f64, f64, f64) {
        let contact = ApertureContactStep::new(lay, y).unwrap();
        let next = step_with_contact_response_and_stop(
            reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
            |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
            |opening, moving| solve_moving_aperture_wave(
                reed.width_m, opening, rho, zc, p_minus, p_m, moving),
            |after, velocity| contact.response(after, velocity)
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string())),
            false,
        ).unwrap();
        (next.outgoing, next.opening, next.velocity)
    }

    // ---- closure stop (MB 3381): a reed cannot pass through its mouthpiece -----------------------

    /// March the reed from rest under a held mouth pressure against a lay, with the closure stop on or off.
    fn march(closure_stop: bool, mouth: f64, steps: usize, lay: &Obstacle) -> Vec<OpeningStep> {
        let reed = reed();
        let (rho, zc, incoming, body, dt) = (1.2, 2.7e7, 0.0, 0.0, 1e-5);
        let (mut y, mut v) = (reed.rest_opening_m, 0.0);
        let mut out = Vec::new();
        for _ in 0..steps {
            let contact = ApertureContactStep::new(lay, y).unwrap();
            let next = step_with_contact_response_and_stop(
                reed, rho, zc, incoming, mouth, y, v, dt, body,
                |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
                |opening, moving| solve_moving_aperture_wave(
                    reed.width_m, opening, rho, zc, incoming, mouth, moving),
                |after, velocity| contact.response(after, velocity)
                    .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string())),
                closure_stop,
            ).unwrap();
            (y, v) = (next.opening, next.velocity);
            out.push(next);
        }
        out
    }

    fn same_bits(a: &OpeningStep, b: &OpeningStep) -> bool {
        [(a.outgoing, b.outgoing), (a.opening, b.opening), (a.velocity, b.velocity),
            (a.mean_velocity, b.mean_velocity)]
            .iter().all(|(x, y)| x.to_bits() == y.to_bits())
    }

    #[test]
    fn g1_a_reed_pressed_through_a_soft_lay_stops_at_zero_opening_and_stays() {
        // A lay far too soft to hold the closing pressure: without the stop the reed goes through.
        let lay = slit_lay(1e3, 2.0).unwrap();
        let free = march(false, 10000.0, 400, &lay);
        let stopped = march(true, 10000.0, 400, &lay);
        let first = free.iter().position(|s| s.opening < 0.0)
            .expect("this setup must press the free reed through zero opening");
        // Until the free reed would have crossed zero, the stop changes nothing, bit for bit.
        for i in 0..first {
            assert!(same_bits(&free[i], &stopped[i]), "step {i} moved before any crossing");
        }
        assert!(stopped.iter().all(|s| s.opening >= 0.0 && s.opening.is_finite()));
        // It lands at exactly zero in the very step the free reed crossed, with zero velocity,
        // and stays shut while the closing pressure holds.
        for s in &stopped[first..] {
            assert_eq!(s.opening, 0.0);
            assert_eq!(s.velocity, 0.0);
        }
    }

    #[test]
    fn g1_the_stop_never_acts_on_a_reed_that_stays_open() {
        let lay = slit_lay(1e8, 2.0).unwrap();
        let free = march(false, 500.0, 400, &lay);
        let stopped = march(true, 500.0, 400, &lay);
        assert!(free.iter().all(|s| s.opening > 0.0), "this setup must keep the reed open");
        assert!(free.iter().zip(&stopped).all(|(a, b)| same_bits(a, b)));
    }

    #[test]
    fn g1_flow_through_a_closed_channel_is_zero() {
        let reed = reed();
        let (rho, zc, incoming, mouth, dt) = (1.2, 2.7e7, 0.0, 10000.0, 1e-5);
        let lay = slit_lay(1e3, 2.0).unwrap();
        // Shut, at rest, with the mouth pressure holding it closed.
        let contact = ApertureContactStep::new(&lay, 0.0).unwrap();
        let step = step_with_contact_response_and_stop(
            reed, rho, zc, incoming, mouth, 0.0, 0.0, dt, 0.0,
            |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
            |opening, moving| solve_moving_aperture_wave(
                reed.width_m, opening, rho, zc, incoming, mouth, moving),
            |after, velocity| contact.response(after, velocity)
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string())),
            true,
        ).unwrap();
        assert_eq!((step.opening, step.velocity, step.mean_velocity), (0.0, 0.0, 0.0));
        let closed_wave = solve_moving_aperture_wave(reed.width_m, 0.0, rho, zc, incoming, mouth, 0.0)
            .unwrap();
        assert_eq!(step.outgoing, closed_wave);
        let drop = mouth - (step.outgoing + incoming);
        assert_eq!(crate::bernoulli_aperture::moving::volume_flow(reed.width_m, 0.0, drop, rho), 0.0);
    }

    #[test]
    fn g0_a_closure_stop_that_would_pull_the_reed_shut_refuses() {
        let reed = reed();
        let (rho, zc, incoming, mouth, dt) = (1.2, 2.7e7, 0.0, 10000.0, 1e-5);
        // A pathological contact reaction, non-zero only in a thin band around zero opening: the free
        // solve ends below zero (outside the band) but at the landing the reaction would pull the reed
        // shut, so the momentum residual there is negative. That must refuse, not publish a state.
        let spike = |after: f64, _velocity: f64| {
            let force = if after.abs() < 1e-9 { 1e6 } else { 0.0 };
            Ok(fs_dcontact::AffineContactResponse { elastic_force: force, force, dissipated_power: 0.0 })
        };
        let run = |closure_stop| step_with_contact_response_and_stop(
            reed, rho, zc, incoming, mouth, 0.0, 0.0, dt, 0.0,
            |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
            |opening, moving| solve_moving_aperture_wave(
                reed.width_m, opening, rho, zc, incoming, mouth, moving),
            spike, closure_stop,
        );
        assert!(run(false).unwrap().opening < 0.0, "the free solve must end below zero");
        match run(true) {
            Err(AcousticRealizeError::Reed { what }) => assert!(what.contains("pull the reed shut"), "{what}"),
            other => panic!("expected the closure stop to refuse, got {:?}", other.map(|s| s.opening)),
        }
    }

    // ---- staple pressure recovery (Almeida, Vergez and Caussé 2007, eqs. 14-17): the blade is driven by (dp)_c --------

    /// The oboe description's reed: rest gap 357 um on 5.6 mm (2.0 mm^2 at rest, inside S_c / sqrt(C_P) = 4.47 mm^2).
    fn oboe_reed() -> BeatingReed {
        BeatingReed {
            rest_opening_m: 3.57e-4,
            width_m: 0.0056,
            closing_pressure_pa: 17720.0,
            blowing_pressure_pa: 7400.0,
            attack_s: 0.008,
            mass_kg: 2.2e-5,
            stiffness_n_m: 3464.0,
            damping_ratio: 0.114,
        }
    }

    /// One free midpoint step against a lay through the viscous slit junction, with the staple recovery on or off. Off is
    /// the junction as the voice has always called it: `characteristic_pressure` handing back a bare `f64`.
    #[allow(clippy::too_many_arguments)]
    fn step_through_the_slit(
        reed: BeatingReed, rho: f64, zc: f64, p_minus: f64, p_m: f64, y: f64, v: f64, dt: f64, u_body: f64,
        lay: &Obstacle, slit: ViscousSlit, recovery: Option<StapleRecovery>,
    ) -> OpeningStep {
        let contact = ApertureContactStep::new(lay, y).unwrap();
        let response = |after: f64, velocity: f64| contact.response(after, velocity)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()));
        let result = match recovery {
            Some(recovery) => step_with_contact_response_and_stop(
                reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
                |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
                |opening, moving| {
                    let junction = slit.characteristic_junction_with_recovery(
                        recovery, reed.width_m, opening, rho, zc, p_minus, p_m, moving)?;
                    Ok(JunctionPressure {
                        outgoing: junction.outgoing_pa,
                        blade_drop_pa: Some(junction.constriction_drop_pa),
                    })
                },
                response, false,
            ),
            None => step_with_contact_response_and_stop(
                reed, rho, zc, p_minus, p_m, y, v, dt, u_body,
                |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
                |opening, moving| slit.characteristic_pressure(
                    reed.width_m, opening, rho, zc, p_minus, p_m, moving),
                response, false,
            ),
        };
        result.unwrap()
    }

    #[test]
    fn g0_a_junction_that_names_no_blade_drop_is_the_original_step_bit_for_bit() {
        let reed = oboe_reed();
        let slit = ViscousSlit::new(1.81e-5, 1e-3).unwrap();
        let lay = slit_lay(1e8, 2.0).unwrap();
        for y in [1e-5, 2e-4, 3.57e-4] {
            for mouth in [-500.0, 500.0, 7400.0] {
                for dt in [1e-6, 2e-5] {
                    // A bare f64 and the same value wrapped without a blade drop take the same arithmetic.
                    let bare = step_through_the_slit(reed, 1.2, 5.9e7, 75.0, mouth, y, 0.3, dt, 1e-7, &lay, slit, None);
                    let contact = ApertureContactStep::new(&lay, y).unwrap();
                    let wrapped = step_with_contact_response_and_stop(
                        reed, 1.2, 5.9e7, 75.0, mouth, y, 0.3, dt, 1e-7,
                        |opening| Ok(opening.max(0.0)), None, OpeningPropagator::Midpoint,
                        |opening, moving| slit.characteristic_pressure(
                            reed.width_m, opening, 1.2, 5.9e7, 75.0, mouth, moving)
                            .map(|outgoing| JunctionPressure { outgoing, blade_drop_pa: None }),
                        |after, velocity| contact.response(after, velocity)
                            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string())),
                        false,
                    ).unwrap();
                    assert_eq!(
                        (bare.outgoing.to_bits(), bare.opening.to_bits(), bare.velocity.to_bits()),
                        (wrapped.outgoing.to_bits(), wrapped.opening.to_bits(), wrapped.velocity.to_bits()),
                    );
                }
            }
        }
    }

    #[test]
    fn g1_staple_blade_force_follows_the_constriction_drop_and_the_step_balances_energy() {
        let reed = oboe_reed();
        let recovery = StapleRecovery::new(0.8, 4e-6).unwrap();
        let slit = ViscousSlit::new(1.81e-5, 1e-3).unwrap();
        let (rho, zc, incoming, body) = (1.2, 5.9e7, 75.0, 1e-7);
        let face = reed_pressure_face(reed);
        let (_, damping) = reed_structural(reed);
        let (mut recovering, mut driven_by_total) = (0, 0);
        for alpha in [2.0, 2.3] {
            for internal_loss in [0.0, 5.0] {
                let lay = slit_lay(1e8, alpha).unwrap().with_internal_loss(internal_loss).unwrap();
                // Independent analytical potential, not the solver's gradient.
                let potential = |opening: f64| 1e8 * (-opening).max(0.0).powf(alpha + 1.0) / (alpha + 1.0);
                let energy = |opening: f64, velocity: f64| {
                    0.5 * reed.mass_kg * velocity * velocity
                        + 0.5 * reed.stiffness_n_m * (opening - reed.rest_opening_m).powi(2)
                        + potential(opening)
                };
                for y in [-5e-5, 1e-5, 2e-4, 3.57e-4, 4.5e-4] {
                    for v in [-2.0, -0.2, 0.1, 1.0] {
                        for dt in [1e-6, 1e-5, 5e-5] {
                            for mouth in [-500.0, 500.0, 7400.0, 15000.0] {
                                let step = step_through_the_slit(
                                    reed, rho, zc, incoming, mouth, y, v, dt, body, &lay, slit, Some(recovery));
                                let (pressure, y1, v1) = (step.outgoing, step.opening, step.velocity);
                                let vm = 0.5 * (v + v1);
                                let drop_total = mouth - (pressure + incoming);
                                let swept = -face * vm;
                                let wave = (pressure - incoming) / zc;
                                let jet = wave - body - swept; // the junction's flow balance
                                let recovered = recovery.recovered_pressure_pa(rho, jet);
                                let drop_c = drop_total + recovered; // (dp)_c = p_m - p_c, eq. 17 rearranged
                                // The blade is in balance under (dp)_c: inertia + spring + damper + face (dp)_c = lay force.
                                let inertia = reed.mass_kg * (v1 - v) / dt;
                                let spring = reed.stiffness_n_m * (0.5 * (y + y1) - reed.rest_opening_m);
                                let viscous = damping * vm;
                                let pressure_force = face * drop_c;
                                let force = inertia + spring + viscous + pressure_force;
                                // The solver stops at a momentum residual of 2e-11 of its own scale, or at the rounding floor of
                                // the pressures it differences (8 eps of their magnitudes, with the bore load z*moving among
                                // them), whichever is larger. In force units that is the most the balance can be off by here.
                                let half = 0.5 * dt;
                                let solver_scale = (reed.mass_kg * vm / half + damping * vm + half * reed.stiffness_n_m * vm).abs()
                                    + (reed.mass_kg * v / half - reed.stiffness_n_m * (y - reed.rest_opening_m)).abs()
                                    + pressure_force.abs() + force.abs();
                                let pressure_magnitude = face.abs()
                                    * (mouth.abs() + incoming.abs() + pressure.abs() + (zc * (body + swept)).abs());
                                let force_tolerance = (2e-11 * solver_scale)
                                    .max(8.0 * f64::EPSILON * (pressure_magnitude + solver_scale))
                                    .max(f64::MIN_POSITIVE);
                                let contact_loss = potential(y) - potential(y1) - dt * force * vm;
                                // Energy with the recovery in it. The slit and staple element dissipate drop_total * jet, which
                                // is not negative; the blade feels p_c while the swept flow enters the bore at p_r, and that
                                // exchange, -recovered * swept, is the only term of either sign.
                                let element = drop_total * jet;
                                let exchange = -recovered * swept;
                                let before = energy(y, v);
                                let after = energy(y1, v1);
                                let supplied = dt * drop_total * (wave - body);
                                let loss = dt * (damping * vm * vm + element + exchange) + contact_loss;
                                let scale = before.abs() + after.abs() + supplied.abs() + loss.abs();
                                // The contact loss is a difference of potentials and a force-times-path term: the force the solver
                                // may leave unresolved (above) is an energy of at most dt * force_tolerance * |vm| there.
                                let tolerance = 2e-11 * scale.max(f64::MIN_POSITIVE) + dt * force_tolerance * vm.abs();
                                assert!((after - before + loss - supplied).abs() <= tolerance,
                                    "y={y} v={v} dt={dt} mouth={mouth}");
                                let flow_scale = wave.abs() + body.abs() + swept.abs();
                                assert!(element >= -1e-9 * drop_total.abs() * flow_scale,
                                    "the slit and staple element created energy: y={y} v={v} dt={dt} mouth={mouth} element={element:e}");
                                assert!(contact_loss >= -tolerance,
                                    "the contact created energy: y={y} v={v} dt={dt} mouth={mouth} {contact_loss:e} vs {tolerance:e}");
                                assert!(force >= -force_tolerance,
                                    "adhesive lay: y={y} v={v} dt={dt} mouth={mouth} force={force:e} tolerance={force_tolerance:e}");
                                if (y1 - y).abs() > 1e-10 * y.abs().max(y1.abs()) {
                                    let elastic = (potential(y) - potential(y1)) / (y1 - y);
                                    let expected = (elastic * (1.0 - internal_loss * vm)).max(0.0);
                                    assert!((force - expected).abs() <= force_tolerance,
                                        "the blade is not driven by (dp)_c: y={y} v={v} dt={dt} mouth={mouth} {force:e} vs {expected:e} tolerance={force_tolerance:e}");
                                    // Where the recovery matters, the total drop would NOT balance the same blade.
                                    if face * recovered > 100.0 * force_tolerance {
                                        recovering += 1;
                                        let by_total = inertia + spring + viscous + face * drop_total;
                                        if (by_total - expected).abs() > force_tolerance {
                                            driven_by_total += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        assert!(recovering > 50, "the grid must exercise the recovery, got {recovering}");
        assert_eq!(driven_by_total, recovering, "the test must be able to tell (dp)_c from the total drop");
    }
}
