//! One-DOF canard control mechanism (bead wf-root-guzez.5.14.1,
//! E4.6b-i). Plan §5.1.3 + canard-mechanics-v1: pilot lever force →
//! lever/cable gain → hinge dynamics
//!
//!   I_c·δ̈ = M_aero + M_pilot − M_friction − M_stop
//!
//! with regularized Coulomb + viscous friction and penetration-spring
//! travel stops at the photo-inferred ±30° (absence-by-verification
//! status carried in the dossier). The aero hinge moment is an INPUT —
//! the caller evaluates fs-wing::hinge::hinge_load on the coupled
//! solution each tick (E4.2b); this module owns only the mechanism.
//!
//! Dossier doctrine honored: the lever gearing is NOT PUBLISHED — the
//! gain here is a tagged MODELING CHOICE input; the hinge-axis chordwise
//! position must sit inside the declared WIDE PRIOR [25, 50]% (point
//! guesses forbidden; admission-gated). Quantitative hinge-moment levels
//! remain Estimated (A7a is the only promotion path).
//!
//! Integration: a bounded, safeguarded scalar solve for midpoint velocity.
//! The stop spring uses its discrete energy gradient when a step crosses a
//! stop; within one contact branch this is ordinary implicit midpoint. Stop
//! damping acts over the fraction of the step spent beyond either stop.
//! Friction, stop loss, and drive work use the same accepted velocity.

use crate::Refusal;
use fs_math::det;

/// Hinge-axis prior band [% chord from LE] (canard-mechanics-v1
/// unknown_priors — admission refuses outside it).
pub const HINGE_AXIS_PRIOR_PCT: (f64, f64) = (25.0, 50.0);

/// Maximum safeguarded midpoint iterations (deterministic).
pub const MIDPOINT_ITERATIONS: u32 = 64;

/// The mechanism parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanardMechanism {
    /// Canard + linkage inertia about the hinge [kg·m²] (Estimated).
    pub inertia_kg_m2: f64,
    /// Pilot force → hinge torque gain [N·m/N] (MODELING CHOICE — the
    /// lever throw/ratio is not published; tag travels in provenance).
    pub lever_gain_nm_per_n: f64,
    /// Coulomb friction level [N·m].
    pub coulomb_nm: f64,
    /// Viscous friction [N·m·s/rad].
    pub viscous_nm_s: f64,
    /// Friction regularization velocity [rad/s] (tanh knee).
    pub friction_reg_rad_s: f64,
    /// Travel stop [rad] (±30° photo-inferred).
    pub stop_rad: f64,
    /// Stop penetration stiffness [N·m/rad].
    pub stop_stiffness_nm_per_rad: f64,
    /// Stop penetration damping [N·m·s/rad].
    pub stop_damping_nm_s: f64,
    /// Hinge-axis chordwise position [% chord from LE] — must sit in
    /// the declared prior band.
    pub hinge_axis_pct_chord: f64,
}

/// The registered v1 mechanism (Estimated values inside declared
/// priors; the diary's 'balanced too near the center' motivates the
/// mid-band axis).
pub const CANARD_MECH_V1: CanardMechanism = CanardMechanism {
    inertia_kg_m2: 6.0,
    lever_gain_nm_per_n: 0.35,
    coulomb_nm: 2.0,
    viscous_nm_s: 0.8,
    friction_reg_rad_s: 0.02,
    stop_rad: 30.0 * core::f64::consts::PI / 180.0,
    stop_stiffness_nm_per_rad: 4000.0,
    stop_damping_nm_s: 80.0,
    hinge_axis_pct_chord: 40.0,
};

/// Mechanism state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechState {
    /// Surface deflection [rad] (+ = nose-up command sense).
    pub delta_rad: f64,
    /// Deflection rate [rad/s].
    pub rate_rad_s: f64,
}

/// One step's receipt (per-step oracles feed on these, never totals).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReceipt {
    /// Friction dissipation this step [J] (must be ≥ 0).
    pub friction_dissipation_j: f64,
    /// Stop damping dissipation this step [J] (must be ≥ 0).
    pub stop_dissipation_j: f64,
    /// Work supplied by the held aero and pilot torques [J].
    pub drive_work_j: f64,
    /// Discrete stop torque [N·m], averaged across any stop crossing.
    pub stop_torque_nm: f64,
    /// Net torque at the accepted midpoint velocity [N·m].
    pub net_torque_nm: f64,
}

impl CanardMechanism {
    /// Admit the mechanism.
    ///
    /// # Errors
    /// `canard-mech-invalid` (non-finite / non-positive constants);
    /// `hinge-axis-outside-prior` (the dossier band, tested at both
    /// edges AND one ulp past each).
    pub fn admit(&self) -> Result<(), Refusal> {
        let pos = [
            self.inertia_kg_m2,
            self.lever_gain_nm_per_n,
            self.friction_reg_rad_s,
            self.stop_rad,
            self.stop_stiffness_nm_per_rad,
        ];
        let nonneg = [self.coulomb_nm, self.viscous_nm_s, self.stop_damping_nm_s];
        if !(pos.iter().all(|v| v.is_finite() && *v > 0.0)
            && nonneg.iter().all(|v| v.is_finite() && *v >= 0.0))
        {
            return Err(Refusal {
                code: "canard-mech-invalid",
                message: format!("{self:?}"),
                ranked_repairs: vec![
                    "positive inertia/gain/reg/stop; non-negative friction".into(),
                ],
            });
        }
        if !self.hinge_axis_pct_chord.is_finite()
            || self.hinge_axis_pct_chord < HINGE_AXIS_PRIOR_PCT.0
            || self.hinge_axis_pct_chord > HINGE_AXIS_PRIOR_PCT.1
        {
            return Err(Refusal {
                code: "hinge-axis-outside-prior",
                message: format!(
                    "hinge axis {}% chord outside the declared prior {:?} — point guesses \
                     outside the band are forbidden (canard-mechanics-v1 doctrine)",
                    self.hinge_axis_pct_chord, HINGE_AXIS_PRIOR_PCT
                ),
                ranked_repairs: vec![
                    "sample inside the prior; widening the prior is a dossier change".into(),
                ],
            });
        }
        Ok(())
    }

    /// Friction, stop resistance, stop damping, and their rate derivative.
    fn resistance(&self, start: f64, rate: f64, dt: f64) -> (f64, f64, f64, f64) {
        let displacement = dt * rate;
        let end = start + displacement;
        let positive = stop_segment(start - self.stop_rad, end - self.stop_rad, displacement);
        let negative = stop_segment(-start - self.stop_rad, -end - self.stop_rad, -displacement);
        let spring = self.stop_stiffness_nm_per_rad * (positive.0 - negative.0);
        let damping = self.stop_damping_nm_s * (positive.1 + negative.1) * rate;
        let tanh = det::tanh(rate / self.friction_reg_rad_s);
        let friction = self.coulomb_nm * tanh + self.viscous_nm_s * rate;
        let friction_slope = self.coulomb_nm * ((1.0 - tanh * tanh) / self.friction_reg_rad_s);
        let spring_slope = self.stop_stiffness_nm_per_rad * dt * (positive.2 + negative.2);
        let damping_slope = if end.abs() > self.stop_rad {
            self.stop_damping_nm_s
        } else {
            0.0
        };
        (
            friction,
            spring + damping,
            damping,
            friction_slope + self.viscous_nm_s + spring_slope + damping_slope,
        )
    }

    /// One passive midpoint/discrete-gradient step. `m_aero_nm` is the coupled-solution
    /// hinge load (E4.2b) held over the step; `pilot_force_n` positive =
    /// pull back = positive pitch command (control-signs-v1).
    ///
    /// # Errors
    /// Admission refusals; `mech-state-invalid` (non-finite state or
    /// inputs, dt ≤ 0); `mech-step-non-finite` for derived overflow;
    /// `mech-step-unresolved` if the bounded scalar solve does not converge.
    pub fn step(
        &self,
        state: MechState,
        m_aero_nm: f64,
        pilot_force_n: f64,
        dt_s: f64,
    ) -> Result<(MechState, StepReceipt), Refusal> {
        self.admit()?;
        if !(state.delta_rad.is_finite()
            && state.rate_rad_s.is_finite()
            && m_aero_nm.is_finite()
            && pilot_force_n.is_finite()
            && dt_s.is_finite()
            && dt_s > 0.0)
        {
            return Err(Refusal {
                code: "mech-state-invalid",
                message: format!("state {state:?}, m_aero {m_aero_nm:?}, dt {dt_s:?}"),
                ranked_repairs: vec!["check the tick inputs".into()],
            });
        }
        // q1=q0+dt*vm, v1=2*vm-v0 reduces the coupled solve to one
        // strictly increasing residual: (2I/dt)(vm-v0)+R(vm)-drive.
        // A fixed number of Picard substitutions is not an implicit solve:
        // its contraction fails exactly at stiff friction/stops.
        let inertia_scale = 2.0 * (self.inertia_kg_m2 / dt_s);
        let drive = m_aero_nm + pilot_force_n * self.lever_gain_nm_per_n;
        let initial_inertia = inertia_scale * state.rate_rad_s;
        let at_zero = self.resistance(state.delta_rad, 0.0, dt_s);
        let residual_zero = at_zero.1 - initial_inertia - drive;
        // R is monotone, so the inertial-only correction from zero is a
        // valid opposite endpoint. No unbounded bracket search is needed.
        let bound = -residual_zero / inertia_scale;
        if ![inertia_scale, initial_inertia, drive, residual_zero, bound]
            .iter()
            .all(|value| value.is_finite())
            || inertia_scale <= 0.0
        {
            return Err(step_refusal("mech-step-non-finite"));
        }
        let (mut lower, mut upper) = (bound.min(0.0), bound.max(0.0));
        let mut rate = 0.0;
        for _ in 0..MIDPOINT_ITERATIONS {
            let (friction, stop, damping, slope) = self.resistance(state.delta_rad, rate, dt_s);
            let inertia = inertia_scale * rate;
            let net = drive - friction - stop;
            let residual = inertia - initial_inertia - net;
            if ![friction, stop, damping, inertia, net, residual]
                .iter()
                .all(|value| value.is_finite())
            {
                return Err(step_refusal("mech-step-non-finite"));
            }
            let scale = inertia
                .abs()
                .max(initial_inertia.abs())
                .max(drive.abs())
                .max(friction.abs())
                .max(stop.abs())
                .max(f64::MIN_POSITIVE);
            if residual.abs() <= 128.0 * f64::EPSILON * scale {
                let displacement = dt_s * rate;
                let next = MechState {
                    delta_rad: state.delta_rad + displacement,
                    rate_rad_s: 2.0 * rate - state.rate_rad_s,
                };
                let receipt = StepReceipt {
                    friction_dissipation_j: friction * displacement,
                    stop_dissipation_j: damping * displacement,
                    drive_work_j: drive * displacement,
                    stop_torque_nm: -stop,
                    net_torque_nm: net,
                };
                if ![
                    next.delta_rad,
                    next.rate_rad_s,
                    receipt.friction_dissipation_j,
                    receipt.stop_dissipation_j,
                    receipt.drive_work_j,
                ]
                .iter()
                .all(|value| value.is_finite())
                {
                    return Err(step_refusal("mech-step-non-finite"));
                }
                return Ok((next, receipt));
            }
            if residual < 0.0 {
                lower = rate;
            } else {
                upper = rate;
            }
            let newton = rate - residual / (inertia_scale + slope);
            rate = if newton.is_finite() && newton >= lower && newton <= upper && newton != rate {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
        }
        Err(step_refusal("mech-step-unresolved"))
    }
}

/// Discrete gradient of 1/2 max(x,0)^2, contact fraction, and its
/// derivative with respect to displacement. Factoring the energy difference
/// avoids subtracting nearly equal squared penetrations. The fractions also
/// integrate the stop damper along the step's linear configuration path.
fn stop_segment(start: f64, end: f64, displacement: f64) -> (f64, f64, f64) {
    if start <= 0.0 && end <= 0.0 {
        (0.0, 0.0, 0.0)
    } else if start > 0.0 && end > 0.0 {
        (0.5 * start + 0.5 * end, 1.0, 0.5)
    } else if start > 0.0 {
        let fraction = (-start / displacement).clamp(0.0, 1.0);
        (0.5 * start * fraction, fraction, 0.5 * fraction * fraction)
    } else {
        let fraction = (end / displacement).clamp(0.0, 1.0);
        (
            0.5 * end * fraction,
            fraction,
            fraction * (1.0 - 0.5 * fraction),
        )
    }
}

fn step_refusal(code: &'static str) -> Refusal {
    Refusal {
        code,
        message: "canard friction/stop step has no finite converged result".into(),
        ranked_repairs: vec![
            "reduce the time step or rescale the declared mechanism inputs".into(),
        ],
    }
}
