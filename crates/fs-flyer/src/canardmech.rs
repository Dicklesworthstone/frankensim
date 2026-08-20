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
//! Integration: deterministic implicit midpoint (stiff-safe for the
//! stop spring), fixed iteration count — bitwise repeatable.

use crate::Refusal;
use fs_math::det;

/// Hinge-axis prior band [% chord from LE] (canard-mechanics-v1
/// unknown_priors — admission refuses outside it).
pub const HINGE_AXIS_PRIOR_PCT: (f64, f64) = (25.0, 50.0);

/// Fixed midpoint iterations (deterministic).
pub const MIDPOINT_ITERATIONS: u32 = 12;

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
    /// Stop torque at the midpoint [N·m] (0 unless penetrating).
    pub stop_torque_nm: f64,
    /// Net torque at the midpoint [N·m].
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

    /// Torque at a state (aero + pilot held constant over the step).
    fn torque(&self, delta: f64, rate: f64, m_aero: f64, pilot_n: f64) -> (f64, f64, f64) {
        let friction =
            self.coulomb_nm * det::tanh(rate / self.friction_reg_rad_s) + self.viscous_nm_s * rate;
        let stop = if delta > self.stop_rad {
            -self.stop_stiffness_nm_per_rad * (delta - self.stop_rad)
                - self.stop_damping_nm_s * rate
        } else if delta < -self.stop_rad {
            -self.stop_stiffness_nm_per_rad * (delta + self.stop_rad)
                - self.stop_damping_nm_s * rate
        } else {
            0.0
        };
        let net = m_aero + pilot_n * self.lever_gain_nm_per_n - friction + stop;
        (net, friction, stop)
    }

    /// One implicit-midpoint step. `m_aero_nm` is the coupled-solution
    /// hinge load (E4.2b) held over the step; `pilot_force_n` positive =
    /// pull back = positive pitch command (control-signs-v1).
    ///
    /// # Errors
    /// Admission refusals; `mech-state-invalid` (non-finite state or
    /// inputs, dt ≤ 0).
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
        // Implicit midpoint: x1 = x0 + dt*f((x0+x1)/2), fixed-point
        // iterated a FIXED count (deterministic; converges for dt well
        // under the stop-spring period, and the fixed count makes the
        // result bitwise stable regardless).
        let (mut d1, mut r1) = (state.delta_rad, state.rate_rad_s);
        let mut mid = (0.0, 0.0, 0.0);
        for _ in 0..MIDPOINT_ITERATIONS {
            let dm = 0.5 * (state.delta_rad + d1);
            let rm = 0.5 * (state.rate_rad_s + r1);
            mid = self.torque(dm, rm, m_aero_nm, pilot_force_n);
            d1 = state.delta_rad + dt_s * rm;
            r1 = state.rate_rad_s + dt_s * mid.0 / self.inertia_kg_m2;
        }
        let rm = 0.5 * (state.rate_rad_s + r1);
        let receipt = StepReceipt {
            friction_dissipation_j: mid.1 * rm * dt_s,
            stop_torque_nm: mid.2,
            net_torque_nm: mid.0,
        };
        Ok((
            MechState {
                delta_rad: d1,
                rate_rad_s: r1,
            },
            receipt,
        ))
    }
}
