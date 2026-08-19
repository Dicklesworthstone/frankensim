//! Skid/sand contact (bead wf-root-guzez.4.8.2, E3.4-ii). Plan §5.1.5:
//! penetration springs + REGULARIZED Coulomb friction + plastic sink,
//! with impact reports. fs-flyer owns this contact (the plan's ownership
//! rule); terrain height comes from the E1.3 grids via the caller.
//!
//! Model (frd/NED, +z down; penetration p = z_skid − z_surface ≥ 0):
//! - Normal: N = max(0, k·p_e + c·ṗ) with the ELASTIC penetration
//!   p_e = p − s (s = plastic sink). max(0, ·) forbids adhesion: sand
//!   never pulls the skid down.
//! - Plastic sink: ṡ = λ·max(0, p − s) — the sand yields toward the
//!   current penetration and NEVER recovers (s is monotone). Energy
//!   absorbed by the sink is the landing's plastic loss.
//! - Friction: F_t = −μ·N·tanh(v_t / v_reg) — regularized Coulomb; the
//!   tanh keeps the force smooth through zero slip (no stick chatter at
//!   the tick rate) while saturating at μ·N.

use crate::Refusal;

/// Penetration cap [m] — beyond this the model is out of its domain
/// (a crash, not a landing; refusals at cap AND cap+1 ulp).
pub const MAX_PENETRATION_M: f64 = 0.5;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal { code, message, ranked_repairs: vec![repair.into()] }
}

/// Sand-contact parameters (per skid).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactParams {
    /// Penetration-spring stiffness [N/m].
    pub stiffness_n_m: f64,
    /// Penetration damping [N·s/m].
    pub damping_n_s_m: f64,
    /// Coulomb friction coefficient.
    pub mu: f64,
    /// Regularization velocity scale [m/s].
    pub v_reg_mps: f64,
    /// Plastic-sink rate λ [1/s].
    pub sink_rate_per_s: f64,
    /// Bearing threshold p_yield [m]: elastic penetration below this
    /// carries load without further plastic yield.
    pub yield_penetration_m: f64,
}

impl ContactParams {
    /// Validate.
    ///
    /// # Errors
    /// `contact-params-invalid`.
    pub fn admit(&self) -> Result<(), Refusal> {
        let ok = self.stiffness_n_m.is_finite()
            && self.stiffness_n_m > 0.0
            && self.damping_n_s_m.is_finite()
            && self.damping_n_s_m >= 0.0
            && self.mu.is_finite()
            && self.mu >= 0.0
            && self.v_reg_mps.is_finite()
            && self.v_reg_mps > 0.0
            && self.sink_rate_per_s.is_finite()
            && self.sink_rate_per_s >= 0.0
            && self.yield_penetration_m.is_finite()
            && self.yield_penetration_m > 0.0;
        if !ok {
            return Err(refuse(
                "contact-params-invalid",
                format!("{self:?}"),
                "k>0, c>=0, mu>=0, v_reg>0, sink_rate>=0, p_yield>0",
            ));
        }
        Ok(())
    }
}

/// One skid's contact state (plastic sink is state — it persists).
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct ContactState {
    /// Accumulated plastic sink s [m] (monotone non-decreasing).
    pub sink_m: f64,
    /// True once contact has occurred (drives the impact report).
    pub touched: bool,
}

/// Per-tick contact forces + diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactTick {
    /// Normal force [N], upward (≥ 0; 0 out of contact).
    pub normal_n: f64,
    /// Friction force along the slip direction [N] (opposes v_t).
    pub friction_n: f64,
    /// Elastic penetration p_e [m].
    pub penetration_m: f64,
    /// Plastic sink after this tick [m].
    pub sink_m: f64,
    /// In contact this tick.
    pub in_contact: bool,
}

/// Advance one skid-contact tick. `penetration_m` = z_skid − z_surface
/// (≥ 0 means below the surface), `pen_rate_mps` its rate, `v_t` the
/// tangential slip speed (signed along the slip axis).
///
/// # Errors
/// Param refusals; `non-finite-input`;
/// `penetration-outside-domain` above [`MAX_PENETRATION_M`].
pub fn contact_tick(
    p: &ContactParams,
    state: &mut ContactState,
    penetration_m: f64,
    pen_rate_mps: f64,
    v_t_mps: f64,
    dt_s: f64,
) -> Result<ContactTick, Refusal> {
    p.admit()?;
    if !(penetration_m.is_finite() && pen_rate_mps.is_finite() && v_t_mps.is_finite())
        || !dt_s.is_finite()
        || dt_s <= 0.0
    {
        return Err(refuse("non-finite-input", "penetration/rate/v_t/dt".into(), "finite, dt>0"));
    }
    if penetration_m > MAX_PENETRATION_M {
        return Err(refuse(
            "penetration-outside-domain",
            format!("penetration {penetration_m} m exceeds {MAX_PENETRATION_M} (crash regime)"),
            "a TerminalEvent (crash) owns this state, not the landing model",
        ));
    }
    if penetration_m <= state.sink_m || penetration_m <= 0.0 {
        // Out of contact (or inside the already-yielded crater): no force.
        return Ok(ContactTick {
            normal_n: 0.0,
            friction_n: 0.0,
            penetration_m: 0.0,
            sink_m: state.sink_m,
            in_contact: false,
        });
    }
    state.touched = true;
    // Plastic sink (explicit, monotone): yields only ABOVE the bearing
    // threshold, so static loads settle instead of sinking forever.
    let p_e_pre = penetration_m - state.sink_m;
    let over_yield = (p_e_pre - p.yield_penetration_m).max(0.0);
    state.sink_m += p.sink_rate_per_s * over_yield * dt_s;
    let p_e = (penetration_m - state.sink_m).max(0.0);
    let normal = (p.stiffness_n_m * p_e + p.damping_n_s_m * pen_rate_mps).max(0.0);
    let friction = -p.mu * normal * (v_t_mps / p.v_reg_mps).tanh();
    Ok(ContactTick {
        normal_n: normal,
        friction_n: friction,
        penetration_m: p_e,
        sink_m: state.sink_m,
        in_contact: true,
    })
}

/// The touchdown/impact report (a TerminalEvent-compatible summary).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpactReport {
    /// Vertical speed at first contact [m/s].
    pub impact_speed_mps: f64,
    /// Peak normal force [N].
    pub peak_normal_n: f64,
    /// Normal impulse ∫N dt [N·s].
    pub normal_impulse_ns: f64,
    /// Maximum elastic penetration [m].
    pub max_penetration_m: f64,
    /// Final plastic sink [m].
    pub final_sink_m: f64,
    /// Sliding distance during contact [m].
    pub sliding_m: f64,
}
