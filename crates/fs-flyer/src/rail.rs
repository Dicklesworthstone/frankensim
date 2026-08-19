//! Launch-rail unilateral constraint + force-based release (bead
//! wf-root-guzez.4.8.1, E3.4-i). Plan §5.1.5: the monorail is a
//! UNILATERAL constrained system — while on the rail the vertical DOF is
//! constrained and the normal reaction N must stay compressive (N ≥ 0,
//! never tensile); release happens when the ADMISSIBLE reaction hits zero
//! AND the free state is separating, sustained for a declared hysteresis
//! tick count. There is NO speed threshold: under a gust the aircraft
//! lifts when the FORCES say so, which the battery's speed-threshold twin
//! proves is a different (wrong) time.
//!
//! Frame: frd/NED (+z down). The rail lies along +x at height z_rail;
//! "upward net force" = negative total z-force. On the rail, the along-x
//! dynamics run free (dolly rolling resistance is a declared load the
//! caller includes); z is pinned and N = max(0, F_down_net) with the
//! complementarity pair (N, gap) both never negative and never both
//! positive.

use crate::Refusal;
use crate::spine::RigidBody;

/// Hysteresis-tick cap (refusals at cap AND cap+1).
pub const MAX_HYSTERESIS_TICKS: u32 = 120;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal { code, message, ranked_repairs: vec![repair.into()] }
}

/// Rail configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RailSpec {
    /// Rail height z (NED, negative = above ground plane origin) [m].
    pub z_rail_m: f64,
    /// Rail length from the start [m] (60 ft = 18.29 for Dec-17).
    pub length_m: f64,
    /// Consecutive separating ticks required before release.
    pub hysteresis_ticks: u32,
}

impl RailSpec {
    /// Validate the spec.
    ///
    /// # Errors
    /// `rail-spec-invalid` (non-finite/zero length; hysteresis above the
    /// cap, tested at cap AND cap+1).
    pub fn admit(&self) -> Result<(), Refusal> {
        if !self.z_rail_m.is_finite() || !self.length_m.is_finite() || self.length_m <= 0.0 {
            return Err(refuse(
                "rail-spec-invalid",
                format!("{self:?}"),
                "finite z, positive length (Dec-17 rail: 18.29 m)",
            ));
        }
        if self.hysteresis_ticks == 0 || self.hysteresis_ticks > MAX_HYSTERESIS_TICKS {
            return Err(refuse(
                "rail-spec-invalid",
                format!("hysteresis {} outside [1, {MAX_HYSTERESIS_TICKS}]", self.hysteresis_ticks),
                "2-6 ticks is the intended band",
            ));
        }
        Ok(())
    }
}

/// The rail phase of a run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RailPhase {
    /// Constrained on the rail (dolly engaged).
    OnRail,
    /// Airborne — the constraint no longer exists (one-way transition).
    Released,
}

/// Per-tick rail report (the V-11a receipt row).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RailTick {
    /// Phase AFTER this tick.
    pub phase: RailPhase,
    /// Applied normal reaction [N] (compressive ≥ 0; 0 once released).
    pub normal_n: f64,
    /// Vertical gap to the rail [m] (0 while on rail).
    pub gap_m: f64,
    /// Along-rail position [m].
    pub x_m: f64,
    /// Along-rail speed [m/s].
    pub vx_mps: f64,
    /// Work done by the reaction this tick [J] (must be 0: N ⊥ motion).
    pub reaction_work_j: f64,
    /// Consecutive separating ticks accumulated.
    pub separating_streak: u32,
    /// End-of-rail overrun (the dolly ran off the end while loaded).
    pub end_of_rail: bool,
}

/// The rail runner: 1-D along-rail state + the unilateral vertical logic.
#[derive(Clone, Debug, PartialEq)]
pub struct RailRun {
    spec: RailSpec,
    phase: RailPhase,
    x_m: f64,
    vx_mps: f64,
    streak: u32,
}

impl RailRun {
    /// Start a run at the rail origin, at rest (the held-on-rail state).
    ///
    /// # Errors
    /// Spec refusals.
    pub fn start(spec: RailSpec) -> Result<RailRun, Refusal> {
        spec.admit()?;
        Ok(RailRun { spec, phase: RailPhase::OnRail, x_m: 0.0, vx_mps: 0.0, streak: 0 })
    }

    /// Current phase.
    #[must_use]
    pub fn phase(&self) -> RailPhase {
        self.phase
    }

    /// Advance one tick while the constraint governs. `force_x_n` is the
    /// net along-rail force (thrust − drag − dolly resistance);
    /// `force_z_n` is the net vertical force INCLUDING weight (+z down,
    /// so weight is +m·g and lift is negative). Returns the tick receipt;
    /// after release the caller hands the state to the free 6-DOF spine.
    ///
    /// # Errors
    /// `non-finite-input`; `rail-already-released` (the transition is
    /// one-way; a released run never re-enters the constraint here —
    /// touchdown is contact's job, E3.4-ii).
    pub fn tick(&mut self, body: &RigidBody, force_x_n: f64, force_z_n: f64, dt_s: f64)
    -> Result<RailTick, Refusal> {
        body.admit()?;
        if !(force_x_n.is_finite() && force_z_n.is_finite() && dt_s.is_finite() && dt_s > 0.0) {
            return Err(refuse(
                "non-finite-input",
                format!("fx {force_x_n}, fz {force_z_n}, dt {dt_s}"),
                "finite forces, positive dt",
            ));
        }
        if self.phase == RailPhase::Released {
            return Err(refuse(
                "rail-already-released",
                "the rail constraint is one-way; the free spine owns the state now".into(),
                "route post-release ticks to spine::step; touchdown is E3.4-ii contact",
            ));
        }
        // Unilateral reaction: the rail can only PUSH (upward, −z). The
        // admissible reaction cancels net downward force; a net upward
        // force cannot be resisted (no tensile rail).
        let normal_n = force_z_n.max(0.0);
        let separating = force_z_n < 0.0;
        self.streak = if separating { self.streak + 1 } else { 0 };
        // Along-rail dynamics (semi-implicit Euler on the 1-D coordinate;
        // the full spine takes over at release).
        let ax = force_x_n / body.mass_kg;
        self.vx_mps += ax * dt_s;
        self.x_m += self.vx_mps * dt_s;
        let end_of_rail = self.x_m > self.spec.length_m;
        if self.streak >= self.spec.hysteresis_ticks {
            self.phase = RailPhase::Released;
        }
        Ok(RailTick {
            phase: self.phase,
            // Once released this tick, the constraint force is zero.
            normal_n: if self.phase == RailPhase::Released { 0.0 } else { normal_n },
            gap_m: 0.0,
            x_m: self.x_m,
            vx_mps: self.vx_mps,
            // Reaction is vertical, motion is horizontal: exactly zero work.
            reaction_work_j: 0.0,
            separating_streak: self.streak,
            end_of_rail,
        })
    }
}
