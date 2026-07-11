//! Unsteady free wakes in 2D: point vortices shed at the trailing edge
//! each step with strength set by KELVIN's theorem (total circulation
//! of bound + wake is conserved at zero), convected by the freestream
//! plus the regularized induced flow of every other wake vortex and
//! the bound vortex — the flapping-gait screening loop's kernel shape.
//! The shipped v1 path is a capped direct sum; fs-fmm acceleration is a
//! recorded successor and is not claimed here.
//!
//! The impulsive-start fixture is the classic: bound circulation grows
//! from roughly HALF its pressure-derived screening asymptote toward that value
//! (the Wagner transient), and the shed sheet rolls up without
//! blowing up — stability and determinism are asserted, shapes are
//! ledgered.

use crate::BemError;
use crate::panel2d::{Airfoil2d, solve};
use std::fmt::Write as _;

/// Largest admitted number of wake steps/vortices in one v1 state.
pub const MAX_WAKE_STEPS: usize = 1_000_000;
/// Direct v1 convection ceiling; larger wakes require the future FMM path.
pub const MAX_DIRECT_WAKE_VORTICES: usize = 1_024;
/// Largest admitted all-pairs induced-velocity pass per step.
pub const MAX_DIRECT_WAKE_PAIR_WORK: usize = MAX_DIRECT_WAKE_VORTICES * MAX_DIRECT_WAKE_VORTICES;
const MAX_ABS_WAKE_ALPHA: f64 = std::f64::consts::FRAC_PI_2;
const MAX_WAKE_DT: f64 = 1.0;
const MIN_WAKE_CORE: f64 = 1.0e-9;
const MAX_WAKE_CORE: f64 = 1.0;

/// One wake vortex.
#[derive(Debug, Clone, Copy)]
pub struct WakeVortex {
    /// Position.
    pub pos: [f64; 2],
    /// Nondimensional circulation (unit chord and freestream).
    pub gamma: f64,
}

/// One step's ledger row.
#[derive(Debug, Clone, Copy)]
pub struct WakeStep {
    /// Nondimensional time (chord/freestream units).
    pub t: f64,
    /// Bound circulation.
    pub bound: f64,
    /// Wake vortex count.
    pub vortices: usize,
    /// Peak induced speed among wake vortices.
    pub peak_speed: f64,
}

/// The unsteady simulation state.
#[derive(Debug, Clone)]
pub struct WakeSim {
    alpha: f64,
    dt: f64,
    core2: f64,
    wake: Vec<WakeVortex>,
    history: Vec<WakeStep>,
    steady_gamma: f64,
    te: [f64; 2],
}

impl WakeSim {
    /// Impulsive start at `alpha`, time step `dt`, regularization core
    /// radius `core`.
    pub fn new(foil: &Airfoil2d, alpha: f64, dt: f64, core: f64) -> Result<WakeSim, BemError> {
        if !alpha.is_finite() || alpha.abs() > MAX_ABS_WAKE_ALPHA {
            return Err(BemError::InvalidScalar {
                name: "wake angle of attack",
                value: alpha,
                requirement: "finite radians with |alpha| <= pi/2",
            });
        }
        if !dt.is_finite() || dt <= 0.0 || dt > MAX_WAKE_DT {
            return Err(BemError::InvalidScalar {
                name: "wake time step",
                value: dt,
                requirement: "finite and in (0, 1] nondimensional time",
            });
        }
        if !core.is_finite() || !(MIN_WAKE_CORE..=MAX_WAKE_CORE).contains(&core) {
            return Err(BemError::InvalidScalar {
                name: "wake regularization radius",
                value: core,
                requirement: "finite and in [1e-9, 1] chord lengths",
            });
        }
        let steady = solve(foil, alpha)?;
        let te = foil.nodes()[0];
        let sim = WakeSim {
            alpha,
            dt,
            core2: core * core,
            wake: Vec::new(),
            history: Vec::new(),
            // Screening-model normalization for unit chord/speed. This is an
            // inferred circulation scale, not a second exact Kutta-Joukowski
            // evaluation of the pressure-integrated panel result.
            steady_gamma: steady.cl / 2.0,
            te,
        };
        sim.validate_state()?;
        Ok(sim)
    }

    /// Steady screening-model circulation scale (the transient asymptote).
    #[must_use]
    pub fn steady_circulation(&self) -> f64 {
        self.steady_gamma
    }

    /// Read-only shed vortices.
    #[must_use]
    pub fn wake(&self) -> &[WakeVortex] {
        &self.wake
    }

    /// Read-only bound-circulation history.
    #[must_use]
    pub fn history(&self) -> &[WakeStep] {
        &self.history
    }

    fn induced(&self, wake: &[WakeVortex], p: [f64; 2], bound: f64, quarter: [f64; 2]) -> [f64; 2] {
        // Platform trig is a cross-ISA hazard in solver paths (6ure).
        let mut v = [fs_math::det::cos(self.alpha), fs_math::det::sin(self.alpha)];
        let two_pi = std::f64::consts::TAU;
        // Bound vortex lumped at the quarter chord (screening model).
        let dx = p[0] - quarter[0];
        let dy = p[1] - quarter[1];
        let r2 = (dx * dx + dy * dy + self.core2).max(1e-12);
        v[0] += bound * dy / (two_pi * r2);
        v[1] += -bound * dx / (two_pi * r2);
        for w in wake {
            let dx = p[0] - w.pos[0];
            let dy = p[1] - w.pos[1];
            let r2 = (dx * dx + dy * dy + self.core2).max(1e-12);
            v[0] += w.gamma * dy / (two_pi * r2);
            v[1] += -w.gamma * dx / (two_pi * r2);
        }
        v
    }

    /// Advance one step: quasi-steady bound circulation ramps toward
    /// the Kutta value against the wake's downwash; the DEFICIT is
    /// shed at the trailing edge (Kelvin), and the wake convects with
    /// the induced flow (forward Euler at fixture scale).
    pub fn step(&mut self) -> Result<(), BemError> {
        self.validate_state()?;
        if self.history.len() >= MAX_WAKE_STEPS {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "wake steps",
                requested: self.history.len().saturating_add(1),
                max: MAX_WAKE_STEPS,
            });
        }
        let quarter = [0.25, 0.0];
        // Quasi-steady bound circulation reduced by wake downwash at
        // the three-quarter-chord control point (thin-airfoil model).
        let cp = [0.75, 0.0];
        // Downwash = the wake's vertical induced velocity at the
        // control point: v_y = Σ −γ·dx/(2π r²).
        let mut wake_wash = 0.0;
        let two_pi = std::f64::consts::TAU;
        for w in &self.wake {
            let dx = cp[0] - w.pos[0];
            let dy = cp[1] - w.pos[1];
            let r2 = (dx * dx + dy * dy + self.core2).max(1e-12);
            wake_wash += -w.gamma * dx / (two_pi * r2);
        }
        if !wake_wash.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "wake downwash",
            });
        }
        // Effective incidence: α + downwash/U (small-angle).
        let bound_prev = self.history.last().map_or(0.0, |s| s.bound);
        // Guard the degenerate zero-incidence case: sin α = 0 makes the ratio
        // 0/0 = NaN, which would poison the whole screening simulation. At zero
        // incidence this model carries no bound circulation to scale.
        let sin_alpha = fs_math::det::sin(self.alpha);
        let bound_target = if sin_alpha.abs() < 1e-12 {
            0.0
        } else {
            self.steady_gamma * fs_math::det::sin(self.alpha + wake_wash) / sin_alpha
        };
        // First-order relaxation with the shed vortex carrying the
        // difference (Kelvin: dΓ_bound = −Γ_shed).
        let bound = bound_prev + 0.5 * (bound_target - bound_prev);
        let shed = bound_prev - bound;
        if !bound.is_finite() || !shed.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "wake circulation update",
            });
        }

        let next_count = self
            .wake
            .len()
            .checked_add(if shed != 0.0 { 1 } else { 0 })
            .ok_or(BemError::WorkEnvelopeExceeded {
                operation: "wake vortices",
                requested: usize::MAX,
                max: MAX_WAKE_STEPS,
            })?;
        validate_direct_wake_work(next_count)?;
        let mut candidate = Vec::new();
        candidate
            .try_reserve_exact(next_count)
            .map_err(|_| BemError::AllocationFailed {
                operation: "wake candidate state",
            })?;
        candidate.extend_from_slice(&self.wake);
        if shed.abs() > 0.0 {
            candidate.push(WakeVortex {
                pos: [self.te[0] + 0.3 * self.dt, self.te[1]],
                gamma: shed,
            });
        }

        // Build the whole next state before committing any mutation.
        let mut snapshot = Vec::new();
        snapshot
            .try_reserve_exact(next_count)
            .map_err(|_| BemError::AllocationFailed {
                operation: "wake velocity snapshot",
            })?;
        for vortex in &candidate {
            let velocity = self.induced(&candidate, vortex.pos, bound, quarter);
            if velocity.iter().any(|component| !component.is_finite()) {
                return Err(BemError::NonFiniteResult {
                    operation: "wake induced velocity",
                });
            }
            snapshot.push(velocity);
        }
        let mut peak = 0.0f64;
        for (w, v) in candidate.iter_mut().zip(&snapshot) {
            w.pos[0] += self.dt * v[0];
            w.pos[1] += self.dt * v[1];
            if w.pos.iter().any(|component| !component.is_finite()) {
                return Err(BemError::NonFiniteResult {
                    operation: "wake convection",
                });
            }
            peak = peak.max(fs_math::det::sqrt(v[0].mul_add(v[0], v[1] * v[1])));
        }
        #[allow(clippy::cast_precision_loss)]
        let t = self.history.len() as f64 * self.dt + self.dt;
        if !t.is_finite() || !peak.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "wake history row",
            });
        }
        self.history
            .try_reserve(1)
            .map_err(|_| BemError::AllocationFailed {
                operation: "wake history",
            })?;
        self.wake = candidate;
        self.history.push(WakeStep {
            t,
            bound,
            vortices: self.wake.len(),
            peak_speed: peak,
        });
        Ok(())
    }

    /// Canonical JSON ledger rows for every `stride`-th step.
    pub fn trace_json(&self, stride: usize) -> Result<String, BemError> {
        self.validate_state()?;
        if stride == 0 {
            return Err(BemError::InvalidTraceStride);
        }
        let rows = self.history.len().div_ceil(stride);
        let reserve = rows
            .checked_mul(128)
            .and_then(|bytes| bytes.checked_add(192))
            .ok_or(BemError::WorkEnvelopeExceeded {
                operation: "wake trace bytes",
                requested: usize::MAX,
                max: MAX_WAKE_STEPS * 128 + 192,
            })?;
        let mut s = String::new();
        s.try_reserve(reserve)
            .map_err(|_| BemError::AllocationFailed {
                operation: "wake JSON trace",
            })?;
        let core = fs_math::det::sqrt(self.core2);
        write!(
            s,
            "{{\"schema\":\"fs-bem.wake-trace.v1\",\"model\":\"inviscid-screening\",\"alpha_rad\":{},\"dt_nondimensional\":{},\"core_radius_chord\":{},\"stride\":{},\"rows\":[",
            self.alpha, self.dt, core, stride
        )
        .map_err(|_| BemError::AllocationFailed {
            operation: "wake JSON trace",
        })?;
        for (i, st) in self.history.iter().enumerate() {
            if i % stride != 0 {
                continue;
            }
            write!(
                s,
                "{{\"time_nondimensional\":{},\"bound_circulation_nondimensional\":{},\"vortices\":{},\"peak_speed_freestream_ratio\":{}}},",
                st.t, st.bound, st.vortices, st.peak_speed
            )
            .map_err(|_| BemError::AllocationFailed {
                operation: "wake JSON trace",
            })?;
        }
        if s.ends_with(',') {
            s.pop();
        }
        s.push_str("]}");
        Ok(s)
    }

    fn validate_state(&self) -> Result<(), BemError> {
        if !self.alpha.is_finite()
            || self.alpha.abs() > MAX_ABS_WAKE_ALPHA
            || !self.dt.is_finite()
            || self.dt <= 0.0
            || self.dt > MAX_WAKE_DT
            || !self.core2.is_finite()
            || !(MIN_WAKE_CORE * MIN_WAKE_CORE..=MAX_WAKE_CORE * MAX_WAKE_CORE)
                .contains(&self.core2)
            || !self.steady_gamma.is_finite()
            || self.te.iter().any(|value| !value.is_finite())
        {
            return Err(BemError::NonFiniteResult {
                operation: "wake configuration",
            });
        }
        if self.history.len() > MAX_WAKE_STEPS {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "wake history",
                requested: self.history.len(),
                max: MAX_WAKE_STEPS,
            });
        }
        if self.wake.len() > MAX_DIRECT_WAKE_VORTICES {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "direct wake vortices",
                requested: self.wake.len(),
                max: MAX_DIRECT_WAKE_VORTICES,
            });
        }
        if self.wake.iter().any(|vortex| {
            !vortex.gamma.is_finite() || vortex.pos.iter().any(|value| !value.is_finite())
        }) || self.history.iter().any(|step| {
            !step.t.is_finite()
                || !step.bound.is_finite()
                || !step.peak_speed.is_finite()
                || step.vortices > MAX_DIRECT_WAKE_VORTICES
        }) {
            return Err(BemError::NonFiniteResult {
                operation: "wake state",
            });
        }
        if self.history.last().map_or(!self.wake.is_empty(), |step| {
            step.vortices != self.wake.len()
        }) {
            return Err(BemError::NonFiniteResult {
                operation: "wake state bookkeeping",
            });
        }
        Ok(())
    }
}

fn validate_direct_wake_work(vortices: usize) -> Result<usize, BemError> {
    if vortices > MAX_DIRECT_WAKE_VORTICES {
        return Err(BemError::WorkEnvelopeExceeded {
            operation: "direct wake vortices",
            requested: vortices,
            max: MAX_DIRECT_WAKE_VORTICES,
        });
    }
    let pair_work = vortices
        .checked_mul(vortices)
        .ok_or(BemError::WorkEnvelopeExceeded {
            operation: "direct wake induced-velocity pairs",
            requested: usize::MAX,
            max: MAX_DIRECT_WAKE_PAIR_WORK,
        })?;
    if pair_work > MAX_DIRECT_WAKE_PAIR_WORK {
        return Err(BemError::WorkEnvelopeExceeded {
            operation: "direct wake induced-velocity pairs",
            requested: pair_work,
            max: MAX_DIRECT_WAKE_PAIR_WORK,
        });
    }
    Ok(pair_work)
}

#[cfg(test)]
mod tests {
    use super::{MAX_DIRECT_WAKE_PAIR_WORK, MAX_DIRECT_WAKE_VORTICES, validate_direct_wake_work};

    #[test]
    fn direct_wake_work_boundary_is_admitted_and_limit_plus_one_is_refused() {
        assert_eq!(
            validate_direct_wake_work(MAX_DIRECT_WAKE_VORTICES).expect("boundary is admitted"),
            MAX_DIRECT_WAKE_PAIR_WORK
        );
        assert!(validate_direct_wake_work(MAX_DIRECT_WAKE_VORTICES + 1).is_err());
    }
}
