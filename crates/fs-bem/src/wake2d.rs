//! Unsteady free wakes in 2D: point vortices shed at the trailing edge
//! each step with strength set by KELVIN's theorem (total circulation
//! of bound + wake is conserved at zero), convected by the freestream
//! plus the regularized induced flow of every other wake vortex and
//! the bound vortex — the flapping-gait screening loop's kernel shape
//! (fs-fmm accelerates the N-body sum at production scale; the fixture
//! battery runs the direct sum).
//!
//! The impulsive-start fixture is the classic: bound circulation grows
//! from roughly HALF its steady value toward the steady Kutta value
//! (the Wagner transient), and the shed sheet rolls up without
//! blowing up — stability and determinism are asserted, shapes are
//! ledgered.

use crate::panel2d::{Airfoil2d, solve};
use std::fmt::Write as _;

/// One wake vortex.
#[derive(Debug, Clone, Copy)]
pub struct WakeVortex {
    /// Position.
    pub pos: [f64; 2],
    /// Circulation.
    pub gamma: f64,
}

/// One step's ledger row.
#[derive(Debug, Clone, Copy)]
pub struct WakeStep {
    /// Time.
    pub t: f64,
    /// Bound circulation.
    pub bound: f64,
    /// Wake vortex count.
    pub vortices: usize,
    /// Peak induced speed among wake vortices.
    pub peak_speed: f64,
}

/// The unsteady simulation state.
pub struct WakeSim {
    alpha: f64,
    dt: f64,
    core2: f64,
    /// The shed wake.
    pub wake: Vec<WakeVortex>,
    /// Bound circulation history.
    pub history: Vec<WakeStep>,
    steady_gamma: f64,
    te: [f64; 2],
}

impl WakeSim {
    /// Impulsive start at `alpha`, time step `dt`, regularization core
    /// radius `core`.
    #[must_use]
    pub fn new(foil: &Airfoil2d, alpha: f64, dt: f64, core: f64) -> WakeSim {
        let steady = solve(foil, alpha);
        let te = foil.nodes[0];
        WakeSim {
            alpha,
            dt,
            core2: core * core,
            wake: Vec::new(),
            history: Vec::new(),
            steady_gamma: steady.cl / 2.0, // Γ_steady = Cl/2 (unit chord/speed)
            te,
        }
    }

    /// Steady-state total circulation (the Wagner asymptote).
    #[must_use]
    pub fn steady_circulation(&self) -> f64 {
        self.steady_gamma
    }

    fn induced(&self, p: [f64; 2], bound: f64, quarter: [f64; 2]) -> [f64; 2] {
        let mut v = [self.alpha.cos(), self.alpha.sin()];
        let two_pi = std::f64::consts::TAU;
        // Bound vortex lumped at the quarter chord (screening model).
        let dx = p[0] - quarter[0];
        let dy = p[1] - quarter[1];
        let r2 = (dx * dx + dy * dy + self.core2).max(1e-12);
        v[0] += bound * dy / (two_pi * r2);
        v[1] += -bound * dx / (two_pi * r2);
        for w in &self.wake {
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
    pub fn step(&mut self) {
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
        // Effective incidence: α + downwash/U (small-angle).
        let bound_prev = self.history.last().map_or(0.0, |s| s.bound);
        let bound_target = self.steady_gamma * (self.alpha + wake_wash).sin() / self.alpha.sin();
        // First-order relaxation with the shed vortex carrying the
        // difference (Kelvin: dΓ_bound = −Γ_shed).
        let bound = bound_prev + 0.5 * (bound_target - bound_prev);
        let shed = bound_prev - bound;
        if shed.abs() > 0.0 {
            self.wake.push(WakeVortex {
                pos: [self.te[0] + 0.3 * self.dt, self.te[1]],
                gamma: shed,
            });
        }
        // Convect.
        let snapshot: Vec<[f64; 2]> = self
            .wake
            .iter()
            .map(|w| self.induced(w.pos, bound, quarter))
            .collect();
        let mut peak = 0.0f64;
        for (w, v) in self.wake.iter_mut().zip(&snapshot) {
            w.pos[0] += self.dt * v[0];
            w.pos[1] += self.dt * v[1];
            peak = peak.max(v[0].hypot(v[1]));
        }
        #[allow(clippy::cast_precision_loss)]
        let t = self.history.len() as f64 * self.dt + self.dt;
        self.history.push(WakeStep {
            t,
            bound,
            vortices: self.wake.len(),
            peak_speed: peak,
        });
    }

    /// Ledger rows for every `stride`-th step.
    #[must_use]
    pub fn trace_json(&self, stride: usize) -> String {
        let mut s = String::from("[");
        for (i, st) in self.history.iter().enumerate() {
            if i % stride != 0 {
                continue;
            }
            let _ = write!(
                s,
                "{{\"t\":{:.3},\"bound\":{:.5},\"n\":{},\"peak\":{:.3}}},",
                st.t, st.bound, st.vortices, st.peak_speed
            );
        }
        let mut s = s.trim_end_matches(',').to_string();
        s.push(']');
        s
    }
}
