//! ReducedAeroelasticWarp (bead wf-root-guzez.5.15, E4.6b0). Plan
//! Round-2 revision: wing-warping is CABLE-DRIVEN TORSION of a compliant
//! structure — a prescribed-kinematic twist overstates control power at
//! speed. This reduced model supports LATERAL-CONTROL claims only
//! (loaded twist, effective control power, slack risk) and explicitly
//! does NOT provide structural margins (no stress, no failure loads —
//! any Validated V-12 lateral claim requires at least this mode).
//!
//! Model: per-strip commanded twist θ_cmd = δw·basis (antisymmetric,
//! control-signs-v1: +δw = +roll = right wing down); the loaded twist
//! solves the per-strip torsion balance
//!
//!   θ_L = θ_cmd − compliance · q·c²·w·clα·(α₀ + θ_L)
//!
//! (closed form — the aeroelastic effectiveness η = 1/(1 + compliance·
//! q·c²·w·clα) emerges q-dependent). Wire-slack bound: the cables
//! cannot push — a load-induced twist DEFICIT beyond the slack bound
//! means a wire has gone slack; the diagnostic fires per strip. The
//! deficit is measured RELATIVE to the zero-command trim washout (the
//! rigging is taut at trim load), so only the command-proportional
//! deficit k/(1+k)·|θ_cmd| carries slack risk.
//! Optional first-order actuation lag with the exact exponential update.

use crate::Refusal;
use fs_math::det;

/// Slack-bound cap [rad] (absurd-input guard).
pub const MAX_SLACK_BOUND_RAD: f64 = 0.5;

/// The reduced aeroelastic warp model.
#[derive(Clone, Debug, PartialEq)]
pub struct ReducedAeroelasticWarp {
    /// Antisymmetric twist basis per wing strip [rad/rad of δw]
    /// (+ on the LEFT half, − on the RIGHT: +δw rolls right-wing-down).
    pub basis: Vec<f64>,
    /// Per-strip torsional compliance [rad/(N·m)] (Estimated).
    pub compliance_rad_per_nm: Vec<f64>,
    /// Strip chords [m].
    pub chord_m: Vec<f64>,
    /// Strip widths [m].
    pub width_m: Vec<f64>,
    /// Section lift slope [1/rad].
    pub cl_alpha: f64,
    /// Wire slack bound [rad]: a warp-relative twist deficit beyond
    /// this means a slack cable.
    pub slack_bound_rad: f64,
    /// Optional first-order actuation lag time constant [s].
    pub lag_tau_s: Option<f64>,
}

/// One strip's loaded-twist line (per-item oracle surface).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StripWarp {
    /// Commanded twist [rad].
    pub commanded_rad: f64,
    /// Loaded (achieved) twist [rad].
    pub loaded_rad: f64,
    /// Aeroelastic effectiveness η ∈ (0, 1].
    pub effectiveness: f64,
    /// Slack margin [rad] (slack_bound − warp-relative deficit;
    /// < 0 ⇒ AT RISK).
    pub slack_margin_rad: f64,
}

/// The loaded-warp evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct LoadedWarpReport {
    /// Per-strip lines.
    pub strips: Vec<StripWarp>,
    /// Strips whose slack margin is negative.
    pub slack_risk_strips: Vec<usize>,
    /// Twist-weighted mean effectiveness.
    pub mean_effectiveness: f64,
}

/// First-order actuation lag state (exact exponential update).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WarpLagState {
    /// Current effective command [rad].
    pub delta_w_rad: f64,
}

impl ReducedAeroelasticWarp {
    /// The registered Wright v1 model: 16 wing strips (two planes of 8),
    /// linear root→tip antisymmetric basis, uniform Estimated compliance
    /// (set for ~0.78 cruise effectiveness — the warp was effective at
    /// Wright speeds; the aeroelastic loss is real but not dominant).
    #[must_use]
    pub fn wright_v1() -> Self {
        let mut basis = Vec::with_capacity(16);
        for _plane in 0..2 {
            for s in 0..8 {
                // Strip centers at y/b = −0.4375 … +0.4375 (8 strips):
                // basis = −y/(b/2) so the LEFT half (y<0) is positive.
                let y_frac = (s as f64 + 0.5) / 8.0 * 2.0 - 1.0; // −0.875..0.875
                basis.push(-y_frac);
            }
        }
        ReducedAeroelasticWarp {
            basis,
            compliance_rad_per_nm: vec![6.0e-5; 16],
            chord_m: vec![1.981; 16],
            width_m: vec![12.29 / 8.0; 16],
            cl_alpha: core::f64::consts::TAU,
            slack_bound_rad: 0.035,
            lag_tau_s: Some(0.12),
        }
    }

    /// Admit the model.
    ///
    /// # Errors
    /// `warp-model-invalid` (length mismatches, non-finite entries,
    /// negative compliance, non-positive chords/widths/slope, slack
    /// bound outside (0, [`MAX_SLACK_BOUND_RAD`]] — tested at cap AND
    /// one ulp past, non-positive lag).
    pub fn admit(&self) -> Result<(), Refusal> {
        let n = self.basis.len();
        let ok = n > 0
            && self.compliance_rad_per_nm.len() == n
            && self.chord_m.len() == n
            && self.width_m.len() == n
            && self.basis.iter().all(|v| v.is_finite() && v.abs() <= 1.0)
            && self
                .compliance_rad_per_nm
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0)
            && self.chord_m.iter().all(|v| v.is_finite() && *v > 0.0)
            && self.width_m.iter().all(|v| v.is_finite() && *v > 0.0)
            && self.cl_alpha.is_finite()
            && self.cl_alpha > 0.0
            && self.slack_bound_rad.is_finite()
            && self.slack_bound_rad > 0.0
            && self.slack_bound_rad <= MAX_SLACK_BOUND_RAD
            && self.lag_tau_s.is_none_or(|t| t.is_finite() && t > 0.0);
        if !ok {
            return Err(Refusal {
                code: "warp-model-invalid",
                message: format!(
                    "n {n}; lens {}/{}/{}; slack {:?}; lag {:?}",
                    self.compliance_rad_per_nm.len(),
                    self.chord_m.len(),
                    self.width_m.len(),
                    self.slack_bound_rad,
                    self.lag_tau_s
                ),
                ranked_repairs: vec![
                    "equal-length finite vectors; |basis| <= 1; compliance >= 0; slack in \
                     (0, 0.5]; positive lag"
                        .into(),
                ],
            });
        }
        Ok(())
    }

    /// Evaluate the loaded warp at a flight state.
    ///
    /// # Errors
    /// Admission refusals; `warp-state-invalid` (non-finite inputs,
    /// non-positive q, |δw| > π/4).
    pub fn evaluate(
        &self,
        delta_w_rad: f64,
        q_pa: f64,
        alpha0_rad: f64,
    ) -> Result<LoadedWarpReport, Refusal> {
        self.admit()?;
        if !(delta_w_rad.is_finite()
            && delta_w_rad.abs() <= core::f64::consts::FRAC_PI_4
            && q_pa.is_finite()
            && q_pa > 0.0
            && alpha0_rad.is_finite())
        {
            return Err(Refusal {
                code: "warp-state-invalid",
                message: format!("dw {delta_w_rad:?}, q {q_pa:?}, alpha {alpha0_rad:?}"),
                ranked_repairs: vec!["stay inside the warp travel and a physical q".into()],
            });
        }
        let mut strips = Vec::with_capacity(self.basis.len());
        let mut risk = Vec::new();
        let (mut eff_sum, mut w_sum) = (0.0f64, 0.0f64);
        for i in 0..self.basis.len() {
            let cmd = delta_w_rad * self.basis[i];
            let k_aero = self.compliance_rad_per_nm[i]
                * q_pa
                * self.chord_m[i]
                * self.chord_m[i]
                * self.width_m[i]
                * self.cl_alpha;
            // θ_L = θ_cmd − k_aero·(α₀ + θ_L)  ⇒  closed form:
            let loaded = (cmd - k_aero * alpha0_rad) / (1.0 + k_aero);
            let effectiveness = 1.0 / (1.0 + k_aero);
            // Slack risk: WARP-RELATIVE deficit only (rigging taut at
            // the trim load; the alpha0 washout exists at zero command).
            let deficit = (k_aero / (1.0 + k_aero)) * cmd.abs();
            let slack_margin = self.slack_bound_rad - deficit;
            if slack_margin < 0.0 {
                risk.push(i);
            }
            let w = cmd.abs();
            eff_sum += effectiveness * w;
            w_sum += w;
            strips.push(StripWarp {
                commanded_rad: cmd,
                loaded_rad: loaded,
                effectiveness,
                slack_margin_rad: slack_margin,
            });
        }
        let mean_effectiveness = if w_sum > 0.0 {
            eff_sum / w_sum
        } else {
            strips.iter().map(|s| s.effectiveness).sum::<f64>() / strips.len() as f64
        };
        Ok(LoadedWarpReport {
            strips,
            slack_risk_strips: risk,
            mean_effectiveness,
        })
    }

    /// Advance the optional actuation lag by the EXACT exponential:
    /// x ← cmd + (x − cmd)·e^(−dt/τ). Without a lag constant this is
    /// the identity on the command.
    ///
    /// # Errors
    /// `warp-state-invalid` (non-finite, dt ≤ 0).
    pub fn lag_step(
        &self,
        state: WarpLagState,
        command_rad: f64,
        dt_s: f64,
    ) -> Result<WarpLagState, Refusal> {
        if !(state.delta_w_rad.is_finite() && command_rad.is_finite() && dt_s > 0.0) {
            return Err(Refusal {
                code: "warp-state-invalid",
                message: format!("state {state:?}, cmd {command_rad:?}, dt {dt_s:?}"),
                ranked_repairs: vec!["check the tick inputs".into()],
            });
        }
        let out = match self.lag_tau_s {
            None => command_rad,
            Some(tau) => command_rad + (state.delta_w_rad - command_rad) * det::exp(-dt_s / tau),
        };
        Ok(WarpLagState { delta_w_rad: out })
    }
}
