//! ReducedLateralBuildUp (bead frankensim-4pa2k; plan §5.1.2 build-up
//! tier under the admitted [`crate::aerowarp::ReducedAeroelasticWarp`]
//! envelope).
//!
//! Scope honesty: this is the REDUCED lateral-control tier. It supports
//! aerodynamic-control and lateral-mode claims (loaded twist, roll and
//! yaw response, adverse-yaw attribution) exactly as the plan sanctions
//! for `ReducedAeroelasticWarp`, and it does NOT provide structural
//! margins, spiral-mode eigenvalues, or validated control power.
//!
//! Model (all reduced constants DECLARED here, never silently tuned):
//! - Loaded twist arrives from `aerowarp` per-strip evaluation; this
//!   module consumes the SIGNED total loaded twist [rad] (positive =
//!   right wing down, control-signs-v1).
//! - Roll: differential lift moment ∝ q·S·b·clα·Θ_L·k_arm against a
//!   strip-theory damping term ∝ q·S·b²·clα·p/(4u); first-order ṗ,
//!   φ̇ = p.
//! - Yaw decomposition (published rows summed by the sim plane; the
//!   view REFUSES a broken split):
//!   - induced-drag differential yaw opposes the warp command when the
//!     rudder linkage is DECOUPLED (the 1901 configuration) — the
//!     historical adverse-yaw sign law;
//!   - rudder yaw is ZERO decoupled, or the declared linkage gain times
//!     the warp command when coupled (1902+);
//!   - profile yaw is reduced yaw damping only (∝ r·b/(2u)).
//! - Envelope: declared attitude/rate caps; a non-finite state or an
//!   attitude-cap crossing refuses with a typed lateral refusal (the
//!   caller closes the run; nothing saturates silently).

use crate::Refusal;

/// Declared roll inertia [kg·m²] (Estimated): strip-model bound for the
/// 340 kg flyer at 12.29 m span, Ix ≈ m·(b/8)²; no source publishes the
/// measured roll inertia.
pub const FLYER_ROLL_INERTIA_KG_M2: f64 = 340.2 * (12.29 / 8.0) * (12.29 / 8.0);

/// Declared yaw inertia [kg·m²] (Estimated): the same strip bound scaled
/// by a b²-dominant axis ratio 4; declared, not sourced.
pub const REDUCED_YAW_INERTIA_KG_M2: f64 = FLYER_ROLL_INERTIA_KG_M2 * 4.0;

/// Declared roll-attitude admission cap [rad]: beyond ~70° the reduced
/// linear lateral model has no claim (post-stall kinematics absent).
pub const PHI_CAP_RAD: f64 = 1.2;

/// Declared yaw-rate cap [rad/s] (absurd-input guard).
pub const R_CAP_RAD_S: f64 = 1.5;

/// Bundle coefficient of the induced-drag differential yaw moment
/// [N·m/rad at the reference dynamic pressure 76.8 Pa ≈ the Dec-17 trim
/// q]. NEGATIVE by the adverse-yaw sign law: with the rudder linkage
/// DECOUPLED, induced-drag yaw opposes the warp roll command. Declared;
/// first-order CL_trim·ΔCL product only.
const K_INDUCED_YAW_NM_PER_RAD_AT_REF_Q: f64 = -42.0;

/// Lateral state of one tick.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct LateralState {
    /// Roll attitude [rad] (+ = right wing down).
    pub phi_rad: f64,
    /// Roll rate [rad/s].
    pub p_rad_s: f64,
    /// Heading [rad] (+ = nose right).
    pub psi_rad: f64,
    /// Yaw rate [rad/s].
    pub r_rad_s: f64,
}

/// Published yaw-moment decomposition for one tick [N·m]. The view's
/// sum-check law: induced + rudder + profile == net, summed in THIS
/// component order so the float result is reproducible bit-for-bit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YawDecompositionNm {
    /// Induced-drag differential yaw moment (adverse when negative for a
    /// positive warp command).
    pub induced_drag_yaw_nm: f64,
    /// Rudder linkage yaw moment (0 when decoupled).
    pub rudder_yaw_nm: f64,
    /// Profile/other yaw moment (reduced damping term).
    pub profile_yaw_nm: f64,
}

impl YawDecompositionNm {
    /// The sim plane's OWN net yaw moment: the exact float sum in the
    /// published component order.
    #[must_use]
    pub fn net(self) -> f64 {
        self.induced_drag_yaw_nm + self.rudder_yaw_nm + self.profile_yaw_nm
    }
}

/// One tick's published lateral row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LateralTick {
    /// Post-integration state.
    pub state: LateralState,
    /// Roll moment applied this tick [N·m] (before integration).
    pub roll_moment_nm: f64,
    /// Yaw-moment decomposition applied this tick [N·m].
    pub yaw: YawDecompositionNm,
    /// Total LOADED twist consumed this tick [rad] (signed).
    pub loaded_twist_rad: f64,
}

/// Rudder linkage configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RudderLinkage {
    /// 1901 configuration: no rudder-to-warp coupling.
    Decoupled,
    /// 1902+ configuration: linked rudder with a declared compensating
    /// gain per rad of warp command (positive = pro-yaw with command,
    /// turning the NET yaw proverse exactly as the closed E7.4b
    /// attribution documents).
    Linked {
        /// Linkage gain [N·m/rad].
        gain_nm_per_rad: f64,
    },
}

/// Reduced lateral model parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LateralModel {
    /// Wing area [m²].
    pub wing_area_m2: f64,
    /// Wing span [m].
    pub wing_span_m: f64,
    /// Section lift slope [1/rad] (shared with the warp model).
    pub cl_alpha: f64,
    /// Differential-lift arm factor [dimensionless]: centroid of the
    /// antisymmetric lift distribution as a fraction of full span
    /// (declared elliptical-like 1/3).
    pub arm_factor: f64,
    /// Rudder linkage mode.
    pub rudder: RudderLinkage,
}

impl LateralModel {
    /// The Wright Flyer reference parameterization (matches the
    /// `aircraft` geometry constants; declared factors above).
    #[must_use]
    pub fn wright_v1(rudder: RudderLinkage) -> Self {
        Self {
            wing_area_m2: 51.0,
            wing_span_m: 12.29,
            cl_alpha: 5.5,
            arm_factor: 1.0 / 3.0,
            rudder,
        }
    }

    /// One 120 Hz lateral step.
    ///
    /// `loaded_twist_rad`: signed total loaded twist from the warp model
    /// (positive = right wing down). `u_mps`: body airspeed. `rho`:
    /// air density.
    ///
    /// # Errors
    /// `lateral-envelope-exceeded` on non-finite inputs/state or a
    /// roll/yaw cap crossing.
    pub fn step(
        &self,
        state: &mut LateralState,
        loaded_twist_rad: f64,
        u_mps: f64,
        rho_kg_m3: f64,
        dt_s: f64,
    ) -> Result<LateralTick, Refusal> {
        if ![loaded_twist_rad, u_mps, rho_kg_m3, dt_s]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(Refusal {
                code: "lateral-envelope-exceeded",
                message: "non-finite lateral input".into(),
                ranked_repairs: vec!["check the warp/atmosphere upstream stages".into()],
            });
        }
        let q_dyn = 0.5 * rho_kg_m3 * u_mps * u_mps;
        let speed_floor = u_mps.max(0.1);
        // Differential-lift roll moment: the signed total twist already
        // carries the antisymmetry of the strip distribution.
        let roll_moment = q_dyn
            * self.wing_area_m2
            * self.wing_span_m
            * self.cl_alpha
            * loaded_twist_rad
            * self.arm_factor
            * 0.5;
        // Strip-theory roll damping: Clp ≈ −clα/4 evaluated at p·b/(2u).
        let roll_damping = -q_dyn
            * self.wing_area_m2
            * self.wing_span_m
            * self.cl_alpha
            * 0.25
            * (state.p_rad_s * self.wing_span_m / (4.0 * speed_floor));
        let p_dot = (roll_moment + roll_damping) / FLYER_ROLL_INERTIA_KG_M2;
        let p_new = state.p_rad_s + p_dot * dt_s;
        let phi_new = state.phi_rad + p_new * dt_s;
        if !phi_new.is_finite() || !p_new.is_finite() || phi_new.abs() > PHI_CAP_RAD {
            return Err(Refusal {
                code: "lateral-envelope-exceeded",
                message: format!(
                    "roll attitude {} rad outside the admitted ±{PHI_CAP_RAD} band",
                    phi_new
                ),
                ranked_repairs: vec![
                    "reduce the warp command".into(),
                    "close the run at the admitted boundary".into(),
                ],
            });
        }

        // Yaw decomposition (published order fixes the float sum).
        // Induced drag rises quadratically in CL; an antisymmetric twist
        // moves each half's CL by ∓ΔCL, so the first-order DIFFERENTIAL
        // drag scales with CL_trim·ΔCL and its moment about the wing
        // carry-through ahead of the CG OPPOSES the roll command for the
        // decoupled 1901 machine (the historical adverse-yaw law).
        let induced = K_INDUCED_YAW_NM_PER_RAD_AT_REF_Q * loaded_twist_rad * (q_dyn / 76.8);
        let rudder = match self.rudder {
            RudderLinkage::Decoupled => 0.0,
            RudderLinkage::Linked { gain_nm_per_rad } => gain_nm_per_rad * loaded_twist_rad,
        };
        // Reduced profile term: yaw damping only at this tier.
        let profile = -q_dyn
            * self.wing_area_m2
            * self.wing_span_m
            * 0.06
            * (state.r_rad_s * self.wing_span_m / (2.0 * speed_floor));
        let yaw = YawDecompositionNm {
            induced_drag_yaw_nm: induced,
            rudder_yaw_nm: rudder,
            profile_yaw_nm: profile,
        };
        let r_dot = yaw.net() / REDUCED_YAW_INERTIA_KG_M2;
        let r_new = state.r_rad_s + r_dot * dt_s;
        if r_new.abs() > R_CAP_RAD_S {
            return Err(Refusal {
                code: "lateral-envelope-exceeded",
                message: format!("yaw rate {r_new} rad/s outside the admitted ±{R_CAP_RAD_S} band"),
                ranked_repairs: vec![
                    "reduce the warp/rudder command".into(),
                    "close the run at the admitted boundary".into(),
                ],
            });
        }
        let psi_new = state.psi_rad + r_new * dt_s;
        if !psi_new.is_finite() || !r_new.is_finite() {
            return Err(Refusal {
                code: "lateral-envelope-exceeded",
                message: "non-finite yaw state".into(),
                ranked_repairs: vec!["check the published decomposition rows".into()],
            });
        }
        *state = LateralState {
            phi_rad: phi_new,
            p_rad_s: p_new,
            psi_rad: psi_new,
            r_rad_s: r_new,
        };
        Ok(LateralTick {
            state: *state,
            roll_moment_nm: roll_moment,
            yaw,
            loaded_twist_rad,
        })
    }
}
