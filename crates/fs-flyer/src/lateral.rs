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
//!   strip-theory damping term ∝ q·S·b²·clα·p/(4u); the held-coefficient
//!   linear rate equation and φ̇ = p are integrated analytically together.
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
use fs_math::{decay, det};

mod gust;
pub use gust::{DistributedLateralLoads, LateralMomentsNm};

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
    /// Step-mean profile/other yaw moment (reduced damping term).
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
    /// Held differential-lift roll moment this tick [N·m], excluding damping.
    pub roll_moment_nm: f64,
    /// Step-mean yaw-moment decomposition applied this tick [N·m].
    pub yaw: YawDecompositionNm,
    /// Total LOADED twist consumed this tick [rad] (signed).
    pub loaded_twist_rad: f64,
}

/// One accepted lateral step with independently resolved physical moments.
/// The original reduced attribution stays separate from the incident loads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForcedLateralTick {
    /// Original control/damping attribution and the resulting accepted state.
    pub reduced: LateralTick,
    /// Additional body-FRD moments, applied once in the original rate update.
    pub incident: LateralMomentsNm,
    /// Step-mean roll moment including control, damping and incident load [N m].
    pub net_roll_moment_nm: f64,
    /// Step-mean yaw moment including control, damping and incident load [N m].
    pub net_yaw_moment_nm: f64,
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

    /// One lateral step, normally called at 120 Hz.
    ///
    /// `loaded_twist_rad`: signed total loaded twist from the warp model
    /// (positive = right wing down). `u_mps`: nonnegative body airspeed.
    /// Density and duration are positive. Aerodynamic coefficients and external
    /// moments are held over the interval; rate-dependent damping evolves with
    /// the actual rate. Both angular rate and its integral are analytic for
    /// that local model, including the zero-damping limit.
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
        self.step_with_moments(
            state,
            loaded_twist_rad,
            u_mps,
            rho_kg_m3,
            dt_s,
            LateralMomentsNm::default(),
        )
        .map(|step| step.reduced)
    }

    /// Advance the same reduced lateral state with resolved incident torques.
    ///
    /// `incident` is the sum of actual distributed panel, profile-drag and
    /// propeller moments, not a second warp or damping law. Both controls and
    /// incident moments enter the same held-coefficient rate/attitude solution.
    /// Returned damping and total moments are interval means, so their impulses
    /// match the accepted rate change to roundoff. Zero moments preserve the
    /// unforced entry's arithmetic. The additional load remains separate from
    /// the original yaw attribution.
    /// This remains a reduced lateral response, not a banked six-DoF path.
    ///
    /// # Errors
    /// As [`Self::step`], also on nonfinite incident moments. Every refusal
    /// leaves the supplied lateral state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_moments(
        &self,
        state: &mut LateralState,
        loaded_twist_rad: f64,
        u_mps: f64,
        rho_kg_m3: f64,
        dt_s: f64,
        incident: LateralMomentsNm,
    ) -> Result<ForcedLateralTick, Refusal> {
        if !incident.is_finite() {
            return Err(Refusal {
                code: "lateral-envelope-exceeded",
                message: "non-finite incident lateral moment".into(),
                ranked_repairs: vec!["retain the finite accepted distributed loads".into()],
            });
        }
        if ![
            loaded_twist_rad,
            u_mps,
            rho_kg_m3,
            dt_s,
            self.wing_area_m2,
            self.wing_span_m,
            self.cl_alpha,
            self.arm_factor,
            state.phi_rad,
            state.p_rad_s,
            state.psi_rad,
            state.r_rad_s,
        ]
        .into_iter()
        .all(f64::is_finite)
            || u_mps < 0.0
            || rho_kg_m3 <= 0.0
            || dt_s <= 0.0
            || self.wing_area_m2 <= 0.0
            || self.wing_span_m <= 0.0
            || self.cl_alpha < 0.0
            || !(0.0..=0.5).contains(&self.arm_factor)
            || state.phi_rad.abs() > PHI_CAP_RAD
            || state.r_rad_s.abs() > R_CAP_RAD_S
            || matches!(self.rudder, RudderLinkage::Linked { gain_nm_per_rad }
                if !gain_nm_per_rad.is_finite())
        {
            return Err(invalid(
                "finite in-envelope state, nonnegative speed/slope, positive density/duration/geometry and a physical arm factor required",
            ));
        }
        let q_dyn = 0.5 * rho_kg_m3 * u_mps * u_mps;
        // q/u = rho*u/2 has a regular zero-speed limit. An arbitrary speed
        // floor changes the declared linear damping law near rest.
        let q_over_speed = 0.5 * rho_kg_m3 * u_mps;
        // Differential-lift roll moment: the signed total twist already
        // carries the antisymmetry of the strip distribution.
        let roll_moment = q_dyn
            * self.wing_area_m2
            * self.wing_span_m
            * self.cl_alpha
            * loaded_twist_rad
            * self.arm_factor
            * 0.5;
        // Preserve the declared reduced damping: clα/4 at p·b/(4u).
        let roll_damping_coefficient = q_over_speed
            * self.wing_area_m2
            * self.wing_span_m
            * self.cl_alpha
            * 0.25
            * (self.wing_span_m / 4.0);
        let roll_forcing = roll_moment + incident.roll_nm;
        let roll = held_rate_step(
            state.p_rad_s,
            roll_forcing,
            roll_damping_coefficient,
            FLYER_ROLL_INERTIA_KG_M2,
            dt_s,
        )?;
        let p_new = roll.endpoint;
        let phi_new = state.phi_rad + roll.mean * dt_s;
        if !phi_new.is_finite() || phi_new.abs() > PHI_CAP_RAD {
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
        // A forced reversal has one interior angular extremum. Checking only
        // the endpoint would admit a coarse step that exits and re-enters the
        // reduced model's roll envelope.
        if state.p_rad_s != 0.0
            && p_new != 0.0
            && state.p_rad_s.is_sign_negative() != p_new.is_sign_negative()
        {
            let initial_time = -state.p_rad_s / (roll_forcing / FLYER_ROLL_INERTIA_KG_M2);
            let z = roll_damping_coefficient / FLYER_ROLL_INERTIA_KG_M2 * initial_time;
            let one_plus_z = 1.0 + z;
            let turn_time = initial_time
                * if one_plus_z == 1.0 {
                    1.0
                } else {
                    det::ln(one_plus_z) / (one_plus_z - 1.0)
                };
            if !turn_time.is_finite() || turn_time <= 0.0 || turn_time >= dt_s {
                return Err(invalid("unrepresentable lateral reversal time"));
            }
            let turn = held_rate_step(
                state.p_rad_s,
                roll_forcing,
                roll_damping_coefficient,
                FLYER_ROLL_INERTIA_KG_M2,
                turn_time,
            )?;
            if (state.phi_rad + turn.mean * turn_time).abs() > PHI_CAP_RAD {
                return Err(invalid("roll reversal crosses the admitted attitude band"));
            }
        }
        let net_roll_moment_nm = roll_forcing - roll_damping_coefficient * roll.mean;

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
        let yaw_damping_coefficient =
            q_over_speed * self.wing_area_m2 * self.wing_span_m * 0.06 * (self.wing_span_m / 2.0);
        let yaw_step = held_rate_step(
            state.r_rad_s,
            (induced + rudder) + incident.yaw_nm,
            yaw_damping_coefficient,
            REDUCED_YAW_INERTIA_KG_M2,
            dt_s,
        )?;
        let yaw = YawDecompositionNm {
            induced_drag_yaw_nm: induced,
            rudder_yaw_nm: rudder,
            profile_yaw_nm: -yaw_damping_coefficient * yaw_step.mean,
        };
        let mut net_yaw_moment_nm = yaw.net();
        if incident.yaw_nm != 0.0 {
            net_yaw_moment_nm += incident.yaw_nm;
        }
        let r_new = yaw_step.endpoint;
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
        let psi_new = state.psi_rad + yaw_step.mean * dt_s;
        if ![psi_new, net_roll_moment_nm, net_yaw_moment_nm]
            .into_iter()
            .all(f64::is_finite)
        {
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
        Ok(ForcedLateralTick {
            reduced: LateralTick {
                state: *state,
                roll_moment_nm: roll_moment,
                yaw,
                loaded_twist_rad,
            },
            incident,
            net_roll_moment_nm,
            net_yaw_moment_nm,
        })
    }
}

fn invalid(message: &str) -> Refusal {
    Refusal {
        code: "lateral-envelope-exceeded",
        message: message.into(),
        ranked_repairs: vec![
            "use finite physical inputs inside the reduced lateral envelope".into(),
        ],
    }
}

struct RateStep {
    endpoint: f64,
    mean: f64,
}

// I w' = torque - damping*w, with all coefficients held. For x=damping*dt/I,
// w1 = exp(-x)*w0 + (torque/I)*dt*phi1(-x), and the mean is
// phi1(-x)*w0 + (torque/I)*dt*phi2(-x). Integrating the mean advances angle
// without the spurious extra half-impulse of an endpoint-Euler drift.
fn held_rate_step(
    initial: f64,
    torque: f64,
    damping: f64,
    inertia: f64,
    dt: f64,
) -> Result<RateStep, Refusal> {
    let rate = damping / inertia;
    let exposure = rate * dt;
    if ![torque, damping, exposure].into_iter().all(f64::is_finite) || damping < 0.0 {
        return Err(invalid("unrepresentable lateral torque or damping"));
    }
    let phi1 = if exposure == 0.0 {
        1.0
    } else {
        decay::decayed_fraction(exposure) / exposure
    };
    let phi2 = if exposure < 1e-3 {
        // Avoid subtracting nearly equal numbers as damping tends to zero.
        0.5 + exposure
            * (-1.0 / 6.0 + exposure * (1.0 / 24.0 + exposure * (-1.0 / 120.0 + exposure / 720.0)))
    } else {
        (1.0 - phi1) / exposure
    };
    let acceleration = torque / inertia;
    let endpoint = decay::attenuate(initial, exposure) + acceleration * (dt * phi1);
    let mean = initial * phi1 + acceleration * (dt * phi2);
    if !endpoint.is_finite() || !mean.is_finite() {
        return Err(invalid("unrepresentable lateral angular motion"));
    }
    Ok(RateStep { endpoint, mean })
}
