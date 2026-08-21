//! fs-airscrew — propulsion (L3). Bead frankensim-wf-root-guzez.5.12.1
//! (E4.5-i, Wright Flyer program). Spec: COMPREHENSIVE_PLAN §5.3.
//!
//! BEMT in the w-FORMULATION (induced axial velocity w, not the a-factor),
//! which stays valid through J = 0 — the Dec-17 held-on-rail state starts
//! at J ≈ 0.7 but static bench anchors (E1.6) sit at J = 0. Per station:
//!
//!   U_ax = V + w,  U_tan = Ω·r (swirl folded into the section closure's
//!   effective angle via the small-swirl approximation — DISCLOSED),
//!   φ = atan2(U_ax, U_tan), α = β − φ, (cl, cd) from the section closure,
//!   blade-element dT vs momentum dT = 4π·r·ρ·U_ax·w·F  (Prandtl F),
//!   fixed-point on w, warm-started from the previous station, bounded,
//!   with a per-station convergence receipt — nonconvergence is a TYPED
//!   refusal (plan law), never a clamped answer.
//!
//! Propeller disks are THIS crate's business: rotor spin-up dynamics
//! I_eq·Ω̇ = Q_engine − ΣQ_prop − Q_drivetrain with the declared engine
//! torque curve. fs-wing never sees a disk (plan ownership).

use fs_airfoil::flat_plate_separated;
use fs_math::det;

/// Station-count cap (refusals at cap AND cap+1).
pub const MAX_STATIONS: usize = 64;
/// Per-station iteration cap.
pub const MAX_STATION_ITERS: u32 = 120;
/// Relative tolerance on the induced-velocity fixed point.
pub const W_TOL: f64 = 1.0e-10;

/// A typed refusal (workspace law).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs.
    pub ranked_repairs: Vec<String>,
}

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// One blade station (from prop-geometry-v1: 1911 calibration table or a
/// declared 1903 reconstruction — the PROVENANCE lives with the caller).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BladeStation {
    /// Radial fraction r/R (0, 1).
    pub r_over_r: f64,
    /// Chord [m].
    pub chord_m: f64,
    /// Geometric pitch angle β from the rotation plane [rad].
    pub beta_rad: f64,
}

/// The rotor.
#[derive(Clone, Debug, PartialEq)]
pub struct Rotor {
    /// Tip radius [m].
    pub radius_m: f64,
    /// Blade count.
    pub n_blades: u32,
    /// Section camber ratio for the closure (Estimated reconstruction).
    pub camber_ratio: f64,
    /// Stations, strictly ascending in r/R.
    pub stations: Vec<BladeStation>,
}

impl Rotor {
    /// Validate the rotor.
    ///
    /// # Errors
    /// `rotor-invalid` (counts/caps at cap AND cap+1, ordering, ranges).
    pub fn admit(&self) -> Result<(), Refusal> {
        if !(self.radius_m > 0.0) || self.n_blades == 0 || self.n_blades > 8 {
            return Err(refuse(
                "rotor-invalid",
                "radius/blades".into(),
                "R>0, 1-8 blades",
            ));
        }
        if self.stations.len() < 3 || self.stations.len() > MAX_STATIONS {
            return Err(refuse(
                "rotor-invalid",
                format!(
                    "{} stations outside [3, {MAX_STATIONS}]",
                    self.stations.len()
                ),
                "the 1911 table has 8-10",
            ));
        }
        let mut prev = 0.0;
        for s in &self.stations {
            let ok = s.r_over_r > prev
                && s.r_over_r < 1.0
                && s.chord_m > 0.0
                && s.beta_rad.is_finite()
                && s.beta_rad > 0.0
                && s.beta_rad < 1.5;
            if !ok {
                return Err(refuse(
                    "rotor-invalid",
                    format!("station at r/R {} invalid", s.r_over_r),
                    "ascending r/R in (0,1); positive chord; beta in (0, 1.5) rad",
                ));
            }
            prev = s.r_over_r;
        }
        Ok(())
    }
}

/// Per-station convergence receipt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StationReceipt {
    /// r/R.
    pub r_over_r: f64,
    /// Converged induced axial velocity w [m/s].
    pub w_mps: f64,
    /// Iterations used.
    pub iterations: u32,
    /// Local effective angle of attack [rad].
    pub alpha_rad: f64,
    /// Prandtl combined tip/root factor at this station.
    pub prandtl_f: f64,
}

/// The BEMT solution.
#[derive(Clone, Debug, PartialEq)]
pub struct BemtSolution {
    /// Thrust [N] (one rotor).
    pub thrust_n: f64,
    /// Torque [N·m].
    pub torque_nm: f64,
    /// CT = T/(ρ n² D⁴), n in rev/s (prop-geometry-v1 convention).
    pub ct: f64,
    /// CQ = Q/(ρ n² D⁵).
    pub cq: f64,
    /// Advance ratio J = V/(nD).
    pub j: f64,
    /// Per-station receipts.
    pub stations: Vec<StationReceipt>,
}

/// Section closure: attached thin-airfoil + flat-plate blend (the same
/// documented shape as the E4.1 datasets; prop tables refine later).
fn section_cl_cd(alpha: f64, camber: f64) -> (f64, f64) {
    let a = alpha.clamp(-1.3, 1.3);
    let attached_cl = 2.0 * core::f64::consts::PI * (a + 2.0 * camber);
    let cd0 = 0.012 + 0.8 * a * a; // declared profile-drag model
    let mag = a.abs();
    if mag <= 0.25 {
        (attached_cl, cd0)
    } else {
        // Blend to the separated plate (V-03's honest post-stall shape).
        let sep = flat_plate_separated(a, 6.0).expect("domain-clamped");
        let (cl_s, cd_s) = fs_airfoil::body_to_wind(sep.cn, sep.ca, a);
        let t = ((mag - 0.25) / 0.2).min(1.0);
        let s = t * t * (3.0 - 2.0 * t);
        (
            attached_cl * (1.0 - s) + cl_s * s,
            cd0 * (1.0 - s) + cd_s.abs() * s,
        )
    }
}

/// Solve one rotor at (axial speed V, rotation Ω). Warm-started station
/// marching (each station starts from its inboard neighbor's w).
///
/// # Errors
/// Rotor admission; `operating-point-invalid`;
/// `station-did-not-converge` naming the station (typed, never clamped).
pub fn bemt_solve(
    rotor: &Rotor,
    rho_kg_m3: f64,
    v_axial_mps: f64,
    omega_rad_s: f64,
) -> Result<BemtSolution, Refusal> {
    rotor.admit()?;
    if !(rho_kg_m3 > 0.0) || !v_axial_mps.is_finite() || v_axial_mps < 0.0 || !(omega_rad_s > 0.0) {
        return Err(refuse(
            "operating-point-invalid",
            format!("rho {rho_kg_m3}, V {v_axial_mps}, omega {omega_rad_s}"),
            "V >= 0 (reverse flight is out of domain), omega > 0",
        ));
    }
    let r_tip = rotor.radius_m;
    let b = f64::from(rotor.n_blades);
    let r_root = rotor.stations[0].r_over_r * 0.8; // hub cutout fraction
    let mut thrust = 0.0;
    let mut torque = 0.0;
    let mut receipts = Vec::with_capacity(rotor.stations.len());
    let mut w_warm = 0.1 * v_axial_mps.max(1.0);
    for (si, st) in rotor.stations.iter().enumerate() {
        let r = st.r_over_r * r_tip;
        let u_tan = omega_rad_s * r;
        let mut w = w_warm;
        let mut converged = false;
        let mut iters = 0;
        let (mut alpha, mut f_prandtl) = (0.0, 1.0);
        let (mut dt_be, mut dq_be) = (0.0, 0.0);
        while iters < MAX_STATION_ITERS {
            iters += 1;
            let u_ax = v_axial_mps + w;
            let phi = det::atan2(u_ax, u_tan);
            alpha = st.beta_rad - phi;
            let (cl, cd) = section_cl_cd(alpha, rotor.camber_ratio);
            let u2 = u_ax * u_ax + u_tan * u_tan;
            let q = 0.5 * rho_kg_m3 * u2;
            let (sin_phi, cos_phi) = (det::sin(phi), det::cos(phi));
            dt_be = q * st.chord_m * (cl * cos_phi - cd * sin_phi) * b;
            dq_be = q * st.chord_m * (cl * sin_phi + cd * cos_phi) * b * r;
            // Prandtl tip + root factors.
            // Prandtl: F = (2/pi)·acos(e^(−f)).
            let f_tip = if sin_phi.abs() > 1e-9 {
                let ft = b / 2.0 * (r_tip - r) / (r * sin_phi.abs());
                2.0 / core::f64::consts::PI * det::exp(-ft).acos()
            } else {
                1.0
            };
            let f_root = if sin_phi.abs() > 1e-9 {
                let fr = b / 2.0 * (r - r_root * r_tip) / (r * sin_phi.abs());
                2.0 / core::f64::consts::PI * det::exp(-fr).acos()
            } else {
                1.0
            };
            f_prandtl = (f_tip * f_root).clamp(1.0e-3, 1.0);
            // Momentum balance for THIS annulus width is applied by the
            // caller-side integration; the fixed point equates the
            // per-length loadings: dT_be = 4π r ρ (V + w) w F.
            let denom = 4.0 * core::f64::consts::PI * r * rho_kg_m3 * f_prandtl;
            let target_prod = dt_be / denom; // (V + w)·w target
            // Solve (V + w_new)·w_new = target_prod for w_new >= 0.
            let disc = v_axial_mps * v_axial_mps + 4.0 * target_prod.max(0.0);
            let w_new = 0.5 * (-v_axial_mps + det::sqrt(disc));
            let step = w_new - w;
            w += 0.5 * step;
            if step.abs() <= W_TOL * (w.abs() + 1e-9) {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(refuse(
                "station-did-not-converge",
                format!(
                    "station {si} (r/R {}) after {MAX_STATION_ITERS} iterations",
                    st.r_over_r
                ),
                "high-loading/static closure limits; the operating point is outside the admitted domain",
            ));
        }
        // Annulus width from station midpoints.
        let r_lo = if si == 0 {
            rotor.stations[0].r_over_r * 0.8
        } else {
            (rotor.stations[si - 1].r_over_r + st.r_over_r) / 2.0
        };
        let r_hi = if si + 1 == rotor.stations.len() {
            1.0
        } else {
            (st.r_over_r + rotor.stations[si + 1].r_over_r) / 2.0
        };
        let dr = (r_hi - r_lo) * r_tip;
        thrust += dt_be * dr;
        torque += dq_be * dr;
        receipts.push(StationReceipt {
            r_over_r: st.r_over_r,
            w_mps: w,
            iterations: iters,
            alpha_rad: alpha,
            prandtl_f: f_prandtl,
        });
        w_warm = w;
    }
    let n_rev = omega_rad_s / (2.0 * core::f64::consts::PI);
    let d = 2.0 * r_tip;
    Ok(BemtSolution {
        thrust_n: thrust,
        torque_nm: torque,
        ct: thrust / (rho_kg_m3 * n_rev * n_rev * d.powi(4)),
        cq: torque / (rho_kg_m3 * n_rev * n_rev * d.powi(5)),
        j: v_axial_mps / (n_rev * d),
        stations: receipts,
    })
}

/// The declared engine torque curve at the PROP shaft (12 hp sustained at
/// 1025 engine rpm through the 23:8 chain; flat-torque approximation
/// DISCLOSED — the real curve is an E1 follow-on).
#[must_use]
pub fn engine_torque_at_prop_nm(omega_prop_rad_s: f64) -> f64 {
    let p_w = 8948.0; // 12 hp sustained (flyer-reference)
    let omega_rated_prop = 1025.0 / 60.0 * 2.0 * core::f64::consts::PI * (8.0 / 23.0);
    // Flat torque up to rated speed, power-limited beyond.
    let q_rated = p_w / omega_rated_prop / 2.0; // per prop (two props share)
    if omega_prop_rad_s <= omega_rated_prop {
        q_rated
    } else {
        p_w / omega_prop_rad_s / 2.0
    }
}

/// One rotor spin-up step: I_eq·Ω̇ = Q_avail − Q_prop (per prop).
///
/// # Errors
/// `rotor-dynamics-invalid`.
pub fn rotor_spinup_step(
    i_eq_kgm2: f64,
    omega_rad_s: f64,
    q_engine_nm: f64,
    q_prop_nm: f64,
    dt_s: f64,
) -> Result<f64, Refusal> {
    if !(i_eq_kgm2 > 0.0) || !omega_rad_s.is_finite() || !dt_s.is_finite() || dt_s <= 0.0 {
        return Err(refuse(
            "rotor-dynamics-invalid",
            "I/omega/dt".into(),
            "positive I and dt",
        ));
    }
    Ok(omega_rad_s + (q_engine_nm - q_prop_nm) / i_eq_kgm2 * dt_s)
}
