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
//!   fixed-point on w, warm-started from the previous station, with a
//!   bounded local bracket finish for slowly contracting iterations,
//!   with a per-station convergence receipt — nonconvergence is a TYPED
//!   refusal (plan law), never a clamped answer.
//!
//! Propeller disks are THIS crate's business: rotor spin-up dynamics
//! I_eq·Ω̇ = Q_engine − ΣQ_prop − Q_drivetrain with the declared engine
//! torque curve. fs-wing never sees a disk (plan ownership).

use fs_airfoil::flat_plate_separated;
use fs_math::det;

pub mod dynamic_inflow;
pub mod oblique;

/// Station-count cap (refusals at cap AND cap+1).
pub const MAX_STATIONS: usize = 64;
/// Per-station force-evaluation cap, including any local bracket finish.
pub const MAX_STATION_ITERS: u32 = 120;
/// Relative tolerance on the induced-velocity fixed point.
pub const W_TOL: f64 = 1.0e-10;
/// Maximum admitted braking induction `a = -w/V` for `V > 0`.
///
/// Ordinary momentum theory requires a high-induction correction beyond
/// about 0.4 (AeroDyn Theory Manual, NREL/TP-500-36881, pp. 6–7). This
/// empirical limit is stricter than the real normal-root bound `w >= -V/2`.
pub const MAX_BRAKING_INDUCTION: f64 = 0.4;

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
    fn station_width_m(&self, index: usize) -> f64 {
        let station = &self.stations[index];
        let low = if index == 0 {
            self.stations[0].r_over_r * 0.8
        } else {
            (self.stations[index - 1].r_over_r + station.r_over_r) / 2.0
        };
        let high = if index + 1 == self.stations.len() {
            1.0
        } else {
            (station.r_over_r + self.stations[index + 1].r_over_r) / 2.0
        };
        (high - low) * self.radius_m
    }

    /// Validate the rotor.
    ///
    /// # Errors
    /// `rotor-invalid` (counts/caps at cap AND cap+1, ordering, ranges).
    pub fn admit(&self) -> Result<(), Refusal> {
        if !self.radius_m.is_finite()
            || !(self.radius_m > 0.0)
            || !self.camber_ratio.is_finite()
            || self.n_blades == 0
            || self.n_blades > 8
        {
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
                && s.chord_m.is_finite()
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
    /// Converged induced axial velocity w [m/s]; negative in braking annuli.
    pub w_mps: f64,
    /// Actual station force evaluations used, including local bracket probes.
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
    // Statement-split (guzez.7.2.1): no fused mul-add.
    let a2 = 0.8 * a * a;
    let cd0 = 0.012 + a2; // declared profile-drag model
    let mag = a.abs();
    if mag <= 0.25 {
        (attached_cl, cd0)
    } else {
        // Blend to the separated plate (V-03's honest post-stall shape).
        let sep = flat_plate_separated(a, 6.0).expect("domain-clamped");
        let (cl_s, cd_s) = fs_airfoil::body_to_wind(sep.cn, sep.ca, a);
        let t = ((mag - 0.25) / 0.2).min(1.0);
        let t2 = t * t;
        let s = t2 * (3.0 - 2.0 * t);
        (
            {
                let u = attached_cl * (1.0 - s);
                let v = cl_s * s;
                u + v
            },
            cd0 * (1.0 - s) + cd_s.abs() * s,
        )
    }
}

/// Solve one rotor at (axial speed V, rotation Ω). Warm-started station
/// marching (each station starts from its inboard neighbor's w).
///
/// # Errors
/// Rotor admission; `operating-point-invalid`;
/// `station-did-not-converge` naming the station (typed, never clamped);
/// `station-wake-domain` for high induction or non-normal throughflow;
/// `station-numerical-invalid` for unrepresentable derived arithmetic.
pub fn bemt_solve(
    rotor: &Rotor,
    rho_kg_m3: f64,
    v_axial_mps: f64,
    omega_rad_s: f64,
) -> Result<BemtSolution, Refusal> {
    rotor.admit()?;
    if !rho_kg_m3.is_finite()
        || !(rho_kg_m3 > 0.0)
        || !v_axial_mps.is_finite()
        || v_axial_mps < 0.0
        || !omega_rad_s.is_finite()
        || !(omega_rad_s > 0.0)
    {
        return Err(refuse(
            "operating-point-invalid",
            format!("rho {rho_kg_m3}, V {v_axial_mps}, omega {omega_rad_s}"),
            "V >= 0 (reverse flight is out of domain), omega > 0",
        ));
    }
    let r_tip = rotor.radius_m;
    let r_root = rotor.stations[0].r_over_r * 0.8; // hub cutout fraction
    let mut thrust = 0.0;
    let mut torque = 0.0;
    let mut receipts = Vec::with_capacity(rotor.stations.len());
    let mut w_warm = 0.1 * v_axial_mps.max(1.0);
    for (si, st) in rotor.stations.iter().enumerate() {
        let r = st.r_over_r * r_tip;
        let u_tan = omega_rad_s * r;
        let problem = StationProblem {
            rotor,
            station: st,
            index: si,
            rho: rho_kg_m3,
            axial: v_axial_mps,
            tangential: u_tan,
            radius: r,
            root_radius: r_root * r_tip,
        };
        let (accepted, iters) = problem.solve(w_warm)?;
        let dr = rotor.station_width_m(si);
        let dthrust = accepted.thrust * dr;
        thrust += dthrust;
        let dtorque = accepted.torque * dr;
        torque += dtorque;
        receipts.push(StationReceipt {
            r_over_r: st.r_over_r,
            w_mps: accepted.w,
            iterations: iters,
            alpha_rad: accepted.alpha,
            prandtl_f: accepted.prandtl,
        });
        w_warm = accepted.w;
    }
    let n_rev = omega_rad_s / (2.0 * core::f64::consts::PI);
    let d = 2.0 * r_tip;
    let solution = BemtSolution {
        thrust_n: thrust,
        torque_nm: torque,
        ct: thrust / (rho_kg_m3 * n_rev * n_rev * d.powi(4)),
        cq: torque / (rho_kg_m3 * n_rev * n_rev * d.powi(5)),
        j: v_axial_mps / (n_rev * d),
        stations: receipts,
    };
    if ![
        solution.thrust_n,
        solution.torque_nm,
        solution.ct,
        solution.cq,
        solution.j,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(refuse(
            "operating-point-invalid",
            "unrepresentable rotor loads or coefficients".into(),
            "supply a finite operating point with representable dimensional and nondimensional loads",
        ));
    }
    Ok(solution)
}

/// One original blade-element/momentum equation. The force law and Prandtl
/// factors are identical for Picard, local probes and the accepted receipt.
struct StationProblem<'a> {
    rotor: &'a Rotor,
    station: &'a BladeStation,
    index: usize,
    rho: f64,
    axial: f64,
    tangential: f64,
    radius: f64,
    root_radius: f64,
}

#[derive(Clone, Copy)]
struct StationSample {
    w: f64,
    alpha: f64,
    prandtl: f64,
    thrust: f64,
    torque: f64,
    target_product: f64,
    momentum_product: f64,
    balance: f64,
    defect: f64,
    normal_target: bool,
}

/// Section forces before imposing a steady annular momentum equation.
/// The transient inflow owner uses this identical force law at its stored flow.
struct BladeLoads {
    alpha: f64,
    sin_phi: f64,
    thrust: f64,
    torque: f64,
    profile_drag: f64,
}

impl StationSample {
    fn tolerance(self) -> f64 {
        W_TOL * (self.w.abs() + 1e-9)
    }

    fn admitted(self, axial: f64) -> bool {
        let tolerance = self.tolerance();
        self.normal_target
            && self.w >= -0.5 * axial
            && self.defect.abs() <= tolerance
            && self.balance.abs() <= tolerance * (axial + 2.0 * self.w + tolerance)
    }
}

impl StationProblem<'_> {
    fn refusal(&self, code: &'static str, detail: &str, repair: &str) -> Refusal {
        refuse(
            code,
            format!(
                "station {} (r/R {}): {detail}",
                self.index, self.station.r_over_r
            ),
            repair,
        )
    }

    fn blade_loads(&self, w: f64) -> Result<BladeLoads, Refusal> {
        let u_ax = self.axial + w;
        if !w.is_finite() || !u_ax.is_finite() || !self.tangential.is_finite() {
            return Err(self.refusal(
                "station-numerical-invalid",
                "unrepresentable station velocity",
                "supply a representable finite rotor operating point",
            ));
        }
        if w < -0.5 * self.axial {
            return Err(self.refusal("station-wake-domain", "trial induction reverses the normal far wake",
                "supply a normal axial-throughflow state; reversed-flow wake closures are not implemented"));
        }
        let phi = det::atan2(u_ax, self.tangential);
        let alpha = self.station.beta_rad - phi;
        let (cl, cd) = section_cl_cd(alpha, self.rotor.camber_ratio);
        let ua2 = u_ax * u_ax;
        let ut2 = self.tangential * self.tangential;
        let u2 = ua2 + ut2;
        let q = 0.5 * self.rho * u2;
        let (sin_phi, cos_phi) = (det::sin(phi), det::cos(phi));
        let b = f64::from(self.rotor.n_blades);
        let lift_ax = cl * cos_phi;
        let drag_ax = cd * sin_phi;
        let thrust = q * self.station.chord_m * (lift_ax - drag_ax) * b;
        let lift_tan = cl * sin_phi;
        let drag_tan = cd * cos_phi;
        let torque = q * self.station.chord_m * (lift_tan + drag_tan) * b * self.radius;
        if ![alpha, sin_phi, thrust, torque]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(self.refusal(
                "station-numerical-invalid",
                "nonfinite blade-element arithmetic",
                "supply a representable finite rotor operating point",
            ));
        }
        Ok(BladeLoads {
            alpha,
            sin_phi,
            thrust,
            torque,
            profile_drag: q * self.station.chord_m * cd * b,
        })
    }

    fn evaluate(&self, w: f64) -> Result<StationSample, Refusal> {
        let BladeLoads {
            alpha,
            sin_phi,
            thrust,
            torque,
            ..
        } = self.blade_loads(w)?;
        let b = f64::from(self.rotor.n_blades);
        let f_tip = if sin_phi.abs() > 1e-9 {
            let ft = b / 2.0 * (self.rotor.radius_m - self.radius) / (self.radius * sin_phi.abs());
            2.0 / core::f64::consts::PI * det::acos(det::exp(-ft))
        } else {
            1.0
        };
        let f_root = if sin_phi.abs() > 1e-9 {
            let fr = b / 2.0 * (self.radius - self.root_radius) / (self.radius * sin_phi.abs());
            2.0 / core::f64::consts::PI * det::acos(det::exp(-fr))
        } else {
            1.0
        };
        let prandtl = (f_tip * f_root).clamp(1.0e-3, 1.0);
        let denominator = 4.0 * core::f64::consts::PI * self.radius * self.rho * prandtl;
        let target_product = thrust / denominator;
        let momentum_product = (self.axial + w) * w;
        let balance = target_product - momentum_product;
        // Signed upper quadratic root, continuous through zero loading.
        // A non-real TRIAL target need not mean this section has no root:
        // its loading changes with induction. Such a trial only proposes a
        // half-retreat toward the normal-flow boundary in solve(); it cannot
        // pass convergence or supply local contraction evidence.
        let disc = self.axial * self.axial + 4.0 * target_product;
        let normal_target = disc >= 0.0;
        let target_w = if normal_target {
            let root_denominator = self.axial + det::sqrt(disc);
            if root_denominator == 0.0 {
                0.0
            } else {
                2.0 * target_product / root_denominator
            }
        } else {
            -0.5 * self.axial
        };
        let defect = target_w - w;
        if ![
            w,
            alpha,
            prandtl,
            thrust,
            torque,
            denominator,
            target_product,
            momentum_product,
            balance,
            disc,
            defect,
        ]
        .into_iter()
        .all(f64::is_finite)
            || denominator <= 0.0
        {
            return Err(self.refusal(
                "station-numerical-invalid",
                "nonfinite station arithmetic",
                "supply a representable finite rotor operating point",
            ));
        }
        Ok(StationSample {
            w,
            alpha,
            prandtl,
            thrust,
            torque,
            target_product,
            momentum_product,
            balance,
            defect,
            normal_target,
        })
    }

    fn accept(
        &self,
        sample: StationSample,
        evaluations: u32,
    ) -> Result<(StationSample, u32), Refusal> {
        // A real normal root is not sufficient empirical validity. Do not
        // publish ordinary momentum theory in the high-induction regime.
        if sample.w < -MAX_BRAKING_INDUCTION * self.axial {
            return Err(self.refusal("station-wake-domain", "braking induction exceeds the ordinary-momentum limit a=0.4",
                "supply an applicable high-induction wake closure; the real normal quadratic branch alone is insufficient"));
        }
        Ok((sample, evaluations))
    }

    fn solve(&self, warm: f64) -> Result<(StationSample, u32), Refusal> {
        let mut w = warm;
        let mut previous = None;
        let mut older = None;
        for evaluations in 1..=MAX_STATION_ITERS {
            let sample = self.evaluate(w)?;
            if sample.admitted(self.axial) {
                return self.accept(sample, evaluations);
            }
            if !sample.normal_target && sample.defect.abs() <= sample.tolerance() {
                return Err(self.refusal("station-did-not-converge", "non-real trial target remains unresolved at the normal-flow boundary",
                    "resolve the station continuation or supply an applicable wake closure; this failed iterate does not establish absence of another physical root"));
            }
            // Preserve the ordinary warm-started path. Reserve the final 24
            // evaluations only when recent contraction predicts it will miss
            // the unchanged tolerance, and a small local bracket is justified.
            if evaluations == MAX_STATION_ITERS - 24
                && let (Some(a), Some(b)) = (older, previous)
                && let Some((probe, slope)) = self.local_probe(a, b, sample, 24)
            {
                return self.finish_local(sample, probe, slope, evaluations);
            }
            older = previous;
            previous = Some(sample);
            w += 0.5 * sample.defect;
        }
        Err(self.refusal("station-did-not-converge", "120 station evaluations exhausted",
            "the unchanged station equation remains unresolved; refine the operating-point continuation or section model"))
    }

    fn local_probe(
        &self,
        a: StationSample,
        b: StationSample,
        c: StationSample,
        remaining: i32,
    ) -> Option<(f64, f64)> {
        let first_move = b.w - a.w;
        let last_move = c.w - b.w;
        let first_slope = (b.balance - a.balance) / first_move;
        let last_slope = (c.balance - b.balance) / last_move;
        let ratio = c.defect.abs() / b.defect.abs();
        if !a.normal_target
            || !b.normal_target
            || !c.normal_target
            || ![first_slope, last_slope, ratio]
                .into_iter()
                .all(f64::is_finite)
            || first_slope >= 0.0
            || last_slope >= 0.0
            || !(0.5..=2.0).contains(&(last_slope / first_slope))
            || !(0.0..0.98).contains(&ratio)
            || c.defect.abs() * ratio.powi(remaining) <= c.tolerance()
            || a.balance.signum() != b.balance.signum()
            || b.balance.signum() != c.balance.signum()
            || last_move.signum() != c.balance.signum()
            || first_move.signum() != last_move.signum()
            || !(a.balance.abs() > b.balance.abs() && b.balance.abs() > c.balance.abs())
        {
            return None;
        }
        // Overshoot the local secant root to obtain a sign bracket, never a
        // global interval spanning other roots of the stalled section law.
        let distance = -2.0 * c.balance / last_slope;
        let probe = c.w + distance;
        let local_scale = (c.w.abs() + self.axial + self.tangential).max(1.0);
        if !probe.is_finite()
            || probe < -0.5 * self.axial
            || probe == c.w
            || distance.abs() > 32.0 * last_move.abs()
            || distance.abs() > 1e-3 * local_scale
        {
            return None;
        }
        Some((probe, last_slope))
    }

    fn finish_local(
        &self,
        start: StationSample,
        probe: f64,
        slope: f64,
        mut evaluations: u32,
    ) -> Result<(StationSample, u32), Refusal> {
        let end = self.evaluate(probe)?;
        evaluations += 1;
        let measured_slope = (end.balance - start.balance) / (end.w - start.w);
        if !measured_slope.is_finite()
            || measured_slope >= 0.0
            || !(0.5..=2.0).contains(&(measured_slope / slope))
        {
            return Err(self.refusal(
                "station-did-not-converge",
                "local branch slope changed during bracket probe",
                "resolve the local stalled-section branch; a global root jump is not admitted",
            ));
        }
        if end.admitted(self.axial) {
            return self.accept(end, evaluations);
        }
        if start.balance.signum() == end.balance.signum() {
            return Err(self.refusal(
                "station-did-not-converge",
                "local continuation did not bracket its root",
                "resolve the warm-started station branch without selecting a distant root",
            ));
        }
        let (mut left, mut right) = if start.w < end.w {
            (start, end)
        } else {
            (end, start)
        };
        while evaluations < MAX_STATION_ITERS {
            let sample = self.evaluate(left.w + 0.5 * (right.w - left.w))?;
            evaluations += 1;
            let roundoff = 64.0
                * f64::EPSILON
                * sample
                    .target_product
                    .abs()
                    .max(sample.momentum_product.abs());
            if sample.balance > left.balance + roundoff || sample.balance < right.balance - roundoff
            {
                return Err(self.refusal("station-did-not-converge", "nonmonotone local station bracket",
                    "resolve the local stalled-section branch; multiple-root selection is not inferred"));
            }
            if sample.admitted(self.axial) {
                return self.accept(sample, evaluations);
            }
            if sample.balance > 0.0 {
                left = sample;
            } else {
                right = sample;
            }
        }
        Err(self.refusal(
            "station-did-not-converge",
            "120 station evaluations exhausted in the local bracket",
            "the local station balance remains unresolved at the declared tolerance",
        ))
    }
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
