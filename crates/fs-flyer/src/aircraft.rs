//! Open-loop integrated aircraft: full force build-up + fixed-control
//! trim continuation (bead wf-root-guzez.5.13.2, E4.6a-ii). Plan §5.1:
//! at a FROZEN control state (canard setting; warp/rudder zero; THROTTLE
//! fixed — the prop speed is a torque-balance STATE, not a control) the
//! build-up assembles, in the FRD body frame
//! (frame-conventions-v1: +x nose, +z down, pitch +M nose-up):
//!
//!   - the coupled multisurface solve (biplane wings + biplane canard in
//!     ONE influence system — canard–wing interference is solved, never
//!     superposed) under the converged prop slipstream (E4.5),
//!   - per-panel Kutta–Joukowsky forces F = ρΓ(V×seg) (freestream-KJ
//!     tier: induced drag is NOT in these forces and is carried as the
//!     declared classical term CL²/(π·AR·e) — documented, Estimated),
//!   - per-strip section quarter-chord moments (thin-airfoil cm₀),
//!   - wing profile drag (declared cd₀), the parasite LEDGER (E4.6a-i),
//!   - thrust at the disk lines, gravity at the CG reference.
//!
//! Trim: damped Newton on r(V, α, δc, ω) = (ΣFx, ΣFz, ΣMy,
//! Q_engine − Q_prop) with θ = α (level flight), deterministic FD
//! Jacobian, typed `trim-not-found` refusal carrying the residual
//! trajectory — never a silent best-effort state.
//!
//! Ownership (AeroEffectOwners): 3-D induction lives in the panel solve;
//! separation in the section closure blend; noncirculatory terms are NOT
//! assembled here (AddedMassOnly — they vanish at a steady trim point
//! anyway); far wake = trailing legs at this tier.

use crate::Refusal;
use crate::dragledger::{DragLedger, wright_ledger_v1};
use crate::propcoupling::{
    CANDIDATE_A, CoupledStep, PropCouplingSolverSpec, PropDisk, coupled_prop_airframe_step,
};
use fs_airscrew::{BladeStation, Rotor};
use fs_blake3::hash_domain;
use fs_math::det;
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec};
use fs_wing::{Panel, SurfaceId, flat_surface};

/// Trim iteration cap.
pub const MAX_TRIM_ITERATIONS: u32 = 40;

/// The fixed-control open-loop design (verified geometry from
/// flyer-reference; Estimated placements carry provenance notes).
#[derive(Clone, Debug, PartialEq)]
pub struct OpenLoopDesign {
    /// Wing span [m] (verified 12.29).
    pub wing_span_m: f64,
    /// Wing chord [m] (verified 1.981).
    pub wing_chord_m: f64,
    /// Wing vertical gap [m] (verified 1.89).
    pub wing_gap_m: f64,
    /// Canard span per plane [m].
    pub canard_span_m: f64,
    /// Canard chord [m].
    pub canard_chord_m: f64,
    /// Canard leading edge x (FRD, forward of the wing LE at 0).
    pub canard_x_le_m: f64,
    /// Canard plane heights (z, FRD down-negative-is-up).
    pub canard_z_m: [f64; 2],
    /// Wing/canard section camber ratio (Estimated 1/20 class).
    pub camber_ratio: f64,
    /// Wing profile drag coefficient (declared Estimated).
    pub cd0_wing: f64,
    /// Oswald-class span-efficiency for the induced term (Estimated).
    pub oswald_e: f64,
    /// CG reference (FRD, from the wing LE origin). Estimated ~30% chord
    /// aft, between the planes (Culick/Jex discussion class).
    pub cg_m: [f64; 3],
    /// Gross mass [kg] (verified 340.2).
    pub gross_mass_kg: f64,
    /// Propeller disks (aft = −x).
    pub disks: [PropDisk; 2],
    /// Rotor definition (E1.6 reconstruction path).
    pub rotor: Rotor,
    /// Parasite ledger.
    pub ledger: DragLedger,
    /// z at which the ledger's parasite drag acts (Estimated near CG).
    pub drag_z_m: f64,
    /// Coupling solver spec.
    pub coupling: PropCouplingSolverSpec,
    /// Spanwise panels per wing plane.
    pub n_span_wing: usize,
    /// Chordwise rows per wing plane (>= 2, plan Tier-A law).
    pub n_chord_wing: usize,
}

/// The registered 1903 rotor stations (E1.6 reconstruction path #2 —
/// same table the fs-airscrew battery pinned).
#[must_use]
pub fn wright_rotor_v1() -> Rotor {
    let deg = core::f64::consts::PI / 180.0;
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: vec![
            BladeStation {
                r_over_r: 0.30,
                chord_m: 0.13,
                beta_rad: 40.0 * deg,
            },
            BladeStation {
                r_over_r: 0.45,
                chord_m: 0.17,
                beta_rad: 30.0 * deg,
            },
            BladeStation {
                r_over_r: 0.60,
                chord_m: 0.20,
                beta_rad: 23.0 * deg,
            },
            BladeStation {
                r_over_r: 0.75,
                chord_m: 0.21,
                beta_rad: 18.5 * deg,
            },
            BladeStation {
                r_over_r: 0.88,
                chord_m: 0.20,
                beta_rad: 15.5 * deg,
            },
            BladeStation {
                r_over_r: 0.96,
                chord_m: 0.16,
                beta_rad: 14.0 * deg,
            },
        ],
    }
}

/// The registered v1 open-loop design.
#[must_use]
pub fn wright_openloop_v1() -> OpenLoopDesign {
    let rpm_prop = 1025.0 * 8.0 / 23.0; // engine 1025 rpm through 23:8 chain
    let omega = rpm_prop / 60.0 * core::f64::consts::TAU;
    OpenLoopDesign {
        wing_span_m: 12.29,
        wing_chord_m: 1.981,
        wing_gap_m: 1.89,
        canard_span_m: 3.66,
        canard_chord_m: 0.61,
        canard_x_le_m: 2.9, // canard qc ~2.23 m ahead of the wing qc
        // Canard planes sit BELOW the lower wing (mounted off the skids;
        // period photos/drawings) — FRD +z down. This is also numerically
        // load-bearing: a canard near-coplanar with a wing plane runs its
        // discrete trailing legs within ~0.1 m of the wing's bound line
        // and the multisurface Picard iteration diverges (measured).
        canard_z_m: [1.05, 0.35],
        camber_ratio: 0.05,
        cd0_wing: 0.015,
        oswald_e: 0.7,
        cg_m: [-0.594, 0.0, -0.6], // 30% chord aft of the wing LE, mid-bay
        gross_mass_kg: 340.2,
        disks: [
            PropDisk {
                center_m: [-2.5, -1.7, -1.0],
                omega_rad_s: omega,
            },
            PropDisk {
                center_m: [-2.5, 1.7, -1.0],
                omega_rad_s: omega,
            },
        ],
        rotor: wright_rotor_v1(),
        ledger: wright_ledger_v1(),
        drag_z_m: -0.6,
        coupling: CANDIDATE_A,
        n_span_wing: 8,
        n_chord_wing: 2,
    }
}

/// One force/moment build-up at a frozen state.
#[derive(Clone, Debug, PartialEq)]
pub struct ForceBuildup {
    /// Net body-frame force [N] (FRD; includes gravity).
    pub force_n: [f64; 3],
    /// Net pitch moment about the CG reference [N·m] (+nose-up).
    pub moment_y_nm: f64,
    /// Wing+canard lift (up-positive) [N].
    pub lift_n: f64,
    /// Total drag (parasite + profile + induced) [N].
    pub drag_n: f64,
    /// L/R thrust [N].
    pub thrust_n: [f64; 2],
    /// Induced-drag term [N] (declared classical formula).
    pub induced_drag_n: f64,
    /// The coupled prop step receipt.
    pub coupled: CoupledStep,
    /// Engine-vs-prop torque imbalance [N·m] (per prop, symmetric):
    /// Q_engine(omega) − mean(Q_prop). Zero at a true fixed-throttle trim.
    pub torque_imbalance_nm: f64,
}

/// Trim result + receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct TrimResult {
    /// Trim airspeed [m/s].
    pub v_mps: f64,
    /// Trim angle of attack (= pitch attitude, level flight) [rad].
    pub alpha_rad: f64,
    /// Trim canard setting [rad].
    pub delta_canard_rad: f64,
    /// Trim prop speed [rad/s] (torque-balance state, not a control).
    pub omega_prop_rad_s: f64,
    /// Final residuals (Fx, Fz, My, Q_imbalance).
    pub residuals: [f64; 4],
    /// Newton iterations used.
    pub iterations: u32,
    /// The build-up at trim.
    pub buildup: ForceBuildup,
    /// Design digest the trim ran under.
    pub design_digest: String,
}

impl OpenLoopDesign {
    /// Design content digest (ModelId ingredient).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut p = Vec::new();
        for v in [
            self.wing_span_m,
            self.wing_chord_m,
            self.wing_gap_m,
            self.canard_span_m,
            self.canard_chord_m,
            self.canard_x_le_m,
            self.canard_z_m[0],
            self.canard_z_m[1],
            self.camber_ratio,
            self.cd0_wing,
            self.oswald_e,
            self.cg_m[0],
            self.cg_m[1],
            self.cg_m[2],
            self.gross_mass_kg,
            self.drag_z_m,
        ] {
            p.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        p.extend_from_slice(self.ledger.digest().as_bytes());
        p.extend_from_slice(self.coupling.digest().as_bytes());
        hash_domain("org.frankensim.fs-flyer.openloop-design.v1", &p).to_hex()
    }

    /// Build panels + strips for the frozen control state. The canard
    /// deflection enters as strip twist (positive δc = positive canard
    /// lift = positive pitch moment, control-signs-v1).
    ///
    /// # Errors
    /// Layout refusals from `flat_surface` pass through.
    pub fn layout(&self, delta_canard_rad: f64) -> Result<(Vec<Panel>, Vec<StripSpec>), Refusal> {
        let map_err = |e: fs_wing::Refusal| Refusal {
            code: e.code,
            message: e.message,
            ranked_repairs: e.ranked_repairs,
        };
        let (ns, nc) = (self.n_span_wing, self.n_chord_wing);
        let mut panels = flat_surface(
            SurfaceId::WingLower,
            self.wing_span_m,
            self.wing_chord_m,
            0.0,
            0.0,
            ns,
            nc,
        )
        .map_err(map_err)?;
        panels.extend(
            flat_surface(
                SurfaceId::WingUpper,
                self.wing_span_m,
                self.wing_chord_m,
                0.0,
                -self.wing_gap_m,
                ns,
                nc,
            )
            .map_err(map_err)?,
        );
        let base_canard = panels.len();
        let ncs = 4; // canard spanwise panels per plane
        panels.extend(
            flat_surface(
                SurfaceId::CanardLower,
                self.canard_span_m,
                self.canard_chord_m,
                self.canard_x_le_m,
                self.canard_z_m[0],
                ncs,
                1,
            )
            .map_err(map_err)?,
        );
        panels.extend(
            flat_surface(
                SurfaceId::CanardUpper,
                self.canard_span_m,
                self.canard_chord_m,
                self.canard_x_le_m,
                self.canard_z_m[1],
                ncs,
                1,
            )
            .map_err(map_err)?,
        );
        let mut strips = Vec::new();
        for plane in 0..2 {
            let base = plane * ns * nc;
            for s in 0..ns {
                strips.push(StripSpec {
                    panel_indices: (0..nc).map(|c| base + c * ns + s).collect(),
                    chord_m: self.wing_chord_m,
                    twist_rad: 0.0,
                });
            }
        }
        for plane in 0..2 {
            let base = base_canard + plane * ncs;
            for s in 0..ncs {
                strips.push(StripSpec {
                    panel_indices: vec![base + s],
                    chord_m: self.canard_chord_m,
                    twist_rad: delta_canard_rad,
                });
            }
        }
        Ok((panels, strips))
    }

    /// Full force/moment build-up at (V, α, δc), θ = α (level flight).
    ///
    /// # Errors
    /// Wing/coupling/ledger refusals pass through; `state-invalid`.
    pub fn force_buildup(
        &self,
        v_mps: f64,
        alpha_rad: f64,
        delta_canard_rad: f64,
        omega_prop_rad_s: f64,
        rho_kg_m3: f64,
    ) -> Result<ForceBuildup, Refusal> {
        if !(v_mps.is_finite()
            && v_mps > 1.0
            && alpha_rad.is_finite()
            && alpha_rad.abs() < 0.6
            && delta_canard_rad.is_finite()
            && delta_canard_rad.abs() < 0.6
            && omega_prop_rad_s.is_finite()
            && (10.0..=120.0).contains(&omega_prop_rad_s)
            && rho_kg_m3.is_finite()
            && rho_kg_m3 > 0.0)
        {
            return Err(Refusal {
                code: "state-invalid",
                message: format!("V {v_mps:?}, alpha {alpha_rad:?}, dc {delta_canard_rad:?}"),
                ranked_repairs: vec!["stay inside the open-loop model domain".into()],
            });
        }
        let map_err = |e: fs_wing::Refusal| Refusal {
            code: e.code,
            message: e.message,
            ranked_repairs: e.ranked_repairs,
        };
        let (panels, strips) = self.layout(delta_canard_rad)?;
        // Body-axes airspeed (u = V cos α, w = V sin α, FRD).
        let fs_v = [
            v_mps * det::cos(alpha_rad),
            0.0,
            v_mps * det::sin(alpha_rad),
        ];
        let op = InfluenceOperator::build(&panels, fs_v, rho_kg_m3).map_err(map_err)?;
        // Section closure: thin airfoil + flat-plate blend (E4.0 window).
        let camber = self.camber_ratio;
        let closure = move |_s: usize, a: f64| -> (f64, StripRegime) {
            let attached = core::f64::consts::TAU * (a + 2.0 * camber);
            let abs = a.abs();
            if abs <= 0.30 {
                (attached, StripRegime::Attached)
            } else if abs < 0.45 {
                let t = (abs - 0.30) / 0.15;
                let s = t * t * (3.0 - 2.0 * t);
                let sep = 1.98 * det::sin(a) * det::cos(a);
                (attached * (1.0 - s) + sep * s, StripRegime::Blended)
            } else {
                (1.98 * det::sin(a) * det::cos(a), StripRegime::Separated)
            }
        };
        // The disks spin at the trim's omega state (fixed control =
        // fixed THROTTLE; the prop speed is set by torque balance).
        let disks = [
            PropDisk {
                omega_rad_s: omega_prop_rad_s,
                ..self.disks[0]
            },
            PropDisk {
                omega_rad_s: omega_prop_rad_s,
                ..self.disks[1]
            },
        ];
        let coupled = coupled_prop_airframe_step(
            &op,
            &panels,
            &strips,
            &closure,
            &self.rotor,
            &disks,
            fs_v,
            rho_kg_m3,
            &self.coupling,
            None,
        )?;
        let q = 0.5 * rho_kg_m3 * v_mps * v_mps;
        let vhat = [fs_v[0] / v_mps, 0.0, fs_v[2] / v_mps];
        // Per-panel Kutta–Joukowsky forces + moments about the CG ref.
        let mut force = [0.0f64; 3];
        let mut m_y = 0.0f64;
        for (j, p) in panels.iter().enumerate() {
            let seg = [p.b[0] - p.a[0], p.b[1] - p.a[1], p.b[2] - p.a[2]];
            let fx = fs_v[1] * seg[2] - fs_v[2] * seg[1];
            let fy = fs_v[2] * seg[0] - fs_v[0] * seg[2];
            let fz = fs_v[0] * seg[1] - fs_v[1] * seg[0];
            let s = rho_kg_m3 * coupled.gamma[j];
            let f = [s * fx, s * fy, s * fz];
            let mid = [
                0.5 * (p.a[0] + p.b[0]) - self.cg_m[0],
                0.5 * (p.a[1] + p.b[1]) - self.cg_m[1],
                0.5 * (p.a[2] + p.b[2]) - self.cg_m[2],
            ];
            force[0] += f[0];
            force[1] += f[1];
            force[2] += f[2];
            m_y += mid[2] * f[0] - mid[0] * f[2];
        }
        let lift = -force[2];
        // Per-strip quarter-chord camber moments (thin airfoil cm0 = −π·camber)
        // + profile drag at the strip line.
        let cm0 = -core::f64::consts::PI * self.camber_ratio;
        let mut d_prof = 0.0;
        for st in &strips {
            let p0 = &panels[st.panel_indices[0]];
            let width = (p0.b[1] - p0.a[1]).abs();
            let c = st.chord_m;
            m_y += q * c * c * width * cm0;
            let dp = q * c * width * self.cd0_wing;
            d_prof += dp;
            let mid = [
                0.5 * (p0.a[0] + p0.b[0]) - self.cg_m[0],
                0.0,
                0.5 * (p0.a[2] + p0.b[2]) - self.cg_m[2],
            ];
            let f = [-dp * vhat[0], 0.0, -dp * vhat[2]];
            force[0] += f[0];
            force[2] += f[2];
            m_y += mid[2] * f[0] - mid[0] * f[2];
        }
        // Induced drag: declared classical term at the CG line (the
        // freestream-KJ panel forces exclude it by construction).
        let s_ref = 2.0 * self.wing_span_m * self.wing_chord_m;
        let ar = self.wing_span_m * self.wing_span_m / (s_ref / 2.0);
        let cl = lift / (q * s_ref);
        let d_ind = q * s_ref * cl * cl / (core::f64::consts::PI * ar * self.oswald_e);
        force[0] -= d_ind * vhat[0];
        force[2] -= d_ind * vhat[2];
        // Parasite ledger at its declared line.
        let ledger_rep = self.ledger.evaluate(rho_kg_m3, v_mps)?;
        let dp = ledger_rep.total_parasite_n;
        force[0] -= dp * vhat[0];
        force[2] -= dp * vhat[2];
        m_y += (self.drag_z_m - self.cg_m[2]) * (-dp * vhat[0]) - 0.0 * (-dp * vhat[2]);
        // Thrust along +x at the disk lines.
        for (k, d) in disks.iter().enumerate() {
            let t = coupled.thrust_n[k];
            force[0] += t;
            m_y += (d.center_m[2] - self.cg_m[2]) * t;
        }
        // Gravity in body axes at θ = α (level flight): mg(−sinθ, 0, cosθ).
        let w = self.gross_mass_kg * 9.80665;
        force[0] -= w * det::sin(alpha_rad);
        force[2] += w * det::cos(alpha_rad);
        let torque_imbalance_nm = fs_airscrew::engine_torque_at_prop_nm(omega_prop_rad_s)
            - 0.5 * (coupled.torque_nm[0] + coupled.torque_nm[1]);
        Ok(ForceBuildup {
            force_n: force,
            moment_y_nm: m_y,
            lift_n: lift,
            drag_n: dp + d_prof + d_ind,
            thrust_n: coupled.thrust_n,
            induced_drag_n: d_ind,
            coupled,
            torque_imbalance_nm,
        })
    }

    /// Fixed-control (fixed-THROTTLE) trim: damped Newton on
    /// (ΣFx, ΣFz, ΣMy, Q_engine−Q_prop) over (V, α, δc, ω),
    /// deterministic FD Jacobian. The prop speed is a trim STATE — the
    /// engine holds torque, not rpm.
    ///
    /// # Errors
    /// `trim-not-found` (cap or domain exit; the message carries the
    /// residual trajectory) — never a silent best-effort state.
    pub fn trim(&self, rho_kg_m3: f64, start: [f64; 4]) -> Result<TrimResult, Refusal> {
        let mut x = start; // (V, alpha, delta_c, omega)
        let mut trail: Vec<String> = Vec::new();
        let tol = [0.5, 0.5, 0.5, 0.5]; // N, N, N·m, N·m
        let mut iterations = 0u32;
        loop {
            let b = self.force_buildup(x[0], x[1], x[2], x[3], rho_kg_m3)?;
            let r = [
                b.force_n[0],
                b.force_n[2],
                b.moment_y_nm,
                b.torque_imbalance_nm,
            ];
            let rn = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]).sqrt();
            trail.push(format!("[{iterations}] |r| {rn:.3}"));
            if r[0].abs() < tol[0]
                && r[1].abs() < tol[1]
                && r[2].abs() < tol[2]
                && r[3].abs() < tol[3]
            {
                return Ok(TrimResult {
                    v_mps: x[0],
                    alpha_rad: x[1],
                    delta_canard_rad: x[2],
                    omega_prop_rad_s: x[3],
                    residuals: r,
                    iterations,
                    buildup: b,
                    design_digest: self.digest(),
                });
            }
            if iterations >= MAX_TRIM_ITERATIONS {
                return Err(Refusal {
                    code: "trim-not-found",
                    message: format!(
                        "cap {MAX_TRIM_ITERATIONS} exhausted; residual trail: {}",
                        trail.join(" ")
                    ),
                    ranked_repairs: vec![
                        "the configuration may be untrimmable at this control state".into(),
                        "widen the start or check control authority (canard area/arm)".into(),
                    ],
                });
            }
            // FD Jacobian (fixed steps — deterministic). A perturbed
            // point whose aero solve refuses (domain edge) retries with
            // the step reversed before giving up.
            let h = [0.05, 0.002, 0.002, 0.25];
            let mut jac = [[0.0f64; 4]; 4];
            for c in 0..4 {
                let mut hh = h[c];
                let mut xp = x;
                xp[c] += hh;
                let bp = match self.force_buildup(xp[0], xp[1], xp[2], xp[3], rho_kg_m3) {
                    Ok(b) => b,
                    Err(_) => {
                        hh = -h[c];
                        let mut xm = x;
                        xm[c] += hh;
                        self.force_buildup(xm[0], xm[1], xm[2], xm[3], rho_kg_m3)?
                    }
                };
                let rp = [
                    bp.force_n[0],
                    bp.force_n[2],
                    bp.moment_y_nm,
                    bp.torque_imbalance_nm,
                ];
                for (row, jrow) in jac.iter_mut().enumerate() {
                    jrow[c] = (rp[row] - r[row]) / hh;
                }
            }
            let dx = solve4(&jac, &[-r[0], -r[1], -r[2], -r[3]]).ok_or_else(|| Refusal {
                code: "trim-not-found",
                message: "singular trim Jacobian".into(),
                ranked_repairs: vec!["control authority degenerate (check canard)".into()],
            })?;
            // Damped step with domain clamps. A candidate whose aero
            // solve refuses is treated as NOT AN IMPROVEMENT (the trim
            // search must survive trial points outside the aero domain),
            // never as a fatal error.
            let mut lambda = 1.0;
            let mut accepted = false;
            for _ in 0..6 {
                let cand = [
                    (x[0] + lambda * dx[0]).clamp(5.0, 25.0),
                    (x[1] + lambda * dx[1]).clamp(-0.3, 0.35),
                    (x[2] + lambda * dx[2]).clamp(-0.5, 0.5),
                    (x[3] + lambda * dx[3]).clamp(15.0, 90.0),
                ];
                if let Ok(bc) = self.force_buildup(cand[0], cand[1], cand[2], cand[3], rho_kg_m3) {
                    let rc = [
                        bc.force_n[0],
                        bc.force_n[2],
                        bc.moment_y_nm,
                        bc.torque_imbalance_nm,
                    ];
                    let rcn =
                        (rc[0] * rc[0] + rc[1] * rc[1] + rc[2] * rc[2] + rc[3] * rc[3]).sqrt();
                    if rcn < rn {
                        x = cand;
                        accepted = true;
                        break;
                    }
                }
                lambda *= 0.5;
            }
            if !accepted {
                return Err(Refusal {
                    code: "trim-not-found",
                    message: format!(
                        "line search stalled at |r| {rn:.3}; trail: {}",
                        trail.join(" ")
                    ),
                    ranked_repairs: vec![
                        "no descent direction — likely no equilibrium in the domain".into(),
                    ],
                });
            }
            iterations += 1;
        }
    }
}

/// Deterministic 4×4 solve: Gaussian elimination, partial pivoting with
/// the FIXED tie rule (strictly-greater magnitude wins; first index on
/// ties).
fn solve4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> Option<[f64; 4]> {
    let mut m = *a;
    let mut v = *b;
    for col in 0..4 {
        let mut piv = col;
        for r in (col + 1)..4 {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if !m[piv][col].is_finite() || m[piv][col].abs() < 1e-12 {
            return None;
        }
        m.swap(col, piv);
        v.swap(col, piv);
        for r in (col + 1)..4 {
            let f = m[r][col] / m[col][col];
            for c in col..4 {
                m[r][c] -= f * m[col][c];
            }
            v[r] -= f * v[col];
        }
    }
    let mut out = [0.0f64; 4];
    for r in (0..4).rev() {
        let mut s = v[r];
        for c in (r + 1)..4 {
            s -= m[r][c] * out[c];
        }
        out[r] = s / m[r][r];
    }
    out.iter().all(|x| x.is_finite()).then_some(out)
}
