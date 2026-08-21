//! Augmented-linearization ENGINE (bead wf-root-guzez.9.2.1, E8.2-i).
//! Plan §8: the linear model of the WHOLE control loop, not just the
//! rigid airframe — rigid 4 (u, w, q, θ) from the E4.6a-iii FD
//! linearization, actuator 2 (canard deflection + rate, friction
//! linearized at the regularized-stiction slope — locally EXACT for the
//! regularized model and documented as such), rotor 1 (torque-balance
//! ω), and, ONLY when an E4.6c pilot is supplied, pilot 4
//! (neuromuscular lag states + a first-order Padé of the reaction
//! delay — a DECLARED approximation carried in the report).
//!
//! Eigenvalues come from a deterministic dense solver (Hessenberg
//! reduction + Francis double-shift QR, fixed iteration budget, typed
//! refusal on nonconvergence). Mode-family attribution is STRUCTURAL:
//! zero one block's couplings and measure which eigenvalues move — the
//! named-attribution machinery the design-diff cards (E8.2-ii) consume.

use crate::Refusal;
use crate::aircraft::OpenLoopDesign;
use crate::canardmech::CanardMechanism;
use crate::longitudinal::{IYY_KG_M2, LongitudinalReport, Pole};
use crate::pilot::PilotWrightModel;

/// Maximum augmented dimension (absurd-input guard).
pub const MAX_DIM: usize = 16;

/// QR iteration budget per eigenvalue (typed refusal past it).
pub const QR_MAX_SWEEPS: usize = 400;

/// Rotor inertia per prop [kg·m²] (blades + shaft, Estimated class).
pub const ROTOR_INERTIA_KG_M2: f64 = 1.6;

/// Block identities for attribution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModeFamily {
    /// Rigid-body longitudinal states.
    Rigid,
    /// Canard actuator states.
    Actuator,
    /// Rotor speed state.
    Rotor,
    /// Pilot states (lag + Padé delay).
    Pilot,
}

/// The assembled engine.
#[derive(Clone, Debug, PartialEq)]
pub struct AugmentedEngine {
    /// The augmented A matrix (row-major, n×n).
    pub a: Vec<Vec<f64>>,
    /// State dimension.
    pub n: usize,
    /// State labels (frozen order).
    pub labels: Vec<&'static str>,
    /// Whether pilot columns are present.
    pub pilot_active: bool,
    /// Declared approximations carried into any downstream claim.
    pub declared_approximations: Vec<&'static str>,
}

/// One labeled eigenvalue.
#[derive(Clone, Debug, PartialEq)]
pub struct LabeledPole {
    /// The eigenvalue.
    pub pole: Pole,
    /// Dominant family by structural attribution.
    pub family: ModeFamily,
    /// Movement magnitude when the family's couplings were frozen
    /// (the attribution strength receipt).
    pub attribution_shift: f64,
}

/// Build the augmented engine about a trim.
///
/// # Errors
/// Trim/linearization refusals pass through — the engine WRAPS a trim
/// failure in `augmented-trim-refused` NAMING the limiting subsystem;
/// `augmented-dim-invalid`.
pub fn build_engine(
    design: &OpenLoopDesign,
    mech: &CanardMechanism,
    pilot: Option<&PilotWrightModel>,
    rigid: &LongitudinalReport,
    control_column: &[f64; 4],
    rho_kg_m3: f64,
) -> Result<AugmentedEngine, Refusal> {
    mech.admit()?;
    let trim = &rigid.trim;
    // Layout: [u w q theta | delta ddelta | omega | (nm p_pade x2? -> nm, pade)]
    let base = 7usize;
    let n = if pilot.is_some() { base + 4 } else { base };
    if n > MAX_DIM {
        return Err(Refusal {
            code: "augmented-dim-invalid",
            message: format!("{n} states"),
            ranked_repairs: vec!["the v1 engine caps at 16 states".into()],
        });
    }
    let mut a = vec![vec![0.0f64; n]; n];
    let mut labels: Vec<&'static str> =
        vec!["u", "w", "q", "theta", "delta_c", "ddelta_c", "omega_prop"];
    // Rigid block + control column (delta_c enters like the trim FD).
    for i in 0..4 {
        for j in 0..4 {
            a[i][j] = rigid.a[i][j];
        }
        a[i][4] = control_column[i];
    }
    // Actuator block: I δ̈ = 0.35·F_pilot − c_eff·δ̇ (+ stop springs are
    // outside the linear domain). Friction slope at the origin:
    // c_eff = viscous + coulomb/reg — the regularized model's local
    // truth (stiction-dominated; documented).
    let c_eff = mech.viscous_nm_s + mech.coulomb_nm / mech.friction_reg_rad_s;
    a[4][5] = 1.0;
    a[5][5] = -c_eff / mech.inertia_kg_m2;
    // Rotor block: I_r ω̇ = Q_engine(ω) − Q_prop(ω, V); FD both.
    let om0 = trim.omega_prop_rad_s;
    let h_om = 0.25;
    let dq_engine = (fs_airscrew::engine_torque_at_prop_nm(om0 + h_om)
        - fs_airscrew::engine_torque_at_prop_nm(om0 - h_om))
        / (2.0 * h_om);
    let qp = |om: f64, v: f64| -> Result<f64, Refusal> {
        fs_airscrew::bemt_solve(&design.rotor, rho_kg_m3, v, om)
            .map(|s| s.torque_nm)
            .map_err(|e| Refusal {
                code: e.code,
                message: e.message,
                ranked_repairs: e.ranked_repairs,
            })
    };
    let v0 = trim.v_mps;
    let dq_om = (qp(om0 + h_om, v0)? - qp(om0 - h_om, v0)?) / (2.0 * h_om);
    let dq_v = (qp(om0, v0 + 0.05)? - qp(om0, v0 - 0.05)?) / 0.1;
    a[6][6] = (dq_engine - dq_om) / ROTOR_INERTIA_KG_M2;
    a[6][0] = -dq_v / ROTOR_INERTIA_KG_M2; // u ~ V at small alpha
    // Rotor -> rigid: thrust sensitivity (two props).
    let dt_om = {
        let tp = |om: f64| -> Result<f64, Refusal> {
            fs_airscrew::bemt_solve(&design.rotor, rho_kg_m3, v0, om)
                .map(|s| s.thrust_n)
                .map_err(|e| Refusal {
                    code: e.code,
                    message: e.message,
                    ranked_repairs: e.ranked_repairs,
                })
        };
        (tp(om0 + h_om)? - tp(om0 - h_om)?) / (2.0 * h_om)
    };
    a[0][6] = 2.0 * dt_om / design.gross_mass_kg;
    // Thrust line below CG: pitch coupling.
    let z_arm = design.disks[0].center_m[2] - design.cg_m[2];
    a[2][6] = 2.0 * dt_om * z_arm / IYY_KG_M2;
    let mut declared: Vec<&'static str> = vec![
        "actuator friction linearized at the regularized stiction slope (viscous + coulomb/reg)",
        "rotor state couples through u only (small-alpha V identification)",
    ];
    let mut pilot_active = false;
    if let Some(p) = pilot {
        p.admit()?;
        pilot_active = true;
        labels.extend(["pilot_nm", "pilot_pade", "pilot_nm_lat", "pilot_pade_lat"]);
        declared.push("pilot reaction delay approximated by first-order Pade (declared)");
        let g = &p.gains;
        let tau_nm = g.neuromuscular_tau_s;
        let t_delay = (g.reaction_ticks as f64) / crate::perception::PERCEPTION_HZ;
        // Longitudinal command u_cmd = -k_theta*theta - k_q*q.
        // Pade-1: pade' = (2/T)(u_cmd - pade); u_delayed = 2*pade - u_cmd.
        // nm' = (u_delayed - nm)/tau_nm.
        // Pilot force F = k_pos*(nm - delta) - k_rate*ddelta enters the
        // actuator: delta'' += 0.35*F/I.
        let (i_nm, i_pd) = (7usize, 8usize);
        let two_over_t = 2.0 / t_delay.max(1.0 / 120.0);
        a[i_pd][3] += two_over_t * (-g.k_theta);
        a[i_pd][2] += two_over_t * (-g.k_q);
        a[i_pd][i_pd] = -two_over_t;
        // u_delayed = 2*pade - u_cmd:
        a[i_nm][i_pd] += 2.0 / tau_nm;
        a[i_nm][3] += -(-g.k_theta) / tau_nm;
        a[i_nm][2] += -(-g.k_q) / tau_nm;
        a[i_nm][i_nm] = -1.0 / tau_nm;
        let gain = mech.lever_gain_nm_per_n / mech.inertia_kg_m2;
        a[5][i_nm] += gain * g.k_lever_pos;
        a[5][4] += -gain * g.k_lever_pos;
        a[5][5] += -gain * g.k_lever_rate;
        // Lateral pilot states are carried but decoupled at this trim
        // (longitudinal engine); they contribute their own poles.
        let (j_nm, j_pd) = (9usize, 10usize);
        a[j_pd][j_pd] = -two_over_t;
        a[j_nm][j_pd] = 2.0 / tau_nm;
        a[j_nm][j_nm] = -1.0 / tau_nm;
    }
    Ok(AugmentedEngine {
        a,
        n,
        labels,
        pilot_active,
        declared_approximations: declared,
    })
}

/// Wrap a trim refusal with the limiting-subsystem name (the typed
/// trim-refusal path of the DONE-WHEN).
#[must_use]
pub fn wrap_trim_refusal(inner: Refusal) -> Refusal {
    let subsystem = if inner.message.contains("residual trail") || inner.message.contains("|r|") {
        "canard-authority (trim moment channel)"
    } else {
        "aerodynamic solve"
    };
    Refusal {
        code: "augmented-trim-refused",
        message: format!("limiting subsystem: {subsystem}; inner: {}", inner.message),
        ranked_repairs: inner.ranked_repairs,
    }
}

/// Eigenvalues of the engine (deterministic Hessenberg + Francis QR).
///
/// # Errors
/// `eig-did-not-converge` (fixed sweep budget exhausted).
pub fn eigenvalues(engine: &AugmentedEngine) -> Result<Vec<Pole>, Refusal> {
    eig_dense(&engine.a)
}

/// Structural mode-family attribution: for each family, zero its
/// coupling rows/columns and match eigenvalues by nearest distance —
/// each pole is labeled with the family whose freezing moved it most.
///
/// # Errors
/// Eigensolver refusals pass through.
pub fn attribute_modes(engine: &AugmentedEngine) -> Result<Vec<LabeledPole>, Refusal> {
    let base = eig_dense(&engine.a)?;
    let families: &[(ModeFamily, std::ops::Range<usize>)] = if engine.pilot_active {
        &[
            (ModeFamily::Rigid, 0..4),
            (ModeFamily::Actuator, 4..6),
            (ModeFamily::Rotor, 6..7),
            (ModeFamily::Pilot, 7..11),
        ]
    } else {
        &[
            (ModeFamily::Rigid, 0..4),
            (ModeFamily::Actuator, 4..6),
            (ModeFamily::Rotor, 6..7),
        ]
    };
    let mut best: Vec<(ModeFamily, f64)> = vec![(ModeFamily::Rigid, -1.0); base.len()];
    for (family, range) in families {
        // Freeze: zero the OFF-DIAGONAL couplings into/out of the block
        // (keep the block's own dynamics — freezing asks 'what if this
        // block were disconnected').
        let mut frozen = engine.a.clone();
        for i in 0..engine.n {
            for j in 0..engine.n {
                let i_in = range.contains(&i);
                let j_in = range.contains(&j);
                if i_in != j_in {
                    frozen[i][j] = 0.0;
                }
            }
        }
        let moved = eig_dense(&frozen)?;
        // Distance from each base pole to the nearest frozen pole.
        for (k, b) in base.iter().enumerate() {
            let d = moved
                .iter()
                .map(|m| ((m.re - b.re).powi(2) + (m.im - b.im).powi(2)).sqrt())
                .fold(f64::INFINITY, f64::min);
            if d > best[k].1 {
                best[k] = (*family, d);
            }
        }
    }
    // Structural-membership fallback (E6.2 hardening): a block that is
    // DECOUPLED by design (the m_aero = 0 actuator tier) barely moves
    // under ANY freeze, so its winner rests on sub-ulp noise and can
    // flip with unrelated code changes (measured 2026-08-21). Below the
    // resolution floor, attribute by which family's ISOLATED block owns
    // the pole instead — exact for decoupled blocks, inert otherwise.
    const SHIFT_FLOOR: f64 = 1.0e-9;
    let mut labeled: Vec<(Pole, ModeFamily, f64)> = base
        .iter()
        .zip(best.iter())
        .map(|(p, (f, s))| (*p, *f, *s))
        .collect();
    for (pole, family, shift) in &mut labeled {
        if *shift >= SHIFT_FLOOR {
            continue;
        }
        let mut best_member: Option<(ModeFamily, f64)> = None;
        for (fam, range) in families {
            let mut isolated = engine.a.clone();
            for i in 0..engine.n {
                for j in 0..engine.n {
                    if !(range.contains(&i) && range.contains(&j)) {
                        isolated[i][j] = 0.0;
                    }
                }
            }
            let iso = eig_dense(&isolated)?;
            let d = iso
                .iter()
                .map(|m| ((m.re - pole.re).powi(2) + (m.im - pole.im).powi(2)).sqrt())
                .fold(f64::INFINITY, f64::min);
            if best_member.is_none_or(|(_, bd)| d < bd) {
                best_member = Some((*fam, d));
            }
        }
        if let Some((fam, d)) = best_member {
            // Only claim membership when the isolated block genuinely
            // owns the pole (an isolated-system pole sits on it).
            if d < 1.0e-6 {
                *family = fam;
            }
        }
    }
    Ok(labeled
        .into_iter()
        .map(|(pole, family, shift)| LabeledPole {
            pole,
            family,
            attribution_shift: shift,
        })
        .collect())
}

/// Dense real eigenvalues: Hessenberg reduction (Householder) followed
/// by Francis double-shift QR with deterministic deflation. Public for
/// battery oracles.
///
/// # Errors
/// `eig-did-not-converge`, `eig-input-invalid`.
pub fn eig_dense(a: &[Vec<f64>]) -> Result<Vec<Pole>, Refusal> {
    let n = a.len();
    if n == 0 || n > MAX_DIM || a.iter().any(|r| r.len() != n) {
        return Err(Refusal {
            code: "eig-input-invalid",
            message: format!("{n}x? matrix"),
            ranked_repairs: vec![format!("square, 1..={MAX_DIM}")],
        });
    }
    if a.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Refusal {
            code: "eig-input-invalid",
            message: "non-finite entry".into(),
            ranked_repairs: vec!["check the assembled couplings".into()],
        });
    }
    let mut h: Vec<Vec<f64>> = a.to_vec();
    // Householder reduction to upper Hessenberg.
    for k in 0..n.saturating_sub(2) {
        let mut x: Vec<f64> = (k + 1..n).map(|i| h[i][k]).collect();
        let alpha = -x[0].signum() * x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if alpha == 0.0 || !alpha.is_finite() {
            continue;
        }
        x[0] -= alpha;
        let vnorm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if vnorm < 1e-300 {
            continue;
        }
        let v: Vec<f64> = x.iter().map(|e| e / vnorm).collect();
        // H = (I-2vv^T) H (I-2vv^T) on the trailing block.
        for col in 0..n {
            let dot: f64 = (0..v.len()).map(|i| v[i] * h[k + 1 + i][col]).sum();
            for i in 0..v.len() {
                h[k + 1 + i][col] -= 2.0 * v[i] * dot;
            }
        }
        for row in 0..n {
            let dot: f64 = (0..v.len()).map(|j| h[row][k + 1 + j] * v[j]).sum();
            for j in 0..v.len() {
                h[row][k + 1 + j] -= 2.0 * v[j] * dot;
            }
        }
    }
    // Francis DOUBLE-SHIFT QR with bulge chasing on the Hessenberg form
    // (real arithmetic converges complex pairs; single real shifts stall
    // on a pair inside a 3-block — measured). Deterministic exceptional
    // shifts every 10 stalled sweeps (the classical trick).
    let mut poles: Vec<Pole> = Vec::with_capacity(n);
    let mut hi = n;
    let mut sweeps_on_block = 0usize;
    let mut total_sweeps = 0usize;
    while hi > 0 {
        if hi == 1 {
            poles.push(Pole {
                re: h[0][0],
                im: 0.0,
            });
            break;
        }
        // Deflation scan: find the highest split point.
        let mut lo = 0;
        for i in (1..hi).rev() {
            let tiny = 1e-13 * (h[i - 1][i - 1].abs() + h[i][i].abs() + 1e-300);
            if h[i][i - 1].abs() < tiny {
                h[i][i - 1] = 0.0;
                lo = i;
                break;
            }
        }
        let size = hi - lo;
        if size == 1 {
            poles.push(Pole {
                re: h[hi - 1][hi - 1],
                im: 0.0,
            });
            hi -= 1;
            sweeps_on_block = 0;
            continue;
        }
        if size == 2 {
            let (p_, q_, r_, s_) = (
                h[hi - 2][hi - 2],
                h[hi - 2][hi - 1],
                h[hi - 1][hi - 2],
                h[hi - 1][hi - 1],
            );
            let tr = p_ + s_;
            let det = p_ * s_ - q_ * r_;
            let disc = tr * tr / 4.0 - det;
            if disc >= 0.0 {
                let root = disc.sqrt();
                poles.push(Pole {
                    re: tr / 2.0 + root,
                    im: 0.0,
                });
                poles.push(Pole {
                    re: tr / 2.0 - root,
                    im: 0.0,
                });
            } else {
                let im = (-disc).sqrt();
                poles.push(Pole { re: tr / 2.0, im });
                poles.push(Pole {
                    re: tr / 2.0,
                    im: -im,
                });
            }
            hi -= 2;
            sweeps_on_block = 0;
            continue;
        }
        if total_sweeps >= QR_MAX_SWEEPS {
            return Err(Refusal {
                code: "eig-did-not-converge",
                message: format!("QR budget {QR_MAX_SWEEPS} exhausted at block size {size}"),
                ranked_repairs: vec!["scaling pathology — inspect the assembled matrix".into()],
            });
        }
        total_sweeps += 1;
        sweeps_on_block += 1;
        // Double-shift parameters from the trailing 2x2 (or the
        // deterministic exceptional shift on stall).
        let (tr, det) = if sweeps_on_block % 11 == 10 {
            let w = h[hi - 1][hi - 2].abs() + h[hi - 2][hi - 3].abs();
            (1.5 * w, w * w)
        } else {
            let (p_, q_, r_, s_) = (
                h[hi - 2][hi - 2],
                h[hi - 2][hi - 1],
                h[hi - 1][hi - 2],
                h[hi - 1][hi - 1],
            );
            (p_ + s_, p_ * s_ - q_ * r_)
        };
        // First column of (H - aI)(H - bI) with a+b = tr, ab = det.
        let mut x = h[lo][lo] * h[lo][lo] + h[lo][lo + 1] * h[lo + 1][lo] - tr * h[lo][lo] + det;
        let mut y = h[lo + 1][lo] * (h[lo][lo] + h[lo + 1][lo + 1] - tr);
        let mut z = if lo + 2 < hi {
            h[lo + 2][lo + 1] * h[lo + 1][lo]
        } else {
            0.0
        };
        // Bulge chase.
        for k in lo..(hi - 2) {
            // Householder for [x, y, z].
            let scale = x.abs() + y.abs() + z.abs();
            if scale < 1e-300 {
                x = h[k + 1][k];
                y = h[k + 2][k];
                z = if k + 3 < hi { h[k + 3][k] } else { 0.0 };
                continue;
            }
            let (xs, ys, zs) = (x / scale, y / scale, z / scale);
            let alpha = -(xs).signum() * (xs * xs + ys * ys + zs * zs).sqrt();
            let mut v = [xs - alpha, ys, zs];
            let vn = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if vn > 1e-300 {
                v = [v[0] / vn, v[1] / vn, v[2] / vn];
                let rows = if k + 3 <= hi { 3 } else { 2 };
                // Left: rows k..k+rows over cols max(lo,k-1)..n... use 0..n
                // (safe; entries outside stay zero-consistent).
                for col in k.saturating_sub(1)..n {
                    let mut dot = 0.0;
                    for (ri, vv) in v.iter().enumerate().take(rows) {
                        dot += vv * h[k + ri][col];
                    }
                    for (ri, vv) in v.iter().enumerate().take(rows) {
                        h[k + ri][col] -= 2.0 * vv * dot;
                    }
                }
                // Right: cols k..k+rows over rows 0..min(hi, k+4).
                let rmax = core::cmp::min(hi, k + 4);
                for row in 0..rmax {
                    let mut dot = 0.0;
                    for (ci, vv) in v.iter().enumerate().take(rows) {
                        dot += h[row][k + ci] * vv;
                    }
                    for (ci, vv) in v.iter().enumerate().take(rows) {
                        h[row][k + ci] -= 2.0 * vv * dot;
                    }
                }
            }
            x = h[k + 1][k];
            y = h[k + 2][k];
            z = if k + 3 < hi { h[k + 3][k] } else { 0.0 };
        }
        // Final 2x1 Givens to restore Hessenberg at the tail.
        let k = hi - 2;
        let (gx, gz) = (h[k][k - 1], h[k + 1][k - 1]);
        if k > lo {
            let rr = gx.hypot(gz);
            if rr > 1e-300 {
                let (c, s_) = (gx / rr, gz / rr);
                for col in (k - 1)..n {
                    let (t1, t2) = (h[k][col], h[k + 1][col]);
                    h[k][col] = c * t1 + s_ * t2;
                    h[k + 1][col] = -s_ * t1 + c * t2;
                }
                for row in 0..hi {
                    let (t1, t2) = (h[row][k], h[row][k + 1]);
                    h[row][k] = c * t1 + s_ * t2;
                    h[row][k + 1] = -s_ * t1 + c * t2;
                }
            }
        }
    }
    // Canonical order (fixed tie rule).
    poles.sort_by(|a, b| {
        b.re.partial_cmp(&a.re)
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(
                a.im.partial_cmp(&b.im)
                    .unwrap_or(core::cmp::Ordering::Equal),
            )
    });
    Ok(poles)
}
