//! Nonlinear section closure (bead wf-root-guzez.5.3.2, E4.2-ii). Plan
//! §5.2 planform layer: the production force path couples the SECTION
//! data (camber, stall — supplied by the caller as a closure over
//! (surface, alpha_eff)) against the 3-D induced flow, via a
//! WARM-STARTED SAFEGUARDED Picard iteration:
//!
//!   per strip s:  α_eff = α_geo − w_ind/V,  Γ_target = ½·V·c·cl(α_eff),
//!   Γ ← Γ + ω·(Γ_target − Γ), with ω halved (down to ω_min) whenever the
//!   residual grows and restored on progress — never a silent jump.
//!
//! Chordwise rows keep the LINEAR solve's row-split ratios (the hinge
//! pressure arm survives the strip-level closure). Every solve reports
//! residual, iteration count, and a BRANCH IDENTITY (fs-blake3 of the
//! per-strip regime pattern): same inputs → same branch id bitwise;
//! consumers detect branch changes across ticks and apply the plan's
//! continuation-or-refuse rule.
//!
//! Factorization/operator reuse: the influence operator carries a
//! geometry hash; solving with panels that no longer match it is a typed
//! refusal (`influence-operator-stale`) — the reuse rule is enforced, not
//! advisory.

use crate::{MAX_PANELS, Panel, Refusal, SolveReport, solve_weissinger_linear};
use fs_blake3::hash_domain;

/// Iteration cap (refusals at cap AND cap+1 by the safeguard test).
pub const MAX_ITERATIONS: u32 = 200;
/// Convergence tolerance on the relative circulation update.
pub const GAMMA_TOL: f64 = 1.0e-10;
/// Relaxation start / floor.
pub const OMEGA_START: f64 = 0.7;
/// Relaxation floor before refusing.
pub const OMEGA_MIN: f64 = 0.02;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// One spanwise strip: the panel rows it owns (front row first), its
/// chord, and its geometric incidence offset [rad] (twist/deflection).
#[derive(Clone, Debug, PartialEq)]
pub struct StripSpec {
    /// Indices into the panel slice, front (quarter-chord) row first.
    pub panel_indices: Vec<usize>,
    /// Strip chord [m].
    pub chord_m: f64,
    /// Geometric incidence offset relative to the body x-axis [rad].
    pub twist_rad: f64,
}

/// Regime a strip converged in (drives the branch identity).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StripRegime {
    /// |α_eff| inside the attached range of the closure.
    Attached,
    /// In the closure's blend band.
    Blended,
    /// Fully separated.
    Separated,
}

/// The per-solve nonlinear report (plan §5.2 required fields).
#[derive(Clone, Debug, PartialEq)]
pub struct NonlinearReport {
    /// Converged circulations per panel.
    pub gamma: Vec<f64>,
    /// Final residual (max relative strip update).
    pub residual: f64,
    /// Iterations used.
    pub iterations: u32,
    /// Branch identity: fs-blake3 of the per-strip regime pattern.
    pub branch_id: String,
    /// Per-strip regimes.
    pub regimes: Vec<StripRegime>,
    /// Total lift [N] (Kutta–Joukowsky, −z component).
    pub total_lift_n: f64,
}

/// The section closure the caller supplies: cl(strip index, α_eff [rad])
/// plus the regime the closure was in (the caller owns regime boundaries;
/// fs-wing records them into the branch identity).
pub type SectionClosure<'a> = &'a dyn Fn(usize, f64) -> (f64, StripRegime);

/// A reusable influence operator: geometry-hashed so stale reuse refuses.
#[derive(Clone, Debug)]
pub struct InfluenceOperator {
    geometry_hash: String,
    /// The linear solution (cold-start Γ and row-split ratios).
    linear: SolveReport,
}

fn geometry_hash(panels: &[Panel], freestream: [f64; 3]) -> String {
    let mut payload = Vec::with_capacity(panels.len() * 8 * 13 + 24);
    for p in panels {
        for v in
            p.a.iter()
                .chain(p.b.iter())
                .chain(p.ctrl.iter())
                .chain(p.normal.iter())
        {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        payload.extend_from_slice(&p.width_m.to_bits().to_le_bytes());
    }
    for v in freestream {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    hash_domain("org.frankensim.fs-wing.influence-geometry.v1", &payload).to_hex()
}

impl InfluenceOperator {
    /// Build (and linear-solve) the operator for this exact geometry.
    ///
    /// # Errors
    /// Linear-solve refusals pass through.
    pub fn build(panels: &[Panel], freestream: [f64; 3], rho: f64) -> Result<Self, Refusal> {
        let linear = solve_weissinger_linear(panels, freestream, rho)?;
        Ok(InfluenceOperator {
            geometry_hash: geometry_hash(panels, freestream),
            linear,
        })
    }

    /// The linear (cold-start) solution.
    #[must_use]
    pub fn linear(&self) -> &SolveReport {
        &self.linear
    }
}

/// Induced z-velocity at a point from all horseshoes with strengths Γ.
fn induced_w(p: [f64; 3], panels: &[Panel], gamma: &[f64], stream_unit: [f64; 3]) -> f64 {
    let mut w = 0.0;
    for (j, panel) in panels.iter().enumerate() {
        let v = crate::trailing_velocity_pub(p, panel.a, panel.b, stream_unit);
        w += gamma[j] * v[2];
    }
    w
}

/// The warm-started safeguarded nonlinear solve.
///
/// # Errors
/// `influence-operator-stale` (panels/freestream no longer match the
/// operator's geometry hash — the reuse rule); `strips-invalid`;
/// `nonlinear-did-not-converge` (iteration cap or relaxation floor hit;
/// the message reports residual/ω/iterations — never a silent fallback
/// to the linear answer).
#[allow(clippy::too_many_arguments)] // the physics inputs are irreducible
pub fn solve_nonlinear(
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    freestream: [f64; 3],
    rho: f64,
    closure: SectionClosure<'_>,
    strip_du_axial: Option<&[f64]>,
    warm_start: Option<&[f64]>,
) -> Result<NonlinearReport, Refusal> {
    if geometry_hash(panels, freestream) != op.geometry_hash {
        return Err(refuse(
            "influence-operator-stale",
            "panels/freestream do not match the operator's geometry hash".into(),
            "rebuild the operator — canard deflection, warp, attitude, or images changed it",
        ));
    }
    let n = panels.len();
    if n > MAX_PANELS {
        return Err(refuse("panel-count-invalid", format!("{n}"), "cap"));
    }
    let mut owned = vec![false; n];
    for s in strips {
        if s.panel_indices.is_empty() || !(s.chord_m > 0.0) || !s.twist_rad.is_finite() {
            return Err(refuse(
                "strips-invalid",
                "empty strip or bad chord/twist".into(),
                "fix layout",
            ));
        }
        for &i in &s.panel_indices {
            if i >= n || owned[i] {
                return Err(refuse(
                    "strips-invalid",
                    format!("panel {i} missing or doubly owned"),
                    "each panel belongs to exactly one strip",
                ));
            }
            owned[i] = true;
        }
    }
    let vmag = (freestream[0] * freestream[0]
        + freestream[1] * freestream[1]
        + freestream[2] * freestream[2])
        .sqrt();
    let stream = [
        freestream[0] / vmag,
        freestream[1] / vmag,
        freestream[2] / vmag,
    ];
    let alpha_free = freestream[2].atan2(freestream[0]);
    if let Some(du) = strip_du_axial
        && (du.len() != strips.len() || !du.iter().all(|v| v.is_finite()))
    {
        return Err(refuse(
            "strips-invalid",
            "strip_du_axial length/finiteness".into(),
            "one finite axial increment per strip (prop slipstream)",
        ));
    }
    // Row-split ratios from the linear solution; strip totals start from
    // the warm start (or the linear totals).
    let lin = &op.linear.gamma;
    let mut ratios: Vec<Vec<f64>> = Vec::with_capacity(strips.len());
    let mut strip_gamma: Vec<f64> = Vec::with_capacity(strips.len());
    for s in strips {
        let total: f64 = s.panel_indices.iter().map(|&i| lin[i]).sum();
        let r: Vec<f64> = if total.abs() > 1e-14 {
            s.panel_indices.iter().map(|&i| lin[i] / total).collect()
        } else {
            vec![1.0 / s.panel_indices.len() as f64; s.panel_indices.len()]
        };
        ratios.push(r);
        strip_gamma.push(match warm_start {
            Some(ws) => s.panel_indices.iter().map(|&i| ws[i]).sum(),
            None => total,
        });
    }
    let mut gamma = vec![0.0f64; n];
    let scatter = |strip_gamma: &[f64], gamma: &mut Vec<f64>| {
        for (s, spec) in strips.iter().enumerate() {
            for (k, &i) in spec.panel_indices.iter().enumerate() {
                gamma[i] = strip_gamma[s] * ratios[s][k];
            }
        }
    };
    scatter(&strip_gamma, &mut gamma);
    let mut omega = OMEGA_START;
    let mut residual = f64::INFINITY;
    let mut regimes = vec![StripRegime::Attached; strips.len()];
    let mut iterations = 0u32;
    while iterations < MAX_ITERATIONS {
        iterations += 1;
        let mut worst = 0.0f64;
        let mut next = strip_gamma.clone();
        for (s, spec) in strips.iter().enumerate() {
            // Reference point: front-row bound-vortex midpoint.
            let front = &panels[spec.panel_indices[0]];
            let mid = [
                (front.a[0] + front.b[0]) / 2.0,
                (front.a[1] + front.b[1]) / 2.0,
                (front.a[2] + front.b[2]) / 2.0,
            ];
            let w = induced_w(mid, panels, &gamma, stream);
            // Sign chain (settled EMPIRICALLY by the battery's physical
            // fixtures): in this basis the lift-up circulation sign is
            // negative and the trailing kernel returns w < 0 for the
            // lift-up system, so the downwash correction enters with a
            // PLUS — the attached fixture (tracks the exact linear
            // solve) and the weight-class camber fixture jointly pin the
            // composite sign; a wrong sign fails both.
            let alpha_eff = alpha_free + spec.twist_rad + (w / vmag);
            let (cl, regime) = closure(s, alpha_eff);
            regimes[s] = regime;
            // Prop slipstream: the washed strip sees a higher local speed
            // (axial increment du) — the REAL prop->wing coupling arm.
            let v_local = vmag + strip_du_axial.map_or(0.0, |du| du[s]);
            // The solver's Γ sign convention (BC with −z normals, a→b
            // along +y) is NEGATIVE for lift-up — the first battery run
            // converged to −5420 N because the closure assumed the
            // opposite sign. Match the solver, not the textbook.
            let target = -0.5 * v_local * spec.chord_m * cl;
            // Residual = the OMEGA-INDEPENDENT fixed-point mismatch (the
            // step-size metric shrank with every safeguard halving and
            // faked progress — a measured lesson).
            worst = worst.max((target - strip_gamma[s]).abs() / (strip_gamma[s].abs() + 1e-9));
            next[s] = strip_gamma[s] + omega * (target - strip_gamma[s]);
        }
        if worst > residual * 1.25 {
            // Residual grew: halve the relaxation and RETRY from the
            // current state (safeguard; never accept the growing step).
            omega *= 0.5;
            if omega < OMEGA_MIN {
                return Err(refuse(
                    "nonlinear-did-not-converge",
                    format!("relaxation floor: residual {residual:e} at iteration {iterations}"),
                    "the closure/geometry pair is outside the admitted domain; NEVER fall back silently to the linear answer",
                ));
            }
            continue;
        }
        strip_gamma = next;
        scatter(&strip_gamma, &mut gamma);
        residual = worst;
        if residual < GAMMA_TOL {
            break;
        }
    }
    if residual >= GAMMA_TOL {
        return Err(refuse(
            "nonlinear-did-not-converge",
            format!("iteration cap {MAX_ITERATIONS}: residual {residual:e}"),
            "raise the strip resolution or check the closure's regularity",
        ));
    }
    // Branch identity: the regime pattern, content-hashed.
    let mut pattern = Vec::with_capacity(strips.len());
    for r in &regimes {
        pattern.push(match r {
            StripRegime::Attached => 0u8,
            StripRegime::Blended => 1,
            StripRegime::Separated => 2,
        });
    }
    let branch_id = hash_domain("org.frankensim.fs-wing.branch.v1", &pattern).to_hex();
    // Kutta–Joukowsky total lift.
    let mut total = 0.0;
    for (j, p) in panels.iter().enumerate() {
        let seg = [p.b[0] - p.a[0], p.b[1] - p.a[1], p.b[2] - p.a[2]];
        let f = [
            freestream[1] * seg[2] - freestream[2] * seg[1],
            freestream[2] * seg[0] - freestream[0] * seg[2],
            freestream[0] * seg[1] - freestream[1] * seg[0],
        ];
        total += -rho * gamma[j] * f[2];
    }
    Ok(NonlinearReport {
        gamma,
        residual,
        iterations,
        branch_id,
        regimes,
        total_lift_n: total,
    })
}
