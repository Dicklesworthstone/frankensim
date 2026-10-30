//! Nonlinear section closure (bead wf-root-guzez.5.3.2, E4.2-ii). Plan
//! §5.2 planform layer: the production force path couples the SECTION
//! data (camber, stall — supplied by the caller as a closure over
//! (surface, alpha_eff)) against the 3-D induced flow, via a
//! WARM-STARTED SAFEGUARDED Picard iteration:
//!
//!   per strip s:  α_eff = α_geo − w_ind/V,  Γ_target = ½·V·c·cl(α_eff),
//!   Γ ← Γ + ω·(Γ_target − Γ), with ω halved (down to ω_min) whenever the
//!   residual repeatedly grows or an opposing residual mode stops contracting.
//!   Halving stays on the current continuation path — never a silent jump.
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
use fs_math::det;

/// Explicit bounded continuation of an accepted section solution.
pub mod continuation;

/// Iteration cap (refusals at cap AND cap+1 by the safeguard test).
// 600, not 200: the full wing+canard multisurface system converges
// linearly at ~0.94/iteration through the cross-coupling (measured in
// E4.6a-ii), needing ~370 iterations to GAMMA_TOL; single-surface cases
// stay well under 100.
pub const MAX_ITERATIONS: u32 = 600;
/// Convergence tolerance on the relative circulation update.
pub const GAMMA_TOL: f64 = 1.0e-10;
/// Relaxation start / floor.
pub const OMEGA_START: f64 = 0.7;
/// Relaxation floor before refusing.
pub const OMEGA_MIN: f64 = 0.02;
// The original row-share fallback floor. Explicit continuation refuses a
// multirow split at this floor instead of treating equal shares as resolved.
const ROW_SUM_FLOOR: f64 = 1e-14;

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
    solve_nonlinear_impl(
        op,
        panels,
        strips,
        freestream,
        rho,
        closure,
        None,
        strip_du_axial,
        warm_start,
        None,
    )
}

/// Nonlinear closure with one incident body-frame relative-air velocity per
/// strip. Axial prop wash adds to the local x component BEFORE computing its
/// angle of attack and chord-plane speed; both the circulation and returned
/// Kutta–Joukowsky lift use that same local velocity. Spanwise velocity does not
/// artificially increase the 2-D section speed. The common `freestream` still
/// defines the frozen straight trailing-leg geometry and linear row ratios.
/// This is a quasi-steady strip model, not chord-resolved unsteady gust memory.
///
/// # Errors
/// As [`solve_nonlinear`]; `incident-flow-invalid` for a count mismatch,
/// non-finite local velocity, or non-positive local chordwise flow after wash.
#[allow(clippy::too_many_arguments)]
pub fn solve_nonlinear_with_incident(
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    freestream: [f64; 3],
    rho: f64,
    closure: SectionClosure<'_>,
    strip_incident_mps: &[[f64; 3]],
    strip_du_axial: Option<&[f64]>,
    warm_start: Option<&[f64]>,
) -> Result<NonlinearReport, Refusal> {
    solve_nonlinear_impl(
        op,
        panels,
        strips,
        freestream,
        rho,
        closure,
        Some(strip_incident_mps),
        strip_du_axial,
        warm_start,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn solve_nonlinear_impl(
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    freestream: [f64; 3],
    rho: f64,
    closure: SectionClosure<'_>,
    strip_incident_mps: Option<&[[f64; 3]]>,
    strip_du_axial: Option<&[f64]>,
    warm_start: Option<&[f64]>,
    mut continuation: Option<&mut continuation::Guard<'_>>,
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
    if warm_start
        .is_some_and(|values| values.len() != n || !values.iter().all(|value| value.is_finite()))
    {
        return Err(refuse(
            "strips-invalid",
            "warm circulation must have one finite value per panel".into(),
            "supply a complete accepted circulation on the same panel ownership",
        ));
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
    let alpha_free = det::atan2(freestream[2], freestream[0]);
    if let Some(du) = strip_du_axial
        && (du.len() != strips.len() || !du.iter().all(|v| v.is_finite()))
    {
        return Err(refuse(
            "strips-invalid",
            "strip_du_axial length/finiteness".into(),
            "one finite axial increment per strip (prop slipstream)",
        ));
    }
    let local_velocity = if let Some(incident) = strip_incident_mps {
        if incident.len() != strips.len() || strips.is_empty() || owned.iter().any(|v| !v) {
            return Err(refuse(
                "incident-flow-invalid",
                "one incident velocity per strip and complete panel ownership are required".into(),
                "sample every owning strip at the same physical tick",
            ));
        }
        let mut local = Vec::with_capacity(strips.len());
        for (s, velocity) in incident.iter().enumerate() {
            let mut velocity = *velocity;
            velocity[0] += strip_du_axial.map_or(0.0, |du| du[s]);
            let speed = (velocity[0] * velocity[0] + velocity[2] * velocity[2]).sqrt();
            if !velocity.iter().all(|v| v.is_finite())
                || velocity[0] <= 0.0
                || !(speed > 0.0 && speed.is_finite())
            {
                return Err(refuse(
                    "incident-flow-invalid",
                    format!("strip {s} has invalid/reversed chordwise flow {velocity:?}"),
                    "this section model requires finite positive chordwise relative flow",
                ));
            }
            local.push(velocity);
        }
        Some(local)
    } else {
        None
    };
    // All incident values and vortex geometry are fixed during this solve.
    // Cache the expensive Biot–Savart coefficients once; the Picard loop
    // retains the original gamma[j] * coefficient[j] accumulation order.
    let mut downwash = Vec::with_capacity(strips.len());
    let mut section_flow = Vec::with_capacity(strips.len());
    let mut panel_strip = vec![0; n];
    for (s, strip) in strips.iter().enumerate() {
        let front = &panels[strip.panel_indices[0]];
        let mid = [
            (front.a[0] + front.b[0]) / 2.0,
            (front.a[1] + front.b[1]) / 2.0,
            (front.a[2] + front.b[2]) / 2.0,
        ];
        downwash.push(
            panels
                .iter()
                .map(|p| crate::trailing_velocity_pub(mid, p.a, p.b, stream)[2])
                .collect::<Vec<_>>(),
        );
        for &panel in &strip.panel_indices {
            panel_strip[panel] = s;
        }
        section_flow.push(match &local_velocity {
            Some(local) => {
                let v = local[s];
                let speed = (v[0] * v[0] + v[2] * v[2]).sqrt();
                (det::atan2(v[2], v[0]), speed, speed)
            }
            None => (
                alpha_free,
                vmag,
                vmag + strip_du_axial.map_or(0.0, |du| du[s]),
            ),
        });
    }
    // Row-split ratios from the linear solution; strip totals start from
    // the warm start (or the linear totals).
    let lin = &op.linear.gamma;
    let mut ratios: Vec<Vec<f64>> = Vec::with_capacity(strips.len());
    let mut strip_gamma: Vec<f64> = Vec::with_capacity(strips.len());
    for s in strips {
        let total: f64 = s.panel_indices.iter().map(|&i| lin[i]).sum();
        let r: Vec<f64> = if total.abs() > ROW_SUM_FLOOR {
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
    if let Some(guard) = &mut continuation {
        guard.prepare(panels, strips, freestream, lin)?;
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
    let mut residual = f64::INFINITY;
    let mut regimes = vec![StripRegime::Attached; strips.len()];
    let mut iterations = 0u32;
    if let Some(guard) = &mut continuation {
        let solved = guard.solve(
            &strip_gamma,
            strips,
            &ratios,
            &downwash,
            &section_flow,
            closure,
            n,
        )?;
        strip_gamma = solved.gamma;
        scatter(&strip_gamma, &mut gamma);
        residual = solved.residual;
        regimes = solved.regimes;
        iterations = solved.iterations;
    } else {
        let mut omega = OMEGA_START;
        let mut growth_streak = 0u32;
        let mut next = vec![0.0; strips.len()];
        let mut previous_mismatch = vec![0.0f64; strips.len()];
        let mut oscillation_streak = 0u32;
        while iterations < MAX_ITERATIONS {
            iterations += 1;
            let mut worst = 0.0f64;
            // Common scaling keeps the residual angle/norm comparison finite
            // without changing the convergence metric or the circulation update.
            let mut comparison_scale = 0.0f64;
            let mut previous_norm_sq = 0.0f64;
            let mut current_norm_sq = 0.0f64;
            let mut direction_dot = 0.0f64;
            for (s, spec) in strips.iter().enumerate() {
                let mut w = 0.0;
                for (j, coefficient) in downwash[s].iter().enumerate() {
                    let term = gamma[j] * coefficient;
                    w += term;
                }
                // Sign chain (settled EMPIRICALLY by the battery's physical
                // fixtures): in this basis the lift-up circulation sign is
                // negative and the trailing kernel returns w < 0 for the
                // lift-up system, so the downwash correction enters with a
                // PLUS — the attached fixture (tracks the exact linear
                // solve) and the weight-class camber fixture jointly pin the
                // composite sign; a wrong sign fails both.
                let (alpha_incident, speed_induction, speed_circulation) = section_flow[s];
                let alpha_eff = alpha_incident + spec.twist_rad + (w / speed_induction);
                let (cl, regime) = closure(s, alpha_eff);
                if !cl.is_finite() {
                    return Err(refuse(
                        "nonlinear-did-not-converge",
                        format!(
                            "closure returned non-finite cl at strip {s} (alpha_eff {alpha_eff})"
                        ),
                        "the closure must be finite over the reachable alpha range; NEVER fall back silently to the linear answer",
                    ));
                }
                regimes[s] = regime;
                // Prop slipstream: the washed strip sees a higher local speed
                // (axial increment du) — the REAL prop->wing coupling arm.
                let v_local = speed_circulation;
                // The solver's Γ sign convention (BC with −z normals, a→b
                // along +y) is NEGATIVE for lift-up — the first battery run
                // converged to −5420 N because the closure assumed the
                // opposite sign. Match the solver, not the textbook.
                let target = -0.5 * v_local * spec.chord_m * cl;
                let mismatch = target - strip_gamma[s];
                if !mismatch.is_finite() {
                    return Err(refuse(
                        "nonlinear-did-not-converge",
                        format!("non-finite circulation mismatch at strip {s}"),
                        "the closure and circulation must remain finite; NEVER fall back silently to the linear answer",
                    ));
                }
                let prior = previous_mismatch[s];
                let component_scale = prior.abs().max(mismatch.abs());
                if component_scale > comparison_scale {
                    let ratio = comparison_scale / component_scale;
                    let squared = ratio * ratio;
                    previous_norm_sq *= squared;
                    current_norm_sq *= squared;
                    direction_dot *= squared;
                    comparison_scale = component_scale;
                }
                if comparison_scale > 0.0 {
                    let old_scaled = prior / comparison_scale;
                    let new_scaled = mismatch / comparison_scale;
                    previous_norm_sq += old_scaled * old_scaled;
                    current_norm_sq += new_scaled * new_scaled;
                    direction_dot += old_scaled * new_scaled;
                }
                previous_mismatch[s] = mismatch;
                // Residual = the OMEGA-INDEPENDENT fixed-point mismatch (the
                // step-size metric shrank with every safeguard halving and
                // faked progress — a measured lesson).
                worst = worst.max((target - strip_gamma[s]).abs() / (strip_gamma[s].abs() + 1e-9));
                // Statement-split (guzez.7.2.1): a + omega*(t-a) is a
                // fusable mul-add; bind the update first.
                let upd = omega * (target - strip_gamma[s]);
                next[s] = strip_gamma[s] + upd;
            }
            // The scalar relative residual can alternate as a strip passes near
            // zero, and an unstable mode may grow less than 25% each time. Both
            // evade the legacy two-growth guard indefinitely. Two successive
            // residual reversals (angle >120 degrees) with no norm contraction
            // identify that oscillatory mode independently of the denominator.
            // The epsilon allowance recognizes a roundoff-level two-cycle; it
            // does not relax GAMMA_TOL or accept any circulation as converged.
            let opposing = direction_dot < 0.0
                && 4.0 * direction_dot * direction_dot > previous_norm_sq * current_norm_sq;
            let noncontracting = previous_norm_sq > 0.0
                && current_norm_sq >= previous_norm_sq * (1.0 - 64.0 * f64::EPSILON);
            if worst >= GAMMA_TOL && opposing && noncontracting {
                oscillation_streak += 1;
            } else {
                oscillation_streak = 0;
            }
            if oscillation_streak >= 2 {
                omega *= 0.5;
                oscillation_streak = 0;
                growth_streak = 0;
                if omega < OMEGA_MIN {
                    return Err(refuse(
                        "nonlinear-did-not-converge",
                        format!(
                            "oscillation relaxation floor: residual {worst:e} at iteration {iterations}"
                        ),
                        "the current branch did not contract within the admitted relaxation range; NEVER fall back silently to the linear answer",
                    ));
                }
                // Reevaluate the same current state with the smaller relaxation;
                // no extrapolation, cold restart, or search for another branch.
                continue;
            }
            if worst > residual * 1.25 {
                growth_streak += 1;
                if growth_streak >= 2 {
                    // Two CONSECUTIVE growths: halve the relaxation and RETRY
                    // from the current state (safeguard; never accept a
                    // persistently growing step). A SINGLE growth is accepted:
                    // in a multisurface system a settling surface legitimately
                    // forces a one-iteration residual spike on the surfaces it
                    // washes (canard -> wing, measured in E4.6a-ii), and
                    // halving on first growth floors omega during a convergent
                    // joint transient. A truly divergent closure grows every
                    // iteration and still hits the floor refusal.
                    omega *= 0.5;
                    growth_streak = 0;
                    if omega < OMEGA_MIN {
                        return Err(refuse(
                            "nonlinear-did-not-converge",
                            format!(
                                "relaxation floor: residual {residual:e} at iteration {iterations}"
                            ),
                            "the closure/geometry pair is outside the admitted domain; NEVER fall back silently to the linear answer",
                        ));
                    }
                    continue;
                }
            } else {
                growth_streak = 0;
            }
            core::mem::swap(&mut strip_gamma, &mut next);
            scatter(&strip_gamma, &mut gamma);
            residual = worst;
            if residual < GAMMA_TOL {
                break;
            }
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
        let freestream = local_velocity
            .as_ref()
            .map_or(freestream, |local| local[panel_strip[j]]);
        let seg = [p.b[0] - p.a[0], p.b[1] - p.a[1], p.b[2] - p.a[2]];
        // Statement-split (guzez.7.2.1): no within-expression
        // mul-sub/mul-add for LLVM to contract into FMA.
        let fx1 = freestream[1] * seg[2];
        let fx2 = freestream[2] * seg[1];
        let fy1 = freestream[2] * seg[0];
        let fy2 = freestream[0] * seg[2];
        let fz1 = freestream[0] * seg[1];
        let fz2 = freestream[1] * seg[0];
        let f = [fx1 - fx2, fy1 - fy2, fz1 - fz2];
        let term = -rho * gamma[j] * f[2];
        total += term;
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
