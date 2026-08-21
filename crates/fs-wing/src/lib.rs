//! fs-wing — coupled multisurface lifting-surface aerodynamics (L3).
//! Bead frankensim-wf-root-guzez.5.3.1 (E4.2-i, Wright Flyer program).
//!
//! Spec: COMPREHENSIVE_PLAN §5.2 (ROUND 6 steady state). E4.2-i ships the
//! WeissingerLLinear machinery: multisurface panel layout (both main
//! wings, both canard planes, vertical surfaces; ≥2 chordwise rows on
//! wing and canard — the hinge pressure arm), horseshoe-vortex influence
//! assembly with the boundary condition applied at the Weissinger-L
//! three-quarter-chord point, and one declared deterministic dense
//! factorization with a condition estimate on every solve.
//!
//! Role boundary (plan law): WeissingerLLinear is an EXACT FIXTURE and an
//! admission-selected debug mode — never the production force path
//! (E4.2-ii's nonlinear closure owns that) and never entered
//! automatically after nonlinear failure. No scalar biplane factor
//! anywhere: the biplane effect EMERGES from the influence matrix, and
//! the battery checks it against the Munk-class trend (a5 verification
//! fixtures) rather than inserting it.
//!
//! Frame: frd-body-v1 (+x forward, +y right, +z down); panels carry unit
//! normals; freestream comes in as a body-frame velocity.

use fs_math::det;

pub mod hinge;
pub mod images;
pub mod nonlinear;
pub mod prescribedwake;
pub mod rom;
pub mod romreduce;

/// Panel-count cap (refusals at cap AND cap+1).
pub const MAX_PANELS: usize = 512;
/// Condition-estimate cap: beyond this the solve refuses (near-singular
/// influence operator — geometry degenerate).
pub const MAX_CONDITION_EST: f64 = 1.0e12;

/// A typed refusal (workspace law).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// Which aircraft surface a panel belongs to (diagnostics + per-surface
/// force reporting; the canard hinge arm needs per-surface rows).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceId {
    /// Lower main wing.
    WingLower,
    /// Upper main wing.
    WingUpper,
    /// Lower canard plane.
    CanardLower,
    /// Upper canard plane.
    CanardUpper,
    /// Vertical surfaces (rudder pair as one lifting surface).
    Vertical,
}

/// One horseshoe panel: bound-vortex segment A→B at the quarter chord,
/// control point at the three-quarter chord (Weissinger-L), unit normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Panel {
    /// Surface this panel belongs to.
    pub surface: SurfaceId,
    /// Bound-vortex endpoint A [m] (left/inboard).
    pub a: [f64; 3],
    /// Bound-vortex endpoint B [m].
    pub b: [f64; 3],
    /// Control point [m] (three-quarter-chord of the strip).
    pub ctrl: [f64; 3],
    /// Unit normal (frd).
    pub normal: [f64; 3],
    /// Panel span width [m] (for force integration).
    pub width_m: f64,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> f64 {
    det::sqrt(dot(a, a))
}

/// Induced velocity at `p` of a UNIT-strength horseshoe: bound segment
/// a→b plus two trailing legs to −∞ along `stream` (unit vector,
/// downstream). Classical Biot–Savart with a hard-core guard.
pub(crate) fn horseshoe_velocity(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    stream: [f64; 3],
) -> [f64; 3] {
    const CORE: f64 = 1.0e-8;
    let seg = |p1: [f64; 3], p2: [f64; 3]| -> [f64; 3] {
        let r1 = sub(p, p1);
        let r2 = sub(p, p2);
        let (n1, n2) = (norm(r1), norm(r2));
        let c = cross3(r1, r2);
        let c2 = dot(c, c);
        if c2 < CORE || n1 < CORE || n2 < CORE {
            return [0.0; 3];
        }
        let r0 = sub(p2, p1);
        let k = (dot(r0, r1) / n1 - dot(r0, r2) / n2) / (4.0 * core::f64::consts::PI * c2);
        [c[0] * k, c[1] * k, c[2] * k]
    };
    // Trailing legs: from far downstream to a, and b to far downstream.
    const FAR: f64 = 1.0e4;
    let a_far = [
        a[0] - stream[0] * FAR,
        a[1] - stream[1] * FAR,
        a[2] - stream[2] * FAR,
    ];
    let b_far = [
        b[0] - stream[0] * FAR,
        b[1] - stream[1] * FAR,
        b[2] - stream[2] * FAR,
    ];
    let v1 = seg(a_far, a);
    let v2 = seg(a, b);
    let v3 = seg(b, b_far);
    [
        v1[0] + v2[0] + v3[0],
        v1[1] + v2[1] + v3[1],
        v1[2] + v2[2] + v3[2],
    ]
}

/// The per-solve report the plan requires.
#[derive(Clone, Debug, PartialEq)]
pub struct SolveReport {
    /// Circulations per panel [m²/s].
    pub gamma: Vec<f64>,
    /// Reciprocal-condition ESTIMATE of the influence matrix (1-norm
    /// power-iteration class; an estimate, honestly labeled).
    pub condition_est: f64,
    /// Per-surface lift (force along −z_freestream-normal) [N].
    pub surface_lift_n: Vec<(SurfaceId, f64)>,
    /// Total lift [N].
    pub total_lift_n: f64,
}

/// The WeissingerLLinear solve: no-penetration at every control point,
/// Γ from one deterministic dense LU (partial pivoting with a FIXED
/// deterministic tie rule), Kutta–Joukowsky strip forces.
///
/// # Errors
/// `panel-count-invalid` (0 or above the cap, tested at cap AND cap+1);
/// `panel-invalid` (non-finite, non-unit normal, zero-width);
/// `freestream-invalid`; `influence-singular` /
/// `influence-ill-conditioned` (condition estimate beyond the cap).
pub fn solve_weissinger_linear(
    panels: &[Panel],
    freestream_mps: [f64; 3],
    rho_kg_m3: f64,
) -> Result<SolveReport, Refusal> {
    let n = panels.len();
    if n == 0 || n > MAX_PANELS {
        return Err(refuse(
            "panel-count-invalid",
            format!("{n} panels outside [1, {MAX_PANELS}]"),
            "the Tier-A target is ~80 panels",
        ));
    }
    let vmag = norm(freestream_mps);
    if !freestream_mps.iter().all(|v| v.is_finite())
        || vmag < 1.0e-6
        || !rho_kg_m3.is_finite()
        || rho_kg_m3 <= 0.0
    {
        return Err(refuse(
            "freestream-invalid",
            format!("V {freestream_mps:?}, rho {rho_kg_m3}"),
            "finite non-zero freestream; positive density",
        ));
    }
    let stream = [
        freestream_mps[0] / vmag,
        freestream_mps[1] / vmag,
        freestream_mps[2] / vmag,
    ];
    for p in panels {
        let finite =
            p.a.iter()
                .chain(p.b.iter())
                .chain(p.ctrl.iter())
                .chain(p.normal.iter())
                .all(|v| v.is_finite());
        let unit = (dot(p.normal, p.normal) - 1.0).abs() < 1.0e-9;
        if !finite || !unit || !(p.width_m > 0.0) {
            return Err(refuse(
                "panel-invalid",
                format!("{:?} panel geometry invalid", p.surface),
                "finite geometry, unit normal, positive width",
            ));
        }
    }
    // Influence matrix: A[i][j] = (velocity at ctrl_i of unit horseshoe j)·n_i.
    let mut a = vec![0.0f64; n * n];
    let mut rhs = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            let v = horseshoe_velocity(panels[i].ctrl, panels[j].a, panels[j].b, stream);
            a[i * n + j] = dot(v, panels[i].normal);
        }
        rhs[i] = -dot(freestream_mps, panels[i].normal);
    }
    // Deterministic dense LU with partial pivoting (fixed tie rule: the
    // smallest row index among equal-magnitude pivots wins).
    let mut lu = a.clone();
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut piv = k;
        let mut best = lu[perm[k] * n + k].abs();
        for r in (k + 1)..n {
            let mag = lu[perm[r] * n + k].abs();
            if mag > best {
                best = mag;
                piv = r;
            }
        }
        if best == 0.0 {
            return Err(refuse(
                "influence-singular",
                format!("zero pivot at column {k}"),
                "degenerate panel geometry (coincident strips?)",
            ));
        }
        perm.swap(k, piv);
        let pk = perm[k];
        for r in (k + 1)..n {
            let pr = perm[r];
            let f = lu[pr * n + k] / lu[pk * n + k];
            lu[pr * n + k] = f;
            for c in (k + 1)..n {
                lu[pr * n + c] -= f * lu[pk * n + c];
            }
        }
    }
    let solve_with = |lu: &[f64], perm: &[usize], b_in: &[f64]| -> Vec<f64> {
        let mut y = vec![0.0f64; n];
        for r in 0..n {
            let mut s = b_in[perm[r]];
            for c in 0..r {
                s -= lu[perm[r] * n + c] * y[c];
            }
            y[r] = s;
        }
        for r in (0..n).rev() {
            let mut s = y[r];
            for c in (r + 1)..n {
                s -= lu[perm[r] * n + c] * y[c];
            }
            y[r] = s / lu[perm[r] * n + r];
        }
        y
    };
    let gamma = solve_with(&lu, &perm, &rhs);
    // Condition ESTIMATE: ||A||_1 · ||A⁻¹ e||_1-ish via one solve on a
    // deterministic ±1 probe (Hager-class single step; an estimate).
    let a_norm1 = (0..n)
        .map(|c| (0..n).map(|r| a[r * n + c].abs()).sum::<f64>())
        .fold(0.0f64, f64::max);
    let probe: Vec<f64> = (0..n)
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    let inv_probe = solve_with(&lu, &perm, &probe);
    let inv_norm_est = inv_probe.iter().map(|v| v.abs()).sum::<f64>() / n as f64;
    let condition_est = a_norm1 * inv_norm_est;
    if !condition_est.is_finite() || condition_est > MAX_CONDITION_EST {
        return Err(refuse(
            "influence-ill-conditioned",
            format!("condition estimate {condition_est:e} beyond {MAX_CONDITION_EST:e}"),
            "geometry near-degenerate; check gap/spacing",
        ));
    }
    // Kutta–Joukowsky per strip: L' = rho·V·Γ·width along the lift
    // direction (freestream × bound-segment direction); report the
    // component opposing +z (frd lift is −z).
    let mut per: Vec<(SurfaceId, f64)> = Vec::new();
    let mut total = 0.0;
    for (j, p) in panels.iter().enumerate() {
        let seg = sub(p.b, p.a);
        let f_vec = cross3(freestream_mps, seg);
        let lift = -rho_kg_m3 * gamma[j] * f_vec[2];
        total += lift;
        match per.iter_mut().find(|(s, _)| *s == p.surface) {
            Some((_, acc)) => *acc += lift,
            None => per.push((p.surface, lift)),
        }
    }
    Ok(SolveReport {
        gamma,
        condition_est,
        surface_lift_n: per,
        total_lift_n: total,
    })
}

/// Build a flat rectangular surface's panels: `n_span` strips × `n_chord`
/// rows at height z, spanning y ∈ [−b/2, b/2], chord along −x from
/// leading edge `x_le` (frd: downstream is −x... the freestream arrives
/// along +x_body from ahead, so chord extends toward −x).
///
/// # Errors
/// `layout-invalid` (zero counts / non-finite dims).
pub fn flat_surface(
    surface: SurfaceId,
    span_m: f64,
    chord_m: f64,
    x_le: f64,
    z_m: f64,
    n_span: usize,
    n_chord: usize,
) -> Result<Vec<Panel>, Refusal> {
    if n_span == 0
        || n_chord == 0
        || !(span_m > 0.0)
        || !(chord_m > 0.0)
        || !x_le.is_finite()
        || !z_m.is_finite()
    {
        return Err(refuse(
            "layout-invalid",
            format!("{surface:?}"),
            "positive dims and counts",
        ));
    }
    let mut out = Vec::with_capacity(n_span * n_chord);
    let dy = span_m / n_span as f64;
    let dx = chord_m / n_chord as f64;
    for c in 0..n_chord {
        let x_qc = x_le - (c as f64 + 0.25) * dx; // quarter chord of the row
        let x_cp = x_le - (c as f64 + 0.75) * dx; // three-quarter (control)
        for s in 0..n_span {
            let y0 = -span_m / 2.0 + s as f64 * dy;
            out.push(Panel {
                surface,
                a: [x_qc, y0, z_m],
                b: [x_qc, y0 + dy, z_m],
                ctrl: [x_cp, y0 + dy / 2.0, z_m],
                normal: [0.0, 0.0, -1.0],
                width_m: dy,
            });
        }
    }
    Ok(out)
}

/// Crate-internal re-export for the nonlinear module (kept out of the
/// public API: the kernel is an implementation detail).
pub(crate) fn horseshoe_velocity_pub(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    stream: [f64; 3],
) -> [f64; 3] {
    horseshoe_velocity(p, a, b, stream)
}

/// Trailing-legs-only induced velocity (no bound segment): the induced-
/// angle probe of the NONLINEAR closure. Classical nonlinear lifting-line
/// evaluates the induced angle from the TRAILING sheet only — a bound
/// segment's near field (especially a strip's own second chordwise row,
/// ~1 m away) is 2-D physics the section data already contains, and
/// including it double-counts (the first battery run refused to converge
/// because of exactly this).
pub(crate) fn trailing_velocity_pub(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    stream: [f64; 3],
) -> [f64; 3] {
    let full = horseshoe_velocity(p, a, b, stream);
    // Subtract the bound-segment contribution (recompute it alone).
    let bound = {
        const CORE: f64 = 1.0e-8;
        let r1 = sub(p, a);
        let r2 = sub(p, b);
        let (n1, n2) = (norm(r1), norm(r2));
        let c = cross3(r1, r2);
        let c2 = dot(c, c);
        if c2 < CORE || n1 < CORE || n2 < CORE {
            [0.0; 3]
        } else {
            let r0 = sub(b, a);
            let k = (dot(r0, r1) / n1 - dot(r0, r2) / n2) / (4.0 * core::f64::consts::PI * c2);
            [c[0] * k, c[1] * k, c[2] * k]
        }
    };
    [full[0] - bound[0], full[1] - bound[1], full[2] - bound[2]]
}

/// Free-air induced velocity of a solved system at any point (the
/// prop-disk inflow probe for the E4.5-ii coupling).
#[must_use]
pub fn induced_velocity_free(
    p: [f64; 3],
    panels: &[Panel],
    gamma: &[f64],
    freestream_mps: [f64; 3],
) -> [f64; 3] {
    let vmag = norm(freestream_mps);
    let stream = [
        freestream_mps[0] / vmag,
        freestream_mps[1] / vmag,
        freestream_mps[2] / vmag,
    ];
    let mut out = [0.0f64; 3];
    for (j, panel) in panels.iter().enumerate() {
        let v = horseshoe_velocity(p, panel.a, panel.b, stream);
        for k in 0..3 {
            out[k] += gamma[j] * v[k];
        }
    }
    out
}
