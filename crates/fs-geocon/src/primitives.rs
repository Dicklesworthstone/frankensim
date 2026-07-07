//! The value/derivative primitives: thickness aggregation, draft
//! angles, envelopes, and certified volume. Smooth forms exist for the
//! optimizer; exact/enclosure forms exist for the ledger — both are
//! returned side by side, never conflated.

use fs_exec::Cx;
use fs_geom::{Chart, Point3, Vec3};
use fs_query::thickness_at;

/// Thickness aggregation over boundary samples.
#[derive(Debug, Clone)]
pub struct ThicknessReport {
    /// Smooth soft-minimum: the C¹ optimizer value (`≥ hard_min`,
    /// converging down as p grows).
    pub soft_min: f64,
    /// The hard sampled minimum (the ledger value).
    pub hard_min: f64,
    /// Indices of samples violating the requirement (LOCALIZATION).
    pub violating: Vec<usize>,
    /// Samples the oracle skipped (medial degeneracies), counted.
    pub skipped: u32,
}

/// Smooth minimum thickness over `samples` with mean p-norm
/// aggregation: `T_soft = (mean(t_i^{-p}))^{-1/p}` — a smooth
/// OVER-approximation of the minimum that converges DOWN to it as `p`
/// grows (exact when samples are uniform, which keeps lever
/// derivatives clean). The optimizer differentiates `soft_min`; the
/// LEDGER/verdict value is `hard_min` and the localized violation
/// list — the two are reported side by side, never conflated.
///
/// # Errors
/// [`fs_query::QueryError`] teaching errors carried through.
pub fn min_thickness_soft(
    chart: &dyn Chart,
    samples: &[Point3],
    required: f64,
    p: f64,
    cx: &Cx<'_>,
) -> Result<ThicknessReport, fs_query::QueryError> {
    let mut inv_sum = 0.0;
    let mut count = 0u32;
    let mut hard_min = f64::INFINITY;
    let mut violating = Vec::new();
    let mut skipped = 0u32;
    for (i, &s) in samples.iter().enumerate() {
        match thickness_at(chart, s, cx) {
            Ok(t) => {
                inv_sum += t.value.powf(-p);
                count += 1;
                hard_min = hard_min.min(t.value);
                if t.value < required {
                    violating.push(i);
                }
            }
            Err(fs_query::QueryError::Cancelled) => {
                return Err(fs_query::QueryError::Cancelled);
            }
            Err(_) => skipped += 1,
        }
    }
    let soft_min = if count == 0 {
        0.0
    } else {
        (inv_sum / f64::from(count)).powf(-1.0 / p)
    };
    Ok(ThicknessReport {
        soft_min,
        hard_min,
        violating,
        skipped,
    })
}

/// Draft-angle assessment against a pull direction.
#[derive(Debug, Clone)]
pub struct DraftReport {
    /// Smooth penalty: mean of squared hinges `max(sinα − n·d, 0)²`
    /// over the assessed samples (C¹ in the normals).
    pub penalty: f64,
    /// EXACT violating regions: sample indices with insufficient draft.
    pub violating: Vec<usize>,
    /// Undercuts (normals pointing AGAINST the pull): worse than mere
    /// low draft; flagged separately.
    pub undercuts: Vec<usize>,
    /// Worst deficit `sinα − n·d` observed.
    pub worst_deficit: f64,
}

/// Assess draft for the mold half pulled along `pull` (unit): surface
/// normals must satisfy `n·pull ≥ sin(min_draft)`. Samples whose
/// normals oppose the pull are undercuts. Samples nearly perpendicular
/// to the pull's mirror-half (`n·pull < −cos_tolerance`) belong to the
/// other mold half and are skipped — the v1 parting model is the plane
/// perpendicular to the pull.
///
/// # Errors
/// [`fs_query::QueryError::NoGradient`] where the chart has no normal.
pub fn draft_violations(
    chart: &dyn Chart,
    samples: &[Point3],
    pull: Vec3,
    min_draft: f64,
    cx: &Cx<'_>,
) -> Result<DraftReport, fs_query::QueryError> {
    let pn = pull.norm().max(1e-300);
    let d = pull.scale(1.0 / pn);
    let sin_a = min_draft.sin();
    let mut penalty = 0.0;
    let mut violating = Vec::new();
    let mut undercuts = Vec::new();
    let mut worst = 0.0f64;
    let mut assessed = 0u32;
    for (i, &s) in samples.iter().enumerate() {
        let sample = chart.eval(s, cx);
        let Some(g) = sample.gradient else {
            return Err(fs_query::QueryError::NoGradient {
                at: [s.x, s.y, s.z],
            });
        };
        let n = g.scale(1.0 / g.norm().max(1e-300));
        let nd = n.dot(d);
        if nd < -0.5 {
            continue; // the other mold half's face
        }
        assessed += 1;
        let deficit = sin_a - nd;
        if deficit > 0.0 {
            if nd < -1e-9 {
                undercuts.push(i);
            } else {
                violating.push(i);
            }
            penalty += deficit * deficit;
            worst = worst.max(deficit);
        }
    }
    penalty /= f64::from(assessed.max(1));
    Ok(DraftReport {
        penalty,
        violating,
        undercuts,
        worst_deficit: worst,
    })
}

/// Envelope containment assessment.
#[derive(Debug, Clone)]
pub struct EnvelopeReport {
    /// Sampled worst signed distance of the design boundary into the
    /// forbidden side (`> 0` means violation).
    pub worst: f64,
    /// Smooth log-sum-exp aggregate: `≥ worst` (conservative), within
    /// `ln(n)/β` of it (the C¹ value the optimizer differentiates).
    pub soft_worst: f64,
    /// Violating sample indices.
    pub violating: Vec<usize>,
}

/// Containment: every design-boundary sample must satisfy
/// `φ_allowed ≤ 0` (inside the allowed region). For keep-outs, pass
/// the keep-out's COMPLEMENT semantics by supplying `flip = true`
/// (violation when the sample is INSIDE the keep-out).
pub fn envelope_violation(
    allowed: &dyn Chart,
    design_boundary: &[Point3],
    beta: f64,
    flip: bool,
    cx: &Cx<'_>,
) -> EnvelopeReport {
    let mut worst = f64::NEG_INFINITY;
    let mut violating = Vec::new();
    // Sum-form log-sum-exp: (1/β)·ln(Σ exp(β·g_i)) ≥ max(g_i) — a
    // CONSERVATIVE smooth upper bound, so driving the soft value to 0
    // drives the true worst to 0 (never stops short).
    let mut acc = 0.0;
    let mut max_g = f64::NEG_INFINITY;
    let gs: Vec<f64> = design_boundary
        .iter()
        .map(|&p| {
            let sd = allowed.eval(p, cx).signed_distance;
            if flip { -sd } else { sd }
        })
        .collect();
    for g in &gs {
        max_g = max_g.max(*g);
    }
    for (i, &g) in gs.iter().enumerate() {
        worst = worst.max(g);
        if g > 0.0 {
            violating.push(i);
        }
        acc += ((g - max_g) * beta).exp();
    }
    let soft_worst = if gs.is_empty() {
        0.0
    } else {
        max_g + acc.ln() / beta
    };
    EnvelopeReport {
        worst,
        soft_worst,
        violating,
    }
}

/// A rigorous volume enclosure.
#[derive(Debug, Clone, Copy)]
pub struct VolumeEnclosure {
    /// Certain lower bound (sure-inside cells).
    pub lo: f64,
    /// Certain upper bound (lower + the uncertainty band).
    pub hi: f64,
    /// Grid step used.
    pub h: f64,
}

/// Certified volume over an EXPLICIT integration domain (fixed
/// independently of design levers, so lever derivatives see the shape
/// change, not grid realignment): sample cell centers at step `h`; cells with
/// `φ ≤ −L·h·√3/2` are SURELY inside, cells with `|φ| < L·h·√3/2` are
/// the uncertainty band. The true volume lies in `[lo, hi]` for
/// 1-Lipschitz-certified charts (the workspace SDF contract).
///
/// # Errors
/// [`fs_query::QueryError::Cancelled`].
pub fn volume_certified(
    chart: &dyn Chart,
    domain: &fs_geom::Aabb,
    h: f64,
    cx: &Cx<'_>,
) -> Result<VolumeEnclosure, fs_query::QueryError> {
    let b = domain;
    let cell = h * h * h;
    let band_half = h * 3.0f64.sqrt() / 2.0;
    let counts = |lo: f64, hi: f64| ((hi - lo) / h).ceil().max(1.0) as u32;
    let (nx, ny, nz) = (
        counts(b.min.x, b.max.x),
        counts(b.min.y, b.max.y),
        counts(b.min.z, b.max.z),
    );
    let mut sure = 0u64;
    let mut band = 0u64;
    for i in 0..nx {
        if cx.checkpoint().is_err() {
            return Err(fs_query::QueryError::Cancelled);
        }
        for j in 0..ny {
            for k in 0..nz {
                let p = Point3::new(
                    b.min.x + (f64::from(i) + 0.5) * h,
                    b.min.y + (f64::from(j) + 0.5) * h,
                    b.min.z + (f64::from(k) + 0.5) * h,
                );
                let sd = chart.eval(p, cx).signed_distance;
                if sd <= -band_half {
                    sure += 1;
                } else if sd < band_half {
                    band += 1;
                }
            }
        }
    }
    Ok(VolumeEnclosure {
        lo: sure as f64 * cell,
        hi: (sure + band) as f64 * cell,
        h,
    })
}

/// Smoothed volume: `Σ h³·σ(−φ/ε)` with the logistic mollifier — the
/// C¹ value whose lever derivative matches the Hadamard shape
/// derivative on fixtures (the battery's validation).
///
/// # Errors
/// [`fs_query::QueryError::Cancelled`].
pub fn volume_smooth(
    chart: &dyn Chart,
    domain: &fs_geom::Aabb,
    h: f64,
    epsilon: f64,
    cx: &Cx<'_>,
) -> Result<f64, fs_query::QueryError> {
    let b = domain;
    let cell = h * h * h;
    let counts = |lo: f64, hi: f64| ((hi - lo) / h).ceil().max(1.0) as u32;
    let (nx, ny, nz) = (
        counts(b.min.x, b.max.x),
        counts(b.min.y, b.max.y),
        counts(b.min.z, b.max.z),
    );
    let mut acc = 0.0;
    for i in 0..nx {
        if cx.checkpoint().is_err() {
            return Err(fs_query::QueryError::Cancelled);
        }
        for j in 0..ny {
            for k in 0..nz {
                let p = Point3::new(
                    b.min.x + (f64::from(i) + 0.5) * h,
                    b.min.y + (f64::from(j) + 0.5) * h,
                    b.min.z + (f64::from(k) + 0.5) * h,
                );
                let sd = chart.eval(p, cx).signed_distance;
                acc += cell / (1.0 + (sd / epsilon).exp());
            }
        }
    }
    Ok(acc)
}
