//! Simplification-equivalence comparison on the ladder seam (IRP M5.5-W2,
//! plan §1.4). It knows nothing about the physics being compared.
//!
//! A detailed and a simplified frozen description are compared on raw
//! calibrated output (complex, in one declared unit) over a predeclared
//! domain. Each side's numerical uncertainty comes from its own M5.5-W1
//! [`LadderReport`] ([`BoundedSolution::from_ladder`]): a GCI where the fit
//! is asymptotic, the flatness threshold where the rungs are flat, a declared
//! bound, or else the rung spread, with every point's [`BoundSource`] kept and
//! the spread-bounded points counted in the report. An observable with no
//! bound at all refuses.
//!
//! **Point bounds.** Each side lies within its numerical bound, so the metric
//! of the difference at point `i` lies in `[lower_i, upper_i]`:
//! - [`Metric::Absolute`]: `|a - b| +- (bound_a + bound_b)`, floored at 0.
//! - [`Metric::LevelDb`]: `|20 log10(|a|/|b|)|` over the moduli intervals;
//!   a modulus interval reaching zero makes `upper_i` infinite.
//!
//! **Sampling term.** The domain is sampled on at least three nested rungs
//! and declared as lines: ordered points along one refined coordinate (for
//! example frequency, one line per hole state). Between samples the upper
//! bound is assumed C2, and the finest rungs are assumed to resolve its
//! narrowest feature. At every sampled local maximum of a line, a concave
//! three-point stencil with second derivative `u''` cannot hide more than
//! `|u''| H^2 / 8` above the samples, `H` being the larger adjacent interval
//! (a parabola's vertex lies within `H/2` of a sample). On rung `k`,
//! `T_k = max(u_j + |u''_j| H_j^2/8) - U_k` over those maxima, with
//! `U_k = max_{i in rung k} upper_i`. The reported term is `T_K`, and it is
//! infinite unless `T_K <= T_{K-1} <= T_{K-2}`: a term that grows under
//! refinement shows an unresolved feature.
//!
//! **Verdict.** `NotEquivalent` when some sampled `lower_i > tolerance`;
//! otherwise `Equivalent` when `U_K + T <= tolerance`; otherwise
//! `Inconclusive`. Both refusals retain the detail (§1.4 rules 4-5). A GCI is
//! an estimate with a safety factor, not a certificate, and the sampling term
//! is a declared extrapolation; neither is a rigorous enclosure.

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, hash_domain};
use fs_math::det;

use super::{AxisReport, AxisVerdict, Coordinate, EnsembleDecl, LadderReport, dimension};
use crate::convergence::ROUNDOFF_ULPS;

/// One solution's raw output and per-point numerical bound.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedSolution {
    /// Content identity of the frozen description that produced the values.
    pub description: ContentHash,
    /// Identity of the M5.5-W1 report the bounds come from.
    pub ladder: ContentHash,
    /// Unit of the values and bounds.
    pub unit: String,
    /// Complex value `(re, im)` per domain point.
    pub values: Vec<(f64, f64)>,
    /// Absolute numerical bound per point, in `unit`.
    pub bounds: Vec<f64>,
    /// Weakest source behind each point's bound.
    pub sources: Vec<BoundSource>,
}

/// Where a numerical bound comes from, strongest first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum BoundSource {
    /// `Unresolved` rungs (flat to roundoff): the flatness threshold
    /// `ROUNDOFF_ULPS eps max|y|`.
    Roundoff,
    /// An `Asymptotic` fit's GCI times the finest value: an estimate with
    /// the safety factor Fs = 1.25, not a certificate.
    Gci,
    /// A bound the producer declared on a non-algebraic axis.
    Declared,
    /// No GCI (pre-asymptotic, scattered or oscillatory rungs): the spread
    /// of all admitted rungs. It bounds the finest rung's error only when
    /// the rungs contract monotonically or vary by roundoff, so reports
    /// count the points that rely on it.
    Spread,
}

/// Bound and source of every observable of one axis, in its report order.
/// `ObservedOnly`, `InsufficientRungs` and `Refused` have none.
fn axis_bounds(axis: &AxisReport) -> Vec<Option<(f64, BoundSource)>> {
    let admitted: Vec<(f64, &super::RungRecord)> = axis
        .rungs
        .iter()
        .enumerate()
        .filter(|(i, _)| !axis.dead_rungs.contains(i) && !axis.failed_rungs.contains(i))
        .filter_map(|(_, r)| {
            let d = dimension(r, &axis.axis)?;
            let h = match axis.axis.coordinate {
                Coordinate::Direct => d.value,
                Coordinate::Reciprocal => 1.0 / d.value,
            };
            Some((h, r))
        })
        .collect();
    let finest = admitted
        .iter()
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, r)| *r);
    axis.qois
        .iter()
        .enumerate()
        .map(|(j, q)| {
            let fine = finest.and_then(|r| r.qoi.get(j)).map(|o| o.value.abs());
            let largest = admitted
                .iter()
                .filter_map(|(_, r)| r.qoi.get(j))
                .map(|o| o.value.abs())
                .fold(0.0, f64::max);
            let (bound, source) = match q.verdict {
                AxisVerdict::Asymptotic { gci: Some(g), .. } => (g * fine?, BoundSource::Gci),
                AxisVerdict::Unresolved => (
                    ROUNDOFF_ULPS * f64::EPSILON * largest,
                    BoundSource::Roundoff,
                ),
                AxisVerdict::Bounded { bound } => (bound, BoundSource::Declared),
                AxisVerdict::Asymptotic { gci: None, .. }
                | AxisVerdict::PreAsymptotic { .. }
                | AxisVerdict::Scattered
                | AxisVerdict::Oscillatory => (q.spread, BoundSource::Spread),
                AxisVerdict::ObservedOnly { .. }
                | AxisVerdict::InsufficientRungs
                | AxisVerdict::Refused { .. } => return None,
            };
            (bound.is_finite() && bound >= 0.0).then_some((bound, source))
        })
        .collect()
}

/// Every observable's bound in one pass: the sum of its per-axis bounds
/// plus its largest measured joint interaction, with the weakest source.
/// `None` when an axis lacks the observable or a bound for it.
fn bound_table(report: &LadderReport) -> BTreeMap<&str, Option<(f64, BoundSource)>> {
    let mut table: BTreeMap<&str, (usize, Option<(f64, BoundSource)>)> = BTreeMap::new();
    for axis in &report.axes {
        for (q, b) in axis.qois.iter().zip(axis_bounds(axis)) {
            let entry = table
                .entry(q.name.as_str())
                .or_insert((0, Some((0.0, BoundSource::Roundoff))));
            entry.0 += 1;
            entry.1 = match (entry.1, b) {
                (Some((total, worst)), Some((v, source))) => Some((total + v, worst.max(source))),
                _ => None,
            };
        }
    }
    let mut interaction: BTreeMap<&str, f64> = BTreeMap::new();
    for x in &report.joint.interactions {
        let e = interaction
            .entry(x.joint_error.name.as_str())
            .or_insert(0.0);
        *e = e.max(x.interaction.abs());
    }
    table
        .into_iter()
        .map(|(name, (count, b))| {
            let extra = interaction.get(name).copied().unwrap_or(0.0);
            let b = b
                .filter(|_| count == report.axes.len())
                .map(|(v, source)| (v + extra, source))
                .filter(|(v, _)| v.is_finite());
            (name, b)
        })
        .collect()
}

/// Numerical bound and its weakest source for the observable `name` in a
/// ladder report: the sum of its per-axis bounds plus the largest measured
/// joint interaction. `None` when any axis lacks a bound for it.
#[must_use]
pub fn qoi_bound(report: &LadderReport, name: &str) -> Option<(f64, BoundSource)> {
    bound_table(report).get(name).copied().flatten()
}

impl BoundedSolution {
    /// Values from the finest joint rung of a solution's own M5.5-W1 report
    /// and bounds from [`qoi_bound`]. Point `i` is the observable pair
    /// `points[i] = (re, im)`, both in `unit`; its bound is the modulus of the
    /// two observable bounds and its source the weaker of theirs.
    ///
    /// # Errors
    /// A missing, duplicated or mis-united observable, or one with no bound.
    pub fn from_ladder(
        description: ContentHash,
        report: &LadderReport,
        unit: &str,
        points: &[(String, String)],
    ) -> Result<Self, String> {
        let finest = report
            .joint
            .rungs
            .last()
            .ok_or("ladder has no joint rungs")?;
        let mut observed: BTreeMap<&str, Vec<&super::Observable>> = BTreeMap::new();
        for q in &finest.qoi {
            observed.entry(q.name.as_str()).or_default().push(q);
        }
        let table = bound_table(report);
        let value = |name: &str| -> Result<(f64, f64, BoundSource), String> {
            let found = observed
                .get(name)
                .ok_or_else(|| format!("observable {name} is absent"))?;
            if found.len() != 1 || found[0].unit != unit {
                return Err(format!("observable {name} is duplicated or not in {unit}"));
            }
            let (bound, source) = table
                .get(name)
                .copied()
                .flatten()
                .ok_or_else(|| format!("observable {name} has no ladder bound"))?;
            Ok((found[0].value, bound, source))
        };
        let mut solution = Self {
            description,
            ladder: report.hash,
            unit: unit.to_owned(),
            values: Vec::with_capacity(points.len()),
            bounds: Vec::with_capacity(points.len()),
            sources: Vec::with_capacity(points.len()),
        };
        for (re, im) in points {
            let ((vr, br, sr), (vi, bi, si)) = (value(re)?, value(im)?);
            solution.values.push((vr, vi));
            solution.bounds.push(det::sqrt(br * br + bi * bi));
            solution.sources.push(sr.max(si));
        }
        Ok(solution)
    }
}

/// How the difference of two complex outputs is measured.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Metric {
    /// `|a - b|`, in the values' unit.
    Absolute,
    /// `|20 log10(|a|/|b|)|` in dB: a level difference of raw outputs.
    LevelDb,
}

/// A predeclared comparison domain, fixed before either side is evaluated.
#[derive(Debug, Clone, PartialEq)]
pub struct EquivalenceDomain {
    /// Coverage: band, states, receivers and, for played cases, the notes,
    /// dynamics, attacks and releases.
    pub description: String,
    /// Nested sampling rungs, coarse to fine, as indices into the points.
    /// Each rung strictly contains the previous; the last holds every point.
    pub rungs: Vec<Vec<usize>>,
    /// Ordered point sequences along the refined coordinate, one per case;
    /// together they hold every point once, and every rung keeps at least
    /// three points of each line.
    pub lines: Vec<Vec<usize>>,
    /// The refined coordinate of each point, strictly increasing along a line.
    pub coordinate: Vec<f64>,
    /// A stochastic source needs a distributional comparison: seeds refuse.
    pub ensemble: EnsembleDecl,
    /// Model-form caveats, carried into the report verbatim.
    pub limitations: Vec<String>,
}

/// Three-way verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum EquivalenceVerdict {
    /// The upper bound, sampling term included, is within tolerance.
    Equivalent {
        /// `U_K + T`.
        max_upper: f64,
    },
    /// A sampled lower bound exceeds the tolerance.
    NotEquivalent {
        /// Largest sampled lower bound.
        max_lower: f64,
        /// Its point index.
        at: usize,
    },
    /// The bounds straddle the tolerance.
    Inconclusive {
        /// Points whose `[lower, upper + T]` contains the tolerance.
        straddling: Vec<usize>,
    },
}

/// Retained comparison, with every per-point bound.
#[derive(Debug, Clone, PartialEq)]
pub struct EquivalenceReport {
    /// The three-way verdict.
    pub verdict: EquivalenceVerdict,
    /// Metric of the difference.
    pub metric: Metric,
    /// Equivalence tolerance in the metric's unit.
    pub tolerance: f64,
    /// Upper bound per point.
    pub upper: Vec<f64>,
    /// Lower bound per point.
    pub lower: Vec<f64>,
    /// Sampled maximum `U_k` per rung.
    pub sampled_maxima: Vec<f64>,
    /// Curvature sampling term `T_k` per rung.
    pub sampling_terms: Vec<f64>,
    /// Reported sampling term: `T_K`, or infinite when the last three
    /// rungs' terms grow.
    pub sampling_term: f64,
    /// Model-form caveats from the domain.
    pub limitations: Vec<String>,
    /// Detailed and simplified description identities.
    pub descriptions: [ContentHash; 2],
    /// Detailed and simplified M5.5-W1 report identities.
    pub ladders: [ContentHash; 2],
    /// Points per side whose bound relies on a rung spread (no GCI).
    pub spread_bounded: [usize; 2],
    /// Content identity of everything above.
    pub hash: ContentHash,
}

impl EquivalenceReport {
    /// Only an `Equivalent` verdict with no model-form caveat authorizes
    /// omitting the detail; any limitation leaves the claim to its consumer.
    #[must_use]
    pub fn authorizes_omission(&self) -> bool {
        matches!(self.verdict, EquivalenceVerdict::Equivalent { .. }) && self.limitations.is_empty()
    }
}

fn modulus(v: (f64, f64)) -> f64 {
    det::sqrt(v.0 * v.0 + v.1 * v.1)
}

fn db(ratio: f64) -> f64 {
    20.0 * det::ln(ratio) / det::ln(10.0)
}

/// `(lower, upper)` of the metric of `a - b` when each side lies within
/// its bound.
fn point_bounds(metric: Metric, a: (f64, f64), ba: f64, b: (f64, f64), bb: f64) -> (f64, f64) {
    match metric {
        Metric::Absolute => {
            let d = modulus((a.0 - b.0, a.1 - b.1));
            ((d - ba - bb).max(0.0), d + ba + bb)
        }
        Metric::LevelDb => {
            let (ma, mb) = (modulus(a), modulus(b));
            let (a_lo, a_hi, b_lo, b_hi) = (ma - ba, ma + ba, mb - bb, mb + bb);
            if a_lo <= 0.0 || b_lo <= 0.0 {
                let lower = if a_lo > b_hi {
                    db(a_lo / b_hi)
                } else if b_lo > a_hi {
                    db(b_lo / a_hi)
                } else {
                    0.0
                };
                return (lower.max(0.0), f64::INFINITY);
            }
            let (lo, hi) = (db(a_lo / b_hi), db(a_hi / b_lo));
            let lower = if lo > 0.0 {
                lo
            } else if hi < 0.0 {
                -hi
            } else {
                0.0
            };
            (lower, lo.abs().max(hi.abs()))
        }
    }
}

/// Largest `u_j + |u''_j| H_j^2 / 8` over the sampled local maxima of one
/// line (`pts` in coordinate order, at least three), from the three-point
/// stencil nearest each maximum; non-finite when any value is.
fn line_reach(pts: &[usize], u: &[f64], x: &[f64]) -> f64 {
    let m = pts.len();
    let mut reach = f64::NEG_INFINITY;
    for j in 0..m {
        let v = u[pts[j]];
        if !v.is_finite() {
            return f64::INFINITY;
        }
        let left = j == 0 || v >= u[pts[j - 1]];
        let right = j + 1 == m || v >= u[pts[j + 1]];
        if !(left && right) {
            continue;
        }
        let b = j.clamp(1, m - 2);
        let (xa, xb, xc) = (x[pts[b - 1]], x[pts[b]], x[pts[b + 1]]);
        let (ua, ub, uc) = (u[pts[b - 1]], u[pts[b]], u[pts[b + 1]]);
        let second = 2.0 * ((uc - ub) / (xc - xb) - (ub - ua) / (xb - xa)) / (xc - xa);
        let gap = |a: usize, c: usize| x[pts[c]] - x[pts[a]];
        let h = match j {
            0 => gap(0, 1),
            _ if j + 1 == m => gap(j - 1, j),
            _ => gap(j - 1, j).max(gap(j, j + 1)),
        };
        reach = reach.max(v + (-second).max(0.0) * h * h / 8.0);
    }
    reach
}

fn valid(s: &BoundedSolution, n: usize) -> bool {
    !s.unit.is_empty()
        && s.values.len() == n
        && s.bounds.len() == n
        && s.values.iter().all(|v| v.0.is_finite() && v.1.is_finite())
        && s.bounds.iter().all(|b| b.is_finite() && *b >= 0.0)
        && s.sources.len() == n
}

/// Compare a detailed and a simplified solution over a predeclared domain.
///
/// # Errors
/// Mismatched or non-finite solutions, a missing bound, fewer than three
/// nested sampling rungs, a last rung that omits a point, lines that do not
/// cover every point once in increasing coordinate or keep fewer than three
/// points on a rung, a declared seed ensemble (pointwise comparison
/// refused), or an invalid tolerance.
pub fn compare_equivalence(
    detailed: &BoundedSolution,
    simplified: &BoundedSolution,
    domain: &EquivalenceDomain,
    metric: Metric,
    tolerance: f64,
) -> Result<EquivalenceReport, String> {
    let n = detailed.values.len();
    if n == 0 || !valid(detailed, n) || !valid(simplified, n) || detailed.unit != simplified.unit {
        return Err("both solutions need one finite value and bound per point, in one unit".into());
    }
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err("equivalence tolerance must be positive and finite".into());
    }
    match &domain.ensemble {
        EnsembleDecl::Deterministic { reason } if !reason.is_empty() => {}
        EnsembleDecl::Deterministic { .. } => {
            return Err("deterministic ensemble declaration needs a reason".into());
        }
        EnsembleDecl::Seeds { .. } => {
            return Err(
                "a seed ensemble needs a distributional comparison, not a pointwise one".into(),
            );
        }
    }
    if domain.description.is_empty() || domain.rungs.len() < 3 {
        return Err("the domain needs a description and at least three sampling rungs".into());
    }
    for (k, rung) in domain.rungs.iter().enumerate() {
        let distinct = rung.iter().enumerate().all(|(i, p)| !rung[..i].contains(p));
        let nested = k == 0
            || (rung.len() > domain.rungs[k - 1].len()
                && domain.rungs[k - 1].iter().all(|p| rung.contains(p)));
        if rung.is_empty() || !distinct || !nested || rung.iter().any(|&p| p >= n) {
            return Err("sampling rungs must be distinct, in range and strictly nested".into());
        }
    }
    if domain.rungs[domain.rungs.len() - 1].len() != n {
        return Err("the finest sampling rung must hold every domain point".into());
    }
    let x = &domain.coordinate;
    if x.len() != n || !x.iter().all(|v| v.is_finite()) {
        return Err("one finite coordinate per point".into());
    }
    let mut seen = vec![false; n];
    for line in &domain.lines {
        for (j, &p) in line.iter().enumerate() {
            if p >= n || seen[p] || (j > 0 && !(x[p] > x[line[j - 1]])) {
                return Err("lines must hold each point once, in increasing coordinate".into());
            }
            seen[p] = true;
        }
    }
    if seen.contains(&false) {
        return Err("the lines must cover every point".into());
    }
    let mut in_rung = vec![vec![false; n]; domain.rungs.len()];
    for (k, rung) in domain.rungs.iter().enumerate() {
        for &p in rung {
            in_rung[k][p] = true;
        }
        if domain
            .lines
            .iter()
            .any(|line| line.iter().filter(|&&p| in_rung[k][p]).count() < 3)
        {
            return Err("every rung must keep at least three points of each line".into());
        }
    }
    let (mut lower, mut upper) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for i in 0..n {
        let (lo, hi) = point_bounds(
            metric,
            detailed.values[i],
            detailed.bounds[i],
            simplified.values[i],
            simplified.bounds[i],
        );
        lower.push(lo);
        upper.push(hi);
    }
    let sampled_maxima: Vec<f64> = domain
        .rungs
        .iter()
        .map(|rung| rung.iter().map(|&i| upper[i]).fold(0.0, f64::max))
        .collect();
    let sampling_terms: Vec<f64> = in_rung
        .iter()
        .zip(&sampled_maxima)
        .map(|(member, &sampled)| {
            let mut reach = sampled;
            for line in &domain.lines {
                let pts: Vec<usize> = line.iter().copied().filter(|&p| member[p]).collect();
                reach = reach.max(line_reach(&pts, &upper, x));
            }
            if reach.is_finite() {
                reach - sampled
            } else {
                f64::INFINITY
            }
        })
        .collect();
    let k = sampled_maxima.len() - 1;
    let settled = sampling_terms[k] <= sampling_terms[k - 1]
        && sampling_terms[k - 1] <= sampling_terms[k - 2];
    let sampling_term = if settled {
        sampling_terms[k]
    } else {
        f64::INFINITY
    };
    let (mut at, mut max_lower) = (0, f64::NEG_INFINITY);
    for (i, &v) in lower.iter().enumerate() {
        if v > max_lower {
            (at, max_lower) = (i, v);
        }
    }
    let max_upper = sampled_maxima[k] + sampling_term;
    let verdict = if max_lower > tolerance {
        EquivalenceVerdict::NotEquivalent { max_lower, at }
    } else if max_upper <= tolerance {
        EquivalenceVerdict::Equivalent { max_upper }
    } else {
        EquivalenceVerdict::Inconclusive {
            straddling: (0..n)
                .filter(|&i| lower[i] <= tolerance && upper[i] + sampling_term > tolerance)
                .collect(),
        }
    };
    let mut report = EquivalenceReport {
        verdict,
        metric,
        tolerance,
        upper,
        lower,
        sampled_maxima,
        sampling_terms,
        sampling_term,
        limitations: domain.limitations.clone(),
        descriptions: [detailed.description, simplified.description],
        ladders: [detailed.ladder, simplified.ladder],
        spread_bounded: [detailed, simplified].map(|s| {
            s.sources
                .iter()
                .filter(|&&b| b == BoundSource::Spread)
                .count()
        }),
        hash: ContentHash([0; 32]),
    };
    let bytes = format!(
        "{:?}{:?}{:?}{:?}{report:?}",
        domain.description, domain.rungs, domain.lines, domain.coordinate
    );
    report.hash = hash_domain("org.frankensim.equivalence.v1", bytes.as_bytes());
    Ok(report)
}
