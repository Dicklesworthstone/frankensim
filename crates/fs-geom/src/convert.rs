//! Certified conversion (plan §7.3's edges, founded here): a functor with
//! a RECEIPT. `Convert<Dst>` produces a [`Certified`] destination chart
//! whose evidence carries the achieved error bound, composable through the
//! Error Ledger — or a structured refusal when the budget is infeasible
//! (admission fails early instead of running a doomed conversion).

use crate::{Aabb, Chart, ChartSample, Differentiability, Point3};
use fs_evidence::{Certified, Evidence, NumericalCertificate, ProvenanceHash};
use fs_exec::Cx;
use core::fmt;

/// The conversion error budget (v1: absolute signed-distance error; the
/// full cost×error Pareto machinery is the Rep Router bead's).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErrBudget {
    /// Maximum tolerated absolute signed-distance error.
    pub abs_sd_error: f64,
}

/// Structured conversion refusal (Decalogue P10: a refusal that teaches).
#[derive(Debug, Clone, PartialEq)]
pub enum ConvertDiag {
    /// The budget cannot be met within this converter's resource cap.
    BudgetInfeasible {
        /// The requested absolute error.
        requested: f64,
        /// The best this converter can achieve at its resolution cap.
        achievable: f64,
        /// Grid resolution the request would need.
        need_resolution: u32,
        /// The converter's per-axis resolution cap.
        cap: u32,
    },
    /// The source declared no (finite) Lipschitz bound, so no rigorous
    /// sampled enclosure exists.
    NoLipschitzBound,
}

impl fmt::Display for ConvertDiag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConvertDiag::BudgetInfeasible {
                requested,
                achievable,
                need_resolution,
                cap,
            } => write!(
                f,
                "conversion refused before running: abs_sd_error {requested} needs a \
                 {need_resolution}^3 grid but the cap is {cap}^3 (achievable: {achievable}). \
                 Fixes (ranked): (1) relax the budget to {achievable}; (2) shrink the region \
                 of interest; (3) wait for the adaptive-sampling chart (rep-sdf bead)"
            ),
            ConvertDiag::NoLipschitzBound => write!(
                f,
                "conversion refused: the source chart claims no certified Lipschitz bound, so \
                 a sampled enclosure would be a guess; use a chart that certifies one"
            ),
        }
    }
}

impl core::error::Error for ConvertDiag {}

/// Certified conversion between representations (plan Appendix B).
pub trait Convert<Dst: Chart>: Chart + Sized {
    /// Convert under `budget`, returning the destination chart WITH its
    /// receipt (achieved bound, provenance chain) — or refuse early.
    ///
    /// # Errors
    /// [`ConvertDiag`] with ranked fixes.
    fn convert(&self, budget: ErrBudget, cx: &Cx<'_>) -> Result<Certified<Dst>, ConvertDiag>;
}

/// A dense sampled-SDF chart with trilinear interpolation: the first
/// concrete conversion TARGET (FrankenVDB-class sparse charts are the
/// rep-sdf bead's; this one exists so conversion receipts are testable
/// end-to-end).
#[derive(Debug, Clone)]
pub struct SampledSdf {
    box_: Aabb,
    n: u32,
    values: Vec<f64>,
    /// Rigorous |sampled - source| bound inside the box.
    bound: f64,
    /// The source's certified Lipschitz constant (outside-box enclosures
    /// lean on it).
    source_lipschitz: f64,
}

impl SampledSdf {
    /// Grid resolution per axis.
    #[must_use]
    pub fn resolution(&self) -> u32 {
        self.n
    }

    /// The rigorous in-box error bound this chart was built with.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.bound
    }

    fn idx(&self, i: u32, j: u32, k: u32) -> usize {
        ((k * self.n + j) * self.n + i) as usize
    }

    fn interp(&self, p: Point3) -> f64 {
        let n = f64::from(self.n - 1);
        let fx = ((p.x - self.box_.min.x) / (self.box_.max.x - self.box_.min.x) * n)
            .clamp(0.0, n);
        let fy = ((p.y - self.box_.min.y) / (self.box_.max.y - self.box_.min.y) * n)
            .clamp(0.0, n);
        let fz = ((p.z - self.box_.min.z) / (self.box_.max.z - self.box_.min.z) * n)
            .clamp(0.0, n);
        let (i0, j0, k0) = (fx as u32, fy as u32, fz as u32);
        let (i1, j1, k1) = (
            (i0 + 1).min(self.n - 1),
            (j0 + 1).min(self.n - 1),
            (k0 + 1).min(self.n - 1),
        );
        let (tx, ty, tz) = (fx - f64::from(i0), fy - f64::from(j0), fz - f64::from(k0));
        let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
        let c00 = lerp(self.values[self.idx(i0, j0, k0)], self.values[self.idx(i1, j0, k0)], tx);
        let c10 = lerp(self.values[self.idx(i0, j1, k0)], self.values[self.idx(i1, j1, k0)], tx);
        let c01 = lerp(self.values[self.idx(i0, j0, k1)], self.values[self.idx(i1, j0, k1)], tx);
        let c11 = lerp(self.values[self.idx(i0, j1, k1)], self.values[self.idx(i1, j1, k1)], tx);
        lerp(lerp(c00, c10, ty), lerp(c01, c11, ty), tz)
    }
}

impl Chart for SampledSdf {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let clamped = Point3::new(
            x.x.clamp(self.box_.min.x, self.box_.max.x),
            x.y.clamp(self.box_.min.y, self.box_.max.y),
            x.z.clamp(self.box_.min.z, self.box_.max.z),
        );
        let dist_out = x.delta_from(clamped).norm();
        let base = self.interp(clamped);
        if dist_out == 0.0 {
            ChartSample {
                signed_distance: base,
                gradient: None,
                lipschitz: None, // interpolant Lipschitz certification is rep-sdf's
                error: NumericalCertificate::enclosure(base - self.bound, base + self.bound),
            }
        } else {
            // Outside the sampled box: value = interp(clamp) + distance.
            // Rigorous for an L-Lipschitz source: true sd ∈
            // [v - bound - (1+L)·dist_out, v + bound] (conservative).
            let v = base + dist_out;
            ChartSample {
                signed_distance: v,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::enclosure(
                    v - self.bound - (1.0 + self.source_lipschitz) * dist_out,
                    v + self.bound,
                ),
            }
        }
    }

    fn support(&self) -> Aabb {
        self.box_
    }

    fn name(&self) -> &'static str {
        "sampled-sdf"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

/// Per-axis resolution cap for the dense sampled target (beyond this the
/// dense grid stops being the right tool — the refusal says so).
pub const SAMPLED_SDF_MAX_RESOLUTION: u32 = 96;

/// Blanket conversion: ANY chart with a certified Lipschitz bound and an
/// exact-or-enclosed error model can be sampled into a [`SampledSdf`] with
/// a rigorous receipt. (The plan's mesh→SDF / F-rep→SDF edges specialize
/// this; the sphere fixture exercises it today.)
impl<C: Chart> Convert<SampledSdf> for C {
    fn convert(&self, budget: ErrBudget, cx: &Cx<'_>) -> Result<Certified<SampledSdf>, ConvertDiag> {
        let box_ = self.support().inflate(budget.abs_sd_error.max(1e-9));
        // Probe the source's Lipschitz claim at the box center.
        let center = Point3::new(
            f64::midpoint(box_.min.x, box_.max.x),
            f64::midpoint(box_.min.y, box_.max.y),
            f64::midpoint(box_.min.z, box_.max.z),
        );
        let lipschitz = match self.eval(center, cx).lipschitz {
            Some(l) if l.is_finite() => l,
            _ => return Err(ConvertDiag::NoLipschitzBound),
        };
        let edge = (box_.max.x - box_.min.x)
            .max(box_.max.y - box_.min.y)
            .max(box_.max.z - box_.min.z);
        // Trilinear error of an L-Lipschitz field sampled at step h is
        // ≤ L·h (conservative; the √3/2 constant is left on the table).
        let h_needed = budget.abs_sd_error / lipschitz.max(f64::MIN_POSITIVE);
        let need_resolution = (edge / h_needed).ceil() as u32 + 1;
        if need_resolution > SAMPLED_SDF_MAX_RESOLUTION {
            let h_cap = edge / f64::from(SAMPLED_SDF_MAX_RESOLUTION - 1);
            return Err(ConvertDiag::BudgetInfeasible {
                requested: budget.abs_sd_error,
                achievable: lipschitz * h_cap,
                need_resolution,
                cap: SAMPLED_SDF_MAX_RESOLUTION,
            });
        }
        let n = need_resolution.max(2);
        let h = edge / f64::from(n - 1);
        let achieved = lipschitz * h;
        let mut values = Vec::with_capacity((n as usize).pow(3));
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let p = Point3::new(
                        box_.min.x + f64::from(i) * (box_.max.x - box_.min.x) / f64::from(n - 1),
                        box_.min.y + f64::from(j) * (box_.max.y - box_.min.y) / f64::from(n - 1),
                        box_.min.z + f64::from(k) * (box_.max.z - box_.min.z) / f64::from(n - 1),
                    );
                    values.push(self.eval(p, cx).signed_distance);
                }
            }
        }
        let chart = SampledSdf {
            box_,
            n,
            values,
            bound: achieved,
            source_lipschitz: lipschitz,
        };
        // The receipt: QoI = achieved bound, enclosed rigorously; the
        // provenance chains source-name → conversion.
        let receipt = Evidence {
            qoi: achieved,
            numerical: NumericalCertificate::enclosure(0.0, achieved),
            statistical: fs_evidence::StatisticalCertificate::None,
            model: fs_evidence::ModelEvidence::none(),
            sensitivity: fs_evidence::SensitivitySummary::default(),
            provenance: ProvenanceHash::chain(
                "convert/sampled-sdf",
                &[ProvenanceHash::of_bytes(self.name().as_bytes())],
            ),
            adjoint_ref: None,
            value: chart,
        };
        receipt.certified().map_err(|e| {
            // Defensive: receipts here are enclosure-grade by construction.
            debug_assert!(false, "conversion receipt must certify: {e}");
            ConvertDiag::NoLipschitzBound
        })
    }
}
