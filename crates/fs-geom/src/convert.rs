//! Evidence-bearing conversion (plan §7.3's edges, founded here): a functor
//! with a RECEIPT. `Convert<Dst>` always exposes the achieved authority and
//! error bound through [`Evidence`]; only a rigorous unclipped global result
//! may promote to [`Certified`]. Weak-source and explicitly clipped local
//! conversions remain plain evidence rather than laundering scope or theorem
//! strength. Infeasible budgets return a structured early refusal.

use crate::{
    Aabb, Axis, Chart, ChartSample, ClippedChart, Differentiability, Point3, SamplingDomain,
    SamplingDomainError, TraceStepClaim,
};
use core::fmt;
use fs_evidence::{Certified, Evidence, NumericalCertificate, NumericalKind, ProvenanceHash};
use fs_exec::Cx;
use fs_ivl::Interval;

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
    /// Cooperative cancellation was observed before publishing a chart or
    /// receipt.
    Cancelled {
        /// Stable conversion stage at which cancellation was observed.
        stage: &'static str,
        /// Grid samples fully validated before cancellation.
        completed_samples: usize,
    },
    /// The requested absolute error must be finite and strictly positive.
    InvalidBudget {
        /// Exact IEEE-754 bits supplied by the caller.
        requested_bits: u64,
    },
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
    /// The source returned a non-finite signed-distance value, so no finite
    /// sampled chart or error receipt can be certified.
    NonFiniteSignedDistance {
        /// Point at which the malformed value was observed.
        point: Point3,
        /// Exact IEEE-754 bits returned by the source chart.
        value_bits: u64,
    },
    /// The requested per-axis grid has more nodes than the admitted floating-
    /// point interval can represent distinctly. Sampling coincident nodes
    /// would invalidate both cell lookup and the reconstruction bound.
    UnrepresentableGrid {
        /// First axis whose nodes could not be made strictly increasing.
        axis: Axis,
        /// First node index that rounded onto a predecessor or endpoint.
        index: u32,
        /// Requested nodes per axis.
        resolution: u32,
        /// Exact admitted lower-endpoint bits.
        min_bits: u64,
        /// Exact admitted upper-endpoint bits.
        max_bits: u64,
    },
    /// The source support and optional clip could not be admitted as a finite
    /// sampling domain.
    SamplingDomain(SamplingDomainError),
    /// Sampling completed, but the source did not carry rigorous authority
    /// relative to the abstract region signed distance. The plain
    /// [`Evidence`] remains available through [`Convert::convert_with_domain`]
    /// or [`Convert::convert_clipped`].
    NoAbstractDistanceClaim {
        /// Weakest source authority observed by the sampler.
        kind: NumericalKind,
    },
}

impl fmt::Display for ConvertDiag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConvertDiag::Cancelled {
                stage,
                completed_samples,
            } => write!(
                f,
                "conversion cancelled during {stage} after {completed_samples} validated grid samples; no chart or receipt was published"
            ),
            ConvertDiag::InvalidBudget { requested_bits } => write!(
                f,
                "conversion refused before sampling: abs_sd_error bits {requested_bits:016x} are not finite and strictly positive; provide a positive finite error budget"
            ),
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
            ConvertDiag::NonFiniteSignedDistance { point, value_bits } => write!(
                f,
                "conversion refused: the source chart returned non-finite signed-distance bits \
                 {value_bits:016x} at ({}, {}, {}); no sampled chart or certified receipt was \
                 published",
                point.x, point.y, point.z
            ),
            ConvertDiag::UnrepresentableGrid {
                axis,
                index,
                resolution,
                ..
            } => write!(
                f,
                "conversion refused before sampling: node {index} of the requested \
                 {resolution}-node grid is not distinctly representable on the {} axis; \
                 widen the finite domain, relax the error budget, or use an adaptive \
                 representation",
                axis.name()
            ),
            ConvertDiag::SamplingDomain(error) => write!(f, "conversion {error}"),
            ConvertDiag::NoAbstractDistanceClaim { kind } => write!(
                f,
                "conversion sampled the source field, but its abstract signed-distance \
                 authority is {kind:?}; keep the result as plain Evidence or provide a \
                 rigorous source chart before requesting Certified output"
            ),
        }
    }
}

impl core::error::Error for ConvertDiag {}

/// Evidence-bearing conversion between representations (plan Appendix B).
/// Global conversion promotes rigorous evidence to [`Certified`]; explicitly
/// clipped conversion remains plain [`Evidence`] because a generic clipped
/// implicit field does not inherit an abstract signed-distance theorem.
pub trait Convert<Dst: Chart + 'static>: Chart + Sized {
    /// Convert over either the source support or its geometric intersection
    /// with `explicit_clip`.
    fn convert_with_domain(
        &self,
        budget: ErrBudget,
        explicit_clip: Option<Aabb>,
        cx: &Cx<'_>,
    ) -> Result<Evidence<Dst>, ConvertDiag>;

    /// Convert under `budget`, returning the destination chart WITH its
    /// receipt (achieved bound, provenance chain) — or refuse early.
    ///
    /// # Errors
    /// [`ConvertDiag`] with ranked fixes.
    fn convert(&self, budget: ErrBudget, cx: &Cx<'_>) -> Result<Certified<Dst>, ConvertDiag> {
        let evidence = self.convert_with_domain(budget, None, cx)?;
        let kind = evidence.numerical.kind;
        evidence
            .certified()
            .map_err(|_| ConvertDiag::NoAbstractDistanceClaim { kind })
    }

    /// Convert the geometric intersection of this chart with an explicit
    /// finite clip AABB.
    ///
    /// # Errors
    /// [`ConvertDiag`] with ranked fixes or a structured sampling-domain
    /// refusal.
    fn convert_clipped(
        &self,
        budget: ErrBudget,
        clip: Aabb,
        cx: &Cx<'_>,
    ) -> Result<Evidence<Dst>, ConvertDiag> {
        self.convert_with_domain(budget, Some(clip), cx)
    }
}

/// A dense sampled-SDF chart with trilinear interpolation: the first
/// concrete conversion TARGET (FrankenVDB-class sparse charts are the
/// rep-sdf bead's; this one exists so conversion receipts are testable
/// end-to-end).
#[derive(Debug, Clone)]
pub struct SampledSdf {
    box_: Aabb,
    n: u32,
    /// Actual representable coordinates used for sampling and cell lookup.
    /// Every axis is strictly increasing and includes the exact box endpoints.
    axis_nodes: [Vec<f64>; 3],
    values: Vec<f64>,
    /// Reconstruction radius relative to the sampled source field. See
    /// [`SampledSdf::nominal_field_bound`] for exactly what it bounds and
    /// [`SampledSdf::nominal_field_bound_kind`] for its authority.
    nominal_field_bound: f64,
    /// Authority of `nominal_field_bound`: `Enclosure` when the source
    /// declared [`TraceStepClaim::ExactDistance`] (whose unit-Lipschitz
    /// theorem is GLOBAL), `Estimate` otherwise — a maximum of per-node
    /// LOCAL Lipschitz values does not bound the slope across a cell.
    nominal_field_bound_kind: NumericalKind,
    /// Weakest source authority observed across the grid, after interpolation
    /// demotes `Exact` to `Enclosure`.
    abstract_distance_kind: NumericalKind,
    /// Total reconstruction-plus-source bound when the source supplied a
    /// finite enclosure or estimate. `None` means honest `NoClaim`.
    abstract_distance_bound: Option<f64>,
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

    /// Actual representable coordinates retained on one sampling axis.
    /// The slice includes both support endpoints and is strictly increasing.
    #[must_use]
    pub fn axis_nodes(&self, axis: Axis) -> &[f64] {
        match axis {
            Axis::X => &self.axis_nodes[0],
            Axis::Y => &self.axis_nodes[1],
            Axis::Z => &self.axis_nodes[2],
        }
    }

    /// Alias for [`Self::nominal_field_bound`] — read that doc before
    /// using the number. It is NOT unconditionally a bound, and it is not
    /// abstract-distance authority: use [`Self::abstract_distance_bound`]
    /// and [`Self::abstract_distance_kind`] before making a
    /// region-distance claim.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.nominal_field_bound
    }

    /// Reconstruction radius relative to the sampled source field:
    /// `L·(largest cell diagonal)` plus an outward interpolation-roundoff
    /// allowance, all outward-rounded.
    ///
    /// It is a RIGOROUS bound on `|eval(p) − source(p)|` over the sampling
    /// box only when the source declared
    /// [`TraceStepClaim::ExactDistance`], because that claim carries a
    /// GLOBAL unit-Lipschitz theorem and the trilinear convex-combination
    /// argument needs a Lipschitz bound valid across the whole cell.
    ///
    /// For every weaker source `L` is the maximum of the per-node
    /// [`ChartSample::lipschitz`] values, which are certified LOCAL to
    /// their query point with no stated radius. A slope spike strictly
    /// between grid nodes is unresolved and CAN exceed this radius, so the
    /// number is an ESTIMATE there. [`Self::nominal_field_bound_kind`]
    /// reports which of the two it is; do not read the bare `f64` out of
    /// that context.
    #[must_use]
    pub fn nominal_field_bound(&self) -> f64 {
        self.nominal_field_bound
    }

    /// Authority of [`Self::nominal_field_bound`]:
    /// [`NumericalKind::Enclosure`] when the source declared
    /// [`TraceStepClaim::ExactDistance`], [`NumericalKind::Estimate`]
    /// otherwise.
    #[must_use]
    pub fn nominal_field_bound_kind(&self) -> NumericalKind {
        self.nominal_field_bound_kind
    }

    /// Weakest abstract signed-distance authority carried by this chart.
    #[must_use]
    pub fn abstract_distance_kind(&self) -> NumericalKind {
        self.abstract_distance_kind
    }

    /// Total error bound relative to abstract region signed distance when the
    /// source supplied a finite enclosure or estimate.
    #[must_use]
    pub fn abstract_distance_bound(&self) -> Option<f64> {
        self.abstract_distance_bound
    }

    fn idx(&self, i: u32, j: u32, k: u32) -> usize {
        ((k * self.n + j) * self.n + i) as usize
    }

    fn interp(&self, p: Point3) -> f64 {
        let (i0, i1, tx) = locate_axis(&self.axis_nodes[0], p.x);
        let (j0, j1, ty) = locate_axis(&self.axis_nodes[1], p.y);
        let (k0, k1, tz) = locate_axis(&self.axis_nodes[2], p.z);
        let (i0, i1, j0, j1, k0, k1) = (
            i0 as u32, i1 as u32, j0 as u32, j1 as u32, k0 as u32, k1 as u32,
        );
        let c00 = stable_lerp(
            self.values[self.idx(i0, j0, k0)],
            self.values[self.idx(i1, j0, k0)],
            tx,
        );
        let c10 = stable_lerp(
            self.values[self.idx(i0, j1, k0)],
            self.values[self.idx(i1, j1, k0)],
            tx,
        );
        let c01 = stable_lerp(
            self.values[self.idx(i0, j0, k1)],
            self.values[self.idx(i1, j0, k1)],
            tx,
        );
        let c11 = stable_lerp(
            self.values[self.idx(i0, j1, k1)],
            self.values[self.idx(i1, j1, k1)],
            tx,
        );
        stable_lerp(stable_lerp(c00, c10, ty), stable_lerp(c01, c11, ty), tz)
    }
}

impl Chart for SampledSdf {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        if !x.x.is_finite() || !x.y.is_finite() || !x.z.is_finite() {
            return ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            };
        }
        let clamped = Point3::new(
            x.x.clamp(self.box_.min.x, self.box_.max.x),
            x.y.clamp(self.box_.min.y, self.box_.max.y),
            x.z.clamp(self.box_.min.z, self.box_.max.z),
        );
        let base = self.interp(clamped);
        let delta = x.delta_from(clamped);
        let Some((dist_out, dist_out_hi)) = finite_norm_with_upper(delta.x, delta.y, delta.z)
        else {
            return ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            };
        };
        let abstract_error = |interval: Interval| {
            if !interval.lo().is_finite() || !interval.hi().is_finite() {
                return NumericalCertificate::no_claim();
            }
            match self.abstract_distance_kind {
                NumericalKind::Exact | NumericalKind::Enclosure => {
                    NumericalCertificate::enclosure(interval.lo(), interval.hi())
                }
                NumericalKind::Estimate => {
                    NumericalCertificate::estimate(interval.lo(), interval.hi())
                }
                NumericalKind::NoClaim => NumericalCertificate::no_claim(),
            }
        };
        let bound = self.abstract_distance_bound.unwrap_or(0.0);
        let at_clamp = centered_interval(base, bound);
        if dist_out == 0.0 {
            ChartSample {
                signed_distance: base,
                gradient: None,
                lipschitz: None, // interpolant Lipschitz certification is rep-sdf's
                error: abstract_error(at_clamp),
            }
        } else {
            // Outside the sampled box, retain the historical positive extension
            // as the nominal value. Authority is attached to the independently
            // derived source interval: f(x) lies within L*d of f(clamp(x)).
            // `dist_out_hi`, the product, and both interval additions are all
            // outward, so endpoint arithmetic cannot narrow the claim.
            let v = base + dist_out;
            if !v.is_finite() {
                return ChartSample {
                    signed_distance: f64::NAN,
                    gradient: None,
                    lipschitz: None,
                    error: NumericalCertificate::no_claim(),
                };
            }
            let excursion = outward_mul_nonnegative(self.source_lipschitz, dist_out_hi);
            let source_at_x = if excursion.is_finite() {
                at_clamp + Interval::new(-excursion, excursion)
            } else {
                Interval::WHOLE
            };
            ChartSample {
                signed_distance: v,
                gradient: None,
                lipschitz: None,
                // Include the published nominal as well as the source theorem;
                // malformed certificates that exclude their nominal fail closed
                // in downstream authority checks.
                error: abstract_error(source_at_x.hull(Interval::point(v))),
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

/// Blanket field sampling into a [`SampledSdf`]. Only the global
/// [`TraceStepClaim::ExactDistance`] theorem can make the interpolation bound
/// rigorous relative to abstract Euclidean distance. Other charts can retain
/// nominal/estimate evidence through [`Convert::convert_with_domain`] or
/// [`Convert::convert_clipped`], but cannot be promoted by [`Convert::convert`].
impl<C: Chart> Convert<SampledSdf> for C {
    #[allow(clippy::too_many_lines)] // One fail-closed admission, sampling, bound, and receipt transaction.
    fn convert_with_domain(
        &self,
        budget: ErrBudget,
        explicit_clip: Option<Aabb>,
        cx: &Cx<'_>,
    ) -> Result<Evidence<SampledSdf>, ConvertDiag> {
        cx.checkpoint().map_err(|_| ConvertDiag::Cancelled {
            stage: "admission",
            completed_samples: 0,
        })?;
        if !budget.abs_sd_error.is_finite() || budget.abs_sd_error <= 0.0 {
            return Err(ConvertDiag::InvalidBudget {
                requested_bits: budget.abs_sd_error.to_bits(),
            });
        }
        let raw_support = self.support();
        SamplingDomain::validate_support(raw_support).map_err(ConvertDiag::SamplingDomain)?;
        let clipped = explicit_clip
            .map(|clip| ClippedChart::new(self, clip))
            .transpose()
            .map_err(ConvertDiag::SamplingDomain)?;
        let source: &dyn Chart = clipped
            .as_ref()
            .map_or(self as &dyn Chart, |chart| chart as &dyn Chart);
        let exact_distance_source = source.trace_step_claim() == TraceStepClaim::ExactDistance;
        let padded_support = source.support().inflate(budget.abs_sd_error.max(1e-9));
        let domain =
            SamplingDomain::resolve(padded_support, None).map_err(ConvertDiag::SamplingDomain)?;
        let box_ = domain.bounds();
        // Probe the source's Lipschitz claim at the box center.
        let center = domain.midpoint();
        let center_sample = source.eval(center, cx);
        cx.checkpoint().map_err(|_| ConvertDiag::Cancelled {
            stage: "source-probe",
            completed_samples: 0,
        })?;
        if !center_sample.signed_distance.is_finite() {
            return Err(ConvertDiag::NonFiniteSignedDistance {
                point: center,
                value_bits: center_sample.signed_distance.to_bits(),
            });
        }
        let sampled_center_lipschitz = match center_sample.lipschitz {
            Some(l) if l.is_finite() && l >= 0.0 => l,
            _ if exact_distance_source => 1.0,
            _ => return Err(ConvertDiag::NoLipschitzBound),
        };
        let lipschitz = if exact_distance_source {
            1.0
        } else {
            sampled_center_lipschitz
        };
        let (center_authority, center_error_radius) =
            sample_abstract_distance_authority(&center_sample);
        // Any trilinear value is a convex combination of the eight corners.
        // Therefore an L-Lipschitz source differs from that combination by at
        // most L times the FULL cell diagonal. Use an outward domain diagonal
        // to select a first resolution; the final bound below is recomputed
        // from the actual representable nodes rather than ideal uniform steps.
        let domain_diagonal = outward_box_diagonal(box_).unwrap_or(f64::INFINITY);
        let full_box_radius = outward_mul_nonnegative(lipschitz, domain_diagonal);
        let need_resolution = resolution_for_radius(full_box_radius, budget.abs_sd_error);
        if need_resolution > SAMPLED_SDF_MAX_RESOLUTION {
            let achievable =
                outward_div_nonnegative(full_box_radius, f64::from(SAMPLED_SDF_MAX_RESOLUTION - 1));
            return Err(ConvertDiag::BudgetInfeasible {
                requested: budget.abs_sd_error,
                achievable,
                need_resolution,
                cap: SAMPLED_SDF_MAX_RESOLUTION,
            });
        }
        // Keep one node of slack when the dense cap allows it. The requested
        // budget is relative to abstract distance, so outward-rounded source
        // certificates (often a few ulps wide) must fit alongside the nominal
        // interpolation error rather than turning every exactly-on-grid budget
        // into a false infeasibility.
        let required_n = need_resolution.max(2);
        let preferred_n = required_n
            .max(2)
            .saturating_add(1)
            .min(SAMPLED_SDF_MAX_RESOLUTION);
        let build_nodes = |n: u32| -> Result<[Vec<f64>; 3], ConvertDiag> {
            Ok([
                build_axis_nodes(Axis::X, box_.min.x, box_.max.x, n)?,
                build_axis_nodes(Axis::Y, box_.min.y, box_.max.y, n)?,
                build_axis_nodes(Axis::Z, box_.min.z, box_.max.z, n)?,
            ])
        };
        // The extra node normally leaves room for outward source and arithmetic
        // radii. It is only slack: if that denser grid cannot be represented,
        // fall back to the mathematically required resolution before refusing.
        let axis_nodes = match build_nodes(preferred_n) {
            Ok(nodes) => nodes,
            Err(_) if preferred_n > required_n => build_nodes(required_n)?,
            Err(error) => return Err(error),
        };
        let n = axis_nodes[0].len() as u32;
        let cell_diagonal = max_cell_diagonal(&axis_nodes).unwrap_or(f64::INFINITY);
        // Exact-distance sources supply the global unit-Lipschitz theorem.
        // For weaker fields, sample-local maxima only parameterize a nominal
        // reconstruction estimate; they never mint Enclosure authority because
        // sub-grid slope spikes remain unresolved.
        let mut l_max = lipschitz;
        let mut abstract_distance_kind = center_authority;
        let mut source_error_radius = center_error_radius.unwrap_or(0.0);
        let mut max_sample_abs = center_sample.signed_distance.abs();
        let mut values = Vec::with_capacity((n as usize).pow(3));
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    cx.checkpoint().map_err(|_| ConvertDiag::Cancelled {
                        stage: "sampling-grid",
                        completed_samples: values.len(),
                    })?;
                    let p = Point3::new(
                        axis_nodes[0][i as usize],
                        axis_nodes[1][j as usize],
                        axis_nodes[2][k as usize],
                    );
                    let sample = source.eval(p, cx);
                    cx.checkpoint().map_err(|_| ConvertDiag::Cancelled {
                        stage: "sampling-grid",
                        completed_samples: values.len(),
                    })?;
                    if !sample.signed_distance.is_finite() {
                        return Err(ConvertDiag::NonFiniteSignedDistance {
                            point: p,
                            value_bits: sample.signed_distance.to_bits(),
                        });
                    }
                    if !exact_distance_source {
                        match sample.lipschitz {
                            Some(l) if l.is_finite() && l >= 0.0 => l_max = l_max.max(l),
                            // A weak chart still needs finite local slope data
                            // for its nominal reconstruction estimate.
                            _ => return Err(ConvertDiag::NoLipschitzBound),
                        }
                    }
                    let (sample_kind, sample_radius) = sample_abstract_distance_authority(&sample);
                    abstract_distance_kind = abstract_distance_kind.max(sample_kind);
                    if let Some(radius) = sample_radius {
                        source_error_radius = source_error_radius.max(radius);
                    }
                    max_sample_abs = max_sample_abs.max(sample.signed_distance.abs());
                    values.push(sample.signed_distance);
                }
            }
        }
        // `cell_diagonal` covers arbitrary convex corner weights against the
        // source field. The second term covers all floating arithmetic in the
        // seven stable lerps, including computed complements. Both terms are
        // finite and outward before being published.
        let reconstruction_radius = outward_mul_nonnegative(l_max, cell_diagonal);
        let interpolation_roundoff = interpolation_roundoff_bound(max_sample_abs);
        let nominal_field_bound =
            outward_add_nonnegative(reconstruction_radius, interpolation_roundoff);
        // If the grid revealed a steeper slope than the center probe assumed,
        // the honest bound may exceed the budget — REFUSE with the true
        // achievable rather than ship a receipt that understates the error.
        if !nominal_field_bound.is_finite() || nominal_field_bound > budget.abs_sd_error {
            let remaining = budget.abs_sd_error - interpolation_roundoff;
            let need_resolution = if remaining.is_finite() && remaining > 0.0 {
                resolution_for_radius(outward_mul_nonnegative(l_max, domain_diagonal), remaining)
            } else {
                u32::MAX
            };
            return Err(ConvertDiag::BudgetInfeasible {
                requested: budget.abs_sd_error,
                achievable: nominal_field_bound,
                need_resolution,
                cap: SAMPLED_SDF_MAX_RESOLUTION,
            });
        }
        // Trilinear reconstruction is approximate even when every source
        // sample was exact. NoClaim absorbs; otherwise the source certificate
        // radius composes with the nominal field reconstruction bound.
        abstract_distance_kind = abstract_distance_kind.max(NumericalKind::Enclosure);
        if !exact_distance_source && abstract_distance_kind != NumericalKind::NoClaim {
            abstract_distance_kind = abstract_distance_kind.max(NumericalKind::Estimate);
        }
        let abstract_distance_bound = if abstract_distance_kind == NumericalKind::NoClaim {
            None
        } else {
            let bound = outward_add_nonnegative(nominal_field_bound, source_error_radius);
            if bound.is_finite() {
                Some(bound)
            } else {
                abstract_distance_kind = NumericalKind::NoClaim;
                None
            }
        };
        if let Some(total_bound) = abstract_distance_bound
            && total_bound > budget.abs_sd_error
        {
            let committed = outward_add_nonnegative(source_error_radius, interpolation_roundoff);
            let remaining = budget.abs_sd_error - committed;
            let need_resolution = if remaining.is_finite() && remaining > 0.0 {
                resolution_for_radius(outward_mul_nonnegative(l_max, domain_diagonal), remaining)
            } else {
                u32::MAX
            };
            return Err(ConvertDiag::BudgetInfeasible {
                requested: budget.abs_sd_error,
                achievable: total_bound,
                need_resolution,
                cap: SAMPLED_SDF_MAX_RESOLUTION,
            });
        }
        let chart = SampledSdf {
            box_,
            n,
            axis_nodes,
            values,
            nominal_field_bound,
            // The reconstruction radius inherits the source's Lipschitz
            // authority: only `ExactDistance` supplies a theorem that holds
            // ACROSS a cell. A maximum over per-node local bounds does not,
            // so it is published as an Estimate.
            nominal_field_bound_kind: if exact_distance_source {
                NumericalKind::Enclosure
            } else {
                NumericalKind::Estimate
            },
            abstract_distance_kind,
            abstract_distance_bound,
            source_lipschitz: l_max,
        };
        // The receipt's QoI is the total abstract-distance bound when one
        // exists, otherwise the finite nominal-field reconstruction bound. Its
        // numerical class prevents the latter from being mistaken for
        // abstract-distance authority.
        let mut provenance_parents = vec![ProvenanceHash::of_bytes(self.name().as_bytes())];
        if let Some(clip) = explicit_clip {
            provenance_parents.push(aabb_provenance(clip));
        }
        let qoi = abstract_distance_bound.unwrap_or(nominal_field_bound);
        let numerical = match abstract_distance_kind {
            NumericalKind::Exact | NumericalKind::Enclosure => {
                NumericalCertificate::enclosure(0.0, qoi)
            }
            NumericalKind::Estimate => NumericalCertificate::estimate(0.0, qoi),
            NumericalKind::NoClaim => NumericalCertificate::no_claim(),
        };
        cx.checkpoint().map_err(|_| ConvertDiag::Cancelled {
            stage: "publication",
            completed_samples: chart.values.len(),
        })?;
        Ok(Evidence {
            qoi,
            numerical,
            statistical: fs_evidence::StatisticalCertificate::None,
            model: fs_evidence::ModelEvidence::none(),
            sensitivity: fs_evidence::SensitivitySummary::default(),
            provenance: ProvenanceHash::chain("convert/sampled-sdf", &provenance_parents),
            adjoint_ref: None,
            value: chart,
        })
    }
}

fn build_axis_nodes(
    axis: Axis,
    min: f64,
    max: f64,
    resolution: u32,
) -> Result<Vec<f64>, ConvertDiag> {
    debug_assert!(resolution >= 2);
    let span = max - min;
    let denominator = f64::from(resolution - 1);
    let mut nodes = Vec::with_capacity(resolution as usize);
    nodes.push(min);
    for index in 1..resolution - 1 {
        // The coordinate need not equal an ideal real-uniform node: the actual
        // representable cells are retained and bounded below. It MUST, however,
        // be a new interior point. Coincident rounded nodes make interpolation
        // and reconstruction geometry undefined and are refused before sampling.
        let t = f64::from(index) / denominator;
        let node = min + span * t;
        if !node.is_finite() || node <= nodes[nodes.len() - 1] || node >= max {
            return Err(ConvertDiag::UnrepresentableGrid {
                axis,
                index,
                resolution,
                min_bits: min.to_bits(),
                max_bits: max.to_bits(),
            });
        }
        nodes.push(node);
    }
    if max <= nodes[nodes.len() - 1] {
        return Err(ConvertDiag::UnrepresentableGrid {
            axis,
            index: resolution - 1,
            resolution,
            min_bits: min.to_bits(),
            max_bits: max.to_bits(),
        });
    }
    nodes.push(max);
    Ok(nodes)
}

fn locate_axis(nodes: &[f64], coordinate: f64) -> (usize, usize, f64) {
    debug_assert!(nodes.len() >= 2);
    let last = nodes.len() - 1;
    if coordinate <= nodes[0] {
        return (0, 1, 0.0);
    }
    if coordinate >= nodes[last] {
        return (last - 1, last, 1.0);
    }
    let upper = nodes.partition_point(|node| *node <= coordinate);
    let lower = upper - 1;
    let width = nodes[upper] - nodes[lower];
    let t = ((coordinate - nodes[lower]) / width).clamp(0.0, 1.0);
    (lower, upper, t)
}

/// Overflow-resistant interpolation. Same-sign endpoints use the exact convex
/// formula `a + (b-a)t`, whose difference cannot overflow. Opposite-sign
/// endpoints use scaled products, whose terms cannot exceed their endpoints.
/// The chart's stored roundoff allowance covers all seven calls.
fn stable_lerp(a: f64, b: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return a;
    }
    if t >= 1.0 {
        return b;
    }
    if a.is_sign_negative() == b.is_sign_negative() {
        a + (b - a) * t
    } else {
        a * (1.0 - t) + b * t
    }
}

fn outward_width(lo: f64, hi: f64) -> Option<f64> {
    let width = hi - lo;
    if !width.is_finite() || width < 0.0 {
        return None;
    }
    if width == 0.0 {
        Some(0.0)
    } else {
        let upper = width.next_up();
        upper.is_finite().then_some(upper)
    }
}

fn outward_add_nonnegative(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        return b;
    }
    if b == 0.0 {
        return a;
    }
    let sum = a + b;
    if sum.is_finite() {
        sum.next_up()
    } else {
        f64::INFINITY
    }
}

fn outward_mul_nonnegative(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        return 0.0;
    }
    let product = a * b;
    if product.is_finite() {
        product.next_up()
    } else {
        f64::INFINITY
    }
}

fn outward_div_nonnegative(numerator: f64, denominator: f64) -> f64 {
    if numerator == 0.0 {
        return 0.0;
    }
    let quotient = numerator / denominator;
    if quotient.is_finite() {
        quotient.next_up()
    } else {
        f64::INFINITY
    }
}

fn outward_norm3(x: f64, y: f64, z: f64) -> Option<f64> {
    let scale = x.max(y).max(z);
    if !scale.is_finite() || scale < 0.0 {
        return None;
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let square = |value: f64| {
        let ratio = outward_div_nonnegative(value, scale);
        outward_mul_nonnegative(ratio, ratio)
    };
    let sum = outward_add_nonnegative(outward_add_nonnegative(square(x), square(y)), square(z));
    let root = sum.sqrt().next_up();
    let norm = outward_mul_nonnegative(scale, root);
    norm.is_finite().then_some(norm)
}

fn outward_box_diagonal(box_: Aabb) -> Option<f64> {
    outward_norm3(
        outward_width(box_.min.x, box_.max.x)?,
        outward_width(box_.min.y, box_.max.y)?,
        outward_width(box_.min.z, box_.max.z)?,
    )
}

fn max_cell_diagonal(axis_nodes: &[Vec<f64>; 3]) -> Option<f64> {
    let mut widths = [0.0_f64; 3];
    for (axis, nodes) in axis_nodes.iter().enumerate() {
        for pair in nodes.windows(2) {
            widths[axis] = widths[axis].max(outward_width(pair[0], pair[1])?);
        }
    }
    outward_norm3(widths[0], widths[1], widths[2])
}

fn resolution_for_radius(full_box_radius: f64, budget: f64) -> u32 {
    let cells = outward_div_nonnegative(full_box_radius, budget).ceil();
    if !cells.is_finite() || cells >= f64::from(u32::MAX - 1) {
        return u32::MAX;
    }
    (cells as u32).saturating_add(1).max(2)
}

/// Seven lerps, at no more than four elementary operations each (including a
/// computed complement), have a forward error comfortably below 64 eps times
/// the largest corner magnitude. The subnormal floor covers absolute rounding
/// when the relative model underflows. Both contributions are nudged outward.
fn interpolation_roundoff_bound(max_abs_value: f64) -> f64 {
    if !max_abs_value.is_finite() {
        return f64::INFINITY;
    }
    let relative = (64.0 * f64::EPSILON).next_up();
    let scaled = outward_mul_nonnegative(max_abs_value, relative);
    let subnormal_floor = 64.0 * f64::from_bits(1);
    scaled.max(subnormal_floor).next_up()
}

fn centered_interval(center: f64, radius: f64) -> Interval {
    if !center.is_finite() || !radius.is_finite() || radius < 0.0 {
        return Interval::WHOLE;
    }
    Interval::point(center) + Interval::new(-radius, radius)
}

fn finite_norm_with_upper(x: f64, y: f64, z: f64) -> Option<(f64, f64)> {
    let (x, y, z) = (x.abs(), y.abs(), z.abs());
    let scale = x.max(y).max(z);
    if !scale.is_finite() {
        return None;
    }
    if scale == 0.0 {
        return Some((0.0, 0.0));
    }
    let nominal = scale
        * ((x / scale) * (x / scale) + (y / scale) * (y / scale) + (z / scale) * (z / scale))
            .sqrt();
    // Each component came from a floating subtraction (`query - clamp`). One
    // outward ULP per nonzero component encloses the corresponding exact real
    // coordinate difference before the norm arithmetic is bounded.
    let component_upper = |value: f64| if value == 0.0 { 0.0 } else { value.next_up() };
    let upper = outward_norm3(component_upper(x), component_upper(y), component_upper(z))?;
    nominal.is_finite().then_some((nominal, upper))
}

/// Validate how strongly one source sample relates its nominal field value to
/// abstract region signed distance. Public certificate fields are deliberately
/// forgeable for composition ergonomics, so malformed finite claims fail
/// closed to `NoClaim` instead of being laundered by sampling.
fn sample_abstract_distance_authority(sample: &ChartSample) -> (NumericalKind, Option<f64>) {
    let certificate = sample.error;
    if certificate.kind == NumericalKind::NoClaim {
        return (NumericalKind::NoClaim, None);
    }
    let finite_ordered = certificate.lo.is_finite()
        && certificate.hi.is_finite()
        && certificate.lo <= certificate.hi;
    let consistent = match certificate.kind {
        NumericalKind::Exact => {
            certificate.lo.to_bits() == sample.signed_distance.to_bits()
                && certificate.hi.to_bits() == sample.signed_distance.to_bits()
        }
        NumericalKind::Enclosure | NumericalKind::Estimate => {
            finite_ordered
                && certificate.lo <= sample.signed_distance
                && sample.signed_distance <= certificate.hi
        }
        NumericalKind::NoClaim => false,
    };
    if !finite_ordered || !consistent {
        return (NumericalKind::NoClaim, None);
    }
    let radius = if certificate.kind == NumericalKind::Exact {
        0.0
    } else {
        outward_width(certificate.lo, sample.signed_distance)
            .and_then(|lo_radius| {
                outward_width(sample.signed_distance, certificate.hi)
                    .map(|hi_radius| lo_radius.max(hi_radius))
            })
            .unwrap_or(f64::INFINITY)
    };
    if radius.is_finite() {
        (certificate.kind, Some(radius))
    } else {
        (NumericalKind::NoClaim, None)
    }
}

fn aabb_provenance(box_: Aabb) -> ProvenanceHash {
    let mut bytes = [0u8; 48];
    for (index, value) in [
        box_.min.x, box_.min.y, box_.min.z, box_.max.x, box_.max.y, box_.max.z,
    ]
    .into_iter()
    .enumerate()
    {
        let start = index * 8;
        bytes[start..start + 8].copy_from_slice(&value.to_bits().to_le_bytes());
    }
    ProvenanceHash::of_bytes(&bytes)
}
