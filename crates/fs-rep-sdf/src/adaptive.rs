//! Adaptively-sampled SDF (plan §7.2): an octree whose leaf cells carry
//! trilinear corner fits, refined where the fit residual against the
//! source exceeds tolerance — the compact representation for smooth
//! shapes with localized features. Residual bounds are MEASURED at probe
//! points and ledgered (an Estimate-grade error model, honestly labeled),
//! then composed with the weakest source-sample authority. A source NoClaim
//! remains NoClaim; interval-verified fits are fs-ivl integration work.

use crate::dense::{SdfBuildError, finite_positive, sample_abstract_distance_authority};
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, ClippedChart, Differentiability, Point3, SamplingDomain};
use std::fmt::Write as _;

/// Deterministic upper bound on octree nodes admitted by an adaptive build.
/// The bound is checked from `max_depth` before source evaluation or allocation.
pub const ADAPTIVE_MAX_NODES: u64 = 1_000_000;

/// One octree node: either a leaf with 8 corner samples or 8 children.
enum Node {
    Leaf { corners: [f64; 8] },
    Branch { children: Box<[Node; 8]> },
}

/// The adaptive chart.
pub struct AdaptiveSdf {
    root: Node,
    box_: Aabb,
    /// Max observed fit residual at probe points (Estimate-grade).
    residual: f64,
    /// Refinement tolerance the build targeted.
    tol: f64,
    /// Cells (leaves) in the tree.
    cells: u64,
    /// Deepest refinement level reached.
    depth: u32,
    /// The source's certified Lipschitz constant (for outside-box math).
    source_lipschitz: f64,
    /// Weakest source authority after the probed fit demotes all finite
    /// authority to Estimate.
    abstract_distance_kind: NumericalKind,
    /// Probed fit band plus maximum source certificate radius. `None` means
    /// honest NoClaim.
    abstract_distance_bound: Option<f64>,
}

struct SourceAuthority {
    kind: NumericalKind,
    max_radius: f64,
}

impl SourceAuthority {
    fn new() -> Self {
        Self {
            kind: NumericalKind::Exact,
            max_radius: 0.0,
        }
    }

    fn observe(&mut self, sample: &ChartSample, point: Point3) -> Result<(), SdfBuildError> {
        if !sample.signed_distance.is_finite() {
            return Err(SdfBuildError::InvalidSample {
                point,
                value_bits: sample.signed_distance.to_bits(),
            });
        }
        let (kind, radius) = sample_abstract_distance_authority(sample);
        self.kind = self.kind.max(kind);
        if let Some(radius) = radius {
            self.max_radius = self.max_radius.max(radius);
        }
        Ok(())
    }
}

/// Build statistics (ledgered: the "fit residual bounds ledgered"
/// requirement).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveStats {
    /// Leaf-cell count.
    pub cells: u64,
    /// Deepest level.
    pub depth: u32,
    /// Max observed probe residual.
    pub residual: f64,
}

impl AdaptiveStats {
    /// Canonical JSON.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(64);
        let _ = write!(
            s,
            "{{\"cells\":{},\"depth\":{},\"residual\":{:.6}}}",
            self.cells, self.depth, self.residual
        );
        s
    }
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
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

fn trilinear(corners: &[f64; 8], t: [f64; 3]) -> f64 {
    let c00 = lerp(corners[0], corners[1], t[0]);
    let c10 = lerp(corners[2], corners[3], t[0]);
    let c01 = lerp(corners[4], corners[5], t[0]);
    let c11 = lerp(corners[6], corners[7], t[0]);
    lerp(lerp(c00, c10, t[1]), lerp(c01, c11, t[1]), t[2])
}

fn corner_point(b: &Aabb, idx: usize) -> Point3 {
    Point3::new(
        if idx & 1 == 0 { b.min.x } else { b.max.x },
        if (idx >> 1) & 1 == 0 {
            b.min.y
        } else {
            b.max.y
        },
        if (idx >> 2) & 1 == 0 {
            b.min.z
        } else {
            b.max.z
        },
    )
}

fn cell_midpoint(b: &Aabb) -> Point3 {
    Point3::new(
        f64::midpoint(b.min.x, b.max.x),
        f64::midpoint(b.min.y, b.max.y),
        f64::midpoint(b.min.z, b.max.z),
    )
}

fn checked_subdivision_midpoint(b: &Aabb) -> Result<Point3, SdfBuildError> {
    let mid = cell_midpoint(b);
    for (axis, (min, max, value)) in [
        (b.min.x, b.max.x, mid.x),
        (b.min.y, b.max.y, mid.y),
        (b.min.z, b.max.z, mid.z),
    ]
    .into_iter()
    .enumerate()
    {
        let strictly_interior = value.is_finite() && min < value && value < max;
        if !strictly_interior {
            return Err(SdfBuildError::AdaptiveSubdivisionUnrepresentable {
                axis,
                min_bits: min.to_bits(),
                max_bits: max.to_bits(),
                midpoint_bits: value.to_bits(),
            });
        }
    }
    Ok(mid)
}

fn octant_at(b: &Aabb, mid: Point3, idx: usize) -> Aabb {
    let min = Point3::new(
        if idx & 1 == 0 { b.min.x } else { mid.x },
        if (idx >> 1) & 1 == 0 { b.min.y } else { mid.y },
        if (idx >> 2) & 1 == 0 { b.min.z } else { mid.z },
    );
    let max = Point3::new(
        if idx & 1 == 0 { mid.x } else { b.max.x },
        if (idx >> 1) & 1 == 0 { mid.y } else { b.max.y },
        if (idx >> 2) & 1 == 0 { mid.z } else { b.max.z },
    );
    Aabb::new(min, max)
}

fn octant(b: &Aabb, idx: usize) -> Aabb {
    octant_at(b, cell_midpoint(b), idx)
}

fn worst_case_octree_nodes(max_depth: u32) -> Option<u128> {
    let mut level = 1u128;
    let mut total = 1u128;
    for _ in 0..max_depth {
        level = level.checked_mul(8)?;
        total = total.checked_add(level)?;
    }
    Some(total)
}

impl AdaptiveSdf {
    /// Build over `source`'s inflated support, splitting cells whose fit
    /// residual (probed at the center and the six face centers) exceeds
    /// `tol`, down to `max_depth`. Polls cancellation per cell.
    ///
    /// # Errors
    /// [`SdfBuildError`] when the tolerance, finite sampling domain, or
    /// deterministic work bound is inadmissible, or when cancellation is
    /// observed mid-build.
    pub fn build(
        source: &dyn Chart,
        tol: f64,
        max_depth: u32,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveSdf, SdfBuildError> {
        let tol = finite_positive(tol, "tol")?;
        let worst_case = worst_case_octree_nodes(max_depth).unwrap_or(u128::MAX);
        if worst_case > u128::from(ADAPTIVE_MAX_NODES) {
            return Err(SdfBuildError::AdaptiveWorkLimit {
                need: worst_case,
                cap: ADAPTIVE_MAX_NODES,
            });
        }
        let support = SamplingDomain::admit(source.support(), None)?.bounds();
        let box_ = SamplingDomain::admit(support.inflate(tol.max(1e-9)), None)?.bounds();
        Self::build_in_domain(source, box_, tol, max_depth, cx)
    }

    /// Build an adaptive SDF of the geometric intersection `source ∩ clip`.
    /// The explicit clip and worst-case octree work are admitted before any
    /// source evaluation or node allocation.
    ///
    /// # Errors
    /// [`SdfBuildError`] under the same conditions as [`Self::build`], plus
    /// an invalid, empty, or degenerate explicit clip.
    pub fn build_clipped(
        source: &dyn Chart,
        tol: f64,
        max_depth: u32,
        clip: Aabb,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveSdf, SdfBuildError> {
        let clipped = ClippedChart::new(source, clip)?;
        Self::build(&clipped, tol, max_depth, cx)
    }

    fn build_in_domain(
        source: &dyn Chart,
        box_: Aabb,
        tol: f64,
        max_depth: u32,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveSdf, SdfBuildError> {
        let center = cell_midpoint(&box_);
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        let center_sample = source.eval(center, cx);
        let mut authority = SourceAuthority::new();
        authority.observe(&center_sample, center)?;
        let lipschitz = center_sample
            .lipschitz
            .filter(|bound| bound.is_finite() && *bound >= 0.0)
            .unwrap_or(1.0);
        let mut cells = 0u64;
        let mut depth_seen = 0u32;
        let mut residual = 0.0f64;
        let root = Self::build_node(
            source,
            &box_,
            tol,
            max_depth,
            0,
            cx,
            &mut cells,
            &mut depth_seen,
            &mut residual,
            &mut authority,
        )?;
        let nominal_field_bound = residual.max(tol);
        let mut abstract_distance_kind = authority.kind.max(NumericalKind::Estimate);
        let abstract_distance_bound = if abstract_distance_kind == NumericalKind::NoClaim {
            None
        } else {
            let bound = nominal_field_bound + authority.max_radius;
            if bound.is_finite() {
                Some(bound)
            } else {
                abstract_distance_kind = NumericalKind::NoClaim;
                None
            }
        };
        Ok(AdaptiveSdf {
            root,
            box_,
            residual,
            tol,
            cells,
            depth: depth_seen,
            source_lipschitz: lipschitz,
            abstract_distance_kind,
            abstract_distance_bound,
        })
    }

    #[allow(clippy::too_many_arguments)] // recursive builder plumbing
    fn build_node(
        source: &dyn Chart,
        b: &Aabb,
        tol: f64,
        max_depth: u32,
        depth: u32,
        cx: &Cx<'_>,
        cells: &mut u64,
        depth_seen: &mut u32,
        residual: &mut f64,
        authority: &mut SourceAuthority,
    ) -> Result<Node, SdfBuildError> {
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        let mut corners = [0.0f64; 8];
        for (i, corner) in corners.iter_mut().enumerate() {
            let point = corner_point(b, i);
            let sample = source.eval(point, cx);
            authority.observe(&sample, point)?;
            *corner = sample.signed_distance;
            cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        }
        // Probe the fit at the center + face centers.
        let mid = cell_midpoint(b);
        let probes = [
            mid,
            Point3::new(b.min.x, mid.y, mid.z),
            Point3::new(b.max.x, mid.y, mid.z),
            Point3::new(mid.x, b.min.y, mid.z),
            Point3::new(mid.x, b.max.y, mid.z),
            Point3::new(mid.x, mid.y, b.min.z),
            Point3::new(mid.x, mid.y, b.max.z),
        ];
        let mut worst = 0.0f64;
        for p in probes {
            let t = [
                (p.x - b.min.x) / (b.max.x - b.min.x),
                (p.y - b.min.y) / (b.max.y - b.min.y),
                (p.z - b.min.z) / (b.max.z - b.min.z),
            ];
            let fit = trilinear(&corners, t);
            if !fit.is_finite() {
                return Err(SdfBuildError::InvalidReconstructionBound { value: fit });
            }
            let sample = source.eval(p, cx);
            authority.observe(&sample, p)?;
            let difference = fit - sample.signed_distance;
            if !difference.is_finite() {
                return Err(SdfBuildError::InvalidReconstructionBound { value: difference });
            }
            let probe_residual = difference.abs();
            if !probe_residual.is_finite() {
                return Err(SdfBuildError::InvalidReconstructionBound {
                    value: probe_residual,
                });
            }
            worst = worst.max(probe_residual);
            cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        }
        if worst <= tol || depth >= max_depth {
            *cells = cells
                .checked_add(1)
                .ok_or(SdfBuildError::AdaptiveWorkLimit {
                    need: u128::MAX,
                    cap: ADAPTIVE_MAX_NODES,
                })?;
            *depth_seen = (*depth_seen).max(depth);
            *residual = residual.max(worst);
            return Ok(Node::Leaf { corners });
        }
        let subdivision_midpoint = checked_subdivision_midpoint(b)?;
        let mut children: Vec<Node> = Vec::with_capacity(8);
        for i in 0..8 {
            children.push(Self::build_node(
                source,
                &octant_at(b, subdivision_midpoint, i),
                tol,
                max_depth,
                depth + 1,
                cx,
                cells,
                depth_seen,
                residual,
                authority,
            )?);
        }
        Ok(Node::Branch {
            children: Box::new(
                children
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("exactly 8 children pushed")),
            ),
        })
    }

    /// Build statistics (ledgered evidence).
    #[must_use]
    pub fn stats(&self) -> AdaptiveStats {
        AdaptiveStats {
            cells: self.cells,
            depth: self.depth,
            residual: self.residual,
        }
    }

    /// Probed fit band relative to the sampled source field.
    #[must_use]
    pub fn nominal_field_bound(&self) -> f64 {
        self.residual.max(self.tol)
    }

    /// Weakest authority relative to abstract region signed distance.
    /// Adaptive probing makes every finite result at most Estimate-grade.
    #[must_use]
    pub fn abstract_distance_kind(&self) -> NumericalKind {
        self.abstract_distance_kind
    }

    /// Probed fit band plus maximum finite source-certificate radius, or
    /// `None` when any observed source sample made no valid claim.
    #[must_use]
    pub fn abstract_distance_bound(&self) -> Option<f64> {
        self.abstract_distance_bound
    }

    fn eval_in_box(&self, p: Point3) -> f64 {
        let mut node = &self.root;
        let mut b = self.box_;
        loop {
            match node {
                Node::Leaf { corners } => {
                    let t = [
                        ((p.x - b.min.x) / (b.max.x - b.min.x)).clamp(0.0, 1.0),
                        ((p.y - b.min.y) / (b.max.y - b.min.y)).clamp(0.0, 1.0),
                        ((p.z - b.min.z) / (b.max.z - b.min.z)).clamp(0.0, 1.0),
                    ];
                    return trilinear(corners, t);
                }
                Node::Branch { children } => {
                    let mid = cell_midpoint(&b);
                    let idx = usize::from(p.x >= mid.x)
                        | (usize::from(p.y >= mid.y) << 1)
                        | (usize::from(p.z >= mid.z) << 2);
                    b = octant(&b, idx);
                    node = &children[idx];
                }
            }
        }
    }
}

impl Chart for AdaptiveSdf {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let clamped = Point3::new(
            x.x.clamp(self.box_.min.x, self.box_.max.x),
            x.y.clamp(self.box_.min.y, self.box_.max.y),
            x.z.clamp(self.box_.min.z, self.box_.max.z),
        );
        let dist_out = x.delta_from(clamped).norm();
        let v = self.eval_in_box(clamped) + dist_out;
        // At best Estimate-grade: the residual is probed, not enclosed
        // (fs-ivl integration promotes this later). A source NoClaim remains
        // absorbing even outside the sampled box.
        let band = self
            .abstract_distance_bound
            .unwrap_or(self.nominal_field_bound());
        let error = match self.abstract_distance_kind {
            NumericalKind::NoClaim => NumericalCertificate::no_claim(),
            NumericalKind::Exact | NumericalKind::Enclosure | NumericalKind::Estimate => {
                let lo = v - band - (1.0 + self.source_lipschitz) * dist_out;
                let hi = v + band;
                if lo.is_finite() && hi.is_finite() && lo <= hi {
                    NumericalCertificate::estimate(lo, hi)
                } else {
                    NumericalCertificate::no_claim()
                }
            }
        };
        ChartSample {
            signed_distance: v,
            gradient: None,
            lipschitz: None,
            error,
        }
    }

    fn support(&self) -> Aabb {
        self.box_
    }

    fn name(&self) -> &'static str {
        "rep-sdf/adaptive"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

#[cfg(test)]
mod tests {
    use super::lerp;

    #[test]
    fn opposite_sign_extreme_lerp_remains_a_finite_convex_value() {
        let low = -f64::MAX;
        let high = f64::MAX;
        assert_eq!(lerp(low, high, 0.0).to_bits(), low.to_bits());
        assert_eq!(lerp(low, high, 1.0).to_bits(), high.to_bits());
        let midpoint = lerp(low, high, 0.5);
        assert!(midpoint.is_finite());
        assert_eq!(midpoint.to_bits(), 0.0f64.to_bits());
    }
}
