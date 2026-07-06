//! Adaptively-sampled SDF (plan §7.2): an octree whose leaf cells carry
//! trilinear corner fits, refined where the fit residual against the
//! source exceeds tolerance — the compact representation for smooth
//! shapes with localized features. Residual bounds are MEASURED at probe
//! points and ledgered (an Estimate-grade error model, honestly labeled —
//! interval-verified fits are fs-ivl integration work, CONTRACT
//! no-claims).

use fs_evidence::NumericalCertificate;
use fs_exec::{Cancelled, Cx};
use fs_geom::{Aabb, Chart, ChartSample, Differentiability, Point3};
use std::fmt::Write as _;

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
    a + (b - a) * t
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
        if (idx >> 1) & 1 == 0 { b.min.y } else { b.max.y },
        if (idx >> 2) & 1 == 0 { b.min.z } else { b.max.z },
    )
}

fn octant(b: &Aabb, idx: usize) -> Aabb {
    let mid = Point3::new(
        f64::midpoint(b.min.x, b.max.x),
        f64::midpoint(b.min.y, b.max.y),
        f64::midpoint(b.min.z, b.max.z),
    );
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

impl AdaptiveSdf {
    /// Build over `source`'s inflated support, splitting cells whose fit
    /// residual (probed at the center and the six face centers) exceeds
    /// `tol`, down to `max_depth`. Polls cancellation per cell.
    ///
    /// # Errors
    /// [`Cancelled`] mid-build.
    pub fn build(
        source: &dyn Chart,
        tol: f64,
        max_depth: u32,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveSdf, Cancelled> {
        let box_ = source.support().inflate(tol.max(1e-9));
        let lipschitz = source
            .eval(
                Point3::new(
                    f64::midpoint(box_.min.x, box_.max.x),
                    f64::midpoint(box_.min.y, box_.max.y),
                    f64::midpoint(box_.min.z, box_.max.z),
                ),
                cx,
            )
            .lipschitz
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
        )?;
        Ok(AdaptiveSdf {
            root,
            box_,
            residual,
            tol,
            cells,
            depth: depth_seen,
            source_lipschitz: lipschitz,
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
    ) -> Result<Node, Cancelled> {
        cx.checkpoint()?;
        let corners: [f64; 8] =
            core::array::from_fn(|i| source.eval(corner_point(b, i), cx).signed_distance);
        // Probe the fit at the center + face centers.
        let mid = Point3::new(
            f64::midpoint(b.min.x, b.max.x),
            f64::midpoint(b.min.y, b.max.y),
            f64::midpoint(b.min.z, b.max.z),
        );
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
            let truth = source.eval(p, cx).signed_distance;
            worst = worst.max((fit - truth).abs());
        }
        if worst <= tol || depth >= max_depth {
            *cells += 1;
            *depth_seen = (*depth_seen).max(depth);
            *residual = residual.max(worst);
            return Ok(Node::Leaf { corners });
        }
        let mut children: Vec<Node> = Vec::with_capacity(8);
        for i in 0..8 {
            children.push(Self::build_node(
                source,
                &octant(b, i),
                tol,
                max_depth,
                depth + 1,
                cx,
                cells,
                depth_seen,
                residual,
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
                    let mid = Point3::new(
                        f64::midpoint(b.min.x, b.max.x),
                        f64::midpoint(b.min.y, b.max.y),
                        f64::midpoint(b.min.z, b.max.z),
                    );
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
        // Estimate-grade: the residual is probed, not enclosed (fs-ivl
        // integration promotes this to a rigorous certificate later).
        let band = self.residual.max(self.tol);
        ChartSample {
            signed_distance: v,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::estimate(
                v - band - (1.0 + self.source_lipschitz) * dist_out,
                v + band,
            ),
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
