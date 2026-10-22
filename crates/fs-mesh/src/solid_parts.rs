//! Closed-surface generators for solid parts and their exact surfaces
//! (IRP M1.3-W4).
//!
//! Two generic generators produce a closed, outward-oriented PLC for
//! [`crate::volumetricize`] together with the exact surface it approximates:
//! - [`profile_bar`]: a bar of declared length and width whose depth follows a
//!   piecewise-cubic Hermite profile `h(x)`, with circular through-holes
//!   across the width;
//! - [`outline_solid`]: a planar outline with polygonal holes, extruded by a
//!   uniform or affine thickness field.
//!
//! Both mesh one planar section with [`crate::planar::mesh_outline`] and
//! extrude it; there is no second surface mesher. Curves are sampled to a
//! declared chord tolerance, a geometry axis of every ladder built on the
//! part. Section coordinates are rounded to a declared power-of-two grid so
//! the planar mesher's boundary midpoints stay exactly collinear (the same
//! convention as fs-plate's dyadic outlines); the receipt records the largest
//! rounding, and on-grid input is unchanged. Every part declares its
//! reference state; a missing one refuses.
//!
//! Frames. Bar: `x` along the length in `[0, L]`, `y` across the width in
//! `[0, W]`, `z` through the depth with the flat face at `z = 0` and the
//! profiled face at `z = h(x)`; hole axes run along `y`. Outline solid: the
//! outline lies in `(x, y)` and the solid spans `0 <= z <= t(x, y)`.
//!
//! Volumetric recovery inserts Steiner points and re-triangulates facets, so
//! curved patches are re-bound through the parent-facet index carried by
//! `source_faces`, never by vertex keys. [`SolidPart::snap`] moves every
//! boundary vertex of a tet complex onto its exact patch;
//! [`SolidPart::refine_snapped`] applies the uniform 1→8 split of
//! [`LabeledTetComplex::refine_uniform`] to the snapped geometry and snaps
//! the new boundary vertices; [`SolidPart::midside_projection`] supplies the
//! curved second-order boundary nodes for TET10 promotion.

use crate::LabeledTetComplex;
use crate::planar::{Loop, PlanarBudget, PlanarError, PlanarMesh, PlanarOptions, PlanarOutline};
use crate::tet_improve::{ImproveError, ImproveOptions, ImproveStats, improve_tets};
use crate::{RegionId, RegionKind, RegionSpec, UnverifiedPlc};
use fs_exec::Cx;
pub use fs_geom::ReferenceState;
use fs_ivl::{Sign, orient3d};
use std::collections::{BTreeMap, BTreeSet};

/// One knot of a piecewise-cubic Hermite depth profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HermiteKnot {
    /// Position along the length [m].
    pub x_m: f64,
    /// Depth at the knot [m].
    pub depth_m: f64,
    /// Depth slope `dh/dx` at the knot [-].
    pub slope: f64,
}

/// `h(x)`: C1 piecewise-cubic Hermite through the knots. The first knot sits
/// at `x = 0` and the last at the bar length.
#[derive(Debug, Clone, PartialEq)]
pub struct DepthProfile {
    /// Strictly increasing knots.
    pub knots: Vec<HermiteKnot>,
}

/// A circular through-hole whose axis runs across the width (along `y`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThroughHole {
    /// Axis position along the length [m].
    pub x_m: f64,
    /// Axis height above the flat face [m].
    pub z_m: f64,
    /// Radius [m].
    pub radius_m: f64,
}

/// A bar described by parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileBarSpec {
    /// Length [m].
    pub length_m: f64,
    /// Width [m].
    pub width_m: f64,
    /// Depth profile `h(x)`.
    pub depth_profile: DepthProfile,
    /// Holes across the width; each must lie strictly inside the material.
    pub through_holes: Vec<ThroughHole>,
    /// Mechanical state the geometry describes; `None` refuses.
    pub reference_state: Option<ReferenceState>,
}

/// Thickness of an extruded outline [m].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThicknessField {
    /// One thickness everywhere.
    Uniform {
        /// Thickness [m].
        thickness_m: f64,
    },
    /// `t(x, y) = at_origin_m + gradient[0] x + gradient[1] y`.
    Affine {
        /// Thickness at the origin [m].
        at_origin_m: f64,
        /// `(dt/dx, dt/dy)` [-].
        gradient: [f64; 2],
    },
}

impl ThicknessField {
    fn at(self, p: [f64; 2]) -> f64 {
        match self {
            Self::Uniform { thickness_m } => thickness_m,
            Self::Affine {
                at_origin_m,
                gradient,
            } => at_origin_m + gradient[0] * p[0] + gradient[1] * p[1],
        }
    }
}

/// A planar outline with holes, extruded along `z`.
#[derive(Debug, Clone, PartialEq)]
pub struct OutlineSolidSpec {
    /// Outer polygon in `(x, y)` [m]; the closing segment is implicit.
    pub outline: Vec<[f64; 2]>,
    /// Hole polygons [m].
    pub holes: Vec<Vec<[f64; 2]>>,
    /// Thickness field; positive at every outline vertex.
    pub thickness: ThicknessField,
    /// Mechanical state the geometry describes; `None` refuses.
    pub reference_state: Option<ReferenceState>,
}

/// Discretization choices and caps of one surface build.
#[derive(Debug, Clone, Copy)]
pub struct SolidPartOptions {
    /// Largest sagitta between a sampled curve and its chords [m].
    pub chord_tolerance_m: f64,
    /// Largest planar edge, also the target wall-layer spacing [m].
    pub maximum_edge_m: f64,
    /// Planar quality target in `[0, 20.7]` degrees.
    pub minimum_angle_deg: f64,
    /// Power-of-two grid [m] the section coordinates are rounded to. Each
    /// halving of a boundary segment consumes one bit below it.
    pub coordinate_grid_m: f64,
    /// Planar mesher caps.
    pub planar_budget: PlanarBudget,
    /// Largest number of closed-surface triangles.
    pub max_surface_triangles: usize,
}

/// Which exact surface a facet approximates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SurfacePatch {
    /// A plane: the PLC is exact there.
    Plane,
    /// The bar's profiled face `z = h(x)`.
    DepthProfile,
    /// The cylinder of through-hole `k`.
    Hole(u32),
}

impl SurfacePatch {
    const fn is_curved(self) -> bool {
        !matches!(self, Self::Plane)
    }
}

/// Refusals of the generators and of the snapping/projection stages.
#[derive(Debug, Clone, PartialEq)]
pub enum SolidPartError {
    /// The part does not declare its reference state.
    ReferenceStateMissing,
    /// A named field is non-finite, non-positive or inconsistent.
    InvalidSpec(&'static str),
    /// The planar section refused (self-intersection, holes, budget, ...).
    Planar(PlanarError),
    /// A cap was exceeded.
    Budget {
        /// Which cap.
        what: &'static str,
        /// Requested count.
        requested: usize,
        /// Cap.
        maximum: usize,
    },
    /// The surface does not enclose positive volume.
    Orientation,
    /// An imported surface is not a closed, consistently oriented 2-manifold
    /// (a leak, a non-manifold edge or an orientation conflict).
    NotClosed {
        /// Stable reason.
        reason: &'static str,
    },
    /// A source-face parent is not a facet of this part.
    UnknownFacet {
        /// Parent facet index.
        facet: u32,
    },
    /// A boundary vertex or edge touches two different curved patches.
    AmbiguousPatch {
        /// Vertex index.
        vertex: u32,
    },
    /// Snapping changed or zeroed the orientation of a tet.
    SnapInverted {
        /// Tet index.
        tet: usize,
    },
    /// Two boundary edges share a midpoint but project differently.
    MidpointCollision,
    /// The quality pass after a snap refused (a tet stayed inverted).
    Improve(ImproveError),
    /// Cancellation was observed.
    Cancelled,
}

impl core::fmt::Display for SolidPartError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MESH-SOLID-PART: {self:?}")
    }
}

impl std::error::Error for SolidPartError {}

impl From<PlanarError> for SolidPartError {
    fn from(value: PlanarError) -> Self {
        match value {
            PlanarError::Cancelled => Self::Cancelled,
            other => Self::Planar(other),
        }
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidPartError> {
    cx.checkpoint().map_err(|_| SolidPartError::Cancelled)
}

fn positive(value: f64, what: &'static str) -> Result<f64, SolidPartError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(SolidPartError::InvalidSpec(what))
    }
}

/// The cubic Hermite interpolant at `t` in `[0, 1]` of an interval of width `dx`.
fn hermite(y0: f64, m0: f64, y1: f64, m1: f64, dx: f64, t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    (2.0 * t3 - 3.0 * t2 + 1.0) * y0
        + (t3 - 2.0 * t2 + t) * dx * m0
        + (-2.0 * t3 + 3.0 * t2) * y1
        + (t3 - t2) * dx * m1
}

impl DepthProfile {
    fn validate(&self, length_m: f64) -> Result<(), SolidPartError> {
        if self.knots.len() < 2 {
            return Err(SolidPartError::InvalidSpec("depth profile needs two knots"));
        }
        for k in &self.knots {
            if !(k.x_m.is_finite() && k.slope.is_finite()) {
                return Err(SolidPartError::InvalidSpec("depth profile knot"));
            }
            positive(k.depth_m, "depth profile knot depth")?;
        }
        if self.knots[0].x_m != 0.0 || self.knots[self.knots.len() - 1].x_m != length_m {
            return Err(SolidPartError::InvalidSpec(
                "depth profile must span exactly [0, length]",
            ));
        }
        if self.knots.windows(2).any(|w| w[1].x_m <= w[0].x_m) {
            return Err(SolidPartError::InvalidSpec(
                "depth profile knots must increase strictly",
            ));
        }
        if !(self.minimum(0.0, length_m) > 0.0) {
            return Err(SolidPartError::InvalidSpec(
                "depth profile must stay positive",
            ));
        }
        Ok(())
    }

    fn interval(&self, x: f64) -> usize {
        let i = self.knots.partition_point(|k| k.x_m <= x);
        i.clamp(1, self.knots.len() - 1) - 1
    }

    /// `h(x)` [m] on `[0, L]`.
    #[must_use]
    pub fn depth(&self, x_m: f64) -> f64 {
        let i = self.interval(x_m);
        let (a, b) = (self.knots[i], self.knots[i + 1]);
        let dx = b.x_m - a.x_m;
        hermite(
            a.depth_m,
            a.slope,
            b.depth_m,
            b.slope,
            dx,
            (x_m - a.x_m) / dx,
        )
    }

    /// `∫_0^L h dx` [m²], exact for the cubic pieces.
    #[must_use]
    pub fn integral(&self) -> f64 {
        self.knots
            .windows(2)
            .map(|w| {
                let dx = w[1].x_m - w[0].x_m;
                dx * (w[0].depth_m + w[1].depth_m) / 2.0
                    + dx * dx * (w[0].slope - w[1].slope) / 12.0
            })
            .sum()
    }

    /// Bound on `|h''|` over interval `i` (linear in `t`, so its endpoints).
    fn curvature_bound(&self, i: usize) -> f64 {
        let (a, b) = (self.knots[i], self.knots[i + 1]);
        let dx = b.x_m - a.x_m;
        let at = |t: f64| {
            ((12.0 * t - 6.0) * a.depth_m
                + (6.0 * t - 4.0) * dx * a.slope
                + (-12.0 * t + 6.0) * b.depth_m
                + (6.0 * t - 2.0) * dx * b.slope)
                / (dx * dx)
        };
        at(0.0).abs().max(at(1.0).abs())
    }

    /// Exact minimum of `h` over `[lo, hi]`: interval ends and the roots of
    /// each cubic piece's derivative.
    fn minimum(&self, lo: f64, hi: f64) -> f64 {
        let mut best = self.depth(lo).min(self.depth(hi));
        for w in self.knots.windows(2) {
            let (a, b) = (w[0], w[1]);
            if b.x_m < lo || a.x_m > hi {
                continue;
            }
            for x in [a.x_m, b.x_m] {
                if (lo..=hi).contains(&x) {
                    best = best.min(self.depth(x));
                }
            }
            // h'(t) dx = A t² + B t + C.
            let dx = b.x_m - a.x_m;
            let (y0, m0, y1, m1) = (a.depth_m, a.slope * dx, b.depth_m, b.slope * dx);
            let qa = 6.0 * y0 + 3.0 * m0 - 6.0 * y1 + 3.0 * m1;
            let qb = -6.0 * y0 - 4.0 * m0 + 6.0 * y1 - 2.0 * m1;
            let qc = m0;
            let mut roots = Vec::new();
            if qa == 0.0 {
                if qb != 0.0 {
                    roots.push(-qc / qb);
                }
            } else {
                let disc = qb * qb - 4.0 * qa * qc;
                if disc >= 0.0 {
                    let s = fs_math::det::sqrt(disc);
                    roots.push((-qb - s) / (2.0 * qa));
                    roots.push((-qb + s) / (2.0 * qa));
                }
            }
            for t in roots {
                let x = a.x_m + t * dx;
                if (0.0..=1.0).contains(&t) && (lo..=hi).contains(&x) {
                    best = best.min(hermite(y0, a.slope, y1, b.slope, dx, t));
                }
            }
        }
        best
    }
}

/// The exact geometry a part's facets approximate.
#[derive(Debug, Clone, PartialEq)]
enum ExactSurface {
    /// Every facet is planar: the PLC is the exact surface.
    Polyhedral,
    /// The bar's profiled face and hole cylinders.
    ProfileBar {
        profile: DepthProfile,
        holes: Vec<ThroughHole>,
    },
}

impl ExactSurface {
    fn project(&self, patch: SurfacePatch, p: [f64; 3]) -> Result<[f64; 3], SolidPartError> {
        match (patch, self) {
            (SurfacePatch::Plane, _) => Ok(p),
            (SurfacePatch::DepthProfile, Self::ProfileBar { profile, .. }) => {
                Ok([p[0], p[1], profile.depth(p[0])])
            }
            (SurfacePatch::Hole(k), Self::ProfileBar { holes, .. }) => {
                let hole = holes
                    .get(k as usize)
                    .ok_or(SolidPartError::InvalidSpec("hole patch index"))?;
                let d = [p[0] - hole.x_m, p[2] - hole.z_m];
                let rho = fs_math::det::sqrt(d[0] * d[0] + d[1] * d[1]);
                if !(rho.is_finite() && rho > 0.0) {
                    return Err(SolidPartError::InvalidSpec("point on a hole axis"));
                }
                let s = hole.radius_m / rho;
                Ok([hole.x_m + s * d[0], p[1], hole.z_m + s * d[1]])
            }
            _ => Err(SolidPartError::InvalidSpec(
                "curved patch on a polyhedral part",
            )),
        }
    }
}

/// What a surface build produced and against what it is checked.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidPartReceipt {
    /// `"profile-bar"`, `"outline-solid"` or `"imported-tessellation"`.
    pub generator: &'static str,
    /// Declared reference state.
    pub reference_state: ReferenceState,
    /// Closed-surface vertices and triangles.
    pub vertices: usize,
    /// Closed-surface triangles.
    pub triangles: usize,
    /// `V - E + F` of the closed surface.
    pub euler_characteristic: i64,
    /// `(2 - χ) / 2`: handles of the closed surface.
    pub genus: i64,
    /// Declared chord tolerance [m] (an import: the declared tessellation
    /// deviation).
    pub chord_tolerance_m: f64,
    /// Declared section coordinate grid [m]; 0 for an import.
    pub coordinate_grid_m: f64,
    /// Largest coordinate change the grid rounding applied [m].
    pub max_quantization_m: f64,
    /// Wall layers per boundary loop (outer first).
    pub wall_layers: Vec<usize>,
    /// Smallest planar-section angle [degrees]; 0 for an import.
    pub planar_minimum_angle_deg: f64,
    /// Divergence-theorem volume of the PLC [m³].
    pub polyhedral_volume_m3: f64,
    /// Volume of the exact solid [m³]; an import's surface of record is its
    /// tessellation, so there it equals the polyhedral volume.
    pub exact_volume_m3: f64,
}

impl SolidPartReceipt {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"generator\":\"{}\",\"reference_state\":\"{}\",\"vertices\":{},\
             \"triangles\":{},\"euler_characteristic\":{},\"genus\":{},\
             \"chord_tolerance_m\":{:e},\"coordinate_grid_m\":{:e},\
             \"max_quantization_m\":{:e},\"wall_layers\":{:?},\
             \"planar_minimum_angle_deg\":{:e},\"polyhedral_volume_m3\":{:e},\
             \"exact_volume_m3\":{:e}}}",
            self.generator,
            self.reference_state.name(),
            self.vertices,
            self.triangles,
            self.euler_characteristic,
            self.genus,
            self.chord_tolerance_m,
            self.coordinate_grid_m,
            self.max_quantization_m,
            self.wall_layers,
            self.planar_minimum_angle_deg,
            self.polyhedral_volume_m3,
            self.exact_volume_m3
        )
    }
}

/// A generated closed surface, its exact geometry and its declarations.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidPart {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    triangle_patch: Vec<SurfacePatch>,
    /// Patch of each unique facet, in [`crate::volumetricize`]'s parent-facet
    /// order (sorted vertex triples).
    facet_patch: Vec<SurfacePatch>,
    seed: [f64; 3],
    surface: ExactSurface,
    receipt: SolidPartReceipt,
}

/// A tet mesh whose boundary vertices lie on a part's exact surface.
#[derive(Debug, Clone, PartialEq)]
pub struct SnappedTetMesh {
    /// Vertex positions [m].
    pub positions: Vec<[f64; 3]>,
    /// Tets, oriented as in the source complex.
    pub tets: Vec<[u32; 4]>,
    /// Boundary faces with their parent PLC facet.
    pub boundary_faces: Vec<([u32; 3], u32)>,
    /// Largest displacement the latest snap applied [m].
    pub max_snap_m: f64,
    /// Sign of `dot(b-a, (c-a) x (d-a))` every tet carries.
    pub sign: f64,
    /// What the quality pass after the latest snap did, when one ran.
    pub improvement: Option<ImproveStats>,
    /// Boundary slivers the latest quality pass peeled.
    pub peeled: usize,
}

/// Exact-surface positions of boundary-edge midside nodes, keyed by the
/// bits of the straight midpoint `f64::midpoint(p[a], p[b])` of the sorted
/// endpoint pair — the point a P1→P2 promotion asks to project.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MidsideProjection {
    targets: BTreeMap<[u64; 3], [f64; 3]>,
    curved: usize,
}

impl MidsideProjection {
    /// The projected position of a boundary-edge midpoint, or `None` when
    /// `p` is not one.
    #[must_use]
    pub fn get(&self, p: [f64; 3]) -> Option<[f64; 3]> {
        self.targets.get(&p.map(f64::to_bits)).copied()
    }

    /// Boundary edges tabulated.
    #[must_use]
    pub fn len(&self) -> usize {
        self.targets.len()
    }

    /// True when no edge is tabulated.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.targets.is_empty()
    }

    /// Edges whose midside node moved onto a curved patch.
    #[must_use]
    pub fn curved_edges(&self) -> usize {
        self.curved
    }
}

impl SolidPart {
    /// Closed-surface vertices [m].
    #[must_use]
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }

    /// Outward-oriented closed-surface triangles.
    #[must_use]
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// Exact patch of each triangle.
    #[must_use]
    pub fn triangle_patches(&self) -> &[SurfacePatch] {
        &self.triangle_patch
    }

    /// Strictly interior seed [m].
    #[must_use]
    pub const fn seed(&self) -> [f64; 3] {
        self.seed
    }

    /// Build facts.
    #[must_use]
    pub const fn receipt(&self) -> &SolidPartReceipt {
        &self.receipt
    }

    /// One solid region holding the whole surface.
    #[must_use]
    pub fn plc(&self, region: RegionId) -> UnverifiedPlc {
        UnverifiedPlc::new(
            self.vertices.clone(),
            vec![RegionSpec {
                id: region,
                kind: RegionKind::Solid,
                seed: self.seed,
                triangles: self.triangles.clone(),
            }],
        )
    }

    /// Project `p` onto `patch` of the exact surface.
    ///
    /// # Errors
    /// A curved patch this part does not have, or a point on a hole axis.
    pub fn project(&self, patch: SurfacePatch, p: [f64; 3]) -> Result<[f64; 3], SolidPartError> {
        self.surface.project(patch, p)
    }

    fn patch_of(&self, parent: u32) -> Result<SurfacePatch, SolidPartError> {
        self.facet_patch
            .get(parent as usize)
            .copied()
            .ok_or(SolidPartError::UnknownFacet { facet: parent })
    }

    /// Curved patch of each boundary vertex at or above `first`.
    fn curved_vertices(
        &self,
        faces: &[([u32; 3], u32)],
        first: u32,
    ) -> Result<BTreeMap<u32, SurfacePatch>, SolidPartError> {
        let mut out = BTreeMap::new();
        for &(face, parent) in faces {
            let patch = self.patch_of(parent)?;
            if !patch.is_curved() {
                continue;
            }
            for v in face {
                if v < first {
                    continue;
                }
                if let Some(old) = out.insert(v, patch) {
                    if old != patch {
                        return Err(SolidPartError::AmbiguousPatch { vertex: v });
                    }
                }
            }
        }
        Ok(out)
    }

    /// Snap `targets`, then either run the quality pass (boundary frozen) or,
    /// without one, refuse any tet the snap inverted.
    #[allow(clippy::too_many_arguments)]
    fn snap_vertices(
        &self,
        positions: &mut [[f64; 3]],
        tets: &mut Vec<[u32; 4]>,
        faces: &mut Vec<([u32; 3], u32)>,
        targets: &BTreeMap<u32, SurfacePatch>,
        sign: f64,
        improve: Option<ImproveOptions>,
        cx: &Cx<'_>,
    ) -> Result<(f64, Option<ImproveStats>, usize), SolidPartError> {
        // Only tets touching a moved vertex can change orientation.
        let mut moved = vec![false; positions.len()];
        let mut new_positions = Vec::with_capacity(targets.len());
        let mut max_snap = 0.0f64;
        for (i, (&v, &patch)) in targets.iter().enumerate() {
            if i % 1024 == 0 {
                checkpoint(cx)?;
            }
            let old = positions[v as usize];
            let new = self.surface.project(patch, old)?;
            let d = [new[0] - old[0], new[1] - old[1], new[2] - old[2]];
            max_snap = max_snap.max(fs_math::det::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
            if new.map(f64::to_bits) != old.map(f64::to_bits) {
                moved[v as usize] = true;
                new_positions.push((v, new));
            }
        }
        let touched: Vec<usize> = tets
            .iter()
            .enumerate()
            .filter(|(_, t)| t.iter().any(|&v| moved[v as usize]))
            .map(|(i, _)| i)
            .collect();
        let before: Vec<Sign> = touched
            .iter()
            .map(|&i| orient3d_of(positions, tets[i]))
            .collect();
        for (v, p) in new_positions {
            positions[v as usize] = p;
        }
        if let Some(options) = improve {
            let peeled = self.peel(
                positions,
                tets,
                faces,
                sign,
                options.target_min_dihedral_deg,
                cx,
            )?;
            let mut frozen = vec![false; positions.len()];
            for (face, _) in faces.iter() {
                for &v in face {
                    frozen[v as usize] = true;
                }
            }
            let boundary: Vec<[u32; 3]> = faces.iter().map(|(f, _)| *f).collect();
            let stats = improve_tets(positions, tets, &frozen, &boundary, sign, options, cx)
                .map_err(|e| match e {
                    ImproveError::Cancelled => SolidPartError::Cancelled,
                    other => SolidPartError::Improve(other),
                })?;
            return Ok((max_snap, Some(stats), peeled));
        }
        for (k, (&tet, sign)) in touched.iter().zip(before).enumerate() {
            if k % 4096 == 0 {
                checkpoint(cx)?;
            }
            let after = orient3d_of(positions, tets[tet]);
            if after == Sign::Zero || after != sign {
                return Err(SolidPartError::SnapInverted { tet });
            }
        }
        Ok((max_snap, None, 0))
    }

    /// Peel boundary slivers: a tet below `target_deg` (or inverted) whose
    /// only boundary faces are two faces of one exact patch sharing an edge
    /// no other tet uses, while its opposite edge is not on the boundary. Its
    /// two interior faces become boundary faces of that patch, so every
    /// boundary vertex stays on the exact surface; the boundary triangulation
    /// changes locally and the topology does not. Repeats until none
    /// qualifies; one round never peels two tets sharing a vertex.
    fn peel(
        &self,
        positions: &[[f64; 3]],
        tets: &mut Vec<[u32; 4]>,
        faces: &mut Vec<([u32; 3], u32)>,
        sign: f64,
        target_deg: f64,
        cx: &Cx<'_>,
    ) -> Result<usize, SolidPartError> {
        let sorted = |mut f: [u32; 3]| {
            f.sort_unstable();
            f
        };
        let mut peeled = 0;
        loop {
            checkpoint(cx)?;
            let boundary: BTreeMap<[u32; 3], u32> = faces
                .iter()
                .map(|(f, parent)| (sorted(*f), *parent))
                .collect();
            let boundary_edges: BTreeSet<[u32; 2]> = boundary
                .keys()
                .flat_map(|f| [[f[0], f[1]], [f[0], f[2]], [f[1], f[2]]])
                .collect();
            let mut incident = vec![Vec::new(); positions.len()];
            for (i, t) in tets.iter().enumerate() {
                for &v in t {
                    incident[v as usize].push(i);
                }
            }
            let mut touched = BTreeSet::new();
            let mut remove = BTreeSet::new();
            let mut dropped = BTreeSet::new();
            let mut added = Vec::new();
            for (i, t) in tets.iter().enumerate() {
                if i % 4096 == 0 {
                    checkpoint(cx)?;
                }
                let p = t.map(|v| positions[v as usize]);
                let volume =
                    sign * dot3(sub3(p[1], p[0]), cross3(sub3(p[2], p[0]), sub3(p[3], p[0])));
                if volume > 0.0 && crate::exude::min_dihedral_deg(&p) >= target_deg {
                    continue;
                }
                let on: Vec<([u32; 3], u32)> = [
                    [t[1], t[2], t[3]],
                    [t[0], t[2], t[3]],
                    [t[0], t[1], t[3]],
                    [t[0], t[1], t[2]],
                ]
                .iter()
                .filter_map(|f| {
                    boundary
                        .get(&sorted(*f))
                        .map(|parent| (sorted(*f), *parent))
                })
                .collect();
                if on.len() != 2 || self.patch_of(on[0].1)? != self.patch_of(on[1].1)? {
                    continue;
                }
                let ab: Vec<u32> = on[0]
                    .0
                    .iter()
                    .copied()
                    .filter(|v| on[1].0.contains(v))
                    .collect();
                let cd: Vec<u32> = t.iter().copied().filter(|v| !ab.contains(v)).collect();
                let (&[a, b], &[c, d]) = (ab.as_slice(), cd.as_slice()) else {
                    continue;
                };
                if boundary_edges.contains(&[c.min(d), c.max(d)])
                    || incident[a as usize]
                        .iter()
                        .filter(|&&u| tets[u].contains(&b))
                        .count()
                        != 1
                    || t.iter().any(|v| touched.contains(v))
                {
                    continue;
                }
                touched.extend(t.iter().copied());
                remove.insert(i);
                dropped.insert(on[0].0);
                dropped.insert(on[1].0);
                added.push((sorted([a, c, d]), on[0].1));
                added.push((sorted([b, c, d]), on[1].1));
            }
            if remove.is_empty() {
                return Ok(peeled);
            }
            peeled += remove.len();
            let mut index = 0;
            tets.retain(|_| {
                let keep = !remove.contains(&index);
                index += 1;
                keep
            });
            faces.retain(|(f, _)| !dropped.contains(&sorted(*f)));
            faces.extend(added);
        }
    }

    /// Move every boundary vertex of `complex` (a volumetricization of this
    /// part's PLC) onto its exact patch. Planar patches are left untouched.
    /// With `improve`, a quality pass that moves only interior vertices and
    /// flips only interior faces follows (volumes depend on the boundary
    /// alone, so it cannot change them); without one, any inversion refuses.
    ///
    /// # Errors
    /// Unknown parent facets, a vertex on two curved patches, a tet whose
    /// orientation the snap flips or zeroes, cancellation.
    pub fn snap(
        &self,
        complex: &LabeledTetComplex,
        improve: Option<ImproveOptions>,
        cx: &Cx<'_>,
    ) -> Result<SnappedTetMesh, SolidPartError> {
        checkpoint(cx)?;
        let mut positions = complex.positions().to_vec();
        let mut tets = complex.tets().to_vec();
        let mut boundary_faces = complex.source_faces().to_vec();
        let sign = tets
            .iter()
            .map(|t| {
                let [a, b, c, d] = t.map(|v| positions[v as usize]);
                dot3(sub3(b, a), cross3(sub3(c, a), sub3(d, a))).signum()
            })
            .find(|s| *s != 0.0)
            .ok_or(SolidPartError::InvalidSpec("complex has no tet volume"))?;
        let targets = self.curved_vertices(&boundary_faces, 0)?;
        let (max_snap_m, improvement, peeled) = self.snap_vertices(
            &mut positions,
            &mut tets,
            &mut boundary_faces,
            &targets,
            sign,
            improve,
            cx,
        )?;
        Ok(SnappedTetMesh {
            positions,
            tets,
            boundary_faces,
            max_snap_m,
            sign,
            improvement,
            peeled,
        })
    }

    /// One uniform 1→8 rung of a snapped mesh: the split of
    /// [`LabeledTetComplex::refine_uniform`], then the new boundary vertices
    /// snapped. Existing vertices keep their bits.
    ///
    /// # Errors
    /// As [`Self::snap`], plus the tet cap.
    pub fn refine_snapped(
        &self,
        mesh: &SnappedTetMesh,
        max_tets: usize,
        improve: Option<ImproveOptions>,
        cx: &Cx<'_>,
    ) -> Result<SnappedTetMesh, SolidPartError> {
        checkpoint(cx)?;
        let requested = mesh.tets.len().saturating_mul(8);
        if requested > max_tets {
            return Err(SolidPartError::Budget {
                what: "refined tets",
                requested,
                maximum: max_tets,
            });
        }
        let split =
            crate::uniform::split_uniform(&mesh.positions, &mesh.tets, &mesh.boundary_faces);
        let first = u32::try_from(mesh.positions.len()).map_err(|_| SolidPartError::Budget {
            what: "vertices",
            requested: mesh.positions.len(),
            maximum: u32::MAX as usize,
        })?;
        let targets = self.curved_vertices(&split.source_faces, first)?;
        let mut positions = split.positions;
        let mut tets = split.tets;
        let mut boundary_faces = split.source_faces;
        let (max_snap_m, improvement, peeled) = self.snap_vertices(
            &mut positions,
            &mut tets,
            &mut boundary_faces,
            &targets,
            mesh.sign,
            improve,
            cx,
        )?;
        Ok(SnappedTetMesh {
            positions,
            tets,
            boundary_faces,
            max_snap_m,
            sign: mesh.sign,
            improvement,
            peeled,
        })
    }

    /// Exact-surface midside nodes of every boundary edge of `mesh`. An edge
    /// on a curved patch projects its straight midpoint onto that patch; an
    /// edge on planes only keeps the midpoint.
    ///
    /// # Errors
    /// Unknown parent facets, an edge on two curved patches, two edges with
    /// one midpoint and different targets, cancellation.
    pub fn midside_projection(
        &self,
        mesh: &SnappedTetMesh,
        cx: &Cx<'_>,
    ) -> Result<MidsideProjection, SolidPartError> {
        checkpoint(cx)?;
        let mut edges: BTreeMap<[u32; 2], Option<SurfacePatch>> = BTreeMap::new();
        for &(face, parent) in &mesh.boundary_faces {
            let patch = self.patch_of(parent)?;
            for (a, b) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                let slot = edges.entry([a.min(b), a.max(b)]).or_insert(None);
                if patch.is_curved() {
                    if slot.is_some_and(|old| old != patch) {
                        return Err(SolidPartError::AmbiguousPatch { vertex: a.min(b) });
                    }
                    *slot = Some(patch);
                }
            }
        }
        let mut out = MidsideProjection::default();
        for (i, ([a, b], patch)) in edges.into_iter().enumerate() {
            if i % 4096 == 0 {
                checkpoint(cx)?;
            }
            let (pa, pb) = (mesh.positions[a as usize], mesh.positions[b as usize]);
            let mid: [f64; 3] = core::array::from_fn(|k| f64::midpoint(pa[k], pb[k]));
            let target = match patch {
                Some(patch) => {
                    out.curved += 1;
                    self.surface.project(patch, mid)?
                }
                None => mid,
            };
            if let Some(old) = out.targets.insert(mid.map(f64::to_bits), target) {
                if old.map(f64::to_bits) != target.map(f64::to_bits) {
                    return Err(SolidPartError::MidpointCollision);
                }
            }
        }
        Ok(out)
    }
}

fn orient3d_of(positions: &[[f64; 3]], t: [u32; 4]) -> Sign {
    let p = |v: u32| positions[v as usize];
    orient3d(p(t[0]), p(t[1]), p(t[2]), p(t[3]))
}

/// The closed surface of an extruded planar section.
struct Extrusion {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    patches: Vec<SurfacePatch>,
    layers: Vec<usize>,
    seed: [f64; 3],
}

fn to_u32(n: usize) -> Result<u32, SolidPartError> {
    u32::try_from(n).map_err(|_| SolidPartError::Budget {
        what: "surface vertices",
        requested: n,
        maximum: u32::MAX as usize,
    })
}

/// Extrude `mesh` (local `(u, v)`) along local `w` from `0` to
/// `height(u, v)`, then map local `(u, v, w)` through the right-handed
/// `frame`. Source segment `s` of the section (in mesh_outline's original
/// order) carries `segment_patch[s]` onto its wall; the two caps are planes.
/// Each boundary loop gets `ceil(max height / mean loop edge)` wall layers.
#[allow(clippy::too_many_arguments)]
fn extrude(
    mesh: &PlanarMesh,
    segment_patch: &[SurfacePatch],
    segment_loop: &[usize],
    height: &dyn Fn([f64; 2]) -> f64,
    maximum_edge_m: f64,
    frame: &dyn Fn([f64; 3]) -> [f64; 3],
    max_triangles: usize,
    cx: &Cx<'_>,
) -> Result<Extrusion, SolidPartError> {
    checkpoint(cx)?;
    let chains = &mesh.receipt.source_segment_chains;
    if chains.len() != segment_patch.len() || chains.len() != segment_loop.len() {
        return Err(SolidPartError::InvalidSpec("segment table mismatch"));
    }
    let n2 = mesh.vertices.len();
    // Counter-clockwise section triangles, and the triangle owning each
    // directed edge.
    let mut tris = Vec::with_capacity(mesh.triangles.len());
    for t in &mesh.triangles {
        let p = t.map(|i| mesh.vertices[i]);
        let cross =
            (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]);
        if cross > 0.0 {
            tris.push(*t);
        } else if cross < 0.0 {
            tris.push([t[0], t[2], t[1]]);
        } else {
            return Err(SolidPartError::InvalidSpec("degenerate section triangle"));
        }
    }
    let directed: BTreeSet<(usize, usize)> = tris
        .iter()
        .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
        .collect();
    let loops = segment_loop.iter().copied().max().map_or(0, |m| m + 1);
    let mut loop_length = vec![0.0f64; loops];
    let mut loop_edges = vec![0usize; loops];
    let mut loop_height = vec![0.0f64; loops];
    let mut loop_vertices: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); loops];
    for (s, chain) in chains.iter().enumerate() {
        let r = segment_loop[s];
        for w in chain.windows(2) {
            let (a, b) = (mesh.vertices[w[0]], mesh.vertices[w[1]]);
            loop_length[r] += fs_math::det::hypot(b[0] - a[0], b[1] - a[1]);
            loop_edges[r] += 1;
        }
        for &v in chain {
            loop_height[r] = loop_height[r].max(height(mesh.vertices[v]));
            loop_vertices[r].insert(v);
        }
    }
    let mut layers = Vec::with_capacity(loops);
    for r in 0..loops {
        if loop_edges[r] == 0 {
            return Err(SolidPartError::InvalidSpec("empty boundary loop"));
        }
        let spacing = (loop_length[r] / loop_edges[r] as f64).min(maximum_edge_m);
        let n = (loop_height[r] / spacing).ceil();
        if !(n.is_finite() && n >= 1.0 && n <= max_triangles as f64) {
            return Err(SolidPartError::Budget {
                what: "wall layers",
                requested: max_triangles.saturating_add(1),
                maximum: max_triangles,
            });
        }
        layers.push(n as usize);
    }
    let requested = loop_edges
        .iter()
        .zip(&layers)
        .try_fold(2usize.saturating_mul(tris.len()), |acc, (e, n)| {
            acc.checked_add(e.checked_mul(*n)?.checked_mul(2)?)
        })
        .unwrap_or(usize::MAX);
    if requested > max_triangles {
        return Err(SolidPartError::Budget {
            what: "surface triangles",
            requested,
            maximum: max_triangles,
        });
    }
    let local = |v: usize, k: usize, n: usize| {
        let p = mesh.vertices[v];
        let h = height(p);
        let w = if k == 0 {
            0.0
        } else if k == n {
            h
        } else {
            h * k as f64 / n as f64
        };
        frame([p[0], p[1], w])
    };
    let mut vertices = Vec::with_capacity(2 * n2);
    for v in 0..n2 {
        vertices.push(local(v, 0, 1));
    }
    for v in 0..n2 {
        vertices.push(local(v, 1, 1));
    }
    // Intermediate wall vertices, loop by loop, vertex by vertex.
    let mut column: BTreeMap<(usize, usize), u32> = BTreeMap::new();
    for r in 0..loops {
        for &v in &loop_vertices[r] {
            for k in 1..layers[r] {
                column.insert((v, k), to_u32(vertices.len())?);
                vertices.push(local(v, k, layers[r]));
            }
        }
    }
    let at = |v: usize, k: usize, n: usize| -> Result<u32, SolidPartError> {
        if k == 0 {
            to_u32(v)
        } else if k == n {
            to_u32(n2 + v)
        } else {
            column
                .get(&(v, k))
                .copied()
                .ok_or(SolidPartError::InvalidSpec("missing wall column"))
        }
    };
    let mut triangles = Vec::with_capacity(requested);
    let mut patches = Vec::with_capacity(requested);
    for t in &tris {
        triangles.push([to_u32(t[0])?, to_u32(t[2])?, to_u32(t[1])?]);
        patches.push(SurfacePatch::Plane);
        triangles.push([to_u32(n2 + t[0])?, to_u32(n2 + t[1])?, to_u32(n2 + t[2])?]);
        patches.push(SurfacePatch::Plane);
    }
    for (s, chain) in chains.iter().enumerate() {
        checkpoint(cx)?;
        let n = layers[segment_loop[s]];
        for w in chain.windows(2) {
            // Material lies left of p→q in the counter-clockwise triangle.
            let (p, q) = if directed.contains(&(w[0], w[1])) {
                (w[0], w[1])
            } else if directed.contains(&(w[1], w[0])) {
                (w[1], w[0])
            } else {
                return Err(SolidPartError::InvalidSpec(
                    "boundary segment has no section triangle",
                ));
            };
            for k in 0..n {
                let (a, b, c, d) = (
                    at(p, k, n)?,
                    at(q, k, n)?,
                    at(q, k + 1, n)?,
                    at(p, k + 1, n)?,
                );
                triangles.push([a, b, c]);
                triangles.push([a, c, d]);
                patches.push(segment_patch[s]);
                patches.push(segment_patch[s]);
            }
        }
    }
    // Seed: centroid of the largest section triangle, at a height fraction no
    // wall layer uses.
    let largest = tris
        .iter()
        .max_by(|a, b| {
            section_area(mesh, **a)
                .partial_cmp(&section_area(mesh, **b))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .ok_or(SolidPartError::InvalidSpec("empty section"))?;
    let c = [0, 1].map(|k| largest.iter().map(|&v| mesh.vertices[v][k]).sum::<f64>() / 3.0);
    let seed = frame([c[0], c[1], height(c) * (core::f64::consts::SQRT_2 - 1.0)]);
    Ok(Extrusion {
        vertices,
        triangles,
        patches,
        layers,
        seed,
    })
}

fn section_area(mesh: &PlanarMesh, t: [usize; 3]) -> f64 {
    let p = t.map(|i| mesh.vertices[i]);
    ((p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0])).abs()
        / 2.0
}

/// Signed divergence-theorem volume of a closed triangle surface.
fn enclosed_volume(vertices: &[[f64; 3]], triangles: &[[u32; 3]]) -> f64 {
    triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| vertices[i as usize]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

fn euler_characteristic(triangles: &[[u32; 3]]) -> i64 {
    let vertices: BTreeSet<u32> = triangles.iter().flatten().copied().collect();
    let edges: BTreeSet<[u32; 2]> = triangles
        .iter()
        .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
        .map(|(a, b)| [a.min(b), a.max(b)])
        .collect();
    vertices.len() as i64 - edges.len() as i64 + triangles.len() as i64
}

/// Discretization facts a build hands to its receipt.
struct BuildFacts {
    chord_tolerance_m: f64,
    coordinate_grid_m: f64,
    planar_minimum_angle_deg: f64,
    max_quantization_m: f64,
    /// `None`: the surface is its own exact geometry.
    exact_volume_m3: Option<f64>,
}

fn assemble_part(
    generator: &'static str,
    reference_state: ReferenceState,
    extrusion: Extrusion,
    surface: ExactSurface,
    facts: BuildFacts,
) -> Result<SolidPart, SolidPartError> {
    let polyhedral_volume_m3 = enclosed_volume(&extrusion.vertices, &extrusion.triangles);
    if !(polyhedral_volume_m3.is_finite() && polyhedral_volume_m3 > 0.0) {
        return Err(SolidPartError::Orientation);
    }
    let chi = euler_characteristic(&extrusion.triangles);
    let mut facets: BTreeMap<[u32; 3], SurfacePatch> = BTreeMap::new();
    for (t, patch) in extrusion.triangles.iter().zip(&extrusion.patches) {
        let mut key = *t;
        key.sort_unstable();
        facets.insert(key, *patch);
    }
    Ok(SolidPart {
        receipt: SolidPartReceipt {
            generator,
            reference_state,
            vertices: extrusion.vertices.len(),
            triangles: extrusion.triangles.len(),
            euler_characteristic: chi,
            genus: (2 - chi) / 2,
            chord_tolerance_m: facts.chord_tolerance_m,
            coordinate_grid_m: facts.coordinate_grid_m,
            max_quantization_m: facts.max_quantization_m,
            wall_layers: extrusion.layers,
            planar_minimum_angle_deg: facts.planar_minimum_angle_deg,
            polyhedral_volume_m3,
            exact_volume_m3: facts.exact_volume_m3.unwrap_or(polyhedral_volume_m3),
        },
        facet_patch: facets.into_values().collect(),
        vertices: extrusion.vertices,
        triangles: extrusion.triangles,
        triangle_patch: extrusion.patches,
        seed: extrusion.seed,
        surface,
    })
}

fn validate_options(options: &SolidPartOptions) -> Result<(), SolidPartError> {
    positive(options.chord_tolerance_m, "chord tolerance")?;
    positive(options.maximum_edge_m, "maximum edge")?;
    let g = positive(options.coordinate_grid_m, "coordinate grid")?;
    // A normal power of two has an all-zero mantissa.
    if !g.is_normal() || g.to_bits() & ((1u64 << 52) - 1) != 0 || g >= options.maximum_edge_m {
        return Err(SolidPartError::InvalidSpec(
            "coordinate grid must be a power of two below the maximum edge",
        ));
    }
    Ok(())
}

/// Round every coordinate to the grid (exact: scaling by a power of two and
/// rounding are exact); returns the largest coordinate change.
fn quantize(rings: &mut [Vec<[f64; 2]>], grid: f64) -> f64 {
    let mut worst = 0.0f64;
    for p in rings.iter_mut().flatten() {
        for c in p.iter_mut() {
            let q = (*c / grid).round() * grid;
            worst = worst.max((q - *c).abs());
            *c = q;
        }
    }
    worst
}

fn planar_options(options: &SolidPartOptions) -> PlanarOptions<'static> {
    PlanarOptions {
        maximum_edge_m: options.maximum_edge_m,
        minimum_angle_deg: options.minimum_angle_deg,
        size_field: None,
        budget: options.planar_budget,
    }
}

/// Straight-segment subdivision count no longer than `maximum_edge_m`.
fn pieces(length: f64, maximum_edge_m: f64, cap: usize) -> Result<usize, SolidPartError> {
    let n = (length / maximum_edge_m).ceil().max(1.0);
    if !(n.is_finite() && n <= cap as f64) {
        return Err(SolidPartError::Budget {
            what: "outline vertices",
            requested: cap.saturating_add(1),
            maximum: cap,
        });
    }
    Ok(n as usize)
}

/// Mesh a bar from its description. The section `{0<=x<=L, 0<=z<=h(x)}`
/// minus the hole disks is meshed in `(x, z)` and extruded across the width.
///
/// # Errors
/// Missing reference state, invalid parameters (including holes that leave
/// the material or overlap), planar refusals, caps, cancellation.
pub fn profile_bar(
    spec: &ProfileBarSpec,
    options: &SolidPartOptions,
    cx: &Cx<'_>,
) -> Result<SolidPart, SolidPartError> {
    let reference_state = spec
        .reference_state
        .ok_or(SolidPartError::ReferenceStateMissing)?;
    checkpoint(cx)?;
    validate_options(options)?;
    let length = positive(spec.length_m, "bar length")?;
    let width = positive(spec.width_m, "bar width")?;
    let profile = &spec.depth_profile;
    profile.validate(length)?;
    let eps = options.chord_tolerance_m;
    let cap = options.planar_budget.max_vertices;
    for (k, hole) in spec.through_holes.iter().enumerate() {
        let r = positive(hole.radius_m, "hole radius")?;
        if !(hole.x_m.is_finite() && hole.z_m.is_finite()) {
            return Err(SolidPartError::InvalidSpec("hole centre"));
        }
        if !(hole.x_m - r > 0.0
            && hole.x_m + r < length
            && hole.z_m - r > 0.0
            && hole.z_m + r < profile.minimum(hole.x_m - r, hole.x_m + r))
        {
            return Err(SolidPartError::InvalidSpec(
                "hole must lie strictly inside the material",
            ));
        }
        for other in &spec.through_holes[..k] {
            let d = fs_math::det::hypot(hole.x_m - other.x_m, hole.z_m - other.z_m);
            if !(d > r + other.radius_m) {
                return Err(SolidPartError::InvalidSpec("holes overlap"));
            }
        }
    }
    // Outer section ring with the patch of the segment starting at each
    // point: flat face, far end, profiled face (right to left), near end.
    let mut ring: Vec<([f64; 2], SurfacePatch)> = Vec::new();
    let n = pieces(length, options.maximum_edge_m, cap)?;
    for i in 0..n {
        ring.push(([length * i as f64 / n as f64, 0.0], SurfacePatch::Plane));
    }
    let end_depth = profile.depth(length);
    let n = pieces(end_depth, options.maximum_edge_m, cap)?;
    for i in 0..n {
        ring.push((
            [length, end_depth * i as f64 / n as f64],
            SurfacePatch::Plane,
        ));
    }
    let mut top = Vec::new();
    for i in 0..profile.knots.len() - 1 {
        let (a, b) = (profile.knots[i].x_m, profile.knots[i + 1].x_m);
        let bound = profile.curvature_bound(i);
        let step = if bound > 0.0 {
            fs_math::det::sqrt(8.0 * eps / bound).min(options.maximum_edge_m)
        } else {
            options.maximum_edge_m
        };
        let n = pieces(b - a, step, cap)?;
        for j in 0..n {
            let x = if j == 0 {
                a
            } else {
                a + (b - a) * j as f64 / n as f64
            };
            top.push(x);
        }
    }
    top.push(length);
    for &x in top.iter().rev().take(top.len() - 1) {
        ring.push(([x, profile.depth(x)], SurfacePatch::DepthProfile));
    }
    let start_depth = profile.depth(0.0);
    let n = pieces(start_depth, options.maximum_edge_m, cap)?;
    for i in 0..n {
        ring.push((
            [0.0, start_depth * (n - i) as f64 / n as f64],
            SurfacePatch::Plane,
        ));
    }
    let mut hole_rings = Vec::with_capacity(spec.through_holes.len());
    for hole in &spec.through_holes {
        let r = hole.radius_m;
        let by_chord = core::f64::consts::PI / fs_math::det::acos((1.0 - eps / r).max(-1.0));
        let by_edge = 2.0 * core::f64::consts::PI * r / options.maximum_edge_m;
        let count = by_chord.max(by_edge).ceil().max(8.0);
        if !(count.is_finite() && count <= cap as f64) {
            return Err(SolidPartError::Budget {
                what: "hole vertices",
                requested: cap.saturating_add(1),
                maximum: cap,
            });
        }
        let count = count as usize;
        hole_rings.push(
            (0..count)
                .map(|j| {
                    let (s, c) = fs_math::det::sin_cos(
                        2.0 * core::f64::consts::PI * j as f64 / count as f64,
                    );
                    [hole.x_m + r * c, hole.z_m + r * s]
                })
                .collect::<Vec<_>>(),
        );
    }
    let mut outer = vec![ring.iter().map(|(p, _)| *p).collect::<Vec<[f64; 2]>>()];
    let mut max_quantization = quantize(&mut outer, options.coordinate_grid_m);
    max_quantization = max_quantization.max(quantize(&mut hole_rings, options.coordinate_grid_m));
    let outer = outer.remove(0);
    let mut segment_patch: Vec<SurfacePatch> = ring.iter().map(|(_, s)| *s).collect();
    let mut segment_loop = vec![0usize; ring.len()];
    for (k, h) in hole_rings.iter().enumerate() {
        let patch = SurfacePatch::Hole(
            u32::try_from(k).map_err(|_| SolidPartError::InvalidSpec("too many holes"))?,
        );
        segment_patch.extend(std::iter::repeat_n(patch, h.len()));
        segment_loop.extend(std::iter::repeat_n(k + 1, h.len()));
    }
    let section = crate::planar::mesh_outline(
        &PlanarOutline {
            outer: Loop::Polyline(&outer),
            holes: hole_rings.iter().map(|h| Loop::Polyline(h)).collect(),
        },
        &planar_options(options),
        cx,
    )?;
    // Local (u, v, w) = (x, z, W - y) is right-handed.
    let extrusion = extrude(
        &section,
        &segment_patch,
        &segment_loop,
        &|_| width,
        options.maximum_edge_m,
        &|p| [p[0], width - p[2], p[1]],
        options.max_surface_triangles,
        cx,
    )?;
    let hole_area: f64 = spec
        .through_holes
        .iter()
        .map(|h| core::f64::consts::PI * h.radius_m * h.radius_m)
        .sum();
    assemble_part(
        "profile-bar",
        reference_state,
        extrusion,
        ExactSurface::ProfileBar {
            profile: profile.clone(),
            holes: spec.through_holes.clone(),
        },
        BuildFacts {
            chord_tolerance_m: options.chord_tolerance_m,
            coordinate_grid_m: options.coordinate_grid_m,
            planar_minimum_angle_deg: section.receipt.minimum_angle_deg,
            max_quantization_m: max_quantization,
            exact_volume_m3: Some(width * (profile.integral() - hole_area)),
        },
    )
}

/// Signed area and first moments `(A, ∫x dA, ∫y dA)` of a ring (shoelace).
fn ring_moments(ring: &[[f64; 2]]) -> [f64; 3] {
    let mut m = [0.0; 3];
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        let cross = a[0] * b[1] - b[0] * a[1];
        m[0] += cross / 2.0;
        m[1] += (a[0] + b[0]) * cross / 6.0;
        m[2] += (a[1] + b[1]) * cross / 6.0;
    }
    m
}

/// Extrude a planar outline with holes by a thickness field. The exact
/// volume `∫ t dA` comes from the polygons' shoelace moments.
///
/// # Errors
/// Missing reference state, non-positive thickness at an outline vertex,
/// planar refusals (self-intersection, holes outside or overlapping), caps,
/// cancellation.
pub fn outline_solid(
    spec: &OutlineSolidSpec,
    options: &SolidPartOptions,
    cx: &Cx<'_>,
) -> Result<SolidPart, SolidPartError> {
    let reference_state = spec
        .reference_state
        .ok_or(SolidPartError::ReferenceStateMissing)?;
    checkpoint(cx)?;
    validate_options(options)?;
    match spec.thickness {
        ThicknessField::Uniform { thickness_m } => {
            positive(thickness_m, "thickness")?;
        }
        ThicknessField::Affine {
            at_origin_m,
            gradient,
        } => {
            if !(at_origin_m.is_finite() && gradient.iter().all(|g| g.is_finite())) {
                return Err(SolidPartError::InvalidSpec("thickness field"));
            }
        }
    }
    // Affine: the minimum over the outline's convex hull sits at an outline
    // vertex, and holes lie inside the outline.
    for p in &spec.outline {
        positive(spec.thickness.at(*p), "thickness at an outline vertex")?;
    }
    let mut rings = vec![spec.outline.clone()];
    rings.extend(spec.holes.iter().cloned());
    let max_quantization = quantize(&mut rings, options.coordinate_grid_m);
    let section = crate::planar::mesh_outline(
        &PlanarOutline {
            outer: Loop::Polyline(&rings[0]),
            holes: rings[1..].iter().map(|h| Loop::Polyline(h)).collect(),
        },
        &planar_options(options),
        cx,
    )?;
    let mut segment_loop = vec![0usize; spec.outline.len()];
    for (k, h) in spec.holes.iter().enumerate() {
        segment_loop.extend(std::iter::repeat_n(k + 1, h.len()));
    }
    let segment_patch = vec![SurfacePatch::Plane; segment_loop.len()];
    let thickness = spec.thickness;
    let extrusion = extrude(
        &section,
        &segment_patch,
        &segment_loop,
        &|p| thickness.at(p),
        options.maximum_edge_m,
        &|p| p,
        options.max_surface_triangles,
        cx,
    )?;
    let outer = ring_moments(&spec.outline);
    let sign = outer[0].signum();
    let mut m = outer.map(|v| v * sign);
    for h in &spec.holes {
        let hm = ring_moments(h);
        let s = hm[0].signum();
        for k in 0..3 {
            m[k] -= hm[k] * s;
        }
    }
    let exact = match thickness {
        ThicknessField::Uniform { thickness_m } => thickness_m * m[0],
        ThicknessField::Affine {
            at_origin_m,
            gradient,
        } => at_origin_m * m[0] + gradient[0] * m[1] + gradient[1] * m[2],
    };
    assemble_part(
        "outline-solid",
        reference_state,
        extrusion,
        ExactSurface::Polyhedral,
        BuildFacts {
            chord_tolerance_m: options.chord_tolerance_m,
            coordinate_grid_m: options.coordinate_grid_m,
            planar_minimum_angle_deg: section.receipt.minimum_angle_deg,
            max_quantization_m: max_quantization,
            exact_volume_m3: Some(exact),
        },
    )
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Generalized winding number of a closed triangle surface at `q`
/// (Van Oosterom–Strackee solid angles): 1 inside, 0 outside.
fn winding_number(vertices: &[[f64; 3]], triangles: &[[u32; 3]], q: [f64; 3]) -> f64 {
    let mut total = 0.0;
    for t in triangles {
        let [a, b, c] = t.map(|i| sub3(vertices[i as usize], q));
        let (la, lb, lc) = (
            fs_math::det::sqrt(dot3(a, a)),
            fs_math::det::sqrt(dot3(b, b)),
            fs_math::det::sqrt(dot3(c, c)),
        );
        let numerator = dot3(a, cross3(b, c));
        let denominator = la * lb * lc + dot3(a, b) * lc + dot3(b, c) * la + dot3(c, a) * lb;
        total += 2.0 * fs_math::det::atan2(numerator, denominator);
    }
    total / (4.0 * core::f64::consts::PI)
}

/// Candidate seeds tried by [`imported_solid`].
const IMPORT_SEED_CANDIDATES: usize = 64;

/// Admit an identified tessellation (for example of a STEP solid) as a
/// closed PLC. The tessellation is the surface of record: every facet is a
/// plane patch, snapping is the identity, and the part is only as accurate
/// as the tessellation, whose declared deviation is recorded as the chord
/// tolerance. Nothing is repaired: a leak, a non-manifold edge or an
/// orientation conflict refuses.
///
/// # Errors
/// Missing reference state, non-finite or out-of-range input, an open or
/// inconsistently oriented surface, non-positive enclosed volume, no
/// interior seed among the candidates, cancellation.
pub fn imported_solid(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    tessellation_deviation_m: f64,
    reference_state: Option<ReferenceState>,
    cx: &Cx<'_>,
) -> Result<SolidPart, SolidPartError> {
    let reference_state = reference_state.ok_or(SolidPartError::ReferenceStateMissing)?;
    checkpoint(cx)?;
    if !(tessellation_deviation_m.is_finite() && tessellation_deviation_m >= 0.0) {
        return Err(SolidPartError::InvalidSpec("tessellation deviation"));
    }
    if vertices.iter().flatten().any(|v| !v.is_finite()) {
        return Err(SolidPartError::InvalidSpec("non-finite vertex"));
    }
    if triangles.len() < 4 {
        return Err(SolidPartError::NotClosed {
            reason: "fewer than four triangles",
        });
    }
    let mut directed: BTreeMap<[u32; 2], u32> = BTreeMap::new();
    for (i, t) in triangles.iter().enumerate() {
        if i % 4096 == 0 {
            checkpoint(cx)?;
        }
        if t.iter().any(|&v| v as usize >= vertices.len()) {
            return Err(SolidPartError::InvalidSpec("triangle vertex out of range"));
        }
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            return Err(SolidPartError::InvalidSpec("degenerate triangle"));
        }
        for e in [[t[0], t[1]], [t[1], t[2]], [t[2], t[0]]] {
            *directed.entry(e).or_insert(0) += 1;
        }
    }
    for (&[a, b], &count) in &directed {
        if count > 1 {
            return Err(SolidPartError::NotClosed {
                reason: "an edge is used twice in one direction (non-manifold or misoriented)",
            });
        }
        match directed.get(&[b, a]) {
            None => {
                return Err(SolidPartError::NotClosed {
                    reason: "an edge has no opposite use (a leak)",
                });
            }
            Some(&back) if back != 1 => {
                return Err(SolidPartError::NotClosed {
                    reason: "an edge has more than two uses (non-manifold)",
                });
            }
            Some(_) => {}
        }
    }
    let surface_triangles = triangles.to_vec();
    if !(enclosed_volume(vertices, &surface_triangles) > 0.0) {
        return Err(SolidPartError::Orientation);
    }
    // Seed: step inward from the centroids of the largest triangles until
    // the winding number says inside.
    let mut order: Vec<(f64, usize)> = triangles
        .iter()
        .enumerate()
        .map(|(i, t)| {
            let [a, b, c] = t.map(|v| vertices[v as usize]);
            let n = cross3(sub3(b, a), sub3(c, a));
            (dot3(n, n), i)
        })
        .collect();
    order.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
    let mut seed = None;
    for &(_, i) in order.iter().take(IMPORT_SEED_CANDIDATES) {
        checkpoint(cx)?;
        let [a, b, c] = triangles[i].map(|v| vertices[v as usize]);
        let n = cross3(sub3(b, a), sub3(c, a));
        let len = fs_math::det::sqrt(dot3(n, n));
        let shortest = [sub3(b, a), sub3(c, b), sub3(a, c)]
            .iter()
            .map(|e| fs_math::det::sqrt(dot3(*e, *e)))
            .fold(f64::INFINITY, f64::min);
        let step = 0.25 * shortest / len;
        let candidate: [f64; 3] =
            core::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0 - step * n[k]);
        if winding_number(vertices, triangles, candidate) > 0.5 {
            seed = Some(candidate);
            break;
        }
    }
    let seed = seed.ok_or(SolidPartError::InvalidSpec("no interior seed found"))?;
    assemble_part(
        "imported-tessellation",
        reference_state,
        Extrusion {
            vertices: vertices.to_vec(),
            patches: vec![SurfacePatch::Plane; surface_triangles.len()],
            triangles: surface_triangles,
            layers: Vec::new(),
            seed,
        },
        ExactSurface::Polyhedral,
        BuildFacts {
            chord_tolerance_m: tessellation_deviation_m,
            coordinate_grid_m: 0.0,
            planar_minimum_angle_deg: 0.0,
            max_quantization_m: 0.0,
            exact_volume_m3: None,
        },
    )
}
