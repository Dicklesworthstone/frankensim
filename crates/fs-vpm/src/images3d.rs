//! Physical-wake flat-plane image integration (bead
//! wf-root-guzez.5.19, E4.7b). Extends the E4.4a
//! FlatPlaneVortexImageExact doctrine (fs-wing::images) to the fs-vpm
//! hybrid wake: every wake element — near filament row, coarsened
//! mid-wake row, far multipole cell — gains an exact mirror across a
//! certified flat plane.
//!
//! Doctrine carried over from V-06a:
//!   - images are BOUNDARY DEVICES, not physical vorticity: they are
//!     built at evaluation time and NEVER enter the physical ledgers
//!     (per-station impulse, Kelvin cell closure, wake digests) — the
//!     battery proves the exclusion by executing the forbidden
//!     materialization and watching the ledger corrupt;
//!   - the image of a vortex segment mirrors its endpoints and
//!     negates its circulation (axial-vector reflection,
//!     ω' = det(R)·R·ω with det(R) = −1);
//!   - image identities are STABLE: a domain-separated hash of the
//!     source identity and the plane, reproducible across runs.
//!
//! Conversion symmetry: mirroring COMMUTES with the near→mid→far
//! conversion chain (coarsening and multipole aggregation), because
//! image cells are built by the SAME `cell_from_segs` arithmetic as
//! real cells. The battery asserts the commutation per station and
//! per cell, never as a totals-only sum.

use crate::farfield3d::{FarCell, FarField, cell_from_segs};
use crate::filament3d::{FilamentWake, Refusal, ShedRow};
use fs_blake3::hash_domain;

/// The certified image plane z = z0 (the wake must sit strictly
/// above; below-plane wake states route to contact, never to images).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImagePlane {
    /// Plane height [m].
    pub z0_m: f64,
}

impl ImagePlane {
    /// Admit the plane against a wake: every node strictly above.
    ///
    /// # Errors
    /// `image-plane-invalid` (non-finite plane, or any wake/line node
    /// at or below it — clearance > 0 admits, clearance = 0 refuses:
    /// the cap AND cap+1 pair is exact here).
    pub fn admit(&self, wake: &FilamentWake) -> Result<(), Refusal> {
        let mut min_z = f64::INFINITY;
        for n in wake
            .rows
            .iter()
            .flat_map(|r| r.nodes.iter())
            .chain(wake.line_nodes.iter())
        {
            min_z = min_z.min(n[2]);
        }
        if self.z0_m.is_finite() && min_z > self.z0_m {
            Ok(())
        } else {
            Err(Refusal {
                code: "image-plane-invalid",
                message: format!("plane z0 {} with min wake z {min_z}", self.z0_m),
                ranked_repairs: vec![
                    "a wake at or below the plane is a contact state, not an image case".into(),
                ],
            })
        }
    }

    /// Mirror a point across the plane.
    #[must_use]
    pub fn mirror(&self, p: [f64; 3]) -> [f64; 3] {
        [p[0], p[1], 2.0 * self.z0_m - p[2]]
    }
}

/// Stable image identity: domain-separated hash of the SOURCE
/// identity and the plane. Reproducible across runs; distinct planes
/// give distinct identities.
#[must_use]
pub fn image_identity(source_digest: &str, plane: ImagePlane) -> String {
    let mut b = source_digest.as_bytes().to_vec();
    b.extend_from_slice(&plane.z0_m.to_bits().to_le_bytes());
    hash_domain("org.frankensim.fs-vpm.image-identity.v1", &b).to_hex()
}

/// The image of the NEAR wake: nodes mirrored, circulations negated
/// (a plain FilamentWake so evaluation reuses the exact segment path
/// — but it is a boundary device the caller must never shed into or
/// feed to the ledgers; `image_system` is the lawful constructor).
#[must_use]
fn mirror_wake(wake: &FilamentWake, plane: ImagePlane) -> FilamentWake {
    let mut img = wake.clone();
    for n in &mut img.line_nodes {
        *n = plane.mirror(*n);
    }
    for row in &mut img.rows {
        for n in &mut row.nodes {
            *n = plane.mirror(*n);
        }
        for g in &mut row.gamma {
            *g = -*g;
        }
    }
    img
}

/// The image of one far cell, built from its mirrored exact segments
/// through the SAME `cell_from_segs` arithmetic as real aggregation
/// (this is what makes conversion symmetry hold).
#[must_use]
pub fn mirror_cell(cell: &FarCell, plane: ImagePlane) -> FarCell {
    let segs: Vec<([f64; 3], [f64; 3], f64)> = cell
        .segs
        .iter()
        .map(|(a, b, g)| (plane.mirror(*a), plane.mirror(*b), -*g))
        .collect();
    cell_from_segs(segs, cell.core2_m2)
}

/// The complete image system for a hybrid wake (near + far).
#[derive(Clone, Debug)]
pub struct ImageSystem {
    /// Plane it was certified against.
    pub plane: ImagePlane,
    /// Stable identity (source digests + plane).
    pub identity: String,
    /// Mirrored near wake (boundary device).
    near_image: FilamentWake,
    /// Mirrored far cells (boundary devices).
    far_image: Vec<FarCell>,
}

impl ImageSystem {
    /// Build the image system. The ONLY lawful constructor: it admits
    /// the plane first and never touches the source wake.
    ///
    /// # Errors
    /// Plane refusals (`image-plane-invalid`).
    pub fn build(near: &FilamentWake, far: &FarField, plane: ImagePlane) -> Result<Self, Refusal> {
        plane.admit(near)?;
        let mut src = near.digest().into_bytes();
        src.extend_from_slice(far.digest().as_bytes());
        let identity = image_identity(core::str::from_utf8(&src).unwrap_or(""), plane);
        Ok(ImageSystem {
            plane,
            identity,
            near_image: mirror_wake(near, plane),
            far_image: far.cells.iter().map(|c| mirror_cell(c, plane)).collect(),
        })
    }

    /// Induced velocity of the IMAGES alone at `p`.
    #[must_use]
    pub fn eval(&self, p: [f64; 3]) -> [f64; 3] {
        let mut v = self.near_image.induced_velocity(p);
        for c in &self.far_image {
            let cv = c.eval(p);
            v[0] += cv[0];
            v[1] += cv[1];
            v[2] += cv[2];
        }
        v
    }

    /// Access the mirrored far cells (battery oracle surface).
    #[must_use]
    pub fn far_image_cells(&self) -> &[FarCell] {
        &self.far_image
    }

    /// Access the mirrored near wake (battery oracle surface).
    #[must_use]
    pub fn near_image_wake(&self) -> &FilamentWake {
        &self.near_image
    }
}

/// Grounded hybrid induced velocity: exact near + far multipoles +
/// their images. The images enter HERE and only here — the physical
/// wake, its ledgers, and its digests are untouched.
#[must_use]
pub fn grounded_hybrid_velocity(
    near: &FilamentWake,
    far: &FarField,
    images: &ImageSystem,
    p: [f64; 3],
) -> [f64; 3] {
    let nv = near.induced_velocity(p);
    let fv = far.eval(p);
    let iv = images.eval(p);
    [
        nv[0] + fv[0] + iv[0],
        nv[1] + fv[1] + iv[1],
        nv[2] + fv[2] + iv[2],
    ]
}

/// The FORBIDDEN materialization (battery falsifier ONLY): append the
/// image rows to the physical wake as if they were real vorticity.
/// This is exactly what the ledger-exclusion doctrine forbids — the
/// battery executes it and proves the per-station impulse ledger
/// corrupts, which is the evidence that the lawful path's exclusion
/// is real and the oracle is sensitive.
#[doc(hidden)]
pub fn unlawful_materialize_images(wake: &mut FilamentWake, plane: ImagePlane) {
    let img = mirror_wake(wake, plane);
    let rows: Vec<ShedRow> = img.rows;
    wake.rows.extend(rows);
}
