//! Source-bound rigid-body mass, balance and inertia without a second dynamics engine.
//!
//! Geometry supplies volume moments in an explicit body frame. Each region's
//! density is resolved from its material card at its full physical state point.
//! The compiler returns the existing `fs-mbd` spatial inertia, with every source
//! bundle retained. Changing material data changes mass, centre of mass, gravity
//! loads and rotational response through the UNCHANGED articulated solvers.
//!
//! This composition boundary is separate from the dependency-light `fs-mbd`
//! numerical kernel: numerical dynamics does not need a material database.
//! Regions must partition a rigid body without overlapping mass. Geometry
//! coverage, rigidity and source accuracy are caller/model obligations, not
//! consequences of a positive inertia tensor. Material uncertainty is retained,
//! not silently converted into a certified inertia bound.

#![deny(unsafe_code)]

/// Complete source-bound mechanism compilation using the existing dynamics engine.
pub mod articulated;
/// Compatible nonuniform deformation with conserved material mass.
pub mod deformation;
/// Principal-frame adaptation for the existing free-body and contact APIs.
pub mod principal;
/// Material-driven buoyancy, displacement and hydrostatic torque.
pub mod hydrostatic;
/// Source-bound creeping-flow resistance and buoyancy on native rigid motion.
pub mod creeping_flow;

/// Source-integrated free expansion with conserved reference material mass.
pub mod thermal;

mod geometry;
pub use geometry::VolumeMoments;

use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_mbd::articulated::{ArticulatedError, SpatialInertia};
use fs_qty::Density;

/// Cold-path work bound; callers may impose a smaller region budget.
pub const MAX_BODY_REGIONS: usize = 4096;

/// Failure never publishes a partially compiled body.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialInertiaError {
    /// Invalid geometry, identity, budget or unrepresentable derived arithmetic.
    Invalid(&'static str),
    /// Original material-card resolution failure.
    Material(MaterialStatePointError),
    /// Existing inertia or dynamics admission refused.
    Dynamics(ArticulatedError),
    /// Cancellation observed between bounded region operations.
    Cancelled,
}
impl core::fmt::Display for MaterialInertiaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "material inertia: {s}"),
            Self::Material(e) => write!(f, "material inertia source: {e}"),
            Self::Dynamics(e) => write!(f, "material inertia dynamics: {e}"),
            Self::Cancelled => f.write_str("material inertia: cancelled"),
        }
    }
}
impl std::error::Error for MaterialInertiaError {}
impl From<MaterialStatePointError> for MaterialInertiaError {
    fn from(e: MaterialStatePointError) -> Self { Self::Material(e) }
}
impl From<ArticulatedError> for MaterialInertiaError {
    fn from(e: ArticulatedError) -> Self { Self::Dynamics(e) }
}

/// Geometry/material assignment for one non-overlapping region, in ONE body frame.
/// Material names never choose density; the complete source query does.
#[derive(Clone, Debug)]
pub struct RegionBinding<'a> {
    /// Stable unique region key; canonical key order fixes the reduction order.
    pub key: &'a str,
    /// Volume, centroid and centroidal volume inertia, in body-frame SI units.
    pub geometry: VolumeMoments,
    /// Immutable source card for this region's density.
    pub material: &'a MaterialCard,
    /// Full physical conditions; different regions may have different temperatures.
    pub point: &'a QueryPoint,
    /// Source selection, including explicit pins when needed.
    pub selection: MaterialPropertySelection,
}

/// One region's geometry, source-resolved density and actual downstream inertia.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedRegion {
    key: String,
    geometry: VolumeMoments,
    material: ResolvedMaterialStatePoint,
    inertia: SpatialInertia,
}
impl ResolvedRegion {
    /// Original region identity.
    pub fn key(&self) -> &str { &self.key }
    /// Geometry at the authored physical state, without inferred expansion.
    pub const fn geometry(&self) -> VolumeMoments { self.geometry }
    /// Exact material state-point bundle and density usage receipt.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    /// Region inertia in the same body frame as the complete body.
    pub const fn inertia(&self) -> SpatialInertia { self.inertia }
}

/// Immutable source-backed body, suitable for the existing `Link::new` API.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBody {
    regions: Vec<ResolvedRegion>,
    inertia: SpatialInertia,
}
impl ResolvedBody {
    /// Existing canonical spatial-inertia type; no copied solver coefficients.
    pub const fn inertia(&self) -> SpatialInertia { self.inertia }
    /// Every original geometry and material bundle in canonical region-key order.
    pub fn regions(&self) -> &[ResolvedRegion] { &self.regions }
}

/// Compile material assignments into one rigid-body inertia.
///
/// Requires only sourced density, not unrelated elastic or yield data. Each
/// source's units, state domain and selection policy are enforced by the shared
/// resolver. An unsupported state or missing density is a refusal, never a
/// representative metal or a manually maintained mass override.
///
/// This creates a NEW physical model. Re-evaluating density at a new temperature
/// with fixed geometry does not conserve the old body's mass; the caller must
/// supply the corresponding geometry when expansion is intended. This API never
/// mutates a live body's momentum or invents work for instantaneous replacement.
///
/// Work is O(n log n) for sorting and O(n) fixed-size resolutions/reductions.
/// Cancellation is checked before work, per region and before publication.
pub fn resolve_body(
    bindings: &[RegionBinding<'_>],
    max_regions: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ResolvedBody, MaterialInertiaError> {
    poll(&mut cancelled)?;
    if max_regions == 0 || max_regions > MAX_BODY_REGIONS
        || bindings.is_empty() || bindings.len() > max_regions {
        return Err(MaterialInertiaError::Invalid("body needs 1..=max_regions regions within the 4096 cap"));
    }
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|b| b.key);
    for (i, binding) in ordered.iter().enumerate() {
        if binding.key.trim().is_empty() || binding.key.len() > 128
            || (i > 0 && ordered[i - 1].key == binding.key) {
            return Err(MaterialInertiaError::Invalid("region keys must be nonblank, unique and at most 128 bytes"));
        }
    }
    let requirement = ScalarPropertyRequirement::try_new(
        DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive,
    )?;
    let mut regions = Vec::with_capacity(ordered.len());
    let mut combined = None;
    for binding in ordered {
        poll(&mut cancelled)?;
        let material = resolve_material_state_point(binding.material, binding.point,
            core::slice::from_ref(&requirement), binding.selection.clone())?;
        let density = material.property(DENSITY_PROPERTY)
            .expect("the required density was resolved").value_si();
        let inertia = binding.geometry.with_density(density)?;
        combined = Some(match combined {
            None => inertia,
            Some(previous) => merge_inertias(previous, inertia)?,
        });
        regions.push(ResolvedRegion { key: binding.key.into(), geometry: binding.geometry, material, inertia });
    }
    let inertia = combined.expect("nonempty region list admitted");
    // A finite centroidal tensor alone need not have a representable spatial
    // matrix at an extreme body-frame origin. Check the actual consumer view.
    inertia.matrix()?;
    poll(&mut cancelled)?;
    Ok(ResolvedBody { regions, inertia })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialInertiaError> {
    if cancelled() { Err(MaterialInertiaError::Cancelled) } else { Ok(()) }
}
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }

// Merge about the centres, never by subtracting enormous origin moments.
// I = Ia + Ib + ma*mb/(ma+mb) * (|cb-ca|² 1 - (cb-ca)(cb-ca)^T).
fn merge_inertias(a: SpatialInertia, b: SpatialInertia) -> Result<SpatialInertia, MaterialInertiaError> {
    let mass = a.mass() + b.mass();
    if !positive(mass) { return Err(MaterialInertiaError::Invalid("combined mass is unrepresentable")); }
    let weight = b.mass() / mass;
    let reduced = a.mass().min(b.mass()) * (a.mass().max(b.mass()) / mass);
    if !positive(weight) || !positive(reduced) {
        return Err(MaterialInertiaError::Invalid("region mass ratio is unrepresentable"));
    }
    let delta = b.center_of_mass() - a.center_of_mass();
    let center = a.center_of_mass() + delta.scale(weight);
    let shift = delta.scale(reduced.sqrt());
    let [x, y, z] = [shift.x, shift.y, shift.z];
    let parallel = [y*y+z*z, -x*y, -x*z, -x*y, x*x+z*z, -y*z, -x*z, -y*z, x*x+y*y];
    let mut tensor = a.inertia_com();
    for (i, entry) in tensor.m.iter_mut().enumerate() {
        *entry += b.inertia_com().m[i] + parallel[i];
    }
    Ok(SpatialInertia::new(mass, center, tensor)?)
}
