//! Source-bound free thermal expansion of an already resolved mass inventory.
//!
//! The shared material owner integrates alpha(T). This adapter applies its
//! engineering free strain to reference geometry and transports the SAME mass
//! into the existing spatial-inertia owner. It never queries a new density on
//! an old volume, silently replaces matter, or updates a live trajectory.
//!
//! This is prescribed isotropic free expansion about explicit fixed points.
//! It is not a constrained thermoelastic solve: different free expansions may
//! be incompatible in a bonded assembly. Stress, heat flow, contact, rigidity,
//! overlap, momentum transfer and work during configuration changes remain the
//! relevant solvers' responsibilities. The shared strain law and source domain
//! bound applicability; no finite-strain constitutive law is introduced here.

/// Material-bound thermal solver temperatures drive conserved geometry and inertia.
pub mod temperature;
/// Shared heat transport and source-bound mechanical configuration.
pub mod heating;

mod articulated;
pub use articulated::{
    ExpandedArticulated, JointMaterialAttachment, LinkThermalExpansion, expand_articulated,
};

use fs_ga::Vec3;
use fs_material::state_point::{IntegratedIsotropicThermalExpansion, ResolvedMaterialStatePoint};
use fs_mbd::articulated::SpatialInertia;

use crate::{MaterialInertiaError, ResolvedBody, VolumeMoments, merge_inertias, poll, positive};

/// One complete, source-bound temperature path for an existing material region.
#[derive(Clone, Copy, Debug)]
pub struct ThermalExpansionBinding<'a> {
    /// Existing unique region key. Every reference region must appear once.
    pub key: &'a str,
    /// Receipt produced by `fs_material::state_point::integrate_isotropic_thermal_expansion`.
    /// Its reference card, material and full state must match the mass inventory.
    pub expansion: &'a IntegratedIsotropicThermalExpansion,
    /// Fixed point of the free expansion, in the original body frame [m].
    /// Choose explicitly; the adapter does not invent mounting constraints.
    pub fixed_point_m: Vec3,
}

/// Deformed geometry and unchanged mass with the actual thermal-source receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedRegion {
    reference_index: usize,
    expansion: IntegratedIsotropicThermalExpansion,
    fixed_point_m: Vec3,
    geometry: VolumeMoments,
    inertia: SpatialInertia,
}
impl ExpandedRegion {
    /// Index in the retained reference body's canonical region order.
    pub const fn reference_index(&self) -> usize { self.reference_index }
    /// Original path/claim/uncertainty receipts, not a synthesized density claim.
    pub const fn expansion(&self) -> &IntegratedIsotropicThermalExpansion { &self.expansion }
    /// Explicit prescribed fixed point in the reference body frame [m].
    pub const fn fixed_point_m(&self) -> Vec3 { self.fixed_point_m }
    /// Expanded volume and complete centroidal volume tensor in the body frame.
    pub const fn geometry(&self) -> VolumeMoments { self.geometry }
    /// Existing canonical mass, centre of mass and rotational inertia.
    pub const fn inertia(&self) -> SpatialInertia { self.inertia }
    /// Density derived by conserving the reference mass on the expanded volume.
    /// This is NOT a newly queried density measurement or a confidence bound.
    pub fn density_kg_m3(&self) -> f64 { self.inertia.mass() / self.geometry.volume_m3() }
    /// Carry an authored reference point with this region's actual free expansion.
    /// The point and result are in the unchanged body frame. Geometric membership
    /// in this region is an explicit caller assignment, not inferred from moments.
    pub fn transport_reference_point(&self, point_m: Vec3) -> Result<Vec3, MaterialInertiaError> {
        transport_point(
            point_m,
            self.fixed_point_m,
            1.0 + self.expansion.free_linear_strain(),
        )
    }
}

/// Immutable configuration of a retained reference mass inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedBody<'a> {
    reference: &'a ResolvedBody,
    regions: Vec<ExpandedRegion>,
    inertia: SpatialInertia,
}
impl<'a> ExpandedBody<'a> {
    /// Original geometry, source density and exact material query receipts.
    pub const fn reference(&self) -> &'a ResolvedBody { self.reference }
    /// Configured regions, in the same canonical order as the reference body.
    pub fn regions(&self) -> &[ExpandedRegion] { &self.regions }
    /// Direct input to the unchanged articulated forward/inverse dynamics APIs.
    pub const fn inertia(&self) -> SpatialInertia { self.inertia }
}

/// Apply source-resolved free thermal strains without creating or destroying mass.
///
/// Each path must begin at the exact card and physical state that supplied its
/// reference region's density. The shared material receipt enforces a single
/// expansion claim along the path. Geometry uses `s = 1 + integral(alpha dT)`,
/// the existing engineering-free-strain convention, not a second alpha(T) law.
/// Volume scales as s^3, central mass inertia as s^2, and each region's mass is
/// copied unchanged. Tensor rotation and parallel-axis composition stay with
/// the existing geometry and dynamics owners.
///
/// Bindings may be listed in any order but must cover all reference regions
/// exactly once. Work is bounded by that body's admitted region count. Source,
/// geometry, density and downstream matrix failures never publish a partial
/// result. Cancellation is checked before work, per region and before return.
/// Call again against the ORIGINAL reference body for a new temperature path;
/// repeated calls do not compound incremental thermal strains.
pub fn expand_body<'a>(
    reference: &'a ResolvedBody,
    bindings: &[ThermalExpansionBinding<'_>],
    mut cancelled: impl FnMut() -> bool,
) -> Result<ExpandedBody<'a>, MaterialInertiaError> {
    poll(&mut cancelled)?;
    if bindings.len() != reference.regions().len() {
        return Err(MaterialInertiaError::Invalid("thermal bindings must cover every reference region exactly once"));
    }
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.key);
    let mut regions = Vec::with_capacity(ordered.len());
    let mut combined = None;
    for (index, (original, binding)) in reference.regions().iter().zip(ordered).enumerate() {
        poll(&mut cancelled)?;
        if binding.key != original.key() {
            return Err(MaterialInertiaError::Invalid("thermal region keys are missing, duplicated or unknown"));
        }
        if !same_reference(original.material(), binding.expansion.reference().resolved()) {
            return Err(MaterialInertiaError::Invalid("thermal path does not start at the exact reference material state"));
        }
        let fixed = binding.fixed_point_m;
        if ![fixed.x, fixed.y, fixed.z].into_iter().all(f64::is_finite) {
            return Err(MaterialInertiaError::Invalid("thermal fixed point must be finite"));
        }
        let stretch = 1.0 + binding.expansion.free_linear_strain();
        if !positive(stretch) {
            return Err(MaterialInertiaError::Invalid("free strain must give a finite positive stretch"));
        }
        let area_scale = stretch * stretch;
        let volume_scale = area_scale * stretch;
        let moment_scale = volume_scale * area_scale;
        if ![area_scale, volume_scale, moment_scale].into_iter().all(positive) {
            return Err(MaterialInertiaError::Invalid("expanded geometry scale is unrepresentable"));
        }
        let old_geometry = original.geometry();
        // Avoid subtracting nearly equal large global coordinates when s~1.
        let center = transport_point(old_geometry.centroid_m(), fixed, stretch)?;
        let mut volume_tensor = old_geometry.inertia_m5();
        for entry in &mut volume_tensor.m { *entry *= moment_scale; }
        let geometry = VolumeMoments::new(old_geometry.volume_m3() * volume_scale, center, volume_tensor)?;
        let mut mass_tensor = original.inertia().inertia_com();
        for entry in &mut mass_tensor.m { *entry *= area_scale; }
        // Do not recompute rho(T)*old_volume or round-trip old_mass/new_volume.
        // The accepted mass inventory is copied bit-for-bit into the new state.
        let inertia = SpatialInertia::new(original.inertia().mass(), center, mass_tensor)?;
        if !positive(inertia.mass() / geometry.volume_m3()) {
            return Err(MaterialInertiaError::Invalid("mass-conserving current density is unrepresentable"));
        }
        combined = Some(match combined {
            None => inertia,
            Some(previous) => merge_inertias(previous, inertia)?,
        });
        regions.push(ExpandedRegion {
            reference_index: index, expansion: (*binding.expansion).clone(),
            fixed_point_m: fixed, geometry, inertia,
        });
    }
    let inertia = combined.expect("a resolved reference body has at least one region");
    inertia.matrix()?;
    poll(&mut cancelled)?;
    Ok(ExpandedBody { reference, regions, inertia })
}

fn transport_point(point: Vec3, fixed: Vec3, stretch: f64) -> Result<Vec3, MaterialInertiaError> {
    if ![point.x, point.y, point.z].into_iter().all(f64::is_finite) {
        return Err(MaterialInertiaError::Invalid(
            "reference material point must be finite",
        ));
    }
    // Preserve exact reference geometry at zero strain and avoid subtracting
    // large fixed-point coordinates merely to reconstruct the original point.
    let result = if stretch == 1.0 {
        point
    } else {
        point + (point - fixed).scale(stretch - 1.0)
    };
    if ![result.x, result.y, result.z]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(MaterialInertiaError::Invalid(
            "expanded material point is unrepresentable",
        ));
    }
    Ok(result)
}

fn same_reference(mass: &ResolvedMaterialStatePoint, expansion: &ResolvedMaterialStatePoint) -> bool {
    mass.card_identity() == expansion.card_identity()
        && mass.material() == expansion.material()
        && mass.axis_quantities() == expansion.axis_quantities()
        && mass.query_point().len() == expansion.query_point().len()
        && mass.query_point().iter().zip(expansion.query_point()).all(|((a, x), (b, y))| {
            a == b && x.to_bits() == y.to_bits()
        })
}
