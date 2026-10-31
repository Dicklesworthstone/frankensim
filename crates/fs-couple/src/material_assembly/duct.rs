//! Source-bound duct walls in the same compiler as other acoustic materials.
//!
//! Property requirements and constitutive reduction come from the wall owner
//! already used by graph and DG wind performances. This path adds no material
//! coefficients, loss fit, gas-state inference or acoustic solver.

pub use crate::bernoulli_aperture::wall::material::source::{
    WallConstitutiveSource, WallMaterialLossSource,
};
use crate::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall};
use fs_scenario::AcousticAssembly;

use super::{MaterialAssemblyError, MaterialSource, missing, physical, resolve};

/// One original wall material, independently prescribed state and local geometry.
/// The same membrane covers all admitted axial segments and chimneys. Regional
/// wall assignments, thermal evolution and axial shell coupling are not inferred.
#[derive(Clone, Debug)]
pub struct DuctWallMaterialBinding<'a> {
    /// Original immutable card, material query and explicit claim selection.
    pub source: MaterialSource<'a>,
    /// Explicit isotropic or principal-aligned constitutive source selection.
    pub constitutive: WallConstitutiveSource,
    /// Radial thickness at this material state [m].
    pub thickness_m: f64,
    /// Physical axial support assumption of the local membrane model.
    pub axial: AxialConstraint,
    /// Positive material applicability interval [Hz], not an audio filter.
    pub band_hz: [f64; 2],
}

pub(super) fn compile(
    template: &AcousticAssembly,
    binding: &DuctWallMaterialBinding<'_>,
) -> Result<LocalMaterialWall, MaterialAssemblyError> {
    let duct = template.duct.as_ref().ok_or_else(|| {
        missing("duct wall", "wall material binding requires a duct description")
    })?;
    if duct.wall.is_some() {
        return Err(missing(
            "duct wall",
            "source-bound duct wall cannot also select authored wall coefficients",
        ));
    }
    let requirements = binding.constitutive.requirements()
        .map_err(|error| physical("duct wall", error))?;
    let material = resolve(&binding.source, &requirements, "duct wall")?;
    let wall = binding.constitutive.resolve(
        material, binding.thickness_m, binding.axial, binding.band_hz,
    ).map_err(|error| physical("duct wall", error))?;
    // The original local wall owner admits both ends of every tapered segment.
    // The acoustic realizer retains topology, tone-hole and numerical admission.
    for segment in &duct.segments {
        for radius in [segment.radius_m, segment.outlet_radius_m] {
            wall.base_patch(radius, segment.length_m)
                .map_err(|error| physical("duct wall", error))?;
        }
    }
    Ok(wall)
}
