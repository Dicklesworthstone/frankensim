//! Source-bound duct walls in the same compiler as other acoustic materials.
//!
//! Property requirements and constitutive reduction come from the wall owner
//! already used by graph and DG wind performances. This path adds no material
//! coefficients, loss fit, gas-state inference or acoustic solver.

pub use crate::bernoulli_aperture::wall::material::source::{
    WallConstitutiveSource, WallMaterialLossSource,
};
use crate::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall};
pub use crate::driving_point::regional::material_wall::MaterialWallLayout;
use fs_scenario::AcousticAssembly;

use super::{MaterialAssemblyError, MaterialSource, missing, physical, resolve};

/// One original wall material, independently prescribed state and local geometry.
/// Thermal evolution and axial shell coupling are not inferred.
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

/// Exact source assignment on the authored geometry, independent of gas regions.
#[derive(Clone, Debug)]
pub enum DuctWallMaterialBindings<'a> {
    /// One original wall material for every admitted segment.
    Uniform(DuctWallMaterialBinding<'a>),
    /// Complete ordered axial partition; no missing-source fallback is allowed.
    Segments(Vec<DuctWallMaterialBinding<'a>>),
}

pub(super) fn compile(
    template: &AcousticAssembly,
    binding: &DuctWallMaterialBindings<'_>,
) -> Result<MaterialWallLayout, MaterialAssemblyError> {
    let duct = template.duct.as_ref().ok_or_else(|| {
        missing("duct wall", "wall material binding requires a duct description")
    })?;
    if duct.wall.is_some() {
        return Err(missing(
            "duct wall",
            "source-bound duct wall cannot also select authored wall coefficients",
        ));
    }
    match binding {
        DuctWallMaterialBindings::Uniform(binding) => {
            let wall = compile_source(binding, "duct wall")?;
            for segment in &duct.segments {
                admit_segment(&wall, segment, "duct wall")?;
            }
            Ok(MaterialWallLayout::Uniform(wall))
        }
        DuctWallMaterialBindings::Segments(bindings) => {
            if bindings.is_empty() || bindings.len() != duct.segments.len()
                || bindings.len() > crate::driving_point::regional::MAX_REGIONAL_SEGMENTS {
                return Err(missing("duct wall", "wall partition needs 1..=256 sources and exactly one per axial segment"));
            }
            if !duct.tone_holes.is_empty() {
                return Err(missing("duct wall", "segment wall sources do not provide cross-region tone-hole coupling"));
            }
            let mut walls = Vec::with_capacity(bindings.len());
            for (index, (binding, segment)) in bindings.iter().zip(&duct.segments).enumerate() {
                let component = format!("duct wall segment {index}");
                let wall = compile_source(binding, &component)?;
                admit_segment(&wall, segment, &component)?;
                walls.push(wall);
            }
            Ok(MaterialWallLayout::Segments(walls))
        }
    }
}

fn compile_source(binding: &DuctWallMaterialBinding<'_>, component: &str)
    -> Result<LocalMaterialWall, MaterialAssemblyError> {
    let requirements = binding.constitutive.requirements()
        .map_err(|error| physical(component, error))?;
    let material = resolve(&binding.source, &requirements, component)?;
    let wall = binding.constitutive.resolve(
        material, binding.thickness_m, binding.axial, binding.band_hz,
    ).map_err(|error| physical(component, error))?;
    Ok(wall)
}

fn admit_segment(wall: &LocalMaterialWall, segment: &fs_scenario::CylinderSegment, component: &str)
    -> Result<(), MaterialAssemblyError> {
    // The original local wall owner admits both ends of every tapered segment.
    // The acoustic realizer retains topology, tone-hole and numerical admission.
    for radius in [segment.radius_m, segment.outlet_radius_m] {
        wall.base_patch(radius, segment.length_m)
            .map_err(|error| physical(component, error))?;
    }
    Ok(())
}
