//! Original full-spectrum film stacks attached to exact native boundaries.

use super::*;

/// An original ordered spectral stack on one existing dielectric boundary.
/// The object's original `Interaction::Dielectric` retains its smooth/GGX
/// roughness and explicit outside/inside media. Its surface must be uncoated:
/// overlapping compact or single-film treatments are not silently replaced.
/// Layers run outside to inside; the actual incident/target bulk indices are
/// evaluated by the native path at every wavelength and physical contact.
#[derive(Clone, Copy, Debug)]
pub struct DielectricStackBinding<'a> {
    pub object_id: u64,
    pub stack: DielectricCoatingStack<'a>,
}

impl<'a> PolarizedScene<'a> {
    /// Replace optional borrowed spectral-stack bindings. Empty input disables
    /// this treatment; geometry, bulk absorption and original sampling streams
    /// remain owned by the existing scene. Preparation stores only bounded
    /// slice descriptors; paths never allocate or fit/resample optical laws.
    pub fn with_dielectric_coating_stacks(
        mut self,
        cx: &Cx<'_>,
        bindings: &[DielectricStackBinding<'a>],
    ) -> Result<Self, SceneError> {
        self.coating_stacks = admit_coating_stacks(cx, &self.bindings, bindings)?;
        Ok(self)
    }
}

pub(super) fn admit_coating_stacks<'a>(
    cx: &Cx<'_>,
    surfaces: &[SurfaceBinding],
    bindings: &[DielectricStackBinding<'a>],
) -> Result<Vec<Option<DielectricCoatingStack<'a>>>, SceneError> {
    cx.checkpoint()?;
    if bindings.len() > surfaces.len() {
        return Err(SceneError::InvalidBindings);
    }
    let mut slots = Vec::new();
    if !bindings.is_empty() {
        slots
            .try_reserve_exact(surfaces.len())
            .map_err(|_| SceneError::AllocationRefused)?;
        slots.resize(surfaces.len(), None);
        for binding in bindings {
            cx.checkpoint()?;
            let index = surfaces
                .binary_search_by_key(&binding.object_id, |surface| surface.object_id)
                .map_err(|_| SceneError::InvalidBindings)?;
            if slots[index].is_some()
                || !matches!(surfaces[index].interaction,
                    Interaction::Dielectric { surface, .. } if surface.coating().is_none())
            {
                return Err(SceneError::InvalidBindings);
            }
            slots[index] = Some(binding.stack);
        }
    }
    cx.checkpoint()?;
    Ok(slots)
}
