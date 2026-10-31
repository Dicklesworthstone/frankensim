//! Attach each exterior cell's appearance to its exact accepted thermal state.

use fs_exec::Cx;
use fs_mbd_material::deformation::thermoelastic::ThermalBody;
use fs_solid::linear3::TetAssemblyBudget;

use super::super::{MaterialGeometryError, poll};
use super::{MaterialRegionMesh, deformed_body_regions};
use crate::material_state::{MaterialOpticalBinding, MaterialThermalBinding, OpticalAppearance};
use crate::tracer::{Primitive, Shape};

pub(super) mod continuation;

/// A supported physical appearance already bound to one accepted material state.
/// Surface finish and gray/diffuse assumptions stay explicit in these bindings.
#[derive(Clone, Copy, Debug)]
pub enum MaterialRegionAppearance<'a> {
    /// Existing opaque conductor observation. Dielectric bindings refuse here
    /// because this exterior partition does not contain internal optical interfaces.
    Optical(&'a MaterialOpticalBinding),
    /// Existing opaque gray/diffuse thermal emission and reflection.
    Thermal(&'a MaterialThermalBinding),
}

/// One exterior region's independently declared surface appearance.
#[derive(Clone, Copy, Debug)]
pub struct MaterialRegionBinding<'a> {
    /// Exact reference-region key. Input order never selects material ownership.
    pub key: &'a str,
    /// Appearance admitted against this cell's current constitutive state.
    pub appearance: MaterialRegionAppearance<'a>,
}

/// Refusal before any candidate primitives are exposed to a scene.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialRegionBindingError {
    /// Original geometry, budget, scale or cancellation refusal.
    Geometry(MaterialGeometryError),
    /// Missing, duplicate, unknown or non-exterior region assignment.
    Bindings(&'static str),
    /// Appearance observes a different material state or a stale temperature.
    StateMismatch {
        /// Original region requiring a fresh binding.
        region: String,
    },
}

impl core::fmt::Display for MaterialRegionBindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry(error) => error.fmt(f),
            Self::Bindings(what) => write!(f, "regional material appearance: {what}"),
            Self::StateMismatch { region } => {
                write!(
                    f,
                    "region {region} appearance differs from its accepted state"
                )
            }
        }
    }
}

impl core::error::Error for MaterialRegionBindingError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Geometry(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MaterialGeometryError> for MaterialRegionBindingError {
    fn from(error: MaterialGeometryError) -> Self {
        Self::Geometry(error)
    }
}

/// Build opaque scene primitives from one accepted body's local state and shape.
///
/// Each binding must match the exact accepted stiffness-state identity of its
/// region, including the card and complete query point. Rebind after the local
/// physical state changes. The original thermal binding owns emission, so the
/// returned primitive never receives a second artistic emission source.
///
/// Exactly one binding is required per exterior region; wholly internal cells
/// do not acquire optical-data requirements. Distinct neighboring temperatures
/// and materials remain distinct. No average temperature, material-name switch,
/// duplicated internal interface, or heat-state update is introduced.
///
/// The caller replaces scene primitives only after the complete result accepts.
/// Original body coordinates are metres in its own frame; use existing instance
/// transforms for placement. Fixed topology and the appearance models' existing
/// physical limits apply. Transparent regional media require internal optical
/// interfaces and medium tracking that this exterior-only observer does not supply;
/// dielectric bindings therefore refuse rather than dropping those interfaces.
///
/// # Errors
/// Refuses geometry/budget/cancellation failures, incomplete or ambiguous region
/// assignments, and stale or foreign accepted material states.
pub fn thermal_body_primitives(
    body: &ThermalBody<'_>,
    bindings: &[MaterialRegionBinding<'_>],
    budget: TetAssemblyBudget,
    maximum_face_entries: usize,
    cx: &Cx<'_>,
) -> Result<Vec<Primitive>, MaterialRegionBindingError> {
    let regions = deformed_body_regions(body.body(), budget, maximum_face_entries, cx)?;
    if bindings.len() != regions.len() {
        return Err(MaterialRegionBindingError::Bindings(
            "exactly one appearance per exterior region is required",
        ));
    }
    let mut ordered = bindings.to_vec();
    ordered.sort_by_key(|binding| binding.key);
    bind_ordered_regions(body, regions, &ordered, cx)
}

// Both fresh and continuing observations consume the same extracted exterior.
// In particular a continuation must not rebuild every region's BVH twice.
fn bind_ordered_regions(
    body: &ThermalBody<'_>,
    regions: Vec<MaterialRegionMesh>,
    bindings: &[MaterialRegionBinding<'_>],
    cx: &Cx<'_>,
) -> Result<Vec<Primitive>, MaterialRegionBindingError> {
    if bindings.len() != regions.len() {
        return Err(MaterialRegionBindingError::Bindings(
            "exactly one appearance per exterior region is required",
        ));
    }
    let mut primitives = Vec::with_capacity(regions.len());
    for (region, binding) in regions.into_iter().zip(bindings) {
        poll(cx)?;
        let key = body.body().reference().regions()[region.region_index].key();
        if binding.key != key {
            return Err(MaterialRegionBindingError::Bindings(
                "appearance keys must match the exact exterior region inventory",
            ));
        }
        let (identity, material) = match binding.appearance {
            MaterialRegionAppearance::Optical(binding) => {
                if matches!(binding.appearance(), OpticalAppearance::Dielectric { .. }) {
                    return Err(MaterialRegionBindingError::Bindings(
                        "regional exterior appearance requires an opaque material",
                    ));
                }
                (binding.accepted_state_identity(), binding.material())
            }
            MaterialRegionAppearance::Thermal(binding) => {
                (binding.accepted_state_identity(), binding.material())
            }
        };
        if identity != body.sources()[region.region_index].stiffness().identity() {
            return Err(MaterialRegionBindingError::StateMismatch {
                region: key.to_owned(),
            });
        }
        primitives.push(Primitive {
            shape: Shape::Mesh(region.mesh),
            material,
            emission: None,
        });
    }
    poll(cx)?;
    Ok(primitives)
}
