//! Conforming fluid interfaces from original material/EOS acoustic states.
//!
//! Each tetrahedron refers to an explicitly supplied source slot. A used source
//! is queried once during preparation and its coherent density/sound-speed pair
//! is frozen in the original heterogeneous pressure operator. Shared pressure
//! DOFs and inverse-density stiffness impose pressure and normal-velocity
//! continuity at interfaces without a separate reflection/transmission model.
//! Gas, liquid, dense-fluid and sourced material-card owners all use the same
//! `FluidAcoustics` interface; their phase and relaxation policies remain theirs.

use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

/// Refusals preserve the original source error and its explicit region slot.
#[derive(Debug)]
pub enum FluidRegionsInteriorError {
    /// Every original tetrahedron must have exactly one region assignment.
    ElementCount { expected: usize, supplied: usize },
    /// A tetrahedron names a source slot which the caller did not supply.
    UnknownRegion {
        tetrahedron: usize,
        region: usize,
        source_count: usize,
    },
    /// An original material/EOS refused its coherent acoustic state.
    Material {
        region: usize,
        source: MaterialError,
    },
    /// The original finite-element owner refused geometry or assembly.
    Interior(AcousticFeError),
    /// Preparation stopped without publishing an operator or partial binding.
    Cancelled,
}

impl core::fmt::Display for FluidRegionsInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ElementCount { expected, supplied } => write!(
                f,
                "fluid regions require {expected} tetrahedron assignments, received {supplied}"
            ),
            Self::UnknownRegion {
                tetrahedron,
                region,
                source_count,
            } => write!(
                f,
                "acoustic tetrahedron {tetrahedron} names fluid region {region}, but only {source_count} sources were supplied"
            ),
            Self::Material { region, source } => write!(f, "fluid region {region}: {source}"),
            Self::Interior(error) => error.fmt(f),
            Self::Cancelled => write!(f, "fluid-region acoustic preparation cancelled"),
        }
    }
}

impl std::error::Error for FluidRegionsInteriorError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material { source, .. } => Some(source),
            Self::Interior(error) => Some(error),
            _ => None,
        }
    }
}

impl From<AcousticFeError> for FluidRegionsInteriorError {
    fn from(error: AcousticFeError) -> Self {
        Self::Interior(error)
    }
}

/// Original fluid sources and their frozen, conforming three-dimensional volume.
///
/// Region boundaries are original tetrahedron faces. The geometry's homogeneous
/// medium supplies no volume coefficients. This binding adds no constitutive law,
/// fictitious reference density, loss model, or phase conversion. It describes a
/// stationary, linear, isotropic fluid with no mean flow. Changing a source after
/// preparation does not update the assembled coefficients: prepare a new binding
/// for a new equilibrium. The borrowed sources retain their provenance/receipts;
/// the snapshots below expose the exact numerical state used in assembly.
pub struct FluidRegionsInterior<'f, 'c, 's> {
    sources: Vec<&'s dyn FluidAcoustics>,
    element_regions: Vec<usize>,
    states: Vec<Option<FluidAcousticState>>,
    operator: HeterogeneousHelmholtzFe<'f, 'c>,
}

impl<'f, 'c, 's> FluidRegionsInterior<'f, 'c, 's> {
    /// Assemble using one explicit source index per original tetrahedron.
    ///
    /// The complete map is checked before any material query. Only used source
    /// slots are sampled, once each, through their original `FluidAcoustics`
    /// owner. Repeated slots share a snapshot across all elements and quadrature
    /// points. The original FE owner admits quadrature and assembles both physical
    /// matrices together; its modes, harmonic drives and transient preparation
    /// reuse those matrices without resampling sources. No interface may cut the
    /// interior of an element in this region representation.
    pub fn new(
        geometry: &'f HelmholtzFe<'c>,
        sources: &[&'s dyn FluidAcoustics],
        element_regions: &[usize],
        quadrature_order: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FluidRegionsInteriorError> {
        if cancelled() {
            return Err(FluidRegionsInteriorError::Cancelled);
        }
        if element_regions.len() != geometry.tetrahedra().len() {
            return Err(FluidRegionsInteriorError::ElementCount {
                expected: geometry.tetrahedra().len(),
                supplied: element_regions.len(),
            });
        }
        for (tetrahedron, &region) in element_regions.iter().enumerate() {
            if cancelled() {
                return Err(FluidRegionsInteriorError::Cancelled);
            }
            if region >= sources.len() {
                return Err(FluidRegionsInteriorError::UnknownRegion {
                    tetrahedron,
                    region,
                    source_count: sources.len(),
                });
            }
        }
        let mut states = vec![None; sources.len()];
        let mut was_cancelled = false;
        let assembled = HeterogeneousHelmholtzFe::new(
            geometry,
            quadrature_order,
            &mut |element, _| {
                let region = element_regions[element];
                let state = match states[region] {
                    Some(state) => state,
                    None => {
                        let state = sources[region].acoustic_state().map_err(|source| {
                            FluidRegionsInteriorError::Material { region, source }
                        })?;
                        states[region] = Some(state);
                        state
                    }
                };
                Ok::<_, FluidRegionsInteriorError>(Medium {
                    density_kg_m3: state.density_kg_m3(),
                    sound_speed_m_s: state.sound_speed_m_s(),
                })
            },
            &mut || {
                was_cancelled |= cancelled();
                was_cancelled
            },
        );
        if was_cancelled {
            return Err(FluidRegionsInteriorError::Cancelled);
        }
        Ok(Self {
            sources: sources.to_vec(),
            element_regions: element_regions.to_vec(),
            states,
            operator: assembled?,
        })
    }

    /// Original borrowed source slots, including explicitly supplied unused ones.
    #[must_use]
    pub fn sources(&self) -> &[&'s dyn FluidAcoustics] {
        &self.sources
    }

    /// Region indices in the original geometry's tetrahedron order.
    #[must_use]
    pub fn element_regions(&self) -> &[usize] {
        &self.element_regions
    }

    /// Exact frozen acoustic snapshot, or `None` for an unused/unknown region.
    /// This never reevaluates an original source.
    #[must_use]
    pub fn region_state(&self, region: usize) -> Option<FluidAcousticState> {
        self.states.get(region).copied().flatten()
    }

    /// Original sparse physical operator and its native solver consumers.
    #[must_use]
    pub fn operator(&self) -> &HeterogeneousHelmholtzFe<'f, 'c> {
        &self.operator
    }
}
