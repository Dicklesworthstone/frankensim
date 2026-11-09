//! Original material membranes coupled to resolved lossy acoustic volumes.
//!
//! Real mass, stiffness, resistance and every positive Maxwell arm enter the
//! existing complete pressure/wall harmonic pencil. Its exp(+i omega t)
//! convention follows directly from those real equations; the source owner's
//! exp(-i omega t) impedance is not copied into this operator. Surface area is
//! supplied exactly once by the original FE trace mass.

use super::{CylindricalMaterialBoundary, MaterialMechanicalError};
use crate::cavity_modal::interior::heterogeneous::porous::{
    AcousticVolumeRegion, PorousRegionsInterior, PorousRegionsInteriorError,
};
use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::MechanicalBoundary;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidControls, HeterogeneousError, PreparedEquivalentFluid,
};

/// Original source or geometry refusals remain distinct from volume errors.
#[derive(Debug)]
pub enum MaterialMechanicalHarmonicError {
    Wall(MaterialMechanicalError),
    Volume(PorousRegionsInteriorError),
    /// One original wall source does not admit the requested frequency.
    Frequency {
        boundary: usize,
        frequency_hz: f64,
        source_band_hz: [f64; 2],
    },
}

impl core::fmt::Display for MaterialMechanicalHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Wall(error) => error.fmt(f),
            Self::Volume(error) => error.fmt(f),
            Self::Frequency {
                boundary,
                frequency_hz,
                source_band_hz,
            } => write!(
                f,
                "material wall {boundary}: {frequency_hz} Hz is outside original band {source_band_hz:?} Hz",
            ),
        }
    }
}
impl std::error::Error for MaterialMechanicalHarmonicError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Wall(error) => Some(error),
            Self::Volume(error) => Some(error),
            Self::Frequency { .. } => None,
        }
    }
}
impl From<MaterialMechanicalError> for MaterialMechanicalHarmonicError {
    fn from(value: MaterialMechanicalError) -> Self {
        Self::Wall(value)
    }
}
impl From<PorousRegionsInteriorError> for MaterialMechanicalHarmonicError {
    fn from(value: PorousRegionsInteriorError) -> Self {
        Self::Volume(value)
    }
}

/// Original wall and volume material sources retained with one coupled solve.
///
/// The wall remains the admitted thin-cylinder local membrane model: there is
/// no shell bending, lateral wall elasticity, prestress or geometry evolution.
/// Porous regions retain their rigid-frame equivalent-fluid applicability.
/// All coefficients are frozen at the requested equilibrium and frequency.
pub struct MaterialMechanicalHarmonic<'f, 'c, 's> {
    boundaries: Vec<CylindricalMaterialBoundary<'s>>,
    wall_band_hz: [f64; 2],
    interior: PorousRegionsInterior<'f, 'c, 's>,
}

impl<'f, 'c, 's> MaterialMechanicalHarmonic<'f, 'c, 's> {
    /// Prepare the full pressure, local membrane and hereditary wall equations.
    ///
    /// Wall geometry/source admission precedes volume material queries. Each
    /// source contributes surface-specific real coefficients and its complete
    /// positive relaxation spectrum, under its original frequency band. No
    /// reciprocal impedance or instantaneous replacement resistance is formed.
    /// The volume owner separately admits fluid/porous domains at this same
    /// frequency. The shared dense budget counts every original coordinate.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        geometry: &'f HelmholtzFe<'c>,
        frequency_hz: f64,
        sources: &[AcousticVolumeRegion<'s>],
        element_regions: &[usize],
        admittances: &[BoundaryAdmittance],
        boundaries: Vec<CylindricalMaterialBoundary<'s>>,
        pressure_release_faces: &[[u32; 3]],
        controls: EquivalentFluidControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialMechanicalHarmonicError> {
        if cancelled() {
            return Err(MaterialMechanicalError::Interior(HeterogeneousError::Cancelled).into());
        }
        let first = boundaries
            .first()
            .ok_or(MaterialMechanicalError::Geometry {
                boundary: 0,
                what: "at least one original material boundary is required",
            })?;
        let mut wall_band_hz = first.source.band_hz();
        let mut numerical = Vec::with_capacity(boundaries.len());
        for (index, boundary) in boundaries.iter().enumerate() {
            if cancelled() {
                return Err(
                    MaterialMechanicalError::Interior(HeterogeneousError::Cancelled).into(),
                );
            }
            let band = boundary.source.band_hz();
            if !frequency_hz.is_finite() || frequency_hz < band[0] || frequency_hz > band[1] {
                return Err(MaterialMechanicalHarmonicError::Frequency {
                    boundary: index,
                    frequency_hz,
                    source_band_hz: band,
                });
            }
            let span = boundary.admit_geometry(index, geometry, cancelled)?;
            let wall = boundary
                .source
                .acoustic_cell_wall(boundary.radius_m, span)
                .map_err(|source| MaterialMechanicalError::Source {
                    boundary: index,
                    source,
                })?;
            numerical.push(MechanicalBoundary {
                faces: boundary.faces.clone(),
                wall: wall.wall,
                relaxation: wall.relaxation,
            });
            wall_band_hz[0] = wall_band_hz[0].max(band[0]);
            wall_band_hz[1] = wall_band_hz[1].min(band[1]);
        }
        let interior = PorousRegionsInterior::new_with_coupled_mechanical_walls(
            geometry,
            frequency_hz,
            sources,
            element_regions,
            admittances,
            &numerical,
            pressure_release_faces,
            controls,
            cancelled,
        )?;
        Ok(Self {
            boundaries,
            wall_band_hz,
            interior,
        })
    }

    /// Original material and geometry bindings in native wall index order.
    #[must_use]
    pub fn boundaries(&self) -> &[CylindricalMaterialBoundary<'s>] {
        &self.boundaries
    }

    /// Intersection of original WALL applicability bands. Volume source
    /// applicability is admitted independently at the preparation frequency.
    #[must_use]
    pub const fn wall_band_hz(&self) -> [f64; 2] {
        self.wall_band_hz
    }

    /// Original source regions, mapping and retained material responses.
    #[must_use]
    pub const fn interior(&self) -> &PorousRegionsInterior<'f, 'c, 's> {
        &self.interior
    }

    /// The complete native solve, wall forces and independent physical losses.
    #[must_use]
    pub const fn prepared(&self) -> &PreparedEquivalentFluid<'f, 'c> {
        self.interior.operator()
    }
}
