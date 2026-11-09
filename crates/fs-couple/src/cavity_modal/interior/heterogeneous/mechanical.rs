//! Original material walls in full spatial acoustic time and harmonic operators.
//!
//! The source law is a thin cylindrical hoop membrane. Its local radial motion
//! approximates the outward normal motion of explicitly admitted mesh facets;
//! it is not a flat-plate law. Every original Maxwell arm reaches the same
//! implicit fluid/wall solve. The FE trace mass supplies actual area exactly once.

use super::AxialGasEmbedding;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::wall::material::{LocalMaterialWall, WallMaterialLoss};
use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::{
    MechanicalBoundary, MechanicalTimeLimits, PreparedMechanicalTime,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{HeterogeneousError, HeterogeneousHelmholtzFe};

mod harmonic;
pub use harmonic::{MaterialMechanicalHarmonic, MaterialMechanicalHarmonicError};

/// One original wall source on a faceted approximation of a cylinder's side.
///
/// Radius and axis must describe the actual reference geometry. The explicit
/// tolerance bounds both relative radial-distance error and the Euclidean
/// difference between the facet's outward unit normal and each vertex's radial
/// unit direction. This admits a controlled faceting approximation, not curved
/// elements. End caps and inward-facing walls refuse. Callers must refine the
/// mesh and tolerance together for spatial accuracy.
///
/// Original thin-wall, axial restraint and hoop/axial material-axis assumptions
/// remain binding. There is no circumferential/axial shell coupling, bending,
/// pressure prestress, exterior loading, changing geometry or thermal evolution.
#[derive(Clone, Debug)]
pub struct CylindricalMaterialBoundary<'s> {
    pub faces: Vec<[u32; 3]>,
    pub source: &'s LocalMaterialWall,
    pub axis: AxialGasEmbedding,
    pub radius_m: f64,
    /// Finite, strictly between zero and one; never inferred from the mesh.
    pub relative_geometry_tolerance: f64,
}

/// Keep source and geometry refusals distinct from the original solver errors.
#[derive(Debug)]
pub enum MaterialMechanicalError {
    Geometry {
        boundary: usize,
        what: &'static str,
    },
    Source {
        boundary: usize,
        source: AcousticRealizeError,
    },
    Interior(HeterogeneousError),
}

impl core::fmt::Display for MaterialMechanicalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry { boundary, what } => write!(f, "material wall {boundary}: {what}"),
            Self::Source { boundary, source } => write!(f, "material wall {boundary}: {source}"),
            Self::Interior(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for MaterialMechanicalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Source { source, .. } => Some(source),
            Self::Interior(error) => Some(error),
            _ => None,
        }
    }
}
impl From<HeterogeneousError> for MaterialMechanicalError {
    fn from(error: HeterogeneousError) -> Self {
        Self::Interior(error)
    }
}

/// The exact source bindings and their one prepared acoustic/wall time owner.
///
/// The original fluid operator remains borrowed, including any gas/material
/// source binding that owns it. Wall source identities, bands and cylindrical
/// geometry remain attached for the lifetime of all borrowed transient states.
pub struct MaterialMechanicalTime<'m, 'f, 'c, 's> {
    boundaries: Vec<CylindricalMaterialBoundary<'s>>,
    band_hz: [f64; 2],
    prepared: PreparedMechanicalTime<'m, 'f, 'c>,
}

impl<'m, 'f, 'c, 's> MaterialMechanicalTime<'m, 'f, 'c, 's> {
    /// Cold preparation from complete original material laws on disjoint faces.
    ///
    /// Each source supplies mass, equilibrium stiffness, viscosity and every
    /// positive hereditary arm through its existing acoustic reduction. The
    /// actual axial span only admits that surface-specific law; its canonical
    /// cylindrical strip area is not applied to the FE boundary a second time.
    /// Original Maxwell clock limits apply to this midpoint solver. They are
    /// local resolution checks, not a global temporal error bound.
    ///
    /// The returned common source band is applicability, not an input filter.
    /// The caller must qualify the excited frequency range and independently
    /// refine time and space. Updating sources requires a new frozen binding;
    /// accepted wall histories are never implicitly transferred to another one.
    pub fn new(
        operator: &'m HeterogeneousHelmholtzFe<'f, 'c>,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        boundaries: Vec<CylindricalMaterialBoundary<'s>>,
        limits: MechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialMechanicalError> {
        if cancelled() {
            return Err(HeterogeneousError::Cancelled.into());
        }
        let first = boundaries
            .first()
            .ok_or(MaterialMechanicalError::Geometry {
                boundary: 0,
                what: "at least one original material boundary is required",
            })?;
        let mut band_hz = first.source.band_hz();
        let mut numerical = Vec::with_capacity(boundaries.len());
        for (index, boundary) in boundaries.iter().enumerate() {
            if cancelled() {
                return Err(HeterogeneousError::Cancelled.into());
            }
            let span = boundary.admit_geometry(index, operator.geometry(), cancelled)?;
            let material_error = |source| MaterialMechanicalError::Source {
                boundary: index,
                source,
            };
            let wall = boundary
                .source
                .acoustic_cell_wall(boundary.radius_m, span)
                .map_err(material_error)?;
            if let WallMaterialLoss::Maxwell(source) = boundary.source.loss() {
                source
                    .validate_local_time_step(wall.wall, duration_s)
                    .map_err(material_error)?;
            }
            let band = boundary.source.band_hz();
            band_hz[0] = band_hz[0].max(band[0]);
            band_hz[1] = band_hz[1].min(band[1]);
            numerical.push(MechanicalBoundary {
                faces: boundary.faces.clone(),
                wall: wall.wall,
                relaxation: wall.relaxation,
            });
        }
        if band_hz[0] >= band_hz[1] {
            return Err(MaterialMechanicalError::Geometry {
                boundary: 0,
                what: "wall sources have no common positive frequency interval",
            });
        }
        let prepared = operator
            .prepare_time_with_walls(duration_s, resistive, &numerical, limits, cancelled)?;
        Ok(Self {
            boundaries,
            band_hz,
            prepared,
        })
    }

    /// Original sources and explicit geometry, in native wall-index order.
    #[must_use]
    pub fn boundaries(&self) -> &[CylindricalMaterialBoundary<'s>] {
        &self.boundaries
    }
    /// Intersection of original source applicability intervals [Hz].
    #[must_use]
    pub const fn band_hz(&self) -> [f64; 2] {
        self.band_hz
    }
    /// The original sparse implicit solver, with all wall histories and work.
    #[must_use]
    pub const fn prepared(&self) -> &PreparedMechanicalTime<'m, 'f, 'c> {
        &self.prepared
    }
}

impl CylindricalMaterialBoundary<'_> {
    fn admit_geometry(
        &self,
        index: usize,
        geometry: &HelmholtzFe<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MaterialMechanicalError> {
        let invalid = |what| MaterialMechanicalError::Geometry {
            boundary: index,
            what,
        };
        let tolerance = self.relative_geometry_tolerance;
        if !self.radius_m.is_finite()
            || self.radius_m <= 0.0
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || tolerance >= 1.0
            || self.faces.is_empty()
        {
            return Err(invalid(
                "positive cylinder radius, nonempty faces and geometry tolerance in (0,1) required",
            ));
        }
        let origin = self.axis.origin_m();
        let direction = self.axis.unit_direction();
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for &face in &self.faces {
            if cancelled() {
                return Err(HeterogeneousError::Cancelled.into());
            }
            let (vertices, normal, _) = geometry
                .boundary_geometry(face)
                .map_err(HeterogeneousError::from)?;
            for vertex in vertices {
                let axial = self
                    .axis
                    .coordinate_m(vertex)
                    .map_err(|_| invalid("cylinder projection is unrepresentable"))?;
                low = low.min(axial);
                high = high.max(axial);
                let radial: [f64; 3] = core::array::from_fn(|i| {
                    (vertex[i] - origin[i] - axial * direction[i]) / self.radius_m
                });
                let radius = fs_math::det::sqrt(radial.iter().map(|r| r * r).sum());
                if !radius.is_finite() || radius <= 0.0 || (radius - 1.0).abs() > tolerance {
                    return Err(invalid(
                        "boundary vertex is outside the declared cylindrical side",
                    ));
                }
                let normal_error = (0..3)
                    .map(|i| {
                        let difference = normal[i] - radial[i] / radius;
                        difference * difference
                    })
                    .sum::<f64>();
                if !normal_error.is_finite() || normal_error > tolerance * tolerance {
                    return Err(invalid(
                        "boundary normal is not the admitted outward radial direction",
                    ));
                }
            }
        }
        let span = high - low;
        if !span.is_finite() || span <= 0.0 {
            return Err(invalid(
                "cylindrical side must have positive finite axial extent",
            ));
        }
        Ok(span)
    }
}
