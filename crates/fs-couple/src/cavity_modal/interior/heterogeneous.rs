//! Frozen source gas in the variable-density three-dimensional pressure equation.
//!
//! The axial coordinate is a physical distance embedded in the original tetrahedral
//! geometry. Temperature and mole fractions are evaluated by `AxialGasField` at
//! each volume quadrature point; density and sound speed are never interpolated
//! independently. The resulting sparse operators can serve many frequencies and
//! drives without repeating material queries. This is a stationary, lossless,
//! isotropic fluid volume with no mean flow or thermal evolution.

use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_material::MaterialError;
use fs_material::gas_field::{AxialGasField, AxialGasSampling};

/// Original source refusals and errors in the physical embedding remain distinct.
#[derive(Debug)]
pub enum AxialGasInteriorError {
    /// Nonfinite geometry, a zero direction, or a nonrepresentable projection.
    Embedding(&'static str),
    /// A mesh vertex lies outside the original field's supplied interval.
    OutsideSource {
        vertex: usize,
        coordinate_m: f64,
        length_m: f64,
    },
    /// A cellwise material jump cuts an element and requires a conforming mesh.
    UnresolvedJump { tetrahedron: usize, station_m: f64 },
    /// Preparation was cancelled before publishing an operator.
    Cancelled,
    /// The original material/phase owner refused a physical sample.
    Material(MaterialError),
    /// The original finite-element owner refused assembly.
    Interior(AcousticFeError),
}

impl core::fmt::Display for AxialGasInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Embedding(what) => write!(f, "axial acoustic embedding: {what}"),
            Self::OutsideSource {
                vertex,
                coordinate_m,
                length_m,
            } => write!(
                f,
                "acoustic mesh vertex {vertex} projects to {coordinate_m} m outside [0, {length_m}] m"
            ),
            Self::UnresolvedJump {
                tetrahedron,
                station_m,
            } => write!(
                f,
                "gas jump at {station_m} m cuts acoustic tetrahedron {tetrahedron}; conform the mesh to the interface"
            ),
            Self::Cancelled => write!(f, "source acoustic preparation cancelled"),
            Self::Material(error) => error.fmt(f),
            Self::Interior(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for AxialGasInteriorError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Interior(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MaterialError> for AxialGasInteriorError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}

impl From<AcousticFeError> for AxialGasInteriorError {
    fn from(error: AcousticFeError) -> Self {
        Self::Interior(error)
    }
}

/// An original gas field extruded along planes normal to an explicit axis.
/// Coordinates use metres: `x = dot(point - origin, unit_direction)`.
/// Normalizing the direction changes no field length or physical wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialGasEmbedding {
    origin_m: [f64; 3],
    unit_direction: [f64; 3],
}

impl AxialGasEmbedding {
    pub fn new(origin_m: [f64; 3], direction: [f64; 3]) -> Result<Self, AxialGasInteriorError> {
        if origin_m
            .iter()
            .chain(&direction)
            .any(|value| !value.is_finite())
        {
            return Err(AxialGasInteriorError::Embedding(
                "origin and direction must be finite",
            ));
        }
        let scale = direction.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        if scale == 0.0 {
            return Err(AxialGasInteriorError::Embedding(
                "direction must be nonzero",
            ));
        }
        let scaled = direction.map(|value| value / scale);
        let norm = fs_math::det::sqrt(scaled.iter().map(|value| value * value).sum());
        Ok(Self {
            origin_m,
            unit_direction: scaled.map(|value| value / norm),
        })
    }

    #[must_use]
    pub const fn origin_m(self) -> [f64; 3] {
        self.origin_m
    }

    #[must_use]
    pub const fn unit_direction(self) -> [f64; 3] {
        self.unit_direction
    }

    /// Project a physical position without clamping or extrapolating a source.
    pub fn coordinate_m(self, point_m: [f64; 3]) -> Result<f64, AxialGasInteriorError> {
        let coordinate = (0..3)
            .map(|axis| (point_m[axis] - self.origin_m[axis]) * self.unit_direction[axis])
            .sum::<f64>();
        if !coordinate.is_finite() {
            return Err(AxialGasInteriorError::Embedding(
                "projected coordinate is not finite",
            ));
        }
        Ok(coordinate)
    }
}

/// Source-retaining preparation of the original spatial pressure space.
///
/// The homogeneous `geometry.medium()` is not a reference fluid for this owner:
/// every volume coefficient comes from `source`. Its pressure traces, boundary
/// normals, element basis and topology are reused unchanged. The numerical owner
/// exposes rigid-wall modes and full harmonic pressure solves with passive Robin
/// admittance and prescribed inward velocity, including the original work balance.
pub struct AxialGasInterior<'f, 'c, 's> {
    source: &'s AxialGasField,
    embedding: AxialGasEmbedding,
    operator: HeterogeneousHelmholtzFe<'f, 'c>,
}

impl<'f, 'c, 's> AxialGasInterior<'f, 'c, 's> {
    /// Assemble once at the caller's volume quadrature order.
    ///
    /// Every mesh vertex must lie in the supplied source interval. A true jump
    /// of a cellwise source must lie on element faces; unresolved jumps refuse
    /// instead of being aliased by quadrature. Continuous primitive gradients
    /// use the requested quadrature and require the caller's mesh/order study.
    pub fn new(
        geometry: &'f HelmholtzFe<'c>,
        source: &'s AxialGasField,
        embedding: AxialGasEmbedding,
        quadrature_order: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AxialGasInteriorError> {
        let mut projected = Vec::with_capacity(geometry.positions().len());
        for (vertex, &point) in geometry.positions().iter().enumerate() {
            if cancelled() {
                return Err(AxialGasInteriorError::Cancelled);
            }
            let coordinate_m = embedding.coordinate_m(point)?;
            if coordinate_m < 0.0 || coordinate_m > source.length_m() {
                return Err(AxialGasInteriorError::OutsideSource {
                    vertex,
                    coordinate_m,
                    length_m: source.length_m(),
                });
            }
            projected.push(coordinate_m);
        }
        if source.sampling() == AxialGasSampling::CellwiseConstant {
            // State i occupies (station[i-1], station[i]]; state 0 is the inlet
            // value only. Equal adjacent primitives do not introduce an interface.
            let jumps: Vec<_> = (1..source.stations_m().len() - 1)
                .filter(|&i| source.primitive_states()[i] != source.primitive_states()[i + 1])
                .map(|i| source.stations_m()[i])
                .collect();
            for (tetrahedron, vertices) in geometry.tetrahedra().iter().enumerate() {
                if cancelled() {
                    return Err(AxialGasInteriorError::Cancelled);
                }
                let coordinates = vertices.map(|vertex| projected[vertex as usize]);
                let low = coordinates.iter().copied().fold(f64::INFINITY, f64::min);
                let high = coordinates
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let next = jumps.partition_point(|station| *station <= low);
                if let Some(&station_m) = jumps.get(next).filter(|station| **station < high) {
                    return Err(AxialGasInteriorError::UnresolvedJump {
                        tetrahedron,
                        station_m,
                    });
                }
            }
        }
        let uniform = source.uniform_state().map(|state| Medium {
            density_kg_m3: state.state().density,
            sound_speed_m_s: state.state().sound_speed,
        });
        let mut sample = |_: usize, point: [f64; 3]| -> Result<Medium, AxialGasInteriorError> {
            if let Some(medium) = uniform {
                return Ok(medium);
            }
            let state = source.state_at(embedding.coordinate_m(point)?)?;
            Ok(Medium {
                density_kg_m3: state.state().density,
                sound_speed_m_s: state.state().sound_speed,
            })
        };
        let operator =
            HeterogeneousHelmholtzFe::new(geometry, quadrature_order, &mut sample, cancelled)?;
        Ok(Self {
            source,
            embedding,
            operator,
        })
    }

    /// Original primitive field, source cards, pressure and phase policy.
    #[must_use]
    pub const fn source(&self) -> &'s AxialGasField {
        self.source
    }

    #[must_use]
    pub const fn embedding(&self) -> AxialGasEmbedding {
        self.embedding
    }

    /// Reuse the prepared spatial operator for modes and frequency-domain drives.
    #[must_use]
    pub const fn operator(&self) -> &HeterogeneousHelmholtzFe<'f, 'c> {
        &self.operator
    }
}
