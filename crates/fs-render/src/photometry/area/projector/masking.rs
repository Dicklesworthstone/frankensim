//! Passive geometric apertures baked into the original represented gobo.
//!
//! Clipping integrates the open part of each output texel. The resulting map
//! is still piecewise constant: its integrated transmission is retained, but
//! its within-cell silhouette is a declared raster approximation. Diffraction,
//! blade thickness, calibration and a total physical error bound are absent.

use core::f64::consts::TAU;
use core::fmt;
use std::sync::Arc;

use fs_exec::{Cancelled, Cx};
use fs_math::det;

use super::{GoboTransmission, MAX_GOBO_TEXELS, SourceError};

/// Maximum number of sides in a declared regular-polygon iris.
pub const MAX_IRIS_BLADES: usize = 64;
/// Maximum independent framing half-planes in one gobo mask.
pub const MAX_FRAMING_SHUTTERS: usize = 4;
const MAX_PLANES: usize = MAX_IRIS_BLADES + MAX_FRAMING_SHUTTERS;
const MAX_VERTICES: usize = 4 + MAX_PLANES;

/// Centered regular polygon with a physical inradius in metres.
///
/// Side normals have angles `rotation + 2*pi*k/blades`; the open region obeys
/// `normal.dot(point) <= inradius` for every side. This is a polygonal iris,
/// not a circle approximation with an implied circular calibration. Zero
/// inradius explicitly closes the aperture.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolygonalIris {
    inradius_m: f64,
    blades: usize,
    rotation_rad: f64,
    rotation: [f64; 2],
}

impl PolygonalIris {
    /// Admit nonnegative finite radius, 3--64 sides, and finite rotation.
    pub fn try_new(inradius_m: f64, blades: usize, rotation_rad: f64) -> Result<Self, SourceError> {
        if !inradius_m.is_finite() || inradius_m < 0.0 || !(3..=MAX_IRIS_BLADES).contains(&blades) {
            return Err(SourceError::InvalidDistribution);
        }
        Ok(Self {
            inradius_m,
            blades,
            rotation_rad,
            rotation: normal(rotation_rad)?,
        })
    }

    /// Perpendicular distance from the optical origin to every side, in metres.
    #[must_use]
    pub const fn inradius_m(&self) -> f64 {
        self.inradius_m
    }

    /// Number of sides in the declared aperture.
    #[must_use]
    pub const fn blades(&self) -> usize {
        self.blades
    }

    /// Angle of the first outward side normal, in radians.
    #[must_use]
    pub const fn rotation_rad(&self) -> f64 {
        self.rotation_rad
    }
}

/// An opaque blade whose open half-plane is `normal.dot(point) <= offset`.
/// The normal is `(cos(angle), sin(angle))` in the original gobo X/Y frame;
/// the signed offset is in metres and can place the blade past the origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FramingShutter {
    normal_angle_rad: f64,
    plane: HalfPlane,
}

impl FramingShutter {
    /// Admit a finite normal angle and signed physical offset.
    pub fn try_new(normal_angle_rad: f64, offset_m: f64) -> Result<Self, SourceError> {
        if !offset_m.is_finite() {
            return Err(SourceError::InvalidDistribution);
        }
        Ok(Self {
            normal_angle_rad,
            plane: HalfPlane {
                normal: normal(normal_angle_rad)?,
                offset_m,
            },
        })
    }

    /// Outward normal angle in radians.
    #[must_use]
    pub const fn normal_angle_rad(&self) -> f64 {
        self.normal_angle_rad
    }

    /// Signed normal offset in metres.
    #[must_use]
    pub const fn offset_m(&self) -> f64 {
        self.plane.offset_m
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct HalfPlane {
    normal: [f64; 2],
    offset_m: f64,
}

/// Immutable intersection of an optional polygonal iris and up to four blades.
/// Clones share the admitted geometry; no physical light-field evaluation or
/// mutable global cache is attached to this value.
#[derive(Debug, Clone, PartialEq)]
pub struct GoboMask {
    data: Arc<MaskData>,
}

#[derive(Debug, PartialEq)]
struct MaskData {
    iris: Option<PolygonalIris>,
    shutters: Vec<FramingShutter>,
    planes: Vec<HalfPlane>,
}

impl GoboMask {
    /// Intersect the supplied apertures in stable iris-side, then blade order.
    pub fn try_new(
        iris: Option<PolygonalIris>,
        shutters: Vec<FramingShutter>,
    ) -> Result<Self, SourceError> {
        if shutters.len() > MAX_FRAMING_SHUTTERS {
            return Err(SourceError::DistributionTooLarge);
        }
        let mut planes = Vec::new();
        planes
            .try_reserve_exact(iris.map_or(0, |value| value.blades) + shutters.len())
            .map_err(|_| SourceError::DistributionTooLarge)?;
        if let Some(iris) = iris {
            for side in 0..iris.blades {
                let local = normal(TAU * side as f64 / iris.blades as f64)?;
                // Rotate a bounded side angle instead of adding it to a huge
                // supplied angle, which could round every side to one normal.
                let rotated = [
                    iris.rotation[0] * local[0] - iris.rotation[1] * local[1],
                    iris.rotation[1] * local[0] + iris.rotation[0] * local[1],
                ];
                planes.push(HalfPlane {
                    normal: unit_normal(rotated)?,
                    offset_m: iris.inradius_m,
                });
            }
        }
        planes.extend(shutters.iter().map(|shutter| shutter.plane));
        Ok(Self {
            data: Arc::new(MaskData {
                iris,
                shutters,
                planes,
            }),
        })
    }

    /// Optional declared polygonal iris.
    #[must_use]
    pub fn iris(&self) -> Option<&PolygonalIris> {
        self.data.iris.as_ref()
    }

    /// Blades in their deterministic clipping order.
    #[must_use]
    pub fn shutters(&self) -> &[FramingShutter] {
        &self.data.shutters
    }

    /// Whether this mask leaves the original gobo unobstructed.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.planes.is_empty()
    }

    /// Transient assembly-cache identity, never an output or semantic identity.
    pub(crate) fn storage_key(&self) -> usize {
        Arc::as_ptr(&self.data) as usize
    }
}

/// Cancellation or an established native source-admission refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GoboMaskError {
    /// The complete replacement map was not published.
    Cancelled,
    /// Invalid, oversized, unallocatable or unrepresentable source data.
    Source(SourceError),
}

impl From<Cancelled> for GoboMaskError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl From<SourceError> for GoboMaskError {
    fn from(error: SourceError) -> Self {
        Self::Source(error)
    }
}

impl fmt::Display for GoboMaskError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "photometric gobo masking refused: {self:?}")
    }
}

impl std::error::Error for GoboMaskError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cancelled => None,
            Self::Source(error) => Some(error),
        }
    }
}

impl GoboTransmission {
    /// Integrate a passive aperture into this map at integer texel subdivisions.
    ///
    /// Physical width/height remain unchanged. Every output cell stays inside
    /// one authored texel and receives its original transmission multiplied by
    /// the fraction left open by all half-planes. This is an area average, not
    /// a center sample. Subdivision refines silhouettes without renormalizing
    /// lost flux; within-cell raster and floating-point errors are unqualified.
    ///
    /// Output is capped by [`MAX_GOBO_TEXELS`], storage is reserved fallibly,
    /// and cancellation is checked per cell. An unrepresentably narrow positive
    /// iris or grid step refuses instead of becoming physical opacity. A closed
    /// iris is valid. Empty masks with unit subdivisions retain the original Arc.
    pub fn masked(
        &self,
        mask: &GoboMask,
        subdivisions: [usize; 2],
        cx: &Cx<'_>,
    ) -> Result<Self, GoboMaskError> {
        cx.checkpoint()?;
        if subdivisions.contains(&0) {
            return Err(SourceError::InvalidDistribution.into());
        }
        let columns = self
            .columns()
            .checked_mul(subdivisions[0])
            .ok_or(SourceError::DistributionTooLarge)?;
        let rows = self
            .rows()
            .checked_mul(subdivisions[1])
            .ok_or(SourceError::DistributionTooLarge)?;
        let count = columns
            .checked_mul(rows)
            .filter(|&value| value <= MAX_GOBO_TEXELS)
            .ok_or(SourceError::DistributionTooLarge)?;
        if mask.is_empty() && subdivisions == [1, 1] {
            return Ok(self.clone());
        }
        let dimensions = [self.width_m(), self.height_m()];
        let resolution = [columns as f64, rows as f64];
        for (size, cells) in dimensions.into_iter().zip(resolution) {
            let step = size / cells;
            if !step.is_finite() || step <= 0.0 || size * 0.5 - step == size * 0.5 {
                return Err(SourceError::NumericalRange.into());
            }
            if let Some(iris) = mask.iris()
                && iris.inradius_m > 0.0
                && iris.inradius_m < step
            {
                let width = iris.inradius_m / step;
                if width == 0.0 || 0.5 + width == 0.5 || 0.5 - width == 0.5 {
                    return Err(SourceError::NumericalRange.into());
                }
            }
        }
        let closed = mask.iris().is_some_and(|iris| iris.inradius_m == 0.0);
        let mut prepared = [CellPlane::EMPTY; MAX_PLANES];
        if !closed {
            for (target, plane) in prepared.iter_mut().zip(&mask.data.planes) {
                cx.checkpoint()?;
                *target = CellPlane::try_new(*plane, dimensions, resolution)?;
            }
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| SourceError::DistributionTooLarge)?;
        for row in 0..rows {
            cx.checkpoint()?;
            let center_y = (row as f64 + 0.5) / resolution[1] - 0.5;
            let original_row = row / subdivisions[1];
            for column in 0..columns {
                cx.checkpoint()?;
                let original =
                    self.data.values[original_row * self.columns() + column / subdivisions[0]];
                let fraction = if closed || original == 0.0 {
                    0.0
                } else {
                    open_fraction(
                        &prepared[..mask.data.planes.len()],
                        [(column as f64 + 0.5) / resolution[0] - 0.5, center_y],
                    )?
                };
                let value = original * fraction;
                if !value.is_finite()
                    || !(0.0..=1.0).contains(&value)
                    || (original > 0.0 && fraction > 0.0 && value == 0.0)
                {
                    return Err(SourceError::NumericalRange.into());
                }
                values.push(value);
            }
        }
        cx.checkpoint()?;
        let result = Self::try_new(dimensions[0], dimensions[1], columns, rows, values)?;
        cx.checkpoint()?;
        Ok(result)
    }
}

fn normal(angle: f64) -> Result<[f64; 2], SourceError> {
    if !angle.is_finite() {
        return Err(SourceError::InvalidDistribution);
    }
    let (sine, cosine) = det::sin_cos(angle);
    unit_normal([cosine, sine])
}

fn unit_normal(value: [f64; 2]) -> Result<[f64; 2], SourceError> {
    let norm = det::sqrt(value[0] * value[0] + value[1] * value[1]);
    if !norm.is_finite() || norm <= 0.0 {
        return Err(SourceError::NumericalRange);
    }
    Ok(value.map(|component| component / norm))
}

#[derive(Clone, Copy)]
struct CellPlane {
    global: [f64; 2],
    local: [f64; 2],
    offset: f64,
}

impl CellPlane {
    const EMPTY: Self = Self {
        global: [0.0; 2],
        local: [0.0; 2],
        offset: 0.0,
    };

    fn try_new(
        plane: HalfPlane,
        dimensions: [f64; 2],
        resolution: [f64; 2],
    ) -> Result<Self, SourceError> {
        let metric = [
            plane.normal[0] * dimensions[0],
            plane.normal[1] * dimensions[1],
        ];
        let scale = metric[0]
            .abs()
            .max(metric[1].abs())
            .max(plane.offset_m.abs());
        if !scale.is_finite() || scale <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        let global = metric.map(|value| value / scale);
        let local = [global[0] / resolution[0], global[1] / resolution[1]];
        let offset = plane.offset_m / scale;
        if (0..2).any(|axis| {
            plane.normal[axis] != 0.0
                && (metric[axis] == 0.0 || global[axis] == 0.0 || local[axis] == 0.0)
        }) || (plane.offset_m != 0.0 && offset == 0.0)
        {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            global,
            local,
            offset,
        })
    }
}

fn open_fraction(planes: &[CellPlane], center: [f64; 2]) -> Result<f64, SourceError> {
    let mut first = [[0.0; 2]; MAX_VERTICES];
    let mut second = [[0.0; 2]; MAX_VERTICES];
    first[..4].copy_from_slice(&[[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]);
    let (mut input, mut output) = (&mut first, &mut second);
    let mut count = 4;
    for plane in planes {
        let offset = plane.offset - plane.global[0] * center[0] - plane.global[1] * center[1];
        let radius = 0.5 * (plane.local[0].abs() + plane.local[1].abs());
        if radius <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        if offset >= radius {
            continue;
        }
        if offset <= -radius {
            return Ok(0.0);
        }
        let scale = plane.local[0].abs().max(plane.local[1].abs());
        let normal = plane.local.map(|value| value / scale);
        count = clip(&input[..count], output, normal, offset / scale)?;
        if count < 3 {
            return Ok(0.0);
        }
        std::mem::swap(&mut input, &mut output);
    }
    let anchor = input[0];
    let (mut sum, mut correction) = (0.0, 0.0);
    for triangle in input[1..count].windows(2) {
        let a = [triangle[0][0] - anchor[0], triangle[0][1] - anchor[1]];
        let b = [triangle[1][0] - anchor[0], triangle[1][1] - anchor[1]];
        let cross = a[0] * b[1] - a[1] * b[0];
        let area = 0.5 * cross;
        if !area.is_finite() || area < 0.0 || (cross > 0.0 && area == 0.0) {
            return Err(SourceError::NumericalRange);
        }
        let next = sum + area;
        correction += if sum >= area {
            (sum - next) + area
        } else {
            (area - next) + sum
        };
        sum = next;
    }
    // The convex clipped polygon lies inside the unit cell. Positive-sum
    // rounding must not turn that construction into an amplifying filter.
    Ok((sum + correction).min(1.0))
}

fn clip(
    input: &[[f64; 2]],
    output: &mut [[f64; 2]; MAX_VERTICES],
    normal: [f64; 2],
    offset: f64,
) -> Result<usize, SourceError> {
    let distance = |point: [f64; 2]| offset - normal[0] * point[0] - normal[1] * point[1];
    let mut previous = input[input.len() - 1];
    let mut previous_distance = distance(previous);
    let mut count = 0;
    for &point in input {
        let point_distance = distance(point);
        if (previous_distance >= 0.0) != (point_distance >= 0.0) {
            let t = previous_distance / (previous_distance - point_distance);
            if !t.is_finite() || !(0.0..=1.0).contains(&t) {
                return Err(SourceError::NumericalRange);
            }
            let intersection = if t <= 0.5 {
                [
                    previous[0] + t * (point[0] - previous[0]),
                    previous[1] + t * (point[1] - previous[1]),
                ]
            } else {
                [
                    point[0] + (1.0 - t) * (previous[0] - point[0]),
                    point[1] + (1.0 - t) * (previous[1] - point[1]),
                ]
            };
            push_vertex(output, &mut count, intersection)?;
        }
        if point_distance >= 0.0 {
            push_vertex(output, &mut count, point)?;
        }
        previous = point;
        previous_distance = point_distance;
    }
    if count > 1 && output[0] == output[count - 1] {
        count -= 1;
    }
    Ok(count)
}

fn push_vertex(
    output: &mut [[f64; 2]; MAX_VERTICES],
    count: &mut usize,
    point: [f64; 2],
) -> Result<(), SourceError> {
    if *count > 0 && output[*count - 1] == point {
        return Ok(());
    }
    if *count == output.len() || point.iter().any(|value| !(-0.5..=0.5).contains(value)) {
        return Err(SourceError::NumericalRange);
    }
    output[*count] = point;
    *count += 1;
    Ok(())
}

#[cfg(test)]
mod tests;
