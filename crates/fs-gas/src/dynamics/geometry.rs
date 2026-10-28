//! Stationary quasi-one-dimensional control-volume geometry.

use super::{Error, MAX_EULER_CELLS, allocated, positive};
use fs_qty::Length;

/// Axial faces, physical face areas, and integrated physical cell volumes.
///
/// The quasi-one-dimensional model assumes smoothly varying area and one
/// transverse equilibrium state per section. This geometry alone does not
/// supply the loss or mixing law of an abrupt junction or side branch.
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh1d {
    positions: Vec<f64>,
    areas: Vec<f64>,
    volumes: Vec<f64>,
    uniform_width: Option<f64>,
    uniform_area: Option<f64>,
}
impl Mesh1d {
    /// Admit a stationary mesh in m, m², and m³.
    ///
    /// Volumes must be integrated from the original area geometry; midpoint
    /// area times width is not exact for a conical cell. The geometry producer
    /// owns that quadrature/model accuracy. All arrays are retained unchanged.
    pub fn new(positions: Vec<f64>, areas: Vec<f64>, volumes: Vec<f64>) -> Result<Self, Error> {
        let n = volumes.len();
        if !(2..=MAX_EULER_CELLS).contains(&n) || positions.len() != n + 1 || areas.len() != n + 1 {
            return Err(Error::Invalid(
                "mesh needs 2..=MAX_EULER_CELLS volumes and one more face",
            ));
        }
        if positions.iter().any(|x| !x.is_finite())
            || positions.windows(2).any(|x| !positive(x[1] - x[0]))
            || areas.iter().chain(&volumes).any(|x| !positive(*x))
            || !positive(positions[n] - positions[0])
            || !positive(volumes.iter().sum())
        {
            return Err(Error::Invalid(
                "mesh geometry must be finite, ordered, and positive",
            ));
        }
        let first_width = positions[1] - positions[0];
        let uniform_width = positions
            .windows(2)
            .all(|x| {
                let width = x[1] - x[0];
                (width - first_width).abs() <= 32.0 * f64::EPSILON * first_width.max(width)
            })
            .then_some(first_width);
        let uniform_area = areas.iter().all(|a| *a == areas[0]).then_some(areas[0]);
        Ok(Self {
            positions,
            areas,
            volumes,
            uniform_width,
            uniform_area,
        })
    }
    /// Constant-area mesh with a declared uniform cell width, starting at zero.
    pub fn uniform(cells: usize, width: Length, area: f64) -> Result<Self, Error> {
        if !(2..=MAX_EULER_CELLS).contains(&cells)
            || !positive(width.value())
            || !positive(area)
            || !positive(width.value() * area)
            || !(width.value() * cells as f64).is_finite()
        {
            return Err(Error::Invalid("uniform mesh dimensions"));
        }
        let mut positions = allocated(cells + 1, 0.0)?;
        for (i, x) in positions.iter_mut().enumerate() {
            *x = width.value() * i as f64;
        }
        let mut mesh = Self::new(
            positions,
            allocated(cells + 1, area)?,
            allocated(cells, width.value() * area)?,
        )?;
        // Retain the declared width rather than a rounded coordinate difference.
        mesh.uniform_width = Some(width.value());
        Ok(mesh)
    }
    /// Face positions in metres, in strictly increasing order.
    pub fn positions(&self) -> &[f64] {
        &self.positions
    }
    /// Physical cross-sectional areas at the faces [m²].
    pub fn face_areas(&self) -> &[f64] {
        &self.areas
    }
    /// Integrated physical volumes [m³].
    pub fn volumes(&self) -> &[f64] {
        &self.volumes
    }
    /// Number of volume averages.
    pub fn cell_count(&self) -> usize {
        self.volumes.len()
    }
    /// Uniform width when all cells share one width, otherwise None.
    pub fn uniform_cell_width(&self) -> Option<Length> {
        self.uniform_width.map(Length::new)
    }
    /// Constant face area when present, otherwise None.
    pub const fn uniform_cross_section(&self) -> Option<f64> {
        self.uniform_area
    }
    pub(super) fn center(&self, i: usize) -> f64 {
        self.positions[i] + 0.5 * (self.positions[i + 1] - self.positions[i])
    }
    pub(super) fn width(&self, i: usize) -> f64 {
        self.positions[i + 1] - self.positions[i]
    }
}
