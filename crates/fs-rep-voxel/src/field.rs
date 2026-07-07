//! Voxel fields on the shared sparse substrate: boolean occupancy,
//! material ids, and SIMP density fractions — with active-set boolean
//! algebra and morphology (dilate/erode/open/close), the manufacturing-
//! constraint primitives. All operations are deterministic (BTreeMap
//! substrate ⇒ sorted iteration).

use crate::VoxelError;
use core::fmt;
use fs_rep_sdf::VdbGrid;

fn clone_grid<T: Copy>(grid: &VdbGrid<T>, background: T) -> VdbGrid<T> {
    // VdbGrid carries no derives; rebuilding from the active set is the
    // substrate-honest clone (deterministic: BTreeMap iteration order).
    let mut out = VdbGrid::new(background);
    for (c, v) in grid.iter_active() {
        out.set(c, v);
    }
    out
}

/// Boolean occupancy on the sparse tree, with world-space placement.
pub struct OccupancyField {
    /// The sparse active set (active = solid).
    pub grid: VdbGrid<bool>,
    /// Voxel edge length (m).
    pub voxel_size: f64,
    /// World position of voxel (0,0,0)'s min corner.
    pub origin: [f64; 3],
}

impl Clone for OccupancyField {
    fn clone(&self) -> Self {
        OccupancyField {
            grid: clone_grid(&self.grid, false),
            voxel_size: self.voxel_size,
            origin: self.origin,
        }
    }
}

impl fmt::Debug for OccupancyField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "OccupancyField {{ active: {}, voxel_size: {} }}",
            self.active(),
            self.voxel_size
        )
    }
}

impl OccupancyField {
    /// An empty field.
    ///
    /// # Errors
    /// [`VoxelError::Parameters`] on a non-positive voxel size.
    pub fn new(voxel_size: f64, origin: [f64; 3]) -> Result<Self, VoxelError> {
        if !(voxel_size.is_finite() && voxel_size > 0.0) {
            return Err(VoxelError::Parameters {
                what: format!("voxel size {voxel_size} must be positive"),
            });
        }
        Ok(OccupancyField {
            grid: VdbGrid::new(false),
            voxel_size,
            origin,
        })
    }

    /// Activate one voxel.
    pub fn set(&mut self, coord: [i32; 3]) {
        self.grid.set(coord, true);
    }

    /// Is a voxel solid?
    #[must_use]
    pub fn is_solid(&self, coord: [i32; 3]) -> bool {
        self.grid.is_active(coord)
    }

    /// The world-space center of a voxel.
    #[must_use]
    pub fn center(&self, coord: [i32; 3]) -> [f64; 3] {
        core::array::from_fn(|k| self.origin[k] + (f64::from(coord[k]) + 0.5) * self.voxel_size)
    }

    /// The voxel containing a world point.
    #[must_use]
    pub fn voxel_of(&self, p: [f64; 3]) -> [i32; 3] {
        core::array::from_fn(|k| {
            #[allow(clippy::cast_possible_truncation)]
            {
                ((p[k] - self.origin[k]) / self.voxel_size).floor() as i32
            }
        })
    }

    /// Active-set union (in place).
    pub fn union(&mut self, other: &OccupancyField) {
        for (c, _) in other.grid.iter_active() {
            self.grid.set(c, true);
        }
    }

    /// Active-set intersection (in place).
    pub fn intersect(&mut self, other: &OccupancyField) {
        let doomed: Vec<[i32; 3]> = self
            .grid
            .iter_active()
            .map(|(c, _)| c)
            .filter(|&c| !other.grid.is_active(c))
            .collect();
        for c in doomed {
            self.grid.deactivate(c);
        }
    }

    /// Active-set subtraction (in place).
    pub fn subtract(&mut self, other: &OccupancyField) {
        let doomed: Vec<[i32; 3]> = self
            .grid
            .iter_active()
            .map(|(c, _)| c)
            .filter(|&c| other.grid.is_active(c))
            .collect();
        for c in doomed {
            self.grid.deactivate(c);
        }
    }

    /// Morphological dilation by `n` voxels (6-connected).
    pub fn dilate(&mut self, n: u32) {
        for _ in 0..n {
            self.grid.dilate();
        }
    }

    /// Morphological erosion by `n` voxels (6-connected).
    pub fn erode(&mut self, n: u32) {
        for _ in 0..n {
            self.grid.erode();
        }
    }

    /// Morphological opening (erode then dilate): removes features
    /// thinner than `2n` voxels — the min-feature-size primitive.
    pub fn open(&mut self, n: u32) {
        self.erode(n);
        self.dilate(n);
    }

    /// Morphological closing (dilate then erode): fills gaps thinner than
    /// `2n` voxels.
    pub fn close(&mut self, n: u32) {
        self.dilate(n);
        self.erode(n);
    }

    /// Active-voxel count.
    #[must_use]
    pub fn active(&self) -> u64 {
        self.grid.active_count()
    }

    /// Active-set statistics as a JSON line (ledger logs).
    #[must_use]
    pub fn stats_json(&self) -> String {
        format!(
            "{{\"kind\":\"occupancy-stats\",\"active\":{},\"voxel_size\":{}}}",
            self.active(),
            self.voxel_size
        )
    }
}

/// Material ids on the sparse tree (0 is reserved for "empty").
pub struct MaterialField {
    /// Material id per active voxel.
    pub grid: VdbGrid<u16>,
}

impl Clone for MaterialField {
    fn clone(&self) -> Self {
        MaterialField {
            grid: clone_grid(&self.grid, 0),
        }
    }
}

impl fmt::Debug for MaterialField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MaterialField {{ active: {} }}", self.grid.active_count())
    }
}

impl Default for MaterialField {
    fn default() -> Self {
        MaterialField::new()
    }
}

impl MaterialField {
    /// An empty multi-material field.
    #[must_use]
    pub fn new() -> Self {
        MaterialField {
            grid: VdbGrid::new(0),
        }
    }

    /// Assign a material (`id > 0`) to a voxel.
    ///
    /// # Errors
    /// [`VoxelError::Parameters`] for id 0 (reserved for empty).
    pub fn assign(&mut self, coord: [i32; 3], id: u16) -> Result<(), VoxelError> {
        if id == 0 {
            return Err(VoxelError::Parameters {
                what: "material id 0 is reserved for empty space".to_string(),
            });
        }
        self.grid.set(coord, id);
        Ok(())
    }

    /// The material at a voxel (0 = empty).
    #[must_use]
    pub fn material(&self, coord: [i32; 3]) -> u16 {
        self.grid.get(coord)
    }

    /// The occupancy footprint of one material.
    ///
    /// # Errors
    /// [`VoxelError::Parameters`] via [`OccupancyField::new`].
    pub fn occupancy_of(
        &self,
        id: u16,
        voxel_size: f64,
        origin: [f64; 3],
    ) -> Result<OccupancyField, VoxelError> {
        let mut occ = OccupancyField::new(voxel_size, origin)?;
        for (c, m) in self.grid.iter_active() {
            if m == id {
                occ.set(c);
            }
        }
        Ok(occ)
    }
}

/// SIMP density fractions in `[0, 1]` on the sparse tree.
pub struct DensityField {
    /// Fraction per active voxel (background 0).
    pub grid: VdbGrid<f32>,
}

impl Clone for DensityField {
    fn clone(&self) -> Self {
        DensityField {
            grid: clone_grid(&self.grid, 0.0),
        }
    }
}

impl fmt::Debug for DensityField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "DensityField {{ active: {} }}", self.grid.active_count())
    }
}

impl Default for DensityField {
    fn default() -> Self {
        DensityField::new()
    }
}

impl DensityField {
    /// An empty density field.
    #[must_use]
    pub fn new() -> Self {
        DensityField {
            grid: VdbGrid::new(0.0),
        }
    }

    /// Set a fraction (clamps are REFUSED, not silent).
    ///
    /// # Errors
    /// [`VoxelError::Parameters`] outside `[0, 1]`.
    pub fn set(&mut self, coord: [i32; 3], fraction: f32) -> Result<(), VoxelError> {
        if !(0.0..=1.0).contains(&fraction) || !fraction.is_finite() {
            return Err(VoxelError::Parameters {
                what: format!("density fraction {fraction} outside [0, 1]"),
            });
        }
        self.grid.set(coord, fraction);
        Ok(())
    }

    /// Threshold to occupancy at a cutoff fraction.
    ///
    /// # Errors
    /// [`VoxelError::Parameters`] via [`OccupancyField::new`].
    pub fn threshold(
        &self,
        cutoff: f32,
        voxel_size: f64,
        origin: [f64; 3],
    ) -> Result<OccupancyField, VoxelError> {
        let mut occ = OccupancyField::new(voxel_size, origin)?;
        for (c, f) in self.grid.iter_active() {
            if f >= cutoff {
                occ.set(c);
            }
        }
        Ok(occ)
    }
}
