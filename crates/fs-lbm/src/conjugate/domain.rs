//! Voxel domain, fluid properties, and solid materials for the conjugate
//! pipeline. Cell `(x, y, z)` occupies `[x dx, (x+1) dx] x [y dx, (y+1) dx] x
//! [z dx, (z+1) dx]` metres and has linear index `(z ny + y) nx + x`.

use super::{ChtError, finite_positive};

/// Occupancy of one voxel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Voxel {
    /// Carries the declared fluid; participates in flow and energy.
    Fluid,
    /// Conducting solid with the given index into the caller's
    /// [`SolidMaterial`] table; participates in energy only.
    Solid(u16),
}

/// Uniform Cartesian voxel domain.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelDomain {
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    voxels: Vec<Voxel>,
}

impl VoxelDomain {
    /// An all-fluid domain of `nx x ny x nz` cubic cells of edge `dx` metres.
    ///
    /// # Errors
    /// [`ChtError::InvalidDomain`] for a zero dimension, an overflowing cell
    /// count, or a non-finite / non-positive spacing.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Result<Self, ChtError> {
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(ChtError::InvalidDomain {
                reason: format!("dimensions must be positive, got {nx}x{ny}x{nz}"),
            });
        }
        if !(dx.is_finite() && dx > 0.0) {
            return Err(ChtError::InvalidDomain {
                reason: format!("cell spacing must be finite and positive, got {dx}"),
            });
        }
        let cells = nx
            .checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .ok_or_else(|| ChtError::InvalidDomain {
                reason: "cell count overflows".into(),
            })?;
        Ok(Self {
            nx,
            ny,
            nz,
            dx,
            voxels: vec![Voxel::Fluid; cells],
        })
    }

    /// A domain whose occupancy is sampled at every cell centre (metres).
    ///
    /// # Errors
    /// As [`VoxelDomain::new`].
    pub fn from_fn(
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        mut occupancy: impl FnMut([f64; 3]) -> Voxel,
    ) -> Result<Self, ChtError> {
        let mut domain = Self::new(nx, ny, nz, dx)?;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let c = domain.index(x, y, z);
                    domain.voxels[c] = occupancy(domain.center(x, y, z));
                }
            }
        }
        Ok(domain)
    }

    /// Cell counts along `x`, `y`, `z`.
    #[must_use]
    pub const fn dims(&self) -> [usize; 3] {
        [self.nx, self.ny, self.nz]
    }

    /// Cell edge length, metres.
    #[must_use]
    pub const fn dx(&self) -> f64 {
        self.dx
    }

    /// Total number of cells.
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.voxels.len()
    }

    /// Number of fluid cells.
    #[must_use]
    pub fn fluid_count(&self) -> usize {
        self.voxels
            .iter()
            .filter(|v| matches!(v, Voxel::Fluid))
            .count()
    }

    /// Linear index of `(x, y, z)`.
    ///
    /// # Panics
    /// If the coordinate is outside the domain.
    #[inline]
    #[must_use]
    pub fn index(&self, x: usize, y: usize, z: usize) -> usize {
        assert!(
            x < self.nx && y < self.ny && z < self.nz,
            "cell ({x},{y},{z}) outside domain"
        );
        (z * self.ny + y) * self.nx + x
    }

    /// Coordinates of linear index `c`.
    #[inline]
    #[must_use]
    pub fn coords(&self, c: usize) -> [usize; 3] {
        let x = c % self.nx;
        let yz = c / self.nx;
        [x, yz % self.ny, yz / self.ny]
    }

    /// Occupancy of `(x, y, z)`.
    #[inline]
    #[must_use]
    pub fn voxel(&self, x: usize, y: usize, z: usize) -> Voxel {
        self.voxels[self.index(x, y, z)]
    }

    /// Occupancy by linear index.
    #[inline]
    #[must_use]
    pub fn voxel_at(&self, c: usize) -> Voxel {
        self.voxels[c]
    }

    /// Whether cell `c` is fluid.
    #[inline]
    #[must_use]
    pub fn is_fluid(&self, c: usize) -> bool {
        matches!(self.voxels[c], Voxel::Fluid)
    }

    /// Overwrite one cell's occupancy.
    pub fn set(&mut self, x: usize, y: usize, z: usize, voxel: Voxel) {
        let c = self.index(x, y, z);
        self.voxels[c] = voxel;
    }

    /// Set every cell whose centre satisfies `inside` to `voxel`; returns the
    /// number of cells written.
    pub fn fill_where(&mut self, voxel: Voxel, mut inside: impl FnMut([f64; 3]) -> bool) -> usize {
        let mut written = 0;
        for c in 0..self.voxels.len() {
            let [x, y, z] = self.coords(c);
            if inside(self.center(x, y, z)) {
                self.voxels[c] = voxel;
                written += 1;
            }
        }
        written
    }

    /// Cell-centre coordinates in metres.
    #[must_use]
    pub fn center(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            (x as f64 + 0.5) * self.dx,
            (y as f64 + 0.5) * self.dx,
            (z as f64 + 0.5) * self.dx,
        ]
    }

    /// Neighbour of `c` across local face `face` (0..6 in [`super::Face3::ALL`]
    /// order), or `None` at the domain boundary.
    #[inline]
    pub(crate) fn neighbor(&self, c: usize, face: usize) -> Option<usize> {
        let [x, y, z] = self.coords(c);
        match face {
            0 => (x > 0).then(|| c - 1),
            1 => (x + 1 < self.nx).then(|| c + 1),
            2 => (y > 0).then(|| c - self.nx),
            3 => (y + 1 < self.ny).then(|| c + self.nx),
            4 => (z > 0).then(|| c - self.nx * self.ny),
            _ => (z + 1 < self.nz).then(|| c + self.nx * self.ny),
        }
    }

    /// Validate every solid material index against a table of `materials`.
    pub(crate) fn check_materials(&self, materials: usize) -> Result<(), ChtError> {
        for (cell, voxel) in self.voxels.iter().enumerate() {
            if let Voxel::Solid(m) = *voxel
                && usize::from(m) >= materials
            {
                return Err(ChtError::UnknownMaterial { cell, material: m });
            }
        }
        Ok(())
    }
}

/// Constant fluid properties at one declared state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluidProperties {
    /// Density, kg/m^3.
    pub density_kg_m3: f64,
    /// Isobaric specific heat, J/(kg K).
    pub specific_heat_j_kg_k: f64,
    /// Thermal conductivity, W/(m K).
    pub conductivity_w_m_k: f64,
    /// Kinematic viscosity, m^2/s.
    pub kinematic_viscosity_m2_s: f64,
}

impl FluidProperties {
    /// Dry air at 300 K and 1 atm (Incropera, DeWitt, Bergman & Lavine,
    /// *Fundamentals of Heat and Mass Transfer*, 6th ed., Table A.4):
    /// `rho = 1.1614`, `c_p = 1007`, `k = 0.0263`, `nu = 15.89e-6`.
    #[must_use]
    pub const fn dry_air_300k() -> Self {
        Self {
            density_kg_m3: 1.1614,
            specific_heat_j_kg_k: 1007.0,
            conductivity_w_m_k: 0.0263,
            kinematic_viscosity_m2_s: 15.89e-6,
        }
    }

    /// Volumetric heat capacity `rho c_p`, J/(m^3 K).
    #[must_use]
    pub fn volumetric_heat_capacity(&self) -> f64 {
        self.density_kg_m3 * self.specific_heat_j_kg_k
    }

    /// Thermal diffusivity `k / (rho c_p)`, m^2/s.
    #[must_use]
    pub fn diffusivity(&self) -> f64 {
        self.conductivity_w_m_k / self.volumetric_heat_capacity()
    }

    /// Prandtl number `nu / alpha`.
    #[must_use]
    pub fn prandtl(&self) -> f64 {
        self.kinematic_viscosity_m2_s / self.diffusivity()
    }

    pub(crate) fn validate(&self) -> Result<(), ChtError> {
        finite_positive("fluid.density_kg_m3", self.density_kg_m3)?;
        finite_positive("fluid.specific_heat_j_kg_k", self.specific_heat_j_kg_k)?;
        finite_positive("fluid.conductivity_w_m_k", self.conductivity_w_m_k)?;
        finite_positive(
            "fluid.kinematic_viscosity_m2_s",
            self.kinematic_viscosity_m2_s,
        )
    }
}

/// One isotropic conducting solid.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidMaterial {
    /// Caller label retained in reports (for example a matdb card id).
    pub label: String,
    /// Thermal conductivity, W/(m K).
    pub conductivity_w_m_k: f64,
}

impl SolidMaterial {
    /// A labelled isotropic solid.
    #[must_use]
    pub fn new(label: impl Into<String>, conductivity_w_m_k: f64) -> Self {
        Self {
            label: label.into(),
            conductivity_w_m_k,
        }
    }
}
