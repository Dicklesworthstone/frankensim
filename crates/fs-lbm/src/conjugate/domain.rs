//! Voxel domain, fluid properties, and solid materials for the conjugate
//! pipeline. Cell `(x, y, z)` has linear index `(z ny + y) nx + x` and
//! occupies `[X_x, X_(x+1)] x [Y_y, Y_(y+1)] x [Z_z, Z_(z+1)]` metres, where
//! the face coordinates are `i dx` on a uniform domain and declared per axis
//! on a graded one ([`VoxelDomain::graded`]: a non-uniform Cartesian grid,
//! fine where the geometry or the boundary layers need it).

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

/// Cartesian voxel domain, uniform or graded per axis.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelDomain {
    nx: usize,
    ny: usize,
    nz: usize,
    /// Uniform spacing (NaN on a graded domain).
    dx: f64,
    /// Face coordinates per axis (`n + 1` each, increasing from 0) on a
    /// graded domain; `None` when uniform.
    faces: Option<Box<[Vec<f64>; 3]>>,
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
            faces: None,
            voxels: vec![Voxel::Fluid; cells],
        })
    }

    /// An all-fluid graded domain from its cell widths along each axis
    /// (metres; faces start at 0). Equal widths along every axis give the
    /// uniform domain.
    ///
    /// # Errors
    /// [`ChtError::InvalidDomain`] for an empty axis, a non-finite or
    /// non-positive width, or an overflowing cell count.
    pub fn graded(widths: [Vec<f64>; 3]) -> Result<Self, ChtError> {
        let dims = [widths[0].len(), widths[1].len(), widths[2].len()];
        for (axis, list) in widths.iter().enumerate() {
            if let Some(w) = list.iter().find(|w| !(w.is_finite() && **w > 0.0)) {
                return Err(ChtError::InvalidDomain {
                    reason: format!(
                        "axis {axis}: cell widths must be finite and positive, got {w}"
                    ),
                });
            }
        }
        let first = widths[0].first().copied().unwrap_or(1.0);
        let uniform = widths.iter().all(|list| list.iter().all(|w| *w == first));
        let mut domain = Self::new(dims[0], dims[1], dims[2], first)?;
        if !uniform {
            domain.dx = f64::NAN;
            domain.faces = Some(Box::new(widths.map(|list| {
                let mut faces = Vec::with_capacity(list.len() + 1);
                faces.push(0.0);
                for w in list {
                    faces.push(faces.last().copied().unwrap_or(0.0) + w);
                }
                faces
            })));
        }
        Ok(domain)
    }

    /// A graded domain whose occupancy is sampled at every cell centre.
    ///
    /// # Errors
    /// As [`VoxelDomain::graded`].
    pub fn graded_from_fn(
        widths: [Vec<f64>; 3],
        mut occupancy: impl FnMut([f64; 3]) -> Voxel,
    ) -> Result<Self, ChtError> {
        let mut domain = Self::graded(widths)?;
        for c in 0..domain.cell_count() {
            let [x, y, z] = domain.coords(c);
            domain.voxels[c] = occupancy(domain.center(x, y, z));
        }
        Ok(domain)
    }

    /// Whether every cell is the same cube.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.faces.is_none()
    }

    /// The uniform spacing, or a refusal naming `what` needs it (the
    /// lattice-Boltzmann paths).
    pub(crate) fn require_uniform(&self, what: &str) -> Result<f64, ChtError> {
        if self.is_uniform() {
            Ok(self.dx)
        } else {
            Err(ChtError::InvalidDomain {
                reason: format!("{what} needs a uniform voxel grid; this domain is graded"),
            })
        }
    }

    /// Width of cell index `i` along `axis`, metres.
    #[inline]
    #[must_use]
    pub fn width(&self, axis: usize, i: usize) -> f64 {
        match &self.faces {
            None => self.dx,
            Some(faces) => faces[axis][i + 1] - faces[axis][i],
        }
    }

    /// Coordinate of face plane `i` (`0..=n`) along `axis`, metres.
    #[inline]
    #[must_use]
    pub fn face_coord(&self, axis: usize, i: usize) -> f64 {
        match &self.faces {
            None => i as f64 * self.dx,
            Some(faces) => faces[axis][i],
        }
    }

    /// Widths of cell `c` along x, y, z, metres.
    #[inline]
    #[must_use]
    pub fn widths(&self, c: usize) -> [f64; 3] {
        let at = self.coords(c);
        [0, 1, 2].map(|a| self.width(a, at[a]))
    }

    /// Area of cell `c`'s faces normal to `axis`, m^2.
    #[inline]
    #[must_use]
    pub fn face_area(&self, c: usize, axis: usize) -> f64 {
        let w = self.widths(c);
        w[(axis + 1) % 3] * w[(axis + 2) % 3]
    }

    /// Volume of cell `c`, m^3.
    #[inline]
    #[must_use]
    pub fn volume(&self, c: usize) -> f64 {
        let w = self.widths(c);
        w[0] * w[1] * w[2]
    }

    /// Smallest cell width over the domain, metres.
    #[must_use]
    pub fn min_width(&self) -> f64 {
        match &self.faces {
            None => self.dx,
            Some(faces) => faces
                .iter()
                .flat_map(|f| f.windows(2).map(|w| w[1] - w[0]))
                .fold(f64::INFINITY, f64::min),
        }
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

    /// Cell edge length of a uniform domain, metres (NaN on a graded one:
    /// geometry there goes through [`Self::width`] and its relatives, so a
    /// spacing assumption that slipped through fails loudly).
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
        match &self.faces {
            None => [
                (x as f64 + 0.5) * self.dx,
                (y as f64 + 0.5) * self.dx,
                (z as f64 + 0.5) * self.dx,
            ],
            Some(faces) => {
                let at = [x, y, z];
                [0, 1, 2].map(|a| 0.5 * (faces[a][at[a]] + faces[a][at[a] + 1]))
            }
        }
    }

    /// Neighbour of `c` across local face `face` (0..6 in [`super::Face3::ALL`]
    /// order), or `None` at the domain boundary. This is geometric adjacency;
    /// the caller checks fluid/solid occupancy separately.
    ///
    /// # Panics
    /// If `c` is outside the domain or `face` is outside `0..6`.
    #[inline]
    #[must_use]
    pub fn neighbor(&self, c: usize, face: usize) -> Option<usize> {
        assert!(c < self.cell_count() && face < 6, "cell/face index outside domain");
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

/// One conducting solid: isotropic, or orthotropic along the grid axes.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidMaterial {
    /// Caller label retained in reports (for example a matdb card id).
    pub label: String,
    /// Thermal conductivity, W/(m K) (the isotropic value; superseded per
    /// axis by `orthotropic_w_m_k` when declared).
    pub conductivity_w_m_k: f64,
    /// Principal conductivities along x, y, z, W/(m K), for orthotropic
    /// solids aligned with the grid (for example a PCB laminate: high
    /// in-plane, low through the board).
    pub orthotropic_w_m_k: Option<[f64; 3]>,
    /// Volumetric heat capacity `rho c`, J/(m^3 K). Required only by
    /// transient marches; steady solves never read it.
    pub volumetric_heat_capacity_j_m3_k: Option<f64>,
    /// Temperature dependence `[(T K, k W/(m K)), ...]` (empty: constant):
    /// piecewise linear in T, constant beyond the ends, scaling every axis
    /// conductivity by `k(T) / conductivity_w_m_k`. Continuous same-material
    /// faces and Dirichlet/Robin half cells integrate the table exactly;
    /// different materials and contacts share the flux from their two
    /// integrated half-cell laws and the declared interface resistance.
    /// The energy solve iterates to consistency.
    pub conductivity_table: Vec<(f64, f64)>,
}

impl SolidMaterial {
    /// A labelled isotropic solid without a declared heat capacity.
    #[must_use]
    pub fn new(label: impl Into<String>, conductivity_w_m_k: f64) -> Self {
        Self {
            label: label.into(),
            conductivity_w_m_k,
            orthotropic_w_m_k: None,
            volumetric_heat_capacity_j_m3_k: None,
            conductivity_table: Vec::new(),
        }
    }

    /// Declare a temperature-dependent conductivity `[(T K, k W/(m K)), ...]`
    /// (see [`Self::conductivity_table`]).
    #[must_use]
    pub fn with_conductivity_table(mut self, points: &[(f64, f64)]) -> Self {
        self.conductivity_table = points.to_vec();
        self
    }

    /// `k(T) / conductivity_w_m_k`: the factor scaling every axis
    /// conductivity at temperature `t` (1 without a table).
    #[must_use]
    pub fn conductivity_factor(&self, t: f64) -> f64 {
        let table = &self.conductivity_table;
        let Some(&(_, k0)) = table.first() else {
            return 1.0;
        };
        let upper = table.partition_point(|&(temperature, _)| temperature < t);
        let k = if upper == 0 {
            k0
        } else if upper == table.len() {
            table[table.len() - 1].1
        } else {
            let (a, b) = (table[upper - 1], table[upper]);
            a.1 + (b.1 - a.1) * (t - a.0) / (b.0 - a.0)
        };
        k / self.conductivity_w_m_k
    }

    /// Mean of the declared conductivity factor over a temperature interval.
    /// Integrating each crossed linear piece exactly avoids subtracting large
    /// Kirchhoff potentials when the two temperatures are nearly equal. The
    /// table retains its declared constant extensions beyond the end knots.
    pub(super) fn conductivity_secant_factor(&self, a: f64, b: f64) -> f64 {
        if a == b || self.conductivity_table.is_empty() {
            return self.conductivity_factor(a);
        }
        let (lo, hi) = (a.min(b), a.max(b));
        let span = hi - lo;
        let mut left = lo;
        let mut left_factor = self.conductivity_factor(left);
        let mut mean = 0.0;
        let first = self
            .conductivity_table
            .partition_point(|&(temperature, _)| temperature <= lo);
        for &(temperature, conductivity) in &self.conductivity_table[first..] {
            if temperature >= hi {
                break;
            }
            let right_factor = conductivity / self.conductivity_w_m_k;
            mean += (temperature - left) / span * (0.5 * left_factor + 0.5 * right_factor);
            left = temperature;
            left_factor = right_factor;
        }
        mean + (hi - left) / span * (0.5 * left_factor + 0.5 * self.conductivity_factor(hi))
    }

    /// Conductivity factor averaged from a cell to its solved Robin surface.
    /// With z the temperature drop from cell toward ambient, the unique root
    /// satisfies integral(0..z, k + d*h) = d*h*abs(Tcell-Tambient). Each table
    /// piece gives a quadratic with a positive derivative. Resistance scaling
    /// keeps its coefficients bounded, and the rationalized root avoids
    /// subtracting nearly equal roots. NaN refuses unrepresentable input via
    /// the energy assembly's existing finite-coefficient check.
    pub(super) fn film_conductivity_factor(
        &self,
        cell: f64,
        ambient: f64,
        axis_k: f64,
        distance: f64,
        h: f64,
    ) -> f64 {
        if cell == ambient || self.conductivity_table.is_empty() {
            return self.conductivity_factor(cell);
        }
        let (lo, hi) = (cell.min(ambient), cell.max(ambient));
        let span = hi - lo;
        let table = &self.conductivity_table;
        let first = table.partition_point(|&(temperature, _)| temperature <= lo);
        let end = table.partition_point(|&(temperature, _)| temperature < hi);
        let cell_factor = self.conductivity_factor(cell);
        let ambient_factor = self.conductivity_factor(ambient);
        let mut maximum = cell_factor.max(ambient_factor);
        for &(_, k) in &table[first..end] {
            maximum = maximum.max(k / self.conductivity_w_m_k);
        }
        let maximum_k = axis_k * maximum;
        let bulk = distance / maximum_k;
        let film = 1.0 / h;
        let resistance_scale = bulk.max(film);
        if !span.is_finite()
            || !maximum_k.is_finite()
            || maximum_k <= 0.0
            || !resistance_scale.is_finite()
            || resistance_scale <= 0.0
        {
            return f64::NAN;
        }
        // alpha*k/kmax + beta is (k+d*h) divided by a common scale.
        let alpha = film / resistance_scale;
        let beta = bulk / resistance_scale;
        if beta == 0.0 {
            return cell_factor;
        }
        let mut remaining = beta;
        let mut start = cell;
        let mut k_start = cell_factor / maximum;
        let mut covered = 0.0;
        let mut mean = 0.0;
        let knots = end - first;
        for piece in 0..=knots {
            let (finish, k_finish) = if piece == knots {
                (ambient, ambient_factor / maximum)
            } else {
                let index = if cell < ambient {
                    first + piece
                } else {
                    end - 1 - piece
                };
                let (temperature, conductivity) = table[index];
                (
                    temperature,
                    conductivity / self.conductivity_w_m_k / maximum,
                )
            };
            let fraction = (finish - start).abs() / span;
            let a = beta + alpha * k_start;
            let b = 0.5 * alpha * (k_finish - k_start);
            let average = beta + 0.5 * alpha * (k_start + k_finish);
            let integral = fraction * average;
            if remaining <= integral || piece == knots {
                // Positive conductivity places the exact root in the closed
                // interval. Clamp only accumulated endpoint roundoff when the
                // film tends to an imposed surface temperature.
                let y = (remaining / fraction).min(average);
                let discriminant = if b < 0.0 {
                    // Equivalent to a*a+4*b*y, but all terms are positive
                    // even when the endpoint slope is much smaller than a.
                    let endpoint = beta + alpha * k_finish;
                    endpoint * endpoint + (-4.0 * b) * (average - y).max(0.0)
                } else {
                    a * a + 4.0 * b * y
                };
                let position = 2.0 * y / (a + discriminant.sqrt());
                if !position.is_finite() {
                    return f64::NAN;
                }
                let position = position.clamp(0.0, 1.0);
                let used = fraction * position;
                let local_mean = k_start + 0.5 * position * (k_finish - k_start);
                // Retain the fractional drop: an absolute surface temperature
                // can round back to the cell even when its mean k differs.
                return maximum
                    * if covered == 0.0 {
                        local_mean
                    } else {
                        let total = covered + used;
                        covered / total * mean + used / total * local_mean
                    };
            }
            if fraction > 0.0 {
                let total = covered + fraction;
                mean = covered / total * mean + fraction / total * (0.5 * k_start + 0.5 * k_finish);
                covered = total;
            }
            remaining -= integral;
            start = finish;
            k_start = k_finish;
        }
        f64::NAN
    }

    /// Declare grid-aligned principal conductivities `[k_x, k_y, k_z]`,
    /// W/(m K).
    #[must_use]
    pub fn with_orthotropic(mut self, conductivity_w_m_k: [f64; 3]) -> Self {
        self.orthotropic_w_m_k = Some(conductivity_w_m_k);
        self
    }

    /// Conductivity along each grid axis, W/(m K).
    #[must_use]
    pub fn axis_conductivity(&self) -> [f64; 3] {
        self.orthotropic_w_m_k
            .unwrap_or([self.conductivity_w_m_k; 3])
    }

    /// Declare the volumetric heat capacity `rho c`, J/(m^3 K).
    #[must_use]
    pub fn with_heat_capacity(mut self, volumetric_heat_capacity_j_m3_k: f64) -> Self {
        self.volumetric_heat_capacity_j_m3_k = Some(volumetric_heat_capacity_j_m3_k);
        self
    }
}
