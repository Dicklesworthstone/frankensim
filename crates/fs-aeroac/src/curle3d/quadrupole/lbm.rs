//! Resolved Reynolds stresses from the original stationary D3Q19 fluid.
//!
//! This adapter samples accepted `BoundaryGrid3` populations, including the
//! native half-force velocity correction. The original `MaterialScaling`
//! supplies reference density, dx and dt. Each selected full fluid cell emits
//! `rho u_i u_j dx^3` in a fixed physical frame; the six components remain
//! separate until the original time history is transformed and propagated.
//!
//! It is a cell-center volume quadrature over a caller-selected box, not an
//! automatic integration of every physical source. Grid index zero is located
//! at the explicitly authored cell center; no half-cell shift is inferred.
//! No periodic images, mean subtraction, filtering, SGS/viscous/thermodynamic
//! stress, moving-frame velocity or partially filled cell volume is supplied.
//! Use a stationary, single-phase, unsaturated voxel domain and refine the
//! source region, lattice and acoustic quadrature independently.

use super::recorded::QuadrupoleWelchStream;
use super::resolved_reynolds_moment_n_m;
use crate::AeroacError;
use fs_lbm::d3q19::BoundaryGrid3;
use fs_lbm::material::MaterialScaling;
use fs_material::fluid_transport::FluidTransport;

/// A stationary source frame and a half-open box of original lattice cells.
#[derive(Clone, Copy, Debug)]
pub struct LbmReynoldsVolumeConfig3 {
    /// Inclusive lower lattice indices [x, y, z].
    pub lower: [usize; 3],
    /// Exclusive upper lattice indices [x, y, z].
    pub upper: [usize; 3],
    /// Physical position of the center of grid cell [0, 0, 0] [m].
    /// Source positions are this point plus R * (dx * [x, y, z]).
    pub cell_zero_center_m: [f64; 3],
    /// Proper orthonormal map from grid Cartesian components to the fixed
    /// source/propagation frame. It rotates both position offsets and velocity.
    /// It does not subtract a mean flow or add a moving frame's velocity.
    pub grid_to_source: [[f64; 3]; 3],
    /// Maximum cells in the selected box, including solids. Bounds work on
    /// every sample and the retained solid-occupancy comparison.
    pub max_cells: usize,
    /// Maximum requested backing-buffer bytes, excluding Vec/allocator metadata.
    /// Includes occupancy, source positions and both moment buffers.
    pub max_workspace_bytes: usize,
}

/// Prepared spatial selection and reusable instantaneous stress observation.
///
/// Source order is z outermost, then y, then x, skipping only the solids
/// recorded during preparation. Each sample verifies dimensions and the
/// complete selected box's occupancy before publishing a new moment/time pair.
/// It does not claim object identity or topology outside that box. The native
/// grid exposes no clock or physical-unit binding: the caller must supply the
/// original scaling, stationary frame and actual accepted solver time.
pub struct PreparedLbmReynoldsVolume3 {
    config: LbmReynoldsVolumeConfig3,
    dimensions: [usize; 3],
    timestep_s: f64,
    density_unit_kg_m3: f64,
    velocity_unit_m_s: f64,
    cell_volume_m3: f64,
    solid: Vec<bool>,
    positions_m: Vec<[f64; 3]>,
    accepted: Vec<[f64; 6]>,
    scratch: Vec<[f64; 6]>,
    accepted_time_s: Option<f64>,
    workspace_bytes: usize,
}

impl PreparedLbmReynoldsVolume3 {
    /// Prepare full-cell geometry from the original material-to-lattice owner.
    ///
    /// This does not read a sample or advance the grid. Reference density is
    /// the physical density corresponding to lattice density one, not a
    /// replacement for each cell's actual density.
    pub fn new<F: FluidTransport + ?Sized>(
        grid: &BoundaryGrid3,
        scaling: &MaterialScaling<'_, F>,
        config: LbmReynoldsVolumeConfig3,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(grid, scaling, config, || false)
    }

    /// Cancellable preparation; no partial observer is returned on refusal.
    pub fn new_cancellable<F: FluidTransport + ?Sized>(
        grid: &BoundaryGrid3,
        scaling: &MaterialScaling<'_, F>,
        config: LbmReynoldsVolumeConfig3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let dimensions = grid.dimensions();
        let cells = config.cell_count(dimensions)?;
        admit_frame(config)?;
        let dx = scaling.cell_size_m();
        let cell_volume_m3 = dx * dx * dx;
        if !cell_volume_m3.is_finite() || cell_volume_m3 <= 0.0 {
            return Err(invalid("D3Q19 source cell volume is unrepresentable"));
        }
        let density_unit_kg_m3 = scaling
            .density_si(1.0)
            .map_err(|_| invalid("D3Q19 source density conversion is unrepresentable"))?;
        let velocity_unit_m_s = scaling
            .velocity_si(1.0)
            .map_err(|_| invalid("D3Q19 source velocity conversion is unrepresentable"))?;
        let mut sources = 0usize;
        for (index, [x, y, z]) in config.cells().enumerate() {
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            sources += usize::from(!grid.is_solid(x, y, z));
        }
        if sources == 0 {
            return Err(invalid("D3Q19 source box contains no full fluid cells"));
        }
        let bytes_per_source =
            core::mem::size_of::<[f64; 3]>() + 2 * core::mem::size_of::<[f64; 6]>();
        let workspace_bytes = sources
            .checked_mul(bytes_per_source)
            .and_then(|bytes| bytes.checked_add(cells * core::mem::size_of::<bool>()))
            .filter(|&bytes| bytes <= config.max_workspace_bytes)
            .ok_or_else(|| invalid("D3Q19 source buffers exceed the declared byte budget"))?;
        let mut solid = allocate(cells, false)?;
        let mut positions_m = allocate(sources, [0.0; 3])?;
        let mut source = 0;
        for (index, cell @ [x, y, z]) in config.cells().enumerate() {
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            solid[index] = grid.is_solid(x, y, z);
            if solid[index] {
                continue;
            }
            let offset = rotate(config.grid_to_source, cell.map(|i| i as f64 * dx));
            let position =
                core::array::from_fn(|axis| config.cell_zero_center_m[axis] + offset[axis]);
            if position.iter().any(|value| !value.is_finite()) {
                return Err(invalid("D3Q19 source position is unrepresentable"));
            }
            positions_m[source] = position;
            source += 1;
        }
        let accepted = allocate(sources, [0.0; 6])?;
        let scratch = allocate(sources, [0.0; 6])?;
        poll(&mut cancelled)?;
        Ok(Self {
            config,
            dimensions,
            timestep_s: scaling.timestep_s(),
            density_unit_kg_m3,
            velocity_unit_m_s,
            cell_volume_m3,
            solid,
            positions_m,
            accepted,
            scratch,
            accepted_time_s: None,
            workspace_bytes,
        })
    }

    /// Fixed positions in the same order as every accepted six-component sample.
    #[must_use]
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }

    /// Original full-cell source count; solids cannot appear/disappear later.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.positions_m.len()
    }

    /// Original lattice timestep in physical seconds. Sampling does not advance
    /// it or infer skipped solver steps; the downstream recorder admits cadence.
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.timestep_s
    }

    /// Physical timestamp of the last successfully published sample.
    #[must_use]
    pub const fn sample_time_s(&self) -> Option<f64> {
        self.accepted_time_s
    }

    /// Accepted instantaneous xx, yy, zz, xy, xz, yz moments [N m].
    /// No zero-valued physical sample is fabricated before the first sample.
    #[must_use]
    pub fn stress_moments_n_m(&self) -> Option<&[[f64; 6]]> {
        self.accepted_time_s.map(|_| self.accepted.as_slice())
    }

    /// Requested backing-buffer payload, excluding allocator metadata.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Feed the last accepted native sample into the original joint recorder.
    ///
    /// The recorder must retain this lattice's original dt, exact geometry and
    /// source order. Its own clock admits the timestamp and refuses skipped or
    /// duplicated observations. A refusal advances neither recorder count nor
    /// segment state, and leaves this accepted sample available for retry.
    pub fn record_sample(
        &self,
        recorder: &mut QuadrupoleWelchStream,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(&mut cancelled)?;
        if recorder.clock().timestep_s().to_bits() != self.timestep_s.to_bits() {
            return Err(invalid("D3Q19 source and recording timesteps differ"));
        }
        let time_s = self
            .accepted_time_s
            .ok_or_else(|| invalid("D3Q19 source has no accepted sample to record"))?;
        recorder.push_sample(&self.positions_m, time_s, &self.accepted, cancelled)
    }

    /// Observe the actual accepted grid at an explicitly supplied physical time.
    /// Times must be finite and strictly increase after the first observation.
    /// The grid is never mutated, and no temporal mean/filter is applied.
    pub fn sample(
        &mut self,
        grid: &BoundaryGrid3,
        time_s: f64,
    ) -> Result<&[[f64; 6]], AeroacError> {
        self.sample_cancellable(grid, time_s, || false)
    }

    /// Allocation-free, cancellable publication of one complete spatial sample.
    ///
    /// Refusal preserves accepted moments and time; reusable scratch may change.
    /// A later, separately called recorder may refuse this accepted observation
    /// without rolling it back. Consumers must record only successful samples.
    pub fn sample_cancellable(
        &mut self,
        grid: &BoundaryGrid3,
        time_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[[f64; 6]], AeroacError> {
        poll(&mut cancelled)?;
        if !time_s.is_finite()
            || self
                .accepted_time_s
                .is_some_and(|previous| time_s <= previous)
        {
            return Err(invalid(
                "D3Q19 source time must be finite and strictly increasing",
            ));
        }
        if grid.dimensions() != self.dimensions {
            return Err(invalid("D3Q19 source grid dimensions changed"));
        }
        let mut source = 0;
        for (index, [x, y, z]) in self.config.cells().enumerate() {
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            let solid = grid.is_solid(x, y, z);
            if solid != self.solid[index] {
                return Err(invalid("D3Q19 source box solid occupancy changed"));
            }
            if solid {
                continue;
            }
            let (lattice_density, lattice_velocity) = grid
                .macroscopic(x, y, z)
                .map_err(|_| invalid("D3Q19 source cell lost positive finite density"))?;
            let density = lattice_density * self.density_unit_kg_m3;
            if !density.is_finite() || density <= 0.0 {
                return Err(invalid("D3Q19 source density must be positive and finite"));
            }
            // The native owner alone supplies the Guo half-force correction.
            // Applying it again here would change both diagonal and shear terms.
            let velocity = rotate(self.config.grid_to_source, lattice_velocity)
                .map(|component| component * self.velocity_unit_m_s);
            self.scratch[source] =
                resolved_reynolds_moment_n_m(density, self.cell_volume_m3, velocity)?;
            source += 1;
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.accepted, &mut self.scratch);
        self.accepted_time_s = Some(time_s);
        Ok(&self.accepted)
    }
}

impl LbmReynoldsVolumeConfig3 {
    fn cell_count(self, dimensions: [usize; 3]) -> Result<usize, AeroacError> {
        let mut cells = 1usize;
        for axis in 0..3 {
            if self.lower[axis] >= self.upper[axis] || self.upper[axis] > dimensions[axis] {
                return Err(invalid(
                    "D3Q19 source box must be nonempty and inside the grid",
                ));
            }
            cells = cells
                .checked_mul(self.upper[axis] - self.lower[axis])
                .ok_or_else(|| invalid("D3Q19 source cell count overflow"))?;
        }
        if cells > self.max_cells {
            return Err(invalid(
                "D3Q19 source box exceeds the declared cell-work budget",
            ));
        }
        Ok(cells)
    }

    fn cells(self) -> impl Iterator<Item = [usize; 3]> {
        (self.lower[2]..self.upper[2]).flat_map(move |z| {
            (self.lower[1]..self.upper[1])
                .flat_map(move |y| (self.lower[0]..self.upper[0]).map(move |x| [x, y, z]))
        })
    }
}

fn admit_frame(config: LbmReynoldsVolumeConfig3) -> Result<(), AeroacError> {
    let r = config.grid_to_source;
    if config
        .cell_zero_center_m
        .iter()
        .any(|value| !value.is_finite())
        || r.iter().flatten().any(|value| !value.is_finite())
    {
        return Err(invalid("D3Q19 source frame must be finite"));
    }
    for i in 0..3 {
        for j in 0..3 {
            let dot = r[i][0] * r[j][0] + r[i][1] * r[j][1] + r[i][2] * r[j][2];
            let target = if i == j { 1.0 } else { 0.0 };
            if (dot - target).abs() > 64.0 * f64::EPSILON {
                return Err(invalid("D3Q19 source frame rotation must be orthonormal"));
            }
        }
    }
    let determinant = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    if (determinant - 1.0).abs() > 128.0 * f64::EPSILON {
        return Err(invalid(
            "D3Q19 source frame rotation must preserve orientation",
        ));
    }
    Ok(())
}

fn rotate(r: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    r.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2])
}

fn allocate<T: Clone>(length: usize, value: T) -> Result<Vec<T>, AeroacError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| invalid("D3Q19 source buffer allocation failed"))?;
    values.resize(length, value);
    Ok(values)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "D3Q19 resolved volume observation",
        })
    } else {
        Ok(())
    }
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

#[cfg(test)]
mod tests;
