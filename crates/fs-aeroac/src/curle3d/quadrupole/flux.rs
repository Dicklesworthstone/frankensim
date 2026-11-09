//! Resolved Lighthill sources from the original compressible solvers.
//!
//! Each selected, stationary Cartesian cell supplies
//! `Q_ij = [rho u_i u_j + ((p-p0)-c0^2(rho-rho0)) delta_ij] dV`.
//! Density, pressure and velocity come from the accepted native state and its
//! original caloric/species owner. No constant-gamma EOS replaces that owner.
//! The explicit acoustic medium supplies rho0 and c0; p0 must describe that
//! same medium. It is not guessed from the (possibly different) source gas.
//!
//! The Euler constructors retain these inviscid terms. Molecular constructors
//! subtract the original configured viscous/SGS stress, including explicit bulk
//! viscosity, reconstructed from the native conservative face tractions. They
//! cannot subsequently accept an inviscid sample that silently drops stress.
//! Applied body-force dipoles, boundary scattering and contributions outside
//! the selected volume are absent. Transform stress histories on their actual
//! accepted clocks before using a frequency-domain observer. Cell averaging
//! loses unresolved velocity products and does not qualify an LES noise model.
//!
//! Reference: Lighthill (1952), DOI 10.1098/rspa.1952.0060; the full tensor and
//! its ambient pressure-density difference appear in NASA report 19770026218,
//! equations (8)-(9):
//! <https://ntrs.nasa.gov/api/citations/19770026218/downloads/19770026218.pdf>.
//! Velocity is absolute in the propagation frame; no mean flow is subtracted.

use fs_flux::compressible::spatial::mixture::{MixtureEuler3d, MixtureGas};
use fs_flux::compressible::spatial::{Boundaries3d, CartesianGrid, Euler3d, FrozenGas, Primitive};
use fs_material::fluid_acoustics::FluidAcousticState;

mod kernel;
mod molecular;
/// Closed stationary permeable-surface mass and momentum sources.
pub mod surface;

/// Explicit stationary mapping into the acoustic propagation frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressibleVolumeConfig3 {
    /// Original ambient density and sound speed, retained without adjustment.
    pub acoustic: FluidAcousticState,
    /// Positive absolute pressure [Pa] of that same acoustic reference medium.
    pub reference_pressure_pa: f64,
    /// Location of the native grid's minimum corner in the acoustic frame [m].
    pub origin_m: [f64; 3],
    /// Proper orthonormal rotation, acoustic components = R * grid components.
    /// This is a stationary rotation, not a translating or rotating frame.
    pub grid_to_acoustic: [[f64; 3]; 3],
}

/// Native source errors remain distinguishable from observation admission.
#[derive(Debug, Clone, PartialEq)]
pub enum CompressibleVolumeError3 {
    Invalid(&'static str),
    Workspace {
        required_bytes: usize,
        maximum_bytes: usize,
    },
    Allocation,
    SourceMismatch,
    GeometryMismatch,
    /// Only strictly later native accepted times may replace a sample.
    SampleClock {
        previous_s: f64,
        candidate_s: f64,
    },
    Frozen(fs_flux::compressible::Error),
    Mixture(fs_flux::compressible::spatial::mixture::Error),
    Molecular(fs_flux::compressible::molecular::MolecularError),
    Dissipation(fs_flux::compressible::spatial::mixture::molecular::DissipationError),
    Cancelled,
}
impl core::fmt::Display for CompressibleVolumeError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "compressible volume acoustic source: {self:?}")
    }
}
impl std::error::Error for CompressibleVolumeError3 {}

#[derive(Debug)]
enum Source {
    Frozen(FrozenGas),
    Mixture(MixtureGas),
}

/// Reusable source extraction with atomic publication and an original clock.
///
/// Cells are a nonempty, strictly increasing list of native x-fast indices.
/// Their source order never changes. Positions use native cell centers and
/// physical volumes; periodic topology does not synthesize periodic images.
/// The numerical budget covers indices, positions and both six-stress buffers
/// (nominally 128 bytes/cell on 64-bit hosts). The cloned immutable EOS metadata
/// is separate, as in the original solver's numerical-workspace accounting.
/// Sampling allocates no adapter buffers and does not mutate the native flow.
#[derive(Debug)]
pub struct PreparedCompressibleVolume3 {
    source: Source,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    config: CompressibleVolumeConfig3,
    cells: Vec<usize>,
    positions: Vec<[f64; 3]>,
    moments: Vec<[f64; 6]>,
    trial: Vec<[f64; 6]>,
    time_s: Option<f64>,
    with_viscous_stress: bool,
}

impl PreparedCompressibleVolume3 {
    /// Checked numerical payload, excluding the original immutable EOS clone.
    pub fn required_workspace_bytes(cells: usize) -> Result<usize, CompressibleVolumeError3> {
        cells
            .checked_mul(size_of::<usize>() + size_of::<[f64; 3]>() + 2 * size_of::<[f64; 6]>())
            .ok_or(CompressibleVolumeError3::Invalid(
                "source workspace overflow",
            ))
    }

    /// Bind the full frozen source and publish its current accepted sample.
    pub fn from_frozen(
        flow: &Euler3d,
        cells: &[usize],
        config: CompressibleVolumeConfig3,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        Self::admit(
            flow.grid(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let mut result = Self::prepare(
            Source::Frozen(flow.gas().clone()),
            *flow.grid(),
            flow.boundaries(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        result.sample_frozen(flow, &mut cancelled)?;
        Ok(result)
    }

    /// Bind every original species, caloric reference and source inversion control.
    /// Initial extraction uses the actual cell composition, including absent species.
    pub fn from_mixture(
        flow: &MixtureEuler3d,
        cells: &[usize],
        config: CompressibleVolumeConfig3,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        Self::admit(
            flow.grid(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let mut result = Self::prepare(
            Source::Mixture(flow.gas().clone()),
            *flow.grid(),
            flow.boundaries(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        result.sample_mixture(flow, &mut cancelled)?;
        Ok(result)
    }

    fn admit(
        grid: &CartesianGrid,
        cells: &[usize],
        config: CompressibleVolumeConfig3,
        maximum: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(cancelled)?;
        if cells.is_empty() {
            return Err(CompressibleVolumeError3::Invalid(
                "nonempty increasing native cell indices required",
            ));
        }
        let required = Self::required_workspace_bytes(cells.len())?;
        admit_bytes(required, maximum)?;
        for (i, &cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            if grid.cell_center(cell).is_none() || i > 0 && cells[i - 1] >= cell {
                return Err(CompressibleVolumeError3::Invalid(
                    "nonempty increasing native cell indices required",
                ));
            }
        }
        if !config.reference_pressure_pa.is_finite()
            || config.reference_pressure_pa <= 0.0
            || config.origin_m.iter().any(|x| !x.is_finite())
            || !kernel::proper_rotation(config.grid_to_acoustic)
        {
            return Err(CompressibleVolumeError3::Invalid(
                "finite origin, proper stationary rotation and positive reference pressure required",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        source: Source,
        grid: CartesianGrid,
        boundaries: Boundaries3d,
        cells: &[usize],
        config: CompressibleVolumeConfig3,
        maximum: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        let mut result = Self {
            source,
            grid,
            boundaries,
            config,
            cells: allocated(cells.len(), 0)?,
            positions: allocated(cells.len(), [0.0; 3])?,
            moments: allocated(cells.len(), [0.0; 6])?,
            trial: allocated(cells.len(), [0.0; 6])?,
            time_s: None,
            with_viscous_stress: false,
        };
        admit_bytes(result.workspace_bytes(), maximum)?;
        result.cells.copy_from_slice(cells);
        for (i, &cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let local = grid
                .cell_center(cell)
                .ok_or(CompressibleVolumeError3::GeometryMismatch)?;
            let rotated = kernel::rotate(config.grid_to_acoustic, local);
            let position = std::array::from_fn(|a| config.origin_m[a] + rotated[a]);
            if position.iter().any(|x| !x.is_finite()) {
                return Err(CompressibleVolumeError3::Invalid(
                    "unrepresentable source position",
                ));
            }
            result.positions[i] = position;
        }
        poll(cancelled)?;
        Ok(result)
    }

    /// Original selected native cell indices in the fixed observation order.
    pub fn cell_indices(&self) -> &[usize] {
        &self.cells
    }
    /// Stationary source positions, directly usable by the quadrupole observer.
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions
    }
    /// Accepted instantaneous xx, yy, zz, xy, xz, yz moments [N m].
    /// They are histories to transform, not already frequency-domain phasors.
    pub fn stress_moments_n_m(&self) -> &[[f64; 6]] {
        &self.moments
    }
    /// Actual accepted native flow time; no caller clock or timestep is substituted.
    pub fn sample_time_s(&self) -> f64 {
        self.time_s.expect("constructors publish an initial sample")
    }
    /// Record this already accepted native sample with its original time/order.
    /// The joint recorder admits a uniform cadence; no resampling or hidden
    /// mean subtraction is performed by this adapter.
    pub fn record_sample(
        &self,
        recorder: &mut super::recorded::QuadrupoleWelchStream,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), crate::AeroacError> {
        recorder.push_sample(
            &self.positions,
            self.sample_time_s(),
            &self.moments,
            cancelled,
        )
    }
    pub const fn config(&self) -> CompressibleVolumeConfig3 {
        self.config
    }
    /// Whether the immutable source selection includes actual configured
    /// viscous/SGS stress. The sampling method must preserve this choice.
    pub const fn includes_viscous_stress(&self) -> bool {
        self.with_viscous_stress
    }
    pub fn workspace_bytes(&self) -> usize {
        self.cells.capacity() * size_of::<usize>()
            + self.positions.capacity() * size_of::<[f64; 3]>()
            + (self.moments.capacity() + self.trial.capacity()) * size_of::<[f64; 6]>()
    }
    /// Full bound original source, not a reduced density/gamma surrogate.
    pub fn frozen_source(&self) -> Option<&FrozenGas> {
        match &self.source {
            Source::Frozen(gas) => Some(gas),
            Source::Mixture(_) => None,
        }
    }
    /// Full bound original species axis and caloric references.
    pub fn mixture_source(&self) -> Option<&MixtureGas> {
        match &self.source {
            Source::Mixture(gas) => Some(gas),
            Source::Frozen(_) => None,
        }
    }

    /// Atomically replace the source sample after a later accepted frozen-gas step.
    pub fn sample_frozen(
        &mut self,
        flow: &Euler3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if self.with_viscous_stress || self.frozen_source() != Some(flow.gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(flow.grid(), flow.boundaries(), flow.time().value())?;
        self.sample_cells(flow.time().value(), &mut cancelled, |cell, cancel| {
            flow.states()[cell]
                .primitive(flow.gas(), cancel)
                .map_err(|error| {
                    if error == fs_flux::compressible::Error::Cancelled {
                        CompressibleVolumeError3::Cancelled
                    } else {
                        CompressibleVolumeError3::Frozen(error)
                    }
                })
        })
    }

    /// Atomically replace the source sample after a later accepted mixture step.
    pub fn sample_mixture(
        &mut self,
        flow: &MixtureEuler3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if self.with_viscous_stress || self.mixture_source() != Some(flow.gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(flow.grid(), flow.boundaries(), flow.time().value())?;
        self.sample_cells(flow.time().value(), &mut cancelled, |cell, cancel| {
            flow.cell(cell)
                .ok_or(CompressibleVolumeError3::GeometryMismatch)?
                .primitive(flow.gas(), cancel)
                .map_err(|error| {
                    if error == fs_flux::compressible::spatial::mixture::Error::Cancelled {
                        CompressibleVolumeError3::Cancelled
                    } else {
                        CompressibleVolumeError3::Mixture(error)
                    }
                })
        })
    }

    fn admit_sample(
        &self,
        grid: &CartesianGrid,
        boundaries: Boundaries3d,
        time: f64,
    ) -> Result<(), CompressibleVolumeError3> {
        if *grid != self.grid || boundaries != self.boundaries {
            return Err(CompressibleVolumeError3::GeometryMismatch);
        }
        if !time.is_finite() || time < 0.0 {
            return Err(CompressibleVolumeError3::Invalid(
                "invalid native accepted time",
            ));
        }
        if let Some(previous_s) = self.time_s
            && time <= previous_s
        {
            return Err(CompressibleVolumeError3::SampleClock {
                previous_s,
                candidate_s: time,
            });
        }
        Ok(())
    }

    fn sample_cells<C: FnMut() -> bool>(
        &mut self,
        time: f64,
        cancelled: &mut C,
        mut cell_primitive: impl FnMut(usize, &mut C) -> Result<Primitive, CompressibleVolumeError3>,
    ) -> Result<(), CompressibleVolumeError3> {
        for (i, &cell) in self.cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let p = cell_primitive(cell, cancelled)?;
            self.trial[i] = kernel::inviscid_moment(
                self.config.acoustic.density_kg_m3(),
                self.config.acoustic.sound_speed_m_s(),
                self.config.reference_pressure_pa,
                p.density,
                p.pressure,
                self.grid.cell_volume(),
                kernel::rotate(self.config.grid_to_acoustic, p.velocity),
            )
            .ok_or(CompressibleVolumeError3::Invalid(
                "unrepresentable inviscid stress moment",
            ))?;
        }
        poll(cancelled)?;
        std::mem::swap(&mut self.moments, &mut self.trial);
        self.time_s = Some(time);
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CompressibleVolumeError3> {
    if cancelled() {
        Err(CompressibleVolumeError3::Cancelled)
    } else {
        Ok(())
    }
}
fn allocated<T: Clone>(n: usize, value: T) -> Result<Vec<T>, CompressibleVolumeError3> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(n)
        .map_err(|_| CompressibleVolumeError3::Allocation)?;
    result.resize(n, value);
    Ok(result)
}
fn admit_bytes(required: usize, maximum: usize) -> Result<(), CompressibleVolumeError3> {
    if required > maximum {
        Err(CompressibleVolumeError3::Workspace {
            required_bytes: required,
            maximum_bytes: maximum,
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
