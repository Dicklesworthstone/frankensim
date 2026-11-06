//! Exact complete modal motion with source-derived membrane strain-rate loss.
//!
//! The original mass-normalized membrane basis fixes inertia, stiffness and
//! geometry. The full material matrix C is retained in q'' + C q' + Lambda q = f.
//! The existing port-Hamiltonian exponential owner prepares the fixed-clock
//! flow and independently integrated positive heat. No modal diagonalization,
//! fitted decay rate, extra prestress loss or algorithmic damping is introduced.
//!
//! State order is [q_0, ..., q_(n-1), qdot_0, ..., qdot_(n-1)]. The n inputs
//! are held generalized forces; the n outputs are interval-mean velocities.
//! Their product is physical work. Motion is linearized at the retained shape;
//! a successful step does not certify a finite-amplitude trajectory, a measured
//! material model, or that a forcing history stays in its declared source band.

use super::{MembraneModalViscosity, MembraneViscosity, MembraneViscosityError};
use crate::membrane::{MembraneError, poll};
use crate::membrane::vibration::{MembraneModes, MembraneVibration};
use fs_phs::{PhsError, PortHamiltonian, PreparedStepRecord, QuadraticStorage, ResistanceWitness};
use fs_phs::linear_exact::{LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, PreparedExactZoh};

/// Complete membrane-mode ceiling for the existing exact owner with one force
/// port per mode: its augmented extent is 2*(2*n + n + 1) <= 256.
pub const MAX_PREPARED_MEMBRANE_MODES: usize = (MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2 - 1) / 3;

/// Original source/geometry refusal or exact propagator refusal.
#[derive(Debug)]
pub enum MembraneDampingTimeError {
    Membrane(MembraneError),
    Viscosity(MembraneViscosityError),
    Propagation(PhsError),
}
impl From<MembraneError> for MembraneDampingTimeError {
    fn from(error: MembraneError) -> Self { Self::Membrane(error) }
}
impl From<MembraneViscosityError> for MembraneDampingTimeError {
    fn from(error: MembraneViscosityError) -> Self { Self::Viscosity(error) }
}
impl From<PhsError> for MembraneDampingTimeError {
    fn from(error: PhsError) -> Self { Self::Propagation(error) }
}
impl core::fmt::Display for MembraneDampingTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Membrane(error) => write!(f, "membrane damping: {error}"),
            Self::Viscosity(error) => write!(f, "membrane viscosity: {error}"),
            Self::Propagation(error) => write!(f, "membrane damping propagation: {error}"),
        }
    }
}
impl std::error::Error for MembraneDampingTimeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Membrane(error) => Some(error),
            Self::Viscosity(error) => Some(error),
            Self::Propagation(error) => Some(error),
        }
    }
}

/// Complete source-bound modal dynamics, prepared once for a fixed duration.
///
/// Source references remain immutable. Callers own the resumable physical state;
/// rejected/cancelled calls leave both caller output buffers unchanged. Successful
/// warm `step_into` calls allocate no owner buffers and compute no exponentials.
pub struct PreparedMembraneDamping<'model, 'patch> {
    viscosity: &'model MembraneViscosity<'patch>,
    modes: &'model MembraneModes<'patch>,
    projected: MembraneModalViscosity,
    flow: PreparedExactZoh<'static>,
    next: Vec<f64>,
    output: Vec<f64>,
}

impl<'model, 'patch> PreparedMembraneDamping<'model, 'patch> {
    /// Prepare all modes of the supplied basis at its retained linearization.
    ///
    /// `budget.max_modes` bounds scalar state coordinates (TWO per membrane
    /// mode), not a truncation count. The complete force-port map also needs
    /// n <= `MAX_PREPARED_MEMBRANE_MODES` (currently 42); no higher modes are
    /// discarded to fit that capacity. `max_projection_terms` bounds the
    /// source-row projection before its allocations. Cancellation is checked
    /// around the original synchronous dense preparation leaves.
    pub fn new(
        viscosity: &'model MembraneViscosity<'patch>,
        modes: &'model MembraneModes<'patch>,
        duration_s: f64,
        budget: LinearSplitBudget,
        max_projection_terms: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneDampingTimeError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(MembraneError::Invalid("positive finite membrane damping duration required").into());
        }
        let count = modes.modes().len();
        let dimension = count.checked_mul(2)
            .ok_or(MembraneError::Budget("membrane damping state extent"))?;
        if count == 0 || count > MAX_PREPARED_MEMBRANE_MODES || dimension > budget.max_modes {
            return Err(MembraneError::Budget("complete membrane damping state exceeds exact-flow budget").into());
        }
        let projected = viscosity.project_modes(modes, max_projection_terms, &mut cancelled)?;
        let square = dimension.checked_mul(dimension)
            .ok_or(MembraneError::Budget("membrane damping matrix extent"))?;
        let mut q = vec![0.0; square];
        let mut j = vec![0.0; square];
        let mut r = vec![0.0; square];
        let mut g = vec![0.0; dimension * count];
        for (i, mode) in modes.modes().iter().enumerate() {
            poll(&mut cancelled)?;
            q[i * dimension + i] = mode.lambda;
            q[(count + i) * dimension + count + i] = 1.0;
            j[i * dimension + count + i] = 1.0;
            j[(count + i) * dimension + i] = -1.0;
            g[(count + i) * count + i] = 1.0;
            for k in 0..count {
                r[(count + i) * dimension + count + k] =
                    projected.matrix_per_s()[i * count + k];
            }
        }
        // Embed the original physical strain-rate squares into velocity effort
        // coordinates. This preserves directional loss and exact null channels
        // through the existing positive-heat owner without an eigenvalue floor.
        let source_rows = projected.dissipation_rows();
        let row_count = source_rows.len() / count;
        let row_extent = row_count.checked_mul(dimension)
            .ok_or(MembraneError::Budget("membrane damping source rows"))?;
        let mut rows = vec![0.0; row_extent];
        for (row, original) in source_rows.chunks_exact(count).enumerate() {
            poll(&mut cancelled)?;
            rows[row * dimension + count..(row + 1) * dimension].copy_from_slice(original);
        }
        let witness = ResistanceWitness::try_from_source_rows(dimension, r.clone(), rows)?;
        poll(&mut cancelled)?;
        let storage = QuadraticStorage::new(q, dimension)?;
        let system = PortHamiltonian::new(dimension, count, j, r, g, Box::new(storage))?
            .with_resistance_witness(witness)?;
        let flow = PreparedExactZoh::new_owned(system, duration_s, budget)?;
        poll(&mut cancelled)?;
        Ok(Self { viscosity, modes, projected, flow,
            next: vec![0.0; dimension], output: vec![0.0; count] })
    }

    /// Bind the SAME accepted nonlinear equilibrium and complete vibration basis.
    pub fn from_vibration(
        viscosity: &'model MembraneViscosity<'patch>,
        vibration: &'model MembraneVibration<'patch>,
        duration_s: f64,
        budget: LinearSplitBudget,
        max_projection_terms: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneDampingTimeError> {
        Self::new(viscosity, vibration.basis(), duration_s, budget, max_projection_terms, cancelled)
    }

    pub fn viscosity(&self) -> &MembraneViscosity<'patch> { self.viscosity }
    pub fn modes(&self) -> &MembraneModes<'patch> { self.modes }
    pub fn projected_viscosity(&self) -> &MembraneModalViscosity { &self.projected }
    pub fn duration_s(&self) -> f64 { self.flow.duration() }
    pub fn mode_count(&self) -> usize { self.output.len() }
    pub fn state_dimension(&self) -> usize { self.next.len() }
    /// Original complete quadratic owner, available for harmonic/port analysis.
    pub fn system(&self) -> &PortHamiltonian { self.flow.system() }

    /// Exact fixed-duration, held-force transition with independent physical heat.
    ///
    /// `record.dissipated` is integrated material heat [J], `record.supplied`
    /// is held-force work [J], and `record.delta_h` is mechanical energy change.
    /// The shared owner enforces their balance; none is numerical damping.
    /// Input forces use the original basis's `project_nodal_force` convention.
    /// Output velocities use its same `reconstruct` convention.
    pub fn step_into(
        &mut self,
        initial: &[f64],
        generalized_force: &[f64],
        endpoint: &mut [f64],
        mean_modal_velocity: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedStepRecord, MembraneDampingTimeError> {
        poll(&mut cancelled)?;
        if endpoint.len() != self.next.len() || mean_modal_velocity.len() != self.output.len() {
            return Err(MembraneError::Invalid("complete membrane damping output buffers required").into());
        }
        let record = self.flow.step_into(initial, generalized_force, &mut self.next, &mut self.output)?;
        poll(&mut cancelled)?;
        endpoint.copy_from_slice(&self.next);
        mean_modal_velocity.copy_from_slice(&self.output);
        Ok(record)
    }
}
