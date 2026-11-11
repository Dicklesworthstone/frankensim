//! Exact held-input evolution of the complete acoustic and causal wall fields.
//!
//! The original consistent volume and surface masses remain in the physical
//! equations. Each wall retains every displacement, velocity and Maxwell memory
//! coordinate, including trace coordinates adjacent to a pressure-release face.
//! The internal invertible change from memory z to extension e = sqrt(k) w - z
//! makes the storage block diagonal. It removes no state or relaxation time.
//!
//! In canonical fluid and wall momenta, the reciprocal coupling is -W v in the
//! fluid and +W p in the wall. The existing pHS exponential/Gramian owner advances
//! this ONE complete operator. Independent same-trajectory observations integrate
//! real Robin heat, wall resistance heat, Maxwell heat and signed fluid-to-wall
//! work. Neither an endpoint energy deficit nor an effective wall admittance is
//! reported as dissipation.
//!
//! Spatial causal-admittance patches may share this same acoustic trajectory.
//! Their original trace pressures obey tau q_dot + q = p, including zero-residue
//! and pressure-release-shared coordinates. They retain their original storage,
//! heat and reversible transfer independently of the mechanical wall accounts.
//! Prescribed pressure tractions may independently drive the original wall
//! coordinates. Their input columns and velocity observations use the SAME
//! physical surface mass, and their signed external work is separate from
//! fluid-to-wall exchange. Warm steps do not reassemble a traction profile.
//!
//! Preparation is dense and explicitly bounded; warm steps reuse cached matrix
//! contractions. Exactness concerns the semidiscrete held-input solution up to
//! numerical exponential error. Spatial convergence, input bandwidth and output
//! sampling remain separate. Wall coefficients are constant on each patch, with
//! the original locally reacting law; this does not add shell bending or infer a
//! causal realization from harmonic admittance samples.

use super::super::super::SpatialRelaxingBoundary;
use super::super::{
    PressureReleaseConstraints,
    exact::{ExactTimeLimits, ExactVelocityPort},
};
use super::relaxing::{self, TraceMemory};
use super::traction::WallTractionProfile;
use super::{
    BoundaryAdmittance, BoundaryVelocitySource, C64, Coo, Csr, Dd, HeterogeneousError,
    HeterogeneousHelmholtzFe, HeterogeneousTimeEnergy, HeterogeneousTimeState, MechanicalBoundary,
    MechanicalTimeEnergy, RelaxationBranch, ResistiveBoundary, SpatialResistiveBoundary,
    VolumeVelocitySource, energy_change, invalid, poll, positive_form,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use fs_la::factor::cholesky;
use fs_math::det;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    prepared::PreparedQuadraticFlow, split_quadratic,
};
use fs_phs::{PortHamiltonian, PreparedStepError, QuadraticStorage};

mod admittance;
use admittance::ExactAdmittance;

/// Admission counts EVERY free fluid and original wall coordinate, not modes.
#[derive(Clone, Copy, Debug)]
pub struct ExactMechanicalTimeLimits {
    /// Two states per free pressure DOF, two per wall trace DOF, and one per
    /// trace DOF per Maxwell arm and spatial-admittance arm. Original sparse forms and full lifted acoustic
    /// vectors are additional. The native augmented-state cap also counts inputs.
    pub max_states: usize,
    /// Native dense cache, byte and multiplication-work admission. Coordinate
    /// maps and physical forms additionally occupy O(max_states^2) storage.
    pub preparation: DenseExponentialBudget,
    /// Relative physical work defect, plus independently evaluated roundoff.
    pub energy_balance_tolerance: f64,
}

impl Default for ExactMechanicalTimeLimits {
    fn default() -> Self {
        let acoustic = ExactTimeLimits::default();
        Self {
            max_states: 2 * acoustic.max_dofs,
            preparation: acoustic.preparation,
            energy_balance_tolerance: acoustic.energy_balance_tolerance,
        }
    }
}

/// A complete original physical state, bound to this wall configuration.
#[derive(Clone)]
pub struct ExactMechanicalTimeState<'p, 'm, 'f, 'c> {
    owner: &'p PreparedExactMechanicalTime<'m, 'f, 'c>,
    acoustic: HeterogeneousTimeState<'m, 'f, 'c>,
    displacements: Vec<Vec<f64>>,
    velocities: Vec<Vec<f64>>,
    memories: Vec<Vec<Vec<f64>>>,
    admittance_memories: Vec<Vec<f64>>,
}

impl core::fmt::Debug for ExactMechanicalTimeState<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ExactMechanicalTimeState")
            .field("acoustic", &self.acoustic)
            .field("wall_patches", &self.displacements.len())
            .finish_non_exhaustive()
    }
}

impl<'m, 'f, 'c> ExactMechanicalTimeState<'_, 'm, 'f, 'c> {
    /// Original acoustic field, directly usable by the existing spatial probes.
    #[must_use]
    pub fn acoustic_state(&self) -> &HeterogeneousTimeState<'m, 'f, 'c> {
        &self.acoustic
    }
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        self.acoustic.pressure_potential_pa_s()
    }
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        self.acoustic.pressure_pa()
    }
    /// Displacement [m], in this patch's original sorted wall_dofs order.
    #[must_use]
    pub fn wall_displacement_m(&self, wall: usize) -> Option<&[f64]> {
        self.displacements.get(wall).map(Vec::as_slice)
    }
    /// Outward normal velocity [m/s], in the same original trace order.
    #[must_use]
    pub fn wall_normal_velocity_m_s(&self, wall: usize) -> Option<&[f64]> {
        self.velocities.get(wall).map(Vec::as_slice)
    }
    /// Original energy-normalized memory z, NOT the internal extension e.
    #[must_use]
    pub fn wall_relaxation_memory(&self, wall: usize, arm: usize) -> Option<&[f64]> {
        self.memories.get(wall)?.get(arm).map(Vec::as_slice)
    }
    /// Original spatial-admittance memory pressure q [Pa], in prepared trace order.
    #[must_use]
    pub fn admittance_memory_pressure_pa(&self, arm: usize) -> Option<&[f64]> {
        self.admittance_memories.get(arm).map(Vec::as_slice)
    }
}

/// Original storage changes, independently integrated heat and reversible work.
#[derive(Clone, Copy, Debug)]
pub struct ExactMechanicalTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    pub acoustic_delta_energy_j: f64,
    pub wall_delta_energy_j: f64,
    pub admittance_delta_energy_j: f64,
    /// All acoustic and applied wall-traction work, counted once.
    pub input_work_j: f64,
    /// Original fluid-source work, summed independently of mechanical forcing.
    /// This retains weak acoustic work without subtracting large external works.
    pub acoustic_input_work_j: f64,
    /// Prepared boundary-source subtotal, already included in input_work_j.
    pub boundary_source_work_j: f64,
    /// Volume-source subtotal, already included in input_work_j.
    pub volume_source_work_j: f64,
    /// Applied pressure-traction work on the wall, already included in input_work_j.
    /// This is external mechanical input, not fluid source work or wall heat.
    pub wall_traction_work_j: f64,
    pub boundary_dissipation_j: f64,
    pub wall_resistance_dissipation_j: f64,
    pub wall_relaxation_dissipation_j: f64,
    /// Signed integral p^T W v dt; negative means the wall returns fluid energy.
    pub wall_work_j: f64,
    /// Independently integrated sum q_j^T Cj q_j [J].
    pub admittance_dissipation_j: f64,
    /// Signed integral sum p^T Cj q_j [J], from fluid to admittance memories.
    pub admittance_work_j: f64,
}

impl ExactMechanicalTimeWork {
    #[must_use]
    pub fn wall_dissipation_j(self) -> f64 {
        self.wall_resistance_dissipation_j + self.wall_relaxation_dissipation_j
    }
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j
            + self.boundary_dissipation_j
            + self.wall_dissipation_j()
            + self.admittance_dissipation_j
            - self.input_work_j
    }
    #[must_use]
    pub fn acoustic_balance_residual_j(self) -> f64 {
        self.acoustic_delta_energy_j
            + self.boundary_dissipation_j
            + self.wall_work_j
            + self.admittance_work_j
            - self.acoustic_input_work_j
    }
    #[must_use]
    pub fn wall_balance_residual_j(self) -> f64 {
        self.wall_delta_energy_j + self.wall_dissipation_j()
            - self.wall_work_j
            - self.wall_traction_work_j
    }
    #[must_use]
    pub fn admittance_balance_residual_j(self) -> f64 {
        self.admittance_delta_energy_j + self.admittance_dissipation_j - self.admittance_work_j
    }
}

/// Checked candidate; a refusal never mutates the accepted state or clock.
#[derive(Clone, Debug)]
pub struct ExactMechanicalTimeStep<'p, 'm, 'f, 'c> {
    pub state: ExactMechanicalTimeState<'p, 'm, 'f, 'c>,
    pub work: ExactMechanicalTimeWork,
    /// Integral of each boundary shape against pressure [Pa m^2 s], with linear
    /// ports first and prepared boundary sources second, in preparation order.
    pub integrated_port_pressure_pa_m2_s: Vec<f64>,
    /// Integral of each original volume source against pressure [Pa s].
    pub integrated_volume_source_pressure_pa_s: Vec<f64>,
    /// Integral of each dimensionless traction shape against outward wall
    /// velocity over its original surface and this interval [m^3]. Multiplying
    /// by its held traction amplitude [Pa] gives signed external work [J].
    pub integrated_wall_traction_displacement_m3: Vec<f64>,
}

pub(in crate::acoustics::helmholtz::heterogeneous::time) struct MassMap {
    pub(in crate::acoustics::helmholtz::heterogeneous::time) scales: Vec<f64>,
    pub(in crate::acoustics::helmholtz::heterogeneous::time) normalized: Vec<f64>,
    pub(in crate::acoustics::helmholtz::heterogeneous::time) inverse: Vec<f64>,
}

struct ExactWall {
    boundary: MechanicalBoundary,
    dofs: Vec<usize>,
    mass: Csr,
    map: MassMap,
    offset: usize,
    inertia_scale: f64,
    stiffness_scale: f64,
    couplings: Vec<f64>,
}

/// Immutable full spatial operator, physical coordinate maps and trajectory cache.
pub struct PreparedExactMechanicalTime<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    duration_s: f64,
    pressure_release: PressureReleaseConstraints,
    potential_scale: Vec<f64>,
    fluid_mass: MassMap,
    walls: Vec<ExactWall>,
    admittance: Vec<ExactAdmittance>,
    states: usize,
    linear_port_count: usize,
    boundary_port_count: usize,
    volume_source_count: usize,
    wall_traction_count: usize,
    flow: PreparedQuadraticFlow,
    limits: ExactMechanicalTimeLimits,
}

fn propagation(error: PreparedStepError) -> HeterogeneousError {
    match error {
        PreparedStepError::Cancelled => HeterogeneousError::Cancelled,
        PreparedStepError::Solver(error) => HeterogeneousError::Propagation(error),
    }
}

fn admit_extent(
    states: usize,
    sources: usize,
    limits: ExactMechanicalTimeLimits,
) -> Result<(), HeterogeneousError> {
    if states > limits.max_states {
        return Err(HeterogeneousError::Limit {
            what: "complete exact acoustic and mechanical wall states",
            requested: states,
            maximum: limits.max_states,
        });
    }
    let augmented = states
        .checked_add(sources)
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| invalid("complete exact wall source extent overflow"))?;
    let maximum = limits
        .preparation
        .max_dimension
        .min(MAX_LINEAR_EXPONENTIAL_BLOCK_DIM);
    if augmented > maximum {
        return Err(HeterogeneousError::Limit {
            what: "complete exact acoustic and wall augmented states",
            requested: augmented,
            maximum,
        });
    }
    Ok(())
}

pub(in crate::acoustics::helmholtz::heterogeneous::time) fn diagonal_scales(
    matrix: &Csr,
    dofs: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, HeterogeneousError> {
    let mut scales = Vec::with_capacity(dofs.len());
    for &row in dofs {
        poll(cancelled)?;
        let (columns, values) = matrix.row(row);
        let diagonal = columns
            .binary_search(&row)
            .ok()
            .map(|i| values[i])
            .unwrap_or(0.0);
        let scale = det::sqrt(diagonal);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(
                invalid("positive original fluid and wall storage diagonals required").into(),
            );
        }
        scales.push(scale);
    }
    Ok(scales)
}

impl MassMap {
    // Restrict the ORIGINAL consistent mass BEFORE inversion. An inverse
    // restricted afterward would give different physical inertia.
    pub(in crate::acoustics::helmholtz::heterogeneous::time) fn new(
        matrix: &Csr,
        dofs: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, HeterogeneousError> {
        let n = dofs.len();
        let scales = diagonal_scales(matrix, dofs, cancelled)?;
        let mut normalized = vec![0.0; n * n];
        for (row, &global) in dofs.iter().enumerate() {
            poll(cancelled)?;
            let (columns, values) = matrix.row(global);
            for (&column, &value) in columns.iter().zip(values) {
                if let Ok(col) = dofs.binary_search(&column) {
                    normalized[row * n + col] = value / scales[row] / scales[col];
                }
            }
        }
        let factor = cholesky(&normalized, n).map_err(HeterogeneousError::MassFactor)?;
        let mut inverse = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for col in 0..n {
            poll(cancelled)?;
            rhs.fill(0.0);
            rhs[col] = 1.0;
            factor.solve(&mut rhs);
            for row in 0..n {
                inverse[row * n + col] = rhs[row];
            }
        }
        // Check every original entry, including the upper triangle. No
        // symmetric projection or diagonal mass replaces the original form.
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..n {
                let mut value = Dd::from_f64(0.0);
                let mut bound = 0.0;
                for k in 0..n {
                    let term =
                        Dd::from_f64(normalized[row * n + k]) * Dd::from_f64(inverse[k * n + col]);
                    value = value + term;
                    bound += term.to_f64().abs();
                }
                let expected: f64 = if row == col { 1.0 } else { 0.0 };
                if !bound.is_finite()
                    || !value.to_f64().is_finite()
                    || (value.to_f64() - expected).abs()
                        > 512.0 * f64::EPSILON * n as f64 * bound.max(expected)
                {
                    return Err(
                        invalid("original fluid or wall mass inverse failed its residual").into(),
                    );
                }
            }
        }
        Ok(Self {
            scales,
            normalized,
            inverse,
        })
    }

    // Original physical field from scaled canonical momentum.
    pub(in crate::acoustics::helmholtz::heterogeneous::time) fn physical_map(
        &self,
        inertia_scale: f64,
    ) -> Vec<f64> {
        let n = self.scales.len();
        (0..n * n)
            .map(|i| self.inverse[i] / inertia_scale / self.scales[i / n])
            .collect()
    }
}

// Original physical form P^T W P, retaining all cross-coordinate terms.
pub(in crate::acoustics::helmholtz::heterogeneous::time) fn pullback(
    matrix: &Csr,
    map: &[f64],
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, HeterogeneousError> {
    let n = matrix.nrows();
    let mut right = vec![0.0; n * n];
    for row in 0..n {
        poll(cancelled)?;
        let (columns, values) = matrix.row(row);
        for col in 0..n {
            let mut value = Dd::from_f64(0.0);
            for (&k, &entry) in columns.iter().zip(values) {
                value = value + Dd::from_f64(entry) * Dd::from_f64(map[k * n + col]);
            }
            right[row * n + col] = value.to_f64();
        }
    }
    let mut form = vec![0.0; n * n];
    for row in 0..n {
        poll(cancelled)?;
        for col in 0..n {
            let mut value = Dd::from_f64(0.0);
            for k in 0..n {
                value = value + Dd::from_f64(map[k * n + row]) * Dd::from_f64(right[k * n + col]);
            }
            form[row * n + col] = scale * value.to_f64();
        }
    }
    Ok(form)
}

pub(in crate::acoustics::helmholtz::heterogeneous::time) fn put_block(
    target: &mut [f64],
    states: usize,
    offset: usize,
    block: &[f64],
    n: usize,
) {
    for row in 0..n {
        target[(offset + row) * states + offset..(offset + row) * states + offset + n]
            .copy_from_slice(&block[row * n..(row + 1) * n]);
    }
}

pub(in crate::acoustics::helmholtz::heterogeneous::time) fn add_skew(
    matrix: &mut [f64],
    states: usize,
    row: usize,
    col: usize,
    value: f64,
) {
    matrix[row * states + col] += value;
    matrix[col * states + row] -= value;
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare one exact held-input trajectory for the complete acoustic field
    /// and the original causal mechanical/Maxwell wall laws.
    ///
    /// Mechanical, real Robin and pressure-release faces must be disjoint.
    /// Adjacent wall patches retain independent complete traces. Driving sources
    /// may overlap mechanical faces, but never pressure-release faces. Every
    /// original wall memory must be supplied when importing a nonzero history.
    ///
    /// Boundary columns are linear ports followed by prepared boundary sources;
    /// their amplitudes are [m/s]. Original volume-source amplitudes are [m^3/s].
    /// All projections and source identity checks reuse their existing owners.
    /// Native dense admission refuses excessive complete state or source counts;
    /// there is no modal truncation or fallback to midpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_walls<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        ports: &[ExactVelocityPort],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactMechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_exact_time_with_walls_and_spatial_boundaries(
            duration_s,
            resistive,
            &[],
            boundaries,
            ports,
            boundary_sources,
            volume_sources,
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Exact full fluid/wall propagation with original spatial resistive traces.
    ///
    /// Pointwise memoryless conductance enters both the complete coupled
    /// dynamics and its trajectory-integrated boundary heat. Causal wall
    /// displacement, velocity, relaxation memories and their independent work
    /// remain unchanged. Spatial, uniform, mechanical and pressure-release
    /// patches must be disjoint. Each spatial law belongs to this exact model;
    /// no material sampling occurs during warm steps.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_walls_and_spatial_boundaries<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        spatial_boundaries: &[&SpatialResistiveBoundary<'_, '_, '_>],
        boundaries: &[MechanicalBoundary],
        ports: &[ExactVelocityPort],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactMechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_exact_time_with_walls_and_relaxing_boundaries(
            duration_s,
            resistive,
            spatial_boundaries,
            boundaries,
            &[],
            ports,
            boundary_sources,
            volume_sources,
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Evolve mechanical/Maxwell walls and causal spatial admittance in ONE
    /// complete held-source trajectory. All physical patches have disjoint faces;
    /// adjacent patches retain independent original trace states at shared DOFs.
    /// Every physical heat and signed exchange is integrated independently.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_walls_and_relaxing_boundaries<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        spatial_boundaries: &[&SpatialResistiveBoundary<'_, '_, '_>],
        boundaries: &[MechanicalBoundary],
        relaxing_boundaries: &[&SpatialRelaxingBoundary<'_, '_, '_>],
        ports: &[ExactVelocityPort],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactMechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_exact_time_with_wall_tractions(
            duration_s,
            resistive,
            spatial_boundaries,
            boundaries,
            relaxing_boundaries,
            ports,
            boundary_sources,
            volume_sources,
            &[],
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Prepare complete fluid/wall evolution with distributed applied tractions.
    ///
    /// Each dimensionless shape uses the selected wall's original GLOBAL
    /// hierarchical trace DOFs. Its held amplitude is pressure [Pa], positive
    /// outward in the same direction as fluid pressure on the wall. The original
    /// consistent surface mass projects both force and conjugate volume velocity;
    /// no prescribed motion, face averaging or alternate structural solve is used.
    /// Several coherent traction profiles may act on the same wall, including
    /// wall coordinates adjacent to a pressure-release face. Inputs follow all
    /// boundary and volume-source columns and count against dense admission.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_wall_tractions<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        spatial_boundaries: &[&SpatialResistiveBoundary<'_, '_, '_>],
        boundaries: &[MechanicalBoundary],
        relaxing_boundaries: &[&SpatialRelaxingBoundary<'_, '_, '_>],
        ports: &[ExactVelocityPort],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        wall_tractions: &[WallTractionProfile],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactMechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        super::super::super::boundary_source::admit_count(boundary_sources.len())?;
        super::super::super::source::admit_count(volume_sources.len())?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || limits.max_states == 0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0
            || limits.energy_balance_tolerance > 1e-6
        {
            return Err(invalid(
                "positive exact wall clock, full-state cap and bounded energy tolerance required",
            )
            .into());
        }
        let boundary_port_count = ports
            .len()
            .checked_add(boundary_sources.len())
            .ok_or_else(|| invalid("exact wall boundary source count overflow"))?;
        let source_count = boundary_port_count
            .checked_add(volume_sources.len())
            .and_then(|count| count.checked_add(wall_tractions.len()))
            .ok_or_else(|| invalid("exact wall source count overflow"))?;
        if pressure_release_faces.is_empty() {
            let states = self
                .ndof()
                .checked_mul(2)
                .ok_or_else(|| invalid("exact acoustic state extent overflow"))?;
            admit_extent(states, source_count, limits)?;
        }
        let pressure_release = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let n = pressure_release.free_dofs.len();
        let mut states = n
            .checked_mul(2)
            .ok_or_else(|| invalid("exact acoustic state extent overflow"))?;
        admit_extent(states, source_count, limits)?;
        let mut seen = pressure_release.faces.clone();
        for faces in resistive
            .iter()
            .map(|b| &b.faces)
            .chain(spatial_boundaries.iter().map(|b| &b.faces))
            .chain(boundaries.iter().map(|b| &b.faces))
            .chain(relaxing_boundaries.iter().map(|b| &b.instantaneous.faces))
        {
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid(
                        "mechanical, causal admittance, real Robin and pressure-release faces must be disjoint",
                    )
                    .into());
                }
            }
        }
        let physical_fields = spatial_boundaries
            .iter()
            .copied()
            .chain(relaxing_boundaries.iter().map(|b| &b.instantaneous))
            .collect::<Vec<_>>();
        let conductance =
            self.resistive_matrix(resistive, &physical_fields, &pressure_release, cancelled)?;
        let conductance = pressure_release.restrict(&conductance, cancelled)?;
        let mut walls = Vec::new();
        for boundary in boundaries {
            poll(cancelled)?;
            fs_phs::wall_specific_impedance(&boundary.wall, 1.0)
                .map_err(HeterogeneousError::Wall)?;
            if boundary.faces.is_empty() {
                return Err(
                    invalid("exact mechanical wall requires original boundary faces").into(),
                );
            }
            // The existing physical face owner retains curved surface metrics
            // and every original pressure trace basis function.
            let (global, _) = self.admittance_matrices(
                &[BoundaryAdmittance {
                    faces: boundary.faces.clone(),
                    admittance_m_pa_s: C64::ONE,
                }],
                cancelled,
            )?;
            let mut dofs = Vec::new();
            for row in 0..self.ndof() {
                poll(cancelled)?;
                let (columns, values) = global.row(row);
                if columns.binary_search(&row).is_ok_and(|i| values[i] > 0.0) {
                    dofs.push(row);
                }
            }
            if dofs.is_empty() {
                return Err(invalid("mechanical wall has no original trace coordinates").into());
            }
            let count = boundary
                .relaxation
                .len()
                .checked_add(2)
                .and_then(|count| count.checked_mul(dofs.len()))
                .ok_or_else(|| invalid("complete wall and memory state extent overflow"))?;
            let offset = states;
            states = states
                .checked_add(count)
                .ok_or_else(|| invalid("complete exact wall state extent overflow"))?;
            admit_extent(states, source_count, limits)?;
            let mut couplings = Vec::with_capacity(boundary.relaxation.len());
            for arm in &boundary.relaxation {
                poll(cancelled)?;
                couplings.push(
                    RelaxationBranch {
                        projection: vec![1.0],
                        stiffness: arm.stiffness_per_area,
                        relaxation_time_s: arm.relaxation_time_s,
                    }
                    .scalar_relaxed_memory(1.0)
                    .map_err(HeterogeneousError::Wall)?,
                );
            }
            let mut local = Coo::new(dofs.len(), dofs.len());
            for (row, &global_row) in dofs.iter().enumerate() {
                poll(cancelled)?;
                let (columns, values) = global.row(global_row);
                for (&col, &value) in columns.iter().zip(values) {
                    if value != 0.0 {
                        let col = dofs
                            .binary_search(&col)
                            .map_err(|_| invalid("incomplete original exact wall trace support"))?;
                        local.push(row, col, value);
                    }
                }
            }
            let mass = local.assemble();
            let map = MassMap::new(&mass, &(0..dofs.len()).collect::<Vec<_>>(), cancelled)?;
            walls.push(ExactWall {
                boundary: boundary.clone(),
                dofs,
                mass,
                map,
                offset,
                inertia_scale: det::sqrt(boundary.wall.surface_density),
                stiffness_scale: det::sqrt(boundary.wall.stiffness_per_area),
                couplings,
            });
        }
        let mut traction_loads = Vec::with_capacity(wall_tractions.len());
        for profile in wall_tractions {
            poll(cancelled)?;
            let wall = walls
                .get(profile.wall)
                .ok_or_else(|| invalid("exact traction references an absent mechanical wall"))?;
            let coefficients = profile.local_coefficients(&wall.dofs, cancelled)?;
            let mut load = Vec::with_capacity(wall.dofs.len());
            for row in 0..wall.dofs.len() {
                poll(cancelled)?;
                let (columns, values) = wall.mass.row(row);
                let mut force = Dd::from_f64(0.0);
                for (&col, &value) in columns.iter().zip(values) {
                    force = force + Dd::from_f64(value) * Dd::from_f64(coefficients[col]);
                }
                let force = force.to_f64();
                if !force.is_finite() {
                    return Err(invalid("original wall traction load is unrepresentable").into());
                }
                load.push(force);
            }
            traction_loads.push((profile.wall, load));
        }
        // The pHS backbone retains the original fluid and all mechanical states.
        // Added raw q memories are NEVER supplied a fictitious invertible residue.
        let base_states = states;
        let traces = relaxing::prepare(
            self,
            relaxing_boundaries,
            limits.max_states - base_states,
            cancelled,
        )?;
        let mut admittance = Vec::with_capacity(traces.len());
        for trace in traces {
            let offset = states;
            states = states
                .checked_add(trace.dofs.len())
                .ok_or_else(|| invalid("complete exact admittance state extent overflow"))?;
            admit_extent(states, source_count, limits)?;
            admittance.push(ExactAdmittance { trace, offset });
        }
        let total_states = states;
        let states = base_states;
        let mut loads = Vec::with_capacity(source_count);
        for port in ports {
            poll(cancelled)?;
            if port.profile.is_empty() {
                return Err(
                    invalid("exact acoustic source requires a nonempty trace profile").into(),
                );
            }
            for &(face, _) in &port.profile {
                let mut key = face;
                key.sort_unstable();
                if pressure_release.faces.contains(&key) {
                    return Err(invalid(
                        "pressure-release and exact velocity-port faces must not overlap",
                    )
                    .into());
                }
            }
            loads.push(self.geometry.boundary_linear_load(&port.profile)?);
        }
        for source in boundary_sources {
            poll(cancelled)?;
            source.validate(self)?;
            for face in &source.faces {
                if pressure_release.faces.contains(face) {
                    return Err(invalid(
                        "pressure-release and prepared boundary-source faces must not overlap",
                    )
                    .into());
                }
            }
            loads.push(source.load.clone());
        }
        for source in volume_sources {
            poll(cancelled)?;
            source.validate(self)?;
            loads.push(source.load.clone());
        }
        let fluid_mass = MassMap::new(&self.mass, &pressure_release.free_dofs, cancelled)?;
        let potential_scale =
            diagonal_scales(&self.stiffness, &pressure_release.free_dofs, cancelled)?;
        let pressure_map = fluid_mass.physical_map(1.0);
        let mut q = vec![0.0; states * states];
        let mut j = vec![0.0; states * states];
        let mut r = vec![0.0; states * states];
        let mut g = vec![0.0; states * source_count];
        for (row, &global_row) in pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let (columns, values) = self.stiffness.row(global_row);
            for (&global_col, &value) in columns.iter().zip(values) {
                let col = pressure_release.original_to_free[global_col];
                if col != usize::MAX {
                    q[row * states + col] = value / potential_scale[row] / potential_scale[col];
                }
            }
            add_skew(
                &mut j,
                states,
                row,
                n + row,
                potential_scale[row] / fluid_mass.scales[row],
            );
            let (columns, values) = conductance.row(row);
            for (&col, &value) in columns.iter().zip(values) {
                r[(n + row) * states + n + col] =
                    value / fluid_mass.scales[row] / fluid_mass.scales[col];
            }
            for (source, load) in loads.iter().enumerate() {
                g[(n + row) * source_count + source] = load[global_row] / fluid_mass.scales[row];
            }
        }
        let traction_offset = boundary_port_count + volume_sources.len();
        for (index, (wall_index, load)) in traction_loads.iter().enumerate() {
            let wall = &walls[*wall_index];
            let velocity_offset = wall.offset + wall.dofs.len();
            for (row, &force) in load.iter().enumerate() {
                poll(cancelled)?;
                g[(velocity_offset + row) * source_count + traction_offset + index] =
                    force / wall.inertia_scale / wall.map.scales[row];
            }
        }
        put_block(&mut q, states, n, &fluid_mass.inverse, n);
        // Four independent physical observations: Robin heat, wall resistance
        // heat, Maxwell heat and signed p^T W v exchange.
        let mut forms = vec![vec![0.0; states * states]; 4];
        put_block(
            &mut forms[0],
            states,
            n,
            &pullback(&conductance, &pressure_map, 1.0, cancelled)?,
            n,
        );
        for wall in &walls {
            poll(cancelled)?;
            let m = wall.dofs.len();
            let w0 = wall.offset;
            let v0 = w0 + m;
            let velocity_map = wall.map.physical_map(wall.inertia_scale);
            put_block(&mut q, states, w0, &wall.map.normalized, m);
            put_block(&mut q, states, v0, &wall.map.inverse, m);
            let resistance =
                wall.boundary.wall.resistance / wall.inertia_scale / wall.inertia_scale;
            for row in 0..m {
                poll(cancelled)?;
                add_skew(
                    &mut j,
                    states,
                    w0 + row,
                    v0 + row,
                    wall.stiffness_scale / wall.inertia_scale,
                );
                for col in 0..m {
                    r[(v0 + row) * states + v0 + col] =
                        resistance * wall.map.normalized[row * m + col];
                }
                let free = pressure_release.original_to_free[wall.dofs[row]];
                if free != usize::MAX {
                    let (columns, values) = wall.mass.row(row);
                    for (&col, &value) in columns.iter().zip(values) {
                        add_skew(
                            &mut j,
                            states,
                            n + free,
                            v0 + col,
                            -value
                                / fluid_mass.scales[free]
                                / wall.inertia_scale
                                / wall.map.scales[col],
                        );
                    }
                }
            }
            put_block(
                &mut forms[1],
                states,
                v0,
                &pullback(
                    &wall.mass,
                    &velocity_map,
                    wall.boundary.wall.resistance,
                    cancelled,
                )?,
                m,
            );
            for (arm, &coupling) in wall.couplings.iter().enumerate() {
                let e0 = w0 + (2 + arm) * m;
                let tau = wall.boundary.relaxation[arm].relaxation_time_s;
                put_block(&mut q, states, e0, &wall.map.normalized, m);
                for row in 0..m {
                    poll(cancelled)?;
                    add_skew(
                        &mut j,
                        states,
                        v0 + row,
                        e0 + row,
                        -coupling / wall.inertia_scale,
                    );
                    for col in 0..m {
                        r[(e0 + row) * states + e0 + col] = wall.map.inverse[row * m + col] / tau;
                        forms[2][(e0 + row) * states + e0 + col] =
                            wall.map.normalized[row * m + col] / tau;
                    }
                }
            }
            // Pull back the ORIGINAL p^T W v, independently of J and R.
            let mut trace_load = vec![0.0; n * m];
            for (row, &global_row) in wall.dofs.iter().enumerate() {
                poll(cancelled)?;
                let free = pressure_release.original_to_free[global_row];
                if free == usize::MAX {
                    continue;
                }
                let (columns, values) = wall.mass.row(row);
                for col in 0..m {
                    let mut value = Dd::from_f64(0.0);
                    for (&k, &entry) in columns.iter().zip(values) {
                        value =
                            value + Dd::from_f64(entry) * Dd::from_f64(velocity_map[k * m + col]);
                    }
                    trace_load[free * m + col] = value.to_f64();
                }
            }
            for row in 0..n {
                poll(cancelled)?;
                for col in 0..m {
                    let mut value = Dd::from_f64(0.0);
                    for k in 0..n {
                        value = value
                            + Dd::from_f64(pressure_map[k * n + row])
                                * Dd::from_f64(trace_load[k * m + col]);
                    }
                    let half = 0.5 * value.to_f64();
                    forms[3][(n + row) * states + v0 + col] += half;
                    forms[3][(v0 + col) * states + n + row] += half;
                }
            }
        }
        let mut observations = Vec::with_capacity(source_count);
        for load in &loads {
            poll(cancelled)?;
            let mut observation = vec![0.0; states];
            for col in 0..n {
                let mut value = Dd::from_f64(0.0);
                for row in 0..n {
                    value = value
                        + Dd::from_f64(load[pressure_release.free_dofs[row]])
                            * Dd::from_f64(pressure_map[row * n + col]);
                }
                observation[n + col] = value.to_f64();
            }
            observations.push(observation);
        }
        // Independently pull back the ORIGINAL shape^T W v into the native
        // velocity coordinates. Its integral is volume displacement, not pressure.
        for (wall_index, load) in &traction_loads {
            poll(cancelled)?;
            let wall = &walls[*wall_index];
            let m = wall.dofs.len();
            let velocity_map = wall.map.physical_map(wall.inertia_scale);
            let mut observation = vec![0.0; states];
            for col in 0..m {
                poll(cancelled)?;
                let mut value = Dd::from_f64(0.0);
                for row in 0..m {
                    value =
                        value + Dd::from_f64(load[row]) * Dd::from_f64(velocity_map[row * m + col]);
                }
                observation[wall.offset + m + col] = value.to_f64();
            }
            observations.push(observation);
        }
        let extension = if admittance.is_empty() {
            None
        } else {
            Some(admittance::extend(
                self,
                duration_s,
                &pressure_release,
                &fluid_mass,
                &pressure_map,
                &admittance,
                states,
                total_states,
                &q,
                &j,
                &r,
                &forms,
                &observations,
                cancelled,
            )?)
        };
        let storage = QuadraticStorage::new(q, states).map_err(HeterogeneousError::Propagation)?;
        let system = PortHamiltonian::new(states, source_count, j, r, g, Box::new(storage))
            .map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: states,
                max_omega_dt: 1.0,
            },
        )
        .map_err(HeterogeneousError::Propagation)?;
        let flow = match extension {
            None => PreparedQuadraticFlow::new_controlled(
                &system,
                &split,
                duration_s,
                &forms,
                &observations,
                limits.preparation,
                &mut *cancelled,
            ),
            Some(extension) => PreparedQuadraticFlow::new_extended_controlled(
                &system,
                &split,
                &extension.generator,
                &extension.scales,
                duration_s,
                &extension.forms,
                &extension.observations,
                limits.preparation,
                &mut *cancelled,
            ),
        }
        .map_err(propagation)?;
        poll(cancelled)?;
        Ok(PreparedExactMechanicalTime {
            model: self,
            duration_s,
            pressure_release,
            potential_scale,
            fluid_mass,
            walls,
            admittance,
            states: total_states,
            linear_port_count: ports.len(),
            boundary_port_count,
            volume_source_count: volume_sources.len(),
            wall_traction_count: wall_tractions.len(),
            flow,
            limits,
        })
    }
}

impl<'m, 'f, 'c> PreparedExactMechanicalTime<'m, 'f, 'c> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Complete evolved extent, including every original wall memory.
    #[must_use]
    pub const fn state_count(&self) -> usize {
        self.states
    }
    #[must_use]
    pub const fn port_count(&self) -> usize {
        self.boundary_port_count
    }
    #[must_use]
    pub const fn volume_source_count(&self) -> usize {
        self.volume_source_count
    }
    /// Number of independently driven original wall-traction shapes.
    #[must_use]
    pub const fn wall_traction_count(&self) -> usize {
        self.wall_traction_count
    }
    #[must_use]
    pub fn wall_count(&self) -> usize {
        self.walls.len()
    }
    /// Active original hierarchical trace coordinates, sorted by global DOF.
    #[must_use]
    pub fn wall_dofs(&self, wall: usize) -> Option<&[usize]> {
        self.walls.get(wall).map(|w| w.dofs.as_slice())
    }
    #[must_use]
    pub fn wall(&self, wall: usize) -> Option<&MechanicalBoundary> {
        self.walls.get(wall).map(|w| &w.boundary)
    }
    #[must_use]
    pub fn admittance_terms(&self) -> usize {
        self.admittance.len()
    }
    #[must_use]
    pub fn admittance_memory_dofs(&self, arm: usize) -> Option<&[usize]> {
        self.admittance.get(arm).map(|a| a.trace.dofs.as_slice())
    }
    #[must_use]
    pub fn admittance_memory_patch_and_term(&self, arm: usize) -> Option<(usize, usize)> {
        self.admittance
            .get(arm)
            .map(|a| (a.trace.patch, a.trace.term))
    }
    #[must_use]
    pub fn admittance_relaxation_time_s(&self, arm: usize) -> Option<f64> {
        self.admittance.get(arm).map(|a| a.trace.tau)
    }
    #[must_use]
    pub fn zero_state(&self) -> ExactMechanicalTimeState<'_, 'm, 'f, 'c> {
        ExactMechanicalTimeState {
            owner: self,
            acoustic: HeterogeneousTimeState {
                model: self.model,
                potential: vec![0.0; self.model.ndof()],
                pressure: vec![0.0; self.model.ndof()],
            },
            displacements: self.walls.iter().map(|w| vec![0.0; w.dofs.len()]).collect(),
            velocities: self.walls.iter().map(|w| vec![0.0; w.dofs.len()]).collect(),
            memories: self
                .walls
                .iter()
                .map(|w| vec![vec![0.0; w.dofs.len()]; w.couplings.len()])
                .collect(),
            admittance_memories: self
                .admittance
                .iter()
                .map(|a| vec![0.0; a.trace.dofs.len()])
                .collect(),
        }
    }
    /// Import an acoustic field with every wall and memory explicitly at rest.
    pub fn state(
        &self,
        potential: &[f64],
        pressure: &[f64],
    ) -> Result<ExactMechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        self.admit_acoustic(potential, pressure)?;
        let mut state = self.zero_state();
        state.acoustic.potential.copy_from_slice(potential);
        state.acoustic.pressure.copy_from_slice(pressure);
        Ok(state)
    }
    /// Import nonhereditary walls. A Maxwell law requires all memory histories.
    pub fn state_with_walls(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
    ) -> Result<ExactMechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        if self.walls.iter().any(|w| !w.couplings.is_empty()) {
            return Err(invalid("complete original Maxwell wall memory is required").into());
        }
        self.state_with_wall_memory(
            potential,
            pressure,
            displacements_m,
            velocities_m_s,
            &vec![Vec::new(); self.walls.len()],
        )
    }
    /// Import original physical coefficients and ALL causal memory. Arrays use
    /// wall, original arm, then wall_dofs order. No history is fitted or reset.
    pub fn state_with_wall_memory(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
        memories: &[Vec<Vec<f64>>],
    ) -> Result<ExactMechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        if !self.admittance.is_empty() {
            return Err(
                invalid("complete original spatial admittance histories are required").into(),
            );
        }
        self.state_with_all_memories(
            potential,
            pressure,
            displacements_m,
            velocities_m_s,
            memories,
            &[],
        )
    }
    /// Import every original physical field and both distinct memory families.
    /// Mechanical memories are energy-normalized z; admittance memories are q
    /// pressures [Pa]. Neither history is reset or pressure-release projected.
    #[allow(clippy::too_many_arguments)]
    pub fn state_with_all_memories(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
        memories: &[Vec<Vec<f64>>],
        admittance_memories: &[Vec<f64>],
    ) -> Result<ExactMechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        self.admit_acoustic(potential, pressure)?;
        if displacements_m.len() != self.walls.len()
            || velocities_m_s.len() != self.walls.len()
            || memories.len() != self.walls.len()
            || self.walls.iter().enumerate().any(|(index, wall)| {
                displacements_m[index].len() != wall.dofs.len()
                    || velocities_m_s[index].len() != wall.dofs.len()
                    || displacements_m[index]
                        .iter()
                        .chain(&velocities_m_s[index])
                        .any(|v| !v.is_finite())
                    || memories[index].len() != wall.couplings.len()
                    || memories[index].iter().any(|arm| {
                        arm.len() != wall.dofs.len() || arm.iter().any(|z| !z.is_finite())
                    })
            })
        {
            return Err(invalid(
                "finite complete original wall displacement, velocity and memory fields required",
            )
            .into());
        }
        if admittance_memories.len() != self.admittance.len()
            || self
                .admittance
                .iter()
                .zip(admittance_memories)
                .any(|(arm, memory)| {
                    memory.len() != arm.trace.dofs.len() || memory.iter().any(|v| !v.is_finite())
                })
        {
            return Err(
                invalid("finite complete original admittance memories are required").into(),
            );
        }
        Ok(ExactMechanicalTimeState {
            owner: self,
            acoustic: HeterogeneousTimeState {
                model: self.model,
                potential: potential.to_vec(),
                pressure: pressure.to_vec(),
            },
            displacements: displacements_m.to_vec(),
            velocities: velocities_m_s.to_vec(),
            memories: memories.to_vec(),
            admittance_memories: admittance_memories.to_vec(),
        })
    }
    fn admit_acoustic(
        &self,
        potential: &[f64],
        pressure: &[f64],
    ) -> Result<(), HeterogeneousError> {
        if potential.len() != self.model.ndof()
            || pressure.len() != self.model.ndof()
            || potential.iter().chain(pressure).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite complete original acoustic potential and pressure required",
            )
            .into());
        }
        for (dof, &fixed) in self.pressure_release.constrained.iter().enumerate() {
            if fixed && (potential[dof] != 0.0 || pressure[dof] != 0.0) {
                return Err(invalid("exact wall pressure-release state requires zero potential and pressure on every trace DOF").into());
            }
        }
        Ok(())
    }
    fn validate(
        &self,
        state: &ExactMechanicalTimeState<'_, '_, '_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self, state.owner) {
            return Err(invalid(
                "exact mechanical state belongs to another prepared wall configuration",
            )
            .into());
        }
        Ok(())
    }
    /// Independently evaluated original sparse fluid and surface storage [J].
    pub fn energy(
        &self,
        state: &ExactMechanicalTimeState<'_, '_, '_, '_>,
    ) -> Result<MechanicalTimeEnergy, HeterogeneousError> {
        self.energy_with_roundoff(state, &mut || false).map(|v| v.0)
    }
    fn energy_with_roundoff(
        &self,
        state: &ExactMechanicalTimeState<'_, '_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(MechanicalTimeEnergy, f64), HeterogeneousError> {
        self.validate(state)?;
        let (kinetic, rk) =
            positive_form(&self.model.stiffness, &state.acoustic.potential, cancelled)?;
        let (compression, rm) =
            positive_form(&self.model.mass, &state.acoustic.pressure, cancelled)?;
        let mut energy = MechanicalTimeEnergy {
            acoustic: HeterogeneousTimeEnergy {
                kinetic_j: 0.5 * kinetic,
                compression_j: 0.5 * compression,
            },
            wall_kinetic_j: 0.0,
            wall_elastic_j: 0.0,
            wall_relaxation_j: 0.0,
            admittance_memory_j: 0.0,
        };
        let mut roundoff = 0.5 * (rk + rm);
        for (index, wall) in self.walls.iter().enumerate() {
            let (kinetic, rk) = positive_form(&wall.mass, &state.velocities[index], cancelled)?;
            let (elastic, re) = positive_form(&wall.mass, &state.displacements[index], cancelled)?;
            let pin = wall.boundary.wall;
            energy.wall_kinetic_j += 0.5 * pin.surface_density * kinetic;
            energy.wall_elastic_j += 0.5 * pin.stiffness_per_area * elastic;
            roundoff += 0.5 * (pin.surface_density * rk + pin.stiffness_per_area * re);
            for (arm, &coupling) in wall.couplings.iter().enumerate() {
                let extension = extensions(
                    coupling,
                    &state.displacements[index],
                    &state.memories[index][arm],
                );
                let (stored, error) = positive_form(&wall.mass, &extension, cancelled)?;
                energy.wall_relaxation_j += 0.5 * stored;
                roundoff += 0.5 * error;
            }
        }
        for (arm, memory) in self.admittance.iter().zip(&state.admittance_memories) {
            let (stored, error) = positive_form(&arm.trace.conductance, memory, cancelled)?;
            energy.admittance_memory_j += 0.5 * arm.trace.tau * stored;
            roundoff += 0.5 * arm.trace.tau * error;
        }
        if !energy.total_j().is_finite() || !roundoff.is_finite() {
            return Err(invalid("exact mechanical acoustic storage is unrepresentable").into());
        }
        Ok((energy, roundoff))
    }
    /// Advance every field with boundary velocities held over this interval.
    /// Amplitudes are linear ports then prepared boundary sources, all [m/s].
    pub fn step<'p>(
        &'p self,
        initial: &ExactMechanicalTimeState<'_, '_, '_, '_>,
        amplitude_m_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactMechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_volume_sources(initial, amplitude_m_s, &[], cancelled)
    }
    /// Advance the SAME coupled trajectory under held boundary and volume inputs.
    /// Every prepared column needs an explicit amplitude, including zero. Work
    /// uses its integrated conjugate pressure, including coherent interference
    /// and negative work; the returned candidate leaves the caller's clock alone.
    pub fn step_with_volume_sources<'p>(
        &'p self,
        initial: &ExactMechanicalTimeState<'_, '_, '_, '_>,
        amplitude_m_s: &[f64],
        volume_flow_m3_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactMechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_wall_tractions(initial, amplitude_m_s, volume_flow_m3_s, &[], cancelled)
    }

    /// Advance the SAME full coupled trajectory with held acoustic inputs and
    /// applied wall tractions [Pa], in preparation order. All prepared amplitudes
    /// must be explicit, including zero. Positive traction acts outward; signed
    /// mechanical work uses its integrated wall-volume-displacement observation.
    /// The original fluid-to-wall work remains an independently integrated
    /// internal exchange. This method does no surface assembly or factorization.
    pub fn step_with_wall_tractions<'p>(
        &'p self,
        initial: &ExactMechanicalTimeState<'_, '_, '_, '_>,
        amplitude_m_s: &[f64],
        volume_flow_m3_s: &[f64],
        traction_pa: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactMechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        self.validate(initial)?;
        if amplitude_m_s.len() != self.boundary_port_count
            || volume_flow_m3_s.len() != self.volume_source_count
            || traction_pa.len() != self.wall_traction_count
            || amplitude_m_s
                .iter()
                .chain(volume_flow_m3_s)
                .chain(traction_pa)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite held boundary, volume-source and wall-traction amplitudes required",
            )
            .into());
        }
        let combined;
        let amplitudes = if self.volume_source_count == 0 && self.wall_traction_count == 0 {
            amplitude_m_s
        } else {
            combined = amplitude_m_s
                .iter()
                .chain(volume_flow_m3_s)
                .chain(traction_pa)
                .copied()
                .collect::<Vec<_>>();
            &combined
        };
        let (initial_energy, initial_roundoff) = self.energy_with_roundoff(initial, cancelled)?;
        let n = self.pressure_release.free_dofs.len();
        let mut x = vec![0.0; self.states];
        for (row, &global) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            x[row] = self.potential_scale[row] * initial.acoustic.potential[global];
            let mut momentum = Dd::from_f64(0.0);
            for col in 0..n {
                momentum = momentum
                    + Dd::from_f64(self.fluid_mass.normalized[row * n + col])
                        * Dd::from_f64(self.fluid_mass.scales[col])
                        * Dd::from_f64(
                            initial.acoustic.pressure[self.pressure_release.free_dofs[col]],
                        );
            }
            x[n + row] = momentum.to_f64();
        }
        for (index, wall) in self.walls.iter().enumerate() {
            let m = wall.dofs.len();
            for row in 0..m {
                poll(cancelled)?;
                x[wall.offset + row] =
                    wall.stiffness_scale * wall.map.scales[row] * initial.displacements[index][row];
                let mut momentum = Dd::from_f64(0.0);
                for col in 0..m {
                    momentum = momentum
                        + Dd::from_f64(wall.map.normalized[row * m + col])
                            * Dd::from_f64(wall.map.scales[col])
                            * Dd::from_f64(initial.velocities[index][col]);
                }
                x[wall.offset + m + row] = wall.inertia_scale * momentum.to_f64();
                for (arm, &coupling) in wall.couplings.iter().enumerate() {
                    x[wall.offset + (2 + arm) * m + row] = wall.map.scales[row]
                        * coupling.mul_add(
                            initial.displacements[index][row],
                            -initial.memories[index][arm][row],
                        );
                }
            }
        }
        for (arm, memory) in self.admittance.iter().zip(&initial.admittance_memories) {
            poll(cancelled)?;
            x[arm.offset..arm.offset + memory.len()].copy_from_slice(memory);
        }
        let result = self
            .flow
            .evaluate(&x, amplitudes)
            .map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let mut next = self.zero_state();
        for (row, &global) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            next.acoustic.potential[global] = result.state[row] / self.potential_scale[row];
            let mut pressure = Dd::from_f64(0.0);
            for col in 0..n {
                pressure = pressure
                    + Dd::from_f64(self.fluid_mass.inverse[row * n + col])
                        * Dd::from_f64(result.state[n + col]);
            }
            next.acoustic.pressure[global] = pressure.to_f64() / self.fluid_mass.scales[row];
        }
        for (index, wall) in self.walls.iter().enumerate() {
            let m = wall.dofs.len();
            for row in 0..m {
                poll(cancelled)?;
                let w =
                    result.state[wall.offset + row] / wall.stiffness_scale / wall.map.scales[row];
                next.displacements[index][row] = w;
                let mut velocity = Dd::from_f64(0.0);
                for col in 0..m {
                    velocity = velocity
                        + Dd::from_f64(wall.map.inverse[row * m + col])
                            * Dd::from_f64(result.state[wall.offset + m + col]);
                }
                next.velocities[index][row] =
                    velocity.to_f64() / wall.inertia_scale / wall.map.scales[row];
                for (arm, &coupling) in wall.couplings.iter().enumerate() {
                    let extension =
                        result.state[wall.offset + (2 + arm) * m + row] / wall.map.scales[row];
                    next.memories[index][arm][row] = coupling.mul_add(w, -extension);
                }
            }
        }
        for (arm, memory) in self.admittance.iter().zip(&mut next.admittance_memories) {
            poll(cancelled)?;
            let count = memory.len();
            memory.copy_from_slice(&result.state[arm.offset..arm.offset + count]);
        }
        // All endpoint coordinates, including zero-energy gauge and relaxed
        // memory directions, must remain finite before publishing a candidate.
        if next
            .acoustic
            .potential
            .iter()
            .chain(&next.acoustic.pressure)
            .chain(next.displacements.iter().flatten())
            .chain(next.velocities.iter().flatten())
            .chain(next.memories.iter().flatten().flatten())
            .chain(next.admittance_memories.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(
                invalid("complete exact acoustic and wall endpoint is unrepresentable").into(),
            );
        }
        let (final_energy, final_roundoff) = self.energy_with_roundoff(&next, cancelled)?;
        let acoustic_change = energy_change(
            &self.model.stiffness,
            &initial.acoustic.potential,
            &next.acoustic.potential,
            cancelled,
        )? + energy_change(
            &self.model.mass,
            &initial.acoustic.pressure,
            &next.acoustic.pressure,
            cancelled,
        )?;
        let mut wall_change = Dd::from_f64(0.0);
        for (index, wall) in self.walls.iter().enumerate() {
            wall_change = wall_change
                + Dd::from_f64(wall.boundary.wall.surface_density)
                    * Dd::from_f64(energy_change(
                        &wall.mass,
                        &initial.velocities[index],
                        &next.velocities[index],
                        cancelled,
                    )?)
                + Dd::from_f64(wall.boundary.wall.stiffness_per_area)
                    * Dd::from_f64(energy_change(
                        &wall.mass,
                        &initial.displacements[index],
                        &next.displacements[index],
                        cancelled,
                    )?);
            for (arm, &coupling) in wall.couplings.iter().enumerate() {
                let before = extensions(
                    coupling,
                    &initial.displacements[index],
                    &initial.memories[index][arm],
                );
                let after = extensions(
                    coupling,
                    &next.displacements[index],
                    &next.memories[index][arm],
                );
                wall_change = wall_change
                    + Dd::from_f64(energy_change(&wall.mass, &before, &after, cancelled)?);
            }
        }
        let mut input = Dd::from_f64(0.0);
        let mut admittance_change = Dd::from_f64(0.0);
        for (index, arm) in self.admittance.iter().enumerate() {
            admittance_change = admittance_change
                + Dd::from_f64(arm.trace.tau)
                    * Dd::from_f64(energy_change(
                        &arm.trace.conductance,
                        &initial.admittance_memories[index],
                        &next.admittance_memories[index],
                        cancelled,
                    )?);
        }
        let source_offset = if self.admittance.is_empty() {
            admittance::BASE_FORMS
        } else {
            admittance::FORM_COUNT
        };
        let mut boundary_source_work = Dd::from_f64(0.0);
        let mut volume_work = Dd::from_f64(0.0);
        let mut acoustic_work = Dd::from_f64(0.0);
        let mut traction_work = Dd::from_f64(0.0);
        let traction_offset = self.boundary_port_count + self.volume_source_count;
        let mut absolute_work = 0.0;
        for (index, (&amplitude, &integral)) in amplitudes
            .iter()
            .zip(&result.integrals[source_offset..])
            .enumerate()
        {
            let work = Dd::from_f64(amplitude) * Dd::from_f64(integral);
            input = input + work;
            if index >= traction_offset {
                traction_work = traction_work + work;
            } else {
                acoustic_work = acoustic_work + work;
                if index >= self.boundary_port_count {
                    volume_work = volume_work + work;
                } else if index >= self.linear_port_count {
                    boundary_source_work = boundary_source_work + work;
                }
            }
            absolute_work += work.to_f64().abs();
        }
        let work = ExactMechanicalTimeWork {
            initial_energy_j: initial_energy.total_j(),
            final_energy_j: final_energy.total_j(),
            delta_energy_j: (Dd::from_f64(acoustic_change) + wall_change + admittance_change)
                .to_f64(),
            acoustic_delta_energy_j: acoustic_change,
            wall_delta_energy_j: wall_change.to_f64(),
            admittance_delta_energy_j: admittance_change.to_f64(),
            input_work_j: input.to_f64(),
            acoustic_input_work_j: acoustic_work.to_f64(),
            boundary_source_work_j: boundary_source_work.to_f64(),
            volume_source_work_j: volume_work.to_f64(),
            wall_traction_work_j: traction_work.to_f64(),
            boundary_dissipation_j: result.integrals[0],
            wall_resistance_dissipation_j: result.integrals[1],
            wall_relaxation_dissipation_j: result.integrals[2],
            wall_work_j: result.integrals[3],
            admittance_dissipation_j: if self.admittance.is_empty() {
                0.0
            } else {
                result.integrals[admittance::HEAT]
            },
            admittance_work_j: if self.admittance.is_empty() {
                0.0
            } else {
                result.integrals[admittance::WORK]
            },
        };
        let physical_scale = work
            .initial_energy_j
            .abs()
            .max(work.final_energy_j.abs())
            .max(work.delta_energy_j.abs())
            .max(work.acoustic_delta_energy_j.abs())
            .max(work.wall_delta_energy_j.abs())
            .max(work.admittance_delta_energy_j.abs())
            .max(work.acoustic_input_work_j.abs())
            .max(work.wall_traction_work_j.abs())
            .max(absolute_work)
            .max(work.boundary_dissipation_j.abs())
            .max(work.wall_dissipation_j().abs())
            .max(work.admittance_dissipation_j.abs())
            .max(work.admittance_work_j.abs())
            .max(work.wall_work_j.abs());
        let limit_j = self.limits.energy_balance_tolerance * physical_scale
            + initial_roundoff
            + final_roundoff
            + 64.0 * f64::EPSILON * absolute_work;
        if !limit_j.is_finite()
            || !work.input_work_j.is_finite()
            || !work.acoustic_input_work_j.is_finite()
            || !work.boundary_source_work_j.is_finite()
            || !work.volume_source_work_j.is_finite()
            || !work.wall_traction_work_j.is_finite()
            || result.integrals.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("exact mechanical acoustic work is unrepresentable").into());
        }
        for residual in [
            work.balance_residual_j(),
            work.acoustic_balance_residual_j(),
            work.wall_balance_residual_j(),
            work.admittance_balance_residual_j(),
        ] {
            if !residual.is_finite() {
                return Err(
                    invalid("exact mechanical acoustic energy balance is unrepresentable").into(),
                );
            }
            if residual.abs() > limit_j {
                return Err(HeterogeneousError::EnergyBalance {
                    observed_j: residual,
                    limit_j,
                });
            }
        }
        for heat in [
            work.boundary_dissipation_j,
            work.wall_resistance_dissipation_j,
            work.wall_relaxation_dissipation_j,
            work.admittance_dissipation_j,
        ] {
            if heat < -limit_j {
                return Err(HeterogeneousError::EnergyBalance {
                    observed_j: heat,
                    limit_j,
                });
            }
        }
        poll(cancelled)?;
        Ok(ExactMechanicalTimeStep {
            state: next,
            work,
            integrated_port_pressure_pa_m2_s: result.integrals
                [source_offset..source_offset + self.boundary_port_count]
                .to_vec(),
            integrated_volume_source_pressure_pa_s: result.integrals
                [source_offset + self.boundary_port_count..source_offset + traction_offset]
                .to_vec(),
            integrated_wall_traction_displacement_m3: result.integrals
                [source_offset + traction_offset..]
                .to_vec(),
        })
    }
}

fn extensions(coupling: f64, displacement: &[f64], memory: &[f64]) -> Vec<f64> {
    displacement
        .iter()
        .zip(memory)
        .map(|(&w, &z)| coupling.mul_add(w, -z))
        .collect()
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod composed_tests;

#[cfg(test)]
mod traction_tests;
