//! Exact held-input evolution of the complete fluid and spatial wall memories.
//!
//! The original equations are `M p_dot + K psi + C0 p + sum Cj qj = b`,
//! `psi_dot = p`, and `tau_j qj_dot + qj = trace(p)`. Canonical fluid
//! coordinates reuse the existing consistent-mass map; every original trace
//! memory remains an independent pressure coordinate, including zero-residue
//! and pressure-release-shared coordinates. No residue form is inverted.
//!
//! The existing complete-generator extension of `PreparedQuadraticFlow`
//! prepares ONE complete operator and same-trajectory Gramian observations. Its
//! added-coordinate scales are numerical only. Physical storage is evaluated
//! from the original fluid forms and `tau_j qj^T Cj qj/2`; the complete
//! generator's corresponding energy identity is checked before publication.
//! Independent forms integrate instantaneous heat, relaxation heat, and
//! signed fluid-to-memory work. No effective midpoint admittance, terminal
//! energy deficit, or endpoint quadrature supplies a physical heat channel.
//!
//! Exactness concerns this semidiscrete held-input model up to numerical
//! exponential error. Dense preparation has explicit state/byte/work limits;
//! subsequent steps reuse its cached contractions. Spatial and material
//! quadrature convergence, source bandwidth and output sampling remain separate.

use super::super::exact::{ExactTimeLimits, ExactVelocityPort};
use super::super::mechanical::exact::{
    MassMap, add_skew, diagonal_scales, pullback, put_block,
};
use super::super::mechanical::balance;
use super::super::PressureReleaseConstraints;
use super::super::super::{boundary_source, source};
use super::{
    BoundaryAdmittance, BoundaryVelocitySource, C64, Coo, Csr, Dd, HeterogeneousError,
    HeterogeneousHelmholtzFe, HeterogeneousTimeState, RelaxingTimeEnergy,
    ResistiveBoundary, SpatialRelaxingBoundary, SpatialResistiveBoundary,
    VolumeVelocitySource, energy_change, invalid, poll, positive_form,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use fs_math::det;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    prepared::PreparedQuadraticFlow, split_quadratic,
};
use fs_phs::{PortHamiltonian, PreparedStepError, QuadraticStorage};

/// Admission counts every free fluid and original trace-memory coordinate.
#[derive(Clone, Copy, Debug)]
pub struct ExactRelaxingTimeLimits {
    /// Two coordinates per free pressure DOF plus all original trace memories.
    /// Original sparse forms, constraint maps and lifted vectors are additional.
    pub max_states: usize,
    /// Native exponential/Gramian dimension, byte and multiplication-work caps.
    /// Physical forms and coordinate maps additionally occupy O(max_states²).
    pub preparation: DenseExponentialBudget,
    pub energy_balance_tolerance: f64,
}

impl Default for ExactRelaxingTimeLimits {
    fn default() -> Self {
        let acoustic = ExactTimeLimits::default();
        Self {
            max_states: 2 * acoustic.max_dofs,
            preparation: acoustic.preparation,
            energy_balance_tolerance: acoustic.energy_balance_tolerance,
        }
    }
}

struct ExactArm {
    patch: usize,
    term: usize,
    dofs: Vec<usize>,
    conductance: Csr,
    tau: f64,
    offset: usize,
}

/// Immutable original operator, coordinate maps and complete trajectory cache.
pub struct PreparedExactRelaxingTime<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    duration_s: f64,
    pressure_release: PressureReleaseConstraints,
    potential_scale: Vec<f64>,
    fluid_mass: MassMap,
    arms: Vec<ExactArm>,
    states: usize,
    linear_port_count: usize,
    boundary_port_count: usize,
    volume_source_count: usize,
    flow: PreparedQuadraticFlow,
    limits: ExactRelaxingTimeLimits,
}

/// Physical fields bound to this exact preparation and its original memories.
#[derive(Clone)]
pub struct ExactRelaxingTimeState<'p, 'm, 'f, 'c> {
    owner: &'p PreparedExactRelaxingTime<'m, 'f, 'c>,
    acoustic: HeterogeneousTimeState<'m, 'f, 'c>,
    memories: Vec<Vec<f64>>,
}

impl core::fmt::Debug for ExactRelaxingTimeState<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ExactRelaxingTimeState")
            .field("acoustic", &self.acoustic)
            .field("relaxation_terms", &self.memories.len())
            .finish_non_exhaustive()
    }
}

impl<'m, 'f, 'c> ExactRelaxingTimeState<'_, 'm, 'f, 'c> {
    pub fn acoustic_state(&self) -> &HeterogeneousTimeState<'m, 'f, 'c> { &self.acoustic }
    pub fn pressure_pa(&self) -> &[f64] { self.acoustic.pressure_pa() }
    pub fn pressure_potential_pa_s(&self) -> &[f64] { self.acoustic.pressure_potential_pa_s() }
    /// Original q_j [Pa], in the prepared arm's `memory_dofs` order.
    pub fn memory_pressure_pa(&self, arm: usize) -> Option<&[f64]> {
        self.memories.get(arm).map(Vec::as_slice)
    }
}

/// Independent physical accounts on the same complete held-input trajectory.
#[derive(Clone, Copy, Debug)]
pub struct ExactRelaxingTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    pub acoustic_delta_energy_j: f64,
    pub memory_delta_energy_j: f64,
    pub input_work_j: f64,
    /// Prepared boundary-source subtotal, already included in input work.
    pub boundary_source_work_j: f64,
    /// Volume-source subtotal, already included in input work.
    pub volume_source_work_j: f64,
    pub instantaneous_dissipation_j: f64,
    pub relaxation_dissipation_j: f64,
    /// Signed integral of sum p^T Cj qj, positive FROM fluid TO memories.
    pub memory_work_j: f64,
}

impl ExactRelaxingTimeWork {
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j + self.instantaneous_dissipation_j
            + self.relaxation_dissipation_j - self.input_work_j
    }
    pub fn acoustic_balance_residual_j(self) -> f64 {
        self.acoustic_delta_energy_j + self.instantaneous_dissipation_j
            + self.memory_work_j - self.input_work_j
    }
    pub fn memory_balance_residual_j(self) -> f64 {
        self.memory_delta_energy_j + self.relaxation_dissipation_j - self.memory_work_j
    }
}

#[derive(Clone, Debug)]
pub struct ExactRelaxingTimeStep<'p, 'm, 'f, 'c> {
    pub state: ExactRelaxingTimeState<'p, 'm, 'f, 'c>,
    pub work: ExactRelaxingTimeWork,
    /// Integral of each boundary shape against p [Pa m² s], linear ports first
    /// and prepared boundary sources second, in preparation order.
    pub integrated_port_pressure_pa_m2_s: Vec<f64>,
    /// Integral of each normalized volume source against p [Pa s].
    pub integrated_volume_source_pressure_pa_s: Vec<f64>,
}

fn propagation(error: PreparedStepError) -> HeterogeneousError {
    match error {
        PreparedStepError::Cancelled => HeterogeneousError::Cancelled,
        PreparedStepError::Solver(error) => HeterogeneousError::Propagation(error),
    }
}

fn admit_extent(states: usize, sources: usize, limits: ExactRelaxingTimeLimits) -> Result<(), HeterogeneousError> {
    if states > limits.max_states {
        return Err(HeterogeneousError::Limit {
            what: "complete exact fluid and spatial wall-memory states",
            requested: states, maximum: limits.max_states,
        });
    }
    let augmented = states.checked_add(sources).and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("exact wall-memory source extent overflow"))?;
    let maximum = limits.preparation.max_dimension.min(MAX_LINEAR_EXPONENTIAL_BLOCK_DIM);
    if augmented > maximum {
        return Err(HeterogeneousError::Limit {
            what: "complete exact wall-memory augmented states", requested: augmented, maximum,
        });
    }
    Ok(())
}

// Audit the COMPLETE generator against original physical storage and heat,
// independently of the backbone pHS metric used for numerical coordinates.
fn admit_energy_identity(
    storage: &[f64], generator: &[f64], heat: &[Vec<f64>], n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), HeterogeneousError> {
    for row in 0..n {
        poll(cancelled)?;
        for col in 0..n {
            let mut sum = Dd::from_f64(0.0);
            let mut scale = 0.0;
            for k in 0..n {
                for term in [
                    Dd::from_f64(storage[row*n+k]) * Dd::from_f64(generator[k*n+col]),
                    Dd::from_f64(generator[k*n+row]) * Dd::from_f64(storage[k*n+col]),
                ] { sum = sum + term; scale += term.to_f64().abs(); }
            }
            for form in &heat[..2] {
                let term = 2.0 * form[row*n+col];
                sum = sum + Dd::from_f64(term); scale += term.abs();
            }
            if !sum.to_f64().is_finite() || !scale.is_finite()
                || sum.to_f64().abs() > 2048.0 * f64::EPSILON * n as f64 * scale {
                return Err(invalid("complete wall-memory generator violates original storage/heat identity").into());
            }
        }
    }
    Ok(())
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare exact complete-state evolution of the ORIGINAL causal wall law.
    ///
    /// Every pressure potential, pressure, and geometric trace memory is
    /// retained. Zero residues remain valid and their memories still evolve.
    /// Held linear/quadratic boundary and volume-source profiles use the same
    /// trajectory for physical work. All physical patches must be disjoint;
    /// sources may share any non-release face. No midpoint condensate is used.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_relaxing_boundaries<'m>(
        &'m self, duration_s: f64, uniform: &[ResistiveBoundary],
        spatial: &[&SpatialResistiveBoundary<'_, '_, '_>],
        relaxing: &[&SpatialRelaxingBoundary<'_, '_, '_>],
        ports: &[ExactVelocityPort], boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]], limits: ExactRelaxingTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactRelaxingTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        boundary_source::admit_count(boundary_sources.len())?;
        source::admit_count(volume_sources.len())?;
        if !duration_s.is_finite() || duration_s <= 0.0 || limits.max_states == 0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0 || limits.energy_balance_tolerance > 1e-6 {
            return Err(invalid("positive exact clock, bounded energy tolerance and complete-state limits required").into());
        }
        let boundary_port_count = ports.len().checked_add(boundary_sources.len())
            .ok_or_else(|| invalid("exact wall-memory boundary source extent overflow"))?;
        let source_count = boundary_port_count.checked_add(volume_sources.len())
            .ok_or_else(|| invalid("exact wall-memory source extent overflow"))?;
        if pressure_release_faces.is_empty() {
            admit_extent(self.ndof().checked_mul(2).ok_or_else(|| invalid("exact fluid state extent overflow"))?, source_count, limits)?;
        }
        let pressure_release = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let n = pressure_release.free_dofs.len();
        let base_states = 2*n;
        let mut states = base_states;
        admit_extent(states, source_count, limits)?;
        let physical_fields = spatial.iter().copied()
            .chain(relaxing.iter().map(|law| &law.instantaneous)).collect::<Vec<_>>();
        let instantaneous = self.resistive_matrix(uniform, &physical_fields, &pressure_release, cancelled)?;
        let instantaneous = pressure_release.restrict(&instantaneous, cancelled)?;
        let mut arms = Vec::new();
        for (patch, law) in relaxing.iter().enumerate() {
            poll(cancelled)?;
            let (trace, _) = self.admittance_matrices(&[BoundaryAdmittance {
                faces: law.faces().to_vec(), admittance_m_pa_s: C64::ONE,
            }], cancelled)?;
            let dofs = (0..self.ndof()).filter(|&row| {
                let (cols, values) = trace.row(row);
                cols.binary_search(&row).is_ok_and(|i| values[i] > 0.0)
            }).collect::<Vec<_>>();
            for (term, original) in law.terms.iter().enumerate() {
                poll(cancelled)?;
                let offset = states;
                states = states.checked_add(dofs.len()).ok_or_else(|| invalid("exact wall-memory state extent overflow"))?;
                admit_extent(states, source_count, limits)?;
                let mut local = Coo::new(dofs.len(), dofs.len());
                for (row, &global) in dofs.iter().enumerate() {
                    let (cols, values) = original.conductance.row(global);
                    for (&col, &value) in cols.iter().zip(values) {
                        if value != 0.0 {
                            let column = dofs.binary_search(&col).map_err(|_| invalid("incomplete exact wall-memory trace support"))?;
                            local.push(row, column, value);
                        }
                    }
                }
                arms.push(ExactArm { patch, term, dofs: dofs.clone(), conductance: local.assemble(), tau: original.relaxation_time_s, offset });
            }
        }
        let mut loads = Vec::with_capacity(source_count);
        for port in ports {
            poll(cancelled)?;
            if port.profile.is_empty() { return Err(invalid("exact acoustic source requires a nonempty trace profile").into()); }
            for &(face, _) in &port.profile {
                let mut key = face; key.sort_unstable();
                if pressure_release.faces.contains(&key) { return Err(invalid("pressure-release and exact velocity-port faces must not overlap").into()); }
            }
            loads.push(self.geometry.boundary_linear_load(&port.profile)?);
        }
        for source in boundary_sources {
            poll(cancelled)?; source.validate(self)?;
            for face in &source.faces {
                if pressure_release.faces.contains(face) { return Err(invalid("pressure-release and prepared boundary-source faces must not overlap").into()); }
            }
            loads.push(source.load.clone());
        }
        for source in volume_sources { poll(cancelled)?; source.validate(self)?; loads.push(source.load.clone()); }
        let fluid_mass = MassMap::new(&self.mass, &pressure_release.free_dofs, cancelled)?;
        let potential_scale = diagonal_scales(&self.stiffness, &pressure_release.free_dofs, cancelled)?;
        let pressure_map = fluid_mass.physical_map(1.0);
        let mut q = vec![0.0; base_states*base_states];
        let mut j = vec![0.0; base_states*base_states];
        let mut r = vec![0.0; base_states*base_states];
        let mut g = vec![0.0; base_states*source_count];
        for (row, &global) in pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let (cols, values) = self.stiffness.row(global);
            for (&col, &value) in cols.iter().zip(values) {
                let col = pressure_release.original_to_free[col];
                if col != usize::MAX { q[row*base_states+col] = value / potential_scale[row] / potential_scale[col]; }
            }
            add_skew(&mut j, base_states, row, n+row, potential_scale[row]/fluid_mass.scales[row]);
            let (cols, values) = instantaneous.row(row);
            for (&col, &value) in cols.iter().zip(values) {
                r[(n+row)*base_states+n+col] = value / fluid_mass.scales[row] / fluid_mass.scales[col];
            }
            for (port, load) in loads.iter().enumerate() { g[(n+row)*source_count+port] = load[global]/fluid_mass.scales[row]; }
        }
        put_block(&mut q, base_states, n, &fluid_mass.inverse, n);
        let mut generator = vec![0.0; states*states];
        let mut storage = vec![0.0; states*states];
        for row in 0..base_states {
            poll(cancelled)?;
            for col in 0..base_states {
                storage[row*states+col] = q[row*base_states+col];
                let mut value = Dd::from_f64(0.0);
                for k in 0..base_states {
                    value = value + (Dd::from_f64(j[row*base_states+k])-Dd::from_f64(r[row*base_states+k])) * Dd::from_f64(q[k*base_states+col]);
                }
                generator[row*states+col] = value.to_f64();
            }
        }
        let mut forms = vec![vec![0.0; states*states]; 3];
        put_block(&mut forms[0], states, n, &pullback(&instantaneous, &pressure_map, 1.0, cancelled)?, n);
        let mut extra_scales = Vec::with_capacity(states-base_states);
        for arm in &arms {
            poll(cancelled)?;
            let m = arm.dofs.len();
            let mut trace = vec![0.0; n*m];
            for (row, &global) in arm.dofs.iter().enumerate() {
                poll(cancelled)?;
                let free = pressure_release.original_to_free[global];
                generator[(arm.offset+row)*states+arm.offset+row] = -1.0/arm.tau;
                if free != usize::MAX {
                    for col in 0..n { generator[(arm.offset+row)*states+n+col] = pressure_map[free*n+col]/arm.tau; }
                }
                let (cols, values) = arm.conductance.row(row);
                for (&col, &value) in cols.iter().zip(values) {
                    forms[1][(arm.offset+row)*states+arm.offset+col] = value;
                    storage[(arm.offset+row)*states+arm.offset+col] = arm.tau*value;
                    if free != usize::MAX {
                        generator[(n+free)*states+arm.offset+col] = -value/fluid_mass.scales[free];
                        trace[free*m+col] = value;
                    }
                }
                // Scale by the original diagonal memory energy when nonzero.
                // A zero residue still has a live driven memory: use a bounded
                // fluid-pressure/time scale solely for numerical conditioning.
                let diagonal = cols.binary_search(&row).ok().map_or(0.0, |i| values[i]);
                let stored_diagonal = arm.tau*diagonal;
                let scale = if diagonal > 0.0 {
                    det::sqrt(stored_diagonal)
                } else {
                    let (mass_cols, mass_values) = self.mass.row(global);
                    let at = mass_cols.binary_search(&global).map_err(|_| invalid("original trace coordinate lacks fluid mass"))?;
                    let scale = arm.tau.max(duration_s);
                    let ratio = (arm.tau/scale) / (arm.tau/scale + duration_s/scale);
                    det::sqrt(mass_values[at]) * det::sqrt(ratio)
                };
                if !scale.is_finite() || scale <= 0.0 { return Err(invalid("exact wall-memory numerical scale is unrepresentable").into()); }
                extra_scales.push(scale);
            }
            for row in 0..n {
                poll(cancelled)?;
                for col in 0..m {
                    let mut value = Dd::from_f64(0.0);
                    for k in 0..n { value = value + Dd::from_f64(pressure_map[k*n+row])*Dd::from_f64(trace[k*m+col]); }
                    let half = 0.5*value.to_f64();
                    forms[2][(n+row)*states+arm.offset+col] = half;
                    forms[2][(arm.offset+col)*states+n+row] = half;
                }
            }
        }
        admit_energy_identity(&storage, &generator, &forms, states, cancelled)?;
        let mut observations = Vec::with_capacity(source_count);
        for load in &loads {
            poll(cancelled)?;
            let mut observation = vec![0.0; states];
            for col in 0..n {
                let mut value = Dd::from_f64(0.0);
                for row in 0..n { value = value + Dd::from_f64(load[pressure_release.free_dofs[row]])*Dd::from_f64(pressure_map[row*n+col]); }
                observation[n+col] = value.to_f64();
            }
            observations.push(observation);
        }
        let storage = QuadraticStorage::new(q, base_states).map_err(HeterogeneousError::Propagation)?;
        let system = PortHamiltonian::new(base_states, source_count, j, r, g, Box::new(storage)).map_err(HeterogeneousError::Propagation)?;
        let split = split_quadratic(&system, LinearSplitBudget { max_modes: base_states, max_omega_dt: 1.0 }).map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let flow = PreparedQuadraticFlow::new_extended_controlled(
            &system, &split, &generator, &extra_scales, duration_s, &forms,
            &observations, limits.preparation, &mut *cancelled,
        ).map_err(propagation)?;
        poll(cancelled)?;
        Ok(PreparedExactRelaxingTime {
            model: self, duration_s, pressure_release, potential_scale, fluid_mass,
            arms, states, linear_port_count: ports.len(), boundary_port_count,
            volume_source_count: volume_sources.len(), flow, limits,
        })
    }
}

impl<'m, 'f, 'c> PreparedExactRelaxingTime<'m, 'f, 'c> {
    pub fn duration_s(&self) -> f64 { self.duration_s }
    pub fn relaxation_terms(&self) -> usize { self.arms.len() }
    pub fn memory_dofs(&self, arm: usize) -> Option<&[usize]> { self.arms.get(arm).map(|a| a.dofs.as_slice()) }
    pub fn memory_patch_and_term(&self, arm: usize) -> Option<(usize, usize)> { self.arms.get(arm).map(|a| (a.patch, a.term)) }
    pub fn memory_relaxation_time_s(&self, arm: usize) -> Option<f64> { self.arms.get(arm).map(|a| a.tau) }

    pub fn state<'p>(&'p self, potential: &[f64], pressure: &[f64]) -> Result<ExactRelaxingTimeState<'p, 'm, 'f, 'c>, HeterogeneousError> {
        let memories = self.arms.iter().map(|arm| vec![0.0; arm.dofs.len()]).collect::<Vec<_>>();
        self.state_with_memories(potential, pressure, &memories.iter().map(Vec::as_slice).collect::<Vec<_>>())
    }

    /// Import every original trace memory explicitly, including coordinates
    /// shared with pressure-release faces. Import never projects a nonzero
    /// pressure/potential release trace to zero.
    pub fn state_with_memories<'p>(
        &'p self, potential: &[f64], pressure: &[f64], memories: &[&[f64]],
    ) -> Result<ExactRelaxingTimeState<'p, 'm, 'f, 'c>, HeterogeneousError> {
        let n = self.model.ndof();
        if potential.len() != n || pressure.len() != n
            || potential.iter().chain(pressure).any(|v| !v.is_finite())
            || memories.len() != self.arms.len()
            || memories.iter().zip(&self.arms).any(|(q, arm)| q.len()!=arm.dofs.len() || q.iter().any(|v| !v.is_finite())) {
            return Err(invalid("finite complete original exact acoustic and memory coordinates required").into());
        }
        for (i, &fixed) in self.pressure_release.constrained.iter().enumerate() {
            if fixed && (potential[i]!=0.0 || pressure[i]!=0.0) { return Err(invalid("exact release trace requires zero potential and pressure").into()); }
        }
        let state = ExactRelaxingTimeState {
            owner: self, acoustic: HeterogeneousTimeState { model: self.model, potential: potential.to_vec(), pressure: pressure.to_vec() },
            memories: memories.iter().map(|q| q.to_vec()).collect(),
        };
        self.energy(&state)?;
        Ok(state)
    }

    fn validate(&self, state: &ExactRelaxingTimeState<'_, '_, '_, '_>) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self, state.owner) { return Err(invalid("exact wall-memory state belongs to another preparation").into()); }
        Ok(())
    }

    pub fn energy(&self, state: &ExactRelaxingTimeState<'_, '_, '_, '_>) -> Result<RelaxingTimeEnergy, HeterogeneousError> {
        self.validate(state)?;
        self.energy_with_roundoff(state, &mut || false).map(|v| v.0)
    }

    fn energy_with_roundoff(
        &self, state: &ExactRelaxingTimeState<'_, '_, '_, '_>, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(RelaxingTimeEnergy, f64, f64), HeterogeneousError> {
        let (kinetic, rk) = positive_form(&self.model.stiffness, &state.acoustic.potential, cancelled)?;
        let (compression, rc) = positive_form(&self.model.mass, &state.acoustic.pressure, cancelled)?;
        let mut memory = Dd::from_f64(0.0);
        let mut memory_roundoff = 0.0;
        for (arm, q) in self.arms.iter().zip(&state.memories) {
            let (stored, roundoff) = positive_form(&arm.conductance, q, cancelled)?;
            memory = memory + Dd::from_f64(0.5)*Dd::from_f64(arm.tau)*Dd::from_f64(stored);
            memory_roundoff += 0.5*arm.tau*roundoff;
        }
        let energy = RelaxingTimeEnergy {
            acoustic: super::HeterogeneousTimeEnergy { kinetic_j: 0.5*kinetic, compression_j: 0.5*compression },
            boundary_memory_j: memory.to_f64(),
        };
        if !energy.total_j().is_finite() || !memory_roundoff.is_finite() { return Err(invalid("exact wall-memory physical storage overflow").into()); }
        Ok((energy, 0.5*(rk+rc), memory_roundoff))
    }

    pub fn step<'p>(
        &'p self, initial: &ExactRelaxingTimeState<'_, '_, '_, '_>, amplitude_m_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactRelaxingTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_volume_sources(initial, amplitude_m_s, &[], cancelled)
    }

    /// Advance one held-input interval, returning independently integrated
    /// boundary/volume work, instantaneous heat, memory heat and exchange.
    /// A refusal never modifies the caller's accepted state.
    pub fn step_with_volume_sources<'p>(
        &'p self, initial: &ExactRelaxingTimeState<'_, '_, '_, '_>, amplitude_m_s: &[f64],
        volume_flow_m3_s: &[f64], cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactRelaxingTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?; self.validate(initial)?;
        if amplitude_m_s.len()!=self.boundary_port_count || volume_flow_m3_s.len()!=self.volume_source_count
            || amplitude_m_s.iter().chain(volume_flow_m3_s).any(|v| !v.is_finite()) {
            return Err(invalid("complete finite exact boundary and volume-source amplitudes required").into());
        }
        let combined;
        let amplitudes = if self.volume_source_count == 0 { amplitude_m_s } else {
            combined = amplitude_m_s.iter().chain(volume_flow_m3_s).copied().collect::<Vec<_>>(); &combined
        };
        let n = self.pressure_release.free_dofs.len();
        let mut x = vec![0.0; self.states];
        for (row, &global) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            x[row] = self.potential_scale[row]*initial.acoustic.potential[global];
            let mut value = Dd::from_f64(0.0);
            for col in 0..n {
                value = value + Dd::from_f64(self.fluid_mass.normalized[row*n+col])
                    * Dd::from_f64(self.fluid_mass.scales[col])
                    * Dd::from_f64(initial.acoustic.pressure[self.pressure_release.free_dofs[col]]);
            }
            x[n+row] = value.to_f64();
        }
        for (arm, memory) in self.arms.iter().zip(&initial.memories) {
            poll(cancelled)?; x[arm.offset..arm.offset+memory.len()].copy_from_slice(memory);
        }
        let result = self.flow.evaluate(&x, amplitudes).map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let mut next = ExactRelaxingTimeState {
            owner: self, acoustic: HeterogeneousTimeState { model: self.model, potential: vec![0.0; self.model.ndof()], pressure: vec![0.0; self.model.ndof()] },
            memories: self.arms.iter().map(|arm| result.state[arm.offset..arm.offset+arm.dofs.len()].to_vec()).collect(),
        };
        for (row, &global) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            next.acoustic.potential[global] = result.state[row]/self.potential_scale[row];
            let mut value = Dd::from_f64(0.0);
            for col in 0..n { value = value + Dd::from_f64(self.fluid_mass.inverse[row*n+col])*Dd::from_f64(result.state[n+col]); }
            next.acoustic.pressure[global] = value.to_f64()/self.fluid_mass.scales[row];
        }
        if next.acoustic.potential.iter().chain(&next.acoustic.pressure)
            .chain(next.memories.iter().flatten()).any(|v| !v.is_finite()) {
            return Err(invalid("exact acoustic or wall-memory endpoint is unrepresentable").into());
        }
        let (initial_energy, initial_fluid_roundoff, initial_memory_roundoff) = self.energy_with_roundoff(initial, cancelled)?;
        let (final_energy, final_fluid_roundoff, final_memory_roundoff) = self.energy_with_roundoff(&next, cancelled)?;
        let acoustic_change = energy_change(&self.model.mass, &initial.acoustic.pressure, &next.acoustic.pressure, cancelled)?
            + energy_change(&self.model.stiffness, &initial.acoustic.potential, &next.acoustic.potential, cancelled)?;
        let mut memory_change = Dd::from_f64(0.0);
        for ((arm, before), after) in self.arms.iter().zip(&initial.memories).zip(&next.memories) {
            memory_change = memory_change + Dd::from_f64(arm.tau)*Dd::from_f64(energy_change(&arm.conductance, before, after, cancelled)?);
        }
        let mut input = Dd::from_f64(0.0);
        let mut boundary_work = Dd::from_f64(0.0);
        let mut volume_work = Dd::from_f64(0.0);
        let mut absolute_work = 0.0;
        for (i, (&amplitude, &integral)) in amplitudes.iter().zip(&result.integrals[3..]).enumerate() {
            let work = Dd::from_f64(amplitude)*Dd::from_f64(integral);
            input = input + work; absolute_work += work.to_f64().abs();
            if i >= self.boundary_port_count { volume_work = volume_work + work; }
            else if i >= self.linear_port_count { boundary_work = boundary_work + work; }
        }
        let work = ExactRelaxingTimeWork {
            initial_energy_j: initial_energy.total_j(), final_energy_j: final_energy.total_j(),
            delta_energy_j: (Dd::from_f64(acoustic_change)+memory_change).to_f64(),
            acoustic_delta_energy_j: acoustic_change, memory_delta_energy_j: memory_change.to_f64(),
            input_work_j: input.to_f64(), boundary_source_work_j: boundary_work.to_f64(), volume_source_work_j: volume_work.to_f64(),
            instantaneous_dissipation_j: result.integrals[0], relaxation_dissipation_j: result.integrals[1], memory_work_j: result.integrals[2],
        };
        if result.integrals.iter().any(|v| !v.is_finite()) || !work.boundary_source_work_j.is_finite() || !work.volume_source_work_j.is_finite() {
            return Err(invalid("exact wall-memory source work or trajectory integral overflow").into());
        }
        let scale = work.initial_energy_j.abs().max(work.final_energy_j.abs()).max(work.delta_energy_j.abs())
            .max(absolute_work).max(work.instantaneous_dissipation_j.abs())
            .max(work.relaxation_dissipation_j.abs()).max(work.memory_work_j.abs());
        let roundoff = initial_fluid_roundoff+final_fluid_roundoff+initial_memory_roundoff+final_memory_roundoff
            +64.0*f64::EPSILON*absolute_work;
        let tolerance = self.limits.energy_balance_tolerance;
        balance(work.acoustic_balance_residual_j(), scale, roundoff, tolerance)?;
        balance(work.memory_balance_residual_j(), scale, roundoff, tolerance)?;
        balance(work.balance_residual_j(), scale, roundoff, tolerance)?;
        let limit = tolerance*scale+roundoff;
        if work.instantaneous_dissipation_j < -limit || work.relaxation_dissipation_j < -limit {
            return Err(HeterogeneousError::EnergyBalance {
                observed_j: work.instantaneous_dissipation_j.min(work.relaxation_dissipation_j), limit_j: limit,
            });
        }
        poll(cancelled)?;
        Ok(ExactRelaxingTimeStep {
            state: next, work,
            integrated_port_pressure_pa_m2_s: result.integrals[3..3+self.boundary_port_count].to_vec(),
            integrated_volume_source_pressure_pa_s: result.integrals[3+self.boundary_port_count..].to_vec(),
        })
    }
}

#[cfg(test)]
mod tests;
