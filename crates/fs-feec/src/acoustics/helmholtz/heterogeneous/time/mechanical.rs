//! Causal locally reacting walls in the original full acoustic trace space.
//!
//! The original [`WallPin`] supplies surface mass `sigma` [kg/m²], stiffness
//! `K` [Pa/m] and resistance `r` [Pa s/m]. Positive displacement `w` and
//! velocity `v` point OUT of the fluid. Optional original Maxwell arms retain
//! energy-normalized memory `z_j`, with `e_j = sqrt(k_j) w - z_j`:
//! `sigma v_dot + r v + K w + sum sqrt(k_j) e_j = p + traction`,
//! `tau_j z_dot_j = e_j`. Each arm stores `e_j^T W e_j/2` and dissipates
//! `e_j^T W e_j/tau_j`; it is not an instantaneous fitted resistance.
//! Every active original boundary coordinate retains both states. With the
//! original face mass `W = integral(phi_i phi_j)`, wall energy is
//! `(sigma v^T W v + K w^T W w)/2`, and fluid feedback is `+W v`.
//! No patch-average piston, diagonal mass or extra area factor is substituted.
//!
//! One midpoint elimination gives
//! `a_j = tau_j/(tau_j+h/2)` condenses each original Maxwell midpoint law.
//! Its force is `a_j sqrt(k_j) (sqrt(k_j)(w_0+h v_mid/2)-z_j0)`.
//! The positive effective admittance enters the ORIGINAL sparse acoustic
//! factor; the history enters its original projected load. The same solve
//! reconstructs all wall states, then checks the uncondensed equations and
//! physical work. This effective admittance is never reported as wall heat.
//!
//! Patches have constant coefficients and disjoint original faces, including
//! disjointness from real Robin and pressure-release patches. Adjacent mechanical patches retain
//! independent states at their shared edges. Imposed velocity may share a
//! wall face. This is the declared local wall law, not shell bending, lateral
//! elastic waves, a material fit, or a moving fluid mesh. Relaxation coefficients
//! and their material validity remain caller-owned. Midpoint is passive but
//! still needs temporal refinement to resolve the declared relaxation times.

use super::super::boundary_admittance::relaxation::midpoint_weights;
use super::super::{
    BoundaryVelocitySource, SpatialRelaxingBoundary, VolumeVelocitySource, source::admit_count,
};
use super::{
    BoundaryAdmittance, C64, Coo, Csr, Dd, HeterogeneousError, HeterogeneousHelmholtzFe,
    HeterogeneousTimeEnergy, HeterogeneousTimeLimits, HeterogeneousTimeState,
    PreparedHeterogeneousTime, RealVelocityTrace, ResistiveBoundary, SpatialResistiveBoundary,
    energy_change, invalid, poll, positive_form, relative,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use fs_phs::acoustic_outlet::AcousticWallRelaxation;
use fs_phs::{RelaxationBranch, WallPin};
use std::collections::BTreeSet;

/// Complete-field exponential evolution with the original causal wall states.
pub mod exact;

/// Applied physical forces on original distributed wall coordinates.
pub mod traction;
use traction::PreparedMechanicalTraction;

mod relaxing;

/// Original surface-specific mechanical law on a complete boundary patch.
#[derive(Clone, Debug)]
pub struct MechanicalBoundary {
    pub faces: Vec<[u32; 3]>,
    pub wall: WallPin,
    /// Complete positive Maxwell spectrum, in the same surface units as the
    /// original acoustic wall owner. Empty means the original Kelvin wall.
    pub relaxation: Vec<AcousticWallRelaxation>,
}

/// Existing acoustic admission plus an explicit full wall-state extent.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalTimeLimits {
    pub acoustic: HeterogeneousTimeLimits,
    /// Sum of active boundary DOFs across patches; each has displacement and
    /// velocity. Patch CSR and O(ndof) work vectors are additional to factor bytes.
    pub max_wall_dofs: usize,
    /// Sum of complete mechanical Maxwell and spatial-admittance memory DOFs.
    /// Maxwell memories use energy-normalized coordinates; admittance memories
    /// use pressure [Pa]. Release-shared and zero-residue coordinates count.
    pub max_wall_relaxation_dofs: usize,
}

impl Default for MechanicalTimeLimits {
    fn default() -> Self {
        Self {
            acoustic: HeterogeneousTimeLimits::default(),
            max_wall_dofs: 4096,
            max_wall_relaxation_dofs: 16384,
        }
    }
}

struct PreparedWall {
    boundary: MechanicalBoundary,
    dofs: Vec<usize>,
    mass: Csr,
    admittance: f64,
    inertia_weight: f64,
    spring_weight: f64,
    relaxation: Vec<PreparedRelaxation>,
}

struct PreparedRelaxation {
    branch: RelaxationBranch,
    coupling: f64,
    history_weight: f64,
}

/// One original acoustic factor with every declared distributed wall state.
pub struct PreparedMechanicalTime<'m, 'f, 'c> {
    condensed: PreparedHeterogeneousTime<'m, 'f, 'c>,
    conductance: Csr,
    walls: Vec<PreparedWall>,
    admittance: Vec<relaxing::TraceMemory>,
    admittance_weights: Vec<(f64, f64)>,
}

/// Candidate fields bound to this exact preparation, including its walls.
#[derive(Clone)]
pub struct MechanicalTimeState<'p, 'm, 'f, 'c> {
    owner: &'p PreparedMechanicalTime<'m, 'f, 'c>,
    acoustic: HeterogeneousTimeState<'m, 'f, 'c>,
    displacements: Vec<Vec<f64>>,
    velocities: Vec<Vec<f64>>,
    // Patch, original arm, original hierarchical trace coordinate.
    memories: Vec<Vec<Vec<f64>>>,
    admittance_memories: Vec<Vec<f64>>,
}

impl core::fmt::Debug for MechanicalTimeState<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("MechanicalTimeState")
            .field("acoustic", &self.acoustic)
            .field("wall_patches", &self.displacements.len())
            .finish_non_exhaustive()
    }
}

impl MechanicalTimeState<'_, '_, '_, '_> {
    /// Borrow the original acoustic state for the existing point and surface
    /// probes. Its spatial/material identity and full hierarchical fields are
    /// preserved; observing it cannot advance or modify the coupled walls.
    #[must_use]
    pub fn acoustic_state(&self) -> &HeterogeneousTimeState<'_, '_, '_> {
        &self.acoustic
    }

    /// Original hierarchical pressure-potential coefficients [Pa s].
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        self.acoustic.pressure_potential_pa_s()
    }
    /// Original hierarchical pressure coefficients [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        self.acoustic.pressure_pa()
    }
    /// Displacement [m] in the patch's original `wall_dofs` order.
    #[must_use]
    pub fn wall_displacement_m(&self, wall: usize) -> Option<&[f64]> {
        self.displacements.get(wall).map(Vec::as_slice)
    }
    /// Outward normal velocity [m/s] in the same original trace order.
    #[must_use]
    pub fn wall_normal_velocity_m_s(&self, wall: usize) -> Option<&[f64]> {
        self.velocities.get(wall).map(Vec::as_slice)
    }
    /// Original energy-normalized Maxwell memory, `z = sqrt(k) w` at rest
    /// after complete relaxation. The arm's stored energy per area is
    /// `(sqrt(k) w-z)^2/2`; this is not a displacement or pressure sample.
    #[must_use]
    pub fn wall_relaxation_memory(&self, wall: usize, arm: usize) -> Option<&[f64]> {
        self.memories.get(wall)?.get(arm).map(Vec::as_slice)
    }

    /// Original spatial-admittance memory q [Pa], in its geometric trace order.
    #[must_use]
    pub fn admittance_memory_pressure_pa(&self, arm: usize) -> Option<&[f64]> {
        self.admittance_memories.get(arm).map(Vec::as_slice)
    }
}

/// Independently evaluated physical acoustic and wall storage [J].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechanicalTimeEnergy {
    pub acoustic: HeterogeneousTimeEnergy,
    pub wall_kinetic_j: f64,
    pub wall_elastic_j: f64,
    pub wall_relaxation_j: f64,
    /// Spatial-admittance storage sum tau_j q_j^T Cj q_j/2.
    pub admittance_memory_j: f64,
}

impl MechanicalTimeEnergy {
    #[must_use]
    pub fn total_j(self) -> f64 {
        self.acoustic.total_j()
            + self.wall_kinetic_j
            + self.wall_elastic_j
            + self.wall_relaxation_j
            + self.admittance_memory_j
    }
}

/// Original source work, reversible exchange and irreversible losses [J].
#[derive(Clone, Copy, Debug)]
pub struct MechanicalTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    pub acoustic_delta_energy_j: f64,
    pub wall_delta_energy_j: f64,
    pub admittance_delta_energy_j: f64,
    /// Total signed fluid-source AND applied wall-traction work, counted once.
    pub input_work_j: f64,
    /// Independently accumulated original fluid-source work. Retained directly
    /// so weak acoustic work survives a much larger mechanical input.
    pub acoustic_input_work_j: f64,
    /// Signed prepared boundary-source subtotal already in `input_work_j`.
    pub boundary_source_work_j: f64,
    /// Signed volume-source subtotal already included in `input_work_j`.
    pub volume_source_work_j: f64,
    /// Signed external wall-traction work, already included in input_work_j.
    /// This is distinct from fluid-to-wall exchange and can be negative.
    pub wall_traction_work_j: f64,
    /// Dissipation of the original REAL Robin patches only.
    pub boundary_dissipation_j: f64,
    /// Original Kelvin resistance PLUS all original Maxwell arm dissipation,
    /// never effective-admittance loss.
    pub wall_dissipation_j: f64,
    /// Maxwell subtotal already included in `wall_dissipation_j`.
    pub wall_relaxation_dissipation_j: f64,
    /// Signed fluid-to-wall exchange `h sum p_mid^T W v_mid`.
    pub wall_work_j: f64,
    /// Original spatial-admittance memory heat, separate from mechanical heat.
    pub admittance_dissipation_j: f64,
    /// Signed fluid-to-admittance transfer h sum p_mid^T Cj qj_mid.
    pub admittance_work_j: f64,
    /// Original uncondensed fluid and wall equations, not the stored factor.
    pub backward_error: f64,
}

impl MechanicalTimeWork {
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j
            + self.boundary_dissipation_j
            + self.wall_dissipation_j
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
        self.wall_delta_energy_j + self.wall_dissipation_j
            - self.wall_work_j
            - self.wall_traction_work_j
    }
    #[must_use]
    pub fn admittance_balance_residual_j(self) -> f64 {
        self.admittance_delta_energy_j + self.admittance_dissipation_j - self.admittance_work_j
    }
}

/// Fully checked candidate; the caller owns acceptance and the clock.
#[derive(Clone, Debug)]
pub struct MechanicalTimeStep<'p, 'm, 'f, 'c> {
    pub state: MechanicalTimeState<'p, 'm, 'f, 'c>,
    pub work: MechanicalTimeWork,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare reciprocal full-trace mechanical walls on the original time owner.
    /// All wall/Robin faces must be disjoint. No material data is inferred.
    pub fn prepare_time_with_walls<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        limits: MechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_time_with_walls_and_pressure_release(
            duration_s,
            resistive,
            boundaries,
            &[],
            limits,
            cancelled,
        )
    }

    /// Prepare causal distributed walls together with exact pressure-release faces.
    /// Release faces impose zero pressure and zero pressure-potential trace on
    /// every original hierarchical trace DOF. Imported states must satisfy both;
    /// a stationary flow potential is never silently projected away.
    ///
    /// Mechanical, real Robin and pressure-release faces must be disjoint.
    /// Adjacent patches may share edges; every mechanical displacement, velocity
    /// and Maxwell memory coordinate is retained independently of the fluid
    /// constraint. Imposed velocity may share a mechanical face but cannot share
    /// a pressure-release face. The acoustic factor budget counts free fluid
    /// DOFs; the original wall and memory budgets still count their complete traces.
    pub fn prepare_time_with_walls_and_pressure_release<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        limits: MechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_time_with_walls_and_spatial_boundaries(
            duration_s,
            resistive,
            &[],
            boundaries,
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Couple complete causal mechanical walls to spatial resistive patches.
    ///
    /// Spatial conductance retains its original pressure trace and contributes
    /// only to boundary heat. Mechanical resistance, relaxation and reversible
    /// fluid/wall work retain their independent physical accounts. Uniform,
    /// spatial, mechanical and pressure-release patches must be disjoint;
    /// velocity sources may share any non-release face. No wall or pressure
    /// coordinate is averaged or discarded.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_walls_and_spatial_boundaries<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        spatial_boundaries: &[&SpatialResistiveBoundary<'_, '_, '_>],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        limits: MechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_time_with_walls_and_relaxing_boundaries(
            duration_s,
            resistive,
            spatial_boundaries,
            boundaries,
            &[],
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Evolve mechanical/Maxwell walls and spatial causal admittance memories
    /// against ONE original acoustic solve. Every physical patch must have
    /// disjoint faces, while independent histories at shared edges remain live.
    /// The memory limit counts both kinds of original boundary history.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_walls_and_relaxing_boundaries<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        spatial_boundaries: &[&SpatialResistiveBoundary<'_, '_, '_>],
        boundaries: &[MechanicalBoundary],
        relaxing: &[&SpatialRelaxingBoundary<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: MechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 || limits.max_wall_dofs == 0 {
            return Err(
                invalid("positive mechanical time step and wall DOF limit required").into(),
            );
        }
        if pressure_release_faces.is_empty() && self.ndof() > limits.acoustic.max_dofs {
            return Err(HeterogeneousError::Limit {
                what: "transient pressure degrees of freedom",
                requested: self.ndof(),
                maximum: limits.acoustic.max_dofs,
            });
        }
        let mut seen = BTreeSet::new();
        for &face in pressure_release_faces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if !seen.insert(key) {
                return Err(invalid("pressure-release faces must be unique").into());
            }
        }
        for faces in resistive
            .iter()
            .map(|b| &b.faces)
            .chain(spatial_boundaries.iter().map(|b| &b.faces))
            .chain(boundaries.iter().map(|b| &b.faces))
            .chain(relaxing.iter().map(|b| &b.instantaneous.faces))
        {
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid(
                        "mechanical, spatial admittance, real Robin and pressure-release faces must be disjoint",
                    )
                    .into());
                }
            }
        }
        let pressure_release = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let physical_fields = spatial_boundaries
            .iter()
            .copied()
            .chain(relaxing.iter().map(|law| &law.instantaneous))
            .collect::<Vec<_>>();
        let conductance =
            self.resistive_matrix(resistive, &physical_fields, &pressure_release, cancelled)?;
        let mut total = 0usize;
        let mut total_memory = 0usize;
        let mut walls = Vec::with_capacity(boundaries.len());
        let mut effective = resistive.to_vec();
        for boundary in boundaries {
            poll(cancelled)?;
            // Reuse the original WallPin admission, not a fitted frequency sample.
            fs_phs::wall_specific_impedance(&boundary.wall, 1.0)
                .map_err(HeterogeneousError::Wall)?;
            let pin = boundary.wall;
            let inertia = 2.0 * (pin.surface_density / duration_s);
            let mut relaxation = Vec::with_capacity(boundary.relaxation.len());
            let mut relaxed_midpoint_stiffness = pin.stiffness_per_area;
            for arm in &boundary.relaxation {
                poll(cancelled)?;
                let branch = RelaxationBranch {
                    projection: vec![1.0],
                    stiffness: arm.stiffness_per_area,
                    relaxation_time_s: arm.relaxation_time_s,
                };
                let coupling = branch
                    .scalar_relaxed_memory(1.0)
                    .map_err(HeterogeneousError::Wall)?;
                // Same bounded weights as the original scalar midpoint owner.
                let scale = arm.relaxation_time_s.max(0.5 * duration_s);
                let tau = arm.relaxation_time_s / scale;
                let half_h = (0.5 * duration_s) / scale;
                let a = tau / (tau + half_h);
                relaxed_midpoint_stiffness += a * arm.stiffness_per_area;
                relaxation.push(PreparedRelaxation {
                    branch,
                    coupling,
                    history_weight: a * coupling,
                });
            }
            let impedance =
                inertia + pin.resistance + (0.5 * duration_s) * relaxed_midpoint_stiffness;
            let admittance = 1.0 / impedance;
            let inertia_weight = inertia * admittance;
            let spring_weight = pin.stiffness_per_area * admittance;
            if ![impedance, admittance, inertia_weight, spring_weight]
                .iter()
                .all(|v| v.is_finite() && *v > 0.0)
            {
                return Err(
                    invalid("mechanical wall midpoint coefficients are unrepresentable").into(),
                );
            }
            for arm in &mut relaxation {
                arm.history_weight *= admittance;
                if !arm.history_weight.is_finite() || arm.history_weight <= 0.0 {
                    return Err(invalid("wall memory coupling is unrepresentable").into());
                }
            }
            // This is the SAME symmetric original face quadrature as Robin loading.
            let (global, _) = self.admittance_matrices(
                &[BoundaryAdmittance {
                    faces: boundary.faces.clone(),
                    admittance_m_pa_s: C64::ONE,
                }],
                cancelled,
            )?;
            let dofs = (0..self.ndof())
                .filter(|&row| {
                    let (columns, values) = global.row(row);
                    columns.binary_search(&row).is_ok_and(|i| values[i] > 0.0)
                })
                .collect::<Vec<_>>();
            total = total
                .checked_add(dofs.len())
                .ok_or_else(|| invalid("mechanical trace extent overflow"))?;
            if total > limits.max_wall_dofs {
                return Err(HeterogeneousError::Limit {
                    what: "distributed wall trace degrees of freedom",
                    requested: total,
                    maximum: limits.max_wall_dofs,
                });
            }
            total_memory = dofs
                .len()
                .checked_mul(relaxation.len())
                .and_then(|count| total_memory.checked_add(count))
                .ok_or_else(|| invalid("mechanical wall memory extent overflow"))?;
            if total_memory > limits.max_wall_relaxation_dofs {
                return Err(HeterogeneousError::Limit {
                    what: "distributed wall relaxation degrees of freedom",
                    requested: total_memory,
                    maximum: limits.max_wall_relaxation_dofs,
                });
            }
            let mut local = Coo::new(dofs.len(), dofs.len());
            for (row, &global_row) in dofs.iter().enumerate() {
                poll(cancelled)?;
                let (cols, entries) = global.row(global_row);
                for (&col, &entry) in cols.iter().zip(entries) {
                    if entry != 0.0 {
                        let col = dofs
                            .binary_search(&col)
                            .map_err(|_| invalid("incomplete original wall trace support"))?;
                        local.push(row, col, entry);
                    }
                }
            }
            effective.push(ResistiveBoundary {
                faces: boundary.faces.clone(),
                admittance_m_pa_s: admittance,
            });
            walls.push(PreparedWall {
                boundary: boundary.clone(),
                dofs,
                mass: local.assemble(),
                admittance,
                inertia_weight,
                spring_weight,
                relaxation,
            });
        }
        let admittance = relaxing::prepare(
            self,
            relaxing,
            limits.max_wall_relaxation_dofs - total_memory,
            cancelled,
        )?;
        let admittance_weights = admittance
            .iter()
            .map(|arm| midpoint_weights(duration_s, arm.tau))
            .collect::<Result<Vec<_>, _>>()?;
        let mut condensed_admittance = Vec::with_capacity(relaxing.len());
        for law in relaxing {
            poll(cancelled)?;
            condensed_admittance.push(law.midpoint_condensate(duration_s, cancelled)?);
        }
        let effective_fields = spatial_boundaries
            .iter()
            .copied()
            .chain(condensed_admittance.iter())
            .collect::<Vec<_>>();
        let condensed = self.prepare_time_with_spatial_boundaries(
            duration_s,
            &effective,
            &effective_fields,
            pressure_release_faces,
            limits.acoustic,
            cancelled,
        )?;
        Ok(PreparedMechanicalTime {
            condensed,
            conductance,
            walls,
            admittance,
            admittance_weights,
        })
    }
}

impl<'m, 'f, 'c> PreparedMechanicalTime<'m, 'f, 'c> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.condensed.duration_s()
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
        self.admittance.get(arm).map(|a| a.dofs.as_slice())
    }
    #[must_use]
    pub fn admittance_memory_patch_and_term(&self, arm: usize) -> Option<(usize, usize)> {
        self.admittance.get(arm).map(|a| (a.patch, a.term))
    }
    #[must_use]
    pub fn admittance_relaxation_time_s(&self, arm: usize) -> Option<f64> {
        self.admittance.get(arm).map(|a| a.tau)
    }
    #[must_use]
    pub fn zero_state(&self) -> MechanicalTimeState<'_, 'm, 'f, 'c> {
        MechanicalTimeState {
            owner: self,
            acoustic: self.condensed.zero_state(),
            displacements: self.walls.iter().map(|w| vec![0.0; w.dofs.len()]).collect(),
            velocities: self.walls.iter().map(|w| vec![0.0; w.dofs.len()]).collect(),
            memories: self
                .walls
                .iter()
                .map(|w| vec![vec![0.0; w.dofs.len()]; w.relaxation.len()])
                .collect(),
            admittance_memories: self
                .admittance
                .iter()
                .map(|a| vec![0.0; a.dofs.len()])
                .collect(),
        }
    }
    /// Import the original acoustic field; every mechanical and admittance
    /// memory starts unexcited. No equilibrium memory is inferred from pressure.
    pub fn state(
        &self,
        potential: &[f64],
        pressure: &[f64],
    ) -> Result<MechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        let acoustic = self.condensed.state(potential, pressure)?;
        let mut state = self.zero_state();
        state.acoustic = acoustic;
        Ok(state)
    }
    /// Import a nonhereditary wall in `wall_dofs` order. A wall with Maxwell
    /// arms requires `state_with_wall_memory`; its past is never inferred.
    pub fn state_with_walls(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
    ) -> Result<MechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        if self.walls.iter().any(|w| !w.relaxation.is_empty()) {
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
    /// Import every original physical history. Memory is indexed by patch,
    /// original arm, then `wall_dofs`, and uses `RelaxationBranch`'s existing
    /// energy-normalized coordinate. Even a zero-displacement wall may retain
    /// prestress in its memory, so omitted histories are refused.
    pub fn state_with_wall_memory(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
        memories: &[Vec<Vec<f64>>],
    ) -> Result<MechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        if !self.admittance.is_empty() {
            return Err(invalid("complete original spatial admittance history is required").into());
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

    /// Import every mechanical and spatial-admittance history explicitly.
    /// Admittance q [Pa] is indexed by flattened patch/term then original trace
    /// DOF. Its history is independent of Maxwell z and of fluid constraints.
    #[allow(clippy::too_many_arguments)]
    pub fn state_with_all_memories(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
        memories: &[Vec<Vec<f64>>],
        admittance_memories: &[Vec<f64>],
    ) -> Result<MechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
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
                    || memories[index].len() != wall.relaxation.len()
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
            || admittance_memories
                .iter()
                .zip(&self.admittance)
                .any(|(q, arm)| q.len() != arm.dofs.len() || q.iter().any(|v| !v.is_finite()))
        {
            return Err(
                invalid("finite complete original spatial admittance memories required").into(),
            );
        }
        let state = MechanicalTimeState {
            owner: self,
            acoustic: self.condensed.state(potential, pressure)?,
            displacements: displacements_m.to_vec(),
            velocities: velocities_m_s.to_vec(),
            memories: memories.to_vec(),
            admittance_memories: admittance_memories.to_vec(),
        };
        self.energy(&state)?;
        Ok(state)
    }
    fn validate_state(
        &self,
        state: &MechanicalTimeState<'_, '_, '_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self, state.owner) {
            return Err(
                invalid("mechanical state belongs to another prepared wall configuration").into(),
            );
        }
        Ok(())
    }
    pub fn energy(
        &self,
        state: &MechanicalTimeState<'_, '_, '_, '_>,
    ) -> Result<MechanicalTimeEnergy, HeterogeneousError> {
        self.validate_state(state)?;
        self.energy_with_roundoff(state, &mut || false).map(|v| v.0)
    }
    fn energy_with_roundoff(
        &self,
        state: &MechanicalTimeState<'_, '_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(MechanicalTimeEnergy, f64, f64, f64), HeterogeneousError> {
        let (acoustic, acoustic_roundoff) = self
            .condensed
            .energy_with_roundoff(&state.acoustic, cancelled)?;
        let mut energy = MechanicalTimeEnergy {
            acoustic,
            wall_kinetic_j: 0.0,
            wall_elastic_j: 0.0,
            wall_relaxation_j: 0.0,
            admittance_memory_j: 0.0,
        };
        let mut roundoff = 0.0;
        for (index, wall) in self.walls.iter().enumerate() {
            let (kinetic, rk) = positive_form(&wall.mass, &state.velocities[index], cancelled)?;
            let (elastic, re) = positive_form(&wall.mass, &state.displacements[index], cancelled)?;
            let pin = wall.boundary.wall;
            energy.wall_kinetic_j += 0.5 * pin.surface_density * kinetic;
            energy.wall_elastic_j += 0.5 * pin.stiffness_per_area * elastic;
            roundoff += 0.5 * (pin.surface_density * rk + pin.stiffness_per_area * re);
            for (arm, memory) in wall.relaxation.iter().zip(&state.memories[index]) {
                poll(cancelled)?;
                let extension = state.displacements[index]
                    .iter()
                    .zip(memory)
                    .map(|(&w, &z)| arm.coupling.mul_add(w, -z))
                    .collect::<Vec<_>>();
                let (stored, error) = positive_form(&wall.mass, &extension, cancelled)?;
                energy.wall_relaxation_j += 0.5 * stored;
                roundoff += 0.5 * error;
            }
        }
        let (admittance_energy, admittance_roundoff) =
            relaxing::energy(&self.admittance, &state.admittance_memories, cancelled)?;
        energy.admittance_memory_j = admittance_energy;
        if !energy.total_j().is_finite() || !roundoff.is_finite() {
            return Err(invalid("mechanical acoustic storage is unrepresentable").into());
        }
        Ok((energy, acoustic_roundoff, roundoff, admittance_roundoff))
    }

    /// Solve fluid and every wall together; no participant is advanced explicitly.
    /// Sources are original inward normal velocities at the step midpoint.
    pub fn step<'p>(
        &'p self,
        initial: &MechanicalTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_volume_sources(initial, inward_velocity, &[], cancelled)
    }

    /// Advance the complete fluid/wall state with original distributed sources.
    /// Boundary velocities [m/s] and volume flows [m^3/s] share the same
    /// midpoint pressure and sparse solve. Sources retain their frozen spatial
    /// loads; no material/source callback, extra factor or wall state is added.
    /// Negative flows extract volume and may do negative signed work.
    pub fn step_with_volume_sources<'p>(
        &'p self,
        initial: &MechanicalTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        sources: &[(&VolumeVelocitySource<'_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_boundary_sources(initial, inward_velocity, &[], sources, cancelled)
    }

    /// Compose prepared quadratic boundary motion, original linear traces and
    /// volume injection with every distributed wall and memory coordinate.
    /// Boundary amplitudes are midpoint velocities [m/s]; volume amplitudes
    /// remain [m^3/s]. Each spatial source is projected once at preparation.
    pub fn step_with_boundary_sources<'p>(
        &'p self,
        initial: &MechanicalTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, f64)],
        sources: &[(&VolumeVelocitySource<'_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_wall_tractions(
            initial,
            inward_velocity,
            boundary_sources,
            sources,
            &[],
            cancelled,
        )
    }

    /// Evolve acoustic forcing AND applied mechanical traction on one coupled
    /// midpoint trajectory. Traction amplitudes [Pa] act outward on their
    /// original wall profiles; pressure feedback still determines the motion.
    /// Repeated/coherent shapes sum before solving. Their signed work is kept
    /// separate from fluid sources and reversible fluid-to-wall transfer.
    /// All amplitudes are midpoint values; their temporal resolution remains
    /// caller-owned. Cancellation or a failed physical balance publishes no state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_wall_tractions<'p>(
        &'p self,
        initial: &MechanicalTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, f64)],
        sources: &[(&VolumeVelocitySource<'_, '_, '_>, f64)],
        tractions: &[(&PreparedMechanicalTraction<'_, '_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        self.validate_state(initial)?;
        let core = &self.condensed;
        core.admit_velocity_faces(inward_velocity, cancelled)?;
        core.admit_boundary_sources(boundary_sources, cancelled)?;
        admit_count(sources.len())?;
        for &(source, amplitude) in sources {
            poll(cancelled)?;
            source.validate(core.model)?;
            if !amplitude.is_finite() {
                return Err(invalid("finite midpoint volume-flow amplitude required").into());
            }
        }
        admit_count(tractions.len())?;
        let mut traction_sums = self
            .walls
            .iter()
            .map(|wall| vec![Dd::ZERO; wall.dofs.len()])
            .collect::<Vec<_>>();
        for &(traction, amplitude) in tractions {
            poll(cancelled)?;
            if !core::ptr::eq(traction.owner, self) || !amplitude.is_finite() {
                return Err(invalid(
                    "traction needs its original wall owner and finite Pa amplitude",
                )
                .into());
            }
            for (i, (sum, &coefficient)) in traction_sums[traction.wall]
                .iter_mut()
                .zip(&traction.coefficients)
                .enumerate()
            {
                if i.is_multiple_of(256) {
                    poll(cancelled)?;
                }
                *sum = *sum + Dd::from_f64(coefficient) * Dd::from_f64(amplitude);
            }
        }
        let traction_fields = traction_sums
            .into_iter()
            .map(|field| field.into_iter().map(Dd::to_f64).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        if traction_fields.iter().flatten().any(|p| !p.is_finite()) {
            return Err(invalid("coherent mechanical traction overflow").into());
        }
        let h = self.duration_s();
        let mut source = core.model.geometry.boundary_linear_load(inward_velocity)?;
        if !sources.is_empty() || !boundary_sources.is_empty() {
            for (row, value) in source.iter_mut().enumerate() {
                poll(cancelled)?;
                let mut sum = Dd::from_f64(*value);
                for &(boundary, amplitude) in boundary_sources {
                    sum = sum + Dd::from_f64(boundary.load[row]) * Dd::from_f64(amplitude);
                }
                for &(volume, amplitude) in sources {
                    sum = sum + Dd::from_f64(volume.load[row]) * Dd::from_f64(amplitude);
                }
                *value = sum.to_f64();
                if !value.is_finite() {
                    return Err(invalid("mechanical acoustic source load overflow").into());
                }
            }
        }
        let mut load = source.clone();
        for (index, wall) in self.walls.iter().enumerate() {
            poll(cancelled)?;
            let history = initial.velocities[index]
                .iter()
                .zip(&initial.displacements[index])
                .enumerate()
                .map(|(i, (&v, &w))| {
                    wall.admittance.mul_add(
                        traction_fields[index][i],
                        wall.history(v, w, &initial.memories[index], i),
                    )
                })
                .collect::<Vec<_>>();
            let mut projected = vec![0.0; wall.dofs.len()];
            wall.mass.spmv(&history, &mut projected);
            for (&dof, value) in wall.dofs.iter().zip(projected) {
                load[dof] -= value;
            }
        }
        relaxing::subtract_history(
            &self.admittance,
            &self.admittance_weights,
            &initial.admittance_memories,
            &mut load,
            cancelled,
        )?;
        // Only this existing owner factors/solves. Its condensed Robin work is
        // an algebraic check, not an observable physical dissipation account.
        let acoustic = core.step_projected(&initial.acoustic, &load, cancelled)?;
        let midpoint_p = initial
            .pressure_pa()
            .iter()
            .zip(acoustic.state.pressure_pa())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let midpoint_psi = initial
            .pressure_potential_pa_s()
            .iter()
            .zip(acoustic.state.pressure_potential_pa_s())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let mut next = MechanicalTimeState {
            owner: self,
            acoustic: acoustic.state,
            displacements: Vec::with_capacity(self.walls.len()),
            velocities: Vec::with_capacity(self.walls.len()),
            memories: Vec::with_capacity(self.walls.len()),
            admittance_memories: Vec::with_capacity(self.admittance.len()),
        };
        let mut wall_load = vec![0.0; source.len()];
        let mut wall_change = 0.0;
        let mut wall_loss = 0.0;
        let mut relaxation_loss = 0.0;
        let mut wall_loss_roundoff = 0.0;
        let mut wall_work = 0.0;
        let mut wall_work_roundoff = 0.0;
        let mut traction_work = Dd::ZERO;
        let mut traction_roundoff = 0.0;
        let mut backward_error: f64 = 0.0;
        for (index, wall) in self.walls.iter().enumerate() {
            poll(cancelled)?;
            let pin = wall.boundary.wall;
            let previous_w = &initial.displacements[index];
            let previous_v = &initial.velocities[index];
            let mut displacement = Vec::with_capacity(wall.dofs.len());
            let mut velocity = Vec::with_capacity(wall.dofs.len());
            let mut middle_v = Vec::with_capacity(wall.dofs.len());
            let mut pressure = Vec::with_capacity(wall.dofs.len());
            let mut memories = vec![Vec::with_capacity(wall.dofs.len()); wall.relaxation.len()];
            let mut middle_extensions =
                vec![Vec::with_capacity(wall.dofs.len()); wall.relaxation.len()];
            for (i, &dof) in wall.dofs.iter().enumerate() {
                let p = midpoint_p[dof];
                let history =
                    wall.history(previous_v[i], previous_w[i], &initial.memories[index], i);
                let applied = traction_fields[index][i];
                let total_pressure = (Dd::from_f64(p) + Dd::from_f64(applied)).to_f64();
                let v = wall.admittance.mul_add(total_pressure, history);
                let end_v = v + (v - previous_v[i]);
                let end_w = h.mul_add(v, previous_w[i]);
                if !end_v.is_finite() || !end_w.is_finite() {
                    return Err(invalid("mechanical wall candidate overflow").into());
                }
                let mid_v = f64::midpoint(previous_v[i], end_v);
                let mid_w = f64::midpoint(previous_w[i], end_w);
                let momentum = Dd::from_f64(pin.surface_density)
                    * (Dd::from_f64(end_v) - Dd::from_f64(previous_v[i]));
                let drag = Dd::from_f64(h) * Dd::from_f64(pin.resistance) * Dd::from_f64(mid_v);
                let spring =
                    Dd::from_f64(h) * Dd::from_f64(pin.stiffness_per_area) * Dd::from_f64(mid_w);
                let force = Dd::from_f64(h) * Dd::from_f64(p);
                let applied_force = Dd::from_f64(h) * Dd::from_f64(applied);
                let mut hereditary = Dd::from_f64(0.0);
                let mut hereditary_scale = 0.0;
                for (j, arm) in wall.relaxation.iter().enumerate() {
                    poll(cancelled)?;
                    let z0 = initial.memories[index][j][i];
                    let (z1, arm_force, _, _) = arm
                        .branch
                        .scalar_midpoint_step(previous_w[i], end_w, z0, h)
                        .map_err(HeterogeneousError::Wall)?;
                    let extension = Dd::from_f64(arm.coupling) * Dd::from_f64(mid_w)
                        - (Dd::from_f64(z0) + Dd::from_f64(z1)) * Dd::from_f64(0.5);
                    let change = Dd::from_f64(arm.branch.relaxation_time_s)
                        * (Dd::from_f64(z1) - Dd::from_f64(z0));
                    let relaxation = Dd::from_f64(h) * extension;
                    let memory_scale = arm.branch.relaxation_time_s * (z0.abs() + z1.abs())
                        + h * (arm.coupling * mid_w.abs() + f64::midpoint(z0, z1).abs());
                    if !memory_scale.is_finite() || !extension.to_f64().is_finite() {
                        return Err(invalid("mechanical memory residual overflow").into());
                    }
                    backward_error = backward_error
                        .max(relative((change - relaxation).to_f64().abs(), memory_scale));
                    let term = Dd::from_f64(h) * Dd::from_f64(arm.coupling) * extension;
                    hereditary = hereditary + term;
                    hereditary_scale += term.to_f64().abs();
                    memories[j].push(z1);
                    // Use the original condensed force to retain weak losses
                    // when subtracting two large midpoint memories would round.
                    middle_extensions[j].push(arm_force / arm.coupling);
                }
                let defect =
                    (momentum + drag + spring + hereditary - force - applied_force).to_f64();
                let scale = hereditary_scale
                    + [momentum, drag, spring, force, applied_force]
                        .iter()
                        .map(|v| v.to_f64().abs())
                        .sum::<f64>();
                let kinematics = (Dd::from_f64(end_w)
                    - Dd::from_f64(previous_w[i])
                    - Dd::from_f64(h) * Dd::from_f64(mid_v))
                .to_f64();
                let kinematic_scale = end_w.abs() + previous_w[i].abs() + h * mid_v.abs();
                if !defect.is_finite()
                    || !scale.is_finite()
                    || !kinematics.is_finite()
                    || !kinematic_scale.is_finite()
                {
                    return Err(invalid("mechanical wall physical residual overflow").into());
                }
                backward_error = backward_error
                    .max(relative(defect.abs(), scale))
                    .max(relative(kinematics.abs(), kinematic_scale));
                displacement.push(end_w);
                velocity.push(end_v);
                middle_v.push(mid_v);
                pressure.push(p);
            }
            let mut projected = vec![0.0; wall.dofs.len()];
            wall.mass.spmv(&middle_v, &mut projected);
            for (&dof, &value) in wall.dofs.iter().zip(&projected) {
                wall_load[dof] += value;
            }
            let (work, rw) = paired_work(h, &pressure, &projected);
            wall_work += work;
            wall_work_roundoff += rw;
            let (work, roundoff) = paired_work(h, &traction_fields[index], &projected);
            traction_work = traction_work + Dd::from_f64(work);
            traction_roundoff += roundoff;
            let (loss, rl) = positive_form(&wall.mass, &middle_v, cancelled)?;
            wall_loss += h * pin.resistance * loss;
            wall_loss_roundoff += h * pin.resistance * rl;
            wall_change += pin.surface_density
                * energy_change(&wall.mass, previous_v, &velocity, cancelled)?
                + pin.stiffness_per_area
                    * energy_change(&wall.mass, previous_w, &displacement, cancelled)?;
            for (j, arm) in wall.relaxation.iter().enumerate() {
                poll(cancelled)?;
                let previous_extension = previous_w
                    .iter()
                    .zip(&initial.memories[index][j])
                    .map(|(&w, &z)| arm.coupling.mul_add(w, -z))
                    .collect::<Vec<_>>();
                let next_extension = displacement
                    .iter()
                    .zip(&memories[j])
                    .map(|(&w, &z)| arm.coupling.mul_add(w, -z))
                    .collect::<Vec<_>>();
                wall_change +=
                    energy_change(&wall.mass, &previous_extension, &next_extension, cancelled)?;
                let (loss, error) = positive_form(&wall.mass, &middle_extensions[j], cancelled)?;
                let loss = (loss / arm.branch.relaxation_time_s) * h;
                wall_loss += loss;
                relaxation_loss += loss;
                wall_loss_roundoff += (error / arm.branch.relaxation_time_s) * h;
            }
            next.displacements.push(displacement);
            next.velocities.push(velocity);
            next.memories.push(memories);
        }
        let admittance = relaxing::midpoint(
            &self.admittance,
            &self.admittance_weights,
            &initial.admittance_memories,
            &midpoint_p,
            h,
            cancelled,
        )?;
        next.admittance_memories = admittance.memories;
        for (combined, &memory) in wall_load.iter_mut().zip(&admittance.outward) {
            *combined = (Dd::from_f64(*combined) + Dd::from_f64(memory)).to_f64();
        }
        backward_error = backward_error.max(admittance.backward_error);
        let (_, fluid_error) = core.physical_residual(
            &initial.acoustic,
            &next.acoustic,
            &midpoint_p,
            &midpoint_psi,
            &self.conductance,
            &source,
            Some(&wall_load),
            cancelled,
        )?;
        backward_error = backward_error.max(fluid_error);
        if backward_error > core.limits.residual_tolerance {
            return Err(HeterogeneousError::Residual {
                observed: backward_error,
                limit: core.limits.residual_tolerance,
            });
        }
        let (
            initial_energy,
            initial_fluid_roundoff,
            initial_wall_roundoff,
            initial_admittance_roundoff,
        ) = self.energy_with_roundoff(initial, cancelled)?;
        let (final_energy, final_fluid_roundoff, final_wall_roundoff, final_admittance_roundoff) =
            self.energy_with_roundoff(&next, cancelled)?;
        let (boundary_loss, boundary_roundoff) =
            positive_form(&self.conductance, &midpoint_p, cancelled)?;
        let (input, input_roundoff) = paired_work(h, &midpoint_p, &source);
        let mut boundary_work = Dd::from_f64(0.0);
        for &(boundary, amplitude) in boundary_sources {
            poll(cancelled)?;
            let pressure = boundary.midpoint_pressure(initial.pressure_pa(), next.pressure_pa());
            boundary_work =
                boundary_work + Dd::from_f64(h) * Dd::from_f64(amplitude) * Dd::from_f64(pressure);
        }
        let boundary_source_work_j = boundary_work.to_f64();
        if !boundary_source_work_j.is_finite() {
            return Err(invalid("mechanical midpoint boundary-source work overflow").into());
        }
        let mut volume_work = Dd::from_f64(0.0);
        for &(volume, amplitude) in sources {
            poll(cancelled)?;
            let pressure = volume.midpoint_pressure(initial.pressure_pa(), next.pressure_pa());
            volume_work =
                volume_work + Dd::from_f64(h) * Dd::from_f64(amplitude) * Dd::from_f64(pressure);
        }
        let volume_source_work_j = volume_work.to_f64();
        if !volume_source_work_j.is_finite() {
            return Err(invalid("mechanical midpoint volume-source work overflow").into());
        }
        let work = MechanicalTimeWork {
            initial_energy_j: initial_energy.total_j(),
            final_energy_j: final_energy.total_j(),
            delta_energy_j: (Dd::from_f64(acoustic.work.delta_energy_j)
                + Dd::from_f64(wall_change)
                + Dd::from_f64(admittance.energy_change))
            .to_f64(),
            acoustic_delta_energy_j: acoustic.work.delta_energy_j,
            wall_delta_energy_j: wall_change,
            admittance_delta_energy_j: admittance.energy_change,
            input_work_j: (Dd::from_f64(input) + traction_work).to_f64(),
            acoustic_input_work_j: input,
            boundary_source_work_j,
            volume_source_work_j,
            wall_traction_work_j: traction_work.to_f64(),
            boundary_dissipation_j: h * boundary_loss,
            wall_dissipation_j: wall_loss,
            wall_relaxation_dissipation_j: relaxation_loss,
            wall_work_j: wall_work,
            admittance_dissipation_j: admittance.heat,
            admittance_work_j: admittance.work,
            backward_error,
        };
        let fluid_roundoff =
            initial_fluid_roundoff + final_fluid_roundoff + h * boundary_roundoff + input_roundoff;
        let wall_roundoff =
            initial_wall_roundoff + final_wall_roundoff + wall_loss_roundoff + traction_roundoff;
        let admittance_roundoff =
            initial_admittance_roundoff + final_admittance_roundoff + admittance.heat_roundoff;
        let tolerance = core.limits.energy_balance_tolerance;
        balance(
            work.acoustic_balance_residual_j(),
            work.acoustic_delta_energy_j.abs()
                + work.boundary_dissipation_j
                + wall_work.abs()
                + admittance.work.abs()
                + input.abs(),
            fluid_roundoff + wall_work_roundoff + admittance.work_roundoff,
            tolerance,
        )?;
        balance(
            work.wall_balance_residual_j(),
            wall_change.abs() + wall_loss + wall_work.abs() + work.wall_traction_work_j.abs(),
            wall_roundoff + wall_work_roundoff,
            tolerance,
        )?;
        balance(
            work.admittance_balance_residual_j(),
            admittance.energy_change.abs() + admittance.heat + admittance.work.abs(),
            admittance_roundoff + admittance.work_roundoff,
            tolerance,
        )?;
        balance(
            work.balance_residual_j(),
            work.delta_energy_j.abs()
                + work.boundary_dissipation_j
                + wall_loss
                + admittance.heat
                + input.abs()
                + work.wall_traction_work_j.abs(),
            fluid_roundoff + wall_roundoff + admittance_roundoff,
            tolerance,
        )?;
        poll(cancelled)?;
        Ok(MechanicalTimeStep { state: next, work })
    }
}

impl PreparedWall {
    fn history(&self, velocity: f64, displacement: f64, memory: &[Vec<f64>], dof: usize) -> f64 {
        let mut history = self
            .inertia_weight
            .mul_add(velocity, -self.spring_weight * displacement);
        for (arm, z) in self.relaxation.iter().zip(memory) {
            let extension = arm.coupling.mul_add(displacement, -z[dof]);
            history = (-arm.history_weight).mul_add(extension, history);
        }
        history
    }
}

pub(super) fn paired_work(h: f64, pressure: &[f64], load: &[f64]) -> (f64, f64) {
    let mut result = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (&p, &b) in pressure.iter().zip(load) {
        let term = Dd::from_f64(h) * Dd::from_f64(p) * Dd::from_f64(b);
        result = result + term;
        absolute += term.to_f64().abs();
    }
    (result.to_f64(), 64.0 * f64::EPSILON * absolute)
}

pub(super) fn balance(
    residual: f64,
    scale: f64,
    roundoff: f64,
    tolerance: f64,
) -> Result<(), HeterogeneousError> {
    let limit = tolerance * scale + roundoff;
    if !residual.is_finite() || !limit.is_finite() {
        return Err(invalid("mechanical acoustic work is unrepresentable").into());
    }
    if residual.abs() > limit {
        return Err(HeterogeneousError::EnergyBalance {
            observed_j: residual,
            limit_j: limit,
        });
    }
    Ok(())
}

#[cfg(test)]
mod boundary_source_tests;
#[cfg(test)]
mod composed_tests;
#[cfg(test)]
mod pressure_release_tests;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod traction_tests;
#[cfg(test)]
mod volume_source_tests;
