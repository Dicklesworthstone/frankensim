//! Full hereditary shell material loss in the existing modal-memory owner.
//!
//! The equilibrium pencil is unchanged. Each projected excess stiffness is
//! factored in its original units and enters the original first-order material
//! coordinates; no scalar damping ratio, sampled viscosity or extra mass is
//! substituted for its frequency-dependent storage and dissipation.

use super::{AcousticRealizeError, ModalPlateViscosity, invalid, owner, poll, validate_mass_basis};
use crate::plate_viscosity::relaxation::{InitialPlateMemory, PlateRelaxationLimits};
use crate::render::schedule::force::coupled::CoupledModalSystem;
use fs_exec::CancelGate;
use fs_phs::RelaxationBranch;
use fs_plate::shell::relaxation::{
    ProjectedShellRelaxationBranch, ShellRelaxationBranch, ShellRelaxationBudget,
    project_relaxation, projection_terms,
};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::{ShellMesh, ShellModel};
use std::collections::BTreeSet;

// The primary shear-deformable shell supplies the same projected spectra and
// existing memory runtime as the CST/DKT model-form rung.
mod mitc;

/// Complete positive material strain channels on the unchanged shell basis.
/// Each row c [1/s] has storage `(c.q-z)^2/2` [J] and heat rate
/// `(c.q-z)^2/tau` [W]. An exact null strain carries no material state.
#[derive(Clone, Debug)]
pub struct ModalShellRelaxation {
    projected: Vec<ProjectedShellRelaxationBranch>,
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
    modes: usize,
    limits: PlateRelaxationLimits,
}

impl ModalShellRelaxation {
    /// Project per-facet passive excess A/B/D spectra on the original
    /// mass-orthonormal CST/DKT modes. Lumped and consistent mass use the SAME
    /// sparse-mass check as instantaneous shell viscosity (entry tolerance
    /// 1e-8); neither modes nor eigenvalues are refitted.
    ///
    /// Supply stiffness increments above the fully relaxed elastic pencil.
    /// Each facet may have different tensor directions and relaxation times.
    /// Only exactly equal times aggregate, before the common compensated PSD
    /// factorization; every resolved positive direction and cross term remains.
    /// Up to 64 original modes and 8192 complete memory states are admitted by
    /// this cold adapter; the selected coupled runtime keeps its own tighter
    /// work, state, initial-energy and instantaneous-frequency limits.
    ///
    /// `max_projection_terms` bounds geometry/strain projection, every original
    /// mass Gram entry, all full factors and later branch embedding. Source
    /// spectra, band support, equilibrium interpretation and mesh/basis
    /// identity remain caller obligations. Bare CST/DKT map restrictions are
    /// those of [`fs_plate::shell::relaxation::project_relaxation`].
    pub fn from_cst_dkt_shell(
        mesh: &ShellMesh,
        spectra: &[Vec<ShellRelaxationBranch>],
        model: &ShellModel,
        modes: &[Vec<f64>],
        limits: PlateRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let n = modes.len();
        let full_dofs = mesh
            .nodes
            .len()
            .checked_mul(6)
            .ok_or_else(|| invalid("shell relaxation coordinate extent overflows"))?;
        if !(1..=64).contains(&limits.max_modes)
            || n == 0
            || n > limits.max_modes
            || limits.max_memory_states > 8192
            || spectra.len() != mesh.tris.len()
            || model.free == 0
            || model.free > full_dofs
            || model.m.nrows() != model.free
            || model.m.ncols() != model.free
            || modes.iter().any(|mode| mode.len() != model.free)
        {
            return Err(invalid(
                "shell relaxation requires its complete bounded original basis and spectrum",
            ));
        }
        let mut branches = 0usize;
        let mut times = BTreeSet::new();
        for section in spectra {
            poll(gate)?;
            branches = branches
                .checked_add(section.len())
                .ok_or_else(|| invalid("shell relaxation source extent overflows"))?;
            if branches > limits.max_projection_terms {
                return Err(invalid("shell relaxation source exceeds its declared work"));
            }
            for branch in section {
                poll(gate)?;
                if !branch.is_zero() {
                    times.insert(branch.relaxation_time_s().to_bits());
                }
            }
        }
        let spatial_work = projection_terms(full_dofs, mesh.tris.len(), n, branches, times.len())
            .ok_or_else(|| invalid("shell relaxation projection work overflows"))?;
        let factor_work = (64 * n * n * n)
            .checked_mul(times.len())
            .ok_or_else(|| invalid("shell relaxation factor work overflows"))?;
        let basis_work = model
            .m
            .nnz()
            .checked_mul(n)
            .and_then(|work| {
                model
                    .free
                    .checked_mul(n * n)
                    .and_then(|gram| work.checked_add(gram))
            })
            .and_then(|work| work.checked_add(factor_work))
            .ok_or_else(|| invalid("shell relaxation basis work overflows"))?;
        if spatial_work
            .checked_add(basis_work)
            .is_none_or(|work| work > limits.max_projection_terms)
        {
            return Err(invalid(
                "complete shell relaxation exceeds its declared projection work",
            ));
        }
        validate_mass_basis(model, modes, gate)?;
        let projected = project_relaxation(
            mesh,
            spectra,
            model,
            modes,
            ShellRelaxationBudget {
                projection: ShellViscosityBudget {
                    max_elements: mesh.tris.len(),
                    max_full_dofs: full_dofs,
                    max_dofs: model.free,
                    max_modes: limits.max_modes,
                    max_work: limits.max_projection_terms - basis_work,
                },
                max_element_branches: branches,
                max_relaxation_times: times.len(),
            },
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            fs_plate::PlateError::Cancelled => AcousticRealizeError::Cancelled,
            other => owner(other),
        })?;
        Self::from_projected(projected, n, limits, gate)
    }

    // Shared by the original CST/DKT and MITC physical projection owners.
    // Both callers admit their original basis and complete setup work before
    // material stiffness is converted into existing Maxwell coordinates.
    fn from_projected(
        projected: Vec<ProjectedShellRelaxationBranch>,
        n: usize,
        limits: PlateRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        let mut rows = Vec::new();
        let mut times = Vec::new();
        for branch in &projected {
            poll(gate)?;
            // Use the common compensated factor on the ORIGINAL tensor.
            // Normalizing/rounding the matrix first can erase a correlated
            // weak direction before its physical strain is ever resolved.
            let factor =
                ModalPlateViscosity::from_psd_modal_matrix(branch.matrix().to_vec(), n, 1.0, gate)?;
            if rows
                .len()
                .checked_add(factor.channel_count())
                .is_none_or(|states| states > limits.max_memory_states)
            {
                return Err(invalid(
                    "complete shell relaxation exceeds its declared memory states",
                ));
            }
            times.extend(std::iter::repeat_n(
                branch.relaxation_time_s(),
                factor.channel_count(),
            ));
            rows.extend(factor.channel_shapes().iter().cloned());
        }
        poll(gate)?;
        Ok(Self {
            projected,
            rows,
            times,
            modes: n,
            limits,
        })
    }

    pub const fn mode_count(&self) -> usize {
        self.modes
    }

    pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] {
        &self.rows
    }

    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.times
    }

    pub fn projected_branches(&self) -> &[ProjectedShellRelaxationBranch] {
        &self.projected
    }

    pub const fn limits(&self) -> PlateRelaxationLimits {
        self.limits
    }

    /// Embed physical memory rows in a complete interleaved mechanical layout
    /// `[q0,v0,q1,v1,...]`. This public bridge also serves existing fs-phs
    /// consumers; all velocity and other-component coefficients remain zero.
    /// `first_mode` is the shell's original contiguous coordinate offset.
    pub fn branches(
        &self,
        first_mode: usize,
        total_modes: usize,
        gate: &CancelGate,
    ) -> Result<Vec<RelaxationBranch>, AcousticRealizeError> {
        poll(gate)?;
        let dimension = total_modes
            .checked_mul(2)
            .ok_or_else(|| invalid("shell relaxation embedding extent overflows"))?;
        if first_mode
            .checked_add(self.modes)
            .is_none_or(|end| end > total_modes)
            || self
                .rows
                .len()
                .checked_mul(dimension)
                .is_none_or(|work| work > self.limits.max_projection_terms)
        {
            return Err(invalid(
                "shell relaxation embedding exceeds its original layout or work",
            ));
        }
        let mut branches = Vec::with_capacity(self.rows.len());
        for (row, &time) in self.rows.iter().zip(&self.times) {
            poll(gate)?;
            let mut projection = vec![0.0; dimension];
            for (mode, &value) in row.iter().enumerate() {
                projection[2 * (first_mode + mode)] = value;
            }
            branches.push(RelaxationBranch {
                projection,
                stiffness: 1.0,
                relaxation_time_s: time,
            });
        }
        poll(gate)?;
        Ok(branches)
    }

    /// Compose this image into the existing complete held-force propagator
    /// before playback. Earlier material histories remain unchanged. Relaxed
    /// initialization sets each new z=c.q; unrelaxed initialization sets z=0.
    /// Every interval's physical heat is reported by the original
    /// `CoupledModalFrame::material_dissipation_j`, independently of endpoint
    /// stored-energy subtraction. No new integrator or heat debit is added.
    ///
    /// The target component must use this original equilibrium basis and zero
    /// intrinsic modal damping. Independently selected viscosity, glue,
    /// support and radiation channels may coexist without being inserted here.
    pub fn attach_to(
        &self,
        system: CoupledModalSystem,
        component: usize,
        initial: InitialPlateMemory,
        gate: &CancelGate,
    ) -> Result<CoupledModalSystem, AcousticRealizeError> {
        poll(gate)?;
        let target = system
            .components()
            .get(component)
            .ok_or_else(|| invalid("shell relaxation component is absent"))?;
        if system.samples_rendered() != 0
            || target.modes().len() != self.modes
            || (0..self.modes).any(|i| target.damping_rate_per_s(i) != 0.0)
        {
            return Err(invalid(
                "shell relaxation needs the original cold component with zero intrinsic damping",
            ));
        }
        let start = system.components()[..component]
            .iter()
            .map(|model| model.modes().len())
            .sum();
        let branches = self.branches(start, system.mode_count(), gate)?;
        let history = match initial {
            InitialPlateMemory::Relaxed => None,
            InitialPlateMemory::Unrelaxed => Some(vec![0.0; branches.len()]),
        };
        system
            .with_additional_material_relaxation(branches, history, gate)
            .map_err(owner)
    }
}

#[cfg(test)]
mod tests;
