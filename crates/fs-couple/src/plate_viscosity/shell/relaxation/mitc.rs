//! Original supported MITC coordinates into the common hereditary runtime.

use super::*;
use fs_plate::shell::mitc::relaxation::{
    ShearShellRelaxationBranch, ShearShellRelaxationBudget,
    project_relaxation as project_mitc_relaxation, projection_terms as mitc_projection_terms,
};
use fs_plate::shell::mitc::{MitcBudget, MitcShellModel};

impl ModalShellRelaxation {
    /// Project independent facet A/B/D and transverse-shear spectra on the
    /// original mass-orthonormal MITC modes. The model retains the exact mesh,
    /// nodal directors, equilibrium bubble condensation, reference inertia and
    /// support frames that assembled its pencil. Its original expansion is
    /// used before physical material work; neither a bare map nor a diagonal
    /// loss approximation replaces it.
    ///
    /// Supply passive excess stiffness above the fully relaxed section. Each
    /// branch keeps its own positive material time; only exactly equal times
    /// aggregate. This is material memory on the original fixed condensed
    /// trial space, not frequency-dependent elimination of independent
    /// director bubbles. Numerical drilling has no material state.
    ///
    /// Up to 1024 original modes and 8192 material states are admitted here,
    /// subject to complete declared work and the selected runtime's tighter
    /// limits. Both original mass choices, including their translation/rotary
    /// coupling, use the same sparse mass Gram check as CST/DKT. Source bands
    /// and reference temperature remain explicit caller data.
    pub fn from_mitc_shell(
        model: &MitcShellModel,
        spectra: &[Vec<ShearShellRelaxationBranch>],
        modes: &[Vec<f64>],
        limits: PlateRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let pencil = model.pencil();
        let n = modes.len();
        let full_dofs = pencil.dof_map.len();
        if !(1..=1024).contains(&limits.max_modes)
            || n == 0
            || n > limits.max_modes
            || limits.max_memory_states > 8192
            || full_dofs % 6 != 0
            || pencil.free == 0
            || pencil.free > full_dofs
            || pencil.m.nrows() != pencil.free
            || pencil.m.ncols() != pencil.free
            || modes.iter().any(|mode| mode.len() != pencil.free)
        {
            return Err(invalid(
                "MITC relaxation requires its complete bounded original basis and spectrum",
            ));
        }
        let mut branches = 0usize;
        let mut times = BTreeSet::new();
        for section in spectra {
            poll(gate)?;
            branches = branches
                .checked_add(section.len())
                .ok_or_else(|| invalid("MITC relaxation source extent overflows"))?;
            if branches > limits.max_projection_terms {
                return Err(invalid("MITC relaxation source exceeds its declared work"));
            }
            for branch in section {
                poll(gate)?;
                if !branch.is_zero() {
                    times.insert(branch.relaxation_time_s().to_bits());
                }
            }
        }
        let nodes = full_dofs / 6;
        let spatial_work = mitc_projection_terms(nodes, spectra.len(), n, branches, times.len())
            .ok_or_else(|| invalid("MITC relaxation projection work overflows"))?;
        let square = n
            .checked_mul(n)
            .ok_or_else(|| invalid("MITC relaxation modal extent overflows"))?;
        let factor_work = square
            .checked_mul(n)
            .and_then(|work| work.checked_mul(64))
            .and_then(|work| work.checked_mul(times.len()))
            .ok_or_else(|| invalid("MITC relaxation factor work overflows"))?;
        let basis_work = pencil
            .m
            .nnz()
            .checked_mul(n)
            .and_then(|work| {
                pencil
                    .free
                    .checked_mul(square)
                    .and_then(|gram| work.checked_add(gram))
            })
            .and_then(|work| work.checked_add(factor_work))
            .ok_or_else(|| invalid("MITC relaxation basis work overflows"))?;
        if spatial_work
            .checked_add(basis_work)
            .is_none_or(|work| work > limits.max_projection_terms)
        {
            return Err(invalid(
                "complete MITC relaxation exceeds its declared projection work",
            ));
        }
        validate_mass_basis(pencil, modes, gate)?;
        let projected = project_mitc_relaxation(
            model,
            spectra,
            modes,
            ShearShellRelaxationBudget {
                mesh: MitcBudget {
                    max_nodes: nodes,
                    max_triangles: spectra.len(),
                },
                max_modes: limits.max_modes,
                max_element_branches: branches,
                max_relaxation_times: times.len(),
                max_projection_terms: limits.max_projection_terms - basis_work,
            },
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            fs_plate::PlateError::Cancelled => AcousticRealizeError::Cancelled,
            other => owner(other),
        })?;
        Self::from_projected(projected, n, limits, gate)
    }
}

#[cfg(test)]
mod tests;
