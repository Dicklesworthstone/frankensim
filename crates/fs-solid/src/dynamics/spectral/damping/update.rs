//! Instantaneous stiffness and viscosity changes on the retained physical mass.
//!
//! Changing the constitutive operators does not erase an existing vibration,
//! create a damping impulse, or reset the work and heat history.

use super::super::{CompleteSpectralSolidError, SpectralMass};
use super::{CompleteSpectralSolid, DampedSpectralSolid, Result, add, bound, invalid, poll};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_sparse::Csr;

/// Accepted zero-duration update on the same original coordinates and mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralParameterUpdate {
    pub time_s: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Unchanged physical kinetic inventory; no momentum impulse is applied.
    pub kinetic_energy_j: f64,
    /// Signed 1/2 u^T (K_after - K_before) u, directly contracted from the
    /// operators. Positive work must be supplied; negative work is returned.
    pub parameter_work_j: f64,
    pub cumulative_parameter_work_j: f64,
    /// Parameter work minus the independently evaluated storage change [J].
    pub balance_error_j: f64,
    /// Force/impulse work + parameter work - heat - (H - H_initial) [J].
    pub cumulative_balance_error_j: f64,
}

impl DampedSpectralSolid {
    /// Replace stiffness, classical viscosity and their complete eigenbasis.
    ///
    /// The original mass, free-DOF map, reference coordinate origin and sample
    /// clock are retained by this owner: none is an input to this operation.
    /// The caller supplies K/C on THAT physical space, with no geometry,
    /// support, equilibrium or inertia remapping. The complete replacement
    /// basis is admitted against the retained original M and supplied K.
    /// Physical displacement and velocity are kept bit for bit, and projected
    /// anew through M into every replacement mode. Mode order, sign, repeated
    /// eigenspaces and elastic/rigid transitions do not identify the state.
    ///
    /// Signed stiffness parameter work is independent of applied force/impulse
    /// work and viscous heat. Changing C alone performs no instantaneous work
    /// or heat. All prior accounts, their compensation, the original initial
    /// energy, accepted step count and clock survive. Subsequent intervals and
    /// impulses close the same cumulative physical balance, including after
    /// switching C to zero.
    ///
    /// `maximum_parameter_work_j` bounds the ABSOLUTE magnitude of this
    /// physical event's signed work. It is not a numerical residual tolerance;
    /// the existing operator, state and work-error controls still apply.
    /// The caller owns the source/sink providing this prescribed parameter
    /// work. This is a piecewise constant physical-parameter event, not a
    /// continuous variable-coefficient integrator, thermodynamic heat law or
    /// transfer of constitutive memory.
    ///
    /// The existing preparation budget covers the live owner plus the entire
    /// unpublished candidate: `128*n^2 + 4096*n` bytes and
    /// `6*n^3 + 80*n^2 + 8192*n + 512*sum(s_i)` scalar-product terms, with
    /// the same oscillator scaling count as ordinary damping preparation.
    /// Cancellation, invalid operators/basis, insufficient work allowance or
    /// any energy gate leaves all accepted state and caches unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn replace_parameters(
        &mut self,
        stiffness_n_per_m: &Csr,
        viscosity_kg_per_s: &Csr,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        maximum_parameter_work_j: f64,
        cx: &Cx<'_>,
    ) -> Result<SpectralParameterUpdate> {
        let (candidate, report) = self.parameter_candidate(
            stiffness_n_per_m,
            viscosity_kg_per_s,
            eigenvalues_rad2_s2,
            basis_columns,
            maximum_parameter_work_j,
            cx,
        )?;
        *self = candidate;
        Ok(report)
    }

    // Source-specific owners must rebuild their caches before publishing the
    // same physical candidate. No old state or account is changed here.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn parameter_candidate(
        &self,
        stiffness_n_per_m: &Csr,
        viscosity_kg_per_s: &Csr,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        maximum_parameter_work_j: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self, SpectralParameterUpdate)> {
        poll(cx)?;
        if !maximum_parameter_work_j.is_finite() || maximum_parameter_work_j < 0.0 {
            return Err(invalid("finite nonnegative parameter-work budget required"));
        }
        let n = self.free_dofs().len() as u128;
        let controls = self.solid.controls;
        bound(
            128 * n * n + 4096 * n,
            controls.maximum_preparation_bytes as u128,
        )?;
        let terms = 6 * n * n * n + 80 * n * n + 8192 * n;
        bound(terms, controls.maximum_preparation_multiply_terms as u128)?;
        let replacement = CompleteSpectralSolid::prepare_with_operators(
            stiffness_n_per_m,
            SpectralMass::Retained(&self.solid.mass),
            &self.solid.free_dofs,
            eigenvalues_rad2_s2,
            basis_columns,
            self.solid.dt_s,
            &self.solid.displacement,
            &self.solid.velocity,
            controls,
            cx,
        )?;
        let mut candidate = Self::prepare_viscosity(replacement, viscosity_kg_per_s, terms, cx)?;
        let coefficient_change = parameter_work(
            &self.solid.stiffness,
            &candidate.solid.stiffness,
            &self.solid.displacement,
            cx,
        )?;
        let work = coefficient_change.to_f64();
        if work.abs() > maximum_parameter_work_j {
            return Err(CompleteSpectralSolidError::Energy {
                defect_j: work,
                allowance_j: maximum_parameter_work_j,
            });
        }
        let (mut accumulated, mut correction) = (
            self.solid.cumulative_parameter_work_j,
            self.solid.parameter_work_correction_j,
        );
        if work != 0.0 {
            add(&mut accumulated, &mut correction, work);
        }
        if !accumulated.is_finite() || !correction.is_finite() {
            return Err(invalid("unrepresentable accumulated parameter work"));
        }

        // Keep the original reference coordinates. Only the coefficient
        // storage at those fixed coordinates changes; subsequent motion uses
        // the new original operator and a polarized reference-to-state change.
        // This account is independent of parameter work at the current u.
        candidate.solid.reference_parameter_change_j = self.solid.reference_parameter_change_j
            + parameter_work(
                &self.solid.stiffness,
                &candidate.solid.stiffness,
                &self.solid.reference_displacement,
                cx,
            )?;
        candidate.solid.reference_displacement =
            super::super::copy(&self.solid.reference_displacement)?;
        candidate.solid.reference_velocity = super::super::copy(&self.solid.reference_velocity)?;
        let storage_before =
            self.solid
                .storage_change(&self.solid.displacement, &self.solid.velocity, cx)?;
        let storage_after = candidate.solid.storage_change(
            &candidate.solid.displacement,
            &candidate.solid.velocity,
            cx,
        )?;

        // Preparation copies physical u/v exactly and checks reconstruction
        // in the new complete basis. Preserve the accepted kinetic bits too.
        candidate.solid.kinetic_j = self.solid.kinetic_j;
        candidate.solid.steps = self.solid.steps;
        candidate.solid.cumulative_work_j = self.solid.cumulative_work_j;
        candidate.solid.work_correction_j = self.solid.work_correction_j;
        candidate.solid.cumulative_parameter_work_j = accumulated;
        candidate.solid.parameter_work_correction_j = correction;
        candidate.cumulative_dissipated_j = self.cumulative_dissipated_j;
        candidate.dissipation_correction_j = self.dissipation_correction_j;

        let report = SpectralParameterUpdate {
            time_s: self.time_s(),
            stored_before_j: self.solid.stored_j,
            stored_after_j: candidate.solid.stored_j,
            kinetic_energy_j: self.solid.kinetic_j,
            parameter_work_j: work,
            cumulative_parameter_work_j: accumulated,
            balance_error_j: (coefficient_change - (storage_after - storage_before)).to_f64(),
            cumulative_balance_error_j: candidate.solid.cumulative_balance_error(
                candidate.solid.cumulative_work_j,
                candidate.solid.work_correction_j,
                candidate.cumulative_dissipated_j,
                candidate.dissipation_correction_j,
                storage_after,
            ),
        };
        for defect_j in [report.balance_error_j, report.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > controls.maximum_work_error_j {
                return Err(CompleteSpectralSolidError::Energy {
                    defect_j,
                    allowance_j: controls.maximum_work_error_j,
                });
            }
        }
        poll(cx)?;
        Ok((candidate, report))
    }
}

fn parameter_work(before: &[f64], after: &[f64], u: &[f64], cx: &Cx<'_>) -> Result<Dd> {
    let n = u.len();
    let mut work = Dd::ZERO;
    for i in 0..n {
        poll(cx)?;
        let mut row = Dd::ZERO;
        for j in 0..n {
            // Difference the original operator entries before contracting.
            // An unchanged large stiffness inventory cannot erase a smaller
            // work exchange on a different component.
            row = row
                + (Dd::from_f64(after[i * n + j]) - Dd::from_f64(before[i * n + j]))
                    * Dd::from_f64(u[j]);
        }
        work = work + Dd::from_f64(0.5 * u[i]) * row;
    }
    if !work.hi.is_finite() || !work.lo.is_finite() || !work.to_f64().is_finite() {
        return Err(invalid("unrepresentable stiffness parameter work"));
    }
    Ok(work)
}
