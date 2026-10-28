//! Explicit changing-coefficient continuation of the original return map.
//!
//! The parameter port supplies the change in the spring potential at the old
//! elastic history. Motion then advances with the new coefficients. This is a
//! first-order, parameter-then-motion split; it is not an inferred thermal law.
//! For temperature-dependent spring free energy, caloric coupling must also
//! account for entropy and U = psi - T * d(psi)/dT. Neither parameter work nor
//! backward-Euler numerical loss may be relabeled as physical friction heat.

use super::*;

impl PartialSlipLaw {
    /// Continue an accepted contact with explicitly updated scalar coefficients.
    ///
    /// The same model/source and ordered contact lineage must be retained. The
    /// two regions keep their original stiffness/capacity fractions, elastic strains,
    /// and kinetic branches; a new regional decomposition requires another law.
    /// Coulomb coefficients and translational/torsional stiffness may change;
    /// the torsional capacity factor remains fixed.
    /// The immutable source identity names the caller's varying-coefficient
    /// model, not a new source selected opportunistically at each step.
    ///
    /// Before the ORIGINAL return map runs, the parameter port supplies
    /// `C = psi(new coefficients, old state) - psi(old coefficients, old state)`.
    /// The returned storage change is measured against the actual old law, so
    /// `delta(psi) + physical heat + numerical loss = wrench work + C`.
    /// Softening can export parameter work. Changed yield capacity may cause
    /// real relaxation/plastic heat even with zero imposed displacement.
    ///
    /// This method does not mutate the supplied checkpoint. Every trial must
    /// reuse the same accepted checkpoint; the owner commits the returned new
    /// law/history AND parameter-work exchange atomically. Ordinary
    /// [`Self::restore_checkpoint`] remains strict about every coefficient.
    /// No source validity, temperature mapping, heat allocation, or material
    /// thermodynamic potential is invented by this source-neutral leaf.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coefficients(
        &self,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        frame: TangentFrame,
        kinematics: PartialSlipKinematics,
        ownership: &GeneralizedWorkOwnership,
        checkpoint: &PartialSlipCheckpoint,
    ) -> Result<PartialSlipStep, PartialSlipError> {
        self.validate()?;
        let old_law = checkpoint.law();
        // Validate the accepted law and full lineage, not a synthetic rebound
        // checkpoint that could silently forget the original stored energy.
        let state = old_law.restore_checkpoint(patch, interface, checkpoint)?;
        for (same, field) in [
            (self.model_id == old_law.model_id, "law_model_id"),
            (self.source_id == old_law.source_id, "law_source_id"),
            (
                self.parameters.partial_slip_hardening_fraction
                    == old_law.parameters.partial_slip_hardening_fraction,
                "region_stiffness_fraction",
            ),
            (
                self.parameters.partial_slip_onset_fraction
                    == old_law.parameters.partial_slip_onset_fraction,
                "region_capacity_fraction",
            ),
            (
                self.parameters.torsional_capacity_factor
                    == old_law.parameters.torsional_capacity_factor,
                "torsional_capacity_factor",
            ),
        ] {
            if !same {
                return Err(PartialSlipError::CheckpointIdentityMismatch { field });
            }
        }
        let mut step = self.advance(patch, interface, frame, kinematics, ownership, &state)?;
        if self.parameters == old_law.parameters {
            // Keep the unchanged-law path bit-identical to the original law.
            return Ok(step);
        }
        let coefficient_work_j = coefficient_work(old_law.parameters, self.parameters, &state)?;
        step.storage_change_j = checked_add(
            step.storage_change_j,
            coefficient_work_j,
            "coefficient_storage_change",
        )?;
        step.generalized_work.coefficient_work_j = coefficient_work_j;
        step.generalized_work.closure_residual_j = checked_add(
            checked_add(
                checked_add(
                    checked_add(
                        step.storage_change_j,
                        step.dissipation.heat_j,
                        "coefficient_physical_change",
                    )?,
                    step.generalized_work.numerical_dissipation_j,
                    "coefficient_accounted_change",
                )?,
                -step.generalized_work.work_into_interface_j,
                "coefficient_wrench_balance",
            )?,
            -coefficient_work_j,
            "coefficient_energy_closure",
        )?;
        Ok(step)
    }
}

fn coefficient_work(
    old: PartialSlipParameters,
    new: PartialSlipParameters,
    state: &PartialSlipState,
) -> Result<f64, PartialSlipError> {
    let hard = old.partial_slip_hardening_fraction;
    let mut work = 0.0;
    for (region, weight) in [1.0 - hard, hard].into_iter().enumerate() {
        let old_stiffness = [
            old.tangential_stiffness_n_per_m * weight,
            old.tangential_stiffness_n_per_m * weight,
            old.torsional_stiffness_nm_per_rad * weight,
        ];
        let new_stiffness = [
            new.tangential_stiffness_n_per_m * weight,
            new.tangential_stiffness_n_per_m * weight,
            new.torsional_stiffness_nm_per_rad * weight,
        ];
        let strain = [
            state.region_elastic_displacement_m[region][0],
            state.region_elastic_displacement_m[region][1],
            state.region_elastic_twist_rad[region],
        ];
        for i in 0..3 {
            let delta = checked_add(
                new_stiffness[i],
                -old_stiffness[i],
                "coefficient_stiffness_change",
            )?;
            // Evaluate the stiffness difference before multiplying by strain;
            // subtracting two large stored energies loses small parameter work.
            let contribution = checked_mul(
                checked_mul(0.5 * delta, strain[i], "coefficient_work")?,
                strain[i],
                "coefficient_work",
            )?;
            work = checked_add(work, contribution, "coefficient_work")?;
        }
    }
    Ok(work)
}
