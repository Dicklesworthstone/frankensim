//! Held-jump coefficient changes without resetting adhesive relaxation memory.
//!
//! The state stores physical Maxwell SPRING strains q, not coefficient-scaled
//! forces. At fixed total strain and q, psi=h/2*(Einf*epsilon^2+sum Ej*qj^2).
//! Changing spring moduli supplies signed parameter work; changing relaxation
//! times alone supplies exactly zero. No physical time, displacement work or
//! dashpot heat occurs in this transfer. Past work/heat and elapsed time remain
//! unchanged. Subsequent exact ramps use the new kinetics and the SAME memory.
//!
//! This is a coefficients-only continuation, not a thermal-expansion, entropy,
//! heat-capacity or heat-partition law. The source owner separately proves that
//! each ordered branch is the same physical internal variable at the new query.

use super::{CohesiveError, IntactViscoelasticBond, ViscoelasticBondPoint, finite, positive};

/// Accepted instantaneous coefficient continuation per original bonded area.
#[derive(Debug, Clone)]
pub struct ViscoelasticBondTemperatureStep<'a> {
    pub point: ViscoelasticBondPoint<'a>,
    /// Signed parameter work [J/m2]; never reported as viscous heat.
    pub parameter_work_j_m2: f64,
    /// Independently evaluated change of spring storage [J/m2].
    pub stored_change_j_m2: f64,
    /// Parameter work minus storage change [J/m2].
    pub defect_j_m2: f64,
}

// Work integrates psi_E dE directly. Storage uses the change of physical spring
// force at fixed strain, with FMA product residuals retaining small differences.
// Unchanged springs contribute exact zero even beside a large stored inventory.
fn spring_change(
    old: f64,
    new: f64,
    strain: f64,
    thickness: f64,
) -> Result<(f64, f64), CohesiveError> {
    if old == new || strain == 0.0 {
        return Ok((0.0, 0.0));
    }
    let physical_gap = finite(thickness * strain)?;
    let work = finite((0.5 * finite((new - old) * strain)?) * physical_gap)?;
    let before = finite(old * strain)?;
    let after = finite(new * strain)?;
    let correction = new.mul_add(strain, -after) - old.mul_add(strain, -before);
    let force_change = finite((after - before) + correction)?;
    let storage = finite((0.5 * force_change) * physical_gap)?;
    Ok((work, storage))
}

impl ViscoelasticBondPoint<'_> {
    /// Continue this original point to explicitly supplied temperature
    /// coefficients at fixed physical jump and fixed physical spring strains.
    /// Thickness, local-use limits and complete branch counts must agree.
    ///
    /// Branch INDEX is its declared identity: this operation neither sorts
    /// relaxation times nor invents a remapping. A source-bound caller MUST first
    /// call `ResolvedViscoelasticBond::require_temperature_successor`, which
    /// additionally checks the original indexed coefficient claims and query.
    /// The numeric law alone cannot establish physical/source branch identity.
    ///
    /// No mechanical work or heat is added to the original Prony ledgers.
    /// The caller owns this signed parameter-work account and must prepare the
    /// successor law for subsequent motion. Every refusal leaves this point
    /// retryable, including an unresolved independent storage/work balance.
    pub fn transfer_temperature<'b>(
        &self,
        law: &'b IntactViscoelasticBond,
        energy_tolerance_j_m2: f64,
    ) -> Result<ViscoelasticBondTemperatureStep<'b>, CohesiveError> {
        positive(energy_tolerance_j_m2)?;
        if self.law.thickness != law.thickness
            || self.law.limits != law.limits
            || (0..2)
                .any(|mode| self.law.spectra[mode].terms.len() != law.spectra[mode].terms.len())
        {
            return Err(CohesiveError::Invalid(
                "bond temperature transfer requires original geometry, bounds and ordered branches",
            ));
        }
        let mut point = ViscoelasticBondPoint {
            law,
            states: self.states.clone(),
            elapsed: self.elapsed,
        };
        let mut parameter_work = 0.0;
        let mut storage_change = 0.0;
        for i in 0..3 {
            let mode = usize::from(i > 0);
            let before = &self.law.spectra[mode];
            let after = &law.spectra[mode];
            let state = &self.states[i];
            let mut same_storage = before.e_inf == after.e_inf;
            let (work, storage) =
                spring_change(before.e_inf, after.e_inf, state.strain, law.thickness)?;
            parameter_work = finite(parameter_work + work)?;
            storage_change = finite(storage_change + storage)?;
            for (j, &q) in state.q.iter().enumerate() {
                let a = before.terms[j].0;
                let b = after.terms[j].0;
                same_storage &= a == b;
                let (work, storage) = spring_change(a, b, q, law.thickness)?;
                parameter_work = finite(parameter_work + work)?;
                storage_change = finite(storage_change + storage)?;
            }
            if !same_storage {
                let mut stress = finite(after.e_inf * state.strain)?;
                for (&(e, _), &q) in after.terms.iter().zip(&state.q) {
                    stress = finite(stress + e * q)?;
                }
                point.states[i].stress = stress;
            }
        }
        finite(point.stored_energy_j_m2())?;
        let defect = finite(parameter_work - storage_change)?;
        if defect.abs() > energy_tolerance_j_m2 {
            return Err(CohesiveError::EnergyDefect(defect));
        }
        Ok(ViscoelasticBondTemperatureStep {
            point,
            parameter_work_j_m2: parameter_work,
            stored_change_j_m2: storage_change,
            defect_j_m2: defect,
        })
    }
}

#[cfg(test)]
mod tests;
