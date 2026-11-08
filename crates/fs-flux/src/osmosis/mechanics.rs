//! Material-owned nonlinear resistance to finite osmotic swelling.
//!
//! A lumped compartment can retain the same convex bulk energy used by solid
//! mechanics. If its initial solvent volume is Vw, fixed non-solvent volume is
//! Vs, and the original material point has volume ratio J0, then its material
//! reference volume is V0=(Vw+Vs)/J0. Subsequent solvent uptake gives
//! J=J0+(V-Vw)/V0, pressure p=p_background+U'(J), and storage V0 U(J).
//! Pressure, tangent, admissible branch and finite energy changes all come
//! from the original material point; no fitted linear spring replaces it.
//!
//! This is an effective homogeneous volumetric restraint. It does not infer
//! shell geometry, shear deformation, polymer mixing, poroelastic flow, damage,
//! or the mechanics of a particular biological membrane. Solute osmotic
//! pressure remains with the original aqueous owner and is applied once.

use super::{OsmoticCompartment, OsmoticNetworkError, OsmoticPressureVolume, nonnegative, positive};
use fs_material::hyper::discrete_gradient::bulk::DiscreteBulkPoint;

#[derive(Clone, Copy, Debug)]
pub(super) struct BulkPressure {
    point: DiscreteBulkPoint,
    reference_material_volume_m3: f64,
}

fn finite(value: f64) -> Result<f64, OsmoticNetworkError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OsmoticNetworkError::Invalid(
            "nonlinear swelling mechanics is not representable",
        ))
    }
}

impl BulkPressure {
    fn at(self, solvent_increment_m3: f64) -> Result<DiscreteBulkPoint, OsmoticNetworkError> {
        let j = self.point.accepted_volume_ratio()
            + solvent_increment_m3 / self.reference_material_volume_m3;
        self.point.rebind(j).map_err(OsmoticNetworkError::Mechanical)
    }

    pub(super) fn pressure_pa(self, solvent_increment_m3: f64) -> Result<f64, OsmoticNetworkError> {
        finite(self.at(solvent_increment_m3)?.endpoint_pressure_pa())
    }

    pub(super) fn tangent_pa_per_m3(
        self,
        solvent_increment_m3: f64,
    ) -> Result<f64, OsmoticNetworkError> {
        let tangent = self
            .at(solvent_increment_m3)?
            .endpoint_tangent_pa()
            .map_err(OsmoticNetworkError::Mechanical)?
            / self.reference_material_volume_m3;
        if positive(tangent) {
            Ok(tangent)
        } else {
            Err(OsmoticNetworkError::Invalid(
                "positive nonlinear swelling stiffness is not representable",
            ))
        }
    }

    pub(super) fn energy_change_j(
        self,
        before_increment_m3: f64,
        after_increment_m3: f64,
    ) -> Result<f64, OsmoticNetworkError> {
        let before = self.at(before_increment_m3)?;
        let after = self.at(after_increment_m3)?;
        let change = before
            .energy_increment_to(&after)
            .map_err(OsmoticNetworkError::Mechanical)?;
        finite(self.reference_material_volume_m3 * change)
    }
}

impl OsmoticPressureVolume {
    /// Retain an original convex material bulk law as a nonlinear swelling load.
    ///
    /// `reference_solvent_volume_m3` is the initial solvent inventory and
    /// `fixed_non_solvent_volume_m3` is an explicit, nonnegative retained volume
    /// (for example the dry material). Their sum maps to the original point's
    /// accepted volume ratio. The background pressure excludes that point's
    /// elastic pressure: a prestrained point retains its original prestress.
    ///
    /// The point preserves the selected quadratic, log-squared or Simo--Taylor
    /// energy and its explicit convex volume domain, including any held natural
    /// reference. Every water trial rechecks that domain. All evaluations are
    /// direct and allocation-free on success; the material's scalar pressure
    /// inverse is never run by this adapter. This supplies no material calibration
    /// or spatial stress-field claim.
    pub fn from_bulk_material(
        reference_solvent_volume_m3: f64,
        fixed_non_solvent_volume_m3: f64,
        background_pressure_pa: f64,
        point: DiscreteBulkPoint,
    ) -> Result<Self, OsmoticNetworkError> {
        if !positive(reference_solvent_volume_m3)
            || !nonnegative(fixed_non_solvent_volume_m3)
            || !background_pressure_pa.is_finite()
        {
            return Err(OsmoticNetworkError::Invalid(
                "bulk swelling needs positive solvent, nonnegative fixed volume and finite background pressure",
            ));
        }
        let total = finite(reference_solvent_volume_m3 + fixed_non_solvent_volume_m3)?;
        let reference_material_volume_m3 = total / point.accepted_volume_ratio();
        if !positive(reference_material_volume_m3) {
            return Err(OsmoticNetworkError::Invalid(
                "bulk swelling reference volume is not representable",
            ));
        }
        let bulk = BulkPressure {
            point,
            reference_material_volume_m3,
        };
        let law = Self {
            reference_volume_m3: reference_solvent_volume_m3,
            reference_pressure_pa: background_pressure_pa,
            elastance_pa_per_m3: bulk.tangent_pa_per_m3(0.0)?,
            bulk: Some(bulk),
        };
        law.pressure_pa(reference_solvent_volume_m3)?;
        Ok(law)
    }

    /// Actual differential pressure/volume stiffness [Pa/m³] at a water trial.
    /// For a material law this is the original positive endpoint tangent,
    /// including the chain rule through the material reference volume.
    pub fn tangent_pa_per_m3(self, volume_m3: f64) -> Result<f64, OsmoticNetworkError> {
        if !positive(volume_m3) {
            return Err(OsmoticNetworkError::Invalid(
                "pressure tangent needs a finite positive solvent volume",
            ));
        }
        match self.bulk {
            Some(bulk) => bulk.tangent_pa_per_m3(volume_m3 - self.reference_volume_m3),
            None => Ok(self.elastance_pa_per_m3),
        }
    }

    /// Mechanical pressure response excluding the arbitrary background datum.
    pub(super) fn response_pressure_pa(self, volume_m3: f64) -> Result<f64, OsmoticNetworkError> {
        if !positive(volume_m3) {
            return Err(OsmoticNetworkError::Invalid(
                "mechanical response needs a finite positive solvent volume",
            ));
        }
        let delta = volume_m3 - self.reference_volume_m3;
        match self.bulk {
            Some(bulk) => bulk.pressure_pa(delta),
            None => finite(self.elastance_pa_per_m3 * delta),
        }
    }
}

/// Exact closed-pair load work, including nonlinear material storage.
/// Common pressure work cancels against the separately checked water balance.
/// A pressure trapezoid would be exact only for affine chamber laws.
pub(super) fn pressure_work_into_solution_j(
    before: &[OsmoticCompartment; 2],
    after_volumes_m3: [f64; 2],
) -> Result<f64, OsmoticNetworkError> {
    let common = before[0].pressure.reference_pressure_pa;
    let mut change = 0.0;
    let mut correction = 0.0;
    for (cell, volume) in before.iter().zip(after_volumes_m3) {
        super::compensated_add(
            &mut change,
            &mut correction,
            cell.pressure.reference_cancelled_potential_change_j(
                cell.solution().solvent_volume_m3(),
                volume,
                common,
            )?,
        );
    }
    finite(-change)
}

#[cfg(test)]
mod tests;
