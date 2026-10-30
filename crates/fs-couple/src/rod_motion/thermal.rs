//! Source-driven temperature changes without resetting a native rod trajectory.
//!
//! The thermal field is externally prescribed. Its fixed-pose change of native
//! strain storage is parameter work; this is not a sensible/latent heat account
//! or a finite-reservoir thermodynamic balance. Subsequent motion uses the same
//! original analytical forces, reference inertia and Hamiltonian splitting.

use std::sync::Arc;

use fs_exec::CancelGate;
use fs_solid::rod::material_field::thermal::source::RodMaterialSource;

use super::{RodMotion, RodMotionError, RodMotionInertia, checkpoint, finite};

/// One accepted-time section-temperature change with zero instantaneous impulse.
#[derive(Clone, Debug, PartialEq)]
pub struct RodThermalRefresh {
    /// Unchanged physical clock [s]; no integration step is consumed.
    pub time_s: f64,
    /// Old source temperatures in segment/layer order [K].
    pub temperatures_before_k: Vec<f64>,
    /// Admitted new source temperatures in the same order [K].
    pub temperatures_after_k: Vec<f64>,
    /// Native section and attached-link strain storage at the original pose [J].
    pub strain_energy_before_j: f64,
    /// Native strain storage at the identical pose with refreshed sections [J].
    pub strain_energy_after_j: f64,
    /// `strain_energy_after_j - strain_energy_before_j` [J]. Positive is energy
    /// entering the mechanical storage through the prescribed thermal parameter.
    /// Neither a heat-capacity estimate nor a complete thermodynamic heat debit.
    pub thermal_parameter_work_j: f64,
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Refresh section E/rho/G and integrated expansion from the original pinned
    /// material source at the accepted clock. Preserve actual pose, momenta,
    /// original reference mass/inertia, supports, body loads and step count.
    /// Fresh native forces replace the previous response atomically.
    ///
    /// Temperature is a supplied body field, not an ambient reservoir value.
    /// Its fixed-pose strain-energy change is reported separately from later
    /// mechanical work. Fixed section geometry/reference displaced volumes and
    /// separately bound tension-link material states remain the active model.
    /// This path adds no heat solver, thermoelastic heat feedback, phase change,
    /// mass transport, damping or cross-section expansion approximation.
    ///
    /// # Errors
    /// Source/reference/geometry mismatch, unsupported temperature paths, native
    /// actual-pose strain admission, arithmetic or cancellation refuse without
    /// changing the accepted trajectory or either set of source observations.
    pub fn rebind_material_temperatures_candidate(
        &self,
        source: &RodMaterialSource<'_>,
        temperatures_k: &[f64],
        gate: &CancelGate,
    ) -> Result<(Self, RodThermalRefresh), RodMotionError> {
        checkpoint(gate)?;
        let temperatures_before_k = self.model.section_temperatures_k()?;
        let strain_energy_before_j = self.model.strain_energy_j(&self.pose, gate)?;
        let model = Arc::new(self.model.rebind_material_temperatures(
            source,
            temperatures_k,
            gate,
        )?);
        let response = model.evaluate(&self.pose, gate)?;
        let strain_energy_after_j = model.strain_energy_j(&self.pose, gate)?;
        let thermal_parameter_work_j = strain_energy_after_j - strain_energy_before_j;
        finite(&[
            strain_energy_before_j,
            strain_energy_after_j,
            thermal_parameter_work_j,
        ])?;
        let report = RodThermalRefresh {
            time_s: self.time_s,
            temperatures_before_k,
            temperatures_after_k: model.section_temperatures_k()?,
            strain_energy_before_j,
            strain_energy_after_j,
            thermal_parameter_work_j,
        };
        let mut candidate = self.clone();
        candidate.model = model;
        candidate.response = response;
        checkpoint(gate)?;
        Ok((candidate, report))
    }
}
