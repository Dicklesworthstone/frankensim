//! Reusable analytic illumination in one original evaluated molecular gas.
//!
//! Material evaluation and source validation occur once before quadrature.
//! Every later sample keeps the same complex k and real physical rho*omega,
//! using the existing plane-wave, volume-source and force-source kernels.
use super::*;
use crate::helmholtz::bc::incident::{IncidentField, IncidentSample};

/// Immutable source-bound incident field at one positive real frequency.
/// Finite emitter declarations remain borrowed; pressure/force/volume-flow
/// phases cannot be changed underneath a prepared source. Source-to-body
/// clearance belongs to the boundary owner, not this free-space sampler.
#[derive(Clone)]
pub struct PreparedMolecularIncident<'i> {
    source: MolecularExterior,
    field: IncidentField<'i>,
    omega: f64,
    wave: DampedWaveNumber,
    rho_omega: f64,
}

impl MolecularExterior {
    /// Validate the complete coherent source and evaluate the original gas
    /// once. Reuse the resulting sampler for boundary integration or finite
    /// receivers without repeating constitutive evaluation at every point.
    ///
    /// This preserves the original harmonic/source domain: a stationary
    /// homogeneous gas, ideal plane waves and compact volume/force sources.
    /// Plane amplitude is pressure at the original coordinate origin. An
    /// absorbing plane wave grows upstream; unrepresentable samples refuse.
    pub fn prepare_incident<'i>(
        &self,
        omega: f64,
        field: IncidentField<'i>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularIncident<'i>, MolecularError> {
        poll(&mut cancelled)?;
        field.validate()?;
        let wave = self.incident_wave(omega)?;
        let rho_omega = self.source.gas().density * omega;
        poll(&mut cancelled)?;
        Ok(PreparedMolecularIncident {
            source: self.clone(),
            field,
            omega,
            wave,
            rho_omega,
        })
    }
}

impl PreparedMolecularIncident<'_> {
    pub fn source(&self) -> &MolecularExterior {
        &self.source
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn wavenumber_per_m(&self) -> C64 {
        self.wave.k
    }
    /// All coherent terms, including plane waves, for whole-query work admission.
    pub fn source_count(&self) -> usize {
        self.field.source_count()
    }
    /// Finite emitter positions only, in their original order. Each validated
    /// field contains at most the existing MAX_INCIDENT_SOURCES terms.
    pub fn source_positions(&self) -> Vec<[f64; 3]> {
        self.field.source_positions()
    }

    /// Bind another complete illumination to this SAME gas and frequency.
    /// Only the new source declaration is validated; original k and rho*omega
    /// are retained without a second constitutive evaluation. This supports
    /// shared boundary factorizations over independent incident cases.
    pub fn with_field<'j>(
        &self,
        field: IncidentField<'j>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularIncident<'j>, HelmholtzError> {
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        field.validate()?;
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(PreparedMolecularIncident {
            source: self.source.clone(),
            field,
            omega: self.omega,
            wave: self.wave,
            rho_omega: self.rho_omega,
        })
    }

    /// Original coherent pressure [Pa] and Cartesian Euler velocity [m/s].
    /// This performs at most the existing bounded source count of analytic
    /// kernel evaluations; consumers integrating many points admit their
    /// aggregate work separately. There is no new material evaluation.
    pub fn pressure_velocity_at(
        &self,
        point_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IncidentSample, HelmholtzError> {
        let out =
            self.field
                .sample_wave(self.wave, self.rho_omega, None, point_m, &mut cancelled)?;
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(out)
    }

    /// Complete point batch in original order, admitted before allocation.
    /// `maximum_source_point_products` bounds source_count*points.len(); it is
    /// a logical work limit, not a timing or spatial-accuracy certificate.
    /// The returned payload is points.len()*size_of::<IncidentSample>() bytes,
    /// excluding allocator metadata. An empty list returns an empty vector
    /// after cancellation and positive-budget admission. Any failed sample
    /// refuses the batch without publishing partial observations.
    pub fn pressure_velocity_at_many(
        &self,
        points_m: &[[f64; 3]],
        maximum_source_point_products: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        let work = self.source_count() as u128 * points_m.len() as u128;
        if maximum_source_point_products == 0
            || work > u128::from(maximum_source_point_products)
            || points_m
                .len()
                .checked_mul(core::mem::size_of::<IncidentSample>())
                .is_none()
        {
            return Err(HelmholtzError::BadParameter {
                what: "molecular incident batch exceeds its source-point or payload extent",
            });
        }
        let mut out = Vec::new();
        out.try_reserve_exact(points_m.len())
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "prepared molecular incident point samples",
            })?;
        for &point in points_m {
            out.push(self.pressure_velocity_at(point, &mut cancelled)?);
        }
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests;
