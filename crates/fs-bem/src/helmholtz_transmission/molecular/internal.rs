//! Finite interior sources through the original molecular transmission factor.
//!
//! The native interior-source equation is unchanged: M delta = -C x_source
//! + [0, (rho_i-rho_o)/rho_o * w_source], with w=i*rho_o*c_eq,o*v_n.
//! Source fields use the original interior complex k, outgoing fields use the
//! original exterior complex k, and both use REAL rho*omega momentum.
//! Reflection is solved directly, never recovered by subtracting total fields.

use super::*;
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::helmholtz::transmission::internal::{
    InternalSourceControls, InternalSourceRegion, sample_reference, solve_correction,
    validate_sources,
};
use crate::near_field::{Geometry, VectorObservation, complex::Prepared};

impl<'s> MolecularFluidInterface<'s> {
    /// Radiate prescribed volume, force and stress sources strictly INSIDE the
    /// original interface. The existing material evaluation and factorization
    /// are reused. Every source is evaluated by the native analytic multipole
    /// owner at this frequency with the interior molecular response.
    ///
    /// Plane waves and unresolved source/interface separations refuse through
    /// the original internal-source admission. The explicit source-point cap
    /// applies before boundary sampling and before each analytic receiver sum.
    /// All sources remain prescribed: no source feedback, motion or finite
    /// source size is supplied by this linear inviscid-momentum model.
    pub fn radiate_internal_field<'p>(
        &'p self,
        sources: IncidentField<'p>,
        controls: InternalSourceControls,
        maximum_source_point_products: u64,
        gate: &CancelGate,
    ) -> Result<MolecularInternalSourceRadiation<'p, 's>, HelmholtzError> {
        poll(gate)?;
        validate_sources(sources, controls)?;
        let n = self.surface.areas().len();
        admit_products(sources.source_count(), n, maximum_source_point_products)?;
        Geometry::new_interior_cancellable(
            self.surface,
            &sources.source_positions(),
            controls.minimum_source_clearance_m,
            gate,
        )?;
        let exterior = self.exterior.gas();
        let interior = self.interior.gas();
        let reference_impedance = exterior.density * exterior.sound_speed;
        let reference = sample_reference(
            self.surface,
            sources,
            self.k_i,
            self.omega * interior.density,
            None,
            reference_impedance,
            gate,
        )?;
        let (correction, relative_backward_residual) = solve_correction(
            &self.contrast,
            self.factor.as_ref(),
            &reference,
            [exterior.density, interior.density],
            controls.maximum_relative_residual,
            gate,
        )?;
        let reflected_velocity: Vec<_> = correction[n..]
            .iter()
            .map(|&w| C64::new(w.im / reference_impedance, -w.re / reference_impedance))
            .collect();
        let total_pressure: Vec<_> = reference[..n]
            .iter()
            .zip(&correction[..n])
            .map(|(&source, &reflected)| source + reflected)
            .collect();
        let total_velocity: Vec<_> = reference[n..]
            .iter()
            .zip(&correction[n..])
            .map(|(&source, &reflected)| {
                let w = source + reflected;
                C64::new(w.im / reference_impedance, -w.re / reference_impedance)
            })
            .collect();
        if total_pressure
            .iter()
            .chain(&total_velocity)
            .chain(&reflected_velocity)
            .any(|&value| !finite(value))
        {
            return Err(bad("nonfinite internal molecular transmission traces"));
        }
        let (_, interface_work_roundoff_interval, interface_work) = evaluate_radiated_power(
            &total_pressure,
            &total_velocity,
            self.surface.areas(),
            Medium::from(exterior),
        );
        if !finite(interface_work.surface_power) || !interface_work.apparent_power.is_finite() {
            return Err(bad("nonfinite internal molecular transmission work"));
        }
        poll(gate)?;
        Ok(MolecularInternalSourceRadiation {
            interface: self,
            sources,
            maximum_source_point_products,
            total_pressure,
            total_velocity,
            reflected_pressure: correction[..n].to_vec(),
            reflected_velocity,
            interface_work,
            interface_work_roundoff_interval,
            relative_backward_residual,
        })
    }
}

/// Original interior emitters, continuous boundary traces, and directly
/// solved regular reflection. The original gas owners remain borrowed through
/// the interface. No independent attenuation or source spectrum is fitted.
pub struct MolecularInternalSourceRadiation<'p, 's> {
    interface: &'p MolecularFluidInterface<'s>,
    sources: IncidentField<'p>,
    maximum_source_point_products: u64,
    total_pressure: Vec<C64>,
    total_velocity: Vec<C64>,
    reflected_pressure: Vec<C64>,
    reflected_velocity: Vec<C64>,
    /// Signed outward work at the ORIGINAL transmitting surface. This volume
    /// contains emitters, so negative work is NOT a source-free absorption
    /// formula. Exterior propagation may further dissipate power: this is not
    /// a lossless far-field radiation claim. Numerical defects remain signed.
    pub interface_work: RadiationPowerDiagnostics,
    /// Panel-dot roundoff only, excluding source, mesh and constitutive errors.
    pub interface_work_roundoff_interval: (f64, f64),
    /// Native corrected-system backward error, excluding spatial/model error.
    pub relative_backward_residual: f64,
}

impl MolecularInternalSourceRadiation<'_, '_> {
    /// Continuous physical pressure, also supporting the exterior field.
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.total_pressure
    }
    /// One actual outward normal velocity shared by the two original gases.
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.total_velocity
    }
    pub fn reflected_pressure_pa(&self) -> &[C64] {
        &self.reflected_pressure
    }
    pub fn reflected_normal_velocity_m_s(&self) -> &[C64] {
        &self.reflected_velocity
    }

    /// Full pressure AND particle velocity in a declared physical region.
    /// Adaptive receivers integrate the same molecular pressure/gradient Green
    /// kernels as the original scattering owner. Direct emitters are added only
    /// in their own interior medium; exactly matched media retain their analytic
    /// fields in both regions after the original geometry admission.
    ///
    /// Quadrature discrepancies cover the fixed boundary representation only,
    /// not source/mesh/material error or relative error near coherent nulls.
    pub fn pressure_velocity_at(
        &self,
        region: InternalSourceRegion,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        self.observe(region, true, points_m, minimum_clearance_m, options, gate)
    }

    /// Observe weak regular interior reflection separately from the direct
    /// source, including the complete analytic velocity at pressure nulls.
    pub fn reflected_pressure_velocity_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        self.observe(
            InternalSourceRegion::Interior,
            false,
            points_m,
            minimum_clearance_m,
            options,
            gate,
        )
    }

    fn observe(
        &self,
        region: InternalSourceRegion,
        total: bool,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        let owner = self.interface;
        let transparent = owner.factor.is_none();
        let interior = region == InternalSourceRegion::Interior;
        let add_source = total && (interior || transparent);
        if add_source {
            admit_products(
                self.sources.source_count(),
                points.len(),
                self.maximum_source_point_products,
            )?;
        }
        options.admit(owner.surface.areas().len(), points.len())?;
        let geometry = if interior {
            Geometry::new_interior_cancellable(owner.surface, points, minimum_clearance_m, gate)?
        } else {
            Geometry::new_iterative(owner.surface, points, minimum_clearance_m, gate)?
        };
        let mut observation = if transparent {
            adaptive::zero_observation(points.len())
        } else {
            let (wave, density, pressure, velocity) = if interior {
                (
                    owner.k_i,
                    owner.interior.gas().density,
                    &self.reflected_pressure,
                    &self.reflected_velocity,
                )
            } else {
                (
                    owner.k_o,
                    owner.exterior.gas().density,
                    &self.total_pressure,
                    &self.total_velocity,
                )
            };
            let prepared = Prepared::new(
                &geometry,
                wave,
                density * owner.omega,
                options.quadrature,
                &mut || gate.is_requested(),
            )?;
            let mut field = prepared.evaluate(pressure, velocity, &mut || gate.is_requested())?;
            if interior {
                // Interior regular reflection is S_i q_i-D_i p. Reverse both
                // p and u once; discrepancy magnitudes remain nonnegative.
                adaptive::reverse_representation(&mut field, gate)?;
            }
            field
        };
        if add_source {
            adaptive::add_incident(
                self.sources,
                owner.k_i,
                owner.interior.gas().density * owner.omega,
                points,
                &mut observation,
                transparent,
                gate,
            )?;
        }
        poll(gate)?;
        Ok(observation)
    }
}

#[cfg(test)]
#[path = "internal/tests.rs"]
mod tests;
