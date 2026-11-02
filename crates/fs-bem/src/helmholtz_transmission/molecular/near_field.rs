//! Adaptive receivers for the original molecular transmission solution.
use super::*;
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::near_field::{Geometry, Observation, VectorObservation, complex::Prepared};

impl MolecularTransmissionField<'_, '_> {
    /// Integrate scattered pressure and particle velocity over the original
    /// retained triangles. One adaptive row bank serves all receivers and
    /// both physical traces. Pressure and each velocity component have their
    /// own propagated quadrature discrepancy; these exclude BEM/material error.
    /// No boundary solve or centroid fallback occurs.
    pub fn scattered_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        self.observe_adaptive(
            false,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )
    }

    /// Adaptive scattered field plus the original analytic incident field.
    /// Original source-point admission applies before quadrature. Discrepancies
    /// describe the boundary integral only, not analytic source/model error or
    /// relative accuracy of a nearly cancelling total field.
    pub fn total_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        admit_products(
            self.incident.source_count(),
            points.len(),
            self.maximum_source_point_products,
        )?;
        let mut result = self.scattered_pressure_velocity_adaptive_at(
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )?;
        add_incident(
            self.incident,
            self.interface.k_o,
            self.interface.exterior.gas().density * self.interface.omega,
            points,
            &mut result,
            self.interface.factor.is_none(),
            gate,
        )?;
        Ok(result)
    }

    /// Adaptive interior S_i*q_i-D_i*p with the ORIGINAL interior momentum
    /// rho_i*omega. The representation sign reverses both paired fields once;
    /// discrepancy magnitudes retain their positive sign. Exact transparency
    /// evaluates the original incident field directly after geometry admission.
    pub fn interior_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        self.observe_adaptive(
            true,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )
    }

    fn observe_adaptive(
        &self,
        interior: bool,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        let owner = self.interface;
        let panels = owner.surface.areas().len();
        admit_products(panels, points.len(), maximum_panel_point_products)?;
        options.admit(panels, points.len())?;
        let transparent = owner.factor.is_none();
        if interior && transparent {
            admit_products(
                self.incident.source_count(),
                points.len(),
                self.maximum_source_point_products,
            )?;
        }
        let geometry = if interior {
            Geometry::new_interior(owner.surface, points, minimum_clearance_m)?
        } else {
            Geometry::new(owner.surface, points, minimum_clearance_m)?
        };
        poll(gate)?;
        if transparent {
            let mut result = zero_observation(points.len());
            if interior {
                add_incident(
                    self.incident,
                    owner.k_i,
                    owner.interior.gas().density * owner.omega,
                    points,
                    &mut result,
                    true,
                    gate,
                )?;
            }
            poll(gate)?;
            return Ok(result);
        }
        let (wave, density, p, v) = if interior {
            (
                owner.k_i,
                owner.interior.gas().density,
                self.total_pressure_pa(),
                self.total_normal_velocity_m_s(),
            )
        } else {
            (
                owner.k_o,
                owner.exterior.gas().density,
                self.scattered_pressure_pa(),
                self.scattered_normal_velocity_m_s(),
            )
        };
        let prepared = Prepared::new(
            &geometry,
            wave,
            density * owner.omega,
            options.quadrature,
            &mut || gate.is_requested(),
        )?;
        let mut result = prepared.evaluate(p, v, &mut || gate.is_requested())?;
        if interior {
            reverse_representation(&mut result, gate)?;
        }
        poll(gate)?;
        Ok(result)
    }
}

pub(super) fn zero_observation(receivers: usize) -> VectorObservation {
    VectorObservation {
        scalar: Observation {
            pressure: vec![C64::ZERO; receivers],
            quadrature_error_estimate_pa: vec![0.0; receivers],
        },
        particle_velocity_m_s: vec![[C64::ZERO; 3]; receivers],
        quadrature_error_estimate_m_s: vec![[0.0; 3]; receivers],
    }
}

pub(super) fn reverse_representation(
    observation: &mut VectorObservation,
    gate: &CancelGate,
) -> Result<(), HelmholtzError> {
    for (p, v) in observation
        .scalar
        .pressure
        .iter_mut()
        .zip(&mut observation.particle_velocity_m_s)
    {
        poll(gate)?;
        *p = -*p;
        for component in v {
            *component = -*component;
        }
    }
    poll(gate)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_incident(
    incident: IncidentField<'_>,
    wave: DampedWaveNumber,
    rho_omega: f64,
    points: &[[f64; 3]],
    observation: &mut VectorObservation,
    transparent: bool,
    gate: &CancelGate,
) -> Result<(), HelmholtzError> {
    for ((&point, p), v) in points
        .iter()
        .zip(&mut observation.scalar.pressure)
        .zip(&mut observation.particle_velocity_m_s)
    {
        poll(gate)?;
        let source =
            incident.sample_wave(wave, rho_omega, None, point, &mut || gate.is_requested())?;
        if transparent {
            // Preserve original incident values and exact zero components.
            *p = source.pressure_pa;
            *v = source.velocity_m_s;
        } else {
            *p = *p + source.pressure_pa;
            for (component, incoming) in v.iter_mut().zip(source.velocity_m_s) {
                *component = *component + incoming;
            }
        }
        if !finite(*p) || v.iter().any(|&component| !finite(component)) {
            return Err(bad("nonfinite adaptive molecular total field"));
        }
    }
    poll(gate)
}

#[cfg(test)]
#[path = "near_field/tests.rs"]
mod tests;
