//! Adaptive paired receivers retain all original molecular-body interactions.
use super::super::adaptive::{add_incident, reverse_representation, zero_observation};
use super::*;
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::near_field::{Geometry, VectorObservation, complex::Prepared};

impl MolecularMultipleField<'_, '_> {
    /// Adaptive exterior Green integration of LOCAL outgoing coefficients.
    /// Complete physical scattered traces include neighbors already and must
    /// not be used here: doing so would count their fields repeatedly.
    /// One prepared global row bank integrates paired p/u with componentwise
    /// quadrature discrepancies, without repeating the coupled boundary solve.
    pub fn scattered_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        let panels = self.owner.surface.areas().len();
        admit_receivers(panels, points, maximum_panel_point_products, gate)?;
        options.admit(panels, points.len())?;
        let geometry = Geometry::new(&self.owner.surface, points, minimum_clearance_m)?;
        self.adaptive_outgoing(&geometry, options, gate)
    }

    /// Total adaptive exterior p/u plus the original coherent incident field.
    /// Discrepancy estimates cover boundary quadrature only; cancellation of
    /// incident and scattered amplitudes may amplify their relative error.
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
        self.adaptive_add_incident(points, &mut result, gate)?;
        Ok(result)
    }

    /// Paired interior p/u in the selected body's original molecular gas.
    /// Active bodies integrate their TOTAL physical traces with the interior
    /// Green sign. A transparent body receives incident plus every active
    /// neighbor, using one global outgoing row bank and ALL-panel admission.
    /// Geometry still requires the receiver inside the selected body and
    /// outside every other body. No artificial interface source is introduced.
    pub fn interior_pressure_velocity_adaptive_at(
        &self,
        body: usize,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        let body = self
            .owner
            .bodies
            .get(body)
            .ok_or_else(|| bad("unknown molecular body"))?;
        let transparent = body.contrast.is_empty();
        let panels = if transparent {
            self.owner.surface.areas().len()
        } else {
            body.original.surface.areas().len()
        };
        admit_receivers(panels, points, maximum_panel_point_products, gate)?;
        options.admit(panels, points.len())?;
        if transparent {
            admit_products(
                self.incident.source_count(),
                points.len(),
                self.maximum_source_point_products,
            )?;
        }
        let geometry = Geometry::new_interior(body.original.surface, points, minimum_clearance_m)?;
        poll(gate)?;
        if transparent {
            drop(geometry);
            // The selected-body admission above identifies the containing
            // component. Combined admission checks clearance against all bodies
            // and exactly one containing component; the same global rows then
            // integrate zero own coefficients plus active exterior neighbors.
            let combined =
                Geometry::new_interior(&self.owner.surface, points, minimum_clearance_m)?;
            let mut result = self.adaptive_outgoing(&combined, options, gate)?;
            self.adaptive_add_incident(points, &mut result, gate)?;
            return Ok(result);
        }
        let prepared = Prepared::new(
            &geometry,
            body.wave,
            body.original.material.gas().density * self.owner.omega,
            options.quadrature,
            &mut || gate.is_requested(),
        )?;
        let range = body.first_panel..body.first_panel + panels;
        let mut result = prepared.evaluate(
            &self.traces.total_pressure[range.clone()],
            &self.traces.total_velocity[range],
            &mut || gate.is_requested(),
        )?;
        reverse_representation(&mut result, gate)?;
        poll(gate)?;
        Ok(result)
    }

    fn adaptive_outgoing(
        &self,
        geometry: &Geometry<'_>,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        if self.owner.factor.is_none() {
            let result = zero_observation(geometry.points_m().len());
            poll(gate)?;
            return Ok(result);
        }
        let prepared = Prepared::new(
            geometry,
            self.owner.wave,
            self.owner.exterior.gas().density * self.owner.omega,
            options.quadrature,
            &mut || gate.is_requested(),
        )?;
        prepared.evaluate(&self.traces.pressure, &self.traces.velocity, &mut || {
            gate.is_requested()
        })
    }

    fn adaptive_add_incident(
        &self,
        points: &[[f64; 3]],
        observation: &mut VectorObservation,
        gate: &CancelGate,
    ) -> Result<(), HelmholtzError> {
        add_incident(
            self.incident,
            self.owner.wave,
            self.owner.exterior.gas().density * self.owner.omega,
            points,
            observation,
            self.owner.factor.is_none(),
            gate,
        )
    }
}

#[cfg(test)]
#[path = "near_field/tests.rs"]
mod tests;
