//! Adaptive receivers preserve local outgoing and total interface traces.
use super::*;
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::near_field::{Geometry, Observation, VectorObservation, complex::Prepared};

impl PorousMultipleField<'_, '_> {
    /// Exterior scattered paired p/u from one global Green row bank. LOCAL
    /// outgoing coefficients prevent repeated integration of neighbors' fields.
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
        molecular::admit_products(panels, points.len(), maximum_panel_point_products)?;
        options.admit(panels, points.len())?;
        let geometry = Geometry::new(&self.owner.surface, points, minimum_clearance_m)?;
        self.adaptive_outgoing(&geometry, options, gate)
    }

    /// Coherent total exterior p/u. Intensity uses Re(p*conj(u))/2 after the
    /// incident and scattered amplitudes are added, retaining interference.
    /// The discrepancy estimates cover quadrature only, not solve/model error.
    pub fn total_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        molecular::admit_products(
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

    /// Local p/u inside the selected original body. Both q=i*omega*rho*u and
    /// u=grad(p)/(i*omega*rho) retain that body's actual complex porous inertia.
    /// Pressure and superficial velocity match the exterior gross-area traces.
    ///
    /// A transparent fluid body receives incident plus every active neighbor.
    /// Its global bank therefore admits ALL panels. Geometry must identify the
    /// selected containing body; wrong-body and exterior receivers are refused.
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
            .ok_or_else(|| bad("unknown porous/fluid body"))?;
        let transparent = body.contrast.is_empty();
        let panels = if transparent {
            self.owner.surface.areas().len()
        } else {
            body.original.surface.areas().len()
        };
        molecular::admit_products(panels, points.len(), maximum_panel_point_products)?;
        options.admit(panels, points.len())?;
        if transparent {
            molecular::admit_products(
                self.incident.source_count(),
                points.len(),
                self.maximum_source_point_products,
            )?;
        }
        let geometry = Geometry::new_interior(body.original.surface, points, minimum_clearance_m)?;
        poll(gate)?;
        if transparent {
            drop(geometry);
            let combined =
                Geometry::new_interior(&self.owner.surface, points, minimum_clearance_m)?;
            let mut result = self.adaptive_outgoing(&combined, options, gate)?;
            self.adaptive_add_incident(points, &mut result, gate)?;
            return Ok(result);
        }
        let prepared = Prepared::new_equivalent_fluid(
            &geometry,
            body.wave,
            body.density.scale(self.owner.omega),
            options.quadrature,
            &mut || gate.is_requested(),
        )?;
        let range = body.first_panel..body.first_panel + panels;
        let mut result = prepared.evaluate(
            &self.traces.total_pressure[range.clone()],
            &self.traces.total_velocity[range],
            &mut || gate.is_requested(),
        )?;
        for (pressure, velocity) in result
            .scalar
            .pressure
            .iter_mut()
            .zip(&mut result.particle_velocity_m_s)
        {
            poll(gate)?;
            *pressure = -*pressure;
            for component in velocity {
                *component = -*component;
            }
        }
        poll(gate)?;
        Ok(result)
    }

    /// Direct heating from the original selected body's coherent p/u:
    /// omega/2 [Im(rho)|u|² + Im(1/K)|p|²], per total body volume.
    /// Lossless fluid bodies return zero; free molecular gas has compression
    /// loss only. A porous body's compression already includes coupled pore
    /// thermal/molecular loss. No boundary residual supplies either channel.
    /// Output buffers are excluded from the original prepared-row byte budget.
    pub fn interior_absorption_adaptive_at(
        &self,
        body: usize,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<PorousAbsorptionObservation, HelmholtzError> {
        let fields = self.interior_pressure_velocity_adaptive_at(
            body,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )?;
        let material = &self.owner.bodies[body];
        super::super::receivers::absorption(
            fields,
            material.density,
            material.bulk,
            self.owner.omega,
            gate,
        )
    }

    fn adaptive_outgoing(
        &self,
        geometry: &Geometry<'_>,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        if self.owner.factor.is_none() {
            let n = geometry.points_m().len();
            let result = VectorObservation {
                scalar: Observation {
                    pressure: vec![C64::ZERO; n],
                    quadrature_error_estimate_pa: vec![0.0; n],
                },
                particle_velocity_m_s: vec![[C64::ZERO; 3]; n],
                quadrature_error_estimate_m_s: vec![[0.0; 3]; n],
            };
            poll(gate)?;
            return Ok(result);
        }
        let wave = DampedWaveNumber {
            k: C64::from_re(self.owner.k_o),
            phase_origin_m: 0.0,
        };
        let prepared = Prepared::new(
            geometry,
            wave,
            self.owner.exterior.density_kg_m3() * self.owner.omega,
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
        for ((&point, pressure), velocity) in points
            .iter()
            .zip(&mut observation.scalar.pressure)
            .zip(&mut observation.particle_velocity_m_s)
        {
            poll(gate)?;
            let source = self.incident.sample_wave(
                self.owner.k_o,
                self.owner.exterior.density_kg_m3() * self.owner.omega,
                None,
                point,
                &mut || gate.is_requested(),
            )?;
            if self.owner.factor.is_none() {
                *pressure = source.pressure_pa;
                *velocity = source.velocity_m_s;
            } else {
                *pressure = *pressure + source.pressure_pa;
                for (component, incoming) in velocity.iter_mut().zip(source.velocity_m_s) {
                    *component = *component + incoming;
                }
            }
            if !finite(*pressure) || velocity.iter().any(|&v| !finite(v)) {
                return Err(bad("nonfinite porous multiple-body total receiver"));
            }
        }
        poll(gate)
    }
}
