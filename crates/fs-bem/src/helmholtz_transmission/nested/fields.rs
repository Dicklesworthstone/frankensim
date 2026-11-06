//! Region-correct pressure and particle velocity from the original observers.
use super::*;
use crate::near_field::{Geometry, Observation, VectorObservation};

impl CoreShellScatteredField<'_, '_> {
    /// Observe the complete field in one declared physical region. Winding
    /// and face/edge distance enforce the region using the original meshes.
    /// Shell observations use BOTH original boundaries with opposite signs,
    /// the shell's actual k/rho, and one shared quadrature-work allowance.
    ///
    /// The existing adaptive triangle owner supplies pressure and analytic
    /// particle velocity together. Returned discrepancy estimates add absolute
    /// contributions from both surfaces; they do not bound mesh/model error
    /// or guarantee relative accuracy of a nearly cancelling total field.
    pub fn pressure_velocity_at(
        &self,
        region: CoreShellRegion,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        options.validate()?;
        let owner = self.interface;
        let outer = match region {
            CoreShellRegion::Exterior => {
                Geometry::new_iterative(owner.surfaces[0], points_m, minimum_clearance_m, gate)?
            }
            CoreShellRegion::Shell | CoreShellRegion::Core => {
                Geometry::new_interior(owner.surfaces[0], points_m, minimum_clearance_m)?
            }
        };
        poll(gate)?;
        let core = match region {
            CoreShellRegion::Core => {
                Geometry::new_interior(owner.surfaces[1], points_m, minimum_clearance_m)?
            }
            CoreShellRegion::Exterior | CoreShellRegion::Shell => {
                Geometry::new_iterative(owner.surfaces[1], points_m, minimum_clearance_m, gate)?
            }
        };
        poll(gate)?;
        let mut remaining = options.maximum_kernel_evaluations;
        let mut result = zero_observation(points_m.len());
        let has_outer = owner.offsets[0].is_some();
        let has_core = owner.offsets[1].is_some();
        let mut add_incident = false;
        match region {
            CoreShellRegion::Exterior => {
                if has_outer || has_core {
                    let boundary = owner.radiating_interface();
                    let geometry = if boundary == 0 { &outer } else { &core };
                    result = self.observe_geometry(
                        geometry,
                        boundary,
                        0,
                        &self.scattered.pressure,
                        &self.scattered.velocity,
                        1.0,
                        options,
                        &mut remaining,
                        gate,
                    )?;
                }
                add_incident = true;
            }
            CoreShellRegion::Shell if has_outer => {
                let trace = self.traces[0].as_ref().expect("active outer interface");
                result = self.observe_geometry(
                    &outer,
                    0,
                    1,
                    &trace.pressure,
                    &trace.velocity,
                    -1.0,
                    options,
                    &mut remaining,
                    gate,
                )?;
                if has_core {
                    let trace = self.traces[1].as_ref().expect("active core interface");
                    let contribution = self.observe_geometry(
                        &core,
                        1,
                        1,
                        &trace.pressure,
                        &trace.velocity,
                        1.0,
                        options,
                        &mut remaining,
                        gate,
                    )?;
                    accumulate(&mut result, contribution)?;
                }
            }
            CoreShellRegion::Shell => {
                // With a transparent outer interface the shell is the original
                // exterior medium; retain its analytic incident field exactly.
                if has_core {
                    result = self.observe_geometry(
                        &core,
                        1,
                        0,
                        &self.scattered.pressure,
                        &self.scattered.velocity,
                        1.0,
                        options,
                        &mut remaining,
                        gate,
                    )?;
                }
                add_incident = true;
            }
            CoreShellRegion::Core if has_core => {
                let trace = self.traces[1].as_ref().expect("active core interface");
                result = self.observe_geometry(
                    &core,
                    1,
                    2,
                    &trace.pressure,
                    &trace.velocity,
                    -1.0,
                    options,
                    &mut remaining,
                    gate,
                )?;
            }
            CoreShellRegion::Core if has_outer => {
                // Removing a matched core boundary leaves one homogeneous
                // interior. Do not reconstruct an artificial second surface.
                let trace = self.traces[0].as_ref().expect("active outer interface");
                result = self.observe_geometry(
                    &outer,
                    0,
                    2,
                    &trace.pressure,
                    &trace.velocity,
                    -1.0,
                    options,
                    &mut remaining,
                    gate,
                )?;
            }
            CoreShellRegion::Core => add_incident = true,
        }
        if add_incident {
            for (index, &point) in points_m.iter().enumerate() {
                poll(gate)?;
                let sample = self.incident.pressure_velocity_at(
                    owner.k[0],
                    medium(owner.states[0]),
                    point,
                )?;
                result.scalar.pressure[index] = result.scalar.pressure[index] + sample.pressure_pa;
                for axis in 0..3 {
                    result.particle_velocity_m_s[index][axis] =
                        result.particle_velocity_m_s[index][axis] + sample.velocity_m_s[axis];
                }
            }
        }
        validate_observation(&result)?;
        poll(gate)?;
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_geometry(
        &self,
        geometry: &Geometry<'_>,
        boundary: usize,
        region: usize,
        pressure: &[C64],
        velocity: &[C64],
        sign: f64,
        mut options: near_field::Options,
        remaining: &mut usize,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        if *remaining < 80 {
            return Err(bad(
                "core-shell total receiver kernel-evaluation budget exhausted",
            ));
        }
        options.maximum_kernel_evaluations = *remaining;
        let owner = self.interface;
        let representation = owner.radiation(
            boundary,
            region,
            pressure.iter().map(|p| p.scale(sign)).collect(),
            velocity.iter().map(|v| v.scale(sign)).collect(),
        )?;
        let prepared = geometry.prepare_velocity_cancellable(
            owner.k[region],
            medium(owner.states[region]),
            options,
            gate,
        )?;
        *remaining -= prepared.kernel_evaluations;
        let result = prepared.evaluate_velocity(&representation)?;
        poll(gate)?;
        Ok(result)
    }
}

fn zero_observation(count: usize) -> VectorObservation {
    VectorObservation {
        scalar: Observation {
            pressure: vec![C64::ZERO; count],
            quadrature_error_estimate_pa: vec![0.0; count],
        },
        particle_velocity_m_s: vec![[C64::ZERO; 3]; count],
        quadrature_error_estimate_m_s: vec![[0.0; 3]; count],
    }
}

fn accumulate(
    target: &mut VectorObservation,
    other: VectorObservation,
) -> Result<(), HelmholtzError> {
    for index in 0..target.scalar.pressure.len() {
        target.scalar.pressure[index] =
            target.scalar.pressure[index] + other.scalar.pressure[index];
        target.scalar.quadrature_error_estimate_pa[index] +=
            other.scalar.quadrature_error_estimate_pa[index];
        for axis in 0..3 {
            target.particle_velocity_m_s[index][axis] = target.particle_velocity_m_s[index][axis]
                + other.particle_velocity_m_s[index][axis];
            target.quadrature_error_estimate_m_s[index][axis] +=
                other.quadrature_error_estimate_m_s[index][axis];
        }
    }
    validate_observation(target)
}

fn validate_observation(field: &VectorObservation) -> Result<(), HelmholtzError> {
    if field
        .scalar
        .pressure
        .iter()
        .chain(field.particle_velocity_m_s.iter().flatten())
        .any(|&value| !finite(value))
        || field
            .scalar
            .quadrature_error_estimate_pa
            .iter()
            .chain(field.quadrature_error_estimate_m_s.iter().flatten())
            .any(|error| !error.is_finite() || *error < 0.0)
    {
        Err(bad("nonfinite core-shell receiver field or discrepancy"))
    } else {
        Ok(())
    }
}
