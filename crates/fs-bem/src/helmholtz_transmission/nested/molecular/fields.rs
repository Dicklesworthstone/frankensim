//! Adaptive complex Green fields on the same physical three-region topology.
use super::*;
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::near_field::{VectorObservation, complex::Prepared};

impl MolecularCoreShellScatteredField<'_, '_> {
    /// Observe continuous pressure and actual particle velocity in a declared
    /// region, with the original molecular k and REAL rho*omega there. Shell
    /// fields combine both boundaries with opposite signs; multiple reflection
    /// is already included in their common solve. No sound-speed refit occurs.
    ///
    /// Both retained meshes enforce region/clearance even when an interface
    /// is transparent. One quadrature allowance covers all contributions;
    /// row-byte admission conservatively includes all contributing surfaces.
    /// Absolute discrepancies add across surfaces and exclude mesh/model error.
    pub fn pressure_velocity_at(
        &self,
        region: CoreShellRegion,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        options.quadrature.validate()?;
        let owner = self.interface;
        let (contributions, add_incident) = fields::region_plan(region, owner.offsets);
        if add_incident {
            admit_products(
                self.incident.source_count(), points_m.len(), self.maximum_source_point_products,
            )?;
        }
        let panels = contributions.iter().flatten()
            .map(|term| owner.surfaces[term.boundary].areas().len()).sum();
        options.admit(panels, points_m.len())?;
        let geometries = fields::region_geometry(
            owner.surfaces, region, points_m, minimum_clearance_m, gate,
        )?;
        let mut remaining = options.quadrature.maximum_kernel_evaluations;
        let mut result = None;
        for contribution in contributions.into_iter().flatten() {
            poll(gate)?;
            let (pressure, velocity) = if contribution.scattered {
                (&self.scattered_pressure, &self.scattered_velocity)
            } else {
                let trace = self.traces[contribution.boundary].as_ref()
                    .expect("active molecular interface");
                (&trace.pressure, &trace.velocity)
            };
            let mut quadrature = options.quadrature;
            quadrature.maximum_kernel_evaluations = remaining;
            let prepared = Prepared::new(
                &geometries[contribution.boundary],
                owner.waves[contribution.region],
                owner.sources[contribution.region].gas().density * owner.omega,
                quadrature,
                &mut || gate.is_requested(),
            )?;
            remaining -= prepared.kernel_evaluations();
            let mut observed = prepared.evaluate(pressure, velocity, &mut || gate.is_requested())?;
            if contribution.sign < 0.0 {
                for value in &mut observed.scalar.pressure {
                    *value = -*value;
                }
                for velocity in &mut observed.particle_velocity_m_s {
                    for value in velocity {
                        *value = -*value;
                    }
                }
            }
            if let Some(total) = &mut result {
                fields::accumulate(total, observed)?;
            } else {
                result = Some(observed);
            }
        }
        let transparent = result.is_none();
        let mut result = result.unwrap_or_else(|| fields::zero_observation(points_m.len()));
        if add_incident {
            for (index, &point) in points_m.iter().enumerate() {
                poll(gate)?;
                let sample = self.incident.sample_wave(
                    owner.waves[0], owner.sources[0].gas().density * owner.omega,
                    None, point, &mut || gate.is_requested(),
                )?;
                if transparent {
                    // Keep exact original fields, including signed zero, when
                    // there is no contrast anywhere along this representation.
                    result.scalar.pressure[index] = sample.pressure_pa;
                    result.particle_velocity_m_s[index] = sample.velocity_m_s;
                } else {
                    result.scalar.pressure[index] =
                        result.scalar.pressure[index] + sample.pressure_pa;
                    for axis in 0..3 {
                        result.particle_velocity_m_s[index][axis] =
                            result.particle_velocity_m_s[index][axis] + sample.velocity_m_s[axis];
                    }
                }
            }
        }
        fields::validate_observation(&result)?;
        poll(gate)?;
        Ok(result)
    }
}
