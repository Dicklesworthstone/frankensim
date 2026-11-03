//! Fixed-pipe reaction wrenches from the original hydraulic trajectory.
//!
//! This is a slender, constant-section, plug-flow control volume, consistent
//! with the original `rho length/area` inertance. Pressure traction uses the
//! two original node pressures separately. Stored and advected linear AND
//! angular momentum determine the load on the enclosing pipe and any contained
//! pressure drive. Drag is already present in those pressures; adding another
//! `loss*area` force would count it twice. No shaft torque or swirl is inferred.
//!
//! Geometry is fixed during the native step. The wrench can drive a mechanical
//! load case or a one-way mechanical step; moving pipes require a relative-flow
//! momentum/energy owner. Junction loads retain the native common static
//! pressure; reducer recovery and mixing/thermal closure remain external.
//! All coordinates use one caller-declared inertial
//! Cartesian frame; no rotation or frame conversion is inferred.

use super::{InertialGraphError, InertialGraphTrial, checkpoint};
use fs_exec::Cx;
use fs_qty::Pressure;

mod pipe;
pub use pipe::{PipeLoad, PipeLoadError, PipeLoadGeometry};
mod junction;
pub use junction::JunctionLoad;
mod assembly;
pub use assembly::HydraulicAssemblyLoads;

impl InertialGraphTrial {
    fn admit_pipe(
        &self,
        cx: &Cx<'_>,
        branch: usize,
        geometry: &PipeLoadGeometry,
    ) -> Result<(), InertialGraphError> {
        checkpoint(cx)?;
        let c = &self.initial.compiled;
        let original = c.branches.get(branch).ok_or(InertialGraphError::Invalid(
            "pipe load requires an original branch",
        ))?;
        if geometry.area_m2 != original.flow_area_m2
            || (geometry.length_m - original.length_m).abs()
                > 64.0 * f64::EPSILON * geometry.length_m.max(original.length_m)
        {
            return Err(InertialGraphError::Invalid(
                "pipe load geometry differs from original inertia",
            ));
        }
        let edge = &c.topology.branches()[branch];
        if let Some(gravity) = &self.gravity {
            if geometry.vertices_m[0] != gravity.positions_m()[edge.from]
                || *geometry.vertices_m.last().unwrap() != gravity.positions_m()[edge.to]
                || gravity.branch_sections().is_some_and(|paths| {
                    let sections = &paths[branch];
                    sections.len() + 1 != geometry.vertices_m.len()
                        || sections.iter().zip(geometry.vertices_m.windows(2)).any(
                            |(section, pair)| {
                                section.start_m != pair[0] || section.end_m != pair[1]
                            },
                        )
                })
            {
                return Err(InertialGraphError::Invalid(
                    "pipe load path differs from original gravity",
                ));
            }
        }
        Ok(())
    }

    /// Resolve the fluid's net load on one original pipe, about `about_m`.
    ///
    /// The supplied path runs from the branch's `from` node to its `to` node.
    /// `exterior_pressure` is the uniform surrounding interval-average pressure in the
    /// same datum as the original nodal pressures. It supplies the exterior
    /// pressure traction; an arbitrary gauge offset cannot create a pipe load.
    /// Original gravity densities are retained after constitutive rechecking.
    /// The path must match original length/area and, when bound, gravity's
    /// original endpoints and every distributed section.
    ///
    /// Flux uses the exact integral of Q(t)^2 for the native affine trajectory,
    /// including both halves of a reversal. This is a time reconstruction of
    /// the midpoint method, not a higher-order fluid solution. The returned
    /// candidate load may be used in coupling trials; publish it only after
    /// the SAME trial passes `accept`'s impulse and energy gates.
    pub fn pipe_load(
        &self,
        cx: &Cx<'_>,
        branch: usize,
        geometry: &PipeLoadGeometry,
        exterior_pressure: Pressure,
        about_m: [f64; 3],
    ) -> Result<PipeLoad, InertialGraphError> {
        self.pipe_load_with_ports(cx, branch, geometry, exterior_pressure, about_m, [true; 2])
    }

    fn pipe_load_with_ports(
        &self,
        cx: &Cx<'_>,
        branch: usize,
        geometry: &PipeLoadGeometry,
        exterior_pressure: Pressure,
        about_m: [f64; 3],
        exposed_ports: [bool; 2],
    ) -> Result<PipeLoad, InertialGraphError> {
        self.admit_pipe(cx, branch, geometry)?;
        let c = &self.initial.compiled;
        let edge = &c.topology.branches()[branch];
        let report = &self.report;
        let rho = c.reference_density.value();
        let result = geometry.evaluate(
            pipe::FlowInterval {
                duration_s: report.duration_s,
                volume_flows_m3_s: [
                    report.initial_flows_m3_s[branch],
                    report.final_flows_m3_s[branch],
                ],
                density_kg_m3: rho,
                pressures_pa: [
                    report.average_pressures_pa[edge.from],
                    report.average_pressures_pa[edge.to],
                ],
                exterior_pressure_pa: exterior_pressure.value(),
                exposed_ports,
                gravity_m_s2: self
                    .gravity
                    .as_ref()
                    .map_or([0.0; 3], |g| g.acceleration_m_s2()),
            },
            about_m,
            |section| {
                self.gravity.as_ref().map_or(rho, |gravity| {
                    gravity.branch_materials().map_or_else(
                        || {
                            gravity.branch_sections().unwrap()[branch][section]
                                .transport
                                .density_kg_m3()
                        },
                        |samples| samples[branch].density_kg_m3(),
                    )
                })
            },
            || checkpoint(cx).is_err(),
        )?;
        checkpoint(cx)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
