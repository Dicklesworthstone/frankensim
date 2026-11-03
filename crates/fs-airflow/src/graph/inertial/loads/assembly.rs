//! Bind point-junction and assembly wrenches to one native hydraulic trial.

use super::super::{InertialGraph, InertialGraphReport};
use super::junction::{self, JunctionLoad, JunctionState, PortInterval, Sum, sum_loads};
use super::{InertialGraphError, InertialGraphTrial, PipeLoad, PipeLoadGeometry, checkpoint};
use fs_exec::Cx;
use fs_qty::Pressure;

/// Per-component reactions and the directly evaluated whole-assembly load.
/// Pipes follow original branch order; junctions follow original node order.
/// Apply either the component loads to their distinct supports or the total
/// to the assembly. Applying both would count the same fluid impulse twice.
#[derive(Clone, Debug)]
pub struct HydraulicAssemblyLoads {
    pub pipes: Vec<PipeLoad>,
    pub junctions: Vec<JunctionLoad>,
    /// Internal pressure and advective endpoint exchanges are excluded before
    /// summation. This avoids subtracting large opposing interface loads to
    /// recover a small net support load; all fluid storage and weight remain.
    pub total: PipeLoad,
    /// Signed sum of the retained per-junction kinetic-energy transports.
    /// Individual terms remain available; their cancellation proves no local
    /// entropy or thermal closure. No thermal energy is deposited by this owner.
    pub junction_kinetic_energy_transport_j: f64,
}

impl InertialGraphTrial {
    /// Reaction on one original zero-storage point junction, including every
    /// incident branch once. Supply geometry in original branch order.
    ///
    /// A branch starting at the node has outward normal equal to its first
    /// tangent and outward flow Q. An ending branch has minus its last tangent
    /// and outward flow -Q. All incident endpoint positions must coincide.
    /// Both endpoint flow balances must close; a fixed-pressure reservoir with
    /// an unrepresented stream therefore cannot be mistaken for a junction.
    ///
    /// The original common STATIC nodal pressure is retained. Unequal areas
    /// change momentum flux and pressure traction; they do not establish a
    /// reducer pressure-recovery or passive mixing/heat law. The returned signed
    /// kinetic-energy transport exposes that unresolved term explicitly.
    pub fn junction_load(
        &self,
        cx: &Cx<'_>,
        node: usize,
        geometries: &[PipeLoadGeometry],
        exterior_pressure: Pressure,
        about_m: [f64; 3],
    ) -> Result<JunctionLoad, InertialGraphError> {
        checkpoint(cx)?;
        let c = &self.initial.compiled;
        if geometries.len() != c.branches.len() || node >= c.free.len() {
            return Err(InertialGraphError::Invalid(
                "original node and one geometry per branch required",
            ));
        }
        let mut position = None;
        let mut ports = Vec::with_capacity(c.topology.adjacency[node].len());
        for &branch in &c.topology.adjacency[node] {
            checkpoint(cx)?;
            let geometry = &geometries[branch];
            self.admit_pipe(cx, branch, geometry)?;
            let at_start = c.topology.branches()[branch].from == node;
            let (point, normal) = geometry.junction_port(at_start);
            if position.is_some_and(|original| original != point) {
                return Err(InertialGraphError::Invalid(
                    "point-junction pipe endpoints must coincide",
                ));
            }
            position = Some(point);
            let sign = if at_start { 1.0 } else { -1.0 };
            ports.push(PortInterval {
                outward_normal: normal,
                area_m2: geometry.area_m2,
                outflows_m3_s: [
                    sign * self.report.initial_flows_m3_s[branch],
                    sign * self.report.final_flows_m3_s[branch],
                ],
            });
        }
        let position_m = position.ok_or(InertialGraphError::Invalid(
            "junction requires an incident pipe",
        ))?;
        let result = junction::evaluate(
            JunctionState {
                node,
                position_m,
                about_m,
                duration_s: self.report.duration_s,
                density_kg_m3: c.reference_density.value(),
                pressure_pa: self.report.average_pressures_pa[node],
                exterior_pressure_pa: exterior_pressure.value(),
            },
            &ports,
            || checkpoint(cx).is_err(),
        )?;
        checkpoint(cx)?;
        Ok(result)
    }

    /// Resolve a complete fixed hydraulic assembly from the SAME trajectory.
    ///
    /// All free nodes are internal junctions. Fixed-pressure nodes remain
    /// external stream cuts unless listed in enclosed_pressure_nodes. A closed
    /// loop's pressure-datum node must be listed explicitly; no reservoir/datum
    /// role is guessed from zero net flow in just one interval. Listed nodes
    /// declare no external stream; endpoint continuity is checked independently.
    ///
    /// Pipe weight and retained-fluid acceleration are included exactly once.
    /// Internal junction pressure/momentum reactions oppose their pipe faces;
    /// the total cancels these interfaces analytically before accumulation.
    /// There is no extra species, enthalpy, gravity-work or heat transaction.
    pub fn assembly_loads(
        &self,
        cx: &Cx<'_>,
        geometries: &[PipeLoadGeometry],
        enclosed_pressure_nodes: &[usize],
        exterior_pressure: Pressure,
        about_m: [f64; 3],
    ) -> Result<HydraulicAssemblyLoads, InertialGraphError> {
        checkpoint(cx)?;
        let c = &self.initial.compiled;
        if geometries.len() != c.branches.len() || enclosed_pressure_nodes.len() > c.free.len() {
            return Err(InertialGraphError::Invalid(
                "complete bounded native assembly geometry required",
            ));
        }
        let mut enclosed: Vec<_> = c.free.iter().map(Option::is_some).collect();
        for &node in enclosed_pressure_nodes {
            checkpoint(cx)?;
            if node >= enclosed.len() || enclosed[node] {
                return Err(InertialGraphError::Invalid(
                    "distinct original pressure-datum nodes required",
                ));
            }
            enclosed[node] = true;
        }
        let mut junctions = Vec::new();
        let mut energy = Sum::default();
        for (node, &inside) in enclosed.iter().enumerate() {
            if inside {
                let load = self.junction_load(cx, node, geometries, exterior_pressure, about_m)?;
                energy.add(load.kinetic_energy_transport_j)?;
                junctions.push(load);
            }
        }
        let mut pipes = Vec::with_capacity(geometries.len());
        let mut exterior_parts = Vec::with_capacity(geometries.len());
        for (branch, geometry) in geometries.iter().enumerate() {
            checkpoint(cx)?;
            let load = self.pipe_load(cx, branch, geometry, exterior_pressure, about_m)?;
            let edge = &c.topology.branches()[branch];
            let exposed = [!enclosed[edge.from], !enclosed[edge.to]];
            exterior_parts.push(if exposed == [true; 2] {
                load
            } else {
                self.pipe_load_with_ports(
                    cx,
                    branch,
                    geometry,
                    exterior_pressure,
                    about_m,
                    exposed,
                )?
            });
            pipes.push(load);
        }
        let total = sum_loads(self.report.duration_s, about_m, &exterior_parts)?;
        checkpoint(cx)?;
        Ok(HydraulicAssemblyLoads {
            pipes,
            junctions,
            total,
            junction_kinetic_energy_transport_j: energy.value()?,
        })
    }

    /// Accept hydraulic history and its mechanically usable assembly loads
    /// together. Geometry, stream closure, cancellation or native impulse/work
    /// failure returns no accepted momentum and no accepted load package.
    #[allow(clippy::too_many_arguments)]
    pub fn accept_with_assembly_loads(
        self,
        cx: &Cx<'_>,
        geometries: &[PipeLoadGeometry],
        enclosed_pressure_nodes: &[usize],
        exterior_pressure: Pressure,
        about_m: [f64; 3],
        impulse_tolerance_pa_s: f64,
        mechanical_tolerance_j: f64,
    ) -> Result<(InertialGraph, InertialGraphReport, HydraulicAssemblyLoads), InertialGraphError>
    {
        let loads = self.assembly_loads(
            cx,
            geometries,
            enclosed_pressure_nodes,
            exterior_pressure,
            about_m,
        )?;
        let (next, report) = self.accept(cx, impulse_tolerance_pa_s, mechanical_tolerance_j)?;
        Ok((next, report, loads))
    }
}
