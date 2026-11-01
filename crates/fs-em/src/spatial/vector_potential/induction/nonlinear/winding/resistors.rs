//! Source-bound ordinary resistors on the retained three-dimensional winding graph.
//!
//! Reassembly changes dissipative coefficients only. Original A, winding
//! currents, native flux/charge, source profiles and terminal ordering remain
//! unchanged; resistance changes add no magnetic or circuit parameter work.

use super::{
    GraphResistorMaterial, Layout, NonlinearWindingPort, Result, SpatialNonlinearWindingCircuit,
    invalid, poll, positive,
};
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_phs::circuit::{Branch, assemble_circuit};
use std::rc::Rc;

impl<C: SpatialConductorMaterial> SpatialNonlinearWindingCircuit<'_, C> {
    /// Sources at the last accepted resistance update, in original branch order.
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.graph_materials
    }

    /// Replace ALL ordinary graph-resistor values with source-resolved material
    /// data, returning an independent owner at the SAME physical time/history.
    /// This explicit source reconfiguration can attach initial material bindings
    /// to an authored graph. A finite thermal owner must separately pin original
    /// sources and admit their owned-temperature refresh across physical steps.
    /// Every subsequent candidate retains these sources beside graph-local heat.
    ///
    /// Previous algebraic entries remain the last native receipt, not a new
    /// equilibrium claim. The next native descriptor step solves actual stage
    /// efforts. No artificial initial solve may project accepted flux or charge.
    /// Cached geometry/profile/Whitney loads are shared without reconstruction.
    pub fn with_graph_resistors(
        &self,
        mut sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = self
            .graph()
            .branches
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Resistor { .. }))
            .count();
        if sources.len() != count || sources.len() > 64 {
            return Err(invalid(
                "every ordinary graph resistor requires one bounded material source",
            ));
        }
        sources.sort_by_key(|source| source.branch);
        let mut graph = self.graph().clone();
        for (index, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            if index > 0 && sources[index - 1].branch == source.branch
                || !source
                    .conductor
                    .source_property()
                    .is_static_at_fixed_external_conditions()
                || !positive(source.conductor.resistance_ohm())
                || !positive(source.conductor.resistance_ohm().recip())
            {
                return Err(invalid(
                    "unique representable static Ohmic graph sources required",
                ));
            }
            let Some((_, _, branch @ Branch::Resistor { .. })) =
                graph.branches.get_mut(source.branch)
            else {
                return Err(invalid(
                    "material source does not name an original graph resistor",
                ));
            };
            *branch = Branch::Resistor {
                ohms: source.conductor.resistance_ohm(),
            };
        }
        if sources == self.graph_materials {
            poll(cancelled)?;
            return Ok(self.clone());
        }
        // Preserve the exact original winding order and private port numbers.
        let mut lowered = graph.clone();
        for (winding, &branch) in self.branch_indices().iter().enumerate() {
            let port = self.layout.external_ports + winding;
            lowered.branches[branch].2 = match self.terminal_ports()[winding] {
                NonlinearWindingPort::CurrentInput => Branch::CurrentSource { port },
                NonlinearWindingPort::VoltageInput { .. } => Branch::VoltageSource { port },
            };
        }
        let dae = assemble_circuit(&lowered)?;
        if dae.flux_index != self.dae().flux_index
            || dae.charge_index != self.dae().charge_index
            || dae.node_potential_index != self.dae().node_potential_index
            || dae.ports != self.dae().ports
            || dae.system.state_dim() != self.dae().system.state_dim()
            || dae.system.hamiltonian(&self.state).to_bits()
                != self.dae().system.hamiltonian(&self.state).to_bits()
        {
            return Err(invalid(
                "resistance-only reassembly changed original circuit storage or ports",
            ));
        }
        let mut candidate = self.clone();
        candidate.layout = Rc::new(Layout {
            graph,
            dae,
            windings: Rc::clone(&self.layout.windings),
            ports: self.layout.ports.clone(),
            external_ports: self.layout.external_ports,
            loads_per_a: Rc::clone(&self.layout.loads_per_a),
        });
        candidate.graph_materials = sources;
        poll(cancelled)?;
        Ok(candidate)
    }
}
