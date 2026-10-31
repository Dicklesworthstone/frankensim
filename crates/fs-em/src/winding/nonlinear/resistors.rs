//! Ordinary source-bound resistors on the original nonlinear winding graph.
//!
//! Reassembly changes only dissipative coefficients. The native field, source
//! profiles, private ports, conserved circuit coordinates and magnetic linkage
//! remain the same; resistance changes add no stored-energy parameter work.

use super::{Layout, NonlinearWindingCircuit, NonlinearWindingError, Result, poll, positive};
use fs_material::conductor::ResolvedConductor;
use fs_phs::circuit::{Branch, assemble_circuit};
use std::rc::Rc;

/// Complete original material/geometry source for one ordinary graph resistor.
#[derive(Clone, Debug, PartialEq)]
pub struct GraphResistorMaterial {
    /// Original graph branch index, not its position among resistor branches.
    pub branch: usize,
    /// Original conductivity/resistivity basis, source bundle and path geometry.
    pub conductor: ResolvedConductor,
}

impl NonlinearWindingCircuit<'_> {
    /// Sources at the last accepted resistance update, in original branch order.
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.graph_materials
    }

    /// Replace ALL ordinary graph-resistor values with resolved material data.
    /// Returns an independent owner at the SAME physical time and stored state.
    /// This explicit source reconfiguration may also attach the first material
    /// bindings to an authored graph. Thermal callers must separately admit the
    /// accepted source state and pin source choices across their physical step.
    /// Retained algebraic entries are the previous native receipt, not new
    /// equilibrium efforts; the subsequent original native descriptor solve
    /// rebuilds the stage efforts under the actual held terminal inputs.
    ///
    /// This supplies no material interpolant, contact law, new circuit solver or
    /// parameter work. The original source must describe a static isotropic
    /// Ohmic scalar at held external conditions. Every later native step retains
    /// these sources beside its existing graph-local heat records.
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
            return Err(NonlinearWindingError::Invalid(
                "every ordinary graph resistor needs one bounded material source",
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
                return Err(NonlinearWindingError::Invalid(
                    "unique representable static Ohmic graph sources required",
                ));
            }
            let Some((_, _, branch @ Branch::Resistor { .. })) =
                graph.branches.get_mut(source.branch)
            else {
                return Err(NonlinearWindingError::Invalid(
                    "material source does not name an original graph resistor",
                ));
            };
            *branch = Branch::Resistor {
                ohms: source.conductor.resistance_ohm(),
            };
        }
        // Retain the exact original winding ordering and private-port numbers.
        let mut lowered = graph.clone();
        for (winding, &branch) in self.branch_indices().iter().enumerate() {
            lowered.branches[branch].2 =
                self.terminal_ports()[winding].branch(self.layout.external_ports + winding);
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
            return Err(NonlinearWindingError::Invalid(
                "resistance-only reassembly changed original circuit storage or ports",
            ));
        }
        let mut candidate = self.clone();
        candidate.layout = Rc::new(Layout {
            graph,
            dae,
            windings: self.layout.windings.clone(),
            terminal_ports: self.layout.terminal_ports.clone(),
            external_ports: self.layout.external_ports,
        });
        candidate.graph_materials = sources;
        poll(cancelled)?;
        Ok(candidate)
    }
}
