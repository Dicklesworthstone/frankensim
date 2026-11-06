//! Source-bound ordinary resistors on the retained three-dimensional winding graph.
//!
//! Reassembly changes dissipative coefficients only. Original A, winding
//! currents, native flux/charge, source profiles and terminal ordering remain
//! unchanged; resistance changes add no magnetic or circuit parameter work.

use super::{GraphResistorMaterial, Result, SpatialNonlinearWindingCircuit};
use fs_material::conductor::spatial::SpatialConductorMaterial;

impl<C: SpatialConductorMaterial> SpatialNonlinearWindingCircuit<'_, C> {
    /// Sources at the last accepted resistance update, in original branch order.
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.terminal.graph_materials
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
        sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Ok(Self {
            terminal: self.terminal.with_graph_resistors(sources, cancelled)?,
            field: self.field.clone(),
        })
    }
}
