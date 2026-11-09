//! Distributed pipe inventory on the unchanged hydraulic graph.
//!
//! Junctions keep their explicit finite mixing cells. Each branch adds an
//! ordered sequence of original thermal control volumes; one signed reference
//! mass flow passes through every segment. Native upwinding reverses donors
//! automatically, including fluid already stored inside a pipe. Only transport
//! links expand: hydraulic pressure/momentum dimensions never change.
//!
//! This is finite-volume axial transport with perfectly mixed cells. Spatial
//! and temporal refinement remain necessary for a resolved residence-time
//! distribution. Empty paths explicitly retain direct junction-to-junction
//! transport; they supply neither pipe storage nor inferred inertia.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_qty::Density;
use fs_thermal::enthalpy::advection::MassBoundary;
use fs_thermal::enthalpy::{CellGeometry, CellId, DepositionSimulation, MaterialId};

use super::{GraphEnthalpyTransport, TransientTransportError, poll};
use crate::graph::LossGraph;
use crate::graph::gravity::GraphGravity;
use crate::graph::inertial::{InertialBranch, InertialGraph};
use fs_thermal::enthalpy::advection::reconstruction::{
    AxialCellStencil, LimitedEnthalpyReconstruction,
};
use std::sync::Arc;

/// Spatial inventory may be refined independently of the bounded hydraulic
/// descriptor. The native thermal owner's own work/storage limits still apply.
pub const MAX_DISTRIBUTED_FLUID_CELLS: usize = 4096;

#[derive(Clone, Debug)]
struct BoundInventory {
    material: MaterialId,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric_identity: Option<[u8; 32]>,
}

/// Fixed fluid topology and actual finite inventory; no time or heat state is
/// stored here. The original thermal owner retains all caloric laws, energies,
/// boundaries and histories. Source-provided caloric identities, when present,
/// are pinned alongside material IDs, complete geometry and conserved masses.
#[derive(Clone, Debug)]
pub struct DistributedGraphInventory {
    node_cells: Vec<CellId>,
    branch_cells: Vec<Vec<CellId>>,
    topology: Vec<(usize, usize, String)>,
    reference_density: Density,
    all_cells: Vec<CellId>,
    /// Expanded thermal indexes and their original hydraulic branch.
    connections: Vec<(usize, usize, usize)>,
    storage: BTreeMap<CellId, BoundInventory>,
    node_positions_m: Vec<Option<[f64; 3]>>,
    branch_volumes_m3: Vec<f64>,
    inertial_branches: Vec<Option<InertialBranch>>,
    inertances_pa_s2_m3: Vec<Option<f64>>,
    maximum_connections: usize,
    reconstruction: Option<Arc<LimitedEnthalpyReconstruction>>,
}

impl DistributedGraphInventory {
    /// Bind branch cells in positive declared-edge order, never CellId order.
    /// Every node and pipe cell is distinct, active, native enthalpy storage
    /// with M=rho0 V. Nonempty pipe paths require positive straight control-
    /// volume lengths and contiguous endpoints; curved pipes are polylines of
    /// those cells. Shared junction endpoints agree with each other and any
    /// original gravity chart. Cross-section is the authored V/length.
    ///
    /// An explicitly empty path remains a direct connection. Its missing
    /// volume/length is never invented or used to derive hydraulic inertance.
    /// Native component caloric/datum admission is reused even at initial rest.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cx: &Cx<'_>,
        graph: &LossGraph,
        reference_density: Density,
        node_cells: &[CellId],
        branch_cells: &[Vec<CellId>],
        thermal: &DepositionSimulation<'_>,
        maximum_connections: usize,
    ) -> Result<Self, TransientTransportError> {
        poll(cx)?;
        if node_cells.len() != graph.node_count()
            || node_cells.is_empty()
            || branch_cells.len() != graph.branches().len()
            || !(1..=1_000_000).contains(&maximum_connections)
            || !positive(reference_density.value())
        {
            return Err(TransientTransportError::Invalid(
                "complete bounded distributed graph inventory required",
            ));
        }
        let pipe_cells = branch_cells
            .iter()
            .try_fold(0usize, |n, path| n.checked_add(path.len()))
            .ok_or(TransientTransportError::Invalid(
                "pipe inventory count overflow",
            ))?;
        let count =
            node_cells
                .len()
                .checked_add(pipe_cells)
                .ok_or(TransientTransportError::Invalid(
                    "fluid inventory count overflow",
                ))?;
        let connections_count = graph.branches().len().checked_add(pipe_cells).ok_or(
            TransientTransportError::Invalid("transport connection count overflow"),
        )?;
        if count > MAX_DISTRIBUTED_FLUID_CELLS || connections_count > maximum_connections {
            return Err(TransientTransportError::Invalid(
                "distributed fluid cell or connection budget exceeded",
            ));
        }
        let mut all_cells = Vec::with_capacity(count);
        all_cells.extend_from_slice(node_cells);
        for path in branch_cells {
            all_cells.extend_from_slice(path);
        }
        let mut storage = BTreeMap::new();
        for &id in &all_cells {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(TransientTransportError::Invalid(
                    "distributed fluid cell is not active",
                ))?;
            let bound = BoundInventory {
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
            };
            if storage.insert(id, bound).is_some() {
                return Err(TransientTransportError::Invalid(
                    "node and pipe inventory cells must be distinct",
                ));
            }
        }
        let mut node_positions_m = graph.gravity().map_or_else(
            || vec![None; node_cells.len()],
            |g| g.positions_m().iter().copied().map(Some).collect(),
        );
        let mut connections = Vec::with_capacity(connections_count);
        let mut inertial_branches = Vec::with_capacity(branch_cells.len());
        let mut inertances = Vec::with_capacity(branch_cells.len());
        let mut volumes = Vec::with_capacity(branch_cells.len());
        let mut offset = node_cells.len();
        for (branch, (edge, path)) in graph.branches().iter().zip(branch_cells).enumerate() {
            poll(cx)?;
            let mut previous_index = edge.from;
            let mut previous_end = None;
            let mut total_length = 0.0;
            let mut total_volume = 0.0;
            let mut length_over_area = 0.0;
            for (ordinal, &id) in path.iter().enumerate() {
                poll(cx)?;
                let geometry = storage[&id].geometry;
                let delta: [f64; 3] =
                    core::array::from_fn(|axis| geometry.end_m[axis] - geometry.start_m[axis]);
                let length = fs_math::det::sqrt(delta.iter().map(|x| x * x).sum());
                let area = geometry.volume_m3 / length;
                if !positive(length)
                    || !positive(area)
                    || previous_end.is_some_and(|end| end != geometry.start_m)
                {
                    return Err(TransientTransportError::Invalid(
                        "pipe cells require contiguous positive length and volume geometry",
                    ));
                }
                if ordinal == 0 {
                    bind_position(&mut node_positions_m[edge.from], geometry.start_m)?;
                }
                previous_end = Some(geometry.end_m);
                total_length += length;
                total_volume += geometry.volume_m3;
                length_over_area += length / area;
                if !positive(total_length) || !positive(total_volume) || !positive(length_over_area)
                {
                    return Err(TransientTransportError::Invalid(
                        "distributed pipe geometry is unrepresentable",
                    ));
                }
                let index = offset + ordinal;
                connections.push((previous_index, index, branch));
                previous_index = index;
            }
            connections.push((previous_index, edge.to, branch));
            offset += path.len();
            if let Some(end) = previous_end {
                bind_position(&mut node_positions_m[edge.to], end)?;
                let equivalent = InertialBranch {
                    length_m: total_length,
                    flow_area_m2: total_length / length_over_area,
                };
                let inertance = reference_density.value() * length_over_area;
                if !positive(equivalent.flow_area_m2)
                    || !positive(inertance)
                    || !positive(1.0 / inertance)
                {
                    return Err(TransientTransportError::Invalid(
                        "distributed pipe inertance is unrepresentable",
                    ));
                }
                inertial_branches.push(Some(equivalent));
                inertances.push(Some(inertance));
            } else {
                inertial_branches.push(None);
                inertances.push(None);
            }
            volumes.push(total_volume);
        }
        let layout = Self {
            node_cells: node_cells.to_vec(),
            branch_cells: branch_cells.to_vec(),
            topology: graph
                .branches()
                .iter()
                .map(|b| (b.from, b.to, b.loss.name.clone()))
                .collect(),
            reference_density,
            all_cells,
            connections,
            storage,
            node_positions_m,
            branch_volumes_m3: volumes,
            inertial_branches,
            inertances_pa_s2_m3: inertances,
            maximum_connections,
            reconstruction: None,
        };
        // The same native reference-mass, storage-basis, component-law and
        // continuity admission as ordinary graph transport; no second law.
        layout.bind(
            cx,
            graph,
            &vec![0.0; branch_cells.len()],
            thermal,
            vec![],
            maximum_connections,
        )?;
        poll(cx)?;
        Ok(layout)
    }

    /// Original hydraulic node order, excluding interior pipe cells.
    pub fn node_cells(&self) -> &[CellId] {
        &self.node_cells
    }
    /// Original hydraulic branch order and positive-direction pipe traversal.
    pub fn branch_cells(&self) -> &[Vec<CellId>] {
        &self.branch_cells
    }

    /// Prepare conservative minmod face enthalpies on original interior cells.
    ///
    /// A stencil contains three consecutive cells of ONE authored pipe path.
    /// Physical nonuniform center distances and half lengths are derived by
    /// the native thermal reconstruction owner. Junction storage, first/last
    /// pipe cells, short paths and direct connections keep donor-cell upwind
    /// transport. No junction gradient or hidden control volume is inferred.
    ///
    /// The shared plan is used by native inertial and quasi-steady endpoint
    /// residuals, including reversal, and by the original explicit timeline.
    /// Linear frozen-flow shortcuts require upwind transport; reconstructed
    /// implicit intervals use the nonlinear endpoint entry.
    pub fn with_limited_reconstruction(
        mut self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        thermal: &DepositionSimulation<'_>,
        maximum_stencils: usize,
    ) -> Result<Self, TransientTransportError> {
        self.admit_storage(cx, graph, thermal)?;
        let count: usize = self
            .branch_cells
            .iter()
            .map(|path| path.len().saturating_sub(2))
            .sum();
        if maximum_stencils > MAX_DISTRIBUTED_FLUID_CELLS || count > maximum_stencils {
            return Err(TransientTransportError::Invalid(
                "bounded original pipe reconstruction stencil count required",
            ));
        }
        let mut stencils = Vec::with_capacity(count);
        for path in &self.branch_cells {
            for cells in path.windows(3) {
                poll(cx)?;
                stencils.push(AxialCellStencil {
                    previous: cells[0],
                    cell: cells[1],
                    next: cells[2],
                });
            }
        }
        let plan =
            LimitedEnthalpyReconstruction::try_new(thermal, &stencils, maximum_stencils, || {
                super::thermal_checkpoint(cx)
            })?;
        self.reconstruction = (!plan.is_empty()).then(|| Arc::new(plan));
        Ok(self)
    }

    pub fn reconstruction(&self) -> Option<&LimitedEnthalpyReconstruction> {
        self.reconstruction.as_deref()
    }
    pub const fn reference_density(&self) -> Density {
        self.reference_density
    }
    pub fn interior_connection_count(&self) -> usize {
        self.connections.len()
    }
    pub fn fluid_cell_count(&self) -> usize {
        self.all_cells.len()
    }
    /// Inferred only from actual path endpoints or the original gravity chart.
    pub fn node_positions_m(&self) -> &[Option<[f64; 3]>] {
        &self.node_positions_m
    }
    pub fn branch_volumes_m3(&self) -> &[f64] {
        &self.branch_volumes_m3
    }
    /// Equivalent geometry preserves rho0 sum(length_i/area_i), including
    /// varying cross-section. Empty direct paths have no derived geometry.
    pub fn inertial_branches(&self) -> &[Option<InertialBranch>] {
        &self.inertial_branches
    }

    /// Admit the same conserved density, junction frame and original pipe path.
    /// Distributed gravity follows exactly one section per stored pipe cell.
    /// Empty direct paths own no interior inventory, so only their endpoints
    /// are checked; no interior geometry or storage is inferred for them.
    pub fn admit_gravity(
        &self,
        cx: &Cx<'_>,
        gravity: &GraphGravity,
    ) -> Result<(), TransientTransportError> {
        poll(cx)?;
        if gravity.reference_transport().density_kg_m3() != self.reference_density.value()
            || gravity.positions_m().len() != self.node_cells.len()
            || gravity.departure_heads_pa().len() != self.branch_cells.len()
            || gravity
                .positions_m()
                .iter()
                .zip(&self.node_positions_m)
                .any(|(actual, expected)| expected.is_some_and(|expected| *actual != expected))
        {
            return Err(TransientTransportError::Invalid(
                "hydraulic gravity and distributed reference inventory disagree",
            ));
        }
        if let Some(branches) = gravity.branch_sections() {
            for (path, sections) in self.branch_cells.iter().zip(branches) {
                poll(cx)?;
                if path.is_empty() {
                    continue;
                }
                if path.len() != sections.len() {
                    return Err(TransientTransportError::Invalid(
                        "gravity must retain one section per stored pipe cell",
                    ));
                }
                for (id, section) in path.iter().zip(sections) {
                    poll(cx)?;
                    let geometry = &self.storage[id].geometry;
                    if section.start_m != geometry.start_m || section.end_m != geometry.end_m {
                        return Err(TransientTransportError::Invalid(
                            "hydraulic gravity section differs from the stored pipe geometry",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Check that each inventoried branch's retained momentum uses the SAME
    /// fluid geometry and rho0. Empty paths keep caller-declared inertia.
    /// This does not attribute junction mixing volumes to branch kinetic energy.
    pub fn admit_momentum(
        &self,
        cx: &Cx<'_>,
        momentum: &InertialGraph,
    ) -> Result<(), TransientTransportError> {
        poll(cx)?;
        if momentum.node_count() != self.node_cells.len()
            || momentum.branch_count() != self.branch_cells.len()
            || momentum.reference_density() != self.reference_density
        {
            return Err(TransientTransportError::Invalid(
                "momentum and distributed inventory disagree",
            ));
        }
        for (expected, &actual) in self
            .inertances_pa_s2_m3
            .iter()
            .zip(momentum.inertances_pa_s2_m3())
        {
            poll(cx)?;
            if expected.is_some_and(|expected| !near(expected, actual)) {
                return Err(TransientTransportError::Invalid(
                    "retained branch inertia differs from the stored pipe geometry",
                ));
            }
        }
        Ok(())
    }

    /// Nominal pipe-only V/|Q| residence times in seconds under a held rate.
    /// None means zero flow; an empty moving path has Some(0). These are not
    /// arrival-time predictions for an evolving/reversing transient.
    pub fn nominal_residence_times_s(
        &self,
        cx: &Cx<'_>,
        signed_volume_flows_m3_s: &[f64],
    ) -> Result<Vec<Option<f64>>, TransientTransportError> {
        self.admit_rates(signed_volume_flows_m3_s)?;
        self.branch_volumes_m3
            .iter()
            .zip(signed_volume_flows_m3_s)
            .map(|(&volume, &flow)| {
                poll(cx)?;
                if flow == 0.0 {
                    return Ok(None);
                }
                let time = volume / flow.abs();
                if !time.is_finite() || (volume != 0.0 && time == 0.0) {
                    return Err(TransientTransportError::Invalid(
                        "pipe residence time is unrepresentable",
                    ));
                }
                Ok(Some(time))
            })
            .collect()
    }

    /// Admit the unchanged physical inventory and graph identities at a thermal
    /// candidate without compiling transport or inventing a time state. This is
    /// also the source-adapter seam for exact native caloric and geometry checks.
    pub fn admit_storage(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), TransientTransportError> {
        poll(cx)?;
        if graph.node_count() != self.node_cells.len()
            || graph.branches().len() != self.topology.len()
            || graph
                .branches()
                .iter()
                .zip(&self.topology)
                .any(|(a, (from, to, name))| a.from != *from || a.to != *to || a.loss.name != *name)
        {
            return Err(TransientTransportError::Invalid(
                "original distributed hydraulic topology changed",
            ));
        }
        if let Some(gravity) = graph.gravity() {
            self.admit_gravity(cx, gravity)?;
        }
        for (&id, bound) in &self.storage {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(TransientTransportError::Invalid(
                    "distributed fluid cell disappeared",
                ))?;
            if cell.material_id != bound.material
                || cell.mass_kg != bound.mass_kg
                || cell.geometry != bound.geometry
                || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
            {
                return Err(TransientTransportError::Invalid(
                    "original distributed storage geometry, mass or caloric identity changed",
                ));
            }
        }
        Ok(())
    }

    /// Compile original signed branch flow into consecutive native mass links.
    /// External stream cells must be original junction cells; an interior
    /// injection changes hydraulic topology and must be authored as a junction.
    ///
    /// Connection IDs enumerate declared branch order, then each branch's
    /// positive-direction segments. External IDs must be disjoint from that
    /// entire expanded range. Both construction and call budgets remain active.
    #[allow(clippy::too_many_arguments)]
    pub fn bind(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        signed_volume_flows_m3_s: &[f64],
        thermal: &DepositionSimulation<'_>,
        boundaries: Vec<MassBoundary>,
        maximum_connections: usize,
    ) -> Result<GraphEnthalpyTransport, TransientTransportError> {
        poll(cx)?;
        self.admit_rates(signed_volume_flows_m3_s)?;
        self.admit_storage(cx, graph, thermal)?;
        for boundary in &boundaries {
            poll(cx)?;
            let id = match boundary {
                MassBoundary::Inlet { cell, .. } | MassBoundary::Outlet { cell, .. } => cell,
            };
            if !self.node_cells.contains(id) {
                return Err(TransientTransportError::Invalid(
                    "external stream must bind an original hydraulic junction",
                ));
            }
        }
        let mut transport = GraphEnthalpyTransport::bind(
            cx,
            self.all_cells.len(),
            graph
                .gravity()
                .map(|g| g.reference_transport().density_kg_m3()),
            self.connections
                .iter()
                .map(|&(a, b, branch)| (a, b, signed_volume_flows_m3_s[branch])),
            self.reference_density,
            &self.all_cells,
            thermal,
            boundaries,
            maximum_connections.min(self.maximum_connections),
        )?;
        // Only the thermal index map expands; this remains the original
        // hydraulic node map exposed by every public graph transport binding.
        transport.node_cells = self.node_cells.clone();
        transport.reconstruction = self.reconstruction.clone();
        if transport.reconstruction.is_some() {
            // Admit the actual signed graph against the prepared original
            // stencil topology, including zero-rate connections at startup.
            transport.endpoint_transport(cx, thermal)?;
        }
        Ok(transport)
    }

    /// Integrate both signs of the original affine branch histories on every
    /// authored pipe segment. The original finite cells and limited face
    /// reconstruction remain the sole enthalpy owners. A reversed segment
    /// adds one native transfer, with ID `interior_connection_count + index`;
    /// both construction and call connection budgets include these transfers.
    /// External inlet and outlet rates each retain their actual swept mass.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_affine(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        initial_m3_s: &[f64],
        mean_m3_s: &[f64],
        final_m3_s: &[f64],
        thermal: &DepositionSimulation<'_>,
        boundaries: Vec<MassBoundary>,
        maximum_connections: usize,
    ) -> Result<GraphEnthalpyTransport, TransientTransportError> {
        let rates = super::affine::admit_rates(cx, graph, initial_m3_s, mean_m3_s, final_m3_s)?;
        let streams = super::affine::node_mass_streams(cx, graph, self.reference_density,
            initial_m3_s, mean_m3_s, final_m3_s)?;
        super::affine::admit_boundaries(cx, &self.node_cells, &streams, &boundaries)?;
        let base = self.bind(cx, graph, mean_m3_s, thermal, boundaries, maximum_connections)?;
        let expanded: Vec<_> = self.connections.iter().map(|&(_, _, branch)| rates[branch]).collect();
        base.with_affine_rates(cx, thermal, self.reference_density, &expanded)
    }

    fn admit_rates(&self, rates: &[f64]) -> Result<(), TransientTransportError> {
        if rates.len() != self.branch_cells.len() || rates.iter().any(|v| !v.is_finite()) {
            Err(TransientTransportError::Invalid(
                "one finite signed rate per original branch required",
            ))
        } else {
            Ok(())
        }
    }
}

fn bind_position(
    slot: &mut Option<[f64; 3]>,
    position: [f64; 3],
) -> Result<(), TransientTransportError> {
    if position.iter().any(|v| !v.is_finite()) || slot.is_some_and(|old| old != position) {
        return Err(TransientTransportError::Invalid(
            "shared hydraulic junction and pipe endpoints disagree",
        ));
    }
    *slot = Some(position);
    Ok(())
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn near(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

#[cfg(test)]
mod tests;
