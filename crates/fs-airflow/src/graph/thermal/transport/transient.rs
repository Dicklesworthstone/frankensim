//! Finite fluid inventory on the original conservative thermal timeline.
//!
//! Hydraulic nodes are explicitly bound to caller-owned thermal control
//! volumes. Their original mass, geometry, caloric law and enthalpy remain the
//! storage authority. Signed reference-volume flow becomes signed mass flow;
//! native upwind transfers move that enthalpy without a second heat integrator.
//! A pressure datum supplies no fluid or heat. Every external mass connection
//! is explicit and every storage cell must balance mass, including datums.
//!
//! Node storage is perfectly mixed; resolving pipe residence time requires
//! caller-declared spatial control volumes. The base binding prescribes flow;
//! [`coupled`] resolves quasi-steady hydraulic and wall feedback at each native
//! implicit endpoint. [`inertial`] couples retained branch momentum to the same
//! thermal residual using native step-average rates. Reference density stays
//! fixed; neither binding adds compressible storage.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_qty::Density;
use fs_thermal::enthalpy::advection::{Advected, AdvectionAudit, MassBoundary, MassFlow};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, TransportStepCoefficients};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, EnergyAudit, MaterialId,
    StepConfig, ThermalStorageBasis,
};

use super::recirculation::{MAX_RECIRCULATION_SUPPLIES, solve_feedback};
use super::{TransportError, poll};
use crate::graph::{GraphSolution, LossGraph};

/// Existing thermal, linear-solve or explicit binding refusal. No partial
/// timeline is returned on failure or cancellation.
#[derive(Debug)]
pub enum TransientTransportError {
    Invalid(&'static str),
    Thermal(DepositionError),
    Transport(TransportError),
}
impl From<DepositionError> for TransientTransportError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<TransportError> for TransientTransportError {
    fn from(error: TransportError) -> Self {
        Self::Transport(error)
    }
}
impl std::fmt::Display for TransientTransportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "transient airflow transport: {self:?}")
    }
}
impl std::error::Error for TransientTransportError {}

/// Original finite-inventory transport with no implicit pressure reservoirs.
/// Interior connection IDs follow the binding's transport order: original
/// branch indices for direct bindings, expanded branch segments for distributed
/// pipe inventory. Caller-supplied boundary IDs must be distinct from that range.
#[derive(Clone)]
pub struct GraphEnthalpyTransport {
    node_cells: Vec<CellId>,
    storage: BTreeMap<CellId, BoundTransportStorage>,
    advection: Advected<()>,
    max_connections: usize,
}

#[derive(Clone)]
struct BoundTransportStorage {
    material: MaterialId,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric_identity: Option<[u8; 32]>,
}

impl GraphEnthalpyTransport {
    /// Bind explicit storage cells to a prescribed signed hydraulic solution.
    /// Cell mass must match reference density times its declared volume within
    /// roundoff. No volume, heat capacity, inlet enthalpy or caloric datum is
    /// inferred. Native advection checks common caloric identity and mass
    /// closure at ALL cells before this returns, regardless of pressure BCs.
    /// Other thermal cells/links may coexist for the general explicit path.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cx: &Cx<'_>,
        flow: &GraphSolution,
        reference_density: Density,
        node_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
        boundaries: Vec<MassBoundary>,
        max_connections: usize,
    ) -> Result<Self, TransientTransportError> {
        Self::bind(
            cx,
            flow.pressures.len(),
            flow.gravity
                .as_ref()
                .map(|g| g.reference_transport().density_kg_m3()),
            flow.branches
                .iter()
                .map(|branch| (branch.from, branch.to, branch.flow.value())),
            reference_density,
            node_cells,
            thermal,
            boundaries,
            max_connections,
        )
    }

    /// Bind signed reference-volume rates directly to the original graph.
    ///
    /// Inertial owners supply their native step-average rates here; no steady
    /// pressure solution is fabricated for an accelerating fluid. The caller
    /// owns momentum and source admission. This binding independently enforces
    /// the same reference inventory, caloric datum and local mass closure as
    /// [`Self::new`], including explicit external streams at pressure nodes.
    #[allow(clippy::too_many_arguments)]
    pub fn from_volume_flows(
        cx: &Cx<'_>,
        graph: &LossGraph,
        signed_volume_flows_m3_s: &[f64],
        reference_density: Density,
        node_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
        boundaries: Vec<MassBoundary>,
        max_connections: usize,
    ) -> Result<Self, TransientTransportError> {
        poll(cx)?;
        if signed_volume_flows_m3_s.len() != graph.branches().len() {
            return Err(TransientTransportError::Invalid(
                "one signed volume flow per original graph branch is required",
            ));
        }
        Self::bind(
            cx,
            graph.node_count(),
            graph
                .gravity()
                .map(|g| g.reference_transport().density_kg_m3()),
            graph
                .branches()
                .iter()
                .zip(signed_volume_flows_m3_s)
                .map(|(branch, &flow)| (branch.from, branch.to, flow)),
            reference_density,
            node_cells,
            thermal,
            boundaries,
            max_connections,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn bind(
        cx: &Cx<'_>,
        node_count: usize,
        gravity_reference_density: Option<f64>,
        branch_flows: impl ExactSizeIterator<Item = (usize, usize, f64)>,
        reference_density: Density,
        node_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
        boundaries: Vec<MassBoundary>,
        max_connections: usize,
    ) -> Result<Self, TransientTransportError> {
        poll(cx)?;
        let rho = reference_density.value();
        if !rho.is_finite()
            || rho <= 0.0
            || node_cells.is_empty()
            || node_cells.len() != node_count
            || branch_flows
                .len()
                .checked_add(boundaries.len())
                .is_none_or(|n| n > max_connections)
            || gravity_reference_density.is_some_and(|value| value != rho)
        {
            return Err(TransientTransportError::Invalid(
                "invalid node, density or connection binding",
            ));
        }
        let mut storage = BTreeMap::new();
        for &id in node_cells {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(TransientTransportError::Invalid(
                    "hydraulic node has no active storage cell",
                ))?;
            let expected = rho * cell.geometry.volume_m3;
            if !expected.is_finite()
                || expected <= 0.0
                || (cell.mass_kg - expected).abs()
                    > 64.0 * f64::EPSILON * cell.mass_kg.max(expected)
                || thermal.material_storage_basis(cell.material_id)?
                    != ThermalStorageBasis::Enthalpy
                || storage
                    .insert(
                        id,
                        BoundTransportStorage {
                            material: cell.material_id,
                            mass_kg: cell.mass_kg,
                            geometry: cell.geometry,
                            caloric_identity: thermal
                                .material_caloric_identity(cell.material_id)?,
                        },
                    )
                    .is_some()
            {
                return Err(TransientTransportError::Invalid(
                    "storage needs distinct enthalpy cells with the declared reference density",
                ));
            }
        }
        for boundary in &boundaries {
            poll(cx)?;
            let id = match boundary {
                MassBoundary::Inlet { cell, .. } | MassBoundary::Outlet { cell, .. } => cell,
            };
            if !storage.contains_key(id) {
                return Err(TransientTransportError::Invalid(
                    "mass boundary is outside the hydraulic storage binding",
                ));
            }
        }
        let mut flows = Vec::with_capacity(branch_flows.len());
        for (index, (from, to, volume_flow)) in branch_flows.enumerate() {
            poll(cx)?;
            if from >= node_cells.len()
                || to >= node_cells.len()
                || from == to
                || !volume_flow.is_finite()
            {
                return Err(TransientTransportError::Invalid(
                    "invalid hydraulic branch endpoints or flow",
                ));
            }
            let kg_s = rho * volume_flow;
            if !kg_s.is_finite() || (volume_flow != 0.0 && kg_s == 0.0) {
                return Err(TransientTransportError::Invalid(
                    "mass flow is not representable",
                ));
            }
            flows.push(MassFlow {
                id: u64::try_from(index)
                    .map_err(|_| TransientTransportError::Invalid("branch ID overflow"))?,
                a: node_cells[from],
                b: node_cells[to],
                kg_s,
            });
        }
        let advection = Advected::new((), flows, max_connections)?
            .with_boundaries(boundaries, max_connections)?;
        // Reuse the actual datum-aware upwind law and its roundoff-level
        // continuity test; a loose hydraulic stopping tolerance is not a leak.
        advection.endpoint_transport(thermal, || thermal_checkpoint(cx))?;
        Ok(Self {
            node_cells: node_cells.to_vec(),
            storage,
            advection,
            max_connections,
        })
    }

    /// Node-to-storage mapping, in hydraulic node order.
    pub fn node_cells(&self) -> &[CellId] {
        &self.node_cells
    }

    /// Compose with the native source/conduction observer instead of creating
    /// a parallel time owner. Its existing budgets, combined CFL restriction,
    /// transactional histories and material domains remain authoritative.
    pub fn wrap_observer<O>(&self, inner: O) -> Result<Advected<O>, DepositionError> {
        Advected::new(inner, self.advection.flows().to_vec(), self.max_connections)?
            .with_boundaries(self.advection.boundaries().to_vec(), self.max_connections)
    }

    /// Advance the original thermal state using its existing conservative
    /// explicit stepper. Native boundary/source/contact heating participates
    /// in the SAME step and combined stability restriction as advection.
    /// The energy and advection audits overlap; do not add their energy twice.
    pub fn advance_to<'m>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        config: StepConfig,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, AdvectionAudit), TransientTransportError>
    {
        self.validate_storage(cx, thermal)?;
        let (next, advected, energy) =
            thermal
                .advance_observed(&self.advection, target_s, config, || thermal_checkpoint(cx))?;
        self.validate_storage(cx, &next)?;
        Ok((next, energy, advected.audit()))
    }

    /// Native endpoint coefficients for a caller's coupled implicit residual.
    /// Pairwise transfers remain internal; only actual external streams enter
    /// sources. This evaluates no timestep and commits no state or history.
    pub fn endpoint_transport(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<TransportStepCoefficients<()>, TransientTransportError> {
        self.validate_storage(cx, thermal)?;
        Ok(self
            .advection
            .endpoint_transport(thermal, || thermal_checkpoint(cx))?)
    }

    /// One backward-Euler advection/fixed-source step on at most 64 cells.
    /// Specific ENTHALPY transport is linear even when the original T(h) law
    /// is nonlinear or has a latent plateau. Solve energy increments with the
    /// exact upwind matrix, then gate the ORIGINAL native local/global residual.
    /// Closed adiabatic loops require no artificial wall anchor or reservoir.
    ///
    /// Active conduction and Robin/reservoir exchange refuse here: use the
    /// general endpoint interface with their nonlinear owner instead. Authored
    /// fixed watts are included. Geometry, inventory and caloric laws stay fixed;
    /// timestep/spatial refinement is still necessary despite implicit stability.
    pub fn advance_advection_implicit<'m>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit), TransientTransportError> {
        self.validate_storage(cx, thermal)?;
        let n = thermal.cells().len();
        if n != self.storage.len()
            || n > MAX_RECIRCULATION_SUPPLIES
            || thermal.links().values().any(|l| l.conductance_w_k != 0.0)
            || thermal
                .cells()
                .values()
                .any(|c| c.boundary.htc_w_m2_k != 0.0 || c.boundary.reservoir.is_some())
            || !cell_tolerance_j.is_finite()
            || cell_tolerance_j < 0.0
            || !total_tolerance_j.is_finite()
            || total_tolerance_j < 0.0
        {
            return Err(TransientTransportError::Invalid(
                "direct implicit advection needs at most 64 mapped cells, fixed sources and no active conduction/Robin exchange",
            ));
        }
        let step = BackwardEulerStep::new(thermal, target_s, || thermal_checkpoint(cx))?;
        let ids: BTreeMap<_, _> = step
            .cell_ids()
            .iter()
            .enumerate()
            .map(|(i, &id)| (id, i))
            .collect();
        let mass: Vec<_> = step
            .cell_ids()
            .iter()
            .map(|id| thermal.cells()[id].mass_kg)
            .collect();
        let dt = target_s - thermal.time_s();
        let mut matrix = vec![vec![0.0; n]; n];
        let mut balance = vec![0.0; n];
        let mut component: Vec<_> = (0..n).collect();
        for (i, row) in matrix.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        for flow in self.advection.flows().iter().filter(|f| f.kg_s != 0.0) {
            poll(cx)?;
            let (up, down) = if flow.kg_s > 0.0 {
                (ids[&flow.a], ids[&flow.b])
            } else {
                (ids[&flow.b], ids[&flow.a])
            };
            let rate = flow.kg_s.abs();
            matrix[up][up] += dt * rate / mass[up];
            matrix[down][up] -= dt * rate / mass[up];
            balance[up] -= rate;
            balance[down] += rate;
            let a = root(&component, up);
            let b = root(&component, down);
            component[a.max(b)] = a.min(b);
        }
        for boundary in self.advection.boundaries() {
            poll(cx)?;
            match *boundary {
                MassBoundary::Inlet { cell, kg_s, .. } => balance[ids[&cell]] += kg_s,
                MassBoundary::Outlet { cell, kg_s, .. } => {
                    let i = ids[&cell];
                    balance[i] -= kg_s;
                    matrix[i][i] += dt * kg_s / mass[i];
                }
            }
        }
        // Native advection subtracts the CURRENT lowest-CellId component
        // enthalpy. Retain its exact derivative even for admitted roundoff in
        // mass balance, rather than multiplying a large formation enthalpy.
        for i in 0..n {
            let reference = root(&component, i);
            matrix[i][reference] += dt * balance[i] / mass[reference];
        }
        let mut energy = vec![0.0; n];
        for attempt in 0..3 {
            poll(cx)?;
            let trial = step.evaluate_with_transport(
                &energy,
                |candidate| {
                    self.advection
                        .endpoint_transport(candidate, || thermal_checkpoint(cx))
                },
                || thermal_checkpoint(cx),
            )?;
            if attempt == 2
                || (trial.maximum_cell_residual_j() <= cell_tolerance_j
                    && trial.audit().residual_j.abs() <= total_tolerance_j)
            {
                let (next, (), audit) =
                    trial.accept(cell_tolerance_j, total_tolerance_j, || {
                        thermal_checkpoint(cx)
                    })?;
                return Ok((next, audit));
            }
            let rhs: Vec<_> = trial.residual_j().iter().map(|r| -r).collect();
            let correction = solve_feedback(cx, &matrix, &rhs, false)?;
            for (value, delta) in energy.iter_mut().zip(correction) {
                *value += delta;
            }
            energy = step.resolvable_increments_j(&energy);
        }
        unreachable!("bounded implicit loop returns its final native acceptance")
    }

    fn validate_storage(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), TransientTransportError> {
        for (&id, bound) in &self.storage {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(TransientTransportError::Invalid(
                    "bound fluid storage disappeared; rebuild the transport binding",
                ))?;
            if cell.material_id != bound.material
                || cell.mass_kg != bound.mass_kg
                || cell.geometry != bound.geometry
                || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
            {
                return Err(TransientTransportError::Invalid(
                    "bound fluid geometry, mass or caloric identity changed; rebuild the transport binding",
                ));
            }
        }
        Ok(())
    }
}

fn root(parents: &[usize], mut node: usize) -> usize {
    while parents[node] != node {
        node = parents[node];
    }
    node
}

fn thermal_checkpoint(cx: &Cx<'_>) -> Result<(), DepositionError> {
    poll(cx).map_err(|_| DepositionError::Cancelled)
}

pub mod coupled;
pub mod distributed;
pub mod implicit;
pub mod inertial;
pub mod material;

#[cfg(test)]
mod tests;
