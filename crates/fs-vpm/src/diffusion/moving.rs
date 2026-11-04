//! Bounded compact diffusion workspace for the original moving owner.
//!
//! Geometry is invalidated by every advection or candidate reset. The shared
//! canonical graph builder and original PSE exchange keep direct-path ordering.
//! This removes unsupported diffusion visits; advection remains direct O(N²).

use super::adaptive::{AdaptiveBudget, admit_count, overflow};
use super::sparse::{CellEntry, SpatialWork, cell_coordinate, heap_sort, visit_candidates};
use super::{PairKernel, PseDiffusion, exchange_fraction_observed, pair_work};
use crate::{AdmittedBudget, Cx, VortexParticle, VpmError, WorkMeter, admit_cap, reserved_vec};

/// Storage and aggregate search limits for a moving compact diffusion window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompactNeighborBudget {
    /// Maximum supported unordered edges at any one moved geometry.
    pub max_edges: usize,
    /// Maximum candidate distance tests over all rebuilds and rejected trials.
    pub max_candidate_visits: u64,
}

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct MovingStats {
    pub(crate) graph_rebuilds: usize,
    pub(crate) candidate_visits: u64,
    pub(crate) fraction_evaluations: u64,
    pub(crate) graph_work_units: u64,
}

struct Workspace {
    cells: Vec<CellEntry>,
    edges: Vec<[usize; 2]>,
    fractions: Vec<f64>,
    edge_capacity: usize,
    neighbors: CompactNeighborBudget,
    budget: AdaptiveBudget,
    geometry_valid: bool,
    stats: MovingStats,
}

enum Method {
    Dense,
    Compact(Workspace),
}

/// Internal policy shared by the original viscous stages and adaptive driver.
pub(crate) struct MovingDiffusion {
    method: Method,
}

fn edge_capacity(n: usize, neighbors: CompactNeighborBudget) -> usize {
    let all_pairs = (n as u128) * (n.saturating_sub(1) as u128) / 2;
    (neighbors.max_edges as u128).min(all_pairs) as usize
}

pub(crate) fn workspace_bytes(
    n: usize,
    neighbors: CompactNeighborBudget,
) -> Result<usize, VpmError> {
    n.checked_mul(core::mem::size_of::<CellEntry>())
        .and_then(|bytes| {
            edge_capacity(n, neighbors)
                .checked_mul(core::mem::size_of::<[usize; 2]>() + core::mem::size_of::<f64>())
                .and_then(|edges| bytes.checked_add(edges))
        })
        .ok_or(overflow("moving compact workspace bytes"))
}

impl MovingDiffusion {
    pub(crate) const fn dense() -> Self {
        Self {
            method: Method::Dense,
        }
    }

    // The complete workspace payload is admitted by the owner before this
    // fallible reservation. Every subsequent push stays within these bounds.
    pub(crate) fn compact(
        n: usize,
        neighbors: CompactNeighborBudget,
        budget: AdaptiveBudget,
    ) -> Result<Self, VpmError> {
        let capacity = edge_capacity(n, neighbors);
        Ok(Self {
            method: Method::Compact(Workspace {
                cells: reserved_vec("pse-moving-cells", n)?,
                edges: reserved_vec("pse-moving-edges", capacity)?,
                fractions: reserved_vec("pse-moving-fractions", capacity)?,
                edge_capacity: capacity,
                neighbors,
                budget,
                geometry_valid: false,
                stats: MovingStats::default(),
            }),
        })
    }

    pub(crate) fn invalidate(&mut self) {
        if let Method::Compact(workspace) = &mut self.method {
            workspace.geometry_valid = false;
        }
    }

    pub(crate) fn stats(&self) -> MovingStats {
        match &self.method {
            Method::Dense => MovingStats::default(),
            Method::Compact(workspace) => workspace.stats,
        }
    }

    pub(crate) fn admit_work(
        &self,
        units: u64,
        meter: &mut WorkMeter,
        accountant: &mut AdmittedBudget<'_>,
    ) -> Result<(), VpmError> {
        match &self.method {
            Method::Dense => Ok(()),
            Method::Compact(workspace) => admit_stage(workspace, 0, 0, units, meter, accountant),
        }
    }

    pub(crate) fn admit_advection(
        &self,
        n: usize,
        completed_pairs: u64,
        meter: &mut WorkMeter,
        accountant: &mut AdmittedBudget<'_>,
    ) -> Result<(), VpmError> {
        match &self.method {
            Method::Dense => Ok(()),
            Method::Compact(workspace) => {
                let visits = pair_work(n, 1, 1)?;
                admit_stage(
                    workspace,
                    completed_pairs,
                    visits,
                    visits,
                    meter,
                    accountant,
                )
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sweep(
        &mut self,
        state: &mut [VortexParticle],
        areas: &[f64],
        law: PseDiffusion,
        dt: f64,
        fraction: f64,
        completed_pairs: &mut u64,
        meter: &mut WorkMeter,
        accountant: &mut AdmittedBudget<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), VpmError> {
        self.sweep_observed(
            state,
            areas,
            law,
            dt,
            fraction,
            completed_pairs,
            meter,
            accountant,
            cx,
            &mut |_, _, _, _| Ok(()),
        )
    }

    // Dense and moving-compact paths observe the SAME original exchange.
    // Cached fractions and canonical forward/reverse ordering are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sweep_observed(
        &mut self,
        state: &mut [VortexParticle],
        areas: &[f64],
        law: PseDiffusion,
        dt: f64,
        fraction: f64,
        completed_pairs: &mut u64,
        meter: &mut WorkMeter,
        accountant: &mut AdmittedBudget<'_>,
        cx: &Cx<'_>,
        observer: &mut impl FnMut(
            usize,
            [VortexParticle; 2],
            [VortexParticle; 2],
            f64,
        ) -> Result<(), VpmError>,
    ) -> Result<(), VpmError> {
        let Method::Compact(workspace) = &mut self.method else {
            return super::symmetric_sweep_observed(
                state,
                areas,
                law,
                dt,
                fraction,
                completed_pairs,
                meter,
                accountant,
                cx,
                observer,
            );
        };
        if !workspace.geometry_valid {
            workspace.rebuild(state, law.kernel_width_m(), meter, accountant, cx)?;
        }
        let edges = u64::try_from(workspace.edges.len())
            .map_err(|_| overflow("moving compact edge count"))?;
        let visits = edges
            .checked_mul(2)
            .ok_or(overflow("moving compact pair visits"))?;
        let units = edges
            .checked_add(visits)
            .ok_or(overflow("moving compact sweep work"))?;
        admit_stage(
            workspace,
            *completed_pairs,
            visits,
            units,
            meter,
            accountant,
        )?;
        let kernel = PairKernel::new(law, dt, fraction);
        workspace.fractions.clear();
        for &[i, j] in &workspace.edges {
            meter.unit(accountant, cx, "fs-vpm.pse-moving-fraction", || {
                let value = kernel.fraction(state[i], state[j], [areas[i], areas[j]]);
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    return Err(VpmError::NonFiniteResult {
                        phase: "fs-vpm.pse-moving-fraction",
                        particle: i,
                    });
                }
                workspace.fractions.push(value);
                Ok(())
            })?;
        }
        workspace.stats.fraction_evaluations = workspace
            .stats
            .fraction_evaluations
            .checked_add(edges)
            .ok_or(overflow("moving compact kernel evaluations"))?;
        let forward = workspace.edges.iter().zip(&workspace.fractions);
        for (&[i, j], &value) in forward.clone().chain(forward.rev()) {
            meter.unit(accountant, cx, "fs-vpm.pse-moving-pair", || {
                *completed_pairs = completed_pairs
                    .checked_add(1)
                    .ok_or(overflow("executed moving compact pairs"))?;
                exchange_fraction_observed(state, areas, value, i, j, observer)
            })?;
        }
        meter.boundary(accountant, cx, "fs-vpm.pse-moving-sweep-complete")
    }
}

fn admit_stage(
    workspace: &Workspace,
    completed_pairs: u64,
    pairs: u64,
    units: u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
) -> Result<(), VpmError> {
    let admitted = (|| {
        let next_pairs = completed_pairs
            .checked_add(pairs)
            .ok_or(overflow("moving compact accumulated pairs"))?;
        admit_count(
            "pair evaluations",
            next_pairs,
            workspace.budget.max_pair_evaluations,
        )?;
        let next_work = meter
            .completed
            .checked_add(units)
            .ok_or(overflow("moving compact accumulated work"))?;
        admit_count("adaptive work", next_work, workspace.budget.max_work_units)
    })();
    if admitted.is_err() {
        meter.flush(accountant, "fs-vpm.pse-moving-refusal")?;
    }
    admitted
}

impl Workspace {
    fn rebuild(
        &mut self,
        state: &[VortexParticle],
        radius: f64,
        meter: &mut WorkMeter,
        accountant: &mut AdmittedBudget<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), VpmError> {
        let before = meter.completed;
        let mut work = SpatialWork {
            cx,
            accountant,
            meter,
            limit: self.budget.max_work_units,
        };
        self.cells.clear();
        self.edges.clear();
        for (particle, value) in state.iter().enumerate() {
            work.unit("fs-vpm.pse-moving-cell", || {
                self.cells.push(CellEntry {
                    cell: [
                        cell_coordinate(value.pos[0], radius, particle)?,
                        cell_coordinate(value.pos[1], radius, particle)?,
                    ],
                    particle,
                });
                Ok(())
            })?;
        }
        heap_sort(&mut self.cells, &mut work)?;
        // Capacity was reserved once; one fill pass suffices at every geometry.
        visit_candidates(
            state,
            &self.cells,
            radius,
            &mut self.stats.candidate_visits,
            self.neighbors.max_candidate_visits,
            &mut work,
            |edge| {
                let next = self
                    .edges
                    .len()
                    .checked_add(1)
                    .ok_or(overflow("moving compact supported edges"))?;
                admit_cap("moving compact edges", next, self.edge_capacity)?;
                self.edges.push(edge);
                Ok(())
            },
        )?;
        heap_sort(&mut self.edges, &mut work)?;
        work.boundary("fs-vpm.pse-moving-graph-complete")?;
        self.stats.graph_work_units = self
            .stats
            .graph_work_units
            .checked_add(work.meter.completed - before)
            .ok_or(overflow("moving compact graph work"))?;
        self.stats.graph_rebuilds = self
            .stats
            .graph_rebuilds
            .checked_add(1)
            .ok_or(overflow("moving compact graph rebuilds"))?;
        self.geometry_valid = true;
        Ok(())
    }
}
