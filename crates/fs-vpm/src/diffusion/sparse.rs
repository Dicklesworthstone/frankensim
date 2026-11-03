//! Reusable spatial graph for fixed-site compact-support PSE.
//!
//! Preparation sorts radius-sized cells and searches a two-cell halo. Each
//! unordered candidate is tested with the dense owner's exact scaled-distance
//! predicate; accepted edges are then sorted by original particle indices.
//! Forward/reverse traversal therefore retains the dense compact law's pair
//! composition, without evaluating absent interactions. Gaussian tails are
//! never truncated: this owner requires an explicitly compact law.
//!
//! Preparation uses O(N log N + C + E log E) work, with C the nearby candidate
//! visits and E the supported edges. Dense clusters can still have quadratic
//! C and E; independent limits refuse that workload. Positions, areas and law
//! are pinned. A moved/remeshed cloud needs another preparation.
//!
//! Each run evaluates the E half-step relaxation fractions once and reuses
//! them across all steps. There are no per-step allocations or repeated
//! pair-kernel evaluations. Shared exceptional-range circulation arithmetic
//! retains its original deterministic fallback. Finite-cloud, quadrature and
//! splitting-error limitations remain those of the declared PSE law.

use super::{
    PairKernel, PseDiffusion, PseKernel, exchange_fraction, scaled_distance2, validate_area,
    validate_inputs,
};
use crate::{
    VortexParticle, VpmBudget, VpmError, VpmRun, WorkMeter, admit_cap, invalid_scalar, reserved_vec,
    validate_particle,
};
use fs_exec::{AdmittedBudget, BudgetConsumption, Cx};

/// Explicit admission limits for constructing one reusable spatial graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PreparationBudget {
    /// Maximum input particles.
    pub max_particles: usize,
    /// Maximum supported unordered edges, independent of candidate density.
    pub max_edges: usize,
    /// Maximum actual distance tests across the count and fill passes.
    pub max_candidate_visits: u64,
    /// Work ceiling admitted against the ambient cost budget before allocation.
    pub max_work_units: u64,
    /// Maximum logical payload of all preparation vectors alive together.
    pub max_live_bytes: usize,
}

impl PreparationBudget {
    /// Declare all preparation limits explicitly.
    #[must_use]
    pub const fn new(
        max_particles: usize,
        max_edges: usize,
        max_candidate_visits: u64,
        max_work_units: u64,
        max_live_bytes: usize,
    ) -> Self {
        Self {
            max_particles,
            max_edges,
            max_candidate_visits,
            max_work_units,
            max_live_bytes,
        }
    }
}

/// Successful preparation's actual work and peak logical storage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PreparationReport {
    /// Distance tests, including rejected nearby pairs and both graph passes.
    /// The fill pass is skipped when no edge is supported.
    pub candidate_visits: u64,
    /// Exact charged input, cell-probe, search, candidate and sort work.
    pub work_units: u64,
    /// Peak vector payload, including temporary cell keys and final edges.
    /// Allocator metadata and capacity rounding are not claimed.
    pub peak_live_bytes: usize,
    /// The admitted work ceiling and actual ambient consumption.
    pub budget: BudgetConsumption,
}

/// Immutable geometry, positive quadrature areas and canonical compact graph.
///
/// Preparation owns three retained vectors (positions, areas, edges) and one
/// temporary cell vector. A count pass precedes the exact edge reservation, so
/// memory scales with actual E rather than the caller's maximum edge allowance.
/// On a 64-bit target the peak logical payload is 48 N + 16 E bytes; retained
/// payload is 24 N + 16 E. Each reservation is fallible and separately admitted.
#[derive(Debug)]
pub struct PreparedPse {
    positions: Vec<[f64; 2]>,
    areas: Vec<f64>,
    law: PseDiffusion,
    edges: Vec<[usize; 2]>,
    preparation: PreparationReport,
}

impl PreparedPse {
    /// Bind one compact law and construct all supported unordered edges.
    ///
    /// Cell coordinates require finite position/radius within +/-2^50.
    /// The two-cell halo covers quotient rounding near cell boundaries; the
    /// final support decision uses exactly the dense owner's r-squared < 1.
    /// Input circulation must be finite but is not retained in the graph.
    ///
    /// Logical work charges one unit for each input validation/copy, neighbor
    /// probe, binary-search comparison, candidate distance test (including
    /// edge count/fill), bounded heap sift, and extraction swap. Cancellation
    /// is observed within the shared 256-unit checkpoint stride. The explicit
    /// max_work_units ceiling is admitted up front, with actual work reported.
    ///
    /// # Errors
    /// Refuses Gaussian laws, invalid/mismatched inputs, unrepresentable cell
    /// coordinates, any local/ambient cap, cancellation or allocation failure.
    /// No partial graph or physical state is published.
    #[allow(clippy::too_many_lines)]
    pub fn prepare_with_cx(
        cx: &Cx<'_>,
        particles: &[VortexParticle],
        areas_m2: &[f64],
        law: PseDiffusion,
        budget: PreparationBudget,
    ) -> Result<Self, VpmError> {
        let n = particles.len();
        validate_inputs(n, areas_m2, 0.0)?;
        if law.kernel() != PseKernel::CompactWendland {
            return Err(invalid_scalar(
                "prepared diffusion requires compact support",
                law.kernel_width_m(),
                None,
            ));
        }
        admit_cap("particles", n, budget.max_particles)?;
        admit_count(
            "preparation work",
            count(n, "preparation particles")?,
            budget.max_work_units,
        )?;
        let base_bytes = payload(
            n,
            core::mem::size_of::<[f64; 2]>()
                + core::mem::size_of::<f64>()
                + core::mem::size_of::<CellEntry>(),
            "preparation input bytes",
        )?;
        admit_cap("logical live bytes", base_bytes, budget.max_live_bytes)?;
        let mut work = PreparationWork::new(cx, budget.max_work_units)?;
        let mut positions = reserved_vec("pse-prepared-positions", n)?;
        let mut areas = reserved_vec("pse-prepared-areas", n)?;
        let mut cells = reserved_vec("pse-prepared-cells", n)?;
        for (i, particle) in particles.iter().copied().enumerate() {
            work.unit("fs-vpm.pse-prepare-input", || {
                validate_particle(particle, i)?;
                validate_area(areas_m2[i], i)?;
                let cell = [
                    cell_coordinate(particle.pos[0], law.kernel_width_m(), i)?,
                    cell_coordinate(particle.pos[1], law.kernel_width_m(), i)?,
                ];
                positions.push(particle.pos);
                areas.push(areas_m2[i]);
                cells.push(CellEntry { cell, particle: i });
                Ok(())
            })?;
        }
        work.boundary("fs-vpm.pse-prepare-input-complete")?;
        heap_sort(&mut cells, &mut work)?;
        work.boundary("fs-vpm.pse-prepare-cells-complete")?;

        let mut candidate_visits = 0;
        let mut edge_count = 0_usize;
        visit_candidates(
            particles,
            &cells,
            law.kernel_width_m(),
            &mut candidate_visits,
            budget.max_candidate_visits,
            &mut work,
            |_| {
                let next = edge_count.checked_add(1).ok_or(VpmError::PlanOverflow {
                    resource: "prepared edges",
                })?;
                admit_cap("prepared edges", next, budget.max_edges)?;
                edge_count = next;
                Ok(())
            },
        )?;
        work.boundary("fs-vpm.pse-prepare-count-complete")?;
        let peak_live_bytes = base_bytes
            .checked_add(payload(
                edge_count,
                core::mem::size_of::<[usize; 2]>(),
                "preparation edge bytes",
            )?)
            .ok_or(VpmError::PlanOverflow {
                resource: "preparation live bytes",
            })?;
        admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
        if edge_count != 0 {
            let visits = candidate_visits
                .checked_mul(2)
                .ok_or(VpmError::PlanOverflow {
                    resource: "preparation candidate visits",
                })?;
            admit_count("preparation candidate visits", visits, budget.max_candidate_visits)?;
        }
        let mut edges = reserved_vec("pse-prepared-edges", edge_count)?;
        if edge_count != 0 {
            visit_candidates(
                particles,
                &cells,
                law.kernel_width_m(),
                &mut candidate_visits,
                budget.max_candidate_visits,
                &mut work,
                |edge| {
                    let next = edges.len().checked_add(1).ok_or(VpmError::PlanOverflow {
                        resource: "prepared edges",
                    })?;
                    admit_cap("prepared edges", next, edge_count)?;
                    edges.push(edge);
                    Ok(())
                },
            )?;
        }
        if edges.len() != edge_count {
            return Err(VpmError::LengthMismatch {
                field: "prepared edge inventory",
                expected: edge_count,
                actual: edges.len(),
            });
        }
        heap_sort(&mut edges, &mut work)?;
        work.boundary("fs-vpm.pse-prepare-publication")?;
        let consumption = work.accountant.consumption();
        if work.meter.completed != consumption.cost_charged {
            return Err(VpmError::AccountingMismatch {
                planned: work.meter.completed,
                completed: work.meter.completed,
                charged: consumption.cost_charged,
            });
        }
        Ok(Self {
            positions,
            areas,
            law,
            edges,
            preparation: PreparationReport {
                candidate_visits,
                work_units: work.meter.completed,
                peak_live_bytes,
                budget: consumption,
            },
        })
    }

    /// Original particle positions; bit patterns must match every run input.
    #[must_use]
    pub fn positions(&self) -> &[[f64; 2]] {
        &self.positions
    }

    /// Original, immutable positive quadrature areas [m²].
    #[must_use]
    pub fn areas_m2(&self) -> &[f64] {
        &self.areas
    }

    /// Original constant-viscosity compact law.
    #[must_use]
    pub const fn law(&self) -> PseDiffusion {
        self.law
    }

    /// Exactly supported unordered pairs, in lexicographic original-index order.
    #[must_use]
    pub fn edges(&self) -> &[[usize; 2]] {
        &self.edges
    }

    /// Actual preparation receipt, separate from all subsequent run budgets.
    #[must_use]
    pub const fn preparation(&self) -> &PreparationReport {
        &self.preparation
    }

    /// Diffuse any finite circulation inventory on these exact fixed sites.
    ///
    /// Each active run charges N + E + S(2E + 1) logical units and 2 E S pair
    /// visits. Its owned payload is N particles plus E f64 fractions; the
    /// immutable borrowed graph is excluded. Zero steps, dt or viscosity skip
    /// coefficient/pair work and fraction allocation, charging N + S units.
    /// Signed-zero position changes also refuse: geometry is pinned bitwise.
    ///
    /// # Errors
    /// Refuses mismatched geometry, nonfinite circulation, invalid dt,
    /// local/ambient budgets, cancellation, allocation failure, or a shared
    /// pair update that cannot be represented. Inputs and graph stay unchanged.
    #[allow(clippy::too_many_lines)]
    pub fn diffuse_with_cx(
        &self,
        cx: &Cx<'_>,
        particles: &[VortexParticle],
        dt: f64,
        steps: usize,
        budget: VpmBudget,
    ) -> Result<VpmRun, VpmError> {
        let n = self.positions.len();
        if particles.len() != n {
            return Err(VpmError::LengthMismatch {
                field: "prepared particles",
                expected: n,
                actual: particles.len(),
            });
        }
        validate_inputs(n, &self.areas, dt)?;
        admit_cap("particles", n, budget.max_particles)?;
        admit_cap("steps", steps, budget.max_steps)?;
        let active = steps != 0 && dt > 0.0 && self.law.kinematic_viscosity_m2_s() > 0.0;
        let fraction_count = if active { self.edges.len() } else { 0 };
        let e = count(fraction_count, "prepared run edges")?;
        let s = count(steps, "prepared run steps")?;
        let pair_evaluations = e
            .checked_mul(2)
            .and_then(|v| v.checked_mul(s))
            .ok_or(VpmError::PlanOverflow {
                resource: "prepared pair evaluations",
            })?;
        admit_count("pair evaluations", pair_evaluations, budget.max_pair_evaluations)?;
        let work_units = count(n, "prepared run particles")?
            .checked_add(e)
            .and_then(|v| v.checked_add(pair_evaluations))
            .and_then(|v| v.checked_add(s))
            .ok_or(VpmError::PlanOverflow {
                resource: "prepared diffusion total work",
            })?;
        let peak_live_bytes = payload(
            n,
            core::mem::size_of::<VortexParticle>(),
            "prepared state bytes",
        )?
        .checked_add(payload(
            fraction_count,
            core::mem::size_of::<f64>(),
            "prepared fraction bytes",
        )?)
        .ok_or(VpmError::PlanOverflow {
            resource: "prepared run live bytes",
        })?;
        admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
        let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
        accountant.checkpoint("fs-vpm.pse-prepared-admission", cx)?;
        let mut state = reserved_vec("pse-prepared-state", n)?;
        let mut fractions = reserved_vec("pse-prepared-fractions", fraction_count)?;
        let mut meter = WorkMeter::default();
        for (i, particle) in particles.iter().copied().enumerate() {
            meter.unit(&mut accountant, cx, "fs-vpm.pse-prepared-input", || {
                validate_particle(particle, i)?;
                for axis in 0..2 {
                    if particle.pos[axis].to_bits() != self.positions[i][axis].to_bits() {
                        let field = if axis == 0 {
                            "prepared position.x"
                        } else {
                            "prepared position.y"
                        };
                        return Err(invalid_scalar(field, particle.pos[axis], Some(i)));
                    }
                }
                state.push(particle);
                Ok(())
            })?;
        }
        meter.boundary(&mut accountant, cx, "fs-vpm.pse-prepared-input-complete")?;
        if active {
            let kernel = PairKernel::new(self.law, dt, 0.5);
            for &[i, j] in &self.edges {
                meter.unit(&mut accountant, cx, "fs-vpm.pse-prepared-fraction", || {
                    let fraction =
                        kernel.fraction(state[i], state[j], [self.areas[i], self.areas[j]]);
                    if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                        return Err(VpmError::NonFiniteResult {
                            phase: "fs-vpm.pse-prepared-fraction",
                            particle: i,
                        });
                    }
                    fractions.push(fraction);
                    Ok(())
                })?;
            }
            meter.boundary(&mut accountant, cx, "fs-vpm.pse-prepared-fractions-complete")?;
        }
        let mut completed_pairs = 0_u64;
        for _ in 0..steps {
            if active {
                let forward = self.edges.iter().zip(&fractions);
                let reverse = self.edges.iter().zip(&fractions).rev();
                for (&[i, j], &fraction) in forward.chain(reverse) {
                    meter.unit(&mut accountant, cx, "fs-vpm.pse-prepared-pair", || {
                        completed_pairs =
                            completed_pairs.checked_add(1).ok_or(VpmError::PlanOverflow {
                                resource: "executed prepared pairs",
                            })?;
                        exchange_fraction(&mut state, &self.areas, fraction, i, j)
                    })?;
                }
            }
            meter.unit(&mut accountant, cx, "fs-vpm.pse-prepared-step", || Ok(()))?;
            meter.boundary(&mut accountant, cx, "fs-vpm.pse-prepared-step-complete")?;
        }
        meter.flush(&mut accountant, "fs-vpm.pse-prepared-finalize")?;
        if meter.completed != work_units
            || completed_pairs != pair_evaluations
            || accountant.consumption().cost_charged != work_units
        {
            return Err(VpmError::AccountingMismatch {
                planned: work_units,
                completed: meter.completed,
                charged: accountant.consumption().cost_charged,
            });
        }
        accountant.checkpoint("fs-vpm.pse-prepared-publication", cx)?;
        Ok(VpmRun {
            particles: state,
            steps_completed: steps,
            pair_evaluations,
            work_units,
            peak_live_bytes,
            budget: accountant.consumption(),
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct CellEntry {
    cell: [i64; 2],
    particle: usize,
}

// Quotient rounding is at most 1/8 cell at this bound. A supported pair's
// component separation is below one radius up to arithmetic roundoff, so its
// floored quotients differ by at most two. The halo also covers negative
// boundaries and quotients that underflow to signed zero. We refuse a larger
// coordinate/radius scale instead of silently dropping physical neighbors.
fn cell_coordinate(position: f64, radius: f64, particle: usize) -> Result<i64, VpmError> {
    const LIMIT: f64 = 1_125_899_906_842_624.0;
    let scaled = position / radius;
    if !scaled.is_finite() || scaled.abs() > LIMIT {
        return Err(invalid_scalar(
            "diffusion scaled position",
            scaled,
            Some(particle),
        ));
    }
    // Integer casts truncate toward zero; the admitted magnitude makes both
    // conversions exact. This computes floor without a platform libm call.
    let truncated = scaled as i64;
    let correction = if scaled < truncated as f64 { 1 } else { 0 };
    Ok(truncated - correction)
}

fn count(value: usize, resource: &'static str) -> Result<u64, VpmError> {
    u64::try_from(value).map_err(|_| VpmError::PlanOverflow { resource })
}

fn payload(elements: usize, bytes: usize, resource: &'static str) -> Result<usize, VpmError> {
    elements
        .checked_mul(bytes)
        .ok_or(VpmError::PlanOverflow { resource })
}

fn admit_count(resource: &'static str, required: u64, maximum: u64) -> Result<(), VpmError> {
    if required > maximum {
        return Err(VpmError::CapExceeded {
            resource,
            required: required.into(),
            maximum: maximum.into(),
        });
    }
    Ok(())
}

// The graph's actual cost is data-dependent. This narrow wrapper admits the
// caller's ceiling, enforces it before every bounded operation, and delegates
// cancellation, deadlines and charging to the existing production WorkMeter.
struct PreparationWork<'cx, 'clock> {
    cx: &'cx Cx<'clock>,
    accountant: AdmittedBudget<'clock>,
    meter: WorkMeter,
    limit: u64,
}

impl<'cx, 'clock> PreparationWork<'cx, 'clock> {
    fn new(cx: &'cx Cx<'clock>, limit: u64) -> Result<Self, VpmError> {
        let mut accountant = AdmittedBudget::admit_ambient(cx, limit)?;
        accountant.checkpoint("fs-vpm.pse-prepare-admission", cx)?;
        Ok(Self {
            cx,
            accountant,
            meter: WorkMeter::default(),
            limit,
        })
    }

    fn unit<T>(
        &mut self,
        phase: &'static str,
        operation: impl FnOnce() -> Result<T, VpmError>,
    ) -> Result<T, VpmError> {
        if self.meter.completed >= self.limit {
            self.meter.flush(&mut self.accountant, phase)?;
            return Err(VpmError::CapExceeded {
                resource: "preparation work",
                required: u128::from(self.meter.completed) + 1,
                maximum: self.limit.into(),
            });
        }
        self.meter.unit(&mut self.accountant, self.cx, phase, operation)
    }

    fn boundary(&mut self, phase: &'static str) -> Result<(), VpmError> {
        self.meter.boundary(&mut self.accountant, self.cx, phase)
    }
}

// Total keys (cell, original index) and unique pair indices make this
// allocation-free heapsort deterministic without relying on sort internals.
// One sift unit performs at most two comparisons and one swap.
fn heap_sort<T: Ord>(values: &mut [T], work: &mut PreparationWork<'_, '_>) -> Result<(), VpmError> {
    let len = values.len();
    for root in (0..len / 2).rev() {
        sift_down(values, root, len, work)?;
    }
    for end in (1..len).rev() {
        work.unit("fs-vpm.pse-prepare-sort", || {
            values.swap(0, end);
            Ok(())
        })?;
        sift_down(values, 0, end, work)?;
    }
    Ok(())
}

fn sift_down<T: Ord>(
    values: &mut [T],
    mut root: usize,
    end: usize,
    work: &mut PreparationWork<'_, '_>,
) -> Result<(), VpmError> {
    while end > 1 && root <= (end - 2) / 2 {
        let next = work.unit("fs-vpm.pse-prepare-sort", || {
            let left = 2 * root + 1;
            let right = left + 1;
            let child = if right < end && values[left] < values[right] {
                right
            } else {
                left
            };
            if values[root] < values[child] {
                values.swap(root, child);
                Ok(Some(child))
            } else {
                Ok(None)
            }
        })?;
        match next {
            Some(child) => root = child,
            None => break,
        }
    }
    Ok(())
}

fn lower_bound(
    cells: &[CellEntry],
    target: CellEntry,
    work: &mut PreparationWork<'_, '_>,
) -> Result<usize, VpmError> {
    let mut low = 0;
    let mut high = cells.len();
    while low < high {
        let mid = low + (high - low) / 2;
        let below = work.unit("fs-vpm.pse-prepare-search", || Ok(cells[mid] < target))?;
        if below {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    Ok(low)
}

fn visit_candidates(
    particles: &[VortexParticle],
    cells: &[CellEntry],
    radius: f64,
    visits: &mut u64,
    max_visits: u64,
    work: &mut PreparationWork<'_, '_>,
    mut accept: impl FnMut([usize; 2]) -> Result<(), VpmError>,
) -> Result<(), VpmError> {
    for source in cells {
        let i = source.particle;
        for dx in -2_i64..=2 {
            for dy in -2_i64..=2 {
                let target = work.unit("fs-vpm.pse-prepare-probe", || {
                    Ok(CellEntry {
                        cell: [source.cell[0] + dx, source.cell[1] + dy],
                        particle: i + 1,
                    })
                })?;
                // Including the original index in the search skips j <= i
                // without scanning dense lower-index neighbors a second time.
                let mut at = lower_bound(cells, target, work)?;
                while at < cells.len() && cells[at].cell == target.cell {
                    let j = cells[at].particle;
                    work.unit("fs-vpm.pse-prepare-candidate", || {
                        let next = visits.checked_add(1).ok_or(VpmError::PlanOverflow {
                            resource: "preparation candidate visits",
                        })?;
                        admit_count("preparation candidate visits", next, max_visits)?;
                        *visits = next;
                        if scaled_distance2(particles[i], particles[j], radius) < 1.0 {
                            accept([i, j])?;
                        }
                        Ok(())
                    })?;
                    at += 1;
                }
            }
        }
    }
    Ok(())
}
