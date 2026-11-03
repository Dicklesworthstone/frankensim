//! Deterministic balanced source tree for the actual continuous finite segments.
//! No source merging, particle substitution, pruning, or physical-state mutation.

use super::{
    PreparedSource, Refusal, WakeAdvectionReport, WakeTreeSpec, WorkMeter, checkpoint, refusal,
};
use crate::farfield3d::expansion::Expansion;
use crate::filament3d::WakeSegment;

struct Node {
    start: usize,
    end: usize,
    children: Option<[usize; 2]>,
    expansion: Expansion,
    minimum_rho2: f64,
}

#[derive(Default)]
pub(super) struct SourceTree {
    order: Vec<usize>,
    nodes: Vec<Node>,
    stack: Vec<usize>,
}

impl SourceTree {
    pub(super) fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub(super) fn rebuild(
        &mut self,
        sources: &[PreparedSource],
        spec: WakeTreeSpec,
        work: &mut WorkMeter,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        checkpoint(cancelled)?;
        work.charge(sources.len() as u64)?;
        self.order.clear();
        self.nodes.clear();
        self.stack.clear();
        self.order
            .try_reserve(sources.len())
            .map_err(|_| allocation())?;
        self.order.extend(
            sources
                .iter()
                .enumerate()
                .filter_map(|(i, s)| (s.geometry.length > 0.0).then_some(i)),
        );
        if self.order.is_empty() {
            return Ok(());
        }
        let leaves = self
            .order
            .len()
            .div_ceil(spec.leaf_capacity)
            .checked_next_power_of_two()
            .ok_or_else(allocation)?;
        let nodes = leaves
            .checked_mul(2)
            .and_then(|v| v.checked_sub(1))
            .ok_or_else(allocation)?;
        self.nodes.try_reserve(nodes).map_err(|_| allocation())?;
        self.stack
            .try_reserve(usize::BITS as usize + 1)
            .map_err(|_| allocation())?;
        self.build_node(
            0,
            self.order.len(),
            sources,
            spec.leaf_capacity,
            work,
            cancelled,
        )?;
        work.charge(self.nodes.len() as u64)?;
        let total_weight = self.nodes[0].expansion.w_abs_sum;
        let log_budget = (0.5 * spec.absolute_velocity_error_mps).ln() - total_weight.ln()
            + core::f64::consts::PI.ln();
        for node in &mut self.nodes {
            // Both retained spatial bounds are at least W*a³/(π*rho⁵).
            // Invert that necessary condition once per stage. Most near nodes
            // can then be opened using only squared distance, before any roots,
            // core-spread remainder or coordinate-conditioning calculations.
            // This gate only rejects approximations; full bounds still accept.
            let threshold = (1.2 * node.expansion.radius_m.ln() - 0.4 * log_budget).exp();
            node.minimum_rho2 = if threshold.is_nan() { 0.0 } else { threshold };
        }
        Ok(())
    }

    fn build_node(
        &mut self,
        start: usize,
        end: usize,
        sources: &[PreparedSource],
        leaf_capacity: usize,
        work: &mut WorkMeter,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize, Refusal> {
        checkpoint(cancelled)?;
        let count = end - start;
        // Moment construction makes bounded source passes. Charge their
        // arithmetic plus a conservative N*ceil(log2 N) sorting allowance.
        let depth = usize::BITS - (count - 1).leading_zeros();
        work.charge(
            (count as u64)
                .checked_mul(16 + u64::from(depth))
                .ok_or_else(allocation)?,
        )?;
        let expansion = Expansion::from_sources(self.order[start..end].iter().map(|&i| {
            let s = &sources[i];
            WakeSegment {
                start: s.geometry.start,
                end: s.geometry.end,
                circulation: s.circulation,
                core2_m2: s.core2_m2,
            }
        }));
        if !expansion.is_finite() {
            return Err(refusal(
                "wake-tree-moment-invalid",
                "tree source moments are not representable",
            ));
        }
        let index = self.nodes.len();
        self.nodes.push(Node {
            start,
            end,
            children: None,
            expansion,
            minimum_rho2: 0.0,
        });
        if count > leaf_capacity {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for &i in &self.order[start..end] {
                let geometry = &sources[i].geometry;
                for c in 0..3 {
                    let mid = 0.5 * geometry.start[c] + 0.5 * geometry.end[c];
                    lo[c] = lo[c].min(mid);
                    hi[c] = hi[c].max(mid);
                }
            }
            let mut axis = 0;
            for c in 1..3 {
                if hi[c] - lo[c] > hi[axis] - lo[axis] {
                    axis = c;
                }
            }
            // A total key gives identical topology even with coincident segment
            // midpoints. No hash seed, random pivots, or scheduling enters it.
            self.order[start..end].sort_unstable_by(|&i, &j| {
                let a = &sources[i].geometry;
                let b = &sources[j].geometry;
                let a = 0.5 * a.start[axis] + 0.5 * a.end[axis];
                let b = 0.5 * b.start[axis] + 0.5 * b.end[axis];
                a.total_cmp(&b).then(i.cmp(&j))
            });
            let middle = start + count / 2;
            let left = self.build_node(start, middle, sources, leaf_capacity, work, cancelled)?;
            let right = self.build_node(middle, end, sources, leaf_capacity, work, cancelled)?;
            self.nodes[index].children = Some([left, right]);
        }
        Ok(index)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn velocity(
        &mut self,
        p: [f64; 3],
        mut velocity: [f64; 3],
        sources: &[PreparedSource],
        spec: WakeTreeSpec,
        work: &mut WorkMeter,
        report: &mut WakeAdvectionReport,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 3], Refusal> {
        if self.nodes.is_empty() {
            return Ok(velocity);
        }
        self.stack.clear();
        self.stack.push(0);
        let total_weight = self.nodes[0].expansion.w_abs_sum;
        let mut correction = [0.0; 3];
        let mut error_sum = 0.0;
        while let Some(index) = self.stack.pop() {
            checkpoint(cancelled)?;
            work.charge(8)?;
            report.visited_cells += 1;
            let node = &self.nodes[index];
            if let Some([left, right]) = node.children {
                let offset: [f64; 3] = core::array::from_fn(|c| p[c] - node.expansion.centroid[c]);
                let rho2 = offset[0] * offset[0]
                    + offset[1] * offset[1]
                    + offset[2] * offset[2]
                    + node.expansion.core2_m2;
                if rho2 < node.minimum_rho2 {
                    self.stack.push(right);
                    self.stack.push(left);
                    continue;
                }
                let bound = node.expansion.admitted_error_bound(p);
                // Unsigned line measure allocates a summable budget across
                // disjoint accepted cells. Half the budget remains as a margin
                // for summing bounds and partition-dependent moment roundoff.
                let local_budget = (0.5 * spec.absolute_velocity_error_mps)
                    * (node.expansion.w_abs_sum / total_weight);
                if bound.is_finite() && bound <= local_budget {
                    work.charge(16)?;
                    let value = node.expansion.velocity(p);
                    accumulate(&mut velocity, &mut correction, value);
                    report.accepted_cells += 1;
                    error_sum += bound;
                    continue;
                }
                self.stack.push(right);
                self.stack.push(left);
            } else {
                work.charge((node.end - node.start) as u64)?;
                report.direct_interactions += (node.end - node.start) as u64;
                for (j, &i) in self.order[node.start..node.end].iter().enumerate() {
                    if j != 0 && j % 256 == 0 {
                        checkpoint(cancelled)?;
                    }
                    let source = &sources[i];
                    let unit = source.unit_velocity(p)?;
                    accumulate(
                        &mut velocity,
                        &mut correction,
                        unit.map(|v| source.circulation * v),
                    );
                }
            }
        }
        if velocity.iter().any(|v| !v.is_finite()) || !error_sum.is_finite() {
            return Err(refusal(
                "wake-advection-nonfinite",
                "tree velocity exceeds the numerical domain",
            ));
        }
        if error_sum > spec.absolute_velocity_error_mps {
            return Err(refusal(
                "wake-tree-error-budget",
                "summed cell bounds exceed the velocity budget",
            ));
        }
        report.max_velocity_error_bound_mps = report.max_velocity_error_bound_mps.max(error_sum);
        Ok(velocity)
    }
}

#[inline]
fn accumulate(sum: &mut [f64; 3], correction: &mut [f64; 3], value: [f64; 3]) {
    for c in 0..3 {
        let increment = value[c] - correction[c];
        let next = sum[c] + increment;
        correction[c] = (next - sum[c]) - increment;
        sum[c] = next;
    }
}

fn allocation() -> Refusal {
    refusal(
        "wake-advection-allocation",
        "tree scratch exceeds the representable allocation domain",
    )
}
