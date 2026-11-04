//! Prepared field queries through the same deterministic source tree used by
//! connected-wake dynamics. Physical source data remain immutable; only the
//! preallocated traversal stack changes between queries.

use super::{
    PreparedSource, Refusal, SourceTree, WakeTreeSpec, WorkMeter, checkpoint, refusal,
};
use crate::filament3d::{
    SegmentGeometry, WakeKinematics, WakeSegment, segment_kinematics_with_core,
};

/// Absolute source-expansion tolerances for one physical field probe.
/// Zero in either channel selects exact original-source-order kinematics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeTreeAccuracy {
    /// Euclidean velocity error [m/s].
    pub absolute_velocity_error_mps: f64,
    /// Frobenius norm of the complete velocity-Jacobian error [1/s].
    pub absolute_gradient_error_per_s: f64,
}

/// Measured work and summed local expansion remainders for one query.
/// These bound source approximation, not finite-difference truncation,
/// accumulated trajectories, or every floating-point operation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WakeTreeQueryReport {
    pub work_used: u64,
    pub visited_cells: u64,
    pub accepted_cells: u64,
    /// Individual retained-source interactions; core lookups count distance tests.
    pub direct_interactions: u64,
    pub velocity_error_bound_mps: f64,
    pub gradient_error_bound_per_s: f64,
}

/// An owned immutable source snapshot with a reusable traversal workspace.
/// Preparation copies source geometry once. Successful queries allocate no
/// source, index, tree, or traversal buffers and never mutate physical sources.
///
/// This uses the advection tree's subdivision and continuous-line moments,
/// rather than a second tree representation. Reflected image queries can reuse
/// the same tree; orthogonal reflection preserves both reported error norms.
pub struct PreparedWakeTree {
    sources: Vec<PreparedSource>,
    tree: SourceTree,
}

impl core::fmt::Debug for PreparedWakeTree {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedWakeTree")
            .field("source_count", &self.source_count())
            .field("node_count", &self.node_count())
            .finish_non_exhaustive()
    }
}

impl PreparedWakeTree {
    /// Conservative requested storage for at most this many input sources,
    /// including source geometry, worst-case binary nodes (leaf capacity one),
    /// source indices, traversal scratch and the owner itself. This excludes
    /// allocator bookkeeping and any excess capacity supplied by the allocator.
    /// No transient source clone or second moment tree is needed.
    #[must_use]
    pub fn required_bytes(max_sources: usize) -> Option<usize> {
        let nodes = if max_sources == 0 {
            0
        } else {
            max_sources
                .checked_next_power_of_two()?
                .checked_mul(2)?
                .checked_sub(1)?
        };
        let sources = max_sources.checked_mul(core::mem::size_of::<PreparedSource>())?;
        let indices = max_sources.checked_mul(core::mem::size_of::<usize>())?;
        let nodes = nodes.checked_mul(core::mem::size_of::<super::Node>())?;
        let stack = (usize::BITS as usize + 1).checked_mul(core::mem::size_of::<usize>())?;
        core::mem::size_of::<Self>()
            .checked_add(sources)?
            .checked_add(indices)?
            .checked_add(nodes)?
            .checked_add(stack)
    }

    /// Admit memory and work before allocation, copy the exact sources, then
    /// build the shared balanced tree. Leaf capacity must be in 1..=256.
    /// Unlike material-advection admission, diagnostic sources may have zero
    /// core; singular target points are masked by the query, not regularized.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        sources: &[WakeSegment],
        leaf_capacity: usize,
        maximum_bytes: usize,
        max_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, WakeTreeQueryReport), Refusal> {
        checkpoint(cancelled)?;
        if !(1..=256).contains(&leaf_capacity) {
            return Err(refusal("wake-field-tree-invalid", "leaf capacity must be in 1..=256"));
        }
        if Self::required_bytes(sources.len()).is_none_or(|bytes| bytes > maximum_bytes) {
            return Err(refusal(
                "wake-field-memory-exhausted",
                "prepared source/tree storage exceeds the explicit byte cap",
            ));
        }
        let mut work = WorkMeter {
            limit: max_work,
            used: 0,
            metered: true,
        };
        work.charge(1)?;
        work.charge(u64::try_from(sources.len()).map_err(|_| allocation())?)?;
        let mut owned = Vec::new();
        owned.try_reserve_exact(sources.len()).map_err(|_| allocation())?;
        for (i, source) in sources.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(cancelled)?;
            }
            let geometry = SegmentGeometry::new(source.start, source.end);
            if !source.circulation.is_finite()
                || source.start.iter().chain(&source.end).any(|v| !v.is_finite())
                || !source.core2_m2.is_finite()
                || source.core2_m2 < 0.0
                || source.core2_m2 > f64::MAX / 16.0
                || !geometry.length.is_finite()
                || (geometry.length == 0.0 && source.start != source.end)
            {
                return Err(refusal(
                    "wake-field-source-invalid",
                    "finite source geometry, circulation and nonnegative representable cores required",
                ));
            }
            owned.push(PreparedSource {
                geometry,
                circulation: source.circulation,
                core2_m2: source.core2_m2,
            });
        }
        let mut tree = SourceTree::default();
        // Query accuracy can change without rebuilding immutable moments.
        // The dynamics-only minimum_rho2 shortcut is not used by field queries.
        tree.rebuild_sources(
            &owned,
            WakeTreeSpec {
                absolute_velocity_error_mps: f64::MAX,
                leaf_capacity,
            },
            true,
            &mut work,
            cancelled,
        )?;
        checkpoint(cancelled)?;
        let report = WakeTreeQueryReport {
            work_used: work.used,
            ..WakeTreeQueryReport::default()
        };
        Ok((Self { sources: owned, tree }, report))
    }

    /// Number of admitted input sources, including zero-strength and degenerate
    /// records. Such records are retained in their original order.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.sources.len()
    }

    #[must_use]
    pub fn node_count(&self) -> usize {
        self.tree.node_count()
    }

    /// Velocity with an absolute Euclidean expansion-error budget.
    /// Zero tolerance sums original sources directly. None denotes a singular
    /// or unrepresentable field; invalid settings, cancellation and exhausted
    /// work return a typed refusal.
    pub fn velocity(
        &mut self,
        p: [f64; 3],
        absolute_error_mps: f64,
        max_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Option<[f64; 3]>, WakeTreeQueryReport), Refusal> {
        let (value, report) = self.query::<false>(
            p,
            WakeTreeAccuracy {
                absolute_velocity_error_mps: absolute_error_mps,
                absolute_gradient_error_per_s: 0.0,
            },
            max_work,
            cancelled,
        )?;
        Ok((value.map(|v| v.velocity), report))
    }

    /// Velocity and the analytic derivative of the same accepted expansion.
    /// Every accepted cell must satisfy both independent error budgets.
    /// The curl error is at most sqrt(2) times the reported Frobenius bound;
    /// divergence error is at most sqrt(3) times that bound.
    pub fn kinematics(
        &mut self,
        p: [f64; 3],
        accuracy: WakeTreeAccuracy,
        max_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Option<WakeKinematics>, WakeTreeQueryReport), Refusal> {
        self.query::<true>(p, accuracy, max_work, cancelled)
    }

    fn query<const GRADIENT: bool>(
        &mut self,
        p: [f64; 3],
        accuracy: WakeTreeAccuracy,
        max_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Option<WakeKinematics>, WakeTreeQueryReport), Refusal> {
        checkpoint(cancelled)?;
        admit_point(p)?;
        if !accuracy.absolute_velocity_error_mps.is_finite()
            || accuracy.absolute_velocity_error_mps < 0.0
            || !accuracy.absolute_gradient_error_per_s.is_finite()
            || accuracy.absolute_gradient_error_per_s < 0.0
        {
            return Err(refusal(
                "wake-field-accuracy-invalid",
                "finite nonnegative velocity and Jacobian tolerances required",
            ));
        }
        let mut work = WorkMeter {
            limit: max_work,
            used: 0,
            metered: true,
        };
        work.charge(1)?;
        let mut report = WakeTreeQueryReport::default();
        let exact = accuracy.absolute_velocity_error_mps == 0.0
            || (GRADIENT && accuracy.absolute_gradient_error_per_s == 0.0);
        let value = self.tree.field_query::<GRADIENT>(
            &self.sources,
            p,
            accuracy,
            exact,
            &mut work,
            &mut report,
            cancelled,
        )?;
        checkpoint(cancelled)?;
        report.work_used = work.used;
        Ok((value, report))
    }

    /// Exact physical core membership: distance to any nonzero-strength segment
    /// is strictly less than radius_multiplier * sqrt(source.core2_m2).
    /// Node balls only exclude impossible matches; leaves use source distances.
    /// Degenerate but nonzero-strength records preserve the direct mask's
    /// point-distance semantics. A zero multiplier always returns false.
    pub fn in_core(
        &mut self,
        p: [f64; 3],
        radius_multiplier: f64,
        max_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(bool, WakeTreeQueryReport), Refusal> {
        checkpoint(cancelled)?;
        admit_point(p)?;
        if !radius_multiplier.is_finite() || radius_multiplier < 0.0 {
            return Err(refusal(
                "wake-field-core-invalid",
                "finite nonnegative core-radius multiplier required",
            ));
        }
        let mut work = WorkMeter {
            limit: max_work,
            used: 0,
            metered: true,
        };
        work.charge(1)?;
        let mut report = WakeTreeQueryReport::default();
        let value = self.tree.field_in_core(
            &self.sources,
            p,
            radius_multiplier,
            &mut work,
            &mut report,
            cancelled,
        )?;
        checkpoint(cancelled)?;
        report.work_used = work.used;
        Ok((value, report))
    }
}

impl SourceTree {
    #[allow(clippy::too_many_arguments)]
    fn field_query<const GRADIENT: bool>(
        &mut self,
        sources: &[PreparedSource],
        p: [f64; 3],
        accuracy: WakeTreeAccuracy,
        exact: bool,
        work: &mut WorkMeter,
        report: &mut WakeTreeQueryReport,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<WakeKinematics>, Refusal> {
        let mut sum = WakeKinematics::default();
        let mut correction = WakeKinematics::default();
        // Circulation-length products may underflow while the source field
        // is still representable. Zero unsigned weight is not proof of no flow.
        let exact = exact
            || self.nodes.first().is_some_and(|node| node.expansion.w_abs_sum == 0.0);
        if exact {
            for (i, source) in sources.iter().enumerate() {
                if i % 256 == 0 {
                    checkpoint(cancelled)?;
                }
                work.charge(if GRADIENT { 4 } else { 1 })?;
                if source.circulation == 0.0 {
                    continue;
                }
                report.direct_interactions += 1;
                let Some(value) = source_field::<GRADIENT>(source, p) else {
                    return Ok(None);
                };
                accumulate_field::<GRADIENT>(&mut sum, &mut correction, value);
            }
            return Ok(finite_field(sum));
        }
        if self.nodes.is_empty() {
            return Ok(Some(sum));
        }
        let total_weight = self.nodes[0].expansion.w_abs_sum;
        self.stack.clear();
        self.stack.push(0);
        while let Some(index) = self.stack.pop() {
            checkpoint(cancelled)?;
            work.charge(8)?;
            report.visited_cells += 1;
            let node = &self.nodes[index];
            if let Some([left, right]) = node.children {
                let fraction = node.expansion.w_abs_sum / total_weight;
                let velocity_budget = (0.5 * accuracy.absolute_velocity_error_mps) * fraction;
                let gradient_budget = (0.5 * accuracy.absolute_gradient_error_per_s) * fraction;
                let velocity_bound = node.expansion.admitted_error_bound(p);
                let gradient_bound = if !GRADIENT {
                    0.0
                } else if velocity_bound.is_finite() && velocity_bound <= velocity_budget {
                    node.expansion.admitted_gradient_error_bound(p)
                } else {
                    f64::INFINITY
                };
                if node.expansion.w_abs_sum > 0.0
                    && velocity_bound.is_finite()
                    && velocity_bound <= velocity_budget
                    && gradient_bound.is_finite()
                    && (!GRADIENT || gradient_bound <= gradient_budget)
                {
                    work.charge(if GRADIENT { 128 } else { 16 })?;
                    let value = if GRADIENT {
                        node.expansion.kinematics(p)
                    } else {
                        WakeKinematics {
                            velocity: node.expansion.velocity(p),
                            gradient: [[0.0; 3]; 3],
                        }
                    };
                    if finite_components(&value) {
                        accumulate_field::<GRADIENT>(&mut sum, &mut correction, value);
                        report.accepted_cells += 1;
                        report.velocity_error_bound_mps += velocity_bound;
                        report.gradient_error_bound_per_s += gradient_bound;
                        continue;
                    }
                    // An unrepresentable polynomial never becomes a false zero:
                    // open the node and let exact sources establish validity.
                }
                self.stack.push(right);
                self.stack.push(left);
            } else {
                for (j, &i) in self.order[node.start..node.end].iter().enumerate() {
                    if j % 256 == 0 {
                        checkpoint(cancelled)?;
                    }
                    work.charge(if GRADIENT { 4 } else { 1 })?;
                    report.direct_interactions += 1;
                    let Some(value) = source_field::<GRADIENT>(&sources[i], p) else {
                        return Ok(None);
                    };
                    accumulate_field::<GRADIENT>(&mut sum, &mut correction, value);
                }
            }
        }
        if !report.velocity_error_bound_mps.is_finite()
            || report.velocity_error_bound_mps > accuracy.absolute_velocity_error_mps
            || !report.gradient_error_bound_per_s.is_finite()
            || (GRADIENT && report.gradient_error_bound_per_s > accuracy.absolute_gradient_error_per_s)
        {
            return Err(refusal(
                "wake-field-error-budget",
                "summed cell remainders exceed the requested source-accuracy budget",
            ));
        }
        Ok(finite_field(sum))
    }

    #[allow(clippy::too_many_arguments)]
    fn field_in_core(
        &mut self,
        sources: &[PreparedSource],
        p: [f64; 3],
        multiplier: f64,
        work: &mut WorkMeter,
        report: &mut WakeTreeQueryReport,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Refusal> {
        if self.nodes.is_empty() || multiplier == 0.0 {
            return Ok(false);
        }
        if !(multiplier * self.nodes[0].expansion.max_source_core2_m2.sqrt()).is_finite() {
            return Err(refusal(
                "wake-field-core-invalid",
                "requested core-distance scale exceeds the numerical domain",
            ));
        }
        self.stack.clear();
        self.stack.push(0);
        while let Some(index) = self.stack.pop() {
            checkpoint(cancelled)?;
            work.charge(8)?;
            report.visited_cells += 1;
            let node = &self.nodes[index];
            let e = &node.expansion;
            let offset: [f64; 3] = core::array::from_fn(|i| p[i] - e.centroid[i]);
            let distance = offset[0].hypot(offset[1]).hypot(offset[2]);
            let core_radius = multiplier * e.max_source_core2_m2.sqrt();
            // Widen exclusion balls for coordinate subtraction and bounding
            // radius roundoff. Uncertain comparisons descend to exact leaves.
            let coordinate_scale = p.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
                + e.centroid.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
                + e.radius_m
                + core_radius;
            let padding = 64.0 * f64::EPSILON * coordinate_scale;
            if distance.is_finite()
                && padding.is_finite()
                && distance - e.radius_m > core_radius + padding
            {
                continue;
            }
            if let Some([left, right]) = node.children {
                self.stack.push(right);
                self.stack.push(left);
            } else {
                for (j, &i) in self.order[node.start..node.end].iter().enumerate() {
                    if j % 256 == 0 {
                        checkpoint(cancelled)?;
                    }
                    work.charge(1)?;
                    report.direct_interactions += 1;
                    let source = &sources[i];
                    let distance = distance_to_source(source, p).ok_or_else(|| {
                        refusal("wake-field-distance-invalid", "source distance is not representable")
                    })?;
                    if distance < multiplier * fs_math::det::sqrt(source.core2_m2) {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }
}

fn source_field<const GRADIENT: bool>(
    source: &PreparedSource,
    p: [f64; 3],
) -> Option<WakeKinematics> {
    if source.circulation == 0.0 || source.geometry.length == 0.0 {
        return Some(WakeKinematics::default());
    }
    let mut value = if GRADIENT {
        // Reuse the prepared owner's numerical separation admission without
        // evaluating the velocity kernel a second time.
        source.admit_separation(p).ok()?;
        let analytic = segment_kinematics_with_core(
            p,
            source.geometry.start,
            source.geometry.end,
            source.core2_m2,
        )?;
        WakeKinematics {
            velocity: analytic.velocity,
            gradient: analytic.gradient,
        }
    } else {
        // The compatibility velocity kernel defines zero on a zero-core axis.
        // Preserve the analytic owner's singular mask for diagnostic probes.
        if source.core2_m2 == 0.0
            && segment_kinematics_with_core(
                p,
                source.geometry.start,
                source.geometry.end,
                0.0,
            )
            .is_none()
        {
            return None;
        }
        WakeKinematics {
            velocity: source.unit_velocity(p).ok()?,
            gradient: [[0.0; 3]; 3],
        }
    };
    for i in 0..3 {
        value.velocity[i] *= source.circulation;
        if GRADIENT {
            for j in 0..3 {
                value.gradient[i][j] *= source.circulation;
            }
        }
    }
    finite_components(&value).then_some(value)
}

fn distance_to_source(source: &PreparedSource, p: [f64; 3]) -> Option<f64> {
    // Preserve the native direct-field mask's arithmetic, including endpoint
    // clamping, degenerate edges and deterministic square root.
    let a = source.geometry.start;
    let ab: [f64; 3] = core::array::from_fn(|i| source.geometry.end[i] - a[i]);
    let ap: [f64; 3] = core::array::from_fn(|i| p[i] - a[i]);
    let length2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let fraction = if length2 == 0.0 {
        0.0
    } else {
        ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / length2).clamp(0.0, 1.0)
    };
    let q: [f64; 3] = core::array::from_fn(|i| a[i] + fraction * ab[i]);
    let residual: [f64; 3] = core::array::from_fn(|i| p[i] - q[i]);
    let distance = fs_math::det::sqrt(
        residual[0] * residual[0] + residual[1] * residual[1] + residual[2] * residual[2],
    );
    distance.is_finite().then_some(distance)
}

fn accumulate_field<const GRADIENT: bool>(
    sum: &mut WakeKinematics,
    correction: &mut WakeKinematics,
    value: WakeKinematics,
) {
    for i in 0..3 {
        kahan(&mut sum.velocity[i], &mut correction.velocity[i], value.velocity[i]);
        if GRADIENT {
            for j in 0..3 {
                kahan(
                    &mut sum.gradient[i][j],
                    &mut correction.gradient[i][j],
                    value.gradient[i][j],
                );
            }
        }
    }
}

#[inline]
fn kahan(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}

fn finite_components(value: &WakeKinematics) -> bool {
    value
        .velocity
        .iter()
        .chain(value.gradient.iter().flatten())
        .all(|v| v.is_finite())
}

fn finite_field(value: WakeKinematics) -> Option<WakeKinematics> {
    (finite_components(&value)
        && value.vorticity().iter().all(|v| v.is_finite())
        && value.divergence().is_finite())
    .then_some(value)
}

fn admit_point(p: [f64; 3]) -> Result<(), Refusal> {
    if p.iter().any(|v| !v.is_finite()) {
        Err(refusal("wake-field-query-invalid", "finite probe coordinates required"))
    } else {
        Ok(())
    }
}

fn allocation() -> Refusal {
    refusal("wake-field-allocation", "prepared field storage is not representable")
}
