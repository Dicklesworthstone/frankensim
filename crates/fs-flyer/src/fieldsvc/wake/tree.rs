//! Prepared wake acceleration through the original field-service return path.
//! One immutable source tree serves physical and reflected probes. Analytic
//! gradients differentiate the accepted expansion; FD probes retain their own
//! explicitly propagated approximation allowance.

use super::{PreparedWakeField, checkpoint};
use crate::fieldsvc::{
    C_GROUND_IMAGES, C_PHYSICAL_WAKE, FD_H, FieldSampleSet, FieldSourceStateV1,
    GridSpec, sample_field_impl,
};
use crate::{Refusal, refuse};
use fs_vpm::filament3d::{
    PreparedWakeTree, WakeKinematics, WakeTreeAccuracy, WakeTreeQueryReport,
};

/// Actual tree work and maximum approximation allowances over a grid request.
/// Preparation is reported separately. Work uses the VPM weighted logical units;
/// source and node visit counts are reported independently. These are not FLOPs.
/// Error allowances concern the retained physical wake, before addition of the
/// atmosphere; they exclude atmospheric model error and FD truncation error.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WakeTreeFieldReport {
    pub work_used: u64,
    pub visited_cells: u64,
    pub accepted_cells: u64,
    pub direct_interactions: u64,
    /// Maximum Euclidean velocity approximation bound at a grid center [m/s].
    pub maximum_velocity_error_bound_mps: f64,
    /// Maximum Frobenius Jacobian approximation bound at a grid center [1/s].
    pub maximum_gradient_error_bound_per_s: f64,
    /// Maximum propagated wake approximation error in the existing FD
    /// divergence dual [1/s]. This does not bound its finite-difference bias.
    pub maximum_fd_divergence_error_bound_per_s: f64,
}

impl WakeTreeFieldReport {
    /// Curl is a linear contraction of the full Jacobian.
    #[must_use]
    pub fn maximum_vorticity_error_bound_per_s(&self) -> f64 {
        core::f64::consts::SQRT_2 * self.maximum_gradient_error_bound_per_s
    }

    /// Trace is a linear contraction of the full Jacobian.
    #[must_use]
    pub fn maximum_analytic_divergence_error_bound_per_s(&self) -> f64 {
        3.0_f64.sqrt() * self.maximum_gradient_error_bound_per_s
    }
}

/// The original field arrays paired with their requested and achieved wake
/// approximation limits. Singular points retain the original invalid mask.
/// Q, lambda2 and speed gradients are evaluated from the returned Jacobian;
/// these scalar nonlinear invariants do not acquire error certificates here.
#[derive(Clone, Debug)]
pub struct BoundedWakeFieldSample {
    pub sample: FieldSampleSet,
    pub accuracy: WakeTreeAccuracy,
    pub report: WakeTreeFieldReport,
}

/// A prepared acceleration of one borrowed, immutable physical field snapshot.
/// Reuse across requests avoids rebuilding sources, moments or traversal scratch.
pub struct PreparedWakeTreeField<'a> {
    source: &'a PreparedWakeField,
    tree: PreparedWakeTree,
    preparation: WakeTreeQueryReport,
}

impl PreparedWakeField {
    /// Prepare the existing VPM tree once. The byte limit covers additional
    /// tree/source/scratch storage; the borrowed exact snapshot is already owned
    /// by the caller. No image source buffer is created. Rejection preserves
    /// this exact snapshot and publishes no partially prepared tree.
    pub fn prepare_tree(
        &self,
        leaf_capacity: usize,
        maximum_bytes: usize,
        maximum_work: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedWakeTreeField<'_>, Refusal> {
        checkpoint(&mut cancelled)?;
        let (tree, preparation) = PreparedWakeTree::prepare(
            &self.segments,
            leaf_capacity,
            maximum_bytes,
            maximum_work,
            &mut cancelled,
        )
        .map_err(tree_refusal)?;
        checkpoint(&mut cancelled)?;
        Ok(PreparedWakeTreeField {
            source: self,
            tree,
            preparation,
        })
    }
}

impl PreparedWakeTreeField<'_> {
    #[must_use]
    pub fn preparation_report(&self) -> WakeTreeQueryReport {
        self.preparation
    }

    #[must_use]
    pub fn node_count(&self) -> usize {
        self.tree.node_count()
    }

    /// Sample the complete existing return set with one request-wide tree work
    /// budget. Physical and image systems share the requested absolute limits.
    ///
    /// FD velocity tolerances also account for division by the existing FD_H:
    /// their summed divergence approximation allowance fits the requested
    /// gradient tolerance. Different accepted tree partitions cannot silently
    /// masquerade as a physical divergence signal. Zero tolerances use exact
    /// source evaluation. Core classification uses exact source distances with
    /// conservative tree pruning, retaining the original mask definition.
    ///
    /// These are retained-source evaluation bounds with the VPM floating-point
    /// allowance, not directed-rounding, trajectory or turbulence certificates.
    /// Cancellation or exhausted work returns no partial field arrays.
    pub fn sample(
        &mut self,
        state: &FieldSourceStateV1,
        grid: &GridSpec,
        component_mask: u32,
        accuracy: WakeTreeAccuracy,
        maximum_work: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BoundedWakeFieldSample, Refusal> {
        checkpoint(&mut cancelled)?;
        grid.admit()?;
        let mut query = TreeFieldQuery::new(&mut self.tree, accuracy, maximum_work)?;
        let sample = sample_field_impl(
            state,
            grid,
            component_mask,
            Some(self.source),
            Some(&mut query),
            &mut cancelled,
        )?;
        Ok(BoundedWakeFieldSample {
            sample,
            accuracy,
            report: query.report,
        })
    }
}

pub(in crate::fieldsvc) struct TreeFieldQuery<'a> {
    tree: &'a mut PreparedWakeTree,
    accuracy: WakeTreeAccuracy,
    maximum_work: u64,
    report: WakeTreeFieldReport,
    current_fd_error: f64,
}

impl<'a> TreeFieldQuery<'a> {
    fn new(
        tree: &'a mut PreparedWakeTree,
        accuracy: WakeTreeAccuracy,
        maximum_work: u64,
    ) -> Result<Self, Refusal> {
        if [
            accuracy.absolute_velocity_error_mps,
            accuracy.absolute_gradient_error_per_s,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(refuse(
                "field-wake-tree-accuracy",
                "wake velocity and gradient limits must be finite and nonnegative".into(),
                "supply absolute velocity [m/s] and Jacobian [1/s] limits",
            ));
        }
        Ok(Self {
            tree,
            accuracy,
            maximum_work,
            report: WakeTreeFieldReport::default(),
            current_fd_error: 0.0,
        })
    }

    fn remaining(&self) -> u64 {
        self.maximum_work - self.report.work_used
    }

    fn absorb(&mut self, report: WakeTreeQueryReport) -> Result<(), Refusal> {
        self.report.work_used = self
            .report
            .work_used
            .checked_add(report.work_used)
            .filter(|&work| work <= self.maximum_work)
            .ok_or_else(|| {
                refuse(
                    "field-wake-budget",
                    "request-wide wake tree work limit exhausted".into(),
                    "reduce the grid or increase the explicit tree work limit",
                )
            })?;
        self.report.visited_cells += report.visited_cells;
        self.report.accepted_cells += report.accepted_cells;
        self.report.direct_interactions += report.direct_interactions;
        Ok(())
    }

    pub(in crate::fieldsvc) fn kinematics(
        &mut self,
        p: [f64; 3],
        mask: u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<WakeKinematics>, Refusal> {
        self.current_fd_error = 0.0;
        let systems = selected_systems(mask);
        let mut out = WakeKinematics::default();
        if systems == 0 {
            return Ok(Some(out));
        }
        let accuracy = WakeTreeAccuracy {
            absolute_velocity_error_mps: self.accuracy.absolute_velocity_error_mps / systems as f64,
            absolute_gradient_error_per_s: self.accuracy.absolute_gradient_error_per_s / systems as f64,
        };
        let mut velocity_error = 0.0;
        let mut gradient_error = 0.0;
        for (bit, reflection) in systems_with_reflection() {
            if mask & bit == 0 {
                continue;
            }
            let q = core::array::from_fn(|i| reflection[i] * p[i]);
            let remaining = self.remaining();
            let (value, report) = self
                .tree
                .kinematics(q, accuracy, remaining, &mut *cancelled)
                .map_err(tree_refusal)?;
            self.absorb(report)?;
            velocity_error += report.velocity_error_bound_mps;
            gradient_error += report.gradient_error_bound_per_s;
            let Some(value) = value else {
                return Ok(None);
            };
            for i in 0..3 {
                out.velocity[i] += reflection[i] * value.velocity[i];
                for j in 0..3 {
                    out.gradient[i][j] += reflection[i] * reflection[j] * value.gradient[i][j];
                }
            }
        }
        self.report.maximum_velocity_error_bound_mps =
            self.report.maximum_velocity_error_bound_mps.max(velocity_error);
        self.report.maximum_gradient_error_bound_per_s =
            self.report.maximum_gradient_error_bound_per_s.max(gradient_error);
        Ok(out
            .velocity
            .iter()
            .chain(out.gradient.iter().flatten())
            .all(|v| v.is_finite())
            .then_some(out))
    }

    pub(in crate::fieldsvc) fn velocity(
        &mut self,
        p: [f64; 3],
        mask: u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<[f64; 3]>, Refusal> {
        let systems = selected_systems(mask);
        let mut out = [0.0; 3];
        if systems == 0 {
            return Ok(Some(out));
        }
        // Six FD probes, each contributing at most error/(2*FD_H).
        let fd_limit = self.accuracy.absolute_gradient_error_per_s * (FD_H / 3.0);
        let tolerance = self.accuracy.absolute_velocity_error_mps.min(fd_limit) / systems as f64;
        let mut velocity_error = 0.0;
        for (bit, reflection) in systems_with_reflection() {
            if mask & bit == 0 {
                continue;
            }
            let q = core::array::from_fn(|i| reflection[i] * p[i]);
            let remaining = self.remaining();
            let (value, report) = self
                .tree
                .velocity(q, tolerance, remaining, &mut *cancelled)
                .map_err(tree_refusal)?;
            self.absorb(report)?;
            velocity_error += report.velocity_error_bound_mps;
            let Some(value) = value else {
                return Ok(None);
            };
            for i in 0..3 {
                out[i] += reflection[i] * value[i];
            }
        }
        self.current_fd_error += velocity_error / (2.0 * FD_H);
        self.report.maximum_fd_divergence_error_bound_per_s = self
            .report
            .maximum_fd_divergence_error_bound_per_s
            .max(self.current_fd_error);
        Ok(out.iter().all(|v| v.is_finite()).then_some(out))
    }

    pub(in crate::fieldsvc) fn in_core(
        &mut self,
        p: [f64; 3],
        mask: u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Refusal> {
        for (bit, reflection) in systems_with_reflection() {
            if mask & bit == 0 {
                continue;
            }
            let q = core::array::from_fn(|i| reflection[i] * p[i]);
            let remaining = self.remaining();
            let (inside, report) = self
                .tree
                .in_core(q, 4.0, remaining, &mut *cancelled)
                .map_err(tree_refusal)?;
            self.absorb(report)?;
            if inside {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

fn selected_systems(mask: u32) -> u32 {
    (mask & (C_PHYSICAL_WAKE | C_GROUND_IMAGES)).count_ones()
}

fn systems_with_reflection() -> [(u32, [f64; 3]); 2] {
    [
        (C_PHYSICAL_WAKE, [1.0; 3]),
        (C_GROUND_IMAGES, [1.0, 1.0, -1.0]),
    ]
}

fn tree_refusal(error: fs_vpm::filament3d::Refusal) -> Refusal {
    refuse(
        error.code,
        error.message,
        "use admitted wake sources and increase the explicit accuracy, memory or work budget",
    )
}
