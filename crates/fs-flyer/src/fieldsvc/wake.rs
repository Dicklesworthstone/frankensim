//! Exact connected-wake sources for the native field service. Preparation
//! retains physical edges once; repeated samples allocate no source buffers.

use super::*;
use fs_vpm::farfield3d::FarField;
use fs_vpm::filament3d::{FilamentWake, WakeKinematics, WakeSegment};

pub(super) mod tree;
pub use tree::{BoundedWakeFieldSample, PreparedWakeTreeField, WakeTreeFieldReport};

/// An immutable field snapshot of a real connected near/far wake.
///
/// The physical component includes the sheet's bound closing line. It replaces
/// the legacy payload's entire prescribed horseshoe; requesting that horseshoe
/// as an additional component refuses instead of counting the bound line twice.
/// Optional images use the field service's ground plane z=0. They remain
/// boundary devices and never enter the physical wake or its impulse ledger.
///
/// Every retained far edge uses its own core plus its elapsed diffusion offset.
/// This exact-source diagnostic does not differentiate a far-field expansion,
/// inherit its velocity-error bound, recover pruned sources, or resolve the
/// dynamics inside a viscous vortex tube. Rebuild after the physical wake moves.
#[derive(Clone, Debug)]
pub struct PreparedWakeField {
    segments: Vec<WakeSegment>,
    images: bool,
    maximum_core_m: f64,
    identity: String,
}

impl PreparedWakeField {
    /// Bind admitted near/far owners at one epoch. `maximum_segments` bounds
    /// the retained physical segment buffer before allocation; images require
    /// no second buffer. Cancellation is polled during source copying and
    /// before publication. All physical sources must be above z=0 for images.
    pub fn new(
        near: &FilamentWake,
        far: &FarField,
        images: bool,
        maximum_segments: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        checkpoint(&mut cancelled)?;
        far.validate_kinematics_sources(near).map_err(|error| {
            refuse(
                "field-wake-invalid",
                format!("{}: {}", error.code, error.message),
                "use admitted connected near/far sources from the same shedding schedule",
            )
        })?;
        let growth = far.near_core_growth_m2().expect("validated core evolution");
        let mut count = 0usize;
        near.visit_segments_with_core_growth(growth, |_| count = count.saturating_add(1));
        for cell in &far.cells {
            checkpoint(&mut cancelled)?;
            count = count.saturating_add(cell.segs.len());
        }
        if count > maximum_segments || count == usize::MAX {
            return Err(refuse(
                "field-wake-budget",
                format!("{count} segments exceed limit {maximum_segments}"),
                "raise the explicit source limit or sample a smaller admitted wake",
            ));
        }
        let mut segments = Vec::new();
        segments.try_reserve_exact(count).map_err(|_| {
            refuse(
                "field-wake-allocation",
                format!("cannot reserve {count} physical segments"),
                "reduce the source limit or release memory",
            )
        })?;
        let mut stopped = false;
        near.visit_segments_with_core_growth(growth, |segment| {
            if !stopped {
                if segments.len() % 256 == 0 && cancelled() {
                    stopped = true;
                } else {
                    segments.push(segment);
                }
            }
        });
        if stopped {
            return Err(cancelled_refusal());
        }
        for cell in &far.cells {
            for ((start, end, circulation), core) in cell.segs.iter().zip(&cell.source_core2_m2) {
                if segments.len() % 256 == 0 {
                    checkpoint(&mut cancelled)?;
                }
                segments.push(WakeSegment {
                    start: *start,
                    end: *end,
                    circulation: *circulation,
                    core2_m2: core + cell.diffusion_offset_m2,
                });
            }
        }
        let mut maximum_core_m: f64 = 0.0;
        for (i, segment) in segments.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(&mut cancelled)?;
            }
            if !segment.circulation.is_finite()
                || !segment.core2_m2.is_finite()
                || segment.core2_m2 < 0.0
            {
                return Err(refuse(
                    "field-wake-invalid",
                    "assembled edge circulation or evolved core is not representable".into(),
                    "use representable physical edge strengths and core widths",
                ));
            }
            if images && (segment.start[2] <= 0.0 || segment.end[2] <= 0.0) {
                return Err(refuse(
                    "field-wake-ground-contact",
                    "physical wake reaches or crosses the ground plane".into(),
                    "resolve ground contact before using flat-plane vortex images",
                ));
            }
            maximum_core_m = maximum_core_m.max(det::sqrt(segment.core2_m2));
        }
        let mut bytes = near.digest().into_bytes();
        bytes.extend_from_slice(far.digest().as_bytes());
        bytes.push(u8::from(images));
        let identity = hash_domain("org.frankensim.wf.physical-wake-field.v1", &bytes).to_hex();
        checkpoint(&mut cancelled)?;
        Ok(Self {
            segments,
            images,
            maximum_core_m,
            identity,
        })
    }

    /// Number of stored physical edges; images do not duplicate storage.
    #[must_use]
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    #[must_use]
    pub fn supported_components(&self) -> u32 {
        C_PHYSICAL_WAKE | if self.images { C_GROUND_IMAGES } else { 0 }
    }

    /// Sample the existing velocity, gradient, curl and tensor-invariant return
    /// set, including an independent finite-difference divergence of the same
    /// exact field. Atmosphere comes from `state`; its prescribed horseshoe is
    /// replaced by this snapshot. The source identity binds both inputs.
    ///
    /// Work is bounded before output allocation by eight source visits per
    /// selected physical/image edge per point: one analytic sample, six velocity
    /// probes and one core-mask pass. The bound counts visits, not FLOPs. Queries
    /// poll every 256 visits and at point boundaries; refusal publishes no field.
    /// Exact sampling costs O(points * retained segments), with no far-field
    /// derivative accuracy or wall-clock performance claim.
    pub fn sample(
        &self,
        state: &FieldSourceStateV1,
        grid: &GridSpec,
        component_mask: u32,
        maximum_source_visits: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FieldSampleSet, Refusal> {
        checkpoint(&mut cancelled)?;
        grid.admit()?;
        let systems = (component_mask & (C_PHYSICAL_WAKE | C_GROUND_IMAGES)).count_ones() as u64;
        let visits = (self.segments.len() as u64)
            .checked_mul(grid.n_points() as u64)
            .and_then(|n| n.checked_mul(systems))
            .and_then(|n| n.checked_mul(8));
        if visits.is_none_or(|n| n > maximum_source_visits) {
            return Err(refuse(
                "field-wake-budget",
                format!("source visits {visits:?} exceed limit {maximum_source_visits}"),
                "reduce sample points or raise the explicit source-visit limit",
            ));
        }
        sample_field_impl(state, grid, component_mask, Some(self), None, &mut cancelled)
    }

    pub(super) fn selected_core_radius(&self, mask: u32) -> f64 {
        if mask & (C_PHYSICAL_WAKE | C_GROUND_IMAGES) != 0 {
            self.maximum_core_m
        } else {
            0.0
        }
    }

    pub(super) fn snapshot_id(&self, state: &FieldSourceStateV1) -> String {
        let mut bytes = state.field_source_snapshot_id().into_bytes();
        bytes.extend_from_slice(self.identity.as_bytes());
        hash_domain("org.frankensim.wf.connected-field-source.v1", &bytes).to_hex()
    }

    pub(super) fn kinematics(
        &self,
        p: [f64; 3],
        mask: u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<WakeKinematics>, Refusal> {
        let mut out = WakeKinematics::default();
        for (bit, reflection) in [
            (C_PHYSICAL_WAKE, [1.0; 3]),
            (C_GROUND_IMAGES, [1.0, 1.0, -1.0]),
        ] {
            if mask & bit == 0 {
                continue;
            }
            let q = core::array::from_fn(|i| reflection[i] * p[i]);
            for (n, source) in self.segments.iter().enumerate() {
                if n % 256 == 0 {
                    checkpoint(cancelled)?;
                }
                if source.circulation == 0.0 {
                    continue;
                }
                let Some(value) =
                    segment_kinematics_with_core(q, source.start, source.end, source.core2_m2)
                else {
                    return Ok(None);
                };
                for i in 0..3 {
                    out.velocity[i] += reflection[i] * source.circulation * value.velocity[i];
                    for j in 0..3 {
                        out.gradient[i][j] += reflection[i]
                            * reflection[j]
                            * source.circulation
                            * value.gradient[i][j];
                    }
                }
            }
        }
        Ok(out
            .velocity
            .iter()
            .chain(out.gradient.iter().flatten())
            .all(|v| v.is_finite())
            .then_some(out))
    }

    pub(super) fn velocity(
        &self,
        p: [f64; 3],
        mask: u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<[f64; 3]>, Refusal> {
        let mut out = [0.0; 3];
        for (bit, reflection) in [
            (C_PHYSICAL_WAKE, [1.0; 3]),
            (C_GROUND_IMAGES, [1.0, 1.0, -1.0]),
        ] {
            if mask & bit == 0 {
                continue;
            }
            let q = core::array::from_fn(|i| reflection[i] * p[i]);
            for (n, source) in self.segments.iter().enumerate() {
                if n % 256 == 0 {
                    checkpoint(cancelled)?;
                }
                if source.circulation == 0.0 {
                    continue;
                }
                // The velocity-only compatibility kernel defines zero on an
                // unregularized filament. A derivative probe there is singular,
                // so retain the analytic owner's refusal for that case.
                if source.core2_m2 == 0.0
                    && segment_kinematics_with_core(q, source.start, source.end, 0.0).is_none()
                {
                    return Ok(None);
                }
                let value =
                    segment_velocity_with_core(q, source.start, source.end, source.core2_m2);
                for i in 0..3 {
                    out[i] += reflection[i] * source.circulation * value[i];
                }
            }
        }
        Ok(out.iter().all(|v| v.is_finite()).then_some(out))
    }

    pub(super) fn in_core(
        &self,
        p: [f64; 3],
        mask: u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Refusal> {
        for (bit, reflection) in [(C_PHYSICAL_WAKE, 1.0), (C_GROUND_IMAGES, -1.0)] {
            if mask & bit == 0 {
                continue;
            }
            let q = [p[0], p[1], reflection * p[2]];
            for (n, source) in self.segments.iter().enumerate() {
                if n % 256 == 0 {
                    checkpoint(cancelled)?;
                }
                if source.circulation != 0.0
                    && dist_to_segment(q, source.start, source.end)
                        < 4.0 * det::sqrt(source.core2_m2)
                {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

fn cancelled_refusal() -> Refusal {
    refuse(
        "field-cancelled",
        "field sampling cancelled before publication".into(),
        "retry with a live cancellation scope",
    )
}

pub(super) fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(cancelled_refusal())
    } else {
        Ok(())
    }
}
