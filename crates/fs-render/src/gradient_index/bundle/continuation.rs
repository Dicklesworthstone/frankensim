//! Retained launch-angle coordinates across smooth segments and interface maps.
//!
//! A continuation never invents a new point-source launch at an interface.
//! It preserves the original launch ray, wavelength and angular basis, and
//! carries arbitrary admitted endpoint columns into the same augmented solver.
//! Past path length, phase/group optical lengths and absorption remain charged.
//! Interface maps are supplied by their separate physical owner; this module
//! has no reflection/refraction equations, Fresnel weights or geometry offsets.

use super::*;

/// Validated endpoint state with its immutable physical path history.
/// Construction cannot relocate a ray or replace its original angular basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleContinuation {
    previous: RayBundlePropagation,
    start: BundleSegmentStart,
    sample: GradientIndexSample,
    interface_events: u32,
}

impl RayBundlePropagation {
    /// Continue the same smooth ray family from its retained endpoint.
    /// A discrete interface must instead use its differentiated mapping owner.
    pub fn continuation(&self) -> BundleContinuation {
        BundleContinuation {
            previous: *self,
            start: BundleSegmentStart {
                ray: self.propagation.ray,
                axis: self.propagation.transverse_axis,
                columns: self.derivatives,
                expected_sample: Some(self.propagation.final_sample),
            },
            sample: self.propagation.final_sample,
            interface_events: self.interface_events,
        }
    }
}

impl BundleContinuation {
    /// Admit a physical owner's differential interface map. Its columns must
    /// remain derivatives with respect to the ORIGINAL two launch angles.
    /// This validates representation invariants, not the supplied interface law.
    /// Position is unchanged: no uncharged offset or projection is permitted.
    pub(crate) fn from_interface(
        previous: &RayBundlePropagation,
        outgoing: Ray,
        transverse_axis: Vec3,
        outgoing_sample: GradientIndexSample,
        derivatives: [AngularDerivative; 2],
    ) -> Result<Self, BundleError> {
        if outgoing.origin != previous.propagation.ray.origin || !finite_point(outgoing.origin) {
            return Err(BundleError::InvalidContinuation);
        }
        let direction = unit(outgoing.dir).ok_or(BundleError::InvalidContinuation)?;
        let length = vec_norm(outgoing.dir);
        if !length.is_finite() || (length - 1.0).abs() > 2e-10 {
            return Err(BundleError::InvalidContinuation);
        }
        let axis =
            transverse(direction, transverse_axis).ok_or(BundleError::InvalidContinuation)?;
        let mut columns = derivatives;
        for column in &mut columns {
            let v = column.direction_per_rad;
            let magnitude = vec_norm(v);
            if !finite_vec(column.position_m_per_rad)
                || !finite_vec(v)
                || !magnitude.is_finite()
                || direction.dot(v).abs() > 2e-10 * magnitude.max(1.0)
            {
                return Err(BundleError::InvalidContinuation);
            }
            // Differentiate the bounded final unit-direction correction; do
            // not silently project a genuinely non-tangent input derivative.
            column.direction_per_rad =
                sub(v, direction.scale(direction.dot(v))).scale(1.0 / length);
            if !finite_vec(column.direction_per_rad) {
                return Err(BundleError::NumericalRange);
            }
        }
        Ok(Self {
            previous: *previous,
            start: BundleSegmentStart {
                ray: Ray {
                    origin: outgoing.origin,
                    dir: direction,
                },
                axis,
                columns,
                expected_sample: Some(outgoing_sample),
            },
            sample: outgoing_sample,
            interface_events: previous
                .interface_events
                .checked_add(1)
                .ok_or(BundleError::NumericalRange)?,
        })
    }

    pub const fn previous(&self) -> &RayBundlePropagation {
        &self.previous
    }
    pub const fn ray(&self) -> Ray {
        self.start.ray
    }
    pub const fn transverse_axis(&self) -> Vec3 {
        self.start.axis
    }
    pub const fn sample(&self) -> GradientIndexSample {
        self.sample
    }
    pub const fn angular_derivatives(&self) -> [AngularDerivative; 2] {
        self.start.columns
    }
    pub const fn interface_events(&self) -> u32 {
        self.interface_events
    }

    /// Retain the validated mapped endpoint immediately, without a field query
    /// or a positive-length propagation. All past integrals stay unchanged.
    pub fn at_endpoint(&self) -> RayBundlePropagation {
        let mut result = self.previous;
        result.propagation.ray = self.start.ray;
        result.propagation.transverse_axis = self.start.axis;
        result.propagation.final_sample = self.sample;
        result.derivatives = self.start.columns;
        result.interface_events = self.interface_events;
        result
    }
}

/// Propagate additional physical arc length from a retained smooth/interface
/// endpoint. Step budgets apply to this call; returned work counts are
/// cumulative, so an orchestrator can enforce a total multi-segment budget.
pub fn propagate_seeded_bundle<F: SecondOrderIndexField + ?Sized>(
    field: &F,
    seed: &BundleContinuation,
    additional_distance_m: f64,
    options: BundleOptions,
) -> Result<RayBundlePropagation, BundleError> {
    propagate_seeded_bundle_with_control(field, seed, additional_distance_m, options, || true)
}

/// Poll cancellation before the first next-field query and every attempted
/// step, including rejected trials. Refusal leaves the reusable seed unchanged.
pub fn propagate_seeded_bundle_with_control<F: SecondOrderIndexField + ?Sized>(
    field: &F,
    seed: &BundleContinuation,
    additional_distance_m: f64,
    options: BundleOptions,
    mut control: impl FnMut() -> bool,
) -> Result<RayBundlePropagation, BundleError> {
    let (segment, derivatives) = propagate_segment(
        field,
        seed.start,
        seed.previous.wavelength_nm,
        additional_distance_m,
        options,
        &mut control,
    )?;
    let mut result = seed.at_endpoint();
    let past = result.propagation;
    let distance_m = checked_sum(past.distance_m, segment.distance_m)?;
    if segment.distance_m > 0.0 && distance_m <= past.distance_m {
        return Err(BundleError::NumericalRange);
    }
    result.propagation = GradientIndexPropagation {
        ray: segment.ray,
        transverse_axis: segment.transverse_axis,
        initial_sample: past.initial_sample,
        final_sample: segment.final_sample,
        distance_m,
        optical_path_length_m: checked_sum(
            past.optical_path_length_m,
            segment.optical_path_length_m,
        )?,
        // An unsupported value at a zero-length endpoint contributes no delay.
        // Missing group information on a positive-length leg remains missing
        // throughout every later continuation.
        group_optical_path_length_m: if segment.distance_m == 0.0 {
            past.group_optical_path_length_m
        } else if past.distance_m == 0.0 {
            segment.group_optical_path_length_m
        } else {
            match (
                past.group_optical_path_length_m,
                segment.group_optical_path_length_m,
            ) {
                (Some(a), Some(b)) => Some(checked_sum(a, b)?),
                _ => None,
            }
        },
        optical_depth: checked_sum(past.optical_depth, segment.optical_depth)?,
        accepted_steps: past
            .accepted_steps
            .checked_add(segment.accepted_steps)
            .ok_or(BundleError::NumericalRange)?,
        rejected_steps: past
            .rejected_steps
            .checked_add(segment.rejected_steps)
            .ok_or(BundleError::NumericalRange)?,
        field_evaluations: past
            .field_evaluations
            .checked_add(segment.field_evaluations)
            .ok_or(BundleError::NumericalRange)?,
        max_normalized_error_estimate: past
            .max_normalized_error_estimate
            .max(segment.max_normalized_error_estimate),
    };
    result.derivatives = derivatives;
    Ok(result)
}

fn checked_sum(a: f64, b: f64) -> Result<f64, BundleError> {
    let sum = a + b;
    if !a.is_finite() || !b.is_finite() || a < 0.0 || b < 0.0 || !sum.is_finite() {
        return Err(BundleError::NumericalRange);
    }
    Ok(sum)
}

#[cfg(test)]
mod tests;
