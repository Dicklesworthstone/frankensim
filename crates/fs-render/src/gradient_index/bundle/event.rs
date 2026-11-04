//! One-pass angular-bundle transport to native geometric boundary events.
//!
//! The central ray, optical integrals and both original launch-angle columns
//! share a single adaptive trajectory and geometry walk. Each Hessian query is
//! admitted by the same free-ball or signed-departure guard as the optical ray.
//! Boundary contacts keep their integrated endpoint and numerical resolution;
//! the neighboring rays' different arrival lengths remain the responsibility
//! of the differentiated interface or receiving-plane projection owner.

use super::*;
use crate::gradient_index::event::{
    BoundaryOracle, BoundaryTermination, EventStepper, GradientDeparture, GradientEventOptions,
    propagate_with_stepper,
};

/// One central integration policy plus absolute angular-column tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleEventOptions {
    pub boundary: GradientEventOptions,
    /// Metres per radian of original orthonormal launch-angle perturbation.
    pub position_derivative_tolerance_m: f64,
    /// Unit-direction change per radian of original launch-angle perturbation.
    pub direction_derivative_tolerance: f64,
}

impl Default for BundleEventOptions {
    fn default() -> Self {
        let bundle = BundleOptions::default();
        Self {
            boundary: GradientEventOptions::default(),
            position_derivative_tolerance_m: bundle.position_derivative_tolerance_m,
            direction_derivative_tolerance: bundle.direction_derivative_tolerance,
        }
    }
}

impl BundleEventOptions {
    fn bundle(self) -> BundleOptions {
        BundleOptions {
            integration: self.boundary.integration,
            position_derivative_tolerance_m: self.position_derivative_tolerance_m,
            direction_derivative_tolerance: self.direction_derivative_tolerance,
        }
    }
}

/// Cumulative retained-source bundle and this call's geometric/work result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleBoundaryPropagation {
    /// Original launch, arbitrary seeded columns, and all accumulated integrals.
    pub bundle: RayBundlePropagation,
    /// This call only: distance, integrals, accepted/rejected steps and queries.
    pub segment: GradientIndexPropagation,
    pub termination: BoundaryTermination,
    /// Actual geometry queries during this call, including signed stage checks.
    pub boundary_evaluations: u32,
    /// Numerical incidence-scaled contact collar, not a certified path bound.
    pub contact_resolution_m: f64,
}

/// Continue a retained angular bundle to the first admitted native boundary or
/// the additional arc-length limit, with no central-ray pretrace/retrace.
#[allow(clippy::too_many_arguments)]
pub fn propagate_seeded_bundle_to_boundary<
    F: SecondOrderIndexField + ?Sized,
    O: BoundaryOracle + ?Sized,
>(
    field: &F,
    oracle: &mut O,
    seed: &BundleContinuation,
    max_distance_m: f64,
    departure: Option<GradientDeparture>,
    options: BundleEventOptions,
) -> Result<BundleBoundaryPropagation, BundleError> {
    propagate_seeded_bundle_to_boundary_with_control(
        field,
        oracle,
        seed,
        max_distance_m,
        departure,
        options,
        || true,
    )
}

/// Poll before the first field query, every geometry query/attempted step, and
/// signed-field departure checks. Cancellation leaves the seed unchanged.
/// Step and query budgets apply to this call; cumulative work remains on the
/// returned bundle and per-call work on `segment` for scene-wide budgeting.
#[allow(clippy::too_many_arguments)]
pub fn propagate_seeded_bundle_to_boundary_with_control<
    F: SecondOrderIndexField + ?Sized,
    O: BoundaryOracle + ?Sized,
>(
    field: &F,
    oracle: &mut O,
    seed: &BundleContinuation,
    max_distance_m: f64,
    departure: Option<GradientDeparture>,
    options: BundleEventOptions,
    control: impl FnMut() -> bool,
) -> Result<BundleBoundaryPropagation, BundleError> {
    propagate_in_coordinates(
        field,
        oracle,
        seed,
        max_distance_m,
        departure,
        options,
        control,
        None,
    )
    .map(|(result, _)| result)
}

/// Internal continuous-plane adapter. The physical owner must establish that
/// n is continuous at the original plane and supply its original one-sided
/// gradients/Hessians on both sides. The oracle still proves finite departure
/// geometry. Returned derivatives always refer to equal TOTAL arc length.
/// The boolean reports an endpoint still inside the departure coordinate collar.
#[allow(clippy::too_many_arguments)]
pub(crate) fn propagate_seeded_bundle_across_continuous_plane_with_control<
    F: SecondOrderIndexField + ?Sized,
    O: BoundaryOracle + ?Sized,
>(
    field: &F,
    oracle: &mut O,
    seed: &BundleContinuation,
    max_distance_m: f64,
    departure: GradientDeparture,
    options: BundleEventOptions,
    control: impl FnMut() -> bool,
) -> Result<(BundleBoundaryPropagation, bool), BundleError> {
    let plane = continuous::PlaneCoordinates::new(departure, options.boundary)?;
    propagate_in_coordinates(
        field,
        oracle,
        seed,
        max_distance_m,
        Some(departure),
        options,
        control,
        Some(plane),
    )
}

#[allow(clippy::too_many_arguments)]
fn propagate_in_coordinates<F: SecondOrderIndexField + ?Sized, O: BoundaryOracle + ?Sized>(
    field: &F,
    oracle: &mut O,
    seed: &BundleContinuation,
    max_distance_m: f64,
    departure: Option<GradientDeparture>,
    options: BundleEventOptions,
    mut control: impl FnMut() -> bool,
    plane: Option<continuous::PlaneCoordinates>,
) -> Result<(BundleBoundaryPropagation, bool), BundleError> {
    let bundle_options = options.bundle();
    bundle_options.validate()?;
    let (event, stepper) = propagate_with_stepper(
        oracle,
        seed.ray(),
        seed.transverse_axis(),
        seed.previous().wavelength_nm(),
        max_distance_m,
        departure,
        options.boundary,
        &mut control,
        |ray, axis| BundleStepper::new(field, seed, ray, axis, bundle_options, plane),
    )?;
    let derivatives = stepper.fixed_arc_columns()?;
    Ok((
        BundleBoundaryPropagation {
            bundle: seed.append_segment(event.propagation, derivatives)?,
            segment: event.propagation,
            termination: event.termination,
            boundary_evaluations: event.boundary_evaluations,
            contact_resolution_m: event.contact_resolution_m,
        },
        stepper.plane.is_some(),
    ))
}

struct BundleStepper<'a, F: SecondOrderIndexField + ?Sized> {
    field: &'a F,
    ray: Ray,
    wavelength_nm: f64,
    options: BundleOptions,
    initial: IndexJet,
    jet: IndexJet,
    uniform: bool,
    track_group: bool,
    columns: [AngularDerivative; 2],
    state: BundleState,
    derivative: BundleState,
    plane: Option<continuous::PlaneCoordinates>,
    clock: [f64; 2],
    clock_derivative: [f64; 2],
}

impl<'a, F: SecondOrderIndexField + ?Sized> BundleStepper<'a, F> {
    fn new(
        field: &'a F,
        seed: &BundleContinuation,
        ray: Ray,
        axis: Vec3,
        options: BundleOptions,
        plane: Option<continuous::PlaneCoordinates>,
    ) -> Result<Self, BundleError> {
        let wavelength_nm = seed.previous().wavelength_nm();
        let initial = field.sample_second_order(ray.origin, wavelength_nm)?;
        if initial.optical != seed.sample() {
            return Err(BundleError::InitialSampleMismatch);
        }
        let uniform = field.is_spatially_uniform();
        if uniform
            && (initial.optical.phase_gradient_per_m != Vec3::new(0.0, 0.0, 0.0)
                || initial.hessian != IndexHessian::ZERO)
        {
            return Err(GradientIndexError::InvalidSample.into());
        }
        // A plane-coordinate collar uses its two arc-length variations even
        // for constant optics. The original uniform shortcut resumes on the
        // next ordinary request, without changing its original fast path.
        let uniform = uniform && plane.is_none();
        let track_group = initial.optical.group_index.is_some();
        let mut columns = seed.angular_derivatives();
        // The shared event entry normalizes the central tangent. Differentiate
        // this bounded representation correction as well as the ODE projections.
        let raw_norm = vec_norm(seed.ray().dir);
        for column in &mut columns {
            let v = column.direction_per_rad;
            column.direction_per_rad = sub(v, ray.dir.scale(ray.dir.dot(v))).scale(1.0 / raw_norm);
        }
        let mut state = [0.0; 24];
        write_vec(&mut state, 3, ray.dir);
        write_vec(&mut state, 6, axis);
        for (i, column) in columns.iter().enumerate() {
            write_vec(&mut state, 12 + 6 * i, column.position_m_per_rad);
            write_vec(&mut state, 15 + 6 * i, column.direction_per_rad);
        }
        let (clock, derivative, clock_derivative) = if let Some(coordinates) = plane {
            let clock = coordinates.enter(&mut state, initial)?;
            let (derivative, clock_derivative) = coordinates.rhs(&state, initial, track_group)?;
            (clock, derivative, clock_derivative)
        } else {
            (
                [0.0; 2],
                bundle_rhs(&state, initial, track_group)?,
                [0.0; 2],
            )
        };
        Ok(Self {
            field,
            ray,
            wavelength_nm,
            options,
            initial,
            jet: initial,
            uniform,
            track_group,
            columns,
            derivative,
            state,
            plane,
            clock,
            clock_derivative,
        })
    }

    fn fixed_arc_columns(&self) -> Result<[AngularDerivative; 2], GradientIndexError> {
        if let Some(plane) = self.plane {
            plane.fixed_arc(&self.state, self.jet, self.clock)
        } else {
            Ok(core::array::from_fn(|i| AngularDerivative {
                position_m_per_rad: read_vec(&self.state, 12 + 6 * i),
                direction_per_rad: read_vec(&self.state, 15 + 6 * i),
            }))
        }
    }
}

impl<F: SecondOrderIndexField + ?Sized> EventStepper for BundleStepper<'_, F> {
    type Trial = BundleTrial;

    fn state(&self) -> State {
        optical_state(&self.state)
    }
    fn initial_sample(&self) -> GradientIndexSample {
        self.initial.optical
    }
    fn final_sample(&self) -> GradientIndexSample {
        self.jet.optical
    }
    fn is_uniform(&self) -> bool {
        self.uniform
    }
    fn trial_state(trial: &BundleTrial) -> State {
        optical_state(&trial.state)
    }
    fn trial(
        &self,
        step: f64,
        new_distance: f64,
        guard: &mut impl FnMut(Point3, Vec3) -> bool,
        evaluations: &mut u32,
    ) -> Result<BundleTrial, GradientIndexError> {
        if !self.uniform {
            return bundle_step_in_coordinates(
                self.field,
                self.ray.origin,
                self.wavelength_nm,
                &self.state,
                self.derivative,
                step,
                self.track_group,
                guard,
                evaluations,
                self.plane
                    .map(|plane| (plane, self.clock, self.clock_derivative)),
            );
        }
        let mut next = self.state;
        write_vec(&mut next, 0, self.ray.dir.scale(new_distance));
        next[9] = self.initial.optical.phase_index * new_distance;
        next[10] = self.initial.optical.group_index.unwrap_or(0.0) * new_distance;
        next[11] = self.initial.optical.absorption_per_m * new_distance;
        for (i, column) in self.columns.iter().enumerate() {
            write_vec(
                &mut next,
                12 + 6 * i,
                sub(
                    column.position_m_per_rad,
                    column.direction_per_rad.scale(-new_distance),
                ),
            );
        }
        let point = self.ray.origin.offset(read_vec(&next, 0));
        if next.iter().any(|value| !value.is_finite()) || !finite_point(point) {
            return Err(GradientIndexError::NumericalRange);
        }
        if !guard(point, self.ray.dir) {
            return Err(GradientIndexError::BoundaryClearanceExceeded);
        }
        Ok(BundleTrial {
            state: next,
            error: [0.0; 24],
            projection_error: 0.0,
            column_projection_error: [0.0; 2],
            jet: self.initial,
            plane_clock: None,
        })
    }
    fn normalized_error(
        &self,
        trial: &BundleTrial,
        step: f64,
        _: GradientIndexOptions,
    ) -> Result<f64, GradientIndexError> {
        let mut error = bundle_error(&self.state, trial, step, self.options)?;
        if let Some(plane) = self.plane {
            error =
                error.max(plane.error(&self.state, self.jet, self.clock, trial, self.options)?);
        }
        Ok(error)
    }
    fn accept(&mut self, trial: BundleTrial) -> Result<(), GradientIndexError> {
        self.state = trial.state;
        self.jet = trial.jet;
        if let Some(plane) = self.plane {
            self.clock = trial
                .plane_clock
                .ok_or(GradientIndexError::InvalidSample)?
                .value;
            if plane.has_left_collar(self.ray.origin.offset(read_vec(&self.state, 0))) {
                let columns = self.fixed_arc_columns()?;
                for (i, column) in columns.iter().enumerate() {
                    write_vec(&mut self.state, 12 + 6 * i, column.position_m_per_rad);
                    write_vec(&mut self.state, 15 + 6 * i, column.direction_per_rad);
                }
                self.plane = None;
                self.clock = [0.0; 2];
            }
        }
        if let Some(plane) = self.plane {
            (self.derivative, self.clock_derivative) =
                plane.rhs(&self.state, self.jet, self.track_group)?;
        } else {
            self.derivative = bundle_rhs(&self.state, self.jet, self.track_group)?;
            self.clock_derivative = [0.0; 2];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
