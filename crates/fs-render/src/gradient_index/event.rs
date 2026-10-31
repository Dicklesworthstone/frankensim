//! Curved propagation to a geometrically localized boundary.
//!
//! A true empty-ball radius bounds safe *arc length*, since displacement from
//! the step origin cannot exceed traveled length. Consequently thin surfaces
//! cannot be skipped merely because a curved segment's chord misses them.
//! Every internal Runge–Kutta stage is also checked against that ball before
//! evaluating the optical field. This extends distance-based stepping (Hart,
//! 1996, <https://doi.org/10.1007/s003710050084>) to the parent ray ODE; nonlinear
//! tracing with ODE/distance coupling is also described by Seyb et al. (2019),
//! <https://doi.org/10.1145/3355089.3356502>.
//!
//! Geometric clearance and upper proximity are separate obligations. A small
//! implicit residual alone never becomes a hit. Contacts are resolved to the
//! caller's metric tolerance with a transverse approach requirement; grazing
//! ambiguity refuses. When a planar patch is supplied, the incidence-scaled
//! unresolved tail must fit inside its geometric exclusion ball. This also
//! refuses near-edge approaches whose tangent can leave the finite surface.
//! The ODE still has local estimates, not a certified global
//! trajectory enclosure, so these results are numerical contacts, not exact
//! crossing certificates. The returned ray stays at its integrated position.
//!
//! Interface continuation uses an explicitly admitted planar patch or a smooth
//! signed-field neighborhood with a bounded cone of real gradients. Smooth
//! departure evaluates the actual signed chart field at internal stages; its
//! tangent plane is never substituted for the interface. The neighborhood
//! must exclude every other surface, including an opposing sheet of the same
//! object. Departure is integrated, never an uncharged position
//! offset. The newly selected material field must support the contact's small
//! numerical collar on both sides; otherwise its ordinary domain refusal is
//! retained. No entire object is ignored to suppress a self-intersection.
//! Explicitly uniform fields use analytic straight segments limited only by
//! geometry, contact departure and remaining distance, without the ODE step cap.

use super::*;

/// Hard cap on geometrical queries made by one boundary propagation request.
pub const MAX_BOUNDARY_QUERIES: u32 = 1_000_000;

/// Oriented side of a boundary whose normal points toward the positive side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundarySide {
    /// Positive/outward side of the oriented surface.
    Positive,
    /// Negative/inward side of the oriented surface.
    Negative,
    /// The query lies on, or cannot distinguish sides of, the surface.
    OnSurface,
}

impl BoundarySide {
    fn sign(self) -> Option<f64> {
        match self {
            Self::Positive => Some(1.0),
            Self::Negative => Some(-1.0),
            Self::OnSurface => None,
        }
    }
}

/// A geometric proximity witness, separate from a ray-intersection claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryWitness {
    /// Stable scene boundary identity.
    pub boundary_id: usize,
    /// Stable local surface/planar-patch identity within that boundary.
    pub feature_id: usize,
    /// Geometric surface representative; propagation does not teleport here.
    pub surface_point: Point3,
    /// Finite unit geometric normal, oriented toward the positive side.
    pub normal: Vec3,
    /// Valid upper bound on query-point distance to this represented surface.
    pub distance_upper_m: f64,
    /// Query-point side of the oriented represented surface.
    pub side: BoundarySide,
    /// Optional radius of an open ball centered on `surface_point` that
    /// contains only this planar patch and no patch edge or other surface.
    /// This must include other facets of the SAME object and all other objects.
    /// Missing support permits arrival but refuses interface departure.
    pub departure_radius_m: Option<f64>,
}

/// Geometry sampled at one integrated point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryQuery {
    /// Lower bound on distance to EVERY boundary in the oracle. Positive
    /// infinity is allowed when the represented geometry is empty.
    pub clearance_m: f64,
    /// Optional nearest-surface upper proximity and orientation witness.
    /// Indeterminate competing surfaces must be refused by the oracle.
    pub witness: Option<BoundaryWitness>,
}

/// Geometry oracle with an explicit per-path mutable work budget.
///
/// Implementations must preserve geometry/snapshot identity for the complete
/// request. A clearance promise covers all represented surfaces; a witness
/// must not promote an unsupported implicit residual to a distance bound.
pub trait BoundaryOracle {
    /// Query geometric exclusion and, when available, upper proximity.
    fn query(&mut self, position: Point3) -> Result<BoundaryQuery, GradientIndexError>;

    /// Prepare one bounded departure. The default preserves planar witnesses;
    /// smooth backends must prove a regular normal cone throughout the ball.
    fn prepare_departure(
        &mut self,
        request: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        DepartureNeighborhood::planar(request.contact.witness)
    }

    /// Optional stronger arrival admission, using the same regularity owner.
    /// The selected side in `request` is the side toward which the ray travels.
    /// A supplied neighborhood must contain the unresolved incident tail.
    fn admit_arrival(
        &mut self,
        _request: DepartureRequest,
    ) -> Result<Option<DepartureNeighborhood>, GradientIndexError> {
        Ok(None)
    }

    /// Enclose the actual signed metric distance during a smooth departure.
    /// A tangent-plane residual is not an implementation of this method.
    fn signed_departure(
        &mut self,
        _contact: BoundaryWitness,
        _position: Point3,
    ) -> Result<SignedBoundaryDistance, GradientIndexError> {
        Err(GradientIndexError::UnsupportedBoundaryDeparture)
    }
}

/// Geometry needed to admit a selected outgoing direction near an interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DepartureRequest {
    pub contact: GradientDeparture,
    pub origin: Point3,
    pub direction: Vec3,
    pub surface_tolerance_m: f64,
    pub min_approach_cosine: f64,
}

/// One immutable local neighborhood, centered at the contact representative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DepartureNeighborhood {
    /// Existing planar patch, excluding edges and every other surface.
    Planar { radius_m: f64 },
    /// The interface's true oriented unit normals throughout this ball lie
    /// within `max_normal_deviation` of `normal_axis`. The signed field is C1
    /// in the ball; no medial axis or opposing sheet may be hidden inside it.
    /// Every other scene object must also be outside the ball.
    Smooth {
        radius_m: f64,
        normal_axis: Vec3,
        max_normal_deviation: f64,
    },
}

impl DepartureNeighborhood {
    pub fn planar(witness: BoundaryWitness) -> Result<Self, GradientIndexError> {
        let radius_m = witness
            .departure_radius_m
            .ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?;
        Ok(Self::Planar { radius_m })
    }
    fn radius(self) -> f64 {
        match self {
            Self::Planar { radius_m } | Self::Smooth { radius_m, .. } => radius_m,
        }
    }
    fn margin(self, witness: BoundaryWitness, sign: f64, direction: Vec3) -> f64 {
        match self {
            Self::Planar { .. } => sign * witness.normal.dot(direction),
            Self::Smooth {
                normal_axis,
                max_normal_deviation,
                ..
            } => cone_margin(normal_axis, max_normal_deviation, sign, direction),
        }
    }
    fn validate(self) -> Result<(), GradientIndexError> {
        if !self.radius().is_finite() || self.radius() <= 0.0 {
            return Err(GradientIndexError::InvalidBoundaryQuery);
        }
        if let Self::Smooth {
            normal_axis,
            max_normal_deviation,
            ..
        } = self
        {
            if !finite_vec(normal_axis)
                || (normal_axis.dot(normal_axis) - 1.0).abs() > 2e-10
                || !max_normal_deviation.is_finite()
                || !(0.0..1.0).contains(&max_normal_deviation)
            {
                return Err(GradientIndexError::InvalidBoundaryQuery);
            }
        }
        Ok(())
    }
}

/// Rigorous signed-distance bounds in world metres; negative is inside.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SignedBoundaryDistance {
    pub lower_m: f64,
    pub upper_m: f64,
}

impl SignedBoundaryDistance {
    fn selected_lower(self, sign: f64) -> Result<f64, GradientIndexError> {
        if !self.lower_m.is_finite() || !self.upper_m.is_finite() || self.lower_m > self.upper_m {
            return Err(GradientIndexError::InvalidBoundaryQuery);
        }
        Ok(if sign > 0.0 {
            self.lower_m
        } else {
            -self.upper_m
        })
    }
}

/// Continuation from a previously resolved numerical interface contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GradientDeparture {
    /// Previous contact. Its oracle admits a planar or smooth neighborhood.
    pub witness: BoundaryWitness,
    /// Side selected by the reflected or transmitted outgoing direction.
    pub side: BoundarySide,
}

/// Accuracy and work controls for one first-boundary walk.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GradientEventOptions {
    /// Existing optical ODE controls; attempted steps include departure/retries.
    pub integration: GradientIndexOptions,
    /// Maximum admitted upper surface-proximity distance, metres. Must be at
    /// least eight times the integration's absolute position tolerance.
    pub surface_tolerance_m: f64,
    /// Smallest admitted absolute normal/tangent cosine for a contact or
    /// departure. Nearly tangent proximity is an unresolved event, not a hit.
    pub min_approach_cosine: f64,
    /// Maximum oracle calls, including neighborhood preparation, arrival
    /// admission and internal signed-departure stage queries.
    pub max_boundary_queries: u32,
}

impl Default for GradientEventOptions {
    fn default() -> Self {
        Self {
            integration: GradientIndexOptions::default(),
            surface_tolerance_m: 1.0e-7,
            min_approach_cosine: 1.0e-5,
            max_boundary_queries: 65_536,
        }
    }
}

impl GradientEventOptions {
    fn validate(self) -> Result<(), GradientIndexError> {
        self.integration.validate()?;
        if !self.surface_tolerance_m.is_finite()
            || self.surface_tolerance_m <= 0.0
            || self.integration.position_tolerance_m > self.surface_tolerance_m / 8.0
            || !self.min_approach_cosine.is_finite()
            || !(0.0..1.0).contains(&self.min_approach_cosine)
            || self.min_approach_cosine == 0.0
            || self.max_boundary_queries == 0
            || self.max_boundary_queries > MAX_BOUNDARY_QUERIES
        {
            return Err(GradientIndexError::InvalidOptions);
        }
        Ok(())
    }
}

/// Why bounded curved propagation stopped.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoundaryTermination {
    /// Numerically localized transverse contact within the declared tolerance.
    Surface(BoundaryWitness),
    /// The requested arc-length budget ended without a resolved contact.
    DistanceLimit,
}

/// A completed curved segment and its geometric termination.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GradientBoundaryPropagation {
    /// Integrated ray, polarization frame, attenuation and optical path lengths.
    pub propagation: GradientIndexPropagation,
    /// Numerical surface contact or completed distance limit.
    pub termination: BoundaryTermination,
    /// Actual number of geometric oracle calls.
    pub boundary_evaluations: u32,
    /// Largest incidence-scaled contact-collar distance encountered, metres.
    /// This reports numerical localization scale, NOT an interval bound on
    /// the true curved tail or on its omitted material-interface contribution.
    pub contact_resolution_m: f64,
}

/// Propagate to the first localized surface or the arc-length budget.
#[allow(clippy::too_many_arguments)]
pub fn propagate_to_boundary<F: GradientIndexField + ?Sized, O: BoundaryOracle + ?Sized>(
    field: &F,
    oracle: &mut O,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    max_distance_m: f64,
    departure: Option<GradientDeparture>,
    options: GradientEventOptions,
) -> Result<GradientBoundaryPropagation, GradientIndexError> {
    propagate_to_boundary_with_control(
        field,
        oracle,
        ray,
        transverse_axis,
        wavelength_nm,
        max_distance_m,
        departure,
        options,
        || true,
    )
}

/// As [`propagate_to_boundary`], polling before geometry queries and every
/// attempted optical step. A `false` control result cancels transactionally.
#[allow(clippy::too_many_arguments)]
pub fn propagate_to_boundary_with_control<
    F: GradientIndexField + ?Sized,
    O: BoundaryOracle + ?Sized,
>(
    field: &F,
    oracle: &mut O,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    max_distance_m: f64,
    mut departure: Option<GradientDeparture>,
    options: GradientEventOptions,
    mut control: impl FnMut() -> bool,
) -> Result<GradientBoundaryPropagation, GradientIndexError> {
    options.validate()?;
    if !max_distance_m.is_finite() || max_distance_m < 0.0 {
        return Err(GradientIndexError::InvalidOptions);
    }
    if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
        return Err(DielectricError::InvalidWavelength.into());
    }
    if !finite_point(ray.origin) {
        return Err(GradientIndexError::InvalidRay);
    }
    let direction = unit(ray.dir).ok_or(GradientIndexError::InvalidRay)?;
    let axis =
        transverse(direction, transverse_axis).ok_or(GradientIndexError::InvalidTransverseAxis)?;
    let mut contact_resolution = 0.0_f64;
    let mut boundary_evaluations = 0;
    let mut departure_neighborhood = None;
    if let Some(previous) = departure {
        validate_witness(previous.witness)?;
        let sign = previous
            .side
            .sign()
            .ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?;
        if !control() {
            return Err(GradientIndexError::Cancelled);
        }
        charge_query(&mut boundary_evaluations, options.max_boundary_queries)?;
        let neighborhood = oracle.prepare_departure(DepartureRequest {
            contact: previous,
            origin: ray.origin,
            direction,
            surface_tolerance_m: options.surface_tolerance_m,
            min_approach_cosine: options.min_approach_cosine,
        })?;
        neighborhood.validate()?;
        let radius = neighborhood.radius();
        let proximity = distance_upper(ray.origin, previous.witness.surface_point);
        let cosine = neighborhood.margin(previous.witness, sign, direction);
        if proximity > options.surface_tolerance_m
            || radius <= 4.0 * options.surface_tolerance_m
            || cosine < options.min_approach_cosine
        {
            return Err(GradientIndexError::UnsupportedBoundaryDeparture);
        }
        contact_resolution = proximity / cosine;
        departure_neighborhood = Some(neighborhood);
    }
    if !control() {
        return Err(GradientIndexError::Cancelled);
    }
    let initial_sample = field.sample(ray.origin, wavelength_nm)?;
    let uniform = field.is_spatially_uniform();
    if uniform && initial_sample.phase_gradient_per_m != Vec3::new(0.0, 0.0, 0.0) {
        return Err(GradientIndexError::InvalidSample);
    }
    let track_group = initial_sample.group_index.is_some();
    let mut state = [0.0; 12];
    put_vec(&mut state, 3, direction);
    put_vec(&mut state, 6, axis);
    let mut final_sample = initial_sample;
    let mut derivative = rhs(&state, initial_sample, track_group)?;
    let mut distance = 0.0;
    let mut step = options.integration.initial_step_m.min(max_distance_m);
    let mut accepted = 0;
    let mut rejected = 0;
    let mut evaluations = 1;
    let mut max_error = 0.0_f64;
    let mut cached_query = None;

    let termination = loop {
        if !control() {
            return Err(GradientIndexError::Cancelled);
        }
        let point = ray.origin.offset(get_vec(&state, 0));
        let tangent = get_vec(&state, 3);
        if !finite_point(point) {
            return Err(GradientIndexError::NumericalRange);
        }
        if uniform {
            step = max_distance_m - distance;
        }
        // Once clear of the numerical collar, ordinary global clearance
        // resumes. The departed object is never removed from the oracle.
        let mut departing = None;
        if let Some(previous) = departure {
            let sign = previous
                .side
                .sign()
                .ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?;
            let neighborhood =
                departure_neighborhood.ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?;
            let signed_distance = match neighborhood {
                DepartureNeighborhood::Planar { .. } => {
                    sign * previous
                        .witness
                        .normal
                        .dot(point.delta_from(previous.witness.surface_point))
                }
                DepartureNeighborhood::Smooth { .. } => {
                    if !control() {
                        return Err(GradientIndexError::Cancelled);
                    }
                    charge_query(&mut boundary_evaluations, options.max_boundary_queries)?;
                    oracle
                        .signed_departure(previous.witness, point)?
                        .selected_lower(sign)?
                }
            };
            if signed_distance > 2.0 * options.surface_tolerance_m {
                departure = None;
                departure_neighborhood = None;
            } else {
                departing = Some((previous, neighborhood, sign, signed_distance));
            }
        }
        let (center, radius, departure_guard) =
            if let Some((previous, neighborhood, sign, signed_distance)) = departing {
                let radius = neighborhood.radius();
                let cosine = neighborhood.margin(previous.witness, sign, tangent);
                if cosine < options.min_approach_cosine {
                    return Err(GradientIndexError::BoundaryResolutionLimit);
                }
                let available = radius - distance_upper(point, previous.witness.surface_point);
                let target = (3.0 * options.surface_tolerance_m - signed_distance) / cosine;
                step = step.min(target).min(0.875 * available);
                (
                    previous.witness.surface_point,
                    radius,
                    Some((previous, neighborhood, sign, signed_distance)),
                )
            } else {
                let query = if let Some(query) = cached_query {
                    query
                } else {
                    charge_query(&mut boundary_evaluations, options.max_boundary_queries)?;
                    let query = oracle.query(point)?;
                    validate_query(query)?;
                    cached_query = Some(query);
                    query
                };
                if let Some(witness) = query.witness
                    && witness.distance_upper_m <= options.surface_tolerance_m
                {
                    let side = witness
                        .side
                        .sign()
                        .ok_or(GradientIndexError::BoundaryResolutionLimit)?;
                    let approach = -side * witness.normal.dot(tangent);
                    if approach >= options.min_approach_cosine {
                        // A small perpendicular gap need not imply contact with
                        // a finite face: at shallow incidence its tangent may
                        // leave the patch before reaching the supporting plane.
                        // The free ball must contain both the current endpoint
                        // and its unresolved tangent extension. This is a local
                        // contact admission check, not a curved-tail certificate.
                        let mut residual = witness.distance_upper_m / approach;
                        if let Some(radius) = witness.departure_radius_m
                            && (distance_upper(point, witness.surface_point) + residual.next_up())
                                .next_up()
                                >= radius
                        {
                            return Err(GradientIndexError::BoundaryResolutionLimit);
                        }
                        if !control() {
                            return Err(GradientIndexError::Cancelled);
                        }
                        charge_query(&mut boundary_evaluations, options.max_boundary_queries)?;
                        if let Some(neighborhood) = oracle.admit_arrival(DepartureRequest {
                            contact: GradientDeparture {
                                witness,
                                side: if side > 0.0 {
                                    BoundarySide::Negative
                                } else {
                                    BoundarySide::Positive
                                },
                            },
                            origin: point,
                            direction: tangent,
                            surface_tolerance_m: options.surface_tolerance_m,
                            min_approach_cosine: options.min_approach_cosine,
                        })? {
                            neighborhood.validate()?;
                            let margin = neighborhood.margin(witness, -side, tangent);
                            if margin < options.min_approach_cosine {
                                return Err(GradientIndexError::BoundaryResolutionLimit);
                            }
                            residual = (witness.distance_upper_m / margin).next_up();
                            if (distance_upper(point, witness.surface_point) + residual).next_up()
                                >= neighborhood.radius()
                            {
                                return Err(GradientIndexError::BoundaryResolutionLimit);
                            }
                        }
                        contact_resolution = contact_resolution.max(residual);
                        break BoundaryTermination::Surface(witness);
                    }
                    if approach.abs() < options.min_approach_cosine {
                        return Err(GradientIndexError::BoundaryResolutionLimit);
                    }
                }
                step = step.min((0.875 * query.clearance_m).next_down());
                (point, query.clearance_m, None)
            };
        if distance == max_distance_m {
            break BoundaryTermination::DistanceLimit;
        }
        if accepted + rejected >= options.integration.max_steps {
            return Err(GradientIndexError::StepBudgetExceeded);
        }
        let remaining = max_distance_m - distance;
        step = step.min(remaining);
        if !step.is_finite() || step <= 0.0 || distance + step == distance {
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        let mut guard_error = None;
        let mut guard = |candidate: Point3, stage_tangent: Vec3| {
            if distance_upper(candidate, center) > radius {
                return false;
            }
            departure_guard.is_none_or(|(previous, neighborhood, sign, initial_signed)| {
                let Some(stage_tangent) = unit(stage_tangent) else {
                    return false;
                };
                if neighborhood.margin(previous.witness, sign, stage_tangent)
                    < options.min_approach_cosine
                {
                    return false;
                }
                let signed = match neighborhood {
                    DepartureNeighborhood::Planar { .. } => Ok(sign
                        * previous
                            .witness
                            .normal
                            .dot(candidate.delta_from(previous.witness.surface_point))),
                    DepartureNeighborhood::Smooth { .. } => {
                        if !control() {
                            Err(GradientIndexError::Cancelled)
                        } else {
                            charge_query(&mut boundary_evaluations, options.max_boundary_queries)
                                .and_then(|()| oracle.signed_departure(previous.witness, candidate))
                                .and_then(|distance| distance.selected_lower(sign))
                        }
                    }
                };
                match signed {
                    Ok(signed) => {
                        signed >= initial_signed - options.integration.position_tolerance_m
                    }
                    Err(error) => {
                        guard_error = Some(error);
                        false
                    }
                }
            })
        };
        let trial = if uniform {
            let new_distance = if step == remaining {
                max_distance_m
            } else {
                distance + step
            };
            let mut next = state;
            put_vec(&mut next, 0, direction.scale(new_distance));
            next[9] = initial_sample.phase_index * new_distance;
            next[10] = initial_sample.group_index.unwrap_or(0.0) * new_distance;
            next[11] = initial_sample.absorption_per_m * new_distance;
            let next_point = ray.origin.offset(get_vec(&next, 0));
            if next.iter().any(|value| !value.is_finite()) || !finite_point(next_point) {
                return Err(GradientIndexError::NumericalRange);
            }
            if guard(next_point, direction) {
                Ok(Trial {
                    state: next,
                    error: [0.0; 12],
                    projection_error: 0.0,
                    sample: initial_sample,
                })
            } else {
                Err(GradientIndexError::BoundaryClearanceExceeded)
            }
        } else {
            trial_step_guarded(
                field,
                ray.origin,
                wavelength_nm,
                &state,
                derivative,
                step,
                track_group,
                guard,
                &mut evaluations,
            )
        };
        // A cancelled or exhausted geometry owner is not a geometric step
        // rejection. Preserve it before retry logic can halve the ODE step.
        if let Some(error) = guard_error {
            return Err(error);
        }
        let trial = match trial {
            Ok(trial) => trial,
            Err(GradientIndexError::BoundaryClearanceExceeded) => {
                rejected += 1;
                if uniform || step <= options.integration.min_step_m {
                    return Err(GradientIndexError::BoundaryResolutionLimit);
                }
                step = (step * 0.5).max(options.integration.min_step_m);
                continue;
            }
            Err(error) => return Err(error),
        };
        let error = normalized_error(&state, &trial, step, options.integration)?;
        let factor = if error == 0.0 {
            5.0
        } else {
            (0.9 * det::pow(error, -0.2)).clamp(0.2, 5.0)
        };
        if error <= 1.0 {
            let next_point = ray.origin.offset(get_vec(&trial.state, 0));
            if next_point == point {
                return Err(GradientIndexError::NumericalRange);
            }
            state = trial.state;
            final_sample = trial.sample;
            derivative = rhs(&state, final_sample, track_group)?;
            distance = if step == remaining {
                max_distance_m
            } else {
                distance + step
            };
            accepted += 1;
            max_error = max_error.max(error);
            step = (step * factor)
                .min(options.integration.max_step_m)
                .max(options.integration.min_step_m);
            cached_query = None;
        } else {
            rejected += 1;
            if step <= options.integration.min_step_m {
                return Err(GradientIndexError::AccuracyNotReached);
            }
            step = (step * factor.min(0.9)).max(options.integration.min_step_m);
        }
    };
    if !contact_resolution.is_finite() {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(GradientBoundaryPropagation {
        propagation: GradientIndexPropagation {
            ray: Ray {
                origin: ray.origin.offset(get_vec(&state, 0)),
                dir: get_vec(&state, 3),
            },
            transverse_axis: get_vec(&state, 6),
            initial_sample,
            final_sample,
            distance_m: distance,
            optical_path_length_m: state[9],
            group_optical_path_length_m: track_group.then_some(state[10]),
            optical_depth: state[11],
            accepted_steps: accepted,
            rejected_steps: rejected,
            field_evaluations: evaluations,
            max_normalized_error_estimate: max_error,
        },
        termination,
        boundary_evaluations,
        contact_resolution_m: contact_resolution,
    })
}

fn charge_query(used: &mut u32, maximum: u32) -> Result<(), GradientIndexError> {
    if *used >= maximum {
        return Err(GradientIndexError::BoundaryQueryBudgetExceeded);
    }
    *used += 1;
    Ok(())
}

pub(crate) fn cone_margin(axis: Vec3, deviation: f64, sign: f64, tangent: Vec3) -> f64 {
    use fs_ivl::Interval;
    let axis = [axis.x, axis.y, axis.z].map(Interval::point);
    let tangent = [tangent.x, tangent.y, tangent.z].map(Interval::point);
    let dot = axis[0] * tangent[0] + axis[1] * tangent[1] + axis[2] * tangent[2];
    let t = tangent.map(Interval::abs);
    let norm = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
    (Interval::point(sign) * dot - Interval::point(deviation) * norm).lo()
}

fn validate_witness(witness: BoundaryWitness) -> Result<(), GradientIndexError> {
    if !finite_point(witness.surface_point)
        || !finite_vec(witness.normal)
        || (witness.normal.dot(witness.normal) - 1.0).abs() > 2.0e-10
        || !witness.distance_upper_m.is_finite()
        || witness.distance_upper_m < 0.0
        || witness
            .departure_radius_m
            .is_some_and(|radius| !radius.is_finite() || radius <= 0.0)
    {
        return Err(GradientIndexError::InvalidBoundaryQuery);
    }
    Ok(())
}

fn validate_query(query: BoundaryQuery) -> Result<(), GradientIndexError> {
    if query.clearance_m.is_nan() || query.clearance_m < 0.0 {
        return Err(GradientIndexError::InvalidBoundaryQuery);
    }
    if let Some(witness) = query.witness {
        validate_witness(witness)?;
        if query.clearance_m > witness.distance_upper_m {
            return Err(GradientIndexError::InvalidBoundaryQuery);
        }
    }
    Ok(())
}

// Outward-rounded upper bound, including subtraction/division/square/sum and
// final norm rounding. Used for geometry safety only, never for optical length.
fn distance_upper(a: Point3, b: Point3) -> f64 {
    if a == b {
        return 0.0;
    }
    let d = [
        (a.x - b.x).abs().next_up(),
        (a.y - b.y).abs().next_up(),
        (a.z - b.z).abs().next_up(),
    ];
    let scale = d[0].max(d[1]).max(d[2]);
    if !scale.is_finite() {
        return f64::INFINITY;
    }
    let mut sum = 0.0;
    for component in d {
        let relative = (component / scale).next_up();
        sum = (sum + (relative * relative).next_up()).next_up();
    }
    (det::sqrt(sum).next_up() * scale).next_up()
}

#[cfg(test)]
mod tests;
